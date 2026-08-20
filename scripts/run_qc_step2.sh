#!/bin/bash
# Step 2 QC gate for the filtered BAMs, plus bigwigs and TSS enrichment.
#
#   A. bamPEFragmentSize    -- fragment-size ladder
#   B. plotFingerprint      -- enrichment without needing peaks or an IgG
#   C. bamCoverage          -- CPM bigwigs, duplicates excluded
#   D. computeMatrix TSS    -- +/-2kb profile, heatmap, numeric TSS enrichment score
#   E. multiBamSummary      -- 5kb-bin correlation + PCA across all samples
#
# Duplicates are FLAGGED not removed in the step-1 BAMs, so every tool here is
# told to skip them (--ignoreDuplicates / samFlagExclude 1024). Blacklist and
# chrM were already removed upstream.
#
# Usage:  ./run_qc_step2.sh [group]        # group defaults to cutrun
#         PHASES=CDE ./run_qc_step2.sh     # run only some phases
#
# Resumable: each phase skips if its output already exists.
#
# The bamCoverage worker below intentionally uses separate `local` statements
# per line: bash expands every assignment word in a combined `local a=.. b=..`
# before performing any of them, so a combined form referencing the first
# variable in the second would silently collide across parallel workers.

set -u
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/../config.sh"
cd "$PROJECT_DIR"

GROUP="${1:-cutrun}"
PHASES="${PHASES:-ABCDE}"
THREADS="${THREADS:-16}"
BW_JOBS="${BW_JOBS:-4}"          # concurrent bamCoverage, each gets THREADS/BW_JOBS

source "$CONDA_SH"
conda activate "$BIGWIG_ENV"
for t in bamPEFragmentSize plotFingerprint bamCoverage computeMatrix plotProfile \
         plotHeatmap multiBamSummary plotCorrelation plotPCA; do
    command -v "$t" >/dev/null || { echo "ERROR: $t missing from $BIGWIG_ENV"; exit 1; }
done

IN_DIR="filtered/${GROUP}"
QC="qc/${GROUP}"
BW="${IN_DIR}/bigwigs"
mkdir -p "$QC" "$BW" "$QC/logs"

mapfile -t SAMPLES < <(ls "${IN_DIR}"/*.final.bam 2>/dev/null \
    | xargs -r -n1 basename | sed 's/\.final\.bam$//' | sort)
[[ ${#SAMPLES[@]} -gt 0 ]] || { echo "ERROR: no final BAMs in $IN_DIR"; exit 1; }
BAMS=(); for s in "${SAMPLES[@]}"; do BAMS+=("${IN_DIR}/${s}.final.bam"); done
echo "group=$GROUP  samples=${#SAMPLES[@]}  phases=$PHASES"

# ---------------------------------------------------------------- A. fragments
if [[ "$PHASES" == *A* ]] && [[ ! -s "$QC/fragment_sizes.txt" ]]; then
    echo "[A] fragment-size distributions"
    bamPEFragmentSize -b "${BAMS[@]}" --samplesLabel "${SAMPLES[@]}" \
        --histogram "$QC/fragment_sizes.png" --maxFragmentLength 1000 \
        --table "$QC/fragment_sizes.txt" -p "$THREADS" \
        > "$QC/logs/fragment_sizes.log" 2>&1 \
      && echo "[A] done" || echo "[A] FAILED (see $QC/logs/fragment_sizes.log)"
fi

# -------------------------------------------------------------- B. fingerprint
# one plot per target parsed from the sample name suffix ("<tumor>-<target>").
# PBRM1..PBRM8 collapse to one target: those suffixes are tube numbers for a
# single PBRM1 target, not 8 distinct targets -- edit the sed pattern below if
# your own naming convention has a similar tube/replicate-number suffix.
if [[ "$PHASES" == *B* ]]; then
    echo "[B] fingerprints per target"
    mapfile -t TARGETS < <(printf '%s\n' "${SAMPLES[@]}" \
        | sed 's/.*-//; s/^PBRM[0-9]*$/PBRM1/' | sort -u)
    for tgt in "${TARGETS[@]}"; do
        outp="$QC/fingerprint_${tgt}.png"
        [[ -s "$outp" ]] && continue
        tb=(); tl=()
        for s in "${SAMPLES[@]}"; do
            n=$(echo "$s" | sed 's/.*-//; s/^PBRM[0-9]*$/PBRM1/')
            [[ "$n" == "$tgt" ]] && { tb+=("${IN_DIR}/${s}.final.bam"); tl+=("$s"); }
        done
        [[ ${#tb[@]} -ge 1 ]] || continue
        plotFingerprint -b "${tb[@]}" --labels "${tl[@]}" \
            --ignoreDuplicates --minMappingQuality 30 \
            --plotFile "$outp" --outQualityMetrics "$QC/fingerprint_${tgt}.txt" \
            -p "$THREADS" > "$QC/logs/fingerprint_${tgt}.log" 2>&1 \
          && echo "[B] $tgt done" || echo "[B] $tgt FAILED"
    done
fi

# ----------------------------------------------------------------- C. bigwigs
if [[ "$PHASES" == *C* ]]; then
    echo "[C] bigwigs (CPM, 10bp bins, dups excluded)"
    export IN_DIR BW QC
    export BW_THREADS=$(( THREADS / BW_JOBS )); [[ $BW_THREADS -lt 1 ]] && BW_THREADS=1
    make_bw() {
        local s="$1"
        local out="${BW}/${s}.cpm.bw"
        [[ -s "$out" ]] && { echo "[C] skip $s"; return 0; }
        bamCoverage -b "${IN_DIR}/${s}.final.bam" -o "${out}.part" \
            --binSize 10 --normalizeUsing CPM --extendReads \
            --ignoreDuplicates --minMappingQuality 30 \
            -p "$BW_THREADS" > "${QC}/logs/bw_${s}.log" 2>&1 \
          && mv "${out}.part" "$out" && echo "[C] $s done" \
          || { rm -f "${out}.part"; echo "[C] $s FAILED"; }
    }
    export -f make_bw
    printf '%s\n' "${SAMPLES[@]}" | xargs -P "$BW_JOBS" -I{} bash -c 'make_bw "$@"' _ {}
fi

# -------------------------------------------------------------------- D. TSS
if [[ "$PHASES" == *D* ]]; then
    TSSBED="$QC/${HOMER_GENOME}_tss.bed"
    HOMER_TSS=$(ls ~/miniconda/envs/"$HOMER_ENV"/share/homer*/data/genomes/"$HOMER_GENOME"/"$HOMER_GENOME".tss 2>/dev/null | head -1)

    if [[ ! -s "$TSSBED" ]]; then
        [[ -s "${HOMER_TSS:-}" ]] || { echo "[D] ERROR: HOMER $HOMER_GENOME.tss not found"; exit 1; }
        echo "[D] building TSS bed from $HOMER_TSS"
        # HOMER <genome>.tss is: <refseq_id> <chr> <start> <end> <strandcode>
        # with every region exactly 4000bp = TSS +/-2000, so the TSS is the
        # CENTER (start+2000), and strandcode is 0=plus, 1=minus.
        awk -F'\t' 'NF>=5{print $4-$3}' "$HOMER_TSS" | sort -n | uniq -c \
            | sort -rn | head -5 > "$QC/logs/tss_widths.txt"
        W=$(head -1 "$QC/logs/tss_widths.txt" | awk '{print $2}')
        echo "[D] dominant region width = ${W}bp (expect 4000)"
        [[ "$W" == "4000" ]] || { echo "[D] ERROR: unexpected width $W"; exit 1; }
        awk -F'\t' 'NF>=5 && $2 ~ /^chr([0-9]+|X|Y)$/ {
              tss = $3 + 2000
              st  = ($5 == 1) ? "-" : "+"
              if (tss < 1) next
              key = $2"\t"tss"\t"st
              if (!(key in seen)) { seen[key]=1
                  printf "%s\t%d\t%d\t%s\t0\t%s\n", $2, tss-1, tss, $1, st }
            }' "$HOMER_TSS" | sort -k1,1 -k2,2n > "${TSSBED}.part"
        n=$(wc -l < "${TSSBED}.part")
        bad=$(awk -F'\t' '$6!="+" && $6!="-"' "${TSSBED}.part" | wc -l)
        if [[ "$n" -lt 10000 || "$bad" -gt 0 ]]; then
            echo "[D] ERROR: TSS bed looks wrong (n=$n, bad_strand=$bad)"; exit 1
        fi
        mv "${TSSBED}.part" "$TSSBED"
        echo "[D] $n unique TSS positions"
    fi

    if [[ ! -s "$QC/tss_matrix.gz" ]]; then
        mapfile -t BWS < <(for s in "${SAMPLES[@]}"; do
            [[ -s "${BW}/${s}.cpm.bw" ]] && echo "${BW}/${s}.cpm.bw"; done)
        mapfile -t BWL < <(for s in "${SAMPLES[@]}"; do
            [[ -s "${BW}/${s}.cpm.bw" ]] && echo "$s"; done)
        echo "[D] computeMatrix over ${#BWS[@]} bigwigs"
        computeMatrix reference-point --referencePoint TSS \
            -S "${BWS[@]}" -R "$TSSBED" --samplesLabel "${BWL[@]}" \
            -b 2000 -a 2000 --binSize 10 --skipZeros --missingDataAsZero \
            -o "$QC/tss_matrix.gz" -p "$THREADS" \
            > "$QC/logs/computeMatrix.log" 2>&1 \
          && echo "[D] matrix done" || { echo "[D] computeMatrix FAILED"; exit 1; }
    fi

    [[ -s "$QC/tss_profile.png" ]] || plotProfile -m "$QC/tss_matrix.gz" \
        -o "$QC/tss_profile.png" --outFileNameData "$QC/tss_profile_data.tab" \
        --perGroup --plotTitle "TSS +/-2kb, CPM" \
        > "$QC/logs/plotProfile.log" 2>&1
    [[ -s "$QC/tss_heatmap.png" ]] || plotHeatmap -m "$QC/tss_matrix.gz" \
        -o "$QC/tss_heatmap.png" --sortUsing mean --sortRegions descend \
        > "$QC/logs/plotHeatmap.log" 2>&1

    # numeric TSS enrichment score = mean signal within +/-250bp of TSS
    # divided by mean signal in the outer 500bp of each flank (background)
    if [[ -s "$QC/tss_profile_data.tab" ]]; then
        awk -F'\t' 'NR<=2{next}
          { nb=NF-2; c=0; b=0; cn=0; bn=0
            for(i=3;i<=NF;i++){ v=$i+0; pos=(i-3)
              if (pos>=175 && pos<225) { c+=v; cn++ }
              else if (pos<50 || pos>=350) { b+=v; bn++ } }
            if(cn&&bn&&b>0) printf "%s\t%.3f\n", $1, (c/cn)/(b/bn)
            else printf "%s\tNA\n", $1 }' \
          "$QC/tss_profile_data.tab" | sort -k2 -rn > "$QC/tss_enrichment_scores.tsv"
        echo "[D] TSS scores -> $QC/tss_enrichment_scores.tsv"
    fi
fi

# ------------------------------------------------------------- E. correlation
if [[ "$PHASES" == *E* ]]; then
    if [[ ! -s "$QC/bins5kb.npz" ]]; then
        echo "[E] multiBamSummary 5kb bins (slowest phase)"
        multiBamSummary bins -b "${BAMS[@]}" --labels "${SAMPLES[@]}" \
            --binSize 5000 --ignoreDuplicates --minMappingQuality 30 \
            -o "$QC/bins5kb.npz" --outRawCounts "$QC/bins5kb_counts.tab" \
            -p "$THREADS" > "$QC/logs/multiBamSummary.log" 2>&1 \
          && echo "[E] counts done" || { echo "[E] multiBamSummary FAILED"; exit 1; }
    fi
    for m in spearman pearson; do
        [[ -s "$QC/corr_${m}.png" ]] || plotCorrelation -in "$QC/bins5kb.npz" \
            -c "$m" -p heatmap --skipZeros --removeOutliers --plotNumbers \
            --colorMap RdYlBu_r -o "$QC/corr_${m}.png" \
            --outFileCorMatrix "$QC/corr_${m}.tab" \
            > "$QC/logs/corr_${m}.log" 2>&1
    done
    [[ -s "$QC/pca.png" ]] || plotPCA -in "$QC/bins5kb.npz" -o "$QC/pca.png" \
        --outFileNameData "$QC/pca.tab" > "$QC/logs/pca.log" 2>&1
    echo "[E] done"
fi

echo "=== step 2 finished for $GROUP -> $QC ==="
