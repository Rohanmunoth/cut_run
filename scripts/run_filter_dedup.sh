#!/bin/bash
# Step 1 post-alignment processing: filter -> mark duplicates -> blacklist.
#
#   a. samtools view -f 2 -F 1804 -q 30, primary chromosomes only
#      (drops chrM + random/Un scaffolds in the same pass; this is the same
#      -f 2 -F 1804 -q 30 threshold ENCODE uses for ATAC-seq BAM filtering)
#   b. picard MarkDuplicates -- MARKS ONLY, does not remove. Downstream steps
#      exclude flagged duplicates with -F 1024 / --ignoreDuplicates.
#   c. bedtools intersect -v against the blacklist
#
# Usage:  ./run_filter_dedup.sh [group]      # group defaults to cutrun
#         JOBS=4 ./run_filter_dedup.sh atac
#         DRY_RUN=1 ./run_filter_dedup.sh
#
# Resumable: a sample whose final BAM already passes samtools quickcheck is
# skipped. Intermediates are written as .part and removed on failure, so an
# interrupted run never leaves a truncated BAM that looks finished.

set -u
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/../config.sh"
cd "$PROJECT_DIR"

GROUP="${1:-cutrun}"
JOBS="${JOBS:-4}"                  # samples in flight
DRY_RUN="${DRY_RUN:-0}"

source "$CONDA_SH"
conda activate "$ALIGN_ENV"
for t in samtools picard bedtools; do
    command -v "$t" >/dev/null || { echo "ERROR: $t not found in $ALIGN_ENV"; exit 1; }
done

export BAM_DIR="bams/${GROUP}"
export OUT_DIR="filtered/${GROUP}"
export MET_DIR="${OUT_DIR}/metrics"
export LOG_DIR="logs/filter_${GROUP}"
export TMP_DIR="${OUT_DIR}/tmp"
export BLACKLIST
export MAPQ=30
export PICARD_MEM="${PICARD_MEM:-4g}"   # JOBS x PICARD_MEM must stay under available RAM;
                                        # ATAC libraries are typically much deeper, budget accordingly
export SORT_THREADS=2

[[ -d "$BAM_DIR" ]]      || { echo "ERROR: $BAM_DIR not found"; exit 1; }
[[ -s "$BLACKLIST" ]]    || { echo "ERROR: blacklist $BLACKLIST missing"; exit 1; }
mkdir -p "$OUT_DIR" "$MET_DIR" "$LOG_DIR" "$TMP_DIR"

# primary chromosomes; everything else (chrM, *_random, chrUn_*) is dropped
export CHROMS=$(echo chr{1..19} chrX chrY)

process_one() {
    local s="$1"
    local in="${BAM_DIR}/${s}.sorted.bam"
    local final="${OUT_DIR}/${s}.final.bam"
    local log="${LOG_DIR}/${s}.log"
    local filt="${TMP_DIR}/${s}.filt.bam"
    local md="${TMP_DIR}/${s}.md.bam"

    if [[ -s "$final" ]] && samtools quickcheck "$final" 2>/dev/null; then
        echo "[skip] $s already done"
        return 0
    fi

    echo "[start] $s"
    {
        echo "=== $s ==="
        set -x

        samtools view -b -f 2 -F 1804 -q "$MAPQ" -@ "$SORT_THREADS" \
            -o "${filt}.part" "$in" $CHROMS && mv "${filt}.part" "$filt" || exit 1

        picard -Xmx"$PICARD_MEM" MarkDuplicates \
            I="$filt" O="${md}.part" M="${MET_DIR}/${s}_dup_metrics.txt" \
            ASSUME_SORTED=true REMOVE_DUPLICATES=false \
            TMP_DIR="$TMP_DIR" VALIDATION_STRINGENCY=LENIENT || exit 1
        mv "${md}.part" "$md"
        rm -f "$filt"

        bedtools intersect -v -abam "$md" -b "$BLACKLIST" > "${final}.part" || exit 1
        mv "${final}.part" "$final"
        rm -f "$md"

        samtools index "$final" || exit 1
    } >>"$log" 2>&1

    if [[ -s "$final" ]] && samtools quickcheck "$final" 2>/dev/null; then
        local kept dup
        kept=$(samtools view -c -F 1024 "$final" 2>/dev/null)
        # metrics are tab-delimited and LIBRARY is "Unknown Library" (contains a
        # space), so FS must be tab or every field shifts by one
        dup=$(awk -F'\t' '/^LIBRARY\t/{getline; print $9; exit}' \
              "${MET_DIR}/${s}_dup_metrics.txt" 2>/dev/null)
        echo "[done] $s  usable_reads=${kept}  dup_frac=${dup:-NA}"
    else
        rm -f "${final}.part" "${md}.part" "${filt}.part"
        echo "[FAIL] $s -- see $log"
    fi
}
export -f process_one

mapfile -t SAMPLES < <(ls "${BAM_DIR}"/*.sorted.bam 2>/dev/null \
    | xargs -r -n1 basename | sed 's/\.sorted\.bam$//' | sort)
echo "group=$GROUP  samples=${#SAMPLES[@]}  jobs=$JOBS  out=$OUT_DIR"

if [[ "$DRY_RUN" == "1" ]]; then
    printf '  would process: %s\n' "${SAMPLES[@]}"
    exit 0
fi

printf '%s\n' "${SAMPLES[@]}" | xargs -P "$JOBS" -I{} bash -c 'process_one "$@"' _ {}

rmdir "$TMP_DIR" 2>/dev/null
echo "=== finished: $(ls "$OUT_DIR"/*.final.bam 2>/dev/null | wc -l)/${#SAMPLES[@]} final BAMs ==="
