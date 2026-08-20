#!/bin/bash
# Step 3: MACS2 peak calling + FRiP for the filtered BAMs.
#
#   -f BAMPE        real fragment coordinates, no shift modelling
#   --nolambda      genome-wide background instead of local lambda (no IgG /
#                   control library case; drop this flag if you have one)
#   --keep-dup all  duplicates retained during calling; set KEEPDUP=1 to collapse
#   --call-summits  deconvolve merged/subpeak structure
#   -B --SPMR       bedGraph pileup, signal-per-million-reads normalised
#   -q 0.01
#   -g $MACS2_GENOME  effective genome size (see config.sh)
#
# Run in $MACS2_ENV. A plain "macs2-env" build has been seen raising
# "ImportError: undefined symbol: __log_finite" -- if that happens, rebuild
# the env from envs/macs2-env2.yml.
#
# FRiP is reported both including and excluding flagged duplicates -- this
# reads filtered/<group>/step1_summary.tsv (built by build_stats.sh) for the
# non-duplicate denominator; run build_stats.sh first or frip_nodup silently
# comes out as zero.
#
# Usage:  ./run_macs2_frip.sh [group]      # default cutrun
#         JOBS=4 KEEPDUP=1 ./run_macs2_frip.sh
#         PHASES=F ./run_macs2_frip.sh     # recompute FRiP only
#         NOBDG=1 ./run_macs2_frip.sh      # drop -B --SPMR (saves disk)
#
# Resumable: a sample with a non-empty _peaks.narrowPeak is skipped.

set -u
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/../config.sh"
cd "$PROJECT_DIR"

GROUP="${1:-cutrun}"
JOBS="${JOBS:-4}"
QVAL="${QVAL:-0.01}"
GSIZE="${GSIZE:-$MACS2_GENOME}"
KEEPDUP="${KEEPDUP:-all}"
NOBDG="${NOBDG:-0}"
PHASES="${PHASES:-PF}"

source "$CONDA_SH"

export IN_DIR="filtered/${GROUP}"
export OUT_DIR="peaks/${GROUP}"
export LOG_DIR="logs/macs2_${GROUP}"
export QVAL GSIZE KEEPDUP NOBDG
mkdir -p "$OUT_DIR" "$LOG_DIR"

mapfile -t SAMPLES < <(ls "${IN_DIR}"/*.final.bam 2>/dev/null \
    | xargs -r -n1 basename | sed 's/\.final\.bam$//' | sort)
[[ ${#SAMPLES[@]} -gt 0 ]] || { echo "ERROR: no final BAMs in $IN_DIR"; exit 1; }
echo "group=$GROUP samples=${#SAMPLES[@]} jobs=$JOBS q=$QVAL g=$GSIZE keep-dup=$KEEPDUP phases=$PHASES"

# ------------------------------------------------------------------ P. peaks
if [[ "$PHASES" == *P* ]]; then
    conda activate "$MACS2_ENV"
    command -v macs2 >/dev/null || { echo "ERROR: macs2 not found in $MACS2_ENV"; exit 1; }
    macs2 --version >/dev/null 2>&1 || { echo "ERROR: macs2 import is broken in $MACS2_ENV"; exit 1; }

    call_one() {
        local s="$1"
        local np="${OUT_DIR}/${s}/${s}_peaks.narrowPeak"
        if [[ -s "$np" ]]; then
            echo "[P] skip $s ($(wc -l < "$np") peaks)"
            return 0
        fi
        mkdir -p "${OUT_DIR}/${s}"
        local bdg=(-B --SPMR)
        [[ "$NOBDG" == "1" ]] && bdg=()
        if macs2 callpeak \
                -t "${IN_DIR}/${s}.final.bam" \
                -f BAMPE -g "$GSIZE" -q "$QVAL" \
                --nolambda --keep-dup "$KEEPDUP" --call-summits \
                "${bdg[@]}" \
                -n "$s" --outdir "${OUT_DIR}/${s}" \
                > "${LOG_DIR}/${s}.log" 2>&1 && [[ -s "$np" ]]
        then
            echo "[P] $s  peaks=$(wc -l < "$np")"
        else
            echo "[P] FAIL $s -- see ${LOG_DIR}/${s}.log"
        fi
    }
    export -f call_one

    printf '%s\n' "${SAMPLES[@]}" | xargs -P "$JOBS" -I{} bash -c 'call_one "$@"' _ {}
    conda deactivate
fi

# ------------------------------------------------------------------- F. FRiP
if [[ "$PHASES" == *F* ]]; then
    conda activate "$ALIGN_ENV"
    for t in samtools bedtools; do
        command -v $t >/dev/null || { echo "ERROR: $t not found"; exit 1; }
    done

    S1="${IN_DIR}/step1_summary.tsv"
    [[ -s "$S1" ]] || echo "WARNING: $S1 missing -- run build_stats.sh $GROUP first, or frip_nodup will read 0"
    SUM="${OUT_DIR}/frip_summary.tsv"
    printf "sample\tn_peaks\tpeak_bp\tpct_genome\ttotal_all\tin_peaks_all\tfrip_all\tusable_nodup\tin_peaks_nodup\tfrip_nodup\n" > "$SUM"

    for s in "${SAMPLES[@]}"; do
        np="${OUT_DIR}/${s}/${s}_peaks.narrowPeak"
        bam="${IN_DIR}/${s}.final.bam"
        [[ -s "$np" ]] || { echo "[F] skip $s (no peaks)"; continue; }

        # merge first: --call-summits emits overlapping entries per summit, which
        # would double-count both bases and reads
        merged="${OUT_DIR}/${s}/${s}_peaks.merged.bed"
        sort -k1,1 -k2,2n "$np" | bedtools merge -i - > "$merged"

        n=$(wc -l < "$merged")
        bp=$(awk '{t+=$3-$2}END{print t+0}' "$merged")
        tot=$(samtools idxstats "$bam" 2>/dev/null | awk '{t+=$3}END{print t+0}')
        nod=$(awk -F'\t' -v s="$s" '$1==s{print $6}' "$S1" 2>/dev/null)
        [[ -n "${nod:-}" ]] || nod=0
        ina=$(samtools view -c -L "$merged" "$bam" 2>/dev/null)
        inn=$(samtools view -c -F 1024 -L "$merged" "$bam" 2>/dev/null)

        awk -v s="$s" -v n="$n" -v bp="$bp" -v tot="$tot" -v ina="$ina" \
            -v nod="$nod" -v inn="$inn" -v gsize="$GENOME_SIZE_BP" 'BEGIN{
            printf "%s\t%d\t%d\t%.3f\t%d\t%d\t%.4f\t%d\t%d\t%.4f\n",
              s, n, bp, 100*bp/gsize, tot, ina, (tot>0?ina/tot:0),
              nod, inn, (nod>0?inn/nod:0)}' >> "$SUM"
        echo "[F] $s peaks=$n frip_all=$(awk -v a="$ina" -v b="$tot" 'BEGIN{printf "%.3f",(b>0?a/b:0)}')"
    done

    echo; echo "=== FRiP summary -> $SUM ==="
    column -t "$SUM"
fi

echo "=== step 3 finished for $GROUP -> $OUT_DIR ==="
