#!/bin/bash
# Builds the three aggregate tables make_report.py reads. None of these are
# produced automatically by the earlier steps (they only print [done] lines to
# stdout) -- this closes that gap so re-running the pipeline never leaves a
# group's report half-populated, or a group's frip_nodup silently reading zero.
#
#   qc/<group>/contig_stats.tsv        raw (pre-filter) chrM / scaffold / primary read counts
#   qc/<group>/align_stats.tsv         bowtie2 pair counts + alignment rate
#   filtered/<group>/step1_summary.tsv pct_kept / dup_frac / usable_est after step 1
#
# Usage:  ./build_stats.sh [group]      # group defaults to cutrun
#
# Idempotent: rebuilds all three files from source artifacts (raw BAM, final
# BAM, Picard metrics, bowtie2 log) every run -- cheap enough to always redo
# rather than track partial state.

set -u
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/../config.sh"
cd "$PROJECT_DIR"

GROUP="${1:-cutrun}"
source "$CONDA_SH"
conda activate "$ALIGN_ENV"
for t in samtools; do
    command -v "$t" >/dev/null || { echo "ERROR: $t not found in $ALIGN_ENV"; exit 1; }
done

RAW_DIR="bams/${GROUP}"
FILT_DIR="filtered/${GROUP}"
QC_DIR="qc/${GROUP}"
mkdir -p "$QC_DIR"

mapfile -t SAMPLES < <(ls "${FILT_DIR}"/*.final.bam 2>/dev/null \
    | xargs -r -n1 basename | sed 's/\.final\.bam$//' | sort)
[[ ${#SAMPLES[@]} -gt 0 ]] || { echo "ERROR: no final BAMs in $FILT_DIR"; exit 1; }
echo "group=$GROUP samples=${#SAMPLES[@]}"

# ---------------------------------------------------------- A. contig_stats.tsv
CONTIG="${QC_DIR}/contig_stats.tsv"
echo -e "sample\traw_total\tchrM\tpct_chrM\tscaffold\tpct_scaffold\tprimary\tpct_primary" > "$CONTIG"
for s in "${SAMPLES[@]}"; do
    bam="${RAW_DIR}/${s}.sorted.bam"
    if [[ ! -s "$bam" ]]; then
        echo "[contig] skip $s (no raw bam at $bam)"
        continue
    fi
    read -r chrm scaf prim <<< "$(samtools idxstats "$bam" | awk '
        $1=="chrM"{m+=$3; next}
        $1 ~ /^chr([0-9]+|X|Y)$/ {p+=$3; next}
        $1!="*"{s+=$3}
        END{printf "%d %d %d", m+0, s+0, p+0}')"
    tot=$((chrm+scaf+prim))
    awk -v s="$s" -v tot="$tot" -v m="$chrm" -v sc="$scaf" -v p="$prim" 'BEGIN{
        printf "%s\t%d\t%d\t%.3f\t%d\t%.3f\t%d\t%.2f\n", s, tot, m, 100*m/tot, sc, 100*sc/tot, p, 100*p/tot}' >> "$CONTIG"
    echo "[contig] $s"
done

# ----------------------------------------------------------- B. align_stats.tsv
ALIGN="${QC_DIR}/align_stats.tsv"
echo -e "sample\ttotal_pairs\talign_rate\tconc_1\tconc_multi" > "$ALIGN"
for s in "${SAMPLES[@]}"; do
    log="logs/${s}_bowtie2.log"
    if [[ ! -s "$log" ]]; then
        echo "[align] skip $s (no bowtie2 log at $log)"
        continue
    fi
    awk -v s="$s" '
        /^[0-9]+ reads/{pairs=$1}
        /aligned concordantly exactly 1 time/{c1=$1}
        /aligned concordantly >1 times/{cm=$1}
        /overall alignment rate/{rate=$1; gsub("%","",rate)}
        END{printf "%s\t%d\t%s\t%d\t%d\n", s, pairs, rate, c1, cm}' "$log" >> "$ALIGN"
    echo "[align] $s"
done

# ------------------------------------------------------- C. step1_summary.tsv
STEP1="${FILT_DIR}/step1_summary.tsv"
echo -e "sample\traw\tfinal\tpct_kept\tdup_frac\tusable_est" > "$STEP1"
for s in "${SAMPLES[@]}"; do
    raw=$(awk -F'\t' -v s="$s" '$1==s{print $2}' "$CONTIG")
    if [[ -z "${raw:-}" || "$raw" -le 0 ]]; then
        echo "[step1] skip $s (no raw count in $CONTIG)"
        continue
    fi
    final_bam="${FILT_DIR}/${s}.final.bam"
    final=$(samtools idxstats "$final_bam" | awk '{t+=$3}END{print t+0}')
    dup=$(awk -F'\t' '/^LIBRARY\t/{getline; print $9; exit}' \
          "${FILT_DIR}/metrics/${s}_dup_metrics.txt" 2>/dev/null)
    usable=$(samtools view -c -F 1024 "$final_bam")
    awk -v s="$s" -v raw="$raw" -v final="$final" -v dup="${dup:-0}" -v usable="$usable" 'BEGIN{
        printf "%s\t%d\t%d\t%.1f\t%s\t%d\n", s, raw, final, 100*final/raw, dup, usable}' >> "$STEP1"
    echo "[step1] $s"
done

echo "=== stats built for $GROUP -> $CONTIG, $ALIGN, $STEP1 ==="
