#!/bin/bash
# Bowtie2 alignment for all trimmed samples in PROJECT_DIR (run once, before
# splitting by group).
#
# JOBS x THREADS should stay near your core count, e.g. 2x8 on 16 vCPUs.

set -u
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/../config.sh"
cd "$PROJECT_DIR"

source "$CONDA_SH"
conda activate "$ALIGN_ENV"
command -v bowtie2  >/dev/null || { echo "ERROR: bowtie2 not found";  exit 1; }
command -v samtools >/dev/null || { echo "ERROR: samtools not found"; exit 1; }

export THREADS="${THREADS:-8}"          # bowtie2 threads per sample
export SORT_THREADS="${SORT_THREADS:-2}"
export SORT_MEM="${SORT_MEM:-2G}"
export BAM_DIR=bams
JOBS="${JOBS:-2}"                       # samples in flight

mkdir -p "$BAM_DIR" logs

align_one() {
    sample="$1"
    r1="trimmed/${sample}_R1_001_val_1.fq.gz"
    r2="trimmed/${sample}_R2_001_val_2.fq.gz"
    out="${BAM_DIR}/${sample}.sorted.bam"

    if [[ -s "$out" ]]; then
        echo "[skip] $sample already aligned"
        return 0
    fi
    if [[ ! -s "$r2" ]]; then
        echo "[FAIL] $sample missing R2 ($r2)"
        return 0
    fi

    echo "[start] $sample"
    # write to .part so an interrupted run never leaves a truncated BAM that
    # the -s check above would mistake for a finished sample
    if bowtie2 -x "$BOWTIE2_INDEX" \
            -1 "$r1" -2 "$r2" \
            --very-sensitive --no-mixed --no-discordant -X 2000 \
            -p "$THREADS" \
            2> "logs/${sample}_bowtie2.log" \
         | samtools sort -@ "$SORT_THREADS" -m "$SORT_MEM" -o "${out}.part" - \
       && mv "${out}.part" "$out" \
       && samtools index "$out"
    then
        echo "[done] $sample  ($(grep 'overall alignment rate' "logs/${sample}_bowtie2.log"))"
    else
        rm -f "${out}.part"
        echo "[FAIL] $sample -- see logs/${sample}_bowtie2.log"
    fi
}
export -f align_one
export BOWTIE2_INDEX

ls trimmed/*_R1_001_val_1.fq.gz \
  | sed 's|^trimmed/||; s|_R1_001_val_1\.fq\.gz$||' \
  | xargs -P "$JOBS" -I{} bash -c 'align_one "$@"' _ {}

echo "All samples finished."
