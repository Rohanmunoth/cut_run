#!/bin/bash
# Trim Galore for all samples in PROJECT_DIR (run once, before splitting by group).
# ATAC (Tn5/Nextera) libraries are trimmed with --nextera; everything else with
# --illumina (TruSeq CUT&RUN adapters).

set -u
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/../config.sh"
cd "$PROJECT_DIR"

source "$CONDA_SH"
conda activate "$ALIGN_ENV"
command -v trim_galore >/dev/null || { echo "ERROR: trim_galore not found in $ALIGN_ENV"; exit 1; }

OUTDIR=trimmed
JOBS="${JOBS:-4}"
mkdir -p "$OUTDIR" logs

run_one() {
  base="$1"
  if [[ "$base" == *ATAC* ]]; then
    adapter=--nextera
  else
    adapter=--illumina
  fi

  if [[ -s "$OUTDIR/${base}_R2_001_val_2.fq.gz" ]]; then
    echo "[skip] $base already trimmed"
    return 0
  fi

  echo "[start] $base ($adapter)"
  trim_galore --paired --fastqc \
    --quality 20 --length 20 "$adapter" --trim-n \
    --output_dir "$OUTDIR" \
    "${base}_R1_001.fastq.gz" "${base}_R2_001.fastq.gz" \
    > "logs/${base}_trim_galore.log" 2>&1 \
    && echo "[done] $base" || echo "[FAIL] $base -- see logs/${base}_trim_galore.log"
}
export -f run_one
export OUTDIR

ls *_R1_001.fastq.gz | sed 's/_R1_001\.fastq\.gz$//' \
  | xargs -P "$JOBS" -I{} bash -c 'run_one "$@"' _ {}

echo "All samples finished."
