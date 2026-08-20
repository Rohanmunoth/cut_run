#!/bin/bash
# Sort samples into assay groups by filename convention:
#
#   atac    *-ATAC          -- Tn5/Nextera-prepped libraries, processed as CUT&RUN
#   cutrun  everything else -- TruSeq ligation CUT&RUN
#
# Moves fastq (root -> fastq/<group>/), trimmed/ and bams/ into per-group
# subdirs. All moves are same-filesystem renames, so no data is copied.
# Idempotent: re-running finds nothing left to move.
#
# DRY_RUN=1 ./organize_by_assay.sh   # print the plan, move nothing

set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/../config.sh"
cd "$PROJECT_DIR"
DRY_RUN="${DRY_RUN:-0}"

# assay group for a sample name -- edit this if your naming convention differs
group_of() {
    case "$1" in
        *-ATAC) echo atac   ;;
        *)      echo cutrun ;;
    esac
}

# derive the canonical sample list from the raw R1 fastqs
mapfile -t SAMPLES < <(ls *_R1_001.fastq.gz 2>/dev/null | sed 's/_R1_001\.fastq\.gz$//' | sort)
if [[ ${#SAMPLES[@]} -eq 0 ]] && [[ -d fastq ]]; then
    mapfile -t SAMPLES < <(ls fastq/*/*_R1_001.fastq.gz 2>/dev/null \
        | xargs -r -n1 basename | sed 's/_R1_001\.fastq\.gz$//' | sort)
fi
echo "samples found: ${#SAMPLES[@]}"

moved=0
move_to() {   # move_to <dest_dir> <file>
    local dest="$1" f="$2"
    [[ -e "$f" ]] || { return 0; }
    if [[ "$DRY_RUN" == "1" ]]; then
        echo "  mv $f -> $dest/"
    else
        mkdir -p "$dest"
        mv -n "$f" "$dest/"
    fi
    moved=$((moved+1))
}

for s in "${SAMPLES[@]}"; do
    g=$(group_of "$s")

    # --- raw fastq + checksums: project root -> fastq/<group>/
    for r in R1 R2; do
        move_to "fastq/$g"  "${s}_${r}_001.fastq.gz"
        move_to "fastq/$g"  "${s}_${r}_001.fastq.gz.md5"
    done

    # --- trim_galore output: trimmed/ -> trimmed/<group>/
    move_to "trimmed/$g" "trimmed/${s}_R1_001_val_1.fq.gz"
    move_to "trimmed/$g" "trimmed/${s}_R2_001_val_2.fq.gz"
    for r in R1 R2; do
        move_to "trimmed/$g" "trimmed/${s}_${r}_001.fastq.gz_trimming_report.txt"
    done
    move_to "trimmed/$g" "trimmed/${s}_R1_001_val_1_fastqc.html"
    move_to "trimmed/$g" "trimmed/${s}_R1_001_val_1_fastqc.zip"
    move_to "trimmed/$g" "trimmed/${s}_R2_001_val_2_fastqc.html"
    move_to "trimmed/$g" "trimmed/${s}_R2_001_val_2_fastqc.zip"

    # --- bams: already split into per-group dirs, so also look there
    for src in "bams/${s}.sorted.bam" "bams/atac/${s}.sorted.bam" "bams/cutrun/${s}.sorted.bam"; do
        [[ -e "$src" ]] || continue
        [[ "$(dirname "$src")" == "bams/$g" ]] && continue   # already correct
        move_to "bams/$g" "$src"
        move_to "bams/$g" "${src}.bai"
    done
done

echo "files moved: $moved  (DRY_RUN=$DRY_RUN)"
