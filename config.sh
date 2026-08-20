#!/bin/bash
# Central config for the CUT&RUN / ATAC pipeline.
# Every script sources this first. Override any value by exporting it in your
# shell before running a script, e.g.:
#   BOWTIE2_INDEX=/data/genomes/mm10/mm10 ./scripts/run_bowtie2.sh

# Directory holding fastq/, trimmed/, bams/, filtered/, qc/, peaks/, logs/ for
# this run. Defaults to wherever you invoke the scripts from.
export PROJECT_DIR="${PROJECT_DIR:-$(pwd)}"

# Reference genome
export BOWTIE2_INDEX="${BOWTIE2_INDEX:-/home/ubuntu/mm10/mm10}"
export BLACKLIST="${BLACKLIST:-/home/ubuntu/mm10-blacklist.bed}"
export MACS2_GENOME="${MACS2_GENOME:-mm}"          # macs2 -g  (mm=1.87e9, hs=2.7e9, ...)
export GENOME_SIZE_BP="${GENOME_SIZE_BP:-1.87e9}"  # must match MACS2_GENOME; used for %-genome-in-peaks
export HOMER_GENOME="${HOMER_GENOME:-mm10}"        # looked up under <homer_env>/share/homer*/data/genomes/<HOMER_GENOME>

# Conda environments (see envs/*.yml -- create with `conda env create -f envs/<name>.yml`)
export CONDA_SH="${CONDA_SH:-/home/ubuntu/miniconda/etc/profile.d/conda.sh}"
export ALIGN_ENV="${ALIGN_ENV:-cutrun_env}"      # trim_galore, fastqc, bowtie2, samtools, picard, bedtools
export MACS2_ENV="${MACS2_ENV:-macs2-env2}"      # NOTE: a plain "macs2-env" build commonly breaks with
                                                  # "ImportError: undefined symbol: __log_finite" -- use
                                                  # the macs2-env2 spec in envs/, not an old macs2-env.
export BIGWIG_ENV="${BIGWIG_ENV:-bigwig_env}"    # deeptools suite
export HOMER_ENV="${HOMER_ENV:-homer_env}"       # only used to source the TSS annotation file
