# CUT&RUN / ATAC pipeline

Reproducible trim -> align -> filter/dedup -> peak-call -> QC -> report pipeline
for paired-end CUT&RUN and Tn5-tagmented ("ATAC") libraries aligned to a single
reference genome. Both assay groups run through the *same* scripts, parameterized
by `[group]` (`cutrun` or `atac`); the only assay-specific branch is the trim
adapter (`--illumina` vs `--nextera`).

## 1. Set up the environments

Four conda environments cover the whole pipeline:

| env file | used for |
|---|---|
| `envs/cutrun_env.yml` | trim_galore, fastqc, bowtie2, samtools, picard, bedtools |
| `envs/macs2-env2.yml` | MACS2 peak calling |
| `envs/bigwig_env.yml` | deeptools (bamCoverage, plotFingerprint, computeMatrix, ...) |
| `envs/homer_env.yml` | HOMER, only used to source a genome's TSS annotation file |

```bash
for f in envs/*.yml; do conda env create -f "$f"; done
```

> A plain `macs2-env` (without the `2`) is known to break with
> `ImportError: undefined symbol: __log_finite` on some systems — always use the
> `macs2-env2` spec above, and `run_macs2_frip.sh` checks `macs2 --version` at
> startup to catch this before calling peaks.

## 2. Point config.sh at your reference files

Edit [config.sh](config.sh) (or export the same variables before running any
script) to match your genome:

```bash
export BOWTIE2_INDEX=/path/to/genome/genome   # bowtie2 index prefix
export BLACKLIST=/path/to/genome-blacklist.bed
export MACS2_GENOME=mm                        # macs2 -g: mm=1.87e9, hs=2.7e9, ...
export GENOME_SIZE_BP=1.87e9                   # must match MACS2_GENOME
export HOMER_GENOME=mm10                       # looked up in <homer_env>/share/homer*/data/genomes/<name>
```

## 3. Lay out the input data

From an empty `PROJECT_DIR` (defaults to wherever you run the scripts from):

```
<PROJECT_DIR>/
├── <sample>_R1_001.fastq.gz
├── <sample>_R2_001.fastq.gz
└── ...
```

Sample names ending in `-ATAC` are routed to the `atac` group; everything else
goes to `cutrun` (edit `group_of()` in `scripts/organize_by_assay.sh` if your
naming convention differs).

## 4. Run it

```bash
cd /path/to/your/project        # becomes PROJECT_DIR
/path/to/cutrun-atac-pipeline/scripts/run_upstream.sh        # once: organize, trim, align
/path/to/cutrun-atac-pipeline/scripts/run_downstream.sh cutrun
/path/to/cutrun-atac-pipeline/scripts/run_downstream.sh atac
```

Or drive each step yourself for more control (every script is independently
resumable/idempotent — skips samples that already have output):

```bash
scripts/organize_by_assay.sh
scripts/run_trim_galore.sh
scripts/run_bowtie2.sh
scripts/run_filter_dedup.sh   cutrun
scripts/build_stats.sh        cutrun
scripts/run_macs2_frip.sh     cutrun
scripts/run_qc_step2.sh       cutrun
python3 scripts/make_report.py cutrun
```

Each per-group step also accepts `JOBS=`, `PHASES=`, etc. — see the comment
block at the top of the script.

## Output layout

```
<PROJECT_DIR>/
├── fastq/<group>/          raw fastqs, sorted by assay
├── trimmed/<group>/        trim_galore output + fastqc
├── bams/<group>/           bowtie2-aligned, sorted BAMs
├── filtered/<group>/       MAPQ/proper-pair/blacklist-filtered, dedup-marked BAMs
│   ├── metrics/            Picard MarkDuplicates metrics
│   ├── bigwigs/             CPM bigwigs (dups excluded)
│   └── step1_summary.tsv    per-sample raw/final/pct_kept/dup_frac/usable_est
├── qc/<group>/             fragment size, fingerprint, TSS, correlation/PCA, report.html
│   ├── align_stats.tsv       bowtie2 pair counts + alignment rate
│   └── contig_stats.tsv      raw chrM / scaffold / primary-chrom read counts
├── peaks/<group>/<sample>/  MACS2 narrowPeak + bedGraph
│   └── frip_summary.tsv      per-sample FRiP, with and without duplicates
└── logs/                   every step's stdout/stderr
```

## Why `build_stats.sh` exists as its own step

`align_stats.tsv`, `contig_stats.tsv` and `step1_summary.tsv` are *not*
side effects of `run_bowtie2.sh` / `run_filter_dedup.sh` (those only print
`[done]` lines to stdout). `run_macs2_frip.sh`'s duplicate-excluded FRiP column
reads `step1_summary.tsv` for its denominator — if that file is missing, it
silently defaults to 0 rather than erroring, and `make_report.py` will show
blank `Total reads` / `Align %` / `Kept %` / `Dup frac` / `Usable reads` cells.
Always run `build_stats.sh <group>` after `run_filter_dedup.sh <group>` and
before `run_macs2_frip.sh <group>` / `make_report.py <group>`.

## Key thresholds (same for both groups)

- Trim: `-q 20 --length 20`, adapter matched to prep (`--illumina`/`--nextera`)
- Align: `--very-sensitive --no-mixed --no-discordant -X 2000`
- Filter: `-f 2 -F 1804 -q 30` (the ENCODE ATAC-seq BAM-filtering standard)
- Dedup: Picard MarkDuplicates, marked not removed; excluded downstream via `-F 1024`
- Peaks: `-f BAMPE -q 0.01 --nolambda --keep-dup all --call-summits` (no shift
  model — appropriate for real paired-end fragments; drop `--nolambda` if you
  have an IgG/input control)
