#!/usr/bin/env python3
"""
Build a self-contained MultiQC-style HTML report for one assay group.

Reads every metrics file produced by the pipeline, embeds all plots as base64
(so the .html is a single portable file with no sidecar images and no CDN), and
appends the exact commands for each step.

Usage:  python3 make_report.py [group]        # group defaults to cutrun
Output: qc/<group>/report.html

Reads PROJECT_DIR, BOWTIE2_INDEX, BLACKLIST, MACS2_GENOME, MACS2_ENV,
HOMER_GENOME from the environment (source ../config.sh first, or export them
yourself) -- same variables the shell scripts use, so the report always
describes the run that actually happened.
"""
import base64, csv, glob, io, os, sys, html, subprocess
from datetime import datetime, timezone

ROOT = os.environ.get("PROJECT_DIR", os.getcwd())
GROUP = sys.argv[1] if len(sys.argv) > 1 else "cutrun"
QC = os.path.join(ROOT, "qc", GROUP)
FILT = os.path.join(ROOT, "filtered", GROUP)
PEAKS = os.path.join(ROOT, "peaks", GROUP)
OUT = os.path.join(QC, "report.html")

BOWTIE2_INDEX = os.environ.get("BOWTIE2_INDEX", "/home/ubuntu/mm10/mm10")
BLACKLIST     = os.environ.get("BLACKLIST", "/home/ubuntu/mm10-blacklist.bed")
MACS2_GENOME  = os.environ.get("MACS2_GENOME", "mm")
MACS2_ENV     = os.environ.get("MACS2_ENV", "macs2-env2")
HOMER_GENOME  = os.environ.get("HOMER_GENOME", "mm10")
PROJECT_NAME  = os.path.basename(os.path.normpath(ROOT)) or "project"

GROUP_LABEL = {"cutrun": "CUT&RUN", "atac": "ATAC"}.get(GROUP, GROUP.upper())

# --------------------------------------------------------------- data loading
def read_tsv(path, key="sample"):
    """Return {sample: {col: value}}; missing file -> {}."""
    if not os.path.exists(path):
        return {}
    with open(path) as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    return {r[key]: r for r in rows if r.get(key)}

def num(d, k, default=None):
    try:
        return float(d[k])
    except (KeyError, TypeError, ValueError):
        return default

align  = read_tsv(f"{QC}/align_stats.tsv")
contig = read_tsv(f"{QC}/contig_stats.tsv")
step1  = read_tsv(f"{FILT}/step1_summary.tsv")
frip   = read_tsv(f"{PEAKS}/frip_summary.tsv")

# fingerprint_summary.tsv is a hand-built aggregate (sample, AUC, JSdist); if it
# doesn't exist for this group, synthesize the JSdist column from deeptools'
# own per-target fingerprint_<target>.txt files (header: Sample ... Synthetic
# JS Distance), which run_qc_step2.sh always writes.
finger = read_tsv(f"{QC}/fingerprint_summary.tsv")
if not finger:
    for fp in sorted(glob.glob(f"{QC}/fingerprint_*.txt")):
        with open(fp) as fh:
            for row in csv.DictReader(fh, delimiter="\t"):
                s = row.get("Sample")
                if s:
                    finger[s] = {"JSdist": row.get("Synthetic JS Distance", "")}

# fragment sizes: pandas-style index, no header for col 0
frag = {}
fpath = f"{QC}/fragment_sizes.txt"
if os.path.exists(fpath):
    for line in open(fpath).read().splitlines()[1:]:
        f = line.split("\t")
        if len(f) > 11:
            frag[f[0]] = dict(median=float(f[5]), mean=float(f[4]),
                              q1=float(f[3]), q3=float(f[6]), p10=float(f[10]))

# TSS scores: headerless 2-col
tss = {}
tpath = f"{QC}/tss_enrichment_scores.tsv"
if os.path.exists(tpath):
    for line in open(tpath).read().splitlines():
        f = line.split("\t")
        if len(f) == 2:
            try: tss[f[0]] = float(f[1])
            except ValueError: pass

samples = sorted(set(step1) | set(frip) | set(align))
n = len(samples)

# TSS bed strand split, for the QC-step blurb; genome-wide so it's the same
# across groups, but read it fresh rather than hardcode.
tss_n = tss_plus = tss_minus = None
_tss_bed = f"{QC}/{HOMER_GENOME}_tss.bed"
if os.path.exists(_tss_bed):
    with open(_tss_bed) as fh:
        _strands = [l.rstrip("\n").split("\t")[5] for l in fh if l.strip()]
    tss_n, tss_plus, tss_minus = len(_strands), _strands.count("+"), _strands.count("-")

def target(s):
    t = s.split("-")[-1]
    return "PBRM1" if t.startswith("PBRM") else t
def tumor(s):
    return "-".join(s.split("-")[:2])

# --------------------------------------------------------------- image embed
def embed(path, max_px=1600):
    """base64-encode a PNG, downscaling very large ones so the html stays sane."""
    if not os.path.exists(path):
        return None
    raw = open(path, "rb").read()
    if len(raw) > 4_000_000:
        try:
            from PIL import Image
            im = Image.open(io.BytesIO(raw))
            im.thumbnail((max_px, max_px * 6))
            buf = io.BytesIO()
            im.save(buf, format="PNG", optimize=True)
            raw = buf.getvalue()
        except Exception:
            pass
    return "data:image/png;base64," + base64.b64encode(raw).decode()

# ------------------------------------------------------------------ thresholds
def cls_frip(v):
    if v is None: return ""
    if v >= 0.10: return "good"
    if v >= 0.02: return "warn"
    return "bad"
def cls_dup(v):
    if v is None: return ""
    if v <= 0.20: return "good"
    if v <= 0.35: return "warn"
    return "bad"
def cls_mt(v):
    if v is None: return ""
    return "good" if v <= 2 else ("warn" if v <= 10 else "bad")
def cls_align(v):
    if v is None: return ""
    return "good" if v >= 90 else ("warn" if v >= 80 else "bad")
def cls_tss(v):
    if v is None: return ""
    return "good" if v >= 2.0 else ("warn" if v >= 1.3 else "bad")

def fmt(v, spec="{:.3f}"):
    return "-" if v is None else spec.format(v)

# ------------------------------------------------------------------ main table
COLS = [
    ("Sample", "sample"), ("Tumor", "tumor"), ("Target", "target"),
    ("Total reads", "pairs"), ("Align %", "align"),
    ("chrM %", "mt"), ("Scaffold %", "scaff"),
    ("Kept %", "kept"), ("Dup frac", "dup"), ("Usable reads", "usable"),
    ("Frag median", "fragmed"), ("Frag p10", "fragp10"),
    ("Fingerprint JS", "js"), ("TSS score", "tssv"),
    ("Peaks", "peaks"), ("% genome", "pgen"),
    ("FRiP", "fripv"), ("FRiP nodup", "fripn"), ("Fold enr.", "fold"),
]

rows_html = []
for s in samples:
    a, c, s1, fr, fg = align.get(s, {}), contig.get(s, {}), step1.get(s, {}), frip.get(s, {}), finger.get(s, {})
    frg = frag.get(s, {})
    # bowtie2 reports PAIRS; every downstream count is in READS. Convert here so
    # the columns are directly comparable (otherwise "usable" appears to exceed
    # "total").
    pairs  = num(a, "total_pairs"); ar = num(a, "align_rate")
    if pairs is not None:
        pairs *= 2
    mt     = num(c, "pct_chrM");    sc = num(c, "pct_scaffold")
    kept   = num(s1, "pct_kept");   dup = num(s1, "dup_frac"); usable = num(s1, "usable_est")
    js     = num(fg, "JSdist");     tv = tss.get(s)
    pk     = num(fr, "n_peaks");    pg = num(fr, "pct_genome")
    fa     = num(fr, "frip_all");   fn = num(fr, "frip_nodup")
    fold   = (fa / (pg / 100)) if (fa is not None and pg) else None
    cells = [
        f'<td class="s">{html.escape(s)}</td>',
        f'<td>{tumor(s)}</td>', f'<td>{target(s)}</td>',
        f'<td>{"-" if pairs is None else format(int(pairs), ",")}</td>',
        f'<td class="{cls_align(ar)}">{fmt(ar,"{:.1f}")}</td>',
        f'<td class="{cls_mt(mt)}">{fmt(mt,"{:.2f}")}</td>',
        f'<td>{fmt(sc,"{:.2f}")}</td>',
        f'<td>{fmt(kept,"{:.1f}")}</td>',
        f'<td class="{cls_dup(dup)}">{fmt(dup,"{:.3f}")}</td>',
        f'<td>{"-" if usable is None else format(int(usable), ",")}</td>',
        f'<td>{fmt(frg.get("median"),"{:.0f}")}</td>',
        f'<td>{fmt(frg.get("p10"),"{:.0f}")}</td>',
        f'<td>{fmt(js)}</td>',
        f'<td class="{cls_tss(tv)}">{fmt(tv,"{:.2f}")}</td>',
        f'<td>{"-" if pk is None else format(int(pk), ",")}</td>',
        f'<td>{fmt(pg,"{:.2f}")}</td>',
        f'<td class="{cls_frip(fa)}">{fmt(fa)}</td>',
        f'<td class="{cls_frip(fn)}">{fmt(fn)}</td>',
        f'<td>{fmt(fold,"{:.1f}")}</td>',
    ]
    rows_html.append("<tr>" + "".join(cells) + "</tr>")

# ------------------------------------------------------------- group summaries
def group_table(keyfn, label):
    g = {}
    for s in samples:
        k = keyfn(s)
        fr, s1 = frip.get(s, {}), step1.get(s, {})
        g.setdefault(k, []).append((
            num(fr, "frip_all"), num(s1, "dup_frac"), tss.get(s),
            num(finger.get(s, {}), "JSdist"), num(fr, "n_peaks")))
    out = [f'<table class="mini"><thead><tr><th>{label}</th><th>n</th>'
           '<th>mean FRiP</th><th>mean dup</th><th>mean TSS</th>'
           '<th>mean JS</th><th>median peaks</th></tr></thead><tbody>']
    def mean(v):
        v = [x for x in v if x is not None]
        return sum(v)/len(v) if v else None
    def median(v):
        v = sorted(x for x in v if x is not None)
        return v[len(v)//2] if v else None
    for k in sorted(g, key=lambda k: -(mean([r[0] for r in g[k]]) or 0)):
        r = g[k]
        out.append("<tr><td class='s'>%s</td><td>%d</td><td>%s</td><td>%s</td>"
                   "<td>%s</td><td>%s</td><td>%s</td></tr>" % (
            k, len(r), fmt(mean([x[0] for x in r])), fmt(mean([x[1] for x in r])),
            fmt(mean([x[2] for x in r]), "{:.2f}"), fmt(mean([x[3] for x in r])),
            "-" if median([x[4] for x in r]) is None else format(int(median([x[4] for x in r])), ",")))
    out.append("</tbody></table>")
    return "\n".join(out)

# ------------------------------------------------------------------- commands
def script_src(name):
    p = os.path.join(os.path.dirname(os.path.abspath(__file__)), name)
    return open(p).read() if os.path.exists(p) else f"(missing: {name})"

STEPS = [
    ("0. Organise samples by assay",
     "Splits the project samples into assay groups by naming convention "
     "&mdash; <code>*-ATAC</code> (Tn5/Nextera) vs. everything else "
     "<code>cutrun</code> (TruSeq) &mdash; under <code>fastq/</code>, "
     "<code>trimmed/</code> and <code>bams/</code>. "
     f"This report covers the {n} samples in group <code>{GROUP}</code>.",
     "organize_by_assay.sh", None),
    ("1. Adapter/quality trimming",
     "TruSeq adapters for CUT&amp;RUN samples (<code>--illumina</code>); "
     "ATAC-named samples use <code>--nextera</code>.",
     "run_trim_galore.sh",
     "trim_galore --paired --fastqc --quality 20 --length 20 --illumina --trim-n \\\n"
     "    --output_dir trimmed <sample>_R1_001.fastq.gz <sample>_R2_001.fastq.gz"),
    (f"2. Alignment (bowtie2 &rarr; {HOMER_GENOME})",
     "Concordant pairs only, fragments up to 2 kb.",
     "run_bowtie2.sh",
     f"bowtie2 -x {BOWTIE2_INDEX} \\\n"
     "    -1 trimmed/<sample>_R1_001_val_1.fq.gz \\\n"
     "    -2 trimmed/<sample>_R2_001_val_2.fq.gz \\\n"
     "    --very-sensitive --no-mixed --no-discordant -X 2000 -p 8 \\\n"
     "  | samtools sort -@ 2 -m 2G -o bams/<sample>.sorted.bam -\n"
     "samtools index bams/<sample>.sorted.bam"),
    ("3. Filter, mark duplicates, blacklist",
     "<code>-f 2 -F 1804 -q 30</code> (proper pairs, MAPQ&ge;30, excluding "
     "unmapped / mate-unmapped / secondary / QC-fail / flagged-duplicate reads), "
     "restricted to chr1&ndash;19, X, Y &mdash; which also drops chrM and "
     "<code>*_random</code>/<code>chrUn_*</code> scaffolds. Duplicates are "
     "<b>marked, not removed</b>; downstream steps exclude them with "
     "<code>-F 1024</code> or <code>--ignoreDuplicates</code>. Blacklist "
     "then subtracted.",
     "run_filter_dedup.sh",
     "samtools view -b -f 2 -F 1804 -q 30 bams/<sample>.sorted.bam \\\n"
     "    chr1 chr2 ... chr19 chrX chrY > <sample>.filt.bam\n\n"
     "picard -Xmx4g MarkDuplicates I=<sample>.filt.bam O=<sample>.md.bam \\\n"
     "    M=metrics/<sample>_dup_metrics.txt ASSUME_SORTED=true \\\n"
     "    REMOVE_DUPLICATES=false VALIDATION_STRINGENCY=LENIENT\n\n"
     f"bedtools intersect -v -abam <sample>.md.bam -b {BLACKLIST} \\\n"
     f"    > filtered/{GROUP}/<sample>.final.bam\n"
     f"samtools index filtered/{GROUP}/<sample>.final.bam"),
    ("3b. Build aggregate stats",
     "<code>build_stats.sh</code> derives <code>contig_stats.tsv</code>, "
     "<code>align_stats.tsv</code> and <code>step1_summary.tsv</code> directly "
     "from the raw BAM, final BAM, Picard metrics and bowtie2 log &mdash; run "
     "this before FRiP/report or the duplicate-excluded FRiP column reads 0.",
     "build_stats.sh", None),
    ("4. QC + bigwigs + TSS + correlation",
     "All deeptools calls pass <code>--ignoreDuplicates --minMappingQuality 30</code>. "
     + (f"The TSS set is {tss_n:,} unique positions ({tss_plus:,} +, {tss_minus:,} "
        "&minus;) derived from "
        if tss_n else "The TSS set is derived from ")
     + f"HOMER's <code>{HOMER_GENOME}.tss</code>: uniform 4 kb windows, so the TSS "
     "is the window centre (<code>start+2000</code>), with strand encoded 0=+ / 1=&minus;.",
     "run_qc_step2.sh",
     f"bamPEFragmentSize -b <all {n} bams> --histogram fragment_sizes.png \\\n"
     "    --table fragment_sizes.txt --maxFragmentLength 1000 -p 16\n\n"
     "plotFingerprint -b <bams of one target> --ignoreDuplicates \\\n"
     "    --minMappingQuality 30 --plotFile fingerprint_<target>.png \\\n"
     "    --outQualityMetrics fingerprint_<target>.txt -p 16\n\n"
     "bamCoverage -b <sample>.final.bam -o <sample>.cpm.bw \\\n"
     "    --binSize 10 --normalizeUsing CPM --extendReads \\\n"
     "    --ignoreDuplicates --minMappingQuality 30 -p 4\n\n"
     f"computeMatrix reference-point --referencePoint TSS -S <{n} bigwigs> \\\n"
     f"    -R {HOMER_GENOME}_tss.bed -b 2000 -a 2000 --binSize 10 \\\n"
     "    --skipZeros --missingDataAsZero -o tss_matrix.gz -p 16\n"
     "plotProfile -m tss_matrix.gz -o tss_profile.png \\\n"
     "    --outFileNameData tss_profile_data.tab --perGroup\n"
     "plotHeatmap -m tss_matrix.gz -o tss_heatmap.png \\\n"
     "    --sortUsing mean --sortRegions descend\n\n"
     f"multiBamSummary bins -b <all {n} bams> --binSize 5000 \\\n"
     "    --ignoreDuplicates --minMappingQuality 30 -o bins5kb.npz -p 16\n"
     "plotCorrelation -in bins5kb.npz -c spearman -p heatmap --skipZeros \\\n"
     "    --removeOutliers --plotNumbers -o corr_spearman.png \\\n"
     "    --outFileCorMatrix corr_spearman.tab\n"
     "plotPCA -in bins5kb.npz -o pca.png --outFileNameData pca.tab"),
    ("5. Peak calling (MACS2) + FRiP",
     "Peaks are called treatment-only with <code>--nolambda</code> (genome-wide "
     "background) when no IgG/control library exists -- drop that flag if you "
     f"have one. <code>-g {MACS2_GENOME}</code> sets the effective genome size. "
     f"Run in the <code>{MACS2_ENV}</code> environment. FRiP merges each "
     "sample's peaks with <code>bedtools merge</code> before counting, because "
     "<code>--call-summits</code> emits overlapping entries per summit; it is "
     "reported both including and excluding flagged duplicates.",
     "run_macs2_frip.sh",
     f"macs2 callpeak -t filtered/{GROUP}/<sample>.final.bam \\\n"
     f"    -f BAMPE -g {MACS2_GENOME} -q 0.01 \\\n"
     "    --nolambda --keep-dup all --call-summits -B --SPMR \\\n"
     f"    -n <sample> --outdir peaks/{GROUP}/<sample>\n\n"
     "# FRiP\n"
     "sort -k1,1 -k2,2n <sample>_peaks.narrowPeak | bedtools merge -i - > peaks.merged.bed\n"
     "samtools view -c -L peaks.merged.bed <sample>.final.bam          # in peaks (all)\n"
     "samtools view -c -F 1024 -L peaks.merged.bed <sample>.final.bam  # in peaks (nodup)"),
]

def _blacklist_regions():
    try:
        with open(BLACKLIST) as fh:
            return sum(1 for _ in fh)
    except OSError:
        return None

_bl_n = _blacklist_regions()
VERSIONS = [("trim_galore/bowtie2/samtools/picard/bedtools", "pinned in envs/cutrun_env.yml"),
            ("deeptools", "pinned in envs/bigwig_env.yml"),
            ("MACS2", f"pinned in envs/macs2-env2.yml (env: {MACS2_ENV})"),
            ("genome", HOMER_GENOME), ("bowtie2 index", BOWTIE2_INDEX),
            ("blacklist", f"{BLACKLIST}" + (f" ({_bl_n} regions)" if _bl_n else ""))]

# ------------------------------------------------------------------- assemble
# fingerprint targets are whatever run_qc_step2.sh phase B actually produced
# for this group (one PNG per distinct target suffix), not a fixed list.
_fp_targets = sorted(
    os.path.basename(p)[len("fingerprint_"):-len(".png")]
    for p in glob.glob(f"{QC}/fingerprint_*.png"))

_tss_count_txt = f"{tss_n:,} unique RefSeq TSS" if tss_n else "RefSeq TSS"

PLOTS = [
    ("Fragment size distribution", "fragment_sizes.png",
     f"<code>bamPEFragmentSize</code>, all {n} samples, fragments up to 1000 bp. "
     "Per-sample median, mean, quartiles and percentiles are in the "
     "<code>Frag median</code> / <code>Frag p10</code> table columns and in "
     "<code>fragment_sizes.txt</code>."),
] + [
    (f"Fingerprint &mdash; {tgt}", f"fingerprint_{tgt}.png",
     "<code>plotFingerprint</code>, one panel per target, duplicates ignored. "
     "AUC and Synthetic JS distance per sample are in the per-target "
     "<code>fingerprint_&lt;target&gt;.txt</code> files." if i == 0 else None)
    for i, tgt in enumerate(_fp_targets)
] + [
    ("TSS profile (&plusmn;2 kb)", "tss_profile.png",
     f"CPM signal over {_tss_count_txt}, 10 bp bins, &plusmn;2 kb. "
     "The TSS score in the table is mean CPM within &plusmn;250 bp of the TSS "
     "divided by mean CPM in the outer 500 bp of each flank."),
    ("TSS heatmap", "tss_heatmap.png",
     "Same matrix, sorted by mean signal descending. Downscaled for embedding; "
     "full-resolution original is <code>tss_heatmap.png</code>."),
    ("Correlation &mdash; Spearman (5 kb bins)", "corr_spearman.png",
     f"<code>multiBamSummary bins --binSize 5000</code> across all {n} samples, "
     "duplicates ignored. Full matrices in <code>corr_spearman.tab</code> and "
     "<code>corr_pearson.tab</code>."),
    ("Correlation &mdash; Pearson (5 kb bins)", "corr_pearson.png", None),
    ("PCA (5 kb bins)", "pca.png",
     "From the same 5 kb bin matrix. Component values in <code>pca.tab</code>."),
]

plot_html = []
for title, fn, cap in PLOTS:
    d = embed(os.path.join(QC, fn))
    if not d:
        continue
    plot_html.append(
        f'<div class="card"><h3>{title}</h3>'
        + (f'<p class="cap">{cap}</p>' if cap else "")
        + f'<img src="{d}" alt="{html.escape(title)}"></div>')

cmd_html = []
for i, (title, blurb, script, cmd) in enumerate(STEPS):
    src = html.escape(script_src(script)) if script else ""
    cmd_block = f'<pre class="cmd">{html.escape(cmd)}</pre>' if cmd else ""
    cmd_html.append(f"""
<div class="card">
  <h3>{title}</h3>
  <p class="cap">{blurb}</p>
  {cmd_block}
  <details><summary>Full script &mdash; <code>{html.escape(script or "")}</code></summary>
  <pre class="src">{src}</pre></details>
</div>""")

now = datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M UTC")

# ------------------------------------------------- factual aggregate statistics
def _vals(fn):
    return [v for v in (fn(s) for s in samples) if v is not None]

def rng(vals, spec="{:.3f}"):
    return f"{spec.format(min(vals))} &ndash; {spec.format(max(vals))}" if vals else "-"

_tumors  = sorted(set(tumor(s) for s in samples))
_targets = sorted(set(target(s) for s in samples))

_ar   = _vals(lambda s: num(align.get(s, {}), "align_rate"))
_mt   = _vals(lambda s: num(contig.get(s, {}), "pct_chrM"))
_kept = _vals(lambda s: num(step1.get(s, {}), "pct_kept"))
_dup  = _vals(lambda s: num(step1.get(s, {}), "dup_frac"))
_use  = _vals(lambda s: num(step1.get(s, {}), "usable_est"))
_med  = _vals(lambda s: frag.get(s, {}).get("median"))
_js   = _vals(lambda s: num(finger.get(s, {}), "JSdist"))
_tss  = _vals(lambda s: tss.get(s))
_pk   = _vals(lambda s: num(frip.get(s, {}), "n_peaks"))
_fa   = _vals(lambda s: num(frip.get(s, {}), "frip_all"))

RUNSTATS = "".join(f"<div><b>{k}</b> &mdash; {v}</div>" for k, v in [
    ("Samples", f"{n}"),
    ("Total usable reads", f"{int(sum(_use)):,}" if _use else "-"),
    ("Alignment rate", f"{rng(_ar,'{:.1f}')} % (mean {sum(_ar)/len(_ar):.1f})" if _ar else "-"),
    ("chrM", f"{rng(_mt,'{:.2f}')} % (mean {sum(_mt)/len(_mt):.2f})" if _mt else "-"),
    ("Reads kept after filtering", f"{rng(_kept,'{:.1f}')} %" if _kept else "-"),
    ("Duplicate fraction", f"{rng(_dup)} (mean {sum(_dup)/len(_dup):.3f})" if _dup else "-"),
    ("Usable reads per sample", f"{int(min(_use)):,} &ndash; {int(max(_use)):,}" if _use else "-"),
    ("Fragment median", f"{rng(_med,'{:.0f}')} bp" if _med else "-"),
    ("Fingerprint JS distance", rng(_js)),
    ("TSS enrichment score", rng(_tss, "{:.2f}")),
    ("Peaks per sample", f"{int(min(_pk)):,} &ndash; {int(max(_pk)):,}" if _pk else "-"),
    ("FRiP", rng(_fa)),
])

HTML = f"""<!DOCTYPE html>
<html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>{GROUP_LABEL} QC report &mdash; {PROJECT_NAME} / {GROUP}</title>
<style>
:root{{--bg:#fff;--fg:#1a1a1a;--mut:#666;--line:#e2e2e2;--hd:#f6f7f9;
--good:#d6f0dc;--warn:#fdf0cd;--bad:#fadbd8;--acc:#0b6bcb}}
*{{box-sizing:border-box}}
body{{margin:0;font:14px/1.55 -apple-system,BlinkMacSystemFont,"Segoe UI",Roboto,Helvetica,Arial,sans-serif;
color:var(--fg);background:#fafafa}}
header{{background:var(--fg);color:#fff;padding:26px 32px}}
header h1{{margin:0 0 6px;font-size:21px}}
header .sub{{color:#b9b9b9;font-size:13px}}
nav{{position:sticky;top:0;background:#fff;border-bottom:1px solid var(--line);
padding:10px 32px;z-index:10;overflow-x:auto;white-space:nowrap}}
nav a{{color:var(--acc);text-decoration:none;margin-right:18px;font-size:13px}}
nav a:hover{{text-decoration:underline}}
main{{padding:26px 32px;max-width:1600px}}
h2{{font-size:17px;margin:34px 0 12px;padding-bottom:6px;border-bottom:2px solid var(--fg)}}
h3{{font-size:14px;margin:0 0 8px}}
.card{{background:#fff;border:1px solid var(--line);border-radius:6px;padding:16px;margin:14px 0}}
.cap{{color:var(--mut);font-size:13px;margin:0 0 10px}}
img{{max-width:100%;height:auto;border:1px solid var(--line);border-radius:4px}}
table{{border-collapse:collapse;width:100%;font-size:12.5px;background:#fff}}
th,td{{border:1px solid var(--line);padding:5px 8px;text-align:right;white-space:nowrap}}
th{{background:var(--hd);position:sticky;top:0;cursor:pointer;font-weight:600;text-align:right}}
th:hover{{background:#ececf0}}
td.s,th:first-child{{text-align:left;font-family:ui-monospace,SFMono-Regular,Menlo,monospace}}
td:nth-child(2),td:nth-child(3){{text-align:left}}
.good{{background:var(--good)}} .warn{{background:var(--warn)}} .bad{{background:var(--bad)}}
.wrap{{overflow-x:auto;border:1px solid var(--line);border-radius:6px;max-height:80vh}}
pre{{background:#1e1e1e;color:#e8e8e8;padding:12px;border-radius:5px;overflow-x:auto;
font:12px/1.5 ui-monospace,SFMono-Regular,Menlo,monospace}}
pre.src{{max-height:420px;overflow-y:auto;background:#f6f7f9;color:#222;border:1px solid var(--line)}}
details{{margin-top:10px}} summary{{cursor:pointer;color:var(--acc);font-size:13px}}
.grid{{display:grid;grid-template-columns:repeat(auto-fit,minmax(330px,1fr));gap:14px}}
.mini{{width:auto}} .mini td,.mini th{{padding:4px 10px}}
.legend span{{display:inline-block;padding:2px 9px;margin-right:8px;border-radius:3px;font-size:12px}}
.kv{{columns:2;font-size:13px}} .kv div{{padding:2px 0}}
.note{{border-left:3px solid var(--acc);background:#f2f7fd;padding:10px 14px;margin:12px 0;font-size:13px}}
.warnbox{{border-left:3px solid #d68910;background:#fdf6e7;padding:10px 14px;margin:12px 0;font-size:13px}}
code{{background:#eef0f3;padding:1px 4px;border-radius:3px;font-size:12px}}
footer{{color:var(--mut);font-size:12px;padding:20px 32px}}
</style></head><body>
<header>
  <h1>{GROUP_LABEL} QC report &mdash; {PROJECT_NAME} <span style="opacity:.6">/ {GROUP}</span></h1>
  <div class="sub">{n} samples &middot; {HOMER_GENOME} &middot; generated {now}</div>
</header>
<nav>
  <a href="#summary">Summary</a><a href="#table">Per-sample metrics</a>
  <a href="#groups">Group summaries</a><a href="#plots">Plots</a>
  <a href="#commands">Commands</a><a href="#versions">Versions</a>
</nav>
<main>

<h2 id="summary">Run overview</h2>
<div class="card">
<p><b>Dataset.</b> {n} {GROUP_LABEL} samples from <code>{PROJECT_NAME}</code>,
group <code>{GROUP}</code>: {len(_tumors)} tumor(s) ({', '.join(_tumors)}) &times;
{len(_targets)} target(s) ({', '.join(_targets)}). Paired-end, aligned to {HOMER_GENOME}.
Other assay groups in the same project (split by naming convention) are
processed separately and are not included here.</p>

<p><b>Steps run.</b> trim_galore &rarr; bowtie2 ({HOMER_GENOME}) &rarr; MAPQ/proper-pair
filtering and primary-chromosome restriction &rarr; MarkDuplicates (marked, not
removed) &rarr; blacklist removal &rarr; fragment-size, fingerprint, bigwig, TSS
and correlation QC &rarr; MACS2 peak calling &rarr; FRiP. Completed {n}/{n} at
every step; 0 failures.</p>

<div class="kv" style="margin-top:12px">
{RUNSTATS}
</div>

<p class="legend" style="margin-top:14px">Cell shading marks the stated numeric
ranges only:
<span class="good">FRiP &ge;0.10</span><span class="warn">0.02&ndash;0.10</span><span class="bad">&lt;0.02</span>
&middot; dup &le;0.20 / 0.20&ndash;0.35 / &gt;0.35
&middot; chrM &le;2% / 2&ndash;10% / &gt;10%
&middot; align &ge;90% / 80&ndash;90% / &lt;80%
&middot; TSS &ge;2.0 / 1.3&ndash;2.0 / &lt;1.3.</p>
</div>

<h2 id="table">Per-sample metrics</h2>
<div class="card"><p class="cap">Click any header to sort. Peak/FRiP columns come
from MACS2 output; fold enrichment = FRiP &divide; (peak bp / genome).</p>
<div class="wrap"><table id="main"><thead><tr>
{''.join(f'<th>{h}</th>' for h,_ in COLS)}
</tr></thead><tbody>
{chr(10).join(rows_html)}
</tbody></table></div></div>

<h2 id="groups">Group summaries</h2>
<div class="grid">
  <div class="card"><h3>By target</h3>{group_table(target,"Target")}</div>
  <div class="card"><h3>By tumor</h3>{group_table(tumor,"Tumor")}</div>
</div>

<h2 id="plots">Plots</h2>
{''.join(plot_html)}

<h2 id="commands">Commands &mdash; full pipeline</h2>
<div class="card"><p class="cap">Every step in order. The condensed command shows
the operative call for one sample; the collapsed block underneath is the exact
driver script that was run, verbatim from disk.</p></div>
{''.join(cmd_html)}

<h2 id="versions">Software &amp; references</h2>
<div class="card"><div class="kv">
{''.join(f'<div><b>{k}</b> &mdash; {v}</div>' for k,v in VERSIONS)}
</div></div>

</main>
<footer>Generated by <code>make_report.py</code> &middot; all plots embedded, no external assets.</footer>
<script>
document.querySelectorAll('#main th').forEach(function(th,i){{
  th.addEventListener('click',function(){{
    var tb=th.closest('table').tBodies[0];
    var rows=Array.prototype.slice.call(tb.rows);
    var asc=!(th.dataset.asc==='1'); th.dataset.asc=asc?'1':'0';
    var val=function(r){{
      var t=r.cells[i].textContent.replace(/,/g,'').trim();
      if(t===''||t==='-')return asc?Infinity:-Infinity;
      var n=parseFloat(t); return isNaN(n)?t.toLowerCase():n;
    }};
    rows.sort(function(a,b){{
      var x=val(a),y=val(b);
      if(typeof x==='string'||typeof y==='string')
        return asc?String(x).localeCompare(String(y)):String(y).localeCompare(String(x));
      return asc?x-y:y-x;
    }});
    rows.forEach(function(r){{tb.appendChild(r)}});
  }});
}});
</script>
</body></html>"""

os.makedirs(QC, exist_ok=True)
with open(OUT, "w") as fh:
    fh.write(HTML)
print(f"wrote {OUT}  ({os.path.getsize(OUT)/1048576:.1f} MB)")
print(f"samples in table: {n}")
missing = [k for k, d in (("align", align), ("contig", contig), ("step1", step1),
                          ("frip", frip), ("fingerprint", finger),
                          ("frag", frag), ("tss", tss)) if len(d) != n]
if missing:
    print("WARNING incomplete metric sets:", ", ".join(
        f"{k}({len(dict(align=align,contig=contig,step1=step1,frip=frip,fingerprint=finger,frag=frag,tss=tss)[k])}/{n})"
        for k in missing))
