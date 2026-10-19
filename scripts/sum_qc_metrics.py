#!/usr/bin/env python3
"""
Summarize per-sample RNA-seq QC metrics from a RustQC / MultiQC output directory
into a single one-row table (<sample>.qcsum.txt).

This is the replacement for the old `sum_qc_metrics.r`, which parsed raw
Picard / STAR / RNA-SeQC / cal_q30 / coverage.stat files. The new RNA-seq QC
workflow (FastQC, STAR, samtools, Qualimap-rnaseq, RSeQC, dupRadar, preseq,
featureCounts, fastp) aggregates most metrics under multiqc_data/, so we read
those parsed tables instead. RSEM (genes.results), ERCC and SIRV spike-in files
are still produced (symlinked in) and parsed directly, exactly as before.

Only Python standard library is used; `samtools` is required on PATH for the
mitochondrial and ribosomal read counts (it reads the marked-duplicates BAM).

Metric provenance / notable changes vs the old workflow:
  * Q30           -> fastp.json summary (combined R1+R2; raw & trimmed) plus
                     per-read MEAN quality from FastQC. The old per-base
                     PCT_R1/R2_Q30 (cal_q30) is no longer produced.
  * RAW/TRIMMED   -> fastp.json (new: not available in the old workflow).
  * EXON/INTRON/INTERGENIC/5'-3' bias/ambiguous -> Qualimap rnaseq
                     (different tool than the old RNA-SeQC/Picard; values are
                     NOT backward-comparable to historical qcsum.txt files).
  * PCT_RIBOSOMAL -> Picard-style: reads overlapping the rRNA interval list,
                     counted against the BAM (featureCounts %rRNA is ~0 because
                     the Ensembl biotype GTF misses the 45S/rDNA repeats that
                     Picard's interval list covers).
  * NUM_JUNCS     -> RSeQC total_splicing_junctions (old was STAR SJ.out.tab
                     junctions with >=5 reads; different definition).
  * DUPLICATION   -> Picard MarkDuplicates *.duplicate.metrics (PERCENT_DUPLICATION).
  * Inputs are located under --indir by file-name suffix (the sample id in the
    name need not match -s exactly); MultiQC tables live in <sample>multiqc_data/.
  * Dropped       -> AVE_COV, PCT_1X/10X/50X_COV, ESTIMATED_LIBRARY_SIZE
                     (no equivalent in the new outputs).
  * Added         -> TIN, inner distance (RSeQC), strandedness, CHRY_PERC /
                     SEX_PRED (chrY expected-count fraction from RSEM, used to
                     predict sex by a simple threshold; XIST_TPM is the
                     second, independent signal, and SEX_PRED is Unclear
                     when the two disagree).
"""

import argparse
import csv
import glob
import json
import math
import os
import re
import subprocess
import sys
import tempfile

NA = "NA"


def warn(msg):
    sys.stderr.write("[sum_qc_metrics] " + msg + "\n")


# --------------------------------------------------------------------------- #
# Generic helpers
# --------------------------------------------------------------------------- #
def read_mqc_table(path):
    """Read a MultiQC tab-delimited table into {sample: {column: value}}."""
    if not path or not os.path.exists(path):
        warn("missing table: %s" % path)
        return {}
    rows = {}
    with open(path, newline="") as fh:
        reader = csv.reader(fh, delimiter="\t")
        header = next(reader, None)
        if not header:
            return {}
        for rec in reader:
            if not rec:
                continue
            rows[rec[0]] = dict(zip(header[1:], rec[1:]))
    return rows


def pick_row(table, sample):
    """Find a sample's row. Tolerates the '.mdup' suffix some tools add; if no
    name matches but the table holds a single sample (the per-sample case),
    use that sole row rather than requiring an exact sample-name match."""
    for key in (sample, sample + ".mdup"):
        if key in table:
            return table[key]
    if len(table) == 1:
        only = next(iter(table))
        warn("sample %r not matched; using sole row %r" % (sample, only))
        return table[only]
    return None


def to_float(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return None


def rnd(x, n=2):
    return round(x, n) if isinstance(x, (int, float)) else NA


# --------------------------------------------------------------------------- #
# fastp (raw / trimmed reads, Q30)
# --------------------------------------------------------------------------- #
def parse_fastp(path):
    out = dict.fromkeys(
        ["RAW_READS", "READS_AFTER_TRIM", "PCT_READS_PASS",
         "Q30_RATE_RAW", "Q30_RATE_TRIM", "MEAN_LEN_AFTER", "GC_PCT"], NA)
    if not path or not os.path.exists(path):
        warn("missing fastp json: %s" % path)
        return out
    with open(path) as fh:
        d = json.load(fh)
    before = d.get("summary", {}).get("before_filtering", {})
    after = d.get("summary", {}).get("after_filtering", {})
    raw = before.get("total_reads")
    trimmed = after.get("total_reads")
    out["RAW_READS"] = int(raw) if raw is not None else NA
    out["READS_AFTER_TRIM"] = int(trimmed) if trimmed is not None else NA
    if raw:
        out["PCT_READS_PASS"] = rnd(100.0 * trimmed / raw)
    if before.get("q30_rate") is not None:
        out["Q30_RATE_RAW"] = rnd(100.0 * before["q30_rate"])
    if after.get("q30_rate") is not None:
        out["Q30_RATE_TRIM"] = rnd(100.0 * after["q30_rate"])
    # mean read length after filtering (average of R1/R2 means)
    lens = [after.get("read1_mean_length"), after.get("read2_mean_length")]
    lens = [x for x in lens if x is not None]
    if lens:
        out["MEAN_LEN_AFTER"] = rnd(sum(lens) / len(lens), 1)
    if after.get("gc_content") is not None:
        out["GC_PCT"] = rnd(100.0 * after["gc_content"])
    return out


# --------------------------------------------------------------------------- #
# STAR alignment counts
# --------------------------------------------------------------------------- #
def parse_star(path, sample):
    out = {"TOTAL_READS": NA, "UNIQ_ALIGNED": NA, "PCT_UNIQ_ALIGNED": NA}
    row = pick_row(read_mqc_table(path), sample)
    if not row:
        warn("no STAR row for %s" % sample)
        return out
    tot = to_float(row.get("total_reads"))
    uniq = to_float(row.get("uniquely_mapped"))
    pct = to_float(row.get("uniquely_mapped_percent"))
    out["TOTAL_READS"] = int(tot) if tot is not None else NA
    out["UNIQ_ALIGNED"] = int(uniq) if uniq is not None else NA
    out["PCT_UNIQ_ALIGNED"] = rnd(pct) if pct is not None else NA
    return out


# --------------------------------------------------------------------------- #
# Qualimap rnaseq_qc_results.txt (exon/intron/intergenic, bias, ambiguous)
# --------------------------------------------------------------------------- #
def parse_qualimap(path):
    out = dict.fromkeys(
        ["PCT_EXON", "PCT_INTRON", "PCT_INTERGENIC",
         "AMBIGUOUS_ALN", "MEDIAN_5PRIME_TO_3PRIME_BIAS"], NA)
    if not path or not os.path.exists(path):
        warn("missing qualimap results: %s" % path)
        return out
    txt = open(path).read()

    def grab_pct(label):
        m = re.search(label + r"\s*=\s*[\d,]+\s*\(([\d.]+)%\)", txt)
        return rnd(float(m.group(1))) if m else NA

    def grab_num(label):
        m = re.search(label + r"\s*=\s*([\d,]+)", txt)
        return float(m.group(1).replace(",", "")) if m else None

    out["PCT_EXON"] = grab_pct("exonic")
    out["PCT_INTRON"] = grab_pct("intronic")
    out["PCT_INTERGENIC"] = grab_pct("intergenic")

    ambiguous = grab_num("ambiguous alignments")
    total = grab_num("total alignments")
    if ambiguous is not None and total:
        out["AMBIGUOUS_ALN"] = rnd(100.0 * ambiguous / total)

    m = re.search(r"5'-3' bias\s*=\s*([\d.]+)", txt)
    if m:
        out["MEDIAN_5PRIME_TO_3PRIME_BIAS"] = rnd(float(m.group(1)))
    return out


# --------------------------------------------------------------------------- #
# Picard MarkDuplicates metrics (duplication)
# --------------------------------------------------------------------------- #
def parse_duplicate_metrics(path):
    """PERCENT_DUPLICATION from a Picard MarkDuplicates metrics file. The data
    is the single row beneath the 'LIBRARY\t...' header line."""
    out = {"PERCENT_DUPLICATION": NA}
    if not path or not os.path.exists(path):
        warn("missing duplicate metrics: %s" % path)
        return out
    with open(path) as fh:
        for line in fh:
            if line.startswith("LIBRARY\t"):
                header = line.rstrip("\n").split("\t")
                data = next(fh, "").rstrip("\n").split("\t")
                row = dict(zip(header, data))
                pd = to_float(row.get("PERCENT_DUPLICATION"))
                out["PERCENT_DUPLICATION"] = rnd(100.0 * pd) if pd is not None else NA
                break
    return out


# --------------------------------------------------------------------------- #
# samtools stats (reads_mapped, for the ribosomal denominator)
# --------------------------------------------------------------------------- #
def get_reads_mapped(path, sample):
    row = pick_row(read_mqc_table(path), sample)
    if not row:
        warn("no samtools stats row for %s (reads_mapped)" % sample)
        return None
    return to_float(row.get("reads_mapped"))


# --------------------------------------------------------------------------- #
# RSeQC: TIN, junctions, strandedness
# --------------------------------------------------------------------------- #
def parse_tin(path, sample):
    out = {"TIN_MEDIAN": NA, "TIN_MEAN": NA}
    row = pick_row(read_mqc_table(path), sample)
    if row:
        out["TIN_MEDIAN"] = rnd(to_float(row.get("TIN(median)")))
        out["TIN_MEAN"] = rnd(to_float(row.get("TIN(mean)")))
    return out


def parse_inner_distance(path, sample):
    """Mean/median inner distance from RSeQC's <sample>.inner_distance_mean.txt
    (columns Name, Mean, Median, sd). Reported as-is; values are negative when
    paired mates overlap, as in short RNA fragments."""
    out = {"INNER_DISTANCE_MEAN": NA, "INNER_DISTANCE_MEDIAN": NA}
    row = pick_row(read_mqc_table(path), sample)
    if row:
        mean = to_float(row.get("Mean"))
        median = to_float(row.get("Median"))
        out["INNER_DISTANCE_MEAN"] = rnd(mean) if mean is not None else NA
        out["INNER_DISTANCE_MEDIAN"] = rnd(median) if median is not None else NA
    return out


def parse_junctions(path, sample):
    out = {"NUM_JUNCS": NA, "KNOWN_JUNCS": NA, "NOVEL_JUNCS": NA}
    row = pick_row(read_mqc_table(path), sample)
    if row:
        for col, key in (("total_splicing_junctions", "NUM_JUNCS"),
                         ("known_splicing_junctions", "KNOWN_JUNCS"),
                         ("novel_splicing_junctions", "NOVEL_JUNCS")):
            v = to_float(row.get(col))
            out[key] = int(v) if v is not None else NA
    return out


def parse_strandedness(path, sample):
    out = {"STRANDEDNESS": NA, "PCT_ANTISENSE": NA}
    row = pick_row(read_mqc_table(path), sample)
    if not row:
        return out
    sense = to_float(row.get("pe_sense"))
    anti = to_float(row.get("pe_antisense"))
    if sense is None or anti is None:
        return out
    out["PCT_ANTISENSE"] = rnd(100.0 * anti)
    if max(sense, anti) < 0.6:
        out["STRANDEDNESS"] = "unstranded"
    elif anti > sense:
        out["STRANDEDNESS"] = "reverse"
    else:
        out["STRANDEDNESS"] = "forward"
    return out


# --------------------------------------------------------------------------- #
# samtools-based counts: MT %, ribosomal %
# --------------------------------------------------------------------------- #
def samtools_view_count(bam, options=None, region=None, threads=1):
    # samtools view [options] <bam> [region]: options precede the BAM, a
    # region (e.g. "MT") must follow it. -@ parallelizes BGZF decompression.
    cmd = ["samtools", "view", "-c", "-F", "0x100", "-@", str(threads)]
    if options:
        cmd += options
    cmd.append(bam)
    if region:
        cmd.append(region)
    try:
        return int(subprocess.check_output(cmd).decode().strip())
    except (subprocess.CalledProcessError, FileNotFoundError, ValueError) as e:
        warn("samtools failed (%s): %s" % (" ".join(cmd), e))
        return None


def calc_mt_pct(bam, star_total_reads, mt_contig="MT"):
    # idxstats reads only the .bai index (no read decompression) -> very fast.
    if not bam or not os.path.exists(bam) or not star_total_reads:
        warn("cannot compute MT_PCT (bam or STAR total missing)")
        return {"MT_PCT": NA}
    try:
        out = subprocess.check_output(["samtools", "idxstats", bam]).decode()
    except (subprocess.CalledProcessError, FileNotFoundError) as e:
        warn("samtools idxstats failed: %s" % e)
        return {"MT_PCT": NA}
    mt = None
    for line in out.splitlines():
        f = line.split("\t")
        if len(f) >= 3 and f[0] == mt_contig:
            mt = int(f[2])  # mapped read-segments on the MT contig
            break
    if mt is None:
        warn("contig %r not found in idxstats" % mt_contig)
        return {"MT_PCT": NA}
    # old definition: MT alignments / (read pairs * 2 mates)
    return {"MT_PCT": rnd(100.0 * mt / (star_total_reads * 2))}


def interval_list_to_bed(interval_list, bed_path):
    """Convert a Picard .interval_list (1-based, inclusive) to BED (0-based)."""
    n = 0
    with open(interval_list) as fin, open(bed_path, "w") as fout:
        for line in fin:
            if line.startswith("@") or not line.strip():
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 3:
                continue
            chrom, start, end = f[0], int(f[1]), int(f[2])
            fout.write("%s\t%d\t%d\n" % (chrom, start - 1, end))
            n += 1
    return n


def calc_ribosomal_pct(bam, interval_list, reads_mapped, threads=1):
    if not bam or not os.path.exists(bam):
        warn("cannot compute PCT_RIBOSOMAL (bam missing)")
        return {"PCT_RIBOSOMAL": NA}
    if not interval_list or not os.path.exists(interval_list):
        warn("cannot compute PCT_RIBOSOMAL (rRNA interval list missing)")
        return {"PCT_RIBOSOMAL": NA}
    if not reads_mapped:
        warn("cannot compute PCT_RIBOSOMAL (reads_mapped missing)")
        return {"PCT_RIBOSOMAL": NA}
    with tempfile.NamedTemporaryFile("w", suffix=".bed", delete=False) as tf:
        bed = tf.name
    try:
        interval_list_to_bed(interval_list, bed)
        ribo = samtools_view_count(bam, options=["-L", bed], threads=threads)
    finally:
        os.unlink(bed)
    if ribo is None:
        return {"PCT_RIBOSOMAL": NA}
    return {"PCT_RIBOSOMAL": rnd(100.0 * ribo / reads_mapped)}


# --------------------------------------------------------------------------- #
# RSEM expression + ERCC + SIRV spike-ins
# --------------------------------------------------------------------------- #
def read_tsv_dicts(path):
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def parse_expression(path):
    out = {"TPM_0": NA, "TPM_1": NA, "TPM_10": NA, "TPM_50": NA}
    if not path or not os.path.exists(path):
        warn("missing RSEM genes.results: %s" % path)
        return out
    tpms = []
    for row in read_tsv_dicts(path):
        gid = row.get("gene_id", "")
        if "ERCC" in gid or "SIRV" in gid:
            continue
        t = to_float(row.get("TPM"))
        if t is not None:
            tpms.append(t)
    out["TPM_0"] = sum(1 for t in tpms if t > 0)
    out["TPM_1"] = sum(1 for t in tpms if t > 1)
    out["TPM_10"] = sum(1 for t in tpms if t > 10)
    out["TPM_50"] = sum(1 for t in tpms if t > 50)
    return out


XIST_ID = "ENSG00000229807"


def parse_chry_sex(path, chry_genes_path, threshold=0.01, xist_threshold=10.0):
    """Sex from two independent signals in the salmon -> RSEM gene table:

      CHRY_PERC = 100 * sum(expected_count of chrY genes) / sum(expected_count
                  of all genes); male evidence when > threshold.
      XIST_TPM  = XIST TPM; female evidence when > xist_threshold.

    SEX_PRED is Male (chrY high, XIST low), Female (chrY low, XIST high), or
    Unclear when the two disagree (both high: e.g. XXY or a mixed sample; both
    low). The chrY list must exclude the PAR genes: Ensembl 115 gives the chrY
    PAR copies their own gene IDs and salmon splits PAR reads between X and Y,
    which puts every female at ~0.03% chrY, above the cutoff. With PAR removed,
    females sit at ~0.0002% and males at ~0.04%. The salmon table has no ERCC/SIRV
    rows, so they are not in the denominator (shifts CHRY_PERC by <3%). gene_id
    version suffixes (ENSG....N) are stripped before matching. CHRY_PERC is NA
    when the chrY list is missing; SEX_PRED is NA unless both signals exist."""
    out = {"CHRY_PERC": NA, "XIST_TPM": NA, "SEX_PRED": NA}
    if not path or not os.path.exists(path):
        warn("missing RSEM genes.results for chrY/sex: %s" % path)
        return out
    chry = set()
    if not chry_genes_path or not os.path.exists(chry_genes_path):
        warn("missing chrY gene list: %s" % chry_genes_path)
    else:
        with open(chry_genes_path) as fh:
            chry = {ln.strip().split(".")[0] for ln in fh if ln.strip()}
        if not chry:
            warn("chrY gene list is empty: %s" % chry_genes_path)
    total = ychr = 0.0
    xist = None
    for row in read_tsv_dicts(path):
        gid = (row.get("gene_id") or "").split(".")[0]
        if gid == XIST_ID:
            xist = to_float(row.get("TPM"))
        ec = to_float(row.get("expected_count"))
        if ec is None:
            continue
        total += ec
        if gid in chry:
            ychr += ec
    perc = None
    if chry and total > 0:
        perc = 100.0 * ychr / total
        out["CHRY_PERC"] = rnd(perc, 4)
    elif chry:
        warn("zero total expected_count; CHRY_PERC = NA")
    if xist is None:
        warn("XIST (%s) not in %s; XIST_TPM = NA" % (XIST_ID, path))
    else:
        out["XIST_TPM"] = rnd(xist)
    if perc is None or xist is None:
        return out
    y_male, x_female = perc > threshold, xist > xist_threshold
    if y_male != x_female:
        out["SEX_PRED"] = "Male" if y_male else "Female"
    else:
        out["SEX_PRED"] = "Unclear"
        warn("sex evidence conflicts: CHRY_PERC=%.4f (cutoff %g), XIST_TPM=%.2f "
             "(cutoff %g)" % (perc, threshold, xist, xist_threshold))
    return out


def parse_globin(path):
    """GLOBIN_PERC = the percentage on the '#Matched' line of
    <sample>.globin_stats.txt (reads mapping to the globin transcript set, e.g.
    `#Matched<TAB>37164567<TAB>14.94934%`). NA when the file or that line is
    missing."""
    out = {"GLOBIN_PERC": NA}
    if not path or not os.path.exists(path):
        warn("missing globin stats: %s" % path)
        return out
    with open(path) as fh:
        for line in fh:
            if line.lstrip("#").lower().startswith("matched"):
                m = re.search(r"([\d.]+)\s*%", line)
                if m:
                    out["GLOBIN_PERC"] = rnd(float(m.group(1)))
                break
    return out


def parse_tissue_pred(path, sample):
    """TISSUE_PRED from <sample>.tissuePred.txt (the tissue-classifier output: a
    header row plus one data row; the 'predicted' column holds the predicted
    tissue, or 'uncertain/other'). NA when the file is missing."""
    out = {"TISSUE_PRED": NA}
    if not path or not os.path.exists(path):
        warn("missing tissue prediction: %s" % path)
        return out
    row = pick_row(read_mqc_table(path), sample)
    if not row:
        warn("no tissue prediction row for %s" % sample)
        return out
    out["TISSUE_PRED"] = row.get("predicted") or NA
    # purity.py fields, present since the purity module was added to predict.py.
    # Older tissuePred files lack them and score NA, so this stays backward
    # compatible with samples processed before the change.
    #
    # Only three of purity.py's eight fields are reported. The dropped five were
    # either near-constant (purity_ref_source is always "agreed" -- the Snakefile
    # passes no --score-as; purity_not_assessed is a property of the reference
    # tissue, not the sample), redundant with what is kept (purity_reference
    # equals the canonicalized submitted tissue whenever the source is "agreed"),
    # or too narrow to earn a column (skin_form, non-null only when a skin panel
    # trips, and already implied by which panel is named in purity_flags).
    #
    # PURITY_REASON is read but NOT emitted: qc_call() needs it to raise the
    # TISSUE_MISMATCH flag. Do not drop it here. Its other values stay recoverable
    # from the reported columns -- no_submitted leaves TISSUE_SUBMITTED as NA, and
    # unsupported_submitted leaves the offending free-text label visible there.
    for src, dst in (("submitted_tissue", "TISSUE_SUBMITTED"),
                     ("purity_call", "PURITY_CALL"),
                     ("purity_flags", "PURITY_FLAGS"),
                     ("purity_reason", "PURITY_REASON")):
        out[dst] = row.get(src) or NA
    return out


# Absolute QC cutoffs, applied identically to all studies. (metric, direction,
# WARN, FAIL); "min" means higher is better. TOTAL_READS counts aligned PAIRS, so
# the 40 M floor is roughly 80 M reads.
#
# GLOBIN_PERC is deliberately absent: high globin carry-over costs a blood library
# effective depth, but the reads that remain are ordinarily fine, so it is a
# covariate to check rather than a gate.
QC_RULES = [
    ("TOTAL_READS",                  "min", 50e6, 40e6),
    ("PCT_READS_PASS",               "min", 95,   90),
    ("Q30_RATE_TRIM",                "min", 90,   85),
    ("PCT_UNIQ_ALIGNED",             "min", 90,   85),
    ("PCT_EXON",                     "min", 70,   65),
    ("PCT_INTERGENIC",               "max", 3,    5),
    ("PCT_RIBOSOMAL",                "max", 2,    5),
    ("TIN_MEDIAN",                   "min", 75,   70),
    ("PERCENT_DUPLICATION",          "max", 60,   75),
    ("MEDIAN_5PRIME_TO_3PRIME_BIAS", "max", 1.4,  1.6),
    ("ERCC_R2",                      "min", 0.85, 0.80),
]


def qc_call(data):
    """{QC_CALL, QC_FLAGS} from the absolute cutoffs above.

    Any FAIL metric -> FAIL; two or more WARN metrics -> REVIEW; else PASS. A
    missing metric scores NA and never fails, which in practice means ERCC_R2
    only applies to spiked-in libraries.

    This is the absolute arm of the QC rule only. The report's second arm -- three
    or more metrics beyond 3.5 robust z from their within-tissue median -- needs a
    cohort to compute that median against, and this script sees one sample. On the
    296-library cohort the two agree on 289 calls, differing only on REVIEWs the
    outlier arm raised alone, so this is a faithful gate and a slightly permissive
    triage. The outlier arm stays in the report, where the cohort exists.
    """
    failed, warned = [], []
    for metric, direction, warn_at, fail_at in QC_RULES:
        value = data.get(metric)
        if value is None or value == NA:
            continue
        try:
            value = float(value)
        except (TypeError, ValueError):
            continue
        if (value < fail_at) if direction == "min" else (value > fail_at):
            failed.append(metric)
        elif (value < warn_at) if direction == "min" else (value > warn_at):
            warned.append(metric)

    flags = []
    if failed:
        flags.append("FAIL:" + ",".join(failed))
    if warned:
        flags.append("WARN:" + ",".join(warned))
    call = "FAIL" if failed else ("REVIEW" if len(warned) >= 2 else "PASS")

    # Submitted and predicted tissue disagreeing is not a sequencing-metric
    # problem, so it does not move QC_CALL -- but it is worth seeing. Across 365
    # libraries in the 2026-08 cohort there were zero true mismatches (the only
    # disagreements were tissues the four-class model has no label for), so this
    # firing means something happened, not that the check is noisy.
    if data.get("PURITY_REASON") == "tissue_mismatch":
        flags.append("TISSUE_MISMATCH:%s/%s" % (data.get("TISSUE_SUBMITTED", "?"),
                                                data.get("TISSUE_PRED", "?")))
    return {"QC_CALL": call, "QC_FLAGS": ";".join(flags)}


def linreg_r2(xs, ys):
    """R^2 of an ordinary least squares fit y ~ x (pure python)."""
    n = len(xs)
    if n < 2:
        return None
    mx, my = sum(xs) / n, sum(ys) / n
    sxx = sum((x - mx) ** 2 for x in xs)
    syy = sum((y - my) ** 2 for y in ys)
    sxy = sum((x - mx) * (y - my) for x, y in zip(xs, ys))
    if sxx == 0 or syy == 0:
        return None
    return (sxy * sxy) / (sxx * syy)


def parse_ercc(path):
    out = {"ERCC_TPM_SUM": NA, "ERCC_TPM_0": NA, "ERCC_TPM_1": NA, "ERCC_R2": NA}
    if not path or not os.path.exists(path):
        warn("missing ercc file: %s" % path)
        return out
    rows = read_tsv_dicts(path)
    tpm = [to_float(r.get("TPM")) or 0.0 for r in rows]
    conc = [to_float(r.get("conc")) for r in rows]
    if sum(tpm) == 0:
        return out
    out["ERCC_TPM_SUM"] = round(sum(tpm))
    out["ERCC_TPM_0"] = sum(1 for t in tpm if t > 0)
    out["ERCC_TPM_1"] = sum(1 for t in tpm if t >= 1)
    xs, ys = [], []
    for t, c in zip(tpm, conc):
        if t and t > 0 and c and c > 0:
            xs.append(math.log2(c))
            ys.append(math.log2(t))
    r2 = linreg_r2(xs, ys)
    out["ERCC_R2"] = rnd(r2) if r2 is not None else NA
    return out


def stdev(vals):
    n = len(vals)
    if n < 2:
        return None
    m = sum(vals) / n
    return math.sqrt(sum((v - m) ** 2 for v in vals) / (n - 1))


def parse_sirv(path):
    """Per SIRV gene: SD of each transcript's % of that gene's total TPM."""
    keys = ["SIRV_TPM_SUM", "SIRV1", "SIRV2", "SIRV3",
            "SIRV4", "SIRV5", "SIRV6", "SIRV7"]
    out = dict.fromkeys(keys, NA)
    if not path or not os.path.exists(path):
        warn("missing sirv file: %s" % path)
        return out
    rows = read_tsv_dicts(path)
    total = sum(to_float(r.get("TPM")) or 0.0 for r in rows)
    if total == 0:
        return out
    out["SIRV_TPM_SUM"] = round(total)
    by_gene = {}
    for r in rows:
        by_gene.setdefault(r.get("gene_id"), []).append(to_float(r.get("TPM")) or 0.0)
    for gene, tpms in by_gene.items():
        gsum = sum(tpms)
        if gsum == 0 or gene not in out:
            continue
        perc = [100.0 * t / gsum for t in tpms]
        sd = stdev(perc)
        out[gene] = rnd(sd) if sd is not None else NA
    return out


# --------------------------------------------------------------------------- #
# Main
# --------------------------------------------------------------------------- #
def build_parser():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("-s", "--samplename", required=True)
    p.add_argument("-i", "--indir", default=".",
                   help="directory to find input QC files in (default: cwd)")
    p.add_argument("--outdir", default=".",
                   help="directory to write the summary into (default: cwd)")
    p.add_argument("-o", "--outname", default=".qcsum.txt",
                   help="suffix appended to the sample name to form the output "
                        "file, i.e. <samplename><outname> (default: "
                        "<samplename>.qcsum.txt)")
    p.add_argument("-t", "--threads", type=int, default=4,
                   help="threads for samtools decompression (default: 4)")
    p.add_argument("--rrna_intervals",
                   default=os.path.expanduser(
                       "/hpf/projects/dig2/acelik_analysis/clinical_pipeline/genomes/hg38/"
                       "LL_rRNA.interval_list"),
                   help="Picard rRNA .interval_list for PCT_RIBOSOMAL "
                        "(default: %(default)s)")
    p.add_argument("--chry_genes",
                   default="/hpf/largeprojects/tgnode/sandbox/yliang_analysis/"
                           "pacbio_longread/annotation_files/"
                           "Homo_sapiens.GRCh38.115.chrYgenes.txt",
                   help="list of chrY Ensembl gene IDs (one per line), PAR "
                        "genes excluded, for CHRY_PERC / SEX_PRED "
                        "(default: %(default)s)")
    p.add_argument("--chry_sex_threshold", type=float, default=0.01,
                   help="CHRY_PERC%% above which chrY counts as male evidence "
                        "(default: %(default)s)")
    p.add_argument("--xist_threshold", type=float, default=10.0,
                   help="XIST TPM above which XIST counts as female evidence; "
                        "SEX_PRED is Unclear when the two signals disagree "
                        "(default: %(default)s)")
    # individual overrides; otherwise files are located under --indir by suffix
    p.add_argument("--bam_file")
    p.add_argument("--duplicate_file")
    p.add_argument("--star_file")
    p.add_argument("--qualimap_file")
    p.add_argument("--samtools_stats_file")
    p.add_argument("--tin_file")
    p.add_argument("--inner_distance_file")
    p.add_argument("--junction_file")
    p.add_argument("--infer_file")
    p.add_argument("--fastp_file")
    p.add_argument("--expr_file")
    p.add_argument("--ercc_file")
    p.add_argument("--sirv_file")
    p.add_argument("--tissue_pred_file")
    p.add_argument("--globin_file")
    return p


def _find(indir, pattern):
    """First path under indir matching a glob pattern (recursive via **)."""
    hits = sorted(glob.glob(os.path.join(indir, pattern), recursive=True))
    if not hits:
        return None
    if len(hits) > 1:
        warn("multiple matches for %r; using %s" % (pattern, hits[0]))
    return hits[0]


def discover_files(opt):
    """Locate inputs under --indir by file-name suffix (sample id is not
    required to match), filling in any options left unset."""
    d = opt.indir
    mqc = _find(d, "*multiqc_data") or os.path.join(d, "multiqc_data")
    found = {
        # top-level, suffix-matched
        "bam_file": _find(d, "*.mdup.bam"),
        "duplicate_file": _find(d, "*.duplicate.metrics"),
        "fastp_file": _find(d, "*fastp.json"),
        # salmon -> RSEM-format gene table (rule salmon2RSEM); the old RSEM
        # *.genes.results is no longer produced
        "expr_file": _find(d, "*_salmon_rsem.results"),
        "ercc_file": _find(d, "*.ercc.txt"),
        "sirv_file": _find(d, "*.sirv.txt"),
        "tissue_pred_file": _find(d, "*_tissuePred.txt"),
        "globin_file": _find(d, "*_globin_stats.txt"),
        # inside the <sample>multiqc_data folder (standard file names)
        "star_file": os.path.join(mqc, "multiqc_star.txt"),
        "samtools_stats_file": os.path.join(mqc, "multiqc_samtools_stats.txt"),
        "tin_file": os.path.join(mqc, "multiqc_rseqc_tin.txt"),
        "junction_file": os.path.join(mqc, "multiqc_rseqc_junction_annotation.txt"),
        "infer_file": os.path.join(mqc, "multiqc_rseqc_infer_experiment.txt"),
        # under the rustqc/ tree, located by suffix anywhere below --indir
        "qualimap_file": _find(d, "**/rnaseq_qc_results.txt"),
        "inner_distance_file": _find(d, "**/*.inner_distance_mean.txt"),
    }
    for k, v in found.items():
        if getattr(opt, k) is None:
            setattr(opt, k, v)


# Final column order of the one-row output table.
COLUMNS = [
    "SAMPLE_ID",
    "RAW_READS", "READS_AFTER_TRIM", "PCT_READS_PASS",
    "Q30_RATE_RAW", "Q30_RATE_TRIM", "MEAN_LEN_AFTER", "GC_PCT",
    "TOTAL_READS", "UNIQ_ALIGNED", "PCT_UNIQ_ALIGNED",
    "PCT_EXON", "PCT_INTRON", "PCT_INTERGENIC", "AMBIGUOUS_ALN",
    "MEDIAN_5PRIME_TO_3PRIME_BIAS", "PCT_RIBOSOMAL",
    "PERCENT_DUPLICATION", "INNER_DISTANCE_MEAN", "INNER_DISTANCE_MEDIAN",
    "TIN_MEDIAN", "TIN_MEAN", "STRANDEDNESS", "PCT_ANTISENSE", "MT_PCT",
    "TPM_0", "TPM_1", "TPM_10", "TPM_50",
    "CHRY_PERC", "XIST_TPM", "SEX_PRED",
    "NUM_JUNCS", "KNOWN_JUNCS", "NOVEL_JUNCS",
    "ERCC_TPM_SUM", "ERCC_TPM_0", "ERCC_TPM_1", "ERCC_R2",
    "SIRV_TPM_SUM", "SIRV1", "SIRV2", "SIRV3", "SIRV4", "SIRV5", "SIRV6", "SIRV7",
    # the sample-identity block: what came out of the library, what the classifier
    # made of it, what the sheet said it was, and whether anything else is in it
    "GLOBIN_PERC", "TISSUE_PRED",
    "TISSUE_SUBMITTED", "PURITY_CALL", "PURITY_FLAGS",
    "QC_CALL", "QC_FLAGS",
]


def main():
    opt = build_parser().parse_args()
    discover_files(opt)
    s = opt.samplename

    data = {"SAMPLE_ID": s}

    data.update(parse_fastp(opt.fastp_file))

    data.update(parse_star(opt.star_file, s))
    # TOTAL_READS (read pairs) is reused below for the MT_PCT denominator
    star_total = data["TOTAL_READS"] if isinstance(data["TOTAL_READS"], int) else None

    data.update(parse_qualimap(opt.qualimap_file))

    data.update(parse_duplicate_metrics(opt.duplicate_file))
    reads_mapped = get_reads_mapped(opt.samtools_stats_file, s)

    data.update(parse_tin(opt.tin_file, s))
    data.update(parse_inner_distance(opt.inner_distance_file, s))
    data.update(parse_junctions(opt.junction_file, s))
    data.update(parse_strandedness(opt.infer_file, s))

    data.update(calc_mt_pct(opt.bam_file, star_total))
    data.update(calc_ribosomal_pct(opt.bam_file, opt.rrna_intervals,
                                   reads_mapped, threads=opt.threads))

    data.update(parse_expression(opt.expr_file))
    data.update(parse_chry_sex(opt.expr_file, opt.chry_genes,
                               opt.chry_sex_threshold, opt.xist_threshold))
    data.update(parse_tissue_pred(opt.tissue_pred_file, s))
    data.update(parse_globin(opt.globin_file))
    data.update(parse_ercc(opt.ercc_file))
    data.update(parse_sirv(opt.sirv_file))

    # last: every metric it reads must already be in `data`
    data.update(qc_call(data))

    out_path = os.path.join(opt.outdir, s + opt.outname)
    with open(out_path, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(COLUMNS)
        w.writerow([data.get(c, NA) for c in COLUMNS])
    warn("wrote %s" % out_path)


if __name__ == "__main__":
    main()
