#!/usr/bin/env python
"""Marker-panel purity scoring for a single RNA-seq library.

The tissue classifier answers "which of the four tissues is this closest to". It
cannot answer "is anything else mixed in": its probabilities saturate at 1.000, and
SelectKBest keeps only genes that separate the four classes, so every keratin --
near zero in all four -- is dropped before the model ever sees it. A library can be
10% someone else's skin and still be called LCL at p = 1.000.

This module measures that independently, by summing fixed marker panels (see
purity_panels.tsv) and comparing each against an absolute TPM threshold for the
reference tissue. Thresholds, denominators and detection limits live in
purity_limits.tsv, NOT in this file: recalibrating is a data change, not a code
change. Regenerate that table with calibrate_purity.py.

Three rules do most of the work:

  Reference tissue    A sample is only scored when the submitted and predicted
                      tissue agree and are one of the four known classes, or when
                      an operator asserts one with --score-as. Anything else is NA
                      with a reason. Forcing an unsupported tissue (skin, colon,
                      cell_line) through a four-class model and then scoring it
                      against that class's thresholds produced the two largest
                      false signals in the 2026-08 cohort.

  Support             A panel only flags when >= MIN_SUPPORT of its genes clear
                      both an absolute floor and a share of the panel's top gene.
                      A bare sum is owned by its loudest member: four Blood
                      libraries scored 150-308 on skin from KRT1 alone, with KRT10
                      near zero, where real epidermis runs those two about 1:1.

  Expected panels     Panels native to a tissue are not scored against it: muscle
                      legitimately contains connective tissue, blood contains every
                      leukocyte lineage. Encoded as status=expected rows in
                      purity_limits.tsv.

Two skin panels rather than one, because they call for different action:
    cornified >> basal   mature skin tissue reached the tube -- a handling error
    basal alone          keratinocytes growing in the dish -- an impure cell line
Reported as skin_form, since a single "skin score" would merge the two and send
you after the wrong problem. Cultured keratinocytes also invert the basal keratins
(KRT17 > KRT14/KRT5, the activated set) where epidermis does the opposite, so that
ratio is used as a second, independent axis.

Public API:
    load_panels(path=None)                  -> {panel: [gene_id, ...]}
    load_limits(path=None)                  -> {(tissue, panel): Limit}
    panel_scores(tpm, panels)               -> {panel: PanelScore}
    canonical_tissue(name)                  -> str or None
    reference_tissue(submitted, predicted, override=None)
                                            -> (tissue|None, source, reason)
    submitted_tissue(genes_results_path)    -> str or None
    purity_call(...)                        -> dict
"""
import collections
import os
import re

HERE = os.path.dirname(os.path.abspath(__file__))
PANEL_FILE = os.path.join(HERE, "purity_panels.tsv")
LIMIT_FILE = os.path.join(HERE, "purity_limits.tsv")

# ---- the four classes the model can represent -------------------------------
# Submitted tissue is free text typed by whoever filled in the sample sheet. Only
# these spellings map into a class; everything else (Skin, SIGMOIDCOLON, Heart,
# cell_line, ...) is deliberately NOT guessed at -- see reference_tissue().
TISSUE_ALIASES = {
    "blood": "Blood",
    "lcl": "LCL",
    "fibroblast": "Fibroblasts",
    "fibroblasts": "Fibroblasts",
    "muscle": "Muscle",
}

# Submitted labels that name no tissue at all, and so are treated as "unknown"
# rather than "unsupported": the prediction is used as the reference instead.
#
# "cell_line" says how the material was grown, not what it is -- GM17398 and
# GM18429 are fibroblast lines, GM27955 is an LCL. That makes it different in kind
# from Skin or SIGMOIDCOLON, which are real tissues the four-class model has no
# label for and must therefore force into the wrong class. A cell line usually IS
# one of the four, so the prediction is the best evidence available and the
# agreement check has nothing to compare it against.
#
# The prediction still has to be usable: below the model's reject threshold it
# arrives as "uncertain/other" and the sample goes back to being unscored.
UNINFORMATIVE_SUBMITTED = {"cell_line", "cellline", "cell line"}

# ---- support rule -----------------------------------------------------------
# A panel gene counts toward support when it clears BOTH an absolute floor and a
# share of the panel's own top gene. The relative arm is what suppresses the
# KRT1-only artifact: at 5%, KRT10=7 against KRT1=172 does not count, while real
# epidermis (KRT1 and KRT10 within ~1.2x) keeps every member.
MIN_SUPPORT = 2
SUPPORT_FLOOR_TPM = 1.0
SUPPORT_FLOOR_SHARE = 0.05

# Panel score below which nothing is actionable anyway: MIN_SUPPORT genes sitting
# exactly on the absolute floor. Thresholds are never calibrated below this.
MIN_REVIEW = 10.0

# blood_globin is only meaningful where bbduk did NOT run. The Snakefile gates that
# on the SUBMITTED tissue, so this is checked against config.yaml, not against the
# prediction: a sample submitted as Blood has had its globin stripped whatever it
# is predicted to be.
GLOBIN_DEPLETED = {"blood"}

# ...and only where the reference genome still carries the globin loci. config.yaml
# points at a globin-MASKED build; config_unmask.yaml does not. On a masked run
# every globin TPM is ~0 whatever is in the tube, so the panel must report NA
# rather than a reassuring 0. predict.py passes this through from the Snakefile.
GLOBIN_PANEL = "blood_globin"

SKIN_FORM_RATIO = 0.5      # cornified/basal above this reads as tissue
SKIN_ACTIVATED_RATIO = 1.0 # KRT17/(KRT5+KRT14) above this reads as culture
# Ensembl 115 IDs, pinned like the panel file.
KRT17, KRT5, KRT14 = "ENSG00000128422", "ENSG00000186081", "ENSG00000186847"

_SAMPLEINFO = re.compile(r"^sampleInfo:\s*$")
_KEY = re.compile(r"^\s{2}(\w+):\s*(.*?)\s*$")

PanelScore = collections.namedtuple("PanelScore", "score support n_genes top_gene")
Limit = collections.namedtuple(
    "Limit", "tissue panel review fail pure_level detection_limit status")


def load_panels(path=None):
    """{panel: [gene_id, ...]} from the pinned TSV."""
    panels = collections.OrderedDict()
    with open(path or PANEL_FILE) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3 or parts[0] == "gene_id":
                continue
            panels.setdefault(parts[2], []).append(parts[0])
    return panels


def _num(s):
    try:
        return float(s)
    except (TypeError, ValueError):
        return None


def load_limits(path=None):
    """{(tissue, panel): Limit} from purity_limits.tsv.

    status is one of:
        active        scored against (review, fail)
        underpowered  scored, but detection_limit is worse than the 2% target --
                      a PASS means "not assessed below detection_limit", not clean
        expected      native to the tissue, never scored against it
        no_reference  no pure-tissue denominator, so no detection limit is known
    """
    limits = {}
    with open(path or LIMIT_FILE) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 7 or f[0] == "reference_tissue":
                continue
            limits[(f[0], f[1])] = Limit(f[0], f[1], _num(f[2]), _num(f[3]),
                                         _num(f[4]), _num(f[5]), f[6])
    return limits


def panel_scores(tpm, panels, human_total=None):
    """{panel: PanelScore} -- summed TPM plus the support count.

    `tpm` is any gene_id -> TPM mapping (a dict or a pandas Series). It is
    renormalised to human genes first: RSEM's TPM column includes the ERCC and
    SIRV spike-ins, which run ~0.03% of a good library but reached 7% in the worst
    one seen, and leaving them in would make a heavily spiked library read as
    cleaner than it is.

    The score stays a plain sum -- that is what makes it convertible to a fraction
    of the transcriptome in purity_call(). Robustness comes from `support`, not
    from changing the aggregate.

    `human_total` short-circuits the full pass over the matrix when the caller has
    already summed the human genes (calibrate_purity.py does, vectorised). Leave it
    None in the pipeline, where one sample is read at a time.
    """
    if human_total is None:
        human = {}
        total = 0.0
        for gene, value in tpm.items():
            gene = str(gene)
            if not gene.startswith("ENSG"):
                continue
            try:
                value = float(value)
            except (TypeError, ValueError):
                continue
            human[gene] = value
            total += value

        def lookup(g):
            return human.get(g, 0.0)
    else:
        total = float(human_total)

        def lookup(g):
            try:
                return float(tpm.get(g, 0.0) or 0.0)
            except (TypeError, ValueError):
                return 0.0

    scale = (1e6 / total) if total > 0 else 0.0

    out = collections.OrderedDict()
    for panel, genes in panels.items():
        values = [lookup(g) * scale for g in genes]
        total_tpm = sum(values)
        top = max(values) if values else 0.0
        cutoff = max(SUPPORT_FLOOR_TPM, SUPPORT_FLOOR_SHARE * top)
        support = sum(1 for v in values if v >= cutoff)
        out[panel] = PanelScore(round(total_tpm, 2), support, len(genes),
                                round(top, 2))
    return out


def canonical_tissue(name):
    """One of the four class names, or None when the label is out of vocabulary."""
    if not name:
        return None
    return TISSUE_ALIASES.get(str(name).strip().strip("'\"").lower())


def reference_tissue(submitted, predicted, override=None):
    """(tissue, source, reason) -- which threshold table applies, and why.

    Scoring a sample against the wrong tissue's thresholds is worse than not
    scoring it: the two loudest signals in the 2026-08 cohort were a skin biopsy
    and a sigmoid colon biopsy, both forced into Fibroblasts by a four-class model
    and then scored as though they were fibroblast cultures.

    An override is an assertion by an operator, not a hint -- it suppresses the
    mismatch check, which is exactly its purpose, and is recorded as such so a
    reviewer can tell which calls were human-asserted.

    ref_source records which of the three routes was taken -- "agreed",
    "predicted" (an uninformative submitted label, see UNINFORMATIVE_SUBMITTED) or
    "override" -- so a reviewer can tell how much the reference is worth. Only
    "agreed" had two independent sources to check against each other.
    """
    if override:
        tissue = canonical_tissue(override)
        if tissue is None:
            return None, "", "unsupported_override"
        return tissue, "override", ""
    if submitted is None:
        return None, "", "no_submitted"
    if str(submitted).strip().strip("'\"").lower() in UNINFORMATIVE_SUBMITTED:
        # No tissue was submitted to agree with, so there is nothing to check the
        # prediction against. Scoring on the prediction alone is weaker evidence
        # than "agreed" and is labelled as such.
        pred = canonical_tissue(predicted)
        if pred is None:
            return None, "", "cell_line_unpredicted"
        return pred, "predicted", ""
    sub = canonical_tissue(submitted)
    if sub is None:
        return None, "", "unsupported_submitted"
    pred = canonical_tissue(predicted)
    if pred is None or sub != pred:
        return None, "", "tissue_mismatch"
    return sub, "agreed", ""


def submitted_tissue(genes_results_path):
    """Tissue as submitted, from the config.yaml the pipeline writes beside the RSEM
    output (rule rsem writes into config["sampleInfo"]["dir"], the same directory).
    None when the file or the field is absent.

    Hand-parsed rather than via PyYAML, which the pipeline container is not
    guaranteed to carry. Anchored to the sampleInfo: block because resources: uses
    the same two-space indent.
    """
    cfg = os.path.join(os.path.dirname(os.path.abspath(genes_results_path)), "config.yaml")
    if not os.path.exists(cfg):
        return None
    inblock = False
    try:
        with open(cfg) as fh:
            for line in fh:
                if _SAMPLEINFO.match(line):
                    inblock = True
                    continue
                if inblock:
                    if line[:1] not in (" ", "\t", "\n"):   # dedent -> block ended
                        break
                    m = _KEY.match(line.rstrip("\n"))
                    if m and m.group(1) == "tissue":
                        return m.group(2).strip("'\"").strip() or None
    except OSError:
        return None
    return None


def skin_form(scores, tpm=None):
    """'tissue' or 'culture' -- which kind of skin signal is present.

    Only meaningful for a sample that actually tripped a skin panel, so callers go
    through purity_call(), which returns None when neither did. A clean culture has
    a little basal keratin and no story to tell about it.

    Two independent axes, and they agreed on all three known cases in the 2026-08
    cohort. The cornified/basal ratio is the primary one. The keratin program is
    the second: epidermis runs KRT5/KRT14 far above KRT17 (3991/4172 vs 452 in the
    one skin library available), while keratinocytes in a flask invert it -- they
    express the activated set KRT6/16/17 and never cornify. The second axis only
    votes when the first is ambiguous, since it needs the basal panel to be lit.
    """
    cornified = scores.get("skin_cornified")
    basal = scores.get("skin_basal")
    cornified = cornified.score if cornified is not None else 0.0
    basal = basal.score if basal is not None else 0.0
    if cornified / max(basal, 1.0) > SKIN_FORM_RATIO:
        return "tissue"
    if tpm is not None:
        k17 = float(tpm.get(KRT17, 0.0) or 0.0)
        k5_14 = float(tpm.get(KRT5, 0.0) or 0.0) + float(tpm.get(KRT14, 0.0) or 0.0)
        if k5_14 > 0 and k17 / k5_14 > SKIN_ACTIVATED_RATIO:
            return "culture"
    return "culture"


def underpowered_panels(reference, limits):
    """Panels whose detection limit is worse than the 2% target for this tissue.

    A PASS on one of these means "not assessed below detection_limit", which is a
    different statement from "clean" and should not be read as one.
    """
    return sorted(p for (t, p), lim in limits.items()
                  if t == reference and lim.status == "underpowered")


def purity_call(scores, submitted=None, predicted=None, override=None,
                limits=None, globin_masked=False, tpm=None):
    """A dict describing one sample's purity.

    Keys: call, reason, reference, ref_source, flags, skin_form, not_assessed.

    call    PASS / REVIEW / FAIL, or NA when no reference tissue could be fixed.
    flags   'panel=value(fraction)' for each panel that tripped, worst first,
            ';'-joined. Ordering is by how far past its own review threshold a
            panel sits, not by raw TPM -- different panels are not on a common
            scale (muscle spends 4.9% of its transcriptome on its markers, T cells
            0.09%, so equal TPM means very different amounts of contaminant).
    """
    limits = load_limits() if limits is None else limits
    reference, source, reason = reference_tissue(submitted, predicted, override)
    if reference is None:
        return {"call": "NA", "reason": reason, "reference": "", "ref_source": "",
                "flags": "", "skin_form": None, "not_assessed": ""}

    depleted = (submitted or "").strip().strip("'\"").lower() in GLOBIN_DEPLETED

    worst, tripped = 0, []
    for panel, ps in scores.items():
        lim = limits.get((reference, panel))
        if lim is None or lim.status in ("expected", "skip"):
            continue
        if lim.review is None or lim.fail is None:
            continue
        if panel == GLOBIN_PANEL and (depleted or globin_masked or submitted is None):
            # Depleted upstream, or masked out of the reference -- either way the
            # number carries no information about what was in the tube.
            continue
        if ps.support < MIN_SUPPORT:
            # One loud gene is not a lineage. Score still reported, cannot flag.
            continue
        level = 2 if ps.score > lim.fail else (1 if ps.score > lim.review else 0)
        if level:
            ratio = ps.score / lim.review if lim.review > 0 else ps.score
            tripped.append((level, ratio, ps.score, panel, lim))
            worst = max(worst, level)

    tripped.sort(key=lambda t: (-t[0], -t[1]))
    parts = []
    for _, _, value, panel, lim in tripped:
        if lim.pure_level:
            parts.append("%s=%g(~%.2g%%)" % (panel, value, 100.0 * value / lim.pure_level))
        else:
            parts.append("%s=%g" % (panel, value))
    flags = ";".join(parts)

    form = (skin_form(scores, tpm)
            if any(t[3].startswith("skin_") for t in tripped) else None)
    return {"call": ("PASS", "REVIEW", "FAIL")[worst],
            "reason": "",
            "reference": reference,
            "ref_source": source,
            "flags": flags,
            "skin_form": form,
            "not_assessed": ",".join(underpowered_panels(reference, limits))}
