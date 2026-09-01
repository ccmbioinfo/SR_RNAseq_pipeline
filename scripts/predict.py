#!/usr/bin/env python
"""Predict tissue for new RNA-seq samples using the trained model.

Usage:
    # one or more RSEM genes.results files (sample name = parent folder, or file stem):
    python predict.py --model model.joblib SAMPLE1.genes.results SAMPLE2.genes.results ...

    # a precomputed TPM matrix (genes x samples; gene_id as first column/index):
    python predict.py --model model.joblib --matrix expr_tpm.parquet

Output: a TSV (to stdout, or --out FILE) with predicted tissue, per-class probabilities,
and an 'uncertain' flag when the top probability is below the model's reject threshold.
Genes the model expects but absent in a sample are filled with 0 (log2(0+1)=0).

Also emits a marker-panel purity score (see purity.py), which catches contamination
the classifier is structurally unable to see -- its probabilities saturate at 1.000
whether or not foreign RNA is present.

Purity is only scored when the tissue is known: the submitted tissue (from the
sample's config.yaml) and the predicted tissue must agree and be one of the four
classes. Otherwise purity_call is NA and purity_reason says why. To score a sample
whose tissue you know but the rule rejects -- an unsupported submitted label, a
prediction you disagree with, or a --matrix run with no config.yaml to read:

    python predict.py --score-as Fibroblasts SAMPLE.genes.results
    python predict.py --matrix expr_tpm.parquet --score-as-map rescore.tsv

where rescore.tsv is sampleID<TAB>tissue. An override is an assertion, not a hint:
it bypasses the agreement check, and is recorded as purity_ref_source=override so a
reviewer can see which calls a human asserted. Do NOT wire --score-as into the
Snakefile from config.yaml -- that would mark every sample as asserted and delete
the check entirely.

--globin-masked tells purity.py that the reference genome has the globin loci
masked out, in which case the globin panel reports nothing rather than a
reassuring zero.
"""
import argparse
import os
import sys
import warnings
warnings.filterwarnings("ignore")
import numpy as np
np.seterr(all="ignore")
import pandas as pd
import joblib

import purity


def tpm_from_genes_results(path):
    df = pd.read_csv(path, sep="\t", usecols=["gene_id", "TPM"])
    if df["gene_id"].duplicated().any():
        df = df.groupby("gene_id", as_index=False)["TPM"].sum()
    s = df.set_index("gene_id")["TPM"]
    # sample name: parent folder if informative, else file stem
    folder = os.path.basename(os.path.dirname(os.path.abspath(path)))
    stem = os.path.basename(path).replace(".genes.results", "").replace(".results", "")
    name = folder if folder and folder not in (".", "") else stem
    s.name = name if name else stem
    return s


def load_matrix(path):
    if path.endswith(".pkl"):
        m = pd.read_pickle(path)
    elif path.endswith(".parquet"):
        m = pd.read_parquet(path)  # requires pyarrow (optional)
    else:
        m = pd.read_csv(path, sep="\t", index_col=0)
    return m  # genes x samples


def load_score_as_map(path):
    """{sampleID: tissue} from a two-column TSV."""
    out = {}
    if not path:
        return out
    with open(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) >= 2 and f[0] and f[0] != "sampleID":
                out[f[0]] = f[1].strip()
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("inputs", nargs="*", help="*.genes.results files")
    ap.add_argument("--model", default=os.path.join(os.path.dirname(__file__), "model.joblib"))
    ap.add_argument("--matrix", help="genes x samples TPM matrix (parquet/tsv)")
    ap.add_argument("--out", help="output TSV (default: stdout)")
    ap.add_argument("--score-as", dest="score_as",
                    help="assert the tissue to score purity against, for every "
                         "sample in this run (Blood|LCL|Fibroblasts|Muscle)")
    ap.add_argument("--score-as-map", dest="score_as_map",
                    help="TSV of sampleID<TAB>tissue, for batch reruns")
    ap.add_argument("--globin-masked", dest="globin_masked", action="store_true",
                    help="reference genome has the globin loci masked; the globin "
                         "panel reports NA instead of a meaningless 0")
    args = ap.parse_args()

    bundle = joblib.load(args.model)
    pipe = bundle["pipeline"]
    genes = bundle["genes"]
    thr = bundle["reject_threshold"]

    sample_paths = {}
    if args.matrix:
        tpm = load_matrix(args.matrix)
    elif args.inputs:
        cols = {}
        for path in args.inputs:
            col = tpm_from_genes_results(path)
            cols[col.name] = col
            sample_paths[col.name] = path   # purity reads config.yaml beside it
        tpm = pd.DataFrame(cols)
    else:
        sys.exit("provide genes.results files or --matrix")

    # Purity is scored on the FULL matrix, before the reindex below: the marker
    # genes are deliberately not among the ones the classifier selects, so they
    # would all be dropped if this ran on the model's gene subset.
    panels = purity.load_panels()
    limits = purity.load_limits()
    scores = {c: purity.panel_scores(tpm[c], panels) for c in tpm.columns}
    submitted = {c: purity.submitted_tissue(sample_paths[c]) if c in sample_paths else None
                 for c in tpm.columns}
    overrides = load_score_as_map(args.score_as_map)

    # align to the genes the model was trained on; missing -> 0
    model_tpm = tpm.reindex(genes).fillna(0.0)
    X = np.log2(model_tpm.T + 1.0)  # samples x genes

    proba = pipe.predict_proba(X.values)
    classes = list(pipe.classes_)
    pred = np.array(classes)[proba.argmax(axis=1)]
    maxp = proba.max(axis=1)
    final = np.where(maxp < thr, "uncertain/other", pred)

    out = pd.DataFrame(proba, columns=[f"p_{c}" for c in classes], index=X.index)
    out.insert(0, "max_prob", maxp.round(3))
    out.insert(0, "predicted", final)
    out.insert(1, "top_class", pred)
    out = out.round(3)
    out.index.name = "sampleID"

    # purity columns, appended after the classifier's own
    calls = {}
    for i, s in enumerate(X.index):
        calls[s] = purity.purity_call(
            scores[s],
            submitted=submitted.get(s),
            predicted=final[i],
            override=args.score_as or overrides.get(s),
            limits=limits,
            globin_masked=args.globin_masked,
            tpm=tpm[s])
    out["submitted_tissue"] = [submitted.get(s) or "NA" for s in X.index]
    out["purity_call"] = [calls[s]["call"] for s in X.index]
    out["purity_reason"] = [calls[s]["reason"] or "NA" for s in X.index]
    out["purity_reference"] = [calls[s]["reference"] or "NA" for s in X.index]
    out["purity_ref_source"] = [calls[s]["ref_source"] or "NA" for s in X.index]
    out["purity_flags"] = [calls[s]["flags"] for s in X.index]
    out["skin_form"] = [calls[s]["skin_form"] or "NA" for s in X.index]
    out["purity_not_assessed"] = [calls[s]["not_assessed"] or "NA" for s in X.index]
    for panel in sorted(panels):
        out[f"panel_{panel}"] = [scores[s][panel].score for s in X.index]
        out[f"support_{panel}"] = [scores[s][panel].support for s in X.index]

    txt = out.to_csv(sep="\t")
    if args.out:
        with open(args.out, "w") as f:
            f.write(txt)
        print(f"wrote {args.out} ({len(out)} samples)")
    else:
        sys.stdout.write(txt)


if __name__ == "__main__":
    main()
