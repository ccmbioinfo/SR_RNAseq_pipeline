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


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("inputs", nargs="*", help="*.genes.results files")
    ap.add_argument("--model", default=os.path.join(os.path.dirname(__file__), "model.joblib"))
    ap.add_argument("--matrix", help="genes x samples TPM matrix (parquet/tsv)")
    ap.add_argument("--out", help="output TSV (default: stdout)")
    args = ap.parse_args()

    bundle = joblib.load(args.model)
    pipe = bundle["pipeline"]
    genes = bundle["genes"]
    thr = bundle["reject_threshold"]

    if args.matrix:
        tpm = load_matrix(args.matrix)
    elif args.inputs:
        tpm = pd.DataFrame({(s := tpm_from_genes_results(p)).name: s for p in args.inputs})
    else:
        sys.exit("provide genes.results files or --matrix")

    # align to the genes the model was trained on; missing -> 0
    tpm = tpm.reindex(genes).fillna(0.0)
    X = np.log2(tpm.T + 1.0)  # samples x genes

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

    txt = out.to_csv(sep="\t")
    if args.out:
        with open(args.out, "w") as f:
            f.write(txt)
        print(f"wrote {args.out} ({len(out)} samples)")
    else:
        sys.stdout.write(txt)


if __name__ == "__main__":
    main()
