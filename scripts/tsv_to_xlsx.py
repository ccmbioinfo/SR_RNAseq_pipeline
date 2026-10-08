#!/usr/bin/env python3
"""tsv_to_xlsx.py — convert a tab-separated table to an .xlsx workbook with a
frozen header row and an autofilter on the header.

Used by parse_outrider.R / parse_fraser.R so the R side needs no Excel-writing R
package (openxlsx): the R scripts write a plain TSV and shell out to this to make
the .xlsx. Reproduces the old openxlsx output (firstRow=TRUE, withFilter=TRUE).

Usage:
    tsv_to_xlsx.py IN.tsv OUT.xlsx [SHEET_NAME]

Requires: pandas, and one of xlsxwriter / openpyxl (whichever is installed).
Missing values in the TSV are the literal string "NA" (as written by R's
write.table) or an empty field; both become blank cells.
"""
import csv
import sys


def main():
    if len(sys.argv) < 3:
        sys.exit("usage: tsv_to_xlsx.py IN.tsv OUT.xlsx [SHEET_NAME]")
    inp, outp = sys.argv[1], sys.argv[2]
    sheet = sys.argv[3] if len(sys.argv) > 3 else "results"

    try:
        import pandas as pd
    except ImportError:
        sys.exit("tsv_to_xlsx.py needs pandas (pip install pandas openpyxl)")

    # Read every column as text first (R wrote quote=FALSE, so disable quoting),
    # treating "NA"/empty as missing.
    df = pd.read_csv(inp, sep="\t", dtype=object, keep_default_na=False,
                     na_values=["NA", ""], quoting=csv.QUOTE_NONE)

    # Coerce a column to numeric only when every non-missing value is a number, so
    # Excel sorts/filters those columns numerically while ids/labels stay text.
    for c in df.columns:
        conv = pd.to_numeric(df[c], errors="coerce")
        nonmiss = df[c].notna()
        if nonmiss.any() and conv[nonmiss].notna().all():
            df[c] = conv

    # Pick whichever xlsx engine is available.
    engine = None
    for e in ("xlsxwriter", "openpyxl"):
        try:
            __import__(e)
            engine = e
            break
        except ImportError:
            continue
    if engine is None:
        sys.exit("tsv_to_xlsx.py needs xlsxwriter or openpyxl "
                 "(pip install xlsxwriter)")

    nrow, ncol = df.shape
    with pd.ExcelWriter(outp, engine=engine) as w:
        df.to_excel(w, index=False, sheet_name=sheet)
        ws = w.sheets[sheet]
        if engine == "xlsxwriter":
            ws.freeze_panes(1, 0)                       # freeze the header row
            if ncol:
                ws.autofilter(0, 0, nrow, ncol - 1)     # header + all data rows
        else:                                           # openpyxl
            ws.freeze_panes = "A2"
            ws.auto_filter.ref = ws.dimensions
    print("tsv_to_xlsx: wrote %s (%d rows x %d cols) via %s"
          % (outp, nrow, ncol, engine))


if __name__ == "__main__":
    main()
