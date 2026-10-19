import argparse as arg
import pandas as pd

def percent(values):
    """
    get a list of values return a list of percentages
    :param values: list of values whatever it may be
    :return: list of percentages
    """
    total=sum(values)
    percentages=[]
    for value in values:
        if total==0:
            percentages.append(0)
        else:
            pct=value/total
            percentages.append(pct)
    return percentages


if __name__ == "__main__":
    parser = arg.ArgumentParser(description='combine detected ERCC/SIRV along with actual values for ERCC just a '
                                            'comparison of actual vs detected concentraion but in different units'
                                            'using TPM for all these measurements')
    parser.add_argument('--gene', type=str, help='unused; kept for the existing Snakefile call', action="store")
    parser.add_argument('--iso', type=str, help='isoform expression values', action="store")
    parser.add_argument('--actual', type=str, help='actual concentrations', action="store")
    parser.add_argument('--samplename', type=str, help='samplename', action="store")
    parser.add_argument('--analysis_path', type=str, help='analysis path', action="store")

    args = parser.parse_args()

    # Both ERCC and SIRV come from the transcript-level salmon quant.sf (--iso).
    # ERCC used to be read from the gene table (--gene), but salmon2RSEM.R keeps
    # only genes in its --genelist, which drops every ERCC- row, so ercc.txt came
    # out header-only. Each ERCC is a single transcript mapped to itself in
    # tx2gene, so its transcript TPM is its gene TPM. --gene is still accepted so
    # the Snakefile call is unchanged, but it is no longer read.
    quant=pd.read_csv(args.iso, header=0, sep="\t")
    iso=quant[quant["Name"].str.contains("SIRV")]
    ercc=quant[quant["Name"].str.startswith("ERCC-")]
    ercc=ercc[["Name", "TPM"]].rename(columns={"Name": "gene_id"})

    actual = pd.read_csv(args.actual, header=0, sep=",")

    # the concentration file keys ERCC by gene_id (ERCC-00130); its transcript_id
    # column holds the GenBank accession
    gene_merged = pd.merge(ercc, actual, how='inner', on='gene_id')

    gene_name=args.analysis_path+"/"+args.samplename+".ercc.txt"
    gene_merged.to_csv(gene_name, sep="\t", index=False, header=True)

    iso_name = args.analysis_path+"/"+args.samplename + ".sirv.txt"
    #iso_merged = pd.merge(iso[["transcript_id", "TPM", "IsoPct"]],
    #                       actual, how='inner', on='transcript_id')
    iso_merged = pd.merge(
    iso[["Name", "TPM", "NumReads"]],
    actual,
    how="inner",
    left_on="Name",
    right_on="transcript_id"
)

    iso_merged.to_csv(iso_name, sep="\t", index=False, header=True)
