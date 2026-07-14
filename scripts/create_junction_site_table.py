import os
import argparse
import pandas as pd


def main():
    parser = argparse.ArgumentParser(
        description="Convert a STAR SJ.out.tab file into junction and splice-site tables."
    )

    parser.add_argument(
        "--junction",
        required=True,
        help="Input STAR SJ.out.tab file"
    )

    parser.add_argument(
        "--output",
        required=True,
        help="Output directory"
    )

    args = parser.parse_args()

    # Create output directory if it doesn't exist
    os.makedirs(args.output, exist_ok=True)

    # Extract sample ID from filename
    # Example: Sample01SJ.out.tab -> Sample01
    basename = os.path.basename(args.junction)
    sampleID = basename.replace("SJ.out.tab", "")

    # Read STAR SJ.out.tab
    # Columns:
    # 0 chromosome
    # 1 intron start
    # 2 intron end
    # 3 strand
    # 4 intron motif
    # 5 annotation
    # 6 uniquely mapped reads
    # 7 multi-mapped reads
    # 8 maximum spliced alignment overhang
    sj = pd.read_csv(args.junction, sep="\t", header=None)

    # Convert STAR strand encoding
    strand_map = {
        0: "*",
        1: "+",
        2: "-"
    }

    ########################################################
    # Output 1: Junction table
    ########################################################

    junctions = pd.DataFrame({
        "seqnames": sj[0],
        "start": sj[1],
        "end": sj[2],
        "width": sj[2] - sj[1] + 1,
        "strand": sj[3].map(strand_map)
    })

    # Add unique read count as a column named after the sample
    junctions[sampleID] = sj[6]

    junction_file = os.path.join(
        args.output,
        f"{sampleID}_SJtable.tsv"
    )

    junctions.to_csv(
        junction_file,
        sep="\t",
        index=False
    )

    ########################################################
    # Output 2: Unique splice sites
    ########################################################

    donor = pd.DataFrame({
        "seqnames": sj[0],
        "start": sj[1] - 1,
        "end": sj[1],
        "width": 2,
        "strand": sj[3].map(strand_map),
        "type": "Donor"
    })

    acceptor = pd.DataFrame({
        "seqnames": sj[0],
        "start": sj[2],
        "end": sj[2] + 1,
        "width": 2,
        "strand": sj[3].map(strand_map),
        "type": "Acceptor"
    })

    splice_sites = pd.concat(
        [donor, acceptor],
        ignore_index=True
    )

    # Remove duplicated splice sites
    splice_sites = (
        splice_sites
        .drop_duplicates()
        .sort_values(["seqnames", "start", "type"])
        .reset_index(drop=True)
    )

    splice_file = os.path.join(
        args.output,
        f"{sampleID}_splice_sites.tsv"
    )

    splice_sites.to_csv(
        splice_file,
        sep="\t",
        index=False
    )

    print(f"Finished processing {sampleID}")
    print(f"Junction table: {junction_file}")
    print(f"Splice site table: {splice_file}")


if __name__ == "__main__":
    main()
