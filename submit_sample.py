#!/usr/bin/env python3

import os
import yaml
import pandas as pd
import argparse
import subprocess


def create_sample_configs(sample_sheet, output_dir, config_yaml, pipeline, resubmit=False):
    # Load template config once
    with open(config_yaml) as f:
        template_config = yaml.safe_load(f)

    # Read TSV
    df = pd.read_csv(
        sample_sheet,
        sep="\t",
        header=None,
        names=["sampleId", "read", "tissue", "vcf", "hpoterm"],
        dtype=str,
    ).fillna("")

    for _, row in df.iterrows():
        sample_id = row["sampleId"]
        sample_dir = os.path.join(output_dir, sample_id)
        os.makedirs(sample_dir, exist_ok=True)

        sample_config_path = os.path.join(sample_dir, "config.yaml")

        #
        # ---------------------------
        # Create config only if not resubmitting
        # ---------------------------
        #
        if not resubmit:

            # Deep copy template config
            config = yaml.safe_load(yaml.safe_dump(template_config))

        #
        # ---------------------------
        # Parse INPUT (BAM / SE / PE)
        # ---------------------------
        #
            raw_input = row["read"].strip()

            if not raw_input:
                raise ValueError(f"Sample {sample_id}: empty input field")

            input_value = ""

            # CASE 1: BAM
            if raw_input.lower().endswith(".bam"):
                input_value = raw_input

            else:
                parts = raw_input.split(",")

                # CASE 2: single-end FASTQ
                if len(parts) == 1:
                    input_value = parts[0].strip()

                # CASE 3: paired-end FASTQ
                elif len(parts) == 2:
                    read1 = ";".join(
                        f.strip() for f in parts[0].split(";") if f.strip()
                    )
                    read2 = ";".join(
                        f.strip() for f in parts[1].split(";") if f.strip()
                    )

                    input_value = ",".join([read1, read2])

                else:
                    raise ValueError(
                        f"Sample {sample_id}: invalid input '{raw_input}'. "
                        "Expected BAM, single-end FASTQ, or paired-end FASTQ"
                    )

            #
            # ---------------------------
            # Parse HPO file
            # ---------------------------
            #
            hpo_terms = []

            if row["hpoterm"]:
                if not os.path.exists(row["hpoterm"]):
                    raise FileNotFoundError(
                        f"Sample {sample_id}: HPO file not found: {row['hpoterm']}"
                    )

                with open(row["hpoterm"]) as f:
                   hpo_terms = [
                       line.strip()
                       for line in f
                       if line.strip()
                   ]

           #
           # ---------------------------
           # Update config
           # ---------------------------
           #
            config["sampleInfo"] = {
                "sampleID": sample_id,
                "input": input_value,
                "tissue": row["tissue"],
                "vcf": row["vcf"],
                "hpoterm": ",".join(hpo_terms),
                "dir": sample_dir,
            }

        #
        # ---------------------------
        # Write config
        # ---------------------------
        #
            with open(sample_config_path, "w") as out:
               yaml.safe_dump(
                   config,
                   out,
                   default_flow_style=False
               )

            print(f"Created {sample_config_path}")
        
        else:
            if not os.path.exists(sample_config_path):
                raise FileNotFoundError(
                    f"Resubmit requested but {sample_config_path} does not exist."
                )

            print(f"Reusing existing {sample_config_path}")

        sbatch_cmd = "sbatch --export pipeline_location={pipe},configfile={config} -J {name} -e {err} -o {out} {pipe}/submit.sh " \
        .format(**{'pipe': pipeline,
                   'config': sample_config_path,
                   "name": sample_id,
                   'dir': sample_dir,
                   'err': sample_dir + "/" + sample_id + "-%j.err",
                   'out': sample_dir + "/" + sample_id + "-%j.out"})
        
        try:
            jobid = subprocess.run(
            sbatch_cmd,
            shell=True,
            cwd=sample_dir,   # run sbatch from sample_dir
            stdout=subprocess.PIPE,
            check=True
        )
        except subprocess.CalledProcessError as err:
            raise err


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Parse sample table file(s) and submit jobs")
    parser.add_argument('-p', '--pipeline_location', required=True,
                        default="/hpf/largeprojects/ccmbio/yliang/clinical_pipeline",
                        help="Folder with snakefile and other scripts")
    parser.add_argument('-c', '--configfile', default=None,
                        help="config.yaml for the pipeline if None (default) will look in pipeline_location for a config.yaml")
    parser.add_argument('-o', '--output_directory', required=True, help="Where the results will end up")
    parser.add_argument('--resubmit', action='store_true', default=False, help="resubmit these samples")
    parser.add_argument("-s", "--samplefile", required=True,
                       help='samplefile, one line per sample col1 samplename col2 reads comma separated col3 tissue quoted with double (") quotes' )
    
    args = parser.parse_args()

    pipeline = args.pipeline_location
    resubmit = args.resubmit
    config = args.configfile

    config_yaml = (
        os.path.abspath(os.path.join(pipeline, "config.yaml"))
        if config is None
        else os.path.abspath(config)
    )

    sample_sheet = args.samplefile
    main_output_dir = args.output_directory

    create_sample_configs(
        sample_sheet=sample_sheet,
        output_dir=main_output_dir,
        config_yaml=config_yaml,
        pipeline=pipeline,
        resubmit=resubmit
    )
