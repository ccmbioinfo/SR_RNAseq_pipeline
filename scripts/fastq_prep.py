import pandas as pd
import sys, os, subprocess, logging, re
import argparse as arg
from snakemake.shell import shell
import yaml

## This script to do the following:

## if the input is paired-fastq --> create symlink in the output folder
## if the input are multiple paired-fastqs --> concatenate to form one fastq per end in the output folder
## if the input is a bam --> convert bam to paired end fastq in the output folder
## if the input is an ora --> convert ora to fastq in the output folder
## if the input are multiple oras --> convert ora to fastq and then concatenate in the output folder 

## the input for this script will be the sample file in the output folder


def create_symlink(src, dest):
    """symlink file"""
    if not os.path.isfile(dest):
        subprocess.check_call(["ln", "-s", src, dest])


def concatenate_fastq(r1, r2, output_dir, sample):
    """
    based on input file suffix, performs
    1. concatenation (more than 1 fastq per end)
    2. create symlink if only one fastq per end
    """
    if len(r1) > 1:  # multiple fastq files per end
        # check_fastq(r1, r2)
        r1_args = " ".join(r1)
        r2_args = " ".join(r2)
        cmd_r1 = [f"cat {r1_args} > {output_dir}/{sample}_R1.fastq.gz"]
        cmd_r2 = [f"cat {r2_args} > {output_dir}/{sample}_R2.fastq.gz"]
        subprocess.run(cmd_r1, shell=True)
        subprocess.run(cmd_r2, shell=True)

    elif ".fastq" in r1[0] or ".fq" in r1[0] and ".fastq" in r2[0] or ".fq" in r2[0]:
        # one fastq per end; symlink
        # check_fastq(r1, r2)
        dest = os.path.join(f"{output_dir}/{sample}_R1.fastq.gz")
        create_symlink(r1[0], dest)
        dest = os.path.join(f"{output_dir}/{sample}_R2.fastq.gz")
        create_symlink(r2[0], dest)
    else:
        print(r1, r2)
        print(f"Input {r1}, {r2} given for {sample} is not handled. Exiting!")
        exit()

def decompress_ora(r1,r2, output_path):
    decompressed_r1_list=[]
    decompressed_r2_list=[]
    for read in r1:
        command = (
            f"orad --ora-reference /opt/orad.2.7.0.linux/oradata "
            f"--path {output_path} "
            f"--threads 10 "
            f"{read}"
        )
        print(f"Command:{command}")
        run_orad=subprocess.run(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, shell=True)
        #decompressed_r1_list.append(f"{os.path.splitext(read)[0]}.gz")
        decompressed_r1_list.append(
            os.path.join(
                output_path,
                f"{os.path.splitext(os.path.basename(read))[0]}.gz"
            )
        )
    
    print(f"R1 List: {decompressed_r1_list}")

    for read in r2:
        command = (
            f"orad --ora-reference /opt/orad.2.7.0.linux/oradata "
            f"--path {output_path} "
            f"--threads 10 "
            f"{read}"
        )
        print(f"Command:{command}")
        run_orad=subprocess.run(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, shell=True)
        #decompressed_r2_list.append(f"{os.path.splitext(read)[0]}.gz")
        decompressed_r2_list.append(
            os.path.join(
                output_path,
                f"{os.path.splitext(os.path.basename(read))[0]}.gz"
            )
        )
    
    return decompressed_r1_list,decompressed_r2_list

def remove_bgzipped_files(decompressed_r1_list,decompressed_r2_list):
    for read in decompressed_r1_list+decompressed_r2_list:
        command=f"rm {read}"
        subprocess.run(command,shell=True)


if __name__ == "__main__":
    parser = arg.ArgumentParser(description='convert ora/bam to fastq if necessary')
    parser.add_argument('--sample_config', type=str, help='sample config file', action="store")
    parser.add_argument('--output', type=str, help='output directory', action="store")  
    args = parser.parse_args()
    
    with open(args.sample_config, "r") as f:
        config = yaml.safe_load(f)

    raw_input = config["sampleInfo"]["input"]
    samples_input = []

    raw_input = raw_input.strip()

    # BAM or single-end FASTQ
    if "," not in raw_input:
        samples_input = [raw_input]

    # paired-end FASTQ
    else:
        parts = raw_input.split(",")

        if len(parts) != 2:
            raise ValueError(f"Invalid paired input: {raw_input}")

        read1 = ";".join(p.strip() for p in parts[0].split(";") if p.strip())
        read2 = ";".join(p.strip() for p in parts[1].split(";") if p.strip())

        samples_input = [read1, read2]

    sampleID=config["sampleInfo"]['sampleID']
    print(samples_input[0])
    
    if samples_input[0].endswith(".bam"):
        print("File type is BAM")
        bam_file = samples_input[0]
        sort_cmd = "samtools sort -n -@ 10 -T {args.output}/{sampleID} {bam_file} |"
        bamtofastq_cmd = "samtools fastq -@ 8 -1 {args.output}/{sampleID}_R1.fastq.gz -2 {args.output}/{sampleID}_R2.fastq.gz -s /dev/null; "
        shell("(" + sort_cmd + bamtofastq_cmd + ")")
        ## conver bam

    else:
        r1, r2 = samples_input[0], samples_input[1]
        r1, r2 = r1.split(";"), r2.split(";")
        print(r1[0])
    
         # check if fastq files exist
        for read in r1 + r2:
            if not os.path.isfile(read):
                print(f"{read} does not exist. Exiting!")
                exit()
    
        #Checking if the fastq file extension is .gz or .ora
        if os.path.splitext(r1[0])[1] == ".gz":
            #If file is gzipped then concatenate the reads
            print(f"Files for {sampleID} are bgzipped.")
            concatenate_fastq(r1, r2, args.output, sampleID)
        elif os.path.splitext(r1[0])[1] == ".ora":
            #If file is ora compressed, then first decompress it
            print(f"Files for {sampleID} are in ora format. Decompressing the files using orad")
            decompressed_r1, decompressed_r2=decompress_ora(r1,r2,args.output)
            if len(decompressed_r1) == 1:
                ## move the decompressed r1 r2 to the sample folder
                r1_args = decompressed_r1[0]
                r2_args = decompressed_r2[0]
                cmd_r1 = [f"mv {r1_args} {args.output}/{sampleID}_R1.fastq.gz"]
                cmd_r2 = [f"mv {r2_args} {args.output}/{sampleID}_R2.fastq.gz"]
                subprocess.run(cmd_r1, shell=True)
                subprocess.run(cmd_r2, shell=True)
            else:
                ##more than 1 r1/r2
                concatenate_fastq(decompressed_r1, decompressed_r2, args.output, sampleID)
                #remove the bgzipped files
                remove_bgzipped_files(decompressed_r1,decompressed_r2)
