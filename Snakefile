import os
import pandas as pd
import math
import re

tissue=config["sampleInfo"]["tissue"]
sampleID=config["sampleInfo"]["sampleID"]
output_dir=config["sampleInfo"]["dir"]
sample_config=os.path.join(config["sampleInfo"]["dir"], "config.yaml")

trimming=True
SCRATCH = os.environ["SCRATCH"]

# Whether the alignment reference has the globin loci masked out. Detected from
# the fasta name rather than configured separately, so the two cannot drift apart.
globin_masked = "maskedglobin" in os.path.basename(config["resources"]["fasta"]).lower()

if tissue.lower() == "blood":
    filter_cmd = """
    bbduk.sh -Xmx24g threads=12 \
    in1={input.r1} in2={input.r2} \
    out1={output.out1} out2={output.out2} \
    ref={input.ref} outm={params.globinread} \
    k=31 stats={output.stats} 2> {params.bbduk_log} || {{ cat {params.bbduk_log} >&2; exit 1; }}
    cat {params.bbduk_log} >&2

    gzip -f {params.globinread}

    # Scratch writes are only acknowledged by the page cache; a failed writeback
    # is reported to nobody once bbduk exits. Force the flush so the error lands
    # here instead of surfacing as corrupt reads in salmon/STAR.
    sync --data {output.out1} {output.out2}

    # bbduk's "Result:" counts reads across BOTH mates. Do not derive this from
    # stats' #Matched: in paired mode both mates are dropped when either matches,
    # so #Matched undercounts removals and the check would misfire every run.
    n_result=$(awk '/^Result:/{{print $2}}' {params.bbduk_log})
    exp=$(( n_result / 2 ))
    n1=$(( $(wc -l < {output.out1}) / 4 ))
    n2=$(( $(wc -l < {output.out2}) / 4 ))
    if [ "$n1" -ne "$exp" ] || [ "$n2" -ne "$exp" ]; then
        echo "ERROR: bbduk output truncated: expected $exp records/mate, got R1=$n1 R2=$n2" >&2
        exit 1
    fi
    """
else:
    filter_cmd = """
    mv {input.r1} {output.out1}
    mv {input.r2} {output.out2}

    echo "NA" > {output.stats}
    echo "NA" > {params.globinread}
    """


#include: "rules/alignment_expr_qc.smk"

rule all:
    input:
        os.path.join(output_dir,f'{sampleID}.qcsum.txt')

if config["sampleInfo"]["vcf"]:
    rule format_vcf:
            input:
                vcf=config["sampleInfo"]["vcf"]
            output:
                het_only_vcf_gz=temp(os.path.join(output_dir,f'{sampleID}.hetsnp.vcf.gz')),
                STAR_vcf=temp(os.path.join(output_dir,f'{sampleID}_STAR.vcf')),
                dnacount_txt=os.path.join(output_dir,"ASE_output",f'{sampleID}DNA_alleleCount.txt')
            params:
                outdir=output_dir,
                sampleID=f'{sampleID}'
            run:
                shell('module load vt/fb0288b; module load tabix; module load bcftools/1.20; scripts/vcf_formatparser.sh {input.vcf} {output.het_only_vcf_gz} {output.dnacount_txt} {output.STAR_vcf}')

rule fastq_prep:
    input: sample_config
    output:
        fastq_r1=temp(os.path.join(output_dir,f'{sampleID}_R1.fastq.gz')),
        fastq_r2=temp(os.path.join(output_dir,f'{sampleID}_R2.fastq.gz'))
    params:
        outdir=output_dir,
        sampleID=f'{sampleID}',
        script=os.path.join(workflow.basedir, "scripts/fastq_prep.py")
    singularity:
        "/hpf/largeprojects/tgnode/sandbox/yliang_analysis/tools/containers/sr_rnaseq.sif"
    shell:
        """
        python {params.script} --sample_config {input} --output {params.outdir}
        """
 
rule trimming:
    input:
        r1=rules.fastq_prep.output.fastq_r1,
        r2=rules.fastq_prep.output.fastq_r2
    output:
        paired1=temp(os.path.join(os.environ["SCRATCH"],f'{sampleID}_R1_trim_paired.fq')),
        paired2=temp(os.path.join(os.environ["SCRATCH"],f'{sampleID}_R2_trim_paired.fq')),
        #paired1=temp(os.path.join(output_dir,f'{sampleID}_R1_trim_paired.fq')),
        #paired2=temp(os.path.join(output_dir,f'{sampleID}_R2_trim_paired.fq')),
        unpaired1=os.path.join(output_dir,f'{sampleID}_unpairR1.fq.gz'),
        unpaired2=os.path.join(output_dir,f'{sampleID}_unpairR2.fq.gz'),
        failed_read=os.path.join(output_dir,f'{sampleID}_failed.fq.gz'),
        json=os.path.join(output_dir,f'{sampleID}fastp.json'),
        html=os.path.join(output_dir,f'{sampleID}fastp.html')
    params:
        outdir=output_dir,
        sampleID=f'{sampleID}',
        scratch=SCRATCH
    threads: 10
    singularity:
        "/hpf/largeprojects/tgnode/sandbox/yliang_analysis/tools/containers/sr_rnaseq.sif"
    shell:
        """
        fastp -w 10 -i {input.r1} -I {input.r2} -o {output.paired1} -O {output.paired2} \
        --unpaired1 {output.unpaired1} --unpaired2 {output.unpaired2} --failed_out {output.failed_read} \
        -h {output.html} -j {output.json}

        # fastp does not fsync. Force the flush so a failed scratch write fails
        # this rule rather than silently corrupting the input to bbduk.
        sync --data {output.paired1} {output.paired2}

        # Cross-check what landed on disk against what fastp says it wrote.
        n_exp=$(grep -A2 '"read1_after_filtering"' {output.json} | grep '"total_reads"' | grep -o '[0-9][0-9]*')
        n1=$(( $(wc -l < {output.paired1}) / 4 ))
        n2=$(( $(wc -l < {output.paired2}) / 4 ))
        if [ "$n1" -ne "$n_exp" ] || [ "$n2" -ne "$n_exp" ]; then
            echo "ERROR: trimmed FASTQ truncated: fastp reported $n_exp reads/mate, found R1=$n1 R2=$n2" >&2
            exit 1
        fi

        avail=$(df -h "{params.scratch}")
        echo "Available space: $avail"
        """

rule filter_globin_read:
    input:
        r1=rules.trimming.output.paired1,
        r2=rules.trimming.output.paired2,
        ref=config["resources"]["bbduk_ref"] 
    output:
        out1=temp(os.path.join(os.environ["SCRATCH"], f"{sampleID}_trim_paired_filtered_R1.fq")),
        out2=temp(os.path.join(os.environ["SCRATCH"], f"{sampleID}_trim_paired_filtered_R2.fq")),
        #out1=os.path.join(output_dir, f"{sampleID}_trim_paired_filtered_R1.fq"),
        #out2=os.path.join(output_dir, f"{sampleID}_trim_paired_filtered_R2.fq"),
        stats=os.path.join(output_dir,f'{sampleID}_globin_stats.txt')
    params:
        outdir=output_dir,
        sampleID=f'{sampleID}',
        globinread=os.path.join(output_dir,f'{sampleID}_globinread.fq'),
        bbduk_log=os.path.join(output_dir,f'{sampleID}_bbduk.log')
    threads: 10
    singularity:
        "/hpf/largeprojects/tgnode/sandbox/yliang_analysis/tools/containers/sr_rnaseq.sif"
    shell:
        filter_cmd


rule salmon_quant:
    input:
        r1=rules.filter_globin_read.output.out1,
        r2=rules.filter_globin_read.output.out2,
        #r1 = os.path.join(output_dir, f"{sampleID}_R1_trim_paired.fq.gz"),
        #r2 = os.path.join(output_dir, f"{sampleID}_R2_trim_paired.fq.gz"),
        index = config["resources"]["salmon_index"] 
    output:
        quant_sf = os.path.join(output_dir, "salmon115", "quant.sf"),
        out_dir = directory(os.path.join(output_dir, "salmon115"))
    params:
        lib_type = "A", 
        extra = "--seqBias --gcBias"
    threads: 16
    log:
        os.path.join(output_dir, "logs", f"{sampleID}_salmon.log")
    singularity:
        "/hpf/largeprojects/tgnode/sandbox/yliang_analysis/tools/containers/sr_rnaseq.sif"
    shell:
        """
        salmon quant \
            -i {input.index} \
            -l {params.lib_type} \
            -1 {input.r1} \
            -2 {input.r2} \
            -p 16 \
            -o {output.out_dir} \
            {params.extra} \
            > {log} 2>&1
        """

rule star:
    input:
        star_input1=rules.filter_globin_read.output.out1 if trimming else rules.fastq_prep.output.fastq_r1,
        star_input2=rules.filter_globin_read.output.out2 if trimming else rules.fastq_prep.output.fastq_r2,
        config=sample_config
    output:
        tx_align=temp(os.path.join(output_dir,f'{sampleID}Aligned.toTranscriptome.out.bam')),
        genome_align=temp(os.path.join(output_dir,f'{sampleID}Aligned.sortedByCoord.out.bam')),
        junc_file=os.path.join(output_dir,f'{sampleID}SJ.out.tab'),
        log_out=temp(os.path.join(output_dir,f'{sampleID}Log.out')),
        log_final=os.path.join(output_dir,f'{sampleID}Log.final.out'),
        s1_out_wig=temp(os.path.join(output_dir,f'{sampleID}Signal.Unique.str1.out.wig')),
        s2_out_wig=temp(os.path.join(output_dir,f'{sampleID}Signal.Unique.str2.out.wig')),
        s1_out_multiwig=temp(os.path.join(output_dir,f'{sampleID}Signal.UniqueMultiple.str1.out.wig')),
        s2_out_multiwig=temp(os.path.join(output_dir,f'{sampleID}Signal.UniqueMultiple.str2.out.wig')),
        progress_out=temp(os.path.join(output_dir,f'{sampleID}Log.progress.out')),
        star_genome=temp(directory(os.path.join(output_dir,f'{sampleID}_STARgenome'))),
        star_pass=temp(directory(os.path.join(output_dir,f'{sampleID}_STARpass1')))
        #star_tmp=directory(os.path.join(output_dir,f'{sampleID}_STARtmp')),
    params:    
        index=config["resources"]["star_index"],
        gtf=config["resources"]["gtf"],
        prefix=os.path.join(output_dir,f'{sampleID}'),
        scratch=SCRATCH,
        outdir=output_dir,
        sampleID=f'{sampleID}',
        wasp=lambda wildcards, input: (
        f"--varVCFfile {config['sampleInfo']['vcf']} --waspOutputMode SAMtag"
        if config.get("sampleInfo", {}).get("vcf")
        else ""),
        muscle_seed=lambda wildcards, input: (
        "--seedPerWindowNmax 30"
        if config.get("sampleInfo", {}).get("tissue").lower() == "muscle"
        else "")
    singularity:
        "/hpf/largeprojects/tgnode/sandbox/yliang_analysis/tools/containers/sr_rnaseq.sif"
    threads: 16
    shell:
        """
        rand=$(openssl rand -hex 6)
        tmpdir={params.scratch}/{params.sampleID}_${{rand}}

        avail=$(df -h "{params.scratch}")
        echo "Available space: $avail"
        echo "Open file limit: $(ulimit -n)"

        STAR --runMode alignReads \
            --outTmpDir "$tmpdir" \
            --runThreadN 16 \
            --outSJfilterDistToOtherSJmin 0 0 0 0 \
            --readFilesIn {input.star_input1} {input.star_input2} \
            --genomeDir {params.index} \
            --outFileNamePrefix {params.prefix} \
            --twopassMode Basic \
            --sjdbGTFfile {params.gtf} \
            --outFilterType BySJout \
            --limitSjdbInsertNsj 1200000 \
            --outSAMstrandField intronMotif \
            --outFilterIntronMotifs None \
            --alignSoftClipAtReferenceEnds Yes \
            --quantMode TranscriptomeSAM GeneCounts \
            --outSAMtype BAM SortedByCoordinate \
            --outSAMattrRGline ID:{params.sampleID} SM:{params.sampleID} PL:ILLUMINA \
            --outSAMattributes All \
            --outSAMunmapped Within \
            --outSAMprimaryFlag AllBestScore \
            --chimSegmentMin 15 \
            --chimJunctionOverhangMin 15 \
            --chimOutType Junctions \
            --chimMainSegmentMultNmax 1 \
            --outWigType wiggle \
            --genomeLoad NoSharedMemory \
            {params.muscle_seed} \
            {params.wasp}
        """

rule create_junction_table:
    input:
        rules.star.output.junc_file
    output:
        junction_table=os.path.join(output_dir,f'{sampleID}_SJtable.tsv'),
        site_table=os.path.join(output_dir,f'{sampleID}_splice_sites.tsv')
    params:
        outdir=output_dir,
        sampleID=f'{sampleID}',
        script=os.path.join(workflow.basedir, "scripts/create_junction_site_table.py")
    singularity:
        "/hpf/largeprojects/tgnode/sandbox/yliang_analysis/tools/containers/sr_rnaseq.sif"
    shell:
        """
        python3 {params.script} \
            --junction {input} \
            --output {params.outdir}
        """


rule rsem:
    input:
        rules.star.output.tx_align
    output:
        gene=os.path.join(output_dir,f'{sampleID}.genes.results'),
        isoform=os.path.join(output_dir,f'{sampleID}.isoforms.results'),
        stat=temp(directory(os.path.join(output_dir,f'{sampleID}.stat')))
    params:
        forward_prob="0",
        index=config["resources"]["rsem_index"],
        max_len="1000",
        scratch=SCRATCH,
        outdir=output_dir,
        sampleID=f'{sampleID}'
    threads: 12
    singularity:
        "/hpf/largeprojects/tgnode/sandbox/yliang_analysis/tools/containers/sr_rnaseq.sif"
    shell:
        """
        rand=$(openssl rand -hex 5)
        tmpdir={params.scratch}/{params.sampleID}_${{rand}}

        avail=$(df -h "{params.scratch}")
        echo "Available space: $avail"

        rsem-calculate-expression \
            --bam \
            --num-threads 12 \
            --fragment-length-max {params.max_len} \
            --no-bam-output \
            --paired-end \
            --forward-prob {params.forward_prob} \
            --temporary-folder "$tmpdir" \
            {input} \
            {params.index} {params.outdir}/{params.sampleID}
        """


rule to_bigwig:
    input:
        s1=rules.star.output.s1_out_wig,
        s2=rules.star.output.s2_out_wig,
        chrom_info=config["resources"]["chrom_info"]
    output:
        bw1=os.path.join(output_dir,f'{sampleID}Signal.Unique.str1.out.bw'),
        bw2=os.path.join(output_dir,f'{sampleID}Signal.Unique.str2.out.bw')
    params:
        outdir=output_dir,
        sampleID=f'{sampleID}'
    threads: 6
    singularity:
        "/hpf/largeprojects/tgnode/sandbox/yliang_analysis/tools/containers/sr_rnaseq.sif"
    shell:
        """
        wigToBigWig {input.s1} {input.chrom_info} {output.bw1}
        wigToBigWig {input.s2} {input.chrom_info} {output.bw2}
        """

rule mark_duplicates:
    input:
        rules.star.output.genome_align
    output:
        bam=os.path.join(output_dir,f'{sampleID}.mdup.bam'),
        metrics=os.path.join(output_dir,f'{sampleID}.duplicate.metrics')
    params:
        outdir=output_dir,
        sampleID=f'{sampleID}'
    threads: 10
    singularity:
        "/hpf/largeprojects/tgnode/sandbox/yliang_analysis/tools/containers/sr_rnaseq.sif"
    shell:
        """
        picard MarkDuplicates -Xmx32g \
        -I {input} -O {output.bam} -M {output.metrics} \
        -ASSUME_SORT_ORDER coordinate 

        samtools index {output.bam}
        """

rule ercc_sirv:
    input:
        gene=rules.rsem.output.gene,
        iso=rules.rsem.output.isoform,
        ercc_conc=config["resources"]["spikeins"]
    params:
        outdir=output_dir,
        sampleID=f'{sampleID}',
        script=os.path.join(workflow.basedir, "scripts/spikeins.py")
    output:
        ercc=os.path.join(output_dir,f'{sampleID}.ercc.txt'),
        sirv=os.path.join(output_dir,f'{sampleID}.sirv.txt')
    threads: 4
    singularity:
        "/hpf/largeprojects/tgnode/sandbox/yliang_analysis/tools/containers/sr_rnaseq.sif"
    shell:
        """
        python3 {params.script} --gene {input.gene} --iso {input.iso} --actual {input.ercc_conc} \
            --samplename {params.sampleID} --analysis_path {params.outdir}
        """


rule fastqc:
    input:
        read_input1 = rules.filter_globin_read.output.out1 if trimming else rules.fastq_prep.output.fastq_r1,
        read_input2 = rules.filter_globin_read.output.out2 if trimming else rules.fastq_prep.output.fastq_r2
    output:
        report1=os.path.join(output_dir, f'{sampleID}_R1_fastqc/fastqc_data.txt'),
        report2=os.path.join(output_dir, f'{sampleID}_R2_fastqc/fastqc_data.txt')
    params:
        outdir=output_dir
    threads: 8
    singularity:
        "/hpf/largeprojects/tgnode/sandbox/yliang_analysis/tools/containers/sr_rnaseq.sif"
    shell:
        """
        fastqc -t 2 -q --extract -o {params.outdir} {input.read_input1} {input.read_input2}

        current_dirname1={params.outdir}/$(basename {input.read_input1} | cut -d '.' -f1)_fastqc
        target_dirname1=$(dirname {output.report1})

        current_dirname2={params.outdir}/$(basename {input.read_input2} | cut -d '.' -f1)_fastqc
        target_dirname2=$(dirname {output.report2})

        if [[ "$current_dirname1" != "$target_dirname1" ]]
        then
          rm -r $target_dirname1
          rm -r $target_dirname2
          mv $current_dirname1 $target_dirname1
          mv $current_dirname2 $target_dirname2
        fi
        """

rule rustqc:
    input:
        rules.mark_duplicates.output.bam
    output:
        rustqc_dir = directory(os.path.join(output_dir, "rustqc"))
    params:
        gtf=config["resources"]["gtf"],
        outdir=output_dir,
        sampleID=f'{sampleID}'
    threads: 12
    singularity:
        "/hpf/largeprojects/tgnode/sandbox/yliang_analysis/tools/containers/rustqc.sif"
    shell:
        """
        cd {params.outdir}

        mkdir rustqc
        
        cd rustqc

        rustqc rna {input} \
        --gtf {params.gtf} \
        -p -s reverse -t 12 --sample-name {params.sampleID}
        """

rule predict_tissue:
    input:
        rsem_output=rules.rsem.output.gene
    output:
        predict_tisse_dir=os.path.join(output_dir, f'{sampleID}_tissuePred.txt')
    params:
        outdir=output_dir,
        sampleID=f'{sampleID}',
        script=os.path.join(workflow.basedir, "scripts/predict.py"),
        # The globin panel only means anything where the reference still carries
        # the globin loci. config.yaml points at a masked build, config_unmask.yaml
        # does not; on a masked run every globin TPM is ~0 whatever was in the
        # tube, and purity.py must say NA rather than report a reassuring 0.
        #
        # Note there is deliberately NO --score-as here. predict.py reads the
        # submitted tissue from config.yaml itself; passing it as an override
        # would mark every sample human-asserted and disable the agreement check.
        globin_masked="--globin-masked" if globin_masked else ""
    threads: 8
    singularity:
        "/hpf/largeprojects/tgnode/sandbox/yliang_analysis/tools/containers/sr_rnaseq.sif"
    shell:
        """
        python {params.script} {params.globin_masked} {input} > {output}

        """

rule multiqc:
    input:
        star_file = rules.star.output.log_final,
        bam_file = rules.mark_duplicates.output.bam,
        duplicate_file = rules.mark_duplicates.output.metrics,
        ercc_file = rules.ercc_sirv.output.ercc,
        sirv_file = rules.ercc_sirv.output.sirv,
        fastqc1=rules.fastqc.output.report1,
        fastqc2=rules.fastqc.output.report2,
        rustqc_dir = os.path.join(output_dir, "rustqc")
    params:
        outdir=output_dir,
        multiqc_prefix=sampleID+'multiqc'
    threads: 8
    singularity:
        "/hpf/largeprojects/tgnode/sandbox/yliang_analysis/tools/containers/sr_rnaseq.sif"
    output:
        multiqc_html="/".join([output_dir, sampleID+"multiqc.html"])
    shell:
        """
        cd {params.outdir}
        mkdir err_out_folder
        find . -name '*.err' -exec mv -t ./err_out_folder/ {{}} +
        find . -name '*.out' -exec mv -t ./err_out_folder/ {{}} +

        multiqc {params.outdir} -o {params.outdir} -n {params.multiqc_prefix}
        """

rule sum_metrics:
    input:
        multiqc_html = os.path.join(output_dir, sampleID+"multiqc.html"),
        salmon_output = rules.salmon_quant.output.quant_sf,
        predict_tissue_output=rules.predict_tissue.output,
        bigwig1 = rules.to_bigwig.output.bw1,
        bigwig2 = rules.to_bigwig.output.bw2,
        junction_table=rules.create_junction_table.output.junction_table,
        junction_site_table=rules.create_junction_table.output.site_table
    output:
        qu_sum=os.path.join(output_dir,f'{sampleID}.qcsum.txt')
    threads: 4
    params:
        outdir=output_dir,
        sampleID=f'{sampleID}',
        script=os.path.join(workflow.basedir, "scripts/sum_qc_metrics.py")
    singularity:
        "/hpf/largeprojects/tgnode/sandbox/yliang_analysis/tools/containers/sr_rnaseq.sif"
    shell:
        """
        python3 {params.script} -s {params.sampleID} -i {params.outdir} --outdir {params.outdir}
        """
    