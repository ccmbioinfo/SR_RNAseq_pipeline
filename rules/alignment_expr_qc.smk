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
        python_modulefile=config["executables"]["python_modulefile"]
    shell:
        """
        module load orad/2.6.1
        module load samtools/1.10
        module load {params.python_modulefile}

        python scripts/fastq_prep.py --sample_config {input} --output {params.outdir}
        """
 
rule trimming:
    input:
        r1=rules.fastq_prep.output.fastq_r1
        r2=rules.fastq_prep.output.fastq_r2
    output:
        paired1=temp(os.path.join(output_dir,f'{sampleID}_R1_trim_paired.fq.gz')),
        paired2=temp(os.path.join(output_dir,f'{sampleID}_R2_trim_paired.fq.gz')),
        unpaired1=os.path.join(output_dir,f'{sampleID}_unpairR1.fq.gz'),
        unpaired2=os.path.join(output_dir,f'{sampleID}_unpairR2.fq.gz'),
        failed_read=os.path.join(output_dir,f'{sampleID}_failed.fq.gz'),
        json=os.path.join(output_dir,f'{sampleID}fastp.json'),
        html=os.path.join(output_dir,f'{sampleID}fastp.html')
    params:
        outdir=output_dir,
        sampleID=f'{sampleID}',
        fastp=config["executables"]["fastp"]
    shell:
        """
        {params.fastp} -w 10 -i {input.r1} -I {input.r2} -o {output.paired1} -O {output.paired2} \
        --unpaired1 {output.unpaired1} --unpaired2 {output.unpaired2} --failed_out {output.failed_read} \
        -h {output.html} -j {output.json}
        """

rule star:
    input:
        star_input1=rules.trimming.output.paired1 if trimming else rules.fastq_prep.output.fastq_r1,
        star_input2=rules.trimming.output.paired2 if trimming else rules.fastq_prep.output.fastq_r2,
        config=sample_config
    output:
        tx_align=temp(os.path.join(output_dir,f'{sampleID}Aligned.toTranscriptome.out.bam')),
        genome_align=temp(os.path.join(output_dir,f'{sampleID}Aligned.out.bam')),
        junc_file=os.path.join(output_dir,f'{sampleID}SJ.out.tab'),
        log_out=temp(os.path.join(output_dir,f'{sampleID}Log.out')),
        log_final=os.path.join(output_dir,f'{sampleID}Log.final.out'),
        progress_out=temp(os.path.join(output_dir,f'{sampleID}Log.progress.out')),
        star_genome=temp(directory(os.path.join(output_dir,f'{sampleID}_STARgenome'))),
        star_tmp=temp(directory(os.path.join(output_dir,f'{sampleID}_STARtmp'))),
        star_pass=temp(directory(os.path.join(output_dir,f'{sampleID}_STARpass1')))
    params:    
        index=config["resources"]["star_index"],
        gtf=config["resources"]["gtf"],
        prefix=os.path.join(output_dir,"f'{sampleID}"),
        outdir=output_dir,
        sampleID=f'{sampleID}',
        STAR_modulefile=config["executables"]["STAR_modulefile"]
    threads: 15
    run:
        import yaml

        with open(input.config) as f:
            sample_cfg = yaml.safe_load(f)

        vcf = sample_cfg.get("sampleInfo", {}).get("vcf")

        wasp_command = ""
        if vcf:
            wasp_command = f"""
                --varVCFfile {vcf}
                --waspOutputMode SAMtag
            """

        shell(f"""
            module load {params.STAR_modulefile}; STAR --runMode alignReads \
                --runThreadN 15 \
                --readFilesCommand zcat \
                --readFilesIn {input.star_input1} {input.star_input2}\
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
                --outSAMtype BAM Unsorted \
                --outSAMattrRGline ID:{params.sampleID}  SM:{params.sampleID}  PL:ILLUMINA \
                --outSAMattributes All \
                --outSAMunmapped Within \
                --outSAMprimaryFlag AllBestScore \
                --chimSegmentMin 15 \
                --chimJunctionOverhangMin 15 \
                --chimOutType Junctions \
                --chimMainSegmentMultNmax 1 \
                --genomeLoad NoSharedMemory --runThreadN {threads} \
                {wasp_command}
        """)
        
        shell(default_command)

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
        outdir=output_dir,
        sampleID=f'{sampleID}',
        RSEM_modulefile=config["executables"]["RSEM_modulefile"]
    threads: 15
    shell:
        """
        module load {params.RSEM_modulefile}

        rsem-calculate-expression \
            --bam \
            --num-threads 15 \
            --fragment-length-max {params.max_len} \
            --no-bam-output \
            --paired-end \
            --forward-prob {params.forward_prob} \
            {input} \
            {params.index} {params.outdir}/{params.sampleID}
        """


rule sort:
    input:
        genome_align=rules.star.output.genome_align
    output:
        genome_align_sorted = temp(os.path.join(output_dir,f'{sampleID}Aligned.sorted.bam'))
    params:
        outdir=output_dir,
        sampleID=f'{sampleID}'
    threads: 8
    shell:
        """
        module load samtools/1.9

        samtools sort -@ 8 -m 6G -o {output.genome_align_sorted} {input.genome_align}
        """

rule mark_duplicates:
    input:
        rules.sort.output.genome_align_sorted
    output:
        bam=os.path.join(output_dir,f'{sampleID}.mdup.bam'),
        metrics=os.path.join(output_dir,f'{sampleID}.duplicate.metrics')
    params:
        picard=config["executables"]["picard"],
        outdir=output_dir,
        sampleID=f'{sampleID}'
    threads: 10
    shell:
        """
        module load samtools/1.9
        module load java/1.8.0_161

        java -Xmx40g -jar {params.picard} MarkDuplicates \
        I={input} O={output.bam} M={output.metrics} \
        ASSUME_SORT_ORDER=coordinate 

        samtools index {output.bam}
        """


rule rnaseq_metrics:
    input:
        rules.mark_duplicates.output.bam
    output:
        metrics=os.path.join(output_dir,f'{sampleID}.rnaseq_metrics.txt')
    params:
        ribo_int=config["resources"]["ribo_int"],
        refflat=config["resources"]["refflat"],
        strand="SECOND_READ_TRANSCRIPTION_STRAND",
        AS="true",
        picard=config["executables"]["picard"],
        outdir=output_dir,
        sampleID=f'{sampleID}'
    threads: 10
    shell:
        """
        module load java/1.8.0_161

        java -Xmx4096m -jar {params.picard} CollectRnaSeqMetrics \
        REF_FLAT={params.refflat} \
        RIBOSOMAL_INTERVALS={params.ribo_int} \
        STRAND_SPECIFICITY= {params.strand} \
        I={input} \
        O={output.metrics} AS={params.AS}
        """

rule ercc_sirv:
    input:
        gene=rules.rsem.output.gene,
        iso=rules.rsem.output.isoform,
        ercc_conc=config["resources"]["spikeins"],
        python_modulefile=config["executables"]["python_modulefile"]
    params:
        rnaseqc = config["executables"]["rnaseqc"],
        outdir=output_dir,
        sampleID=f'{sampleID}'
    output:
        ercc=os.path.join(output_dir,f'{sampleID}.ercc.txt'),
        sirv=os.path.join(output_dir,f'{sampleID}.sirv.txt')
    threads: 10
    shell:
        """
        module load {input.python_modulefile}

        source python_venv/bin/activate

        python3 scripts/spikeins.py --gene {input.gene} --iso {input.iso} --actual {input.ercc_conc} \
            --samplename {params.sampleID} --analysis_path {params.outdir}
        """

rule rnaseqc:
    input:
        bam = rules.mark_duplicates.output.bam,
        collapsed_gtf = config["resources"]["collapsed_gtf"]
    output:
        os.path.join(output_dir,f'{sampleID}.metrics.tsv')
    params:
        rnaseqc = config["executables"]["rnaseqc"],
        outdir=output_dir,
        sampleID=f'{sampleID}'
    threads: 10
    shell:
        """
        {params.rnaseqc} {input.collapsed_gtf} {input.bam} {params.outdir} -s {params.sampleID}
        rm -f {params.outdir}/*gct
        """

rule insert_size_metrics:
    input:
        rules.mark_duplicates.output.bam
    output:
        metrics=os.path.join(output_dir,f'{sampleID}.insert_metrics.txt'),
        plot=temp(os.path.join(output_dir,f'{sampleID}.insert_metrics.pdf')),
    params:
        picard=config["executables"]["picard"]
    threads: 10
    shell:
        """
        module load java/1.8.0_161
        module load R/4.0.3
 
        java -Xmx4096m -jar {params.picard} CollectInsertSizeMetrics \
        I={input} O={output.metrics} H={output.plot}
        """

rule fastqc:
    input:
        read_input1 = rules.trimming.output.paired1 if trimming else rules.fastq_prep.output.fastq_r1,
        read_input2 = rules.trimming.output.paired2 if trimming else rules.fastq_prep.output.fastq_r2
    output:
        report1=os.path.join(output_dir, f'{sampleID}_R1_fastqc/fastqc_data.txt'),
        report2=os.path.join(output_dir, f'{sampleID}_R2_fastqc/fastqc_data.txt')
    params:
        outdir=output_dir,
    threads: 2
    shell:
        """
        module load fastqc/0.11.5

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

rule alignment_metrics:
    input:
        rules.mark_duplicates.output.bam
    output: os.path.join(output_dir,f'{sampleID}.alignment_metrics.txt')
    params:
        picard=config["executables"]["picard"],
        fasta=config["resources"]["fasta"],
        outdir=output_dir,
        sampleID=f'{sampleID}'
    threads: 10
    shell:
        """
        module load java/1.8.0_161

        java -Xmx4096m -jar {params.picard} CollectAlignmentSummaryMetrics \
        I={input} R={params.fasta} O={output}        
        """
