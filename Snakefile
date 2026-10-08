import os
import pandas as pd
import math
import re

tissue=config["sampleInfo"]["tissue"]
sampleID=config["sampleInfo"]["sampleID"]
output_dir=config["sampleInfo"]["dir"]
sample_config=os.path.join(config["sampleInfo"]["dir"], "config.yaml")

# FRASER2 single-sample-vs-frozen-cohort splicing outlier calling. The control
# cohort is a fixed, read-only bundle (built by fraser_cohort.R build/add) at
# <cohort_root>/<tissue>/fraser_objs -- never written to by a per-sample run, so
# any number of cases can be scored against it concurrently and every case is
# measured against a panel that does not move under it. A tissue with no such
# bundle yet is skipped (not failed) rather than blocking the sample.
fraser_enabled_cfg = config["script_parameters"].get("FRASER", False)
fraser_cohort_dir = os.path.join(config["fraser"]["cohort_root"], tissue, "fraser_objs") if fraser_enabled_cfg else None
FRASER_ENABLED = bool(fraser_enabled_cfg) and os.path.isdir(fraser_cohort_dir)
if fraser_enabled_cfg and not FRASER_ENABLED:
    print(f"NOTE: no FRASER control cohort for tissue '{tissue}' at {fraser_cohort_dir} "
          f"-- skipping FRASER for this sample.")

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

#outrider check
outrider = config["script_parameters"]["OUTRIDER"]

VCF_ENABLED = bool(config["sampleInfo"].get("vcf"))
#include: "rules/alignment_expr_qc.smk"

rule all:
    input:
        os.path.join(output_dir,f'{sampleID}.qcsum.txt'),
        *(
            [os.path.join(output_dir, f"{sampleID}_OUTRIDER_results_annotated.tsv")]
            if outrider
            else []
        ),
        *(
            [os.path.join(output_dir, f"{sampleID}_FRASER_results_annotated.filtered.tsv")]
            if FRASER_ENABLED
            else []
        )


if VCF_ENABLED:
    rule vcf_filter:
        input:
            vcf=config["sampleInfo"]["vcf"]
        output:
            vcf_filter_result=os.path.join(output_dir,f'model_variants.tsv')
        threads: 8
        params:
            outdir=output_dir,
            sampleID=f'{sampleID}',
            script=os.path.join(workflow.basedir, "scripts/trio_rare_filter.py"),
            # Explicit genomic-proband -> RNAseq-sample linkage. sampleInfo.sampleID
            # (RNAseq ID, e.g. SKS_GS_092_03) and the VCF's own sample naming (e.g.
            # DRAGEN quad ids like CS4C-070111_EXP_0001) are NOT the same string and
            # nothing else in the pipeline cross-checks them -- see
            # trio-rare-filter-script memory. Without vcf_proband_id/ped, trio_rare_filter.py
            # falls back to its own *_01/*_02/*_03 suffix autodetection (unchanged
            # from before), which silently fails to find a proband on VCFs that
            # don't follow that convention (e.g. DRAGEN quad ids) and can't resolve
            # siblings for quad mode at all.
            extra=lambda wc: " ".join(filter(None, [
                f"-p {config['sampleInfo']['vcf_proband_id']}"
                if config.get("sampleInfo", {}).get("vcf_proband_id") else "",
                f"--ped {config['sampleInfo']['ped']}"
                if config.get("sampleInfo", {}).get("ped") else "",
            ]))
        singularity:
            config["resources"]["container"]
        shell:
            """
            python3 {params.script} -v {input} -o {params.outdir} {params.extra}
            """

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
        config["resources"]["container"]
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
        config["resources"]["container"]
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
        config["resources"]["container"]
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
        config["resources"]["container"]
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

rule salmon2RSEM:
    input:
        quant_sf = os.path.join(output_dir, "salmon115", "quant.sf"),
        salmon_dir = directory(os.path.join(output_dir, "salmon115"))
    output:
        salmon2RSEM_output = os.path.join(output_dir, f'{sampleID}_salmon_rsem.results')
    params:
        outdir=output_dir,
        sampleID=f'{sampleID}',
        tx2gene=config["resources"]["tx2gene"],
        script=os.path.join(workflow.basedir, "scripts/salmon2RSEM.R")
    singularity:
        config["resources"]["container"]
    shell:
        """
        Rscript --vanilla {params.script} \
        --tx2gene {params.tx2gene} \
        --salmon_dir {input.salmon_dir} \
        --output_dir {params.outdir}
        """

if VCF_ENABLED:
    rule normalize_vcf_for_star:
        # STAR's WASP matches VCF #CHROM against the genome fasta's contig
        # names literally. This genome has no "chr" prefix (MT, not chrM);
        # clinical VCFs do. Rename map is built from chrName.txt so it
        # tracks the genome build instead of assuming a fixed convention.
        input:
            vcf=config["sampleInfo"]["vcf"],
            star_chrnames=os.path.join(config["resources"]["star_index"], "chrName.txt")
        output:
            vcf=temp(os.path.join(output_dir, f'{sampleID}_star_wasp.vcf')),
            chr_map=temp(os.path.join(output_dir, f'{sampleID}_star_chr_map.tsv'))
        singularity:
            config["resources"]["container"]
        shell:
            """
            awk 'BEGIN{{OFS="\\t"}} {{
                n=$1
                if (n == "MT") {{ print "chrM", "MT"; print "chrMT", "MT" }}
                else {{ print "chr"n, n }}
            }}' {input.star_chrnames} > {output.chr_map}

            bcftools annotate --rename-chrs {output.chr_map} {input.vcf} \
                -Ov -o {output.vcf}
            """

rule star:
    input:
        star_input1=rules.filter_globin_read.output.out1 if trimming else rules.fastq_prep.output.fastq_r1,
        star_input2=rules.filter_globin_read.output.out2 if trimming else rules.fastq_prep.output.fastq_r2,
        config=sample_config,
        wasp_vcf=rules.normalize_vcf_for_star.output.vcf if VCF_ENABLED else []
    output:
        tx_align=temp(os.path.join(output_dir,f'{sampleID}Aligned.toTranscriptome.out.bam')),
        genome_align=temp(os.path.join(output_dir,f'{sampleID}Aligned.out.bam')),
        junc_file=os.path.join(output_dir,f'{sampleID}SJ.out.tab'),
        log_out=temp(os.path.join(output_dir,f'{sampleID}Log.out')),
        log_final=os.path.join(output_dir,f'{sampleID}Log.final.out'),
        progress_out=temp(os.path.join(output_dir,f'{sampleID}Log.progress.out')),
        star_genome=temp(directory(os.path.join(output_dir,f'{sampleID}_STARgenome'))),
        star_pass=temp(directory(os.path.join(output_dir,f'{sampleID}_STARpass1')))
        #s1_out_wig=temp(os.path.join(output_dir,f'{sampleID}Signal.Unique.str1.out.wig')),
        #s2_out_wig=temp(os.path.join(output_dir,f'{sampleID}Signal.Unique.str2.out.wig')),
        #s1_out_multiwig=temp(os.path.join(output_dir,f'{sampleID}Signal.UniqueMultiple.str1.out.wig')),
        #s2_out_multiwig=temp(os.path.join(output_dir,f'{sampleID}Signal.UniqueMultiple.str2.out.wig')),
        #star_tmp=directory(os.path.join(output_dir,f'{sampleID}_STARtmp')),
    params:    
        index=config["resources"]["star_index"],
        gtf=config["resources"]["gtf"],
        prefix=os.path.join(output_dir,f'{sampleID}'),
        scratch=SCRATCH,
        outdir=output_dir,
        sampleID=f'{sampleID}',
        wasp=lambda wildcards, input: (
        f"--varVCFfile {input.wasp_vcf} --waspOutputMode SAMtag"
        if VCF_ENABLED
        else ""),
        muscle_seed=lambda wildcards, input: (
        "--seedPerWindowNmax 30"
        if config.get("sampleInfo", {}).get("tissue").lower() == "muscle"
        else "")
    singularity:
        config["resources"]["container"]
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
            --outSAMtype BAM Unsorted \
            --outSAMattrRGline ID:{params.sampleID} SM:{params.sampleID} PL:ILLUMINA \
            --outSAMattributes All \
            --outSAMunmapped Within \
            --outSAMprimaryFlag AllBestScore \
            --chimSegmentMin 15 \
            --chimJunctionOverhangMin 15 \
            --chimOutType Junctions \
            --chimMainSegmentMultNmax 1 \
            --genomeLoad NoSharedMemory \
            {params.muscle_seed} \
            {params.wasp}
        
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
    singularity:
        config["resources"]["container"]
    shell:
        """
        samtools sort -@ 8 -m 6G -o {output.genome_align_sorted} {input.genome_align}
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
        config["resources"]["container"]
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
        config["resources"]["container"]
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


#rule to_bigwig:
#    input:
#        s1=rules.star.output.s1_out_wig,
#        s2=rules.star.output.s2_out_wig,
#        chrom_info=config["resources"]["chrom_info"]
#    output:
#        bw1=os.path.join(output_dir,f'{sampleID}Signal.Unique.str1.out.bw'),
#        bw2=os.path.join(output_dir,f'{sampleID}Signal.Unique.str2.out.bw')
#    params:
#        outdir=output_dir,
#        sampleID=f'{sampleID}'
#    threads: 6
#    singularity:
#        "/hpf/largeprojects/tgnode/sandbox/yliang_analysis/tools/containers/sr_rnaseq.sif"
#    shell:
#        """
#        wigToBigWig {input.s1} {input.chrom_info} {output.bw1}
#        wigToBigWig {input.s2} {input.chrom_info} {output.bw2}
#        """

rule mark_duplicates:
    input:
        rules.sort.output.genome_align_sorted
    output:
        bam=os.path.join(output_dir,f'{sampleID}.mdup.bam'),
        # Was previously written by `samtools index` below but never declared,
        # so nothing downstream could depend on it. FRASER's BAM strand-recovery
        # step does indexed region reads and needs it tracked as a real output.
        bai=os.path.join(output_dir,f'{sampleID}.mdup.bam.bai'),
        metrics=os.path.join(output_dir,f'{sampleID}.duplicate.metrics')
    params:
        outdir=output_dir,
        sampleID=f'{sampleID}'
    threads: 10
    singularity:
        config["resources"]["container"]
    shell:
        """
        picard MarkDuplicates -Xmx32g \
        -I {input} -O {output.bam} -M {output.metrics} \
        -ASSUME_SORT_ORDER coordinate

        samtools index {output.bam}
        """

rule bam_to_bigwig:
    input:
        bam=rules.mark_duplicates.output.bam
    output:
        bw1=os.path.join(output_dir,f'{sampleID}Signal.Unique.plus.bw'),
        bw2=os.path.join(output_dir,f'{sampleID}Signal.Unique.minus.bw')
    params:
        outdir=output_dir,
        sampleID=f'{sampleID}'
    threads: 10
    singularity:
        config["resources"]["container"]"
    shell:
        """
        # Biological + strand
        bamCoverage \
        -b {input} \
        -o {output.bw1} \
        --filterRNAstrand reverse \
        --normalizeUsing CPM \
        --minMappingQuality 255 \
        --binSize 1 -p 10

        # Biological - strand
        bamCoverage \
        -b {input} \
        -o {output.bw2} \
        --filterRNAstrand forward \
        --normalizeUsing CPM \
        --minMappingQuality 255 \
        --binSize 1 -p 10
        """

if FRASER_ENABLED:
    fraser_workdir = os.path.join(output_dir, "fraser_work")

    rule fraser_add:
        # Adds this one sample's counts to a COPY of the frozen cohort's counts
        # object under fraser_workdir -- the cohort bundle itself is never
        # written to (--cohort-dir is read-only), so any number of samples can
        # run this concurrently against the same tissue panel. Old samples are
        # never re-read from BAM: fraser_cohort.R resolves their split-count
        # files and non-split caches straight out of the cohort dir (the latter
        # via symlinks into fraser_workdir/cache, since FRASER only knows how to
        # find its non-split cache under its own workingDir).
        input:
            sj=rules.star.output.junc_file,
            bam=rules.mark_duplicates.output.bam,
            bai=rules.mark_duplicates.output.bai
        output:
            touch(os.path.join(fraser_workdir, ".add.done"))
        params:
            cohort_dir=fraser_cohort_dir,
            workdir=fraser_workdir,
            sampleID=f'{sampleID}',
            script=os.path.join(workflow.basedir, "scripts/fraser_cohort.R")
        threads: 8
        singularity:
            config["resources"]["container"]
        shell:
            # --vanilla is mandatory here, not cosmetic: --use-singularity mounts
            # the host $HOME, so a plain Rscript picks up host ~/.Renviron /
            # ~/.Rprofile, which point R_LIBS_USER at a host package library
            # whose XML.so links against a host libxml2 absent from the
            # container -- library(FRASER) then fails on dyn.load. --vanilla
            # skips those startup files so R sees only the image's packages.
            """
            Rscript --vanilla {params.script} add \
                --cohort-dir {params.cohort_dir} \
                --workdir    {params.workdir} \
                --sample-id  {params.sampleID} \
                --bam        {input.bam} \
                --sj         {input.sj} \
                --threads    {threads}
            """

    rule fraser_fit:
        # Fits the cohort+sample counts object (q pinned per tissue in config,
        # not re-estimated per case -- with a frozen panel and one varying
        # sample, OHT re-estimation would be the only remaining source of
        # cross-sample drift) and writes this sample's result tables straight
        # into output_dir, flat and sample-prefixed via --out-dir +
        # --results-samples, so they sit alongside the rest of the pipeline's
        # per-sample deliverables with no collision risk.
        input:
            rules.fraser_add.output
        output:
            results=os.path.join(output_dir, f'{sampleID}_results.tsv'),
            results_all=os.path.join(output_dir, f'{sampleID}_results_all.tsv.gz'),
            # Gene-level table intentionally NOT requested here: fraser_cohort.R
            # only writes it behind --gene-table, which this rule doesn't pass
            # (that extraction step has OOM-killed cohort-wide fits at 64G even
            # at all=FALSE). Declaring it as an output with no flag to produce
            # it made Snakemake treat every otherwise-successful fit as failed,
            # which then ran this rule's own `rm -rf savedObjects/cohort`
            # cleanup and broke the next fraser_add/fraser_fit cycle.
            private_novel=os.path.join(output_dir, f'{sampleID}_private_novel.tsv'),
            cohort_private_novel=os.path.join(output_dir, f'{sampleID}_cohort_private_novel.tsv'),
            # The FITTED fds (annotated, filtered, PCA-corrected) -- kept on
            # purpose by the cleanup step below, specifically so fraser_parse
            # can point --fds-file at it for splice-impact annotation
            # (potentialImpact/causesFrameshift/UTR_overlap need the fds
            # itself, not just the results table: classifying exon-skip vs
            # truncation needs neighbouring junctions' counts). Declaring it
            # here -- rather than leaving it an untracked side effect of this
            # rule's shell script -- makes fraser_parse's dependency on it
            # explicit in the DAG. Path verified 2026-09-03 against a real
            # fitted object: loadFraserDataSet(file=<this>) loads with
            # hgnc_symbol/annotatedJunction present.
            fds_fit=os.path.join(fraser_workdir, "savedObjects", "cohort_fit", "fds-object.RDS")
        params:
            cohort_dir=fraser_cohort_dir,
            workdir=fraser_workdir,
            out_dir=output_dir,
            sampleID=f'{sampleID}',
            q=config["fraser"]["q"][tissue],
            # TESTING threshold (filterExpressionAndVariability), not the
            # frozen cohort's DISCOVERY threshold (--sj-min-count, pinned at
            # build time and unrelated) -- changed from the fraser_cohort.R
            # default (20) to 10, 2026-09-22. A fit at a different --min-count
            # than the cohort was pinned with is now a WARNING, not a hard
            # error (fraser_cohort.R change, same day), so this proceeds even
            # though the 09182026 cohort bundles were themselves fit at 20.
            mincount=10,
            script=os.path.join(workflow.basedir, "scripts/fraser_cohort.R")
        threads: 12
        singularity:
            config["resources"]["container"]
        shell:
            """
            Rscript --vanilla {params.script} fit \
                --cohort-dir {params.cohort_dir} \
                --workdir    {params.workdir} \
                --results-samples {params.sampleID} \
                --out-dir    {params.out_dir} \
                --q          {params.q} \
                --min-count  {params.mincount} \
                --threads    {threads}

            # savedObjects/cohort is the re-assembled COUNTS object -- cheap to
            # regenerate (fraser_add just re-links the cohort's frozen caches,
            # no BAM re-read) and it scales with cohort size, not with this one
            # sample (~130-150MB at the current LCL cohort size). Prune it once
            # fit has succeeded and the result tables are on disk.
            # savedObjects/cohort_fit (the actual PCA fit, expensive to redo)
            # and the symlinked non-split cache are kept -- a later
            # annotate_splice_impact.R pass on this sample needs the fitted fds.
            rm -rf {params.workdir}/savedObjects/cohort
            """

    rule fraser_parse:
        # Mirrors outrider_parse: takes fraser_fit's ALL-junctions table (not
        # the outliers-only one -- parse_fraser.R's nominal-p rescue gate needs
        # the sub-threshold rows all=FALSE would have dropped, same reason
        # fraser_cohort.R's own closing message points at this exact file) and
        # writes the gene/HPO/variant-annotated, triaged table. Same HPO/VCF
        # sentinel convention as outrider_parse: "NULL" (a literal string, not
        # Python None) when neither is configured -- parse_fraser.R checks
        # file.exists() on it, so a missing sentinel file just degrades to
        # NA/empty columns rather than erroring.
        input:
            results_all=rules.fraser_fit.output.results_all,
            fds_fit=rules.fraser_fit.output.fds_fit,
            # Per-sample private-novel table (already scoped to this sample --
            # NOT the cohort-wide cohort_private_novel output, which is bigger
            # and would need parse_fraser.R's -s-based re-filtering for no
            # benefit here). Added 2026-09-22 alongside parse_fraser.R's new
            # -N mode, so the same gene annotation / blacklist /
            # potential_artifact / HPO overlap the filtered junctions get is
            # also attached to the junctions FRASER never got to test.
            private_novel=rules.fraser_fit.output.private_novel,
            filter_variant=lambda wc: (
                rules.vcf_filter.output.vcf_filter_result
                if VCF_ENABLED
                else []
            )
        output:
            full=os.path.join(output_dir, f'{sampleID}_FRASER_results_annotated.tsv.gz'),
            filtered=os.path.join(output_dir, f'{sampleID}_FRASER_results_annotated.filtered.tsv'),
            private_novel_annotated=os.path.join(output_dir, f'{sampleID}_FRASER_private_novel_annotated.tsv')
        threads: 4
        params:
            outdir=output_dir,
            sampleID=f'{sampleID}',
            hpo=config["sampleInfo"]["hpoterm"] or "NULL",
            v=lambda wc: (
                rules.vcf_filter.output.vcf_filter_result
                if VCF_ENABLED
                else "NULL"
            ),
            txdb=config["fraser"]["txdb"],
            # -o is a BASE path, not a real output file: parse_fraser.R strips
            # its extension and writes <base>.tsv.gz / <base>.filtered.tsv
            # (an .xlsx workbook is also written if -o ends .xlsx, via a
            # separate tsv_to_xlsx.py converter this pipeline doesn't carry --
            # deliberately left off here to keep the R dep set at
            # tidyverse+optparse, same as outrider_parse).
            out=os.path.join(output_dir, f'{sampleID}_FRASER_results_annotated.tsv'),
            script=os.path.join(workflow.basedir, "scripts/parse_fraser.R")
        singularity:
            config["resources"]["container"]
        shell:
            # --txdb + --fds-file turn on splice-impact annotation
            # (potentialImpact/causesFrameshift/UTR_overlap), applied to the
            # FILTERED rows only, matching DROP -- this is what pulls in
            # fraser_impact_annotation.R (sourced lazily, sibling to this
            # script) and needs the FITTED fds, not the counts object (see the
            # comment on fraser_fit's fds_fit output).
            # -N/--private-novel-out run independently of -r/-s's own table --
            # see parse_fraser.R's "--private-novel" option / the
            # "Annotating the private-novel table" section of README_fraser.md.
            """
            Rscript --vanilla {params.script} -r {input.results_all} \
            -s {params.sampleID} -p {params.hpo} \
            -o {params.out} -V {params.v} \
            --txdb {params.txdb} --fds-file {input.fds_fit} \
            -N {input.private_novel} --private-novel-out {output.private_novel_annotated}
            """

rule ercc_sirv:
    input:
        gene=rules.salmon2RSEM.output.salmon2RSEM_output,
        iso=rules.salmon_quant.output.quant_sf,
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
        config["resources"]["container"]
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
        config["resources"]["container"]
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
        config["resources"]["rustqc_container"]
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
        #rsem_output=rules.rsem.output.gene
        rules.salmon2RSEM.output.salmon2RSEM_output
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
        config["resources"]["container"]
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
        bigwig1 = rules.bam_to_bigwig.output.bw1,
        bigwig2 = rules.bam_to_bigwig.output.bw2,
        rustqc_dir = os.path.join(output_dir, "rustqc")
    params:
        outdir=output_dir,
        multiqc_prefix=sampleID+'multiqc'
    threads: 8
    singularity:
        config["resources"]["container"]
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
        junction_table=rules.create_junction_table.output.junction_table,
        junction_site_table=rules.create_junction_table.output.site_table,
        bigwig1 = rules.bam_to_bigwig.output.bw1,
        bigwig2 = rules.bam_to_bigwig.output.bw2,
    output:
        qu_sum=os.path.join(output_dir,f'{sampleID}.qcsum.txt')
    threads: 4
    params:
        outdir=output_dir,
        sampleID=f'{sampleID}',
        script=os.path.join(workflow.basedir, "scripts/sum_qc_metrics.py")
    singularity:
        config["resources"]["container"]
    shell:
        """
        python3 {params.script} -s {params.sampleID} -i {params.outdir} --outdir {params.outdir}
        """

if outrider:
    rule run_outrider:
        input:
            gene=os.path.join(output_dir, "salmon115", "quant.sf")
        output:
            outrider_result=os.path.join(output_dir,f'{sampleID}_OUTRIDER_results.rds')
        threads: 10
        params:
            outdir=output_dir,
            sampleID=f'{sampleID}',
            out_dir=output_dir,
            cohort_salmon=config["outrider"]["cohort_root"][tissue],
            tx2gene=config["resources"]["tx2gene"],
            txdb=config["outrider"]["txdb"],
            q=config["outrider"]["q"][tissue],
            script=os.path.join(workflow.basedir, "scripts/salmon_OUTRIDER.R")
        singularity:
                config["resources"]["container"]
        shell:
                """
                Rscript --vanilla {params.script} -s {input} \
                    -c {params.cohort_salmon} --txdb {params.txdb} \
                    --q_value {params.q} --tx2gene {params.tx2gene} \
                    -o {output.outrider_result}
                """


    rule outrider_parse:
        input:
            outrider_rds=rules.run_outrider.output.outrider_result,
            filter_variant=lambda wc: (
                rules.vcf_filter.output.vcf_filter_result
                if VCF_ENABLED
                else []
            )
        output:
            outrider_anno_result=os.path.join(output_dir,f'{sampleID}_OUTRIDER_results_annotated.tsv')
        threads: 4
        params:
            outdir=output_dir,
            sampleID=f'{sampleID}',
            tissue=tissue,
            hpo=config["sampleInfo"]["hpoterm"] or "NULL",
            v=lambda wc: (
                rules.vcf_filter.output.vcf_filter_result
                if VCF_ENABLED
                else "NULL"
            ),
            script=os.path.join(workflow.basedir, "scripts/parse_outrider.R")
        singularity:
                config["resources"]["container"]
        shell:
                """
                Rscript --vanilla {params.script} -r {input.outrider_rds} \
                -s {params.sampleID} -p {params.hpo} \
                -o {output.outrider_anno_result} -V {params.v}
                """