# Major update -- 2026/Oct/07

* We are uisng salmon as the main quantification tool (RSEM currently unavilable in 115). Added script salmon2RSEM.R to make the gene output looks like RSEM. Also updated salmon_OUTRIDER.R and spikeins.py to take the RSEM-like salmon output.
* Updates in Snakefile:
    - Due to the temp space issue, STAR no more create sorted bam. Added new rule to create sorted bam and bigwig file
    - Officially added FRASER and OUTRIDER in the pipeline. Use config file to turn it ON/OFF
* Both OUTRIDER and FRASER have pre-defined Q values for different cohort tissues. Q value from selected TG internal cohort using grid-search OUTRIDER::estimateBestQ
* A lof of updates in FRASER:
    - rewrite strand recovery, add control-evidence privacy fix, private-novel annotation
    - fraser_nonsplit_from_sjtab.R (09-18): replace pair-based col4=0 strand recovery with a donor-base/SAM-flag vote; drop the --max-region-reads blind-skip gate for a per-junction time budget + multi-sample retry.
    - fraser_cohort.R (09-18→22): track strand_attempts.tsv.gz and add control_evidence to distinguish "tested absent" from "never tested" controls; pin strand_method as a critical cohort param; expose --quantile-min-expression; add opt-in --private-into-counts; downgrade --min-count mismatch from stop() to warning().
    - parse_fraser.R (09-21→24): read deltaPsi_raw's cohort median directly off --fds-file/--fds-dir independent of --txdb; switch IG/TR artifact detection to a gene-symbol regex; add pvalue_rank; rename gene/hgncSymbol output columns; add --private-novel/--private-novel-out to annotate fraser_cohort.R's private-novel table via the same shared annotation functions as the main filtered table.
* future plans:
    - Arriba fusion
    - New HPO-annotation

# Major update -- 2026/Sep/01

* updated qc sum script, tissue predict, quality check
* salmon2RSEM.R included, but not added into the pipeline yet
* IMPORTANT: salmon is now using v115 cdna based index
* update in major snakefile, hopefully dont reduce the required tmp space
* removed unused config yaml files and requirements file, since we dont use masked genome and venv

# First version -- 2026/July/14

#how to run

setup_python.sh to set up the venv

check the samples.tsv for submission template


source /hpf/largeprojects/tgnode/sandbox/rnaseq_repo/SR_RNAseq_pipeline_env/miniforge3/etc/profile.d/conda.sh
conda activate new_SR_RNAseq_pipeline

#test sample
python /hpf/largeprojects/tgnode/sandbox/yliang_analysis/SR_RNAseq_pipeline/submit_sample.py -p /hpf/largeprojects/tgnode/sandbox/yliang_analysis/SR_RNAseq_pipeline -c /hpf/largeprojects/tgnode/sandbox/yliang_analysis/SR_RNAseq_pipeline/config_unmask.yaml -s /hpf/largeprojects/tgnode/sandbox/yliang_analysis/SR_RNAseq_pipeline/samples.tsv -o /hpf/largeprojects/tgnode/sandbox/yliang_analysis/temp/test

#real samples (tmp)
python /hpf/largeprojects/tgnode/sandbox/yliang_analysis/SR_RNAseq_pipeline/submit_sample.py -p /hpf/largeprojects/tgnode/sandbox/yliang_analysis/SR_RNAseq_pipeline -c /hpf/largeprojects/tgnode/sandbox/yliang_analysis/SR_RNAseq_pipeline/config_unmask.yaml -s /hpf/largeprojects/tgnode/sandbox/rnaseq_repo/TG_RNAseq/test.txt -o /hpf/largeprojects/tgnode/sandbox/rnaseq_repo/TG_RNAseq/


can run at any place as long as you module load and activate the venv
```
sampleID\tfastq-path\ttissue\tvcf\thpo_term
```
