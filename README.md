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
