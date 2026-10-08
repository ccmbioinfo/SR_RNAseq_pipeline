#!/usr/bin/env bash
# Parse OUTRIDER + FRASER control-cohort results for every sample in every
# tissue under a control-cohort root, and drop the annotated tables into that
# sample's OWN project folder (TG_RNAseq / ComplexCare / new_GenomeClinic).
# No variant/VCF join here (no -V) -- just OUTRIDER+FRASER, gene annotation,
# and HPO overlap (HPO columns left NA when a sample has no hpoterm/ file).
#
# ONE SLURM JOB PER SAMPLE: this script only submits and exits; it does not do
# the work itself. Samples are independent (each re-reads the same tissue-wide
# results tables), so they run in parallel and a failure isolates to one
# sample instead of the batch. Within a sample's job, OUTRIDER and FRASER run
# as two sequential steps and one failing does not block the other.
# Pattern copied from rnaseq_repo/test/LCL/run_parse_all.sh.
#
# BUILT TO SURVIVE A COHORT REBUILD: tissues are auto-discovered under
# --cohort-root (anything with highQual_samples_<Tissue>.txt +
# fraser_objs/cohort_results_all.tsv.gz + a *.rds is a tissue; anything else,
# e.g. an empty tissue dir or the jobs/archive dirs, is silently skipped), and
# each sample's project folder is auto-discovered under --repo-root rather
# than assumed. So re-pointing --cohort-root at a later cohort snapshot (e.g.
# .../control_cohort/Oct2026) needs no other edits. Already-annotated samples
# are skipped unless FORCE=1 (gate = the expected output files already exist).
#
# ⭐ GOTCHA: at least one tissue's OUTRIDER rds has been seen with sampleID set
# to the FULL source path rather than the bare sample name for some rows (e.g.
# Sept2026 Blood: 134/175 are "/hpf/.../ComplexCare/21-8881", not "21-8881";
# Fibroblasts/LCL in the same cohort snapshot did NOT have this). parse_outrider.R
# -s does an exact match and hard-errors otherwise, so this script resolves the
# literal sampleID per tissue via resolve_outrider_sampleids.R (cached, rebuilt
# whenever the rds is newer than the cache) instead of assuming the highQual
# list's ID is what's in the table. FRASER's cohort_results_all.tsv.gz has
# always used the plain ID in spot checks, so no resolution is done there.
#
# Usage:
#   ./submit_parse_control_cohort.sh                       # every tissue, every sample
#   ./submit_parse_control_cohort.sh --tissue Blood         # just one tissue
#   ./submit_parse_control_cohort.sh --tissue Blood 24-13843 20-10360   # just these samples
#   ./submit_parse_control_cohort.sh --cohort-root /hpf/.../control_cohort/Oct2026
#   DRY=1 ./submit_parse_control_cohort.sh                  # print job scripts, submit nothing
#   FORCE=1 ./submit_parse_control_cohort.sh                # re-run even if outputs exist
set -uo pipefail

SCRIPTS=${SCRIPTS:-/hpf/largeprojects/tgnode/sandbox/yliang_analysis/SR_RNAseq_pipeline/scripts}
COHORT_ROOT=${COHORT_ROOT:-/hpf/largeprojects/tgnode/sandbox/rnaseq_repo/control_cohort/Sept2026}
REPO_ROOT=${REPO_ROOT:-/hpf/largeprojects/tgnode/sandbox/rnaseq_repo}
HPO_DIR=${HPO_DIR:-$REPO_ROOT/hpoterm}
PROJECT_DIRS=(TG_RNAseq ComplexCare new_GenomeClinic)

# Both parsers pull in OUTRIDER/FRASER/Bioconductor, so run inside the
# container. --vanilla is required, not optional: singularity mounts $HOME,
# and the host ~/.Renviron puts a host R_LIBS_USER ahead of the image's own
# library, which breaks library(OUTRIDER)/library(FRASER) with a libxml2
# error. --cleanenv alone does not fix it.
SIF=${SIF:-/hpf/largeprojects/tgnode/sandbox/yliang_analysis/tools/containers/sr_rnaseq.sif}

# Per-job resources. GUESSES, not measured: each job reads the tissue-wide
# OUTRIDER rds and the tissue-wide FRASER cohort_results_all.tsv.gz into R.
# The largest tissue seen so far (Sept2026 Blood, 175 samples) is a 1.2 GB
# gzip / ~44M rows / ~8-9 GB uncompressed FRASER table read via read.delim, so
# this is sized for that, not the smaller tissues. Right-size per tissue after
# the first batch with `seff <jobid>` / `sacct -j <jobid> -o MaxRSS,Elapsed`.
CPUS=${CPUS:-1}
MEM=${MEM:-16G}
TIME=${TIME:-2:00:00}

want=("$@")
tissue_filter=()
while [[ ${#want[@]} -gt 0 ]]; do
    case "${want[0]}" in
        --tissue)
            tissue_filter+=("${want[1]}")
            want=("${want[@]:2}")
            ;;
        --cohort-root)
            COHORT_ROOT="${want[1]}"
            want=("${want[@]:2}")
            ;;
        --repo-root)
            REPO_ROOT="${want[1]}"
            HPO_DIR=${HPO_DIR_OVERRIDE:-$REPO_ROOT/hpoterm}
            want=("${want[@]:2}")
            ;;
        *)
            break
            ;;
    esac
done
sample_filter=("${want[@]}")   # remaining positional args = sample IDs to restrict to

keep_tissue() {
    [[ ${#tissue_filter[@]} -eq 0 ]] && return 0
    local t=$1 f
    for f in "${tissue_filter[@]}"; do [[ "$t" == "$f" ]] && return 0; done
    return 1
}
keep_sample() {
    [[ ${#sample_filter[@]} -eq 0 ]] && return 0
    local s=$1 f
    for f in "${sample_filter[@]}"; do [[ "$s" == "$f" ]] && return 0; done
    return 1
}

[[ -d "$COHORT_ROOT" ]] || { echo "!! cohort root not found: $COHORT_ROOT" >&2; exit 1; }
[[ -d "$REPO_ROOT" ]]   || { echo "!! repo root not found: $REPO_ROOT" >&2; exit 1; }

jobids=()
n_skipped_done=0
n_unresolved=0
shopt -s nullglob

for tissue_dir in "$COHORT_ROOT"/*/; do
    tissue_dir="${tissue_dir%/}"
    tissue=$(basename "$tissue_dir")
    keep_tissue "$tissue" || continue

    samples_file="$tissue_dir/highQual_samples_${tissue}.txt"
    fraser_file="$tissue_dir/fraser_objs/cohort_results_all.tsv.gz"
    rds_candidates=("$tissue_dir"/*.rds)

    if [[ ! -f "$samples_file" || ! -f "$fraser_file" || ${#rds_candidates[@]} -eq 0 ]]; then
        [[ ${#tissue_filter[@]} -gt 0 ]] && \
            echo ">> WARN: $tissue is missing highQual_samples_${tissue}.txt, fraser_objs/cohort_results_all.tsv.gz, or a *.rds -- skipping"
        continue
    fi
    if [[ ${#rds_candidates[@]} -gt 1 ]]; then
        echo ">> WARN: $tissue has ${#rds_candidates[@]} *.rds files, using ${rds_candidates[0]}: ${rds_candidates[*]}"
    fi
    rds="${rds_candidates[0]}"

    echo "=== $tissue  (rds: $(basename "$rds"), samples: $(basename "$samples_file")) ==="

    # --- resolve OUTRIDER's literal sampleID per plain sample ID (cached) ---
    map_file="$tissue_dir/.outrider_sampleid_map.tsv"
    if [[ ! -f "$map_file" || "$rds" -nt "$map_file" ]]; then
        echo "   resolving OUTRIDER sampleIDs -> $map_file"
        singularity exec --cleanenv --bind /hpf/largeprojects "$SIF" \
            Rscript --vanilla "$SCRIPTS/resolve_outrider_sampleids.R" "$rds" "$map_file" \
            || { echo "!! failed to resolve OUTRIDER sampleIDs for $tissue -- skipping tissue" >&2; continue; }
    fi

    STAMP=$(date +%Y%m%d)
    LOGDIR="${tissue_dir}/jobs/parse_${STAMP}"
    mkdir -p "$LOGDIR"

    while IFS= read -r sample; do
        [[ -z "$sample" ]] && continue
        keep_sample "$sample" || continue

        # --- resolve the sample's project folder ---
        proj=""
        for p in "${PROJECT_DIRS[@]}"; do
            if [[ -d "$REPO_ROOT/$p/$sample" ]]; then
                if [[ -n "$proj" ]]; then
                    echo ">> WARN: $sample found in >1 project dir ($proj AND $p); using $proj"
                else
                    proj="$p"
                fi
            fi
        done
        if [[ -z "$proj" ]]; then
            echo ">> WARN: $sample ($tissue) not found under TG_RNAseq/ComplexCare/new_GenomeClinic -- skipping"
            n_unresolved=$((n_unresolved + 1))
            continue
        fi
        sample_dir="$REPO_ROOT/$proj/$sample"

        # --- HPO file (optional; leave HPO columns NA if absent) ---
        hpo="$HPO_DIR/${sample}_hpo.txt"
        hpo_arg=""
        [[ -f "$hpo" ]] && hpo_arg="    -p '$hpo'"

        # --- literal OUTRIDER sampleID (may differ from $sample; see header) ---
        outrider_sid=$(awk -F'\t' -v s="$sample" '$1==s{print $2; found=1} END{if(!found) exit 1}' "$map_file")
        have_outrider=$?

        out_outrider="$sample_dir/${sample}_${tissue}_OUTRIDER_results_annotated.xlsx"
        out_fraser_base="$sample_dir/${sample}_${tissue}_FRASER_results_annotated"
        out_fraser="${out_fraser_base}.tsv"          # parse_fraser.R writes <base>.tsv.gz + <base>.filtered.tsv

        # --- idempotency: skip if already parsed (unless FORCE=1) ---
        if [[ -z "${FORCE:-}" && -f "$out_outrider" && -f "${out_fraser_base}.filtered.tsv" ]]; then
            n_skipped_done=$((n_skipped_done + 1))
            continue
        fi

        outrider_note=""
        [[ $have_outrider -ne 0 ]] && outrider_note=", no OUTRIDER match"
        echo ">> $sample  ($tissue, project: $proj${hpo_arg:+, hpo: yes}${outrider_note})"

        body=$(cat <<EOF
#!/bin/bash
#SBATCH --job-name=parseCC_${tissue}_${sample}
#SBATCH --cpus-per-task=${CPUS}
#SBATCH --mem=${MEM}
#SBATCH --time=${TIME}
#SBATCH --constraint=AlmaLinux9
#SBATCH --output=${LOGDIR}/parse_${sample}_%j.log
#SBATCH --error=${LOGDIR}/parse_${sample}_%j.err
set -uo pipefail
echo "job \${SLURM_JOB_ID:-NA} on \$(hostname) -- $tissue/$sample -- \$(date)"

EOF
)
        if [[ $have_outrider -eq 0 ]]; then
            body+=$(cat <<EOF

cmd_outrider=(
    singularity exec --cleanenv --bind /hpf/largeprojects '${SIF}'
    Rscript --vanilla '${SCRIPTS}/parse_outrider.R'
    -r '${rds}'
    -s '${outrider_sid}'
${hpo_arg}
    -o '${out_outrider}'
)
printf 'running:'; printf ' %q' "\${cmd_outrider[@]}"; printf '\n'
"\${cmd_outrider[@]}" || echo ">> OUTRIDER FAILED: $sample"
EOF
)
        else
            body+=$(cat <<EOF

echo ">> SKIPPING OUTRIDER: $sample has no sampleID match in $(basename "$rds")"
EOF
)
        fi

        body+=$(cat <<EOF

cmd_fraser=(
    singularity exec --cleanenv --bind /hpf/largeprojects '${SIF}'
    Rscript --vanilla '${SCRIPTS}/parse_fraser.R'
    -r '${fraser_file}'
    -s '${sample}'
${hpo_arg}
    -o '${out_fraser}'
)
printf 'running:'; printf ' %q' "\${cmd_fraser[@]}"; printf '\n'
"\${cmd_fraser[@]}" || echo ">> FRASER FAILED: $sample"

echo "done -- \$(date)"
EOF
)

        if [[ -n "${DRY:-}" ]]; then
            printf '%s\n---\n' "$body"
            continue
        fi
        jid=$(sbatch --parsable <<<"$body") \
            && { echo "   submitted job $jid"; jobids+=("$jid"); } \
            || echo ">> SUBMIT FAILED: $sample"
    done < "$samples_file"
done

echo
if [[ -n "${DRY:-}" ]]; then
    echo ">> dry run, nothing submitted"
elif [[ ${#jobids[@]} -gt 0 ]]; then
    echo ">> submitted ${#jobids[@]} job(s); skipped $n_skipped_done already-parsed, $n_unresolved unresolved (no project folder)"
    echo "   watch:  squeue -u \$USER"
    echo "   logs :  <tissue>/jobs/parse_${STAMP:-<date>}/parse_<sample>_<jobid>.{log,err}"
    echo "   after:  sacct -j $(IFS=,; echo "${jobids[*]}") -o JobID,JobName%30,State,Elapsed,MaxRSS"
else
    echo ">> nothing submitted (skipped $n_skipped_done already-parsed, $n_unresolved unresolved)"
fi
