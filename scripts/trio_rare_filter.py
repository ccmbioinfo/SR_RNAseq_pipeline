#!/usr/bin/env python3
"""
trio_rare_filter.py — rare-disease inheritance filtering on a VEP + vcfanno
annotated trio VCF (Care4Rare / crg2 style).

Writes one VCF per inheritance model into an output directory:
    de_novo.vcf.gz         proband has ALT, both parents hom-ref
    recessive_hom.vcf.gz   proband hom-ALT, both parents het (carriers)
    x_linked.vcf.gz        chrX: proband has ALT, mother het carrier, father hom-ref
    comp_het.vcf.gz        genes with >=1 maternally-inherited AND >=1 paternally-
                           inherited het in the proband (candidate compound-het)

Rarity = population AF below --af (gnomAD FAF95 if the VCF carries it, else
gnomAD popmax) plus CoLoRSdb AF when present; variants MISSING from those
databases (i.e. novel) are KEPT. Multiallelic sites are split first.

As a by-product it writes model_variants.tsv, a flat per-variant table for the
downstream RNA-outlier join. Besides the four inheritance models, that table also
carries a 'high_impact' rescue set — proband variants that are HIGH/MODERATE VEP
impact or ClinVar Pathogenic/Likely_pathogenic, kept regardless of inheritance
model (still within the rare set) — and a trailing 'inheritance' column giving the
pattern inferred from the trio genotypes for every row.

Heavy lifting is done by bcftools; this script orchestrates it and does the
gene-level compound-het pairing in Python.

Example
-------
    python3 trio_rare_filter.py -v DSK007.coding.vep.vcfanno.vcf.gz
    python3 trio_rare_filter.py -v in.vcf.gz -p KID -m MUM -f DAD -a 0.001 -o out/

Requires: bcftools (with the split-vep plugin) and tabix on PATH.
"""
import argparse
import os
import re
import shutil
import subprocess
import sys
from collections import defaultdict

# VEP consequence severity, most -> least severe (rank = index). Used to collapse
# the per-transcript split-vep rows in model_variants.tsv down to the single
# worst-consequence row per (model, variant, gene).
VEP_ORDER = [
    "transcript_ablation", "splice_acceptor_variant", "splice_donor_variant",
    "stop_gained", "frameshift_variant", "stop_lost", "start_lost",
    "transcript_amplification", "feature_elongation", "feature_truncation",
    "inframe_insertion", "inframe_deletion", "missense_variant",
    "protein_altering_variant", "splice_donor_5th_base_variant",
    "splice_region_variant", "splice_donor_region_variant",
    "splice_polypyrimidine_tract_variant", "incomplete_terminal_codon_variant",
    "start_retained_variant", "stop_retained_variant", "synonymous_variant",
    "coding_sequence_variant", "mature_miRNA_variant", "5_prime_UTR_variant",
    "3_prime_UTR_variant", "non_coding_transcript_exon_variant", "intron_variant",
    "NMD_transcript_variant", "non_coding_transcript_variant",
    "coding_transcript_variant", "upstream_gene_variant", "downstream_gene_variant",
    "TFBS_ablation", "TFBS_amplification", "TF_binding_site_variant",
    "regulatory_region_ablation", "regulatory_region_amplification",
    "regulatory_region_variant", "intergenic_variant", "sequence_variant",
]
VEP_RANK = {c: i for i, c in enumerate(VEP_ORDER)}


def csq_rank(conseq):
    """Rank of the most severe term in an '&'-joined consequence string."""
    return min((VEP_RANK.get(t, len(VEP_ORDER)) for t in conseq.split("&")),
               default=len(VEP_ORDER))


def geno(gtstr):
    """Collapse a VCF GT string to RR / RA / AA / NA (hom-ref / het / hom-alt /
    missing). Post-`norm -m-any` sites are biallelic, so any non-'0' allele is the
    ALT; a haploid call (e.g. a male chrX '1') that carries ALT collapses to AA."""
    alleles = [a for a in re.split(r"[/|]", gtstr) if a not in (".", "")]
    if not alleles:
        return "NA"
    has_alt = any(a != "0" for a in alleles)
    has_ref = "0" in alleles
    if has_alt and has_ref:
        return "RA"
    if has_alt:
        return "AA"
    return "RR"


def classify_inheritance(chrom, gp, gm, gd):
    """Infer an inheritance pattern for a proband variant from the trio genotypes
    (proband/mother/father, each RR/RA/AA/NA as returned by geno()). One of:
    de_novo, maternal, paternal, biparental_het, recessive_hom, x_linked, hom_alt,
    uncertain. In singleton mode (no parents) the caller emits '.' instead."""
    if gp in ("NA", "RR"):
        return "uncertain"           # proband doesn't confidently carry the ALT
    if gm == "NA" or gd == "NA":
        return "uncertain"           # a parent is uncalled -> can't segregate
    isX = chrom in ("chrX", "X")
    mom_alt = gm in ("RA", "AA")
    dad_alt = gd in ("RA", "AA")
    if gp == "AA":                   # proband hom-ALT (or male-X hemizygous)
        if isX and mom_alt and not dad_alt:
            return "x_linked"
        if gm == "RA" and gd == "RA":
            return "recessive_hom"
        return "hom_alt"
    # gp == "RA": proband het
    if not mom_alt and not dad_alt:
        return "de_novo"
    if isX and mom_alt and not dad_alt:
        return "x_linked"
    if mom_alt and not dad_alt:
        return "maternal"
    if dad_alt and not mom_alt:
        return "paternal"
    return "biparental_het"          # both parents carry -> ambiguous (comp-het side)


def classify_inheritance_duo(chrom, gp, g_present, present_is_mother):
    """Duo counterpart of classify_inheritance() for a proband + ONE parent
    (`g_present`, RR/RA/AA/NA as returned by geno()). Labels are deliberately
    non-committal where the missing parent's genotype would be needed to
    disambiguate (e.g. a proband het not carried by the present parent could be
    de novo OR inherited from the missing parent -- both read the same here).
    One of: denovo_or_missing_parent, maternal_transmitted, paternal_transmitted,
    x_linked_candidate, hom_alt_untransmitted_present_parent,
    <parent>_and_unknown_parent, uncertain."""
    if gp in ("NA", "RR") or g_present == "NA":
        return "uncertain"
    isX = chrom in ("chrX", "X")
    present_alt = g_present in ("RA", "AA")
    parent_label = "maternal" if present_is_mother else "paternal"
    if gp == "AA":                   # proband hom-ALT (or hemizygous)
        if isX and present_is_mother and present_alt:
            return "x_linked_candidate"
        if present_alt:
            # present parent is a carrier; the 2nd allele's origin needs the
            # missing parent to confirm biallelic (recessive) inheritance
            return f"{parent_label}_and_unknown_parent"
        return "hom_alt_untransmitted_present_parent"
    # gp == "RA": proband het
    if present_alt:
        return f"{parent_label}_transmitted"
    return "denovo_or_missing_parent"


def die(msg):
    sys.exit(f"ERROR: {msg}")


def run(cmd, capture=False):
    """Run a shell command; abort on failure. Return stdout if capture=True."""
    r = subprocess.run(cmd, shell=True, text=True,
                       stdout=subprocess.PIPE if capture else None)
    if r.returncode != 0:
        die(f"command failed ({r.returncode}): {cmd}")
    return r.stdout if capture else None


def samples(vcf):
    return run(f"bcftools query -l {vcf}", capture=True).split()


def autodetect(names, suffix):
    """Pick the sample whose name ends in _0<suffix> or _<suffix>."""
    for n in names:
        if n.endswith(f"_0{suffix}") or n.endswith(f"_{suffix}"):
            return n
    return None


def parse_ped(path):
    """Parse a 6+ column PED/peddy pedigree file (family, sample, dad, mom, sex,
    phenotype[, ...extra peddy columns ignored]). '0' parent ids -> None.
    Returns {sample_id: (family_id, dad_id, mom_id, sex, phenotype)}."""
    ped = {}
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line or line.startswith("#"):
                continue
            f = line.split("\t") if "\t" in line else line.split()
            if len(f) < 6:
                continue
            fam, iid, dad, mom, sex, phen = f[0], f[1], f[2], f[3], f[4], f[5]
            ped[iid] = (fam, dad if dad != "0" else None,
                       mom if mom != "0" else None, sex, phen)
    return ped


def parse_vcf_pedigree(vcf):
    """Parse DRAGEN-style '##PEDIGREE=<Child=X,Mother=Y,Father=Z>' header lines as
    a pedigree source alternative to --ped. Sex/phenotype are unknown from this
    source. Returns the same shape as parse_ped()."""
    hdr = run(f"bcftools view -h {vcf}", capture=True)
    ped = {}
    for line in hdr.splitlines():
        if not line.startswith("##PEDIGREE=<"):
            continue
        body = line[len("##PEDIGREE=<"):].rstrip(">")
        kv = dict(item.split("=", 1) for item in body.split(",") if "=" in item)
        child = kv.get("Child")
        if child:
            ped[child] = (None, kv.get("Father"), kv.get("Mother"), None, None)
    return ped


def resolve_family(names, ped, proband):
    """Resolve mother/father/siblings for `proband` from a parsed pedigree dict,
    restricted to samples actually present in the VCF (`names`). A sibling is any
    other pedigree entry sharing BOTH parent ids with the proband (both non-null).
    Returns (mother_id_or_None, father_id_or_None, [sibling_ids])."""
    info = ped.get(proband)
    if not info:
        return None, None, []
    _, dad_id, mom_id, _, _ = info
    mother = mom_id if mom_id in names else None
    father = dad_id if dad_id in names else None
    siblings = []
    if dad_id and mom_id:
        for s in names:
            if s == proband:
                continue
            sinfo = ped.get(s)
            if not sinfo:
                continue
            _, sdad, smom, _, _ = sinfo
            if sdad == dad_id and smom == mom_id:
                siblings.append(s)
    return mother, father, siblings


def chrom_x(vcf):
    hdr = run(f"bcftools view -h {vcf}", capture=True)
    for line in hdr.splitlines():
        if line.startswith("##contig"):
            # ##contig=<ID=chrX,...>
            cid = line.split("ID=", 1)[1].split(",")[0].split(">")[0]
            if cid in ("chrX", "X"):
                return cid
    return "chrX"


def info_ids(vcf):
    """Set of INFO field IDs declared in the VCF header."""
    hdr = run(f"bcftools view -h {vcf}", capture=True)
    return {line.split("ID=", 1)[1].split(",")[0]
            for line in hdr.splitlines() if line.startswith("##INFO=<")}


def format_ids(vcf):
    """Set of FORMAT (per-sample) field IDs declared in the VCF header."""
    hdr = run(f"bcftools view -h {vcf}", capture=True)
    return {line.split("ID=", 1)[1].split(",")[0]
            for line in hdr.splitlines() if line.startswith("##FORMAT=<")}


def csq_subfields(vcf):
    """Names of the VEP CSQ subfields exposed by +split-vep (empty if none)."""
    r = subprocess.run(f"bcftools +split-vep -l {vcf}", shell=True, text=True,
                       stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)
    subs = set()
    for ln in (r.stdout or "").splitlines():
        parts = ln.replace("\t", " ").split()
        if len(parts) >= 2:            # "<idx> <Name> [Type ...]"
            subs.add(parts[1])
    return subs


def write_variant_tsv(outdir, pro, mom, dad, sources, af_field, available,
                      siblings=None, ped=None):
    """Write model_variants.tsv: a flat, VEP-expanded (split-vep -d), model-tagged
    per-variant table for downstream joins (e.g. parse_outrider.R). One row per
    (model, variant, gene) tuple, keeping the worst-consequence transcript. The
    'gt' column is the proband genotype; a trailing 'inheritance' column is the
    pattern inferred from the trio genotypes (mom/dad), the duo genotype (one
    parent present), or '.' in singleton mode.
    Columns:
        model chrom pos ref alt gt gene consequence impact hgvsc hgvsp
        popmax cadd spliceai clinvar revel inheritance
    Annotation fields that the VCF lacks are emitted as '.' so the schema is
    stable across cohorts (TG vs DECODER field-name differences).

    When `siblings` is given (quad mode, or duo + siblings) every row gains a
    'sib_<ID>_gt' column per sibling (raw GT) plus a trailing 'sib_concordant'
    column: a ';'-joined list of siblings that BOTH share the proband's
    carrier/hom-alt state AND are marked affected (phenotype '2') in `ped`.
    This is annotation only -- it never filters or drops a row -- so the base
    (no-sibling) schema above is unchanged for existing consumers."""
    def pick(cands):
        for c in cands:
            if c and c in available:
                return f"%{c}"
        return "."
    cols = {
        "popmax":   ["gnomad_af_grpmax", "gnomad_af_popmax", af_field],
        "cadd":     ["CADD_phred", "CADD_PHRED", "cadd_phred"],
        "spliceai": ["spliceai_score", "SpliceAI", "spliceai"],
        # clinvar_pathogenic (direct CLNSIG concat) first: the C4R/DECODER VCFs
        # declare a derived clinvar_sig that is usually empty, so preferring it
        # blanked the column even when clinvar_pathogenic held the P/LP call.
        "clinvar":  ["clinvar_pathogenic", "clinvar_sig", "CLNSIG", "clinvar_status"],
        "revel":    ["REVEL_score", "REVEL", "revel"],
    }
    tail = r"\t".join(pick(cols[k])
                      for k in ("popmax", "cadd", "spliceai", "clinvar", "revel"))
    # HGVS nomenclature (transcript-specific; split-vep -d keeps it aligned to the
    # per-transcript consequence). '.' when the VCF's VEP config didn't emit it.
    hgvs = pick(["HGVSc"]) + r"\t" + pick(["HGVSp"])
    # [%GT ] emits every sample's GT (space-separated) in one field so we can read
    # the trio and infer inheritance; the proband GT is picked out for the output
    # 'gt' column. -d -> one row per transcript so the join can take the
    # gene-specific worst consequence.
    fmt = (r"%CHROM\t%POS\t%REF\t%ALT\t[%GT ]\t%Gene\t%Consequence\t%IMPACT\t"
           + hgvs + r"\t" + tail + r"\n")
    trio = mom is not None and dad is not None
    duo = (mom is not None) != (dad is not None)   # exactly one parent present
    present_is_mother = mom is not None if duo else None
    ip = _SIDX[pro]
    im = _SIDX[mom] if mom is not None else None
    idd = _SIDX[dad] if dad is not None else None
    sib_info = [(s, _SIDX[s]) for s in (siblings or []) if s in _SIDX]
    ped = ped or {}
    # split-vep -d emits one row per transcript; keep only the worst-consequence
    # row per (model, variant, gene) so the table isn't inflated with redundant
    # per-transcript lines (the downstream join uses exactly this granularity).
    best = {}   # (model, chrom, pos, ref, alt, gene) -> (rank, [fields...])
    order = []  # preserve first-seen order for stable output
    for label, vcf in sources:
        if not os.path.exists(vcf):
            continue
        cmd = f"bcftools +split-vep -d -f '{fmt}' {vcf}"
        p = subprocess.Popen(cmd, shell=True, text=True,
                             stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        for line in p.stdout:
            f = line.rstrip("\n").split("\t")
            if len(f) < 15:
                continue
            gts = f[4].split()                       # per-sample GTs, VCF order
            gt_pro = gts[ip] if ip < len(gts) else "."
            f[4] = gt_pro                            # 'gt' column = proband only
            gp = geno(gt_pro)
            if trio:
                gt_mom = gts[im] if im < len(gts) else "."
                gt_dad = gts[idd] if idd < len(gts) else "."
                inh = classify_inheritance(f[0], gp, geno(gt_mom), geno(gt_dad))
            elif duo:
                pidx = im if present_is_mother else idd
                gt_present = gts[pidx] if pidx < len(gts) else "."
                inh = classify_inheritance_duo(f[0], gp, geno(gt_present),
                                               present_is_mother)
            else:
                inh = "."
            sib_fields, concordant = [], []
            for sname, sidx in sib_info:
                gt_sib = gts[sidx] if sidx < len(gts) else "."
                sib_fields.append(gt_sib)
                if gp in ("RA", "AA"):
                    gsib = geno(gt_sib)
                    same_carrier = ((gp == "AA" and gsib == "AA") or
                                   (gp == "RA" and gsib in ("RA", "AA")))
                    phen = ped.get(sname, (None, None, None, None, None))[4]
                    if same_carrier and phen == "2":
                        concordant.append(sname)
            row = [label] + f + [inh]
            if sib_info:
                row += sib_fields + [";".join(concordant) or "."]
            key = (label, f[0], f[1], f[2], f[3], f[5])   # model,chrom,pos,ref,alt,gene
            r = csq_rank(f[6])
            cur = best.get(key)
            if cur is None:
                order.append(key)
                best[key] = (r, row)
            elif r < cur[0]:
                best[key] = (r, row)
        _, err = p.communicate()
        if p.returncode != 0:
            sys.stderr.write(err)
            die(f"split-vep failed on {vcf}")
    # Normalise the gene namespace. On a VEP --merged VCF, split-vep -d emits BOTH
    # an Ensembl row (Gene=ENSG..., Feature=ENST...) and a RefSeq row
    # (Gene=<Entrez integer>, Feature=NM_/NR_...) per transcript, so %Gene mixes
    # ENSG and Entrez ids. Every downstream join (OUTRIDER/FRASER geneID,
    # gene_coords) is Ensembl-based, so a numeric Entrez gene never matches and is
    # a redundant duplicate whenever the same variant also carries an ENSG row.
    # Drop those, but keep a numeric gene when it's the ONLY annotation the variant
    # has (RefSeq-only genes absent from Ensembl) so no variant loses its gene.
    ensg_variants = {key[1:5] for key in order if key[5].startswith("ENSG")}
    dropped = 0
    kept_order = []
    for key in order:
        if key[5].isdigit() and key[1:5] in ensg_variants:
            dropped += 1
            continue
        kept_order.append(key)
    if dropped:
        print(f"  gene namespace  dropped {dropped} redundant RefSeq/Entrez rows "
              f"(kept Ensembl equivalents)")
    order = kept_order
    tsv = os.path.join(outdir, "model_variants.tsv")
    header = ("model\tchrom\tpos\tref\talt\tgt\tgene\tconsequence\timpact\t"
             "hgvsc\thgvsp\tpopmax\tcadd\tspliceai\tclinvar\trevel\tinheritance")
    if sib_info:
        header += "\t" + "\t".join(f"sib_{s}_gt" for s, _ in sib_info) + "\tsib_concordant"
    with open(tsv, "w") as fh:
        fh.write(header + "\n")
        for key in order:
            fh.write("\t".join(best[key][1]) + "\n")
    print(f"  variant tsv    {tsv}  ({len(order)} rows)")
    return tsv


def build_high_impact(work, pro, info, outdir):
    """Build a temp VCF of proband-carried rare variants worth keeping regardless
    of inheritance model: VEP IMPACT HIGH or MODERATE, OR ClinVar Pathogenic/
    Likely_pathogenic. These feed model_variants.tsv tagged 'high_impact' (each
    row's real segregation shows up in the inheritance column). Operates on the
    already-rare, PASS working file, so no AF/rarity bypass. Returns the temp path
    (caller removes it after the TSV is written)."""
    # ClinVar P/LP across whichever significance field(s) the VCF declares. '~' is
    # a case-sensitive regex, so "Pathogenic"/"Likely_pathogenic" match P/LP and
    # its combos (e.g. "Pathogenic,_risk_factor") but NOT the lowercase
    # "Conflicting_..._pathogenicity". Several DECODER VCFs declare an (empty)
    # clinvar_sig alongside a populated clinvar_pathogenic, so OR across all present.
    clin_fields = [c for c in ("clinvar_sig", "clinvar_pathogenic", "CLNSIG")
                   if c in info]
    keep = 'IMPACT="HIGH" || IMPACT="MODERATE"'
    for c in clin_fields:
        keep += f' || {c}~"Pathogenic" || {c}~"Likely_pathogenic"'
    out = os.path.join(outdir, ".high_impact.vcf.gz")
    # split-vep -s worst -c IMPACT materialises the worst-consequence IMPACT into
    # INFO so the second view can OR it with the site-level ClinVar field(s); then
    # keep only variants the proband actually carries.
    run(f'bcftools +split-vep {work} -s worst -c IMPACT -Ou '
        f'| bcftools view -i \'({keep}) && {gt(pro, "alt")}\' -Oz -o {out}')
    run(f"tabix -f -p vcf {out}")
    return out


_SIDX = {}  # sample name -> 0-based column index; filled in main()


def gt(sample, code):
    """A bcftools per-sample genotype test by 0-based column index,
    e.g. gt('KID','alt') -> GT[0]="alt". bcftools GT[...] subscripts take an
    index, not a sample name, so we translate via _SIDX."""
    return f'GT[{_SIDX[sample]}]="{code}"'


def ad_alt(sample, op, val):
    """A bcftools alt-allele-depth test for one sample. AD is Number=R
    (ref,alt), so allele subscript 1 is the ALT depth:
    ad_alt('KID','>=',3) -> FMT/AD[0:1]>=3."""
    return f'FMT/AD[{_SIDX[sample]}:1]{op}{val}'


def count(vcf):
    out = run(f"bcftools view -H {vcf} | wc -l", capture=True)
    return int(out.strip())


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("-v", "--vcf", required=True, help="annotated trio VCF (.vcf.gz)")
    ap.add_argument("-o", "--outdir", help="output dir (default: <vcf_dir>/trio_filter)")
    ap.add_argument("-p", "--proband", help="proband sample id (else auto *_03)")
    ap.add_argument("-m", "--mother", help="mother sample id (else auto *_02)")
    ap.add_argument("-f", "--father", help="father sample id (else auto *_01)")
    ap.add_argument("--ped", help="PED/peddy pedigree file (family, sample, dad, "
                    "mom, sex, phenotype) to resolve mother/father/siblings for "
                    "--proband. Enables duo (one parent) and quad (parents + "
                    "sibling(s)) modes. Overrides the *_01/*_02/*_03 suffix "
                    "autodetect for role resolution; -m/-f still override the "
                    "pedigree.")
    ap.add_argument("--from-vcf-pedigree", action="store_true",
                    help="also resolve pedigree from the VCF's own "
                         "'##PEDIGREE=<Child=,Mother=,Father=>' header lines "
                         "(e.g. DRAGEN joint-genotyped VCFs); --ped takes "
                         "priority for any sample present in both")
    ap.add_argument("--siblings", help="comma-separated sibling sample ids; "
                    "overrides pedigree-derived siblings")
    ap.add_argument("--no-siblings", action="store_true",
                    help="ignore pedigree-derived siblings (force plain "
                         "trio/duo behaviour even if the VCF has more samples)")
    ap.add_argument("-a", "--af", type=float, default=0.01,
                    help="population-AF cutoff; novel variants kept (default 0.01)")
    ap.add_argument("--faf-field", default="gnomad_fafmax_faf95_max",
                    help="preferred INFO field for the population-AF filter "
                         "(FAF95). Used if present in the header; otherwise the "
                         "script falls back to --popmax-field "
                         "(default: gnomad_fafmax_faf95_max)")
    ap.add_argument("--popmax-field", default="gnomad_af_popmax",
                    help="fallback population-AF INFO field when --faf-field is "
                         "absent (default: gnomad_af_popmax)")
    ap.add_argument("--min-alt", type=int, default=3,
                    help="min proband ALT-read depth. In trio mode this gates "
                         "de_novo calls (and both parents must have 0 ALT reads, "
                         "guarding against parental dropout); in --singleton mode "
                         "it gates every kept variant. 0 disables (default 3)")
    ap.add_argument("--singleton", action="store_true",
                    help="proband-only VCF (no parents): skip the inheritance "
                         "models and emit singleton.vcf.gz = rare PASS variants "
                         "the proband carries, passing basic genotype quality "
                         "(--min-alt, --min-gq). Auto-enabled for a 1-sample VCF")
    ap.add_argument("--min-gq", type=int, default=20,
                    help="min proband genotype quality (FORMAT/GQ) in --singleton "
                         "mode; applied only if the VCF carries GQ. 0 disables "
                         "(default 20)")
    ap.add_argument("-r", "--ref", help="reference FASTA (enables left-align in norm)")
    ap.add_argument("-i", "--impact", default="all",
                    help="VEP IMPACT levels to keep (comma-sep), e.g. "
                         "'HIGH,MODERATE'; 'all' disables the coding-impact filter "
                         "(default: all). A value other than 'all' requires a CSQ "
                         "field.")
    ap.add_argument("--gene-field", default=None,
                    help="INFO field to group compound-het by when the VCF has "
                         "no CSQ (e.g. an SnpEff/gene-symbol tag). With CSQ the "
                         "VEP SYMBOL is used automatically; without CSQ and "
                         "without this, comp_het is skipped")
    ap.add_argument("--no-variant-tsv", action="store_true",
                    help="skip writing model_variants.tsv, the VEP-expanded, "
                         "model-tagged per-variant table that parse_outrider.R "
                         "consumes (written by default when the VCF has CSQ)")
    o = ap.parse_args()

    for tool in ("bcftools", "tabix"):
        if not shutil.which(tool):
            die(f"{tool} not found on PATH")
    if not os.path.isfile(o.vcf):
        die(f"VCF not found: {o.vcf}")

    info = info_ids(o.vcf)
    has_csq = "CSQ" in info
    # split-vep (VEP IMPACT filter + gene SYMBOL for comp-het) is only needed
    # when the VCF carries CSQ; only then must the plugin be loadable.
    if has_csq:
        # Functional test: 'split-vep -l' lists CSQ subfields and exits 0 when
        # the plugin loads. Don't probe with '-h': plugins print help but exit
        # non-zero by convention, so an exit-code check there always false-fails.
        if subprocess.run(f"bcftools +split-vep -l {o.vcf}", shell=True,
                          stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL).returncode != 0:
            die("VCF has a CSQ field but the bcftools split-vep plugin is not "
                "loadable (set BCFTOOLS_PLUGINS)")

    names = samples(o.vcf)
    _SIDX.update({n: i for i, n in enumerate(names)})

    # --- pedigree sources (optional): --ped takes priority over --from-vcf-pedigree
    #     for any sample present in both. Neither given -> old suffix autodetect,
    #     byte-identical to pre-pedigree behaviour. ---
    ped = {}
    if o.from_vcf_pedigree:
        ped.update(parse_vcf_pedigree(o.vcf))
    if o.ped:
        ped.update(parse_ped(o.ped))

    pro = o.proband or autodetect(names, 3) or (names[0] if len(names) == 1 else None)

    siblings = []
    if ped and pro in ped and not o.singleton and len(names) > 1:
        ped_mom, ped_dad, ped_sibs = resolve_family(names, ped, pro)
        mom = o.mother or ped_mom
        dad = o.father or ped_dad
        siblings = ped_sibs
    else:
        mom = o.mother or autodetect(names, 2)
        dad = o.father or autodetect(names, 1)

    if o.siblings:
        siblings = [s.strip() for s in o.siblings.split(",") if s.strip()]
    if o.no_siblings:
        siblings = []
    # never double-count a resolved parent (or the proband) as a sibling
    siblings = [s for s in siblings if s in names and s not in (pro, mom, dad)]

    # a 1-sample VCF (or --singleton) can't support any segregation model: fall
    # back to a proband-only AF + basic-quality filter.
    singleton = o.singleton or len(names) == 1
    if singleton:
        if not pro:
            die(f"could not determine proband; pass -p. Samples: {names}")
        if pro not in names:
            die(f"proband '{pro}' not in VCF. Samples: {names}")
        mom = dad = None
        siblings = []
    else:
        if not pro:
            die(f"could not determine proband; pass -p or --ped. Samples: {names}")
        if pro not in names:
            die(f"proband '{pro}' not in VCF. Samples: {names}")
        if not mom and not dad:
            die("could not determine mother or father; pass -m/-f or --ped/"
                "--from-vcf-pedigree (or --singleton for a proband-only VCF). "
                f"Samples: {names}")
        for role, val in (("mother", mom), ("father", dad)):
            if val and val not in names:
                die(f"{role} '{val}' not in VCF. Samples: {names}")

    mode = ("singleton" if singleton else
            "duo" if not (mom and dad) else
            "quad" if siblings else
            "trio")

    outdir = o.outdir or os.path.join(os.path.dirname(os.path.abspath(o.vcf)), "trio_filter")
    os.makedirs(outdir, exist_ok=True)
    xchr = chrom_x(o.vcf)

    impacts = [x.strip() for x in o.impact.split(",") if x.strip()]
    keep_impact = not (len(impacts) == 1 and impacts[0].lower() == "all")
    if keep_impact and not has_csq:
        die("--impact filtering needs a VEP CSQ field, which this VCF lacks; "
            "re-run with --impact all (or VEP-annotate the VCF first)")

    # --- pick the population-AF field: FAF95 preferred, then a popmax/grpmax
    #     fallback (gnomAD v4 renamed popmax -> grpmax). Add extra cohort-AF
    #     fields only when present, so the script also works on VCFs that lack
    #     them (e.g. short-read sets without CoLoRSdb). ---
    af_candidates, seen = [], set()
    for c in (o.faf_field, o.popmax_field, "gnomad_af_grpmax", "gnomad_af_popmax"):
        if c and c not in seen:
            seen.add(c)
            af_candidates.append(c)
    af_field = next((f for f in af_candidates if f in info), None)
    if not af_field:
        die("no population-AF field found (tried: " + ", ".join(af_candidates)
            + "); name one with --faf-field/--popmax-field")
    extra_af = [f for f in ("CoLoRSdb_AF",) if f in info]
    af_fields = [af_field] + extra_af
    excl = " || ".join(f"{f}>={o.af}" for f in af_fields)

    gene_field = ("SYMBOL" if has_csq
                  else o.gene_field if (o.gene_field and o.gene_field in info)
                  else None)

    has_gq = "GQ" in format_ids(o.vcf)

    # fields available for the downstream variant TSV = CSQ subfields + INFO tags
    variant_tsv = has_csq and not o.no_variant_tsv
    available = (csq_subfields(o.vcf) | info) if variant_tsv else set()

    print("== trio_rare_filter ==")
    print(f"  VCF     : {o.vcf}")
    if singleton:
        print(f"  mode    : SINGLETON (proband-only; no inheritance models)")
        print(f"  proband : {pro}")
    else:
        print(f"  mode    : {mode.upper()}")
        print(f"  proband : {pro}   mother: {mom or '(absent)'}   "
              f"father: {dad or '(absent)'}"
              + (f"   siblings: {', '.join(siblings)}" if siblings else ""))
        if mode == "duo":
            present_parent = mom or dad
            missing = "father" if mom else "mother"
            print(f"  duo     : only {'mother' if mom else 'father'} present "
                  f"({present_parent}); segregation from the missing {missing} "
                  f"cannot be confirmed -> de_novo/comp_het are NOT computed as "
                  f"in trio mode (see 'candidate' model tags below)")
    print(f"  AF < {o.af} ({', '.join(af_fields)}; novel kept)   chrX: {xchr}")
    print(f"  CSQ     : {'yes' if has_csq else 'no'}   "
          f"comp_het gene field: {gene_field or 'n/a (comp_het skipped)'}")
    print(f"  impact  : {'|'.join(impacts) if keep_impact else 'all (no impact filter)'}")
    if singleton:
        q = [f"ALT>={o.min_alt}" if o.min_alt > 0 else None,
             f"GQ>={o.min_gq}" if (o.min_gq > 0 and has_gq) else None]
        print(f"  quality : proband carries ALT"
              + (", " + ", ".join(x for x in q if x) if any(q) else " (no depth/GQ guard)"))
    elif mode in ("trio", "quad"):
        print(f"  de_novo : proband ALT>={o.min_alt}, parents ALT==0"
              if o.min_alt > 0 else "  de_novo : no depth guard (--min-alt 0)")
    if variant_tsv:
        print("  variants: model_variants.tsv (VEP-expanded, for parse_outrider.R)")
        if not singleton:
            print("            + 'high_impact' rescue (HIGH/MODERATE or ClinVar "
                  f"P/LP, any inheritance); {mode} inheritance column"
                  + (f"; sib_* columns for {', '.join(siblings)}" if siblings else ""))
    elif has_csq:
        print("  variants: skipped (--no-variant-tsv)")
    print(f"  outdir  : {outdir}")

    # --- 1. prefilter: split multiallelics, keep PASS, drop common variants,
    #        (optionally) keep only the requested VEP IMPACT levels, and carry
    #        one INFO/SYMBOL (worst consequence) for gene grouping. In trio mode
    #        this one working file feeds every model. In singleton mode there are
    #        no models, so the per-proband genotype-quality gate is folded in here
    #        and rare.vcf.gz IS the final output. ---
    norm = f"bcftools norm -m-any {'-f ' + o.ref if o.ref else ''} {o.vcf} -Ou"
    work = os.path.join(outdir, "rare.vcf.gz")
    stages = [norm, f'bcftools view -f PASS,. -e "{excl}" -Ou']
    if has_csq:
        # split-vep pulls the worst-consequence SYMBOL into INFO (for gene-level
        # comp-het grouping) and applies the VEP IMPACT filter when requested.
        sv_filter = ""
        if keep_impact:
            sv_filter = " -i '" + " || ".join(f'IMPACT="{x}"' for x in impacts) + "'"
        stages.append(f"bcftools +split-vep - -s worst -c SYMBOL{sv_filter} -Ou")
    if singleton:
        # proband-only: no segregation to test. Require the proband to carry the
        # ALT and pass basic genotype quality (ALT depth, GQ) right here, so
        # rare.vcf.gz needs no separate singleton output.
        sexpr = gt(pro, "alt")
        if o.min_alt > 0:
            sexpr += f' && {ad_alt(pro, ">=", o.min_alt)}'
        if o.min_gq > 0 and has_gq:
            sexpr += f' && FMT/GQ[{_SIDX[pro]}]>={o.min_gq}'
        stages.append(f"bcftools view -i '{sexpr}' -Ou")
    # last stage writes the compressed, indexed working file
    stages[-1] = stages[-1][:-len(" -Ou")] + f" -Oz -o {work}"
    run(" | ".join(stages))
    run(f"tabix -f -p vcf {work}")

    if singleton:
        if variant_tsv:
            write_variant_tsv(outdir, pro, None, None,
                              [("singleton", work)], af_field, available)
        print(f">> done. Single output: {work}")
        return

    print(">> applying models")

    def emit(name, expr, region=None):
        out = os.path.join(outdir, f"{name}.vcf.gz")
        reg = f"-r {region} " if region else ""
        run(f"bcftools view {reg}-i '{expr}' {work} -Oz -o {out}")
        run(f"tabix -f -p vcf {out}")
        print(f"  {name:<14} {count(out)} variants")

    def candidates(expr):
        """gene -> set of 'chrom\\tpos' for proband sites matching `expr`."""
        q = (f"bcftools view -i '{expr}' {work} | "
             + rf"bcftools query -f '%CHROM\t%POS\t%INFO/{gene_field}\n'")
        genes = defaultdict(set)
        for line in run(q, capture=True).splitlines():
            chrom, pos, gene = line.split("\t")
            if gene and gene != ".":
                genes[gene].add(f"{chrom}\t{pos}")
        return genes

    def write_gene_vcf(ch_vcf, genes_file, ch_genes, site_map, label):
        """Write comp-het-style outputs: a .genes list and a sites VCF."""
        with open(genes_file, "w") as fh:
            fh.write("\n".join(ch_genes) + ("\n" if ch_genes else ""))
        sites = sorted({s for g in ch_genes for s in site_map[g]},
                       key=lambda x: (x.split("\t")[0], int(x.split("\t")[1])))
        if sites:
            targets = os.path.join(outdir, ".ch_sites.tsv")
            with open(targets, "w") as fh:
                fh.write("\n".join(sites) + "\n")
            run(f"bcftools view -T {targets} {work} -Oz -o {ch_vcf}")
            os.remove(targets)
        else:
            run(f"bcftools view -h {work} -Oz -o {ch_vcf}")   # header only
        run(f"tabix -f -p vcf {ch_vcf}")
        print(f"  {label:<14} {count(ch_vcf)} variants in {len(ch_genes)} genes")

    if mode in ("trio", "quad"):
        denovo = f'{gt(pro,"alt")} && {gt(mom,"RR")} && {gt(dad,"RR")}'
        if o.min_alt > 0:
            # proband must have real ALT support and BOTH parents zero ALT
            # reads, so a low-coverage hom-ref call in a parent can't
            # masquerade as de novo
            denovo += (f' && {ad_alt(pro,">=",o.min_alt)}'
                       f' && {ad_alt(mom,"==",0)} && {ad_alt(dad,"==",0)}')
        emit("de_novo",       denovo)
        emit("recessive_hom", f'{gt(pro,"AA")} && {gt(mom,"RA")} && {gt(dad,"RA")}')
        emit("x_linked",      f'{gt(pro,"alt")} && {gt(mom,"RA")} && {gt(dad,"RR")}', xchr)
        model_names = ["de_novo", "recessive_hom", "x_linked"]

        # --- compound het: pair candidate hets by gene across the two
        #     parents. Needs a per-variant gene label (VEP SYMBOL, or
        #     --gene-field); without one we can't group, so comp_het is
        #     skipped. ---
        ch_vcf = os.path.join(outdir, "comp_het.vcf.gz")
        if gene_field:
            mat = candidates(f'{gt(pro,"RA")} && {gt(mom,"alt")} && {gt(dad,"RR")}')
            pat = candidates(f'{gt(pro,"RA")} && {gt(dad,"alt")} && {gt(mom,"RR")}')
            ch_genes = sorted(set(mat) & set(pat))   # gene het from BOTH parents
            merged = defaultdict(set)
            for g in ch_genes:
                merged[g] = mat[g] | pat[g]
            write_gene_vcf(ch_vcf, os.path.join(outdir, "comp_het.genes"),
                          ch_genes, merged, "comp_het")
        else:
            run(f"bcftools view -h {work} -Oz -o {ch_vcf}")   # header only
            run(f"tabix -f -p vcf {ch_vcf}")
            print(f"  {'comp_het':<14} skipped (no gene field; pass "
                  "--gene-field or use a CSQ-annotated VCF)")
        model_names.append("comp_het")

    else:  # duo: exactly one parent present
        present_parent, present_is_mother = (mom, True) if mom else (dad, False)
        parent_label = "maternal" if present_is_mother else "paternal"

        # candidate de novo: NOT distinguishable from inheritance via the
        # missing parent, so it is explicitly named/tagged as such rather
        # than asserted as de_novo the way trio mode can.
        cand_denovo = f'{gt(pro,"alt")} && {gt(present_parent,"RR")}'
        if o.min_alt > 0:
            cand_denovo += (f' && {ad_alt(pro,">=",o.min_alt)}'
                            f' && {ad_alt(present_parent,"==",0)}')
        emit("denovo_or_untransmitted", cand_denovo)
        emit(f"{parent_label}_transmitted",
             f'{gt(pro,"alt")} && {gt(present_parent,"alt")}')
        model_names = ["denovo_or_untransmitted", f"{parent_label}_transmitted"]

        if present_is_mother:
            # a son/daughter's chrX ALT + carrier mother is still informative
            # without the father; the trio model's father==RR check is a
            # specificity add-on we simply can't apply here.
            emit("x_linked_candidate", f'{gt(pro,"alt")} && {gt(mom,"RA")}', xchr)
            model_names.append("x_linked_candidate")
        else:
            print(f"  {'x_linked_candidate':<14} skipped (needs mother present)")

        # --- true compound-het trans-phasing needs BOTH parents; with one
        #     parent we can only flag genes with >=2 rare proband hets as a
        #     weaker candidate set requiring manual/2nd-parent follow-up. ---
        ch_vcf = os.path.join(outdir, "multi_het_candidate.vcf.gz")
        if gene_field:
            het = candidates(f'{gt(pro,"RA")}')
            ch_genes = sorted(g for g, sites in het.items() if len(sites) >= 2)
            write_gene_vcf(ch_vcf, os.path.join(outdir, "multi_het_candidate.genes"),
                          ch_genes, het, "multi_het_cand")
            print("                 (duo mode: trans phase unconfirmed -- "
                  "needs the missing parent or read-backed phasing)")
        else:
            run(f"bcftools view -h {work} -Oz -o {ch_vcf}")   # header only
            run(f"tabix -f -p vcf {ch_vcf}")
            print(f"  {'multi_het_cand':<14} skipped (no gene field; pass "
                  "--gene-field or use a CSQ-annotated VCF)")
        model_names.append("multi_het_candidate")

    if variant_tsv:
        src = [(m, os.path.join(outdir, f"{m}.vcf.gz")) for m in model_names]
        # rescue: proband variants that are HIGH/MODERATE impact or ClinVar P/LP,
        # kept regardless of inheritance model (built as a temp VCF, folded into
        # the TSV tagged 'high_impact', then removed — no standalone VCF output).
        hi = build_high_impact(work, pro, info, outdir)
        src.append(("high_impact", hi))
        write_variant_tsv(outdir, pro, mom, dad, src, af_field, available,
                          siblings=siblings, ped=ped)
        for f in (hi, hi + ".tbi"):
            if os.path.exists(f):
                os.remove(f)

    print(f">> done. Outputs in {outdir}/")


if __name__ == "__main__":
    main()
