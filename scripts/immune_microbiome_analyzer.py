"""
Immune System & Microbiome Host Genetics Analyzer
Analyzes genetic variants related to:
  A) Immune System Genetics — HLA tagging SNPs, complement, innate immunity
  B) Microbiome Host Genetic Predictors — FUT2, ABO, VDR, NOD2, LCT

Evidence levels:
  Strong   — replicated across multiple large cohorts (N>10,000), consistent effects
  Moderate — consistent evidence, 2-5 replication studies, smaller N
  Emerging — preliminary, single studies or conflicting results

References:
- Bonder et al. (2016) Nat Genet — FUT2 secretor status & microbiome
- Hugot et al. (2001) Nature — NOD2 and Crohn's disease
- Goodrich et al. (2016) Cell — LCT and microbiome
- Wang et al. (2020) Nat Commun — VDR and microbiome diversity
- Hosomi et al. (2022) Cell Host Microbe — ABO and gut microbiome

DISCLAIMER: This is NOT a medical diagnostic tool. HLA typing from tagging
SNPs is indicative only. Microbiome associations are probabilistic; environment,
diet, antibiotic history, and geography dominate microbiome composition.
"""

from typing import Dict, List, Tuple, Any, Optional
from dataclasses import dataclass, field


# ─────────────────────────────────────────────────────────────────────────────
# COMPLEMENT MAP (for strand-flip handling)
# ─────────────────────────────────────────────────────────────────────────────

COMPLEMENT = {"A": "T", "T": "A", "C": "G", "G": "C"}


def complement_allele(allele: str) -> str:
    """Return the complement of a single nucleotide."""
    return COMPLEMENT.get(allele.upper(), allele.upper())


def complement_genotype(genotype: str) -> str:
    """Return the complement of a genotype string."""
    return "".join(COMPLEMENT.get(c, c) for c in genotype.upper())


def count_allele(genotype: str, allele: str, try_complement: bool = True) -> int:
    """
    Count occurrences of an allele in a genotype, optionally handling complement strands.

    Args:
        genotype: Two-character genotype string (e.g. "AG")
        allele: Single nucleotide to count (e.g. "A")
        try_complement: If True, falls back to complement strand matching when the
            allele is not found. Set to False for SNPs where the ref/alt alleles
            are C/G or A/T on the same strand (not a strand flip).
    """
    genotype = genotype.upper().replace("-", "")
    allele = allele.upper()

    count = genotype.count(allele)

    if try_complement and count == 0:
        comp = complement_allele(allele)
        if comp != allele:
            count = genotype.count(comp)

    return min(count, 2)


# SNPs where ref/alt are C/G or A/T on the same strand (NOT a strand flip).
# For these, we must NOT use complement matching.
NO_COMPLEMENT_SNPS = {
    "rs2066845",   # NOD2 G908R: ref=G, alt=C (same strand, not complement)
    "rs2230199",   # C3 R102G: ref=C, alt=G (same strand)
    "rs9273363",   # HLA-DQB1: ref=A, alt=C (not complement issue but safe)
}


# ─────────────────────────────────────────────────────────────────────────────
# SNP DATA CLASSES
# ─────────────────────────────────────────────────────────────────────────────

@dataclass
class ImmuneSNP:
    rsid: str
    gene: str
    function: str
    risk_allele: str
    evidence: str
    category: str  # "hla" or "innate"
    association: str
    population_freq: str
    interpretations: Dict[str, str]  # genotype -> interpretation text
    allele_effects: Dict[str, str] = field(default_factory=dict)


@dataclass
class MicrobiomeSNP:
    rsid: str
    gene: str
    function: str
    risk_allele: str
    evidence: str
    category: str  # "fut2", "abo", "vdr", "nod2", "lct"
    interpretations: Dict[str, str]


# ─────────────────────────────────────────────────────────────────────────────
# SECTION A — IMMUNE SYSTEM SNP DATABASE
# ─────────────────────────────────────────────────────────────────────────────

IMMUNE_SNPS: Dict[str, ImmuneSNP] = {
    # ── HLA TAGGING SNPs (6) ──
    "rs2395185": ImmuneSNP(
        rsid="rs2395185",
        gene="HLA-DRB1",
        function="HLA-DRB1*13 haplotype indicator",
        risk_allele="T",
        evidence="Moderate",
        category="hla",
        association="HLA-DRB1*13 haplotype indicator — reduced susceptibility to autoimmune conditions",
        population_freq="~15-30% in Europeans",
        interpretations={
            "TT": "Homozygous for the T allele tagging HLA-DRB1*13. This haplotype is associated with reduced susceptibility to several autoimmune conditions and enhanced response to Hepatitis B vaccination.",
            "GT": "Heterozygous carrier of the HLA-DRB1*13 tagging allele. One copy of the protective haplotype marker.",
            "TG": "Heterozygous carrier of the HLA-DRB1*13 tagging allele. One copy of the protective haplotype marker.",
            "GG": "No copies of the DRB1*13 tagging allele. This is the most common genotype in European populations.",
        },
    ),
    "rs9268858": ImmuneSNP(
        rsid="rs9268858",
        gene="HLA-DRB1",
        function="HLA-A*01 associated region",
        risk_allele="G",
        evidence="Moderate",
        category="hla",
        association="HLA-A*01 supertype — pharmacogenomic relevance (abacavir sensitivity)",
        population_freq="~25% in Europeans",
        interpretations={
            "GG": "Homozygous for the G allele tagging HLA-A*01 supertype. This is associated with abacavir hypersensitivity risk and certain immune response patterns.",
            "AG": "Heterozygous carrier of HLA-A*01 tagging allele. One copy of this immune marker.",
            "GA": "Heterozygous carrier of HLA-A*01 tagging allele. One copy of this immune marker.",
            "AA": "No copies of the HLA-A*01 tagging allele. Typical genotype.",
        },
    ),
    "rs9273363": ImmuneSNP(
        rsid="rs9273363",
        gene="HLA-DQB1",
        function="HLA-DQB1*06:02 — narcolepsy risk marker",
        risk_allele="C",
        evidence="Strong",
        category="hla",
        association="HLA-DQB1*06:02 tagging SNP — strongest genetic risk factor for narcolepsy with cataplexy",
        population_freq="~25% general population; >90% of narcolepsy patients",
        interpretations={
            "CC": "Homozygous for the C allele strongly tagging HLA-DQB1*06:02. Present in >90% of narcolepsy patients, but DQB1*06:02 alone is not sufficient — environmental triggers are needed.",
            "AC": "Heterozygous carrier of the DQB1*06:02 tagging allele. One copy of this narcolepsy susceptibility marker.",
            "CA": "Heterozygous carrier of the DQB1*06:02 tagging allele. One copy of this narcolepsy susceptibility marker.",
            "AA": "No copies of the DQB1*06:02 tagging allele. Lower narcolepsy genetic susceptibility at this locus.",
        },
    ),
    "rs2187668": ImmuneSNP(
        rsid="rs2187668",
        gene="HLA-DQ2.5",
        function="HLA-DQ2.5 (DQA1*05 + DQB1*02) — celiac disease risk",
        risk_allele="C",
        evidence="Strong",
        category="hla",
        association="Celiac susceptibility marker — present in ~95% of celiac patients",
        population_freq="~25-30% of Europeans; ~95% of celiac patients",
        interpretations={
            "CC": "Homozygous for HLA-DQ2.5 tagging allele. Present in 95% of celiac patients, but only ~3% of carriers develop disease. Homozygosity significantly increases absolute risk.",
            "CT": "Heterozygous carrier of HLA-DQ2.5 tagging allele. Celiac susceptibility marker present but risk remains low without environmental triggers.",
            "TC": "Heterozygous carrier of HLA-DQ2.5 tagging allele. Celiac susceptibility marker present but risk remains low without environmental triggers.",
            "TT": "No copies of the HLA-DQ2.5 tagging allele. Very low celiac disease susceptibility at this locus.",
        },
    ),
    "rs7775228": ImmuneSNP(
        rsid="rs7775228",
        gene="HLA-DQ8",
        function="HLA-DQ8 (DQA1*03 + DQB1*03:02) — celiac/T1D risk",
        risk_allele="C",
        evidence="Moderate",
        category="hla",
        association="HLA-DQ8 tagging SNP — second major celiac-risk haplotype, also type 1 diabetes risk",
        population_freq="~15-20% in Europeans",
        interpretations={
            "CC": "Homozygous for HLA-DQ8 tagging allele. DQ8 is the second major celiac-risk haplotype (~5-10% of celiac patients who lack DQ2.5 carry DQ8). Also confers risk for type 1 diabetes.",
            "CT": "Heterozygous carrier of the DQ8 tagging allele. One copy of this celiac/T1D susceptibility marker.",
            "TC": "Heterozygous carrier of the DQ8 tagging allele. One copy of this celiac/T1D susceptibility marker.",
            "TT": "No copies of the HLA-DQ8 tagging allele at this locus.",
        },
    ),
    "rs2647012": ImmuneSNP(
        rsid="rs2647012",
        gene="HLA-DPB1",
        function="HLA-DPB1 — ankylosing spondylitis / immune signalling",
        risk_allele="A",
        evidence="Moderate",
        category="hla",
        association="HLA-DPB1 region — ankylosing spondylitis risk independent of HLA-B27",
        population_freq="~50-60% in Europeans",
        interpretations={
            "AA": "Homozygous for the A allele tagging certain HLA-DPB1 alleles associated with ankylosing spondylitis risk. Effect size is modest compared to HLA-B27.",
            "AG": "Heterozygous carrier. One copy of the DPB1 region tagging allele.",
            "GA": "Heterozygous carrier. One copy of the DPB1 region tagging allele.",
            "GG": "No copies of the HLA-DPB1 tagging allele at this marker.",
        },
    ),

    # ── COMPLEMENT & INNATE IMMUNITY (4) ──
    "rs4151667": ImmuneSNP(
        rsid="rs4151667",
        gene="CFB",
        function="Complement Factor B — alternative pathway activation",
        risk_allele="A",
        evidence="Moderate",
        category="innate",
        association="Complement pathway activation — AMD protective factor",
        population_freq="~5-10% A allele in Europeans",
        interpretations={
            "AA": "Reduced complement activity. Protective for age-related macular degeneration (AMD). Possibly lower inflammatory tone.",
            "AT": "Intermediate complement activity. One copy of the reduced-activity allele.",
            "TA": "Intermediate complement activity. One copy of the reduced-activity allele.",
            "TT": "Normal/higher complement activity. Typical for most people.",
        },
        allele_effects={
            "AA": "Reduced complement activity — protective for AMD, possibly lower inflammatory tone",
            "AT": "Intermediate complement activity",
            "TA": "Intermediate complement activity",
            "TT": "Normal/higher complement activity — typical for most people",
        },
    ),
    "rs2230199": ImmuneSNP(
        rsid="rs2230199",
        gene="C3",
        function="Complement Component 3 — R102G variant",
        risk_allele="G",
        evidence="Moderate",
        category="innate",
        association="Higher serum C3 levels / complement activation",
        population_freq="~30-35% G allele in Europeans",
        interpretations={
            "GG": "Higher C3 levels. Stronger complement opsonisation of pathogens but modestly higher AMD risk signal.",
            "CG": "Intermediate C3 levels. One copy of the higher-activity allele.",
            "GC": "Intermediate C3 levels. One copy of the higher-activity allele.",
            "CC": "Lower C3 levels. Ancestral/reference state.",
        },
        allele_effects={
            "GG": "Higher C3 levels — stronger complement; higher AMD risk signal",
            "CG": "Intermediate — one copy of higher-activity allele",
            "GC": "Intermediate — one copy of higher-activity allele",
            "CC": "Lower C3 levels — ancestral / reference state",
        },
    ),
    "rs3131379": ImmuneSNP(
        rsid="rs3131379",
        gene="MICA",
        function="MHC class I chain-related protein A — NK cell activation",
        risk_allele="T",
        evidence="Moderate",
        category="innate",
        association="NK cell activation, MHC immunogenetics",
        population_freq="Varies by population",
        interpretations={
            "TT": "Variant MICA allele. Altered NK cell activation signal. Associated with differential susceptibility to autoimmune diseases.",
            "GT": "Heterozygous. Intermediate NK cell activation profile.",
            "TG": "Heterozygous. Intermediate NK cell activation profile.",
            "GG": "Reference MICA allele. Standard NK cell activation response.",
        },
        allele_effects={
            "TT": "Variant MICA allele — altered NK cell activation signal",
            "GT": "Heterozygous — intermediate",
            "TG": "Heterozygous — intermediate",
            "GG": "Reference MICA allele",
        },
    ),
    "rs2523987": ImmuneSNP(
        rsid="rs2523987",
        gene="HLA-A/B",
        function="HLA-A/B intergenic region — MHC diversity marker",
        risk_allele="A",
        evidence="Emerging",
        category="innate",
        association="MHC immunogenetic diversity marker",
        population_freq="Varies by population",
        interpretations={
            "AA": "Variant at this MHC diversity marker. Contributes to HLA diversity assessment.",
            "AG": "Heterozygous at this MHC marker. Two different MHC haplotypes represented.",
            "GA": "Heterozygous at this MHC marker. Two different MHC haplotypes represented.",
            "GG": "Homozygous reference-region allele at this MHC marker.",
        },
        allele_effects={
            "AA": "Variant at this MHC marker",
            "AG": "Heterozygous — two MHC haplotypes",
            "GA": "Heterozygous — two MHC haplotypes",
            "GG": "Homozygous reference-region allele",
        },
    ),
}


# ─────────────────────────────────────────────────────────────────────────────
# SECTION B — MICROBIOME HOST FACTOR SNP DATABASE
# ─────────────────────────────────────────────────────────────────────────────

MICROBIOME_SNPS: Dict[str, MicrobiomeSNP] = {
    # ── FUT2 Secretor Status ──
    "rs601338": MicrobiomeSNP(
        rsid="rs601338",
        gene="FUT2",
        function="Fucosyltransferase 2 — secretor status (W154X)",
        risk_allele="A",
        evidence="Strong",
        category="fut2",
        interpretations={
            "GG": "Secretor (GG). Full FUT2 activity. Gut mucosa is well-coated with H-antigen fucose residues, promoting Bifidobacterium colonization. Higher susceptibility to Norovirus GII and H. pylori attachment.",
            "GA": "Secretor (GA). One functional FUT2 allele is sufficient for full secretor phenotype. Microbiome profile similar to GG with normal Bifidobacterium colonization.",
            "AG": "Secretor (AG). One functional FUT2 allele is sufficient for full secretor phenotype. Microbiome profile similar to GG with normal Bifidobacterium colonization.",
            "AA": "Non-secretor (AA). FUT2 enzyme completely absent. Dramatically lower Bifidobacterium (20-50x reduction). Different microbial community structure. Lower susceptibility to Norovirus GII but potentially higher risk of Crohn's disease.",
        },
    ),

    # ── ABO Blood Group SNPs ──
    "rs8176719": MicrobiomeSNP(
        rsid="rs8176719",
        gene="ABO",
        function="ABO glycosyltransferase — O vs non-O determination",
        risk_allele="T",
        evidence="Moderate",
        category="abo",
        interpretations={
            "CC": "Non-O on both chromosomes (likely A/A, A/B, or B/B).",
            "CT": "Likely A/O or B/O (non-O/O heterozygote).",
            "TC": "Likely A/O or B/O (non-O/O heterozygote).",
            "TT": "Likely O/O blood group.",
        },
    ),
    "rs505922": MicrobiomeSNP(
        rsid="rs505922",
        gene="ABO",
        function="ABO locus broad marker — C=non-O",
        risk_allele="C",
        evidence="Moderate",
        category="abo",
        interpretations={
            "CC": "Non-O ABO type marker (both alleles).",
            "CT": "One non-O allele detected.",
            "TC": "One non-O allele detected.",
            "TT": "O-type marker at this locus.",
        },
    ),
    "rs8176746": MicrobiomeSNP(
        rsid="rs8176746",
        gene="ABO",
        function="ABO A vs B differentiation (Tyr316Ser)",
        risk_allele="T",
        evidence="Moderate",
        category="abo",
        interpretations={
            "CC": "A-type or O-type alleles (no B allele detected).",
            "CT": "Possible AB or admixed A/B — both A and B alleles present.",
            "TC": "Possible AB or admixed A/B — both A and B alleles present.",
            "TT": "B-type or O-type alleles (no A allele at this marker).",
        },
    ),
    "rs8176747": MicrobiomeSNP(
        rsid="rs8176747",
        gene="ABO",
        function="ABO refinement marker",
        risk_allele="T",
        evidence="Moderate",
        category="abo",
        interpretations={
            "CC": "Reference allele at this ABO refinement marker.",
            "CT": "Heterozygous at this ABO refinement marker.",
            "TC": "Heterozygous at this ABO refinement marker.",
            "TT": "Variant allele at this ABO refinement marker.",
        },
    ),
    "rs1053878": MicrobiomeSNP(
        rsid="rs1053878",
        gene="ABO",
        function="ABO A1 vs A2 subtype refinement",
        risk_allele="G",
        evidence="Moderate",
        category="abo",
        interpretations={
            "AA": "Reference allele — typical A1 subtype if A-type.",
            "AG": "Heterozygous — possible A2 subtype indicator.",
            "GA": "Heterozygous — possible A2 subtype indicator.",
            "GG": "Variant allele — may indicate A2 subtype.",
        },
    ),

    # ── VDR (Vitamin D Receptor) ──
    "rs2228570": MicrobiomeSNP(
        rsid="rs2228570",
        gene="VDR",
        function="VDR FokI (C>T) — vitamin D receptor activity",
        risk_allele="T",
        evidence="Moderate",
        category="vdr",
        interpretations={
            "CC": "FF genotype — most active VDR form. Strongest vitamin D response, supporting antimicrobial peptide (defensin) production in the gut.",
            "CT": "Ff genotype — intermediate VDR activity. Moderate defensin production.",
            "TC": "Ff genotype — intermediate VDR activity. Moderate defensin production.",
            "TT": "ff genotype — less active VDR. Reduced intestinal defensin production may affect gut microbiome diversity.",
        },
    ),
    "rs731236": MicrobiomeSNP(
        rsid="rs731236",
        gene="VDR",
        function="VDR TaqI (T>C) — vitamin D receptor expression",
        risk_allele="C",
        evidence="Moderate",
        category="vdr",
        interpretations={
            "TT": "TT genotype — reference. Normal VDR expression and activity.",
            "TC": "Tt genotype — intermediate VDR expression.",
            "CT": "Tt genotype — intermediate VDR expression.",
            "CC": "tt genotype — associated with lower VDR activity in some studies. May reduce antimicrobial peptide production.",
        },
    ),
    "rs1544410": MicrobiomeSNP(
        rsid="rs1544410",
        gene="VDR",
        function="VDR BsmI (G>A) — vitamin D receptor expression",
        risk_allele="A",
        evidence="Moderate",
        category="vdr",
        interpretations={
            "GG": "bb genotype — reference. Normal VDR expression.",
            "GA": "Bb genotype — intermediate VDR expression.",
            "AG": "Bb genotype — intermediate VDR expression.",
            "AA": "BB genotype — lower VDR expression in some studies. Associated with reduced microbiome diversity in some cohorts.",
        },
    ),

    # ── NOD2 (Bacterial Sensing) ──
    "rs2066844": MicrobiomeSNP(
        rsid="rs2066844",
        gene="NOD2",
        function="NOD2 R702W (Arg702Trp) — bacterial peptidoglycan sensing",
        risk_allele="T",
        evidence="Moderate",
        category="nod2",
        interpretations={
            "CC": "No R702W variant. Normal NOD2 function and bacterial peptidoglycan sensing in Paneth cells.",
            "CT": "Heterozygous R702W. Mildly reduced NF-kB activation from bacterial sensing. Most carriers never develop Crohn's disease.",
            "TC": "Heterozygous R702W. Mildly reduced NF-kB activation from bacterial sensing. Most carriers never develop Crohn's disease.",
            "TT": "Homozygous R702W. Substantially impaired bacterial sensing. Reduced defensin secretion may alter gut microbiome composition.",
        },
    ),
    "rs2066845": MicrobiomeSNP(
        rsid="rs2066845",
        gene="NOD2",
        function="NOD2 G908R (Gly908Arg) — LRR domain function",
        risk_allele="C",
        evidence="Moderate",
        category="nod2",
        interpretations={
            "GG": "No G908R variant. Normal NOD2 LRR domain function.",
            "GC": "Heterozygous G908R. Impaired LRR domain function. Modest increase in Crohn's disease susceptibility.",
            "CG": "Heterozygous G908R. Impaired LRR domain function. Modest increase in Crohn's disease susceptibility.",
            "CC": "Homozygous G908R. Rare. Severely impaired NOD2 function.",
        },
    ),

    # ── LCT (Lactase Persistence) ──
    "rs4988235": MicrobiomeSNP(
        rsid="rs4988235",
        gene="LCT",
        function="Lactase persistence (MCM6 enhancer, -13910 C>T)",
        risk_allele="G",
        evidence="Strong",
        category="lct",
        interpretations={
            "AA": "Lactase persistent (AA). Produces lactase throughout life. When consuming dairy, lactose is digested in the small intestine, limiting substrate for colonic fermenters.",
            "AG": "Lactase persistent (AG). One LP allele is sufficient (dominant). Lactase production maintained into adulthood.",
            "GA": "Lactase persistent (GA). One LP allele is sufficient (dominant). Lactase production maintained into adulthood.",
            "GG": "Lactase non-persistent (GG). Ancestral genotype. Lactase production declines after weaning. Undigested lactose reaches the colon, fueling fermentation by Bifidobacterium and altering short-chain fatty acid production.",
        },
    ),
}


# ─────────────────────────────────────────────────────────────────────────────
# ANALYSIS RESULT DATACLASSES
# ─────────────────────────────────────────────────────────────────────────────

@dataclass
class ImmuneVariantResult:
    rsid: str
    gene: str
    genotype: str
    risk_copies: int
    association: str
    evidence: str
    interpretation: str
    note: str = ""


@dataclass
class MicrobiomeVariantResult:
    rsid: str
    gene: str
    genotype: str
    status: str
    interpretation: str
    evidence: str
    microbiome_impact: str  # "high", "moderate", "low"


@dataclass
class ImmuneMicrobiomeResult:
    # Immune section
    hla_markers: List[ImmuneVariantResult]
    innate_markers: List[ImmuneVariantResult]
    hla_diversity_score: float  # 0.0 - 1.0
    hla_diversity_label: str
    hla_diversity_interpretation: str
    immune_snps_typed: int
    immune_snps_missing: List[str]

    # Microbiome section
    secretor_status: str  # "Secretor", "Non-secretor", "Unknown"
    secretor_result: Optional[MicrobiomeVariantResult]
    blood_group: str
    blood_group_confidence: str
    blood_group_interpretation: str
    abo_snps_found: Dict[str, str]  # rsid -> genotype
    vdr_variants: List[MicrobiomeVariantResult]
    vdr_interpretation: str
    nod2_variants: List[MicrobiomeVariantResult]
    nod2_status: str
    nod2_interpretation: str
    lactase_result: Optional[MicrobiomeVariantResult]
    microbiome_snps_typed: int
    microbiome_snps_missing: List[str]

    # Summary
    key_findings: List[str]
    overall_profile: str


# ─────────────────────────────────────────────────────────────────────────────
# HELPER: GENOTYPE LOOKUP WITH COMPLEMENT STRAND HANDLING
# ─────────────────────────────────────────────────────────────────────────────

def _get_interpretation(snp_interpretations: Dict[str, str], genotype: str) -> str:
    """Look up interpretation, trying complement strand if needed."""
    gt = genotype.upper().replace("-", "")
    if gt in snp_interpretations:
        return snp_interpretations[gt]
    # Try reversed genotype (e.g. AG -> GA)
    gt_rev = gt[::-1]
    if gt_rev in snp_interpretations:
        return snp_interpretations[gt_rev]
    # Try complement
    gt_comp = complement_genotype(gt)
    if gt_comp in snp_interpretations:
        return snp_interpretations[gt_comp]
    # Try complement reversed
    gt_comp_rev = gt_comp[::-1]
    if gt_comp_rev in snp_interpretations:
        return snp_interpretations[gt_comp_rev]
    return f"Genotype {gt} detected. Interpretation not available for this combination."


# ─────────────────────────────────────────────────────────────────────────────
# ABO BLOOD GROUP INFERENCE
# ─────────────────────────────────────────────────────────────────────────────

def _infer_blood_group(variants: Dict[str, Tuple[str, str, str]]) -> Tuple[str, str, str]:
    """
    Infer ABO blood group from available SNPs.
    Returns (group, confidence, interpretation).
    """
    abo_snps = {}
    for rsid in ["rs8176719", "rs505922", "rs8176746", "rs8176747", "rs1053878"]:
        if rsid.lower() in variants:
            _, _, gt = variants[rsid.lower()]
            abo_snps[rsid] = gt.upper().replace("-", "")

    if not abo_snps:
        return "Undetermined", "low", "No ABO SNPs found in dataset. Cannot infer blood group."

    inferred = "Undetermined"
    confidence = "low"

    # rs8176746: CC = A or O; TT = B or O; CT/TC = AB-like
    gt_746 = abo_snps.get("rs8176746")
    if gt_746:
        if gt_746 == "CC":
            inferred = "Likely A or O"
            confidence = "moderate"
        elif gt_746 == "TT":
            inferred = "Likely B or O"
            confidence = "moderate"
        elif gt_746 in ("CT", "TC"):
            inferred = "Possible AB"
            confidence = "low"

    # rs8176719: refines O vs non-O
    gt_719 = abo_snps.get("rs8176719")
    if gt_719:
        if gt_719 == "TT":
            inferred = "Likely O"
            confidence = "moderate"
        elif gt_719 == "CC":
            if gt_746 == "CC":
                inferred = "Likely A"
                confidence = "moderate"
            elif gt_746 == "TT":
                inferred = "Likely B"
                confidence = "moderate"
            elif gt_746 in ("CT", "TC"):
                inferred = "Likely AB"
                confidence = "moderate"
            else:
                inferred = "Likely A, B, or AB"
                confidence = "low"
        elif gt_719 in ("CT", "TC"):
            if gt_746 == "CC":
                inferred = "Likely A/O"
                confidence = "moderate"
            elif gt_746 == "TT":
                inferred = "Likely B/O"
                confidence = "moderate"
            else:
                inferred = "Likely A/O or B/O"
                confidence = "moderate"

    # rs505922 as fallback
    if inferred == "Undetermined":
        gt_922 = abo_snps.get("rs505922")
        if gt_922:
            if gt_922 == "TT":
                inferred = "Likely O"
                confidence = "low"
            elif gt_922 == "CC":
                inferred = "Likely non-O"
                confidence = "low"
            elif gt_922 in ("CT", "TC"):
                inferred = "Possible O carrier"
                confidence = "low"

    interp_map = {
        "Likely O": "Blood group O is associated with higher Prevotella copri abundance in some cohorts and different Bacteroides community patterns.",
        "Likely A": "Blood group A is associated with specific patterns of bacterial adhesion to intestinal epithelium and higher Ruminococcaceae abundance.",
        "Likely B": "Blood group B shows an intermediate microbiome profile, similar to O in some studies.",
        "Likely AB": "Blood group AB combines A and B antigens on intestinal epithelium, potentially supporting the widest range of bacterial adhesion patterns.",
        "Likely A/O": "Likely blood group A (carrying one O allele). Associated with specific bacterial adhesion patterns on intestinal epithelium.",
        "Likely B/O": "Likely blood group B (carrying one O allele). Intermediate microbiome profile in most studies.",
    }
    interpretation = interp_map.get(inferred,
        "ABO blood group shapes the glycan scaffold on intestinal mucosa, influencing bacterial colonization. SNP array typing is imprecise for ABO.")

    return inferred, confidence, interpretation


# ─────────────────────────────────────────────────────────────────────────────
# HLA DIVERSITY SCORE
# ─────────────────────────────────────────────────────────────────────────────

def _calculate_hla_diversity(hla_results: List[ImmuneVariantResult]) -> Tuple[float, str, str]:
    """
    Calculate HLA diversity score based on heterozygosity across HLA tagging SNPs.
    Returns (score, label, interpretation).
    """
    if not hla_results:
        return 0.0, "Unknown", "No HLA markers typed."

    typed_count = 0
    heterozygous_count = 0

    for marker in hla_results:
        gt = marker.genotype.upper().replace("-", "")
        if len(gt) == 2 and gt != "NN":
            typed_count += 1
            if gt[0] != gt[1]:
                heterozygous_count += 1

    if typed_count == 0:
        return 0.0, "Unknown", "No HLA markers could be assessed."

    score = round(heterozygous_count / typed_count, 2)

    if score >= 0.80:
        label = "High"
        interp = ("Your HLA heterozygosity is high, suggesting a broad pathogen recognition "
                  "repertoire. High HLA diversity is generally advantageous for immune "
                  "defense against a wide range of pathogens.")
    elif score >= 0.50:
        label = "Above Average"
        interp = ("Your HLA heterozygosity is above the European average, suggesting a broader "
                  "pathogen recognition repertoire. This generally supports diverse immune responses.")
    elif score >= 0.30:
        label = "Average"
        interp = ("Your HLA heterozygosity is within the typical range for European populations. "
                  "Your immune repertoire diversity is normal.")
    else:
        label = "Below Average"
        interp = ("Your HLA heterozygosity is below the European average. This may indicate "
                  "more homogeneous HLA alleles, which is common in certain ancestral backgrounds. "
                  "It does not indicate immune deficiency.")

    return score, label, interp


# ─────────────────────────────────────────────────────────────────────────────
# CORE ANALYSIS
# ─────────────────────────────────────────────────────────────────────────────

def analyze_immune_microbiome(
    variants: Dict[str, Tuple[str, str, str]]
) -> ImmuneMicrobiomeResult:
    """
    Analyze genetic variants for immune system and microbiome host factors.

    Args:
        variants: Dictionary mapping rsID (lowercase) to (chromosome, position, genotype)

    Returns:
        ImmuneMicrobiomeResult with full immune and microbiome analysis
    """
    # ── SECTION A: IMMUNE SYSTEM ──
    hla_markers: List[ImmuneVariantResult] = []
    innate_markers: List[ImmuneVariantResult] = []
    immune_missing: List[str] = []
    immune_typed = 0

    for rsid, snp in IMMUNE_SNPS.items():
        rsid_lower = rsid.lower()
        if rsid_lower in variants:
            immune_typed += 1
            _chrom, _pos, genotype = variants[rsid_lower]
            genotype = genotype.upper().replace("-", "")
            use_comp = rsid not in NO_COMPLEMENT_SNPS
            risk_copies = count_allele(genotype, snp.risk_allele, try_complement=use_comp)
            interpretation = _get_interpretation(snp.interpretations, genotype)

            # Build note for HLA markers
            note = ""
            if snp.category == "hla":
                if "celiac" in snp.association.lower():
                    note = "Present in 95% of celiac patients but only 3% of carriers develop disease"
                elif "narcolepsy" in snp.association.lower():
                    note = "DQB1*06:02 alone is not sufficient for narcolepsy — environmental triggers are needed"
                elif "DRB1*13" in snp.association:
                    note = "Proxy SNP — true HLA typing requires specialized assays"

            result = ImmuneVariantResult(
                rsid=rsid,
                gene=snp.gene,
                genotype=genotype,
                risk_copies=risk_copies,
                association=snp.association,
                evidence=snp.evidence,
                interpretation=interpretation,
                note=note,
            )

            if snp.category == "hla":
                hla_markers.append(result)
            else:
                innate_markers.append(result)
        else:
            immune_missing.append(f"{rsid} ({snp.gene})")

    # Calculate HLA diversity
    hla_diversity_score, hla_diversity_label, hla_diversity_interp = _calculate_hla_diversity(hla_markers)

    # ── SECTION B: MICROBIOME ──
    microbiome_typed = 0
    microbiome_missing: List[str] = []
    key_findings: List[str] = []

    # B1: FUT2 Secretor Status
    secretor_result = None
    secretor_status = "Unknown"
    fut2_snp = MICROBIOME_SNPS["rs601338"]
    if "rs601338" in variants:
        microbiome_typed += 1
        _chrom, _pos, genotype = variants["rs601338"]
        genotype = genotype.upper().replace("-", "")
        interpretation = _get_interpretation(fut2_snp.interpretations, genotype)

        g_count = count_allele(genotype, "G")
        a_count = count_allele(genotype, "A")

        if a_count == 2:
            secretor_status = "Non-secretor"
            status_str = "Non-secretor (AA)"
            impact = "high"
            key_findings.append("Non-secretor status (FUT2 AA) — dramatically different Bifidobacterium levels")
        elif g_count >= 1:
            secretor_status = "Secretor"
            status_str = f"Secretor ({genotype})"
            impact = "high"
            key_findings.append(f"Secretor status (FUT2 {genotype}) — promotes Bifidobacterium colonization")
        else:
            secretor_status = "Unknown"
            status_str = f"Ambiguous ({genotype})"
            impact = "low"

        secretor_result = MicrobiomeVariantResult(
            rsid="rs601338",
            gene="FUT2",
            genotype=genotype,
            status=status_str,
            interpretation=interpretation,
            evidence="Strong",
            microbiome_impact=impact,
        )
    else:
        microbiome_missing.append("rs601338 (FUT2)")

    # B2: ABO Blood Group
    abo_snps_found: Dict[str, str] = {}
    abo_rsids = ["rs8176719", "rs505922", "rs8176746", "rs8176747", "rs1053878"]
    for rsid in abo_rsids:
        rsid_lower = rsid.lower()
        if rsid_lower in variants:
            microbiome_typed += 1
            _, _, gt = variants[rsid_lower]
            abo_snps_found[rsid] = gt.upper().replace("-", "")
        else:
            microbiome_missing.append(f"{rsid} (ABO)")

    blood_group, bg_confidence, bg_interpretation = _infer_blood_group(variants)

    # B3: VDR
    vdr_variants: List[MicrobiomeVariantResult] = []
    vdr_rsids = ["rs2228570", "rs731236", "rs1544410"]
    vdr_risk_count = 0
    vdr_typed = 0
    for rsid in vdr_rsids:
        rsid_lower = rsid.lower()
        snp = MICROBIOME_SNPS[rsid]
        if rsid_lower in variants:
            microbiome_typed += 1
            vdr_typed += 1
            _chrom, _pos, genotype = variants[rsid_lower]
            genotype = genotype.upper().replace("-", "")
            use_comp = rsid not in NO_COMPLEMENT_SNPS
            risk_copies = count_allele(genotype, snp.risk_allele, try_complement=use_comp)
            vdr_risk_count += risk_copies
            interpretation = _get_interpretation(snp.interpretations, genotype)

            vdr_variants.append(MicrobiomeVariantResult(
                rsid=rsid,
                gene="VDR",
                genotype=genotype,
                status=f"{risk_copies} risk allele(s)",
                interpretation=interpretation,
                evidence="Moderate",
                microbiome_impact="moderate" if risk_copies > 0 else "low",
            ))
        else:
            microbiome_missing.append(f"{rsid} (VDR)")

    # VDR overall interpretation
    if vdr_typed == 0:
        vdr_interpretation = "No VDR variants typed. Cannot assess vitamin D receptor function."
    elif vdr_risk_count == 0:
        vdr_interpretation = ("Your VDR variants suggest normal vitamin D receptor function, "
                              "supporting antimicrobial peptide production and healthy gut microbiome diversity.")
    elif vdr_risk_count <= 2:
        vdr_interpretation = ("Your VDR variants suggest mildly reduced vitamin D receptor activity. "
                              "Adequate vitamin D intake may help compensate, supporting defensin production and microbiome diversity.")
    else:
        vdr_interpretation = ("Your VDR variants suggest reduced vitamin D receptor activity across multiple markers. "
                              "This may affect antimicrobial peptide production and gut microbiome diversity. "
                              "Ensuring adequate vitamin D levels through supplementation may be particularly beneficial.")

    # B4: NOD2
    nod2_variants: List[MicrobiomeVariantResult] = []
    nod2_rsids = ["rs2066844", "rs2066845"]
    nod2_risk_count = 0
    nod2_typed = 0
    for rsid in nod2_rsids:
        rsid_lower = rsid.lower()
        snp = MICROBIOME_SNPS[rsid]
        if rsid_lower in variants:
            microbiome_typed += 1
            nod2_typed += 1
            _chrom, _pos, genotype = variants[rsid_lower]
            genotype = genotype.upper().replace("-", "")
            use_comp = rsid not in NO_COMPLEMENT_SNPS
            risk_copies = count_allele(genotype, snp.risk_allele, try_complement=use_comp)
            nod2_risk_count += risk_copies
            interpretation = _get_interpretation(snp.interpretations, genotype)

            nod2_variants.append(MicrobiomeVariantResult(
                rsid=rsid,
                gene="NOD2",
                genotype=genotype,
                status="Normal" if risk_copies == 0 else f"{risk_copies} risk allele(s)",
                interpretation=interpretation,
                evidence="Moderate",
                microbiome_impact="moderate" if risk_copies > 0 else "low",
            ))
        else:
            microbiome_missing.append(f"{rsid} (NOD2)")

    # NOD2 overall
    if nod2_typed == 0:
        nod2_status = "Not typed"
        nod2_interpretation = "No NOD2 variants typed. Cannot assess bacterial sensing function."
    elif nod2_risk_count == 0:
        nod2_status = "Normal"
        nod2_interpretation = ("No NOD2 risk variants detected. Normal bacterial sensing in "
                               "intestinal epithelium. Paneth cell defensin production expected to be normal.")
    elif nod2_risk_count == 1:
        nod2_status = "One risk allele"
        nod2_interpretation = ("One NOD2 risk variant allele detected. Mildly reduced bacterial "
                               "peptidoglycan sensing possible. Most carriers never develop Crohn's disease. "
                               "Gut microbiome composition may be subtly altered.")
        key_findings.append("NOD2 heterozygous variant — mildly reduced bacterial sensing")
    else:
        nod2_status = "Multiple risk alleles"
        nod2_interpretation = (f"{nod2_risk_count} NOD2 risk variant alleles detected. "
                               "Substantially impaired NOD2 signaling may affect gut immunity "
                               "and microbiome composition. This does NOT diagnose Crohn's disease.")
        key_findings.append(f"NOD2: {nod2_risk_count} risk alleles — impaired bacterial sensing")

    # B5: LCT Lactase
    lactase_result = None
    lct_snp = MICROBIOME_SNPS["rs4988235"]
    if "rs4988235" in variants:
        microbiome_typed += 1
        _chrom, _pos, genotype = variants["rs4988235"]
        genotype = genotype.upper().replace("-", "")
        interpretation = _get_interpretation(lct_snp.interpretations, genotype)

        a_count = count_allele(genotype, "A")
        if a_count >= 1:
            lct_status = "Persistent"
            impact = "moderate"
            key_findings.append(f"Lactase persistent ({genotype})")
        else:
            lct_status = "Non-persistent"
            impact = "high"
            key_findings.append("Lactase non-persistent — undigested lactose fuels colonic fermentation")

        lactase_result = MicrobiomeVariantResult(
            rsid="rs4988235",
            gene="LCT",
            genotype=genotype,
            status=lct_status,
            interpretation=interpretation,
            evidence="Strong",
            microbiome_impact=impact,
        )
    else:
        microbiome_missing.append("rs4988235 (LCT)")

    # HLA diversity finding
    if hla_diversity_score >= 0.50:
        key_findings.insert(0, f"High HLA heterozygosity (score: {hla_diversity_score})")
    elif hla_markers:
        key_findings.insert(0, f"HLA heterozygosity score: {hla_diversity_score}")

    # Build overall microbiome profile
    profile_parts = []
    if secretor_status == "Secretor":
        profile_parts.append("secretor status promoting Bifidobacterium colonization")
    elif secretor_status == "Non-secretor":
        profile_parts.append("non-secretor status with reduced Bifidobacterium attachment sites")
    if lactase_result:
        if lactase_result.status == "Non-persistent":
            profile_parts.append("lactase non-persistence driving colonic lactose fermentation")
        else:
            profile_parts.append("lactase persistence allowing dairy-derived substrate absorption")
    if nod2_risk_count > 0:
        profile_parts.append("altered NOD2 bacterial sensing")
    if vdr_risk_count > 2:
        profile_parts.append("reduced VDR-mediated antimicrobial defense")

    if profile_parts:
        overall_profile = (
            "Your genetic profile suggests a microbiome shaped by " +
            ", ".join(profile_parts) + ". "
            "Remember that diet, environment, and lifestyle are the dominant factors "
            "determining gut microbiome composition (80-90% of variance), "
            "while genetics contributes approximately 10-20%."
        )
    else:
        overall_profile = (
            "Limited microbiome-relevant SNPs were typed. Your gut microbiome is primarily "
            "shaped by diet, environment, antibiotic history, and lifestyle factors, "
            "which account for 80-90% of microbiome variance."
        )

    return ImmuneMicrobiomeResult(
        hla_markers=hla_markers,
        innate_markers=innate_markers,
        hla_diversity_score=hla_diversity_score,
        hla_diversity_label=hla_diversity_label,
        hla_diversity_interpretation=hla_diversity_interp,
        immune_snps_typed=immune_typed,
        immune_snps_missing=immune_missing,
        secretor_status=secretor_status,
        secretor_result=secretor_result,
        blood_group=blood_group,
        blood_group_confidence=bg_confidence,
        blood_group_interpretation=bg_interpretation,
        abo_snps_found=abo_snps_found,
        vdr_variants=vdr_variants,
        vdr_interpretation=vdr_interpretation,
        nod2_variants=nod2_variants,
        nod2_status=nod2_status,
        nod2_interpretation=nod2_interpretation,
        lactase_result=lactase_result,
        microbiome_snps_typed=microbiome_typed,
        microbiome_snps_missing=microbiome_missing,
        key_findings=key_findings,
        overall_profile=overall_profile,
    )


# ─────────────────────────────────────────────────────────────────────────────
# JSON OUTPUT
# ─────────────────────────────────────────────────────────────────────────────

def generate_immune_microbiome_json(result: ImmuneMicrobiomeResult) -> dict:
    """
    Generate a JSON-serializable dict for the frontend immune/microbiome report.

    Args:
        result: ImmuneMicrobiomeResult from analyze_immune_microbiome()

    Returns:
        Dict matching the immune_microbiome report JSON schema
    """
    # Build HLA diversity markers
    hla_markers_json = []
    for m in result.hla_markers:
        hla_markers_json.append({
            "rsid": m.rsid,
            "gene": m.gene,
            "genotype": m.genotype,
            "riskCopies": m.risk_copies,
            "association": m.association,
            "evidence": m.evidence,
            "interpretation": m.interpretation,
            "note": m.note,
        })

    # Build innate immunity markers
    innate_markers_json = []
    for m in result.innate_markers:
        innate_markers_json.append({
            "rsid": m.rsid,
            "gene": m.gene,
            "genotype": m.genotype,
            "status": m.interpretation,
            "evidence": m.evidence,
        })

    # Build secretor status
    secretor_json = None
    if result.secretor_result:
        secretor_json = {
            "gene": "FUT2",
            "rsid": result.secretor_result.rsid,
            "genotype": result.secretor_result.genotype,
            "status": result.secretor_status,
            "interpretation": result.secretor_result.interpretation,
            "microbiomeImpact": result.secretor_result.microbiome_impact,
        }

    # Build blood group
    blood_group_json = {
        "inferred": result.blood_group,
        "confidence": result.blood_group_confidence,
        "interpretation": result.blood_group_interpretation,
        "microbiomeImpact": "low",
        "snpsFound": result.abo_snps_found,
    }

    # Build VDR section
    vdr_variants_json = []
    for v in result.vdr_variants:
        vdr_variants_json.append({
            "rsid": v.rsid,
            "gene": v.gene,
            "genotype": v.genotype,
            "status": v.status,
            "interpretation": v.interpretation,
            "evidence": v.evidence,
        })

    vdr_json = {
        "gene": "VDR",
        "variants": vdr_variants_json,
        "interpretation": result.vdr_interpretation,
        "microbiomeImpact": "moderate" if any(v.microbiome_impact != "low" for v in result.vdr_variants) else "low",
    }

    # Build NOD2 section
    nod2_variants_json = []
    for v in result.nod2_variants:
        nod2_variants_json.append({
            "rsid": v.rsid,
            "gene": v.gene,
            "genotype": v.genotype,
            "status": v.status,
            "interpretation": v.interpretation,
            "evidence": v.evidence,
        })

    nod2_json = {
        "gene": "NOD2",
        "variants": nod2_variants_json,
        "status": result.nod2_status,
        "interpretation": result.nod2_interpretation,
        "microbiomeImpact": "moderate" if result.nod2_status != "Normal" else "low",
    }

    # Build lactase section
    lactase_json = None
    if result.lactase_result:
        lactase_json = {
            "gene": "LCT",
            "rsid": result.lactase_result.rsid,
            "genotype": result.lactase_result.genotype,
            "status": result.lactase_result.status,
            "interpretation": result.lactase_result.interpretation,
            "microbiomeImpact": result.lactase_result.microbiome_impact,
        }

    return {
        "reportType": "immune_microbiome",
        "version": "1.0",
        "summary": {
            "hlaDiversityScore": result.hla_diversity_score,
            "hlaDiversityLabel": result.hla_diversity_label,
            "secretorStatus": result.secretor_status,
            "bloodGroup": result.blood_group,
            "immuneSnpsTyped": result.immune_snps_typed,
            "microbiomeSnpsTyped": result.microbiome_snps_typed,
            "keyFindings": result.key_findings,
        },
        "immuneSystem": {
            "title": "Immune System Genetics",
            "description": (
                "Your inherited immune repertoire based on HLA and innate immunity markers. "
                "The immune system is the most genetically diverse system in the human body. "
                "HLA genes on chromosome 6 are the most polymorphic genes known, determining "
                "which pathogens your immune system can recognize."
            ),
            "hlaDiversity": {
                "score": result.hla_diversity_score,
                "maxScore": 1.0,
                "label": result.hla_diversity_label,
                "interpretation": result.hla_diversity_interpretation,
                "markers": hla_markers_json,
            },
            "innateImmunity": innate_markers_json,
        },
        "microbiome": {
            "title": "Microbiome Genetic Predisposition",
            "description": (
                "How your genome influences your gut microbiome composition. "
                "Host genetics explains roughly 10-20% of gut microbiome variance. "
                "The remaining 80-90% is determined by diet, environment, antibiotic "
                "history, birth mode, geography, and life stage."
            ),
            "secretorStatus": secretor_json,
            "bloodGroup": blood_group_json,
            "vitaminD": vdr_json,
            "bacterialSensing": nod2_json,
            "lactase": lactase_json,
            "overallProfile": result.overall_profile,
        },
    }
