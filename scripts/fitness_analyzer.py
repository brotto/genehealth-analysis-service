"""
Fitness Analyzer
Analyzes genetic variants related to athletic performance across four dimensions:
Power/Sprint, Endurance, Recovery, and Injury Risk.

References:
- Yang et al. (2003) AJHG - ACTN3 R577X
- Jones et al. (2002) J Physiol - ACE I/D
- Ahmetov & Fedotovskaya (2015) Biol Sport
- September et al. (2009) Br J Sports Med - COL5A1
- Posthumus et al. (2009) AJSM - COL1A1, MMP3
- Fishman et al. (1998) J Clin Invest - IL6
- Eynon et al. (2011) J Sci Med Sport - PPARGC1A
"""

from typing import Dict, List, Tuple, Any
from dataclasses import dataclass, field


# ─────────────────────────────────────────────────────────────────────────────
# COMPLEMENT MAP (for strand-flip handling)
# ─────────────────────────────────────────────────────────────────────────────

COMPLEMENT = {"A": "T", "T": "A", "C": "G", "G": "C"}


def complement_allele(allele: str) -> str:
    """Return the complement of a single nucleotide."""
    return COMPLEMENT.get(allele.upper(), allele.upper())


def count_allele(genotype: str, allele: str) -> int:
    """
    Count occurrences of an allele in a genotype, handling complement strands.
    For ambiguous SNPs (A/T or C/G), checks both the allele and its complement.
    """
    genotype = genotype.upper().replace("-", "")
    allele = allele.upper()
    comp = complement_allele(allele)

    count = genotype.count(allele)

    # If the allele is A/T or C/G (ambiguous), complement is the same pair
    # so we should NOT double-count. Only count complement if it differs
    # from the original allele and the genotype doesn't contain the original.
    if count == 0 and comp != allele:
        count = genotype.count(comp)

    return min(count, 2)


# ─────────────────────────────────────────────────────────────────────────────
# DIMENSION DEFINITIONS
# ─────────────────────────────────────────────────────────────────────────────

DIMENSIONS = ["Power/Sprint", "Endurance", "Recovery", "Injury Risk"]

DIMENSION_EMOJI = {
    "Power/Sprint": "\U0001f4aa",   # 💪
    "Endurance": "\u26a1",          # ⚡
    "Recovery": "\U0001fa79",       # 🩹
    "Injury Risk": "\U0001f9b4",    # 🦴
}


# ─────────────────────────────────────────────────────────────────────────────
# FITNESS VARIANT DATABASE
# ─────────────────────────────────────────────────────────────────────────────

@dataclass
class FitnessSNP:
    rsid: str
    gene: str
    dimension: str
    trait: str
    risk_allele: str        # allele associated with LESS favorable outcome
    favorable_allele: str   # allele associated with MORE favorable outcome
    effect: str             # description when favorable allele is present
    evidence: str           # "Strong", "Moderate", or "Preliminary"


FITNESS_VARIANTS: Dict[str, FitnessSNP] = {
    # ── POWER / SPRINT ──
    "rs1815739": FitnessSNP(
        rsid="rs1815739", gene="ACTN3", dimension="Power/Sprint",
        trait="Fast-twitch muscle fiber composition",
        risk_allele="T", favorable_allele="C",
        effect="RR genotype (CC) associated with elite sprint/power performance; "
               "alpha-actinin-3 protein present in fast-twitch fibers",
        evidence="Strong",
    ),
    "rs17602729": FitnessSNP(
        rsid="rs17602729", gene="AMPD1", dimension="Power/Sprint",
        trait="Muscle energy metabolism (AMP deaminase)",
        risk_allele="A", favorable_allele="G",
        effect="Normal AMP deaminase activity supports high-intensity "
               "anaerobic energy production",
        evidence="Moderate",
    ),
    "rs699": FitnessSNP(
        rsid="rs699", gene="AGT", dimension="Power/Sprint",
        trait="Angiotensinogen & muscle hypertrophy",
        risk_allele="G", favorable_allele="A",
        effect="M235T variant (A allele) linked to greater strength gains "
               "and muscle hypertrophy response to training",
        evidence="Moderate",
    ),
    "rs1799983": FitnessSNP(
        rsid="rs1799983", gene="NOS3", dimension="Power/Sprint",
        trait="Nitric oxide production & blood flow",
        risk_allele="T", favorable_allele="G",
        effect="G allele (Glu298) supports normal nitric oxide synthase activity, "
               "improving blood flow and power output",
        evidence="Moderate",
    ),

    # ── ENDURANCE ──
    "rs4341": FitnessSNP(
        rsid="rs4341", gene="ACE", dimension="Endurance",
        trait="ACE insertion/deletion (I/D proxy)",
        risk_allele="G", favorable_allele="A",
        effect="I allele (A) associated with endurance performance, "
               "lower ACE activity, and improved metabolic efficiency",
        evidence="Strong",
    ),
    "rs8192678": FitnessSNP(
        rsid="rs8192678", gene="PPARGC1A", dimension="Endurance",
        trait="Mitochondrial biogenesis (PGC-1alpha)",
        risk_allele="A", favorable_allele="G",
        effect="Gly482 (G allele) associated with higher VO2max and "
               "enhanced mitochondrial function",
        evidence="Strong",
    ),
    "rs4253778": FitnessSNP(
        rsid="rs4253778", gene="PPARA", dimension="Endurance",
        trait="Fatty acid oxidation capacity",
        risk_allele="C", favorable_allele="G",
        effect="G allele linked to enhanced fat oxidation and endurance "
               "capacity through PPAR-alpha pathway",
        evidence="Moderate",
    ),
    "rs2010963": FitnessSNP(
        rsid="rs2010963", gene="VEGFA", dimension="Endurance",
        trait="Angiogenesis & capillary density",
        risk_allele="C", favorable_allele="G",
        effect="G allele associated with higher VEGF expression, "
               "improved angiogenesis and oxygen delivery to muscles",
        evidence="Moderate",
    ),
    "rs1042713": FitnessSNP(
        rsid="rs1042713", gene="ADRB2", dimension="Endurance",
        trait="Beta-2 adrenergic receptor (bronchodilation)",
        risk_allele="A", favorable_allele="G",
        effect="Gly16 (G allele) associated with improved bronchodilation "
               "and aerobic capacity under exercise",
        evidence="Moderate",
    ),
    "rs11549465": FitnessSNP(
        rsid="rs11549465", gene="HIF1A", dimension="Endurance",
        trait="Hypoxia-inducible factor (oxygen sensing)",
        risk_allele="C", favorable_allele="T",
        effect="Pro582Ser (T allele) associated with improved oxygen "
               "sensing and aerobic performance at altitude",
        evidence="Moderate",
    ),
    "rs1867785": FitnessSNP(
        rsid="rs1867785", gene="EPAS1", dimension="Endurance",
        trait="Endurance performance factor (HIF-2alpha)",
        risk_allele="A", favorable_allele="G",
        effect="G allele linked to favorable endurance phenotypes "
               "through HIF-2alpha signaling pathway",
        evidence="Preliminary",
    ),

    # ── RECOVERY ──
    "rs1800795": FitnessSNP(
        rsid="rs1800795", gene="IL6", dimension="Recovery",
        trait="Post-exercise inflammation (Interleukin-6)",
        risk_allele="G", favorable_allele="C",
        effect="C allele (-174C) associated with lower basal IL-6 levels "
               "and reduced post-exercise inflammation",
        evidence="Strong",
    ),
    "rs1205": FitnessSNP(
        rsid="rs1205", gene="CRP", dimension="Recovery",
        trait="Baseline inflammation (C-reactive protein)",
        risk_allele="C", favorable_allele="T",
        effect="T allele associated with lower basal CRP levels "
               "and faster recovery from intense exercise",
        evidence="Moderate",
    ),
    "rs1800629": FitnessSNP(
        rsid="rs1800629", gene="TNF", dimension="Recovery",
        trait="TNF-alpha inflammatory response",
        risk_allele="A", favorable_allele="G",
        effect="G allele (-308G) associated with lower TNF-alpha levels "
               "and reduced exercise-induced muscle damage",
        evidence="Strong",
    ),
    "rs4880": FitnessSNP(
        rsid="rs4880", gene="SOD2", dimension="Recovery",
        trait="Oxidative stress defense (superoxide dismutase)",
        risk_allele="T", favorable_allele="C",
        effect="Ala16 (C allele) variant targets mitochondria efficiently, "
               "providing better protection against exercise-induced oxidative stress",
        evidence="Moderate",
    ),
    "rs1800896": FitnessSNP(
        rsid="rs1800896", gene="IL10", dimension="Recovery",
        trait="Anti-inflammatory response (Interleukin-10)",
        risk_allele="T", favorable_allele="C",
        effect="C allele (-1082G) associated with higher IL-10 expression "
               "and enhanced anti-inflammatory recovery response",
        evidence="Moderate",
    ),

    # ── INJURY RISK ──
    "rs12722": FitnessSNP(
        rsid="rs12722", gene="COL5A1", dimension="Injury Risk",
        trait="Tendon & ligament collagen structure",
        risk_allele="T", favorable_allele="C",
        effect="C allele associated with greater range of motion and "
               "reduced risk of Achilles tendinopathy and ACL injuries",
        evidence="Strong",
    ),
    "rs1800012": FitnessSNP(
        rsid="rs1800012", gene="COL1A1", dimension="Injury Risk",
        trait="Type I collagen & bone/tendon strength",
        risk_allele="T", favorable_allele="G",
        effect="G allele (Sp1 site) associated with normal collagen production "
               "and reduced risk of stress fractures and tendon injuries",
        evidence="Strong",
    ),
    "rs679620": FitnessSNP(
        rsid="rs679620", gene="MMP3", dimension="Injury Risk",
        trait="Extracellular matrix remodeling",
        risk_allele="A", favorable_allele="G",
        effect="G allele associated with balanced MMP-3 activity "
               "and reduced tendon degeneration risk",
        evidence="Moderate",
    ),
    "rs143383": FitnessSNP(
        rsid="rs143383", gene="GDF5", dimension="Injury Risk",
        trait="Joint development & osteoarthritis risk",
        risk_allele="T", favorable_allele="C",
        effect="C allele associated with higher GDF-5 expression, "
               "better joint cartilage maintenance, and lower osteoarthritis risk",
        evidence="Strong",
    ),
    "rs4789932": FitnessSNP(
        rsid="rs4789932", gene="TIMP2", dimension="Injury Risk",
        trait="Tissue inhibitor of metalloproteinases",
        risk_allele="C", favorable_allele="T",
        effect="T allele associated with balanced tissue remodeling "
               "and reduced susceptibility to tendon injuries",
        evidence="Preliminary",
    ),

    # ── ADDITIONAL VARIANTS ──
    "rs1049434": FitnessSNP(
        rsid="rs1049434", gene="MCT1", dimension="Power/Sprint",
        trait="Lactate transport capacity",
        risk_allele="T", favorable_allele="A",
        effect="A allele (Asp490) supports efficient lactate clearance "
               "from muscles during high-intensity exercise",
        evidence="Moderate",
    ),
    "rs7460": FitnessSNP(
        rsid="rs7460", gene="BDKRB2", dimension="Endurance",
        trait="Bradykinin receptor (vascular efficiency)",
        risk_allele="C", favorable_allele="T",
        effect="T allele (-9 variant) associated with higher bradykinin "
               "receptor expression and improved endurance capacity",
        evidence="Preliminary",
    ),
    "rs1800169": FitnessSNP(
        rsid="rs1800169", gene="CNTF", dimension="Recovery",
        trait="Ciliary neurotrophic factor (muscle regeneration)",
        risk_allele="A", favorable_allele="G",
        effect="G allele supports normal CNTF expression, aiding "
               "motor neuron survival and muscle regeneration after exercise",
        evidence="Preliminary",
    ),
    "rs970547": FitnessSNP(
        rsid="rs970547", gene="COL12A1", dimension="Injury Risk",
        trait="Type XII collagen (ligament integrity)",
        risk_allele="T", favorable_allele="C",
        effect="C allele associated with normal collagen XII function "
               "and reduced susceptibility to ACL rupture",
        evidence="Moderate",
    ),
}


# ─────────────────────────────────────────────────────────────────────────────
# ANALYSIS RESULT DATACLASSES
# ─────────────────────────────────────────────────────────────────────────────

@dataclass
class FitnessVariantResult:
    rsid: str
    gene: str
    dimension: str
    trait: str
    genotype: str
    favorable_allele_count: int  # 0, 1, or 2
    effect: str
    interpretation: str
    evidence: str


@dataclass
class FitnessAnalysisResult:
    total_checked: int
    found: int
    not_found: int
    dimension_scores: Dict[str, float]          # dimension -> 0-100 score
    dimension_max_scores: Dict[str, float]       # dimension -> max possible
    results_by_dimension: Dict[str, List[FitnessVariantResult]]
    missing_rsids: List[str]


# ─────────────────────────────────────────────────────────────────────────────
# CORE ANALYSIS
# ─────────────────────────────────────────────────────────────────────────────

def analyze_fitness(
    variants: Dict[str, Tuple[str, str, str]]
) -> FitnessAnalysisResult:
    """
    Analyze genetic variants for fitness-related traits.

    Args:
        variants: Dictionary mapping rsID (lowercase) to (chromosome, position, genotype)

    Returns:
        FitnessAnalysisResult with dimension scores and per-variant results
    """
    results_by_dimension: Dict[str, List[FitnessVariantResult]] = {
        dim: [] for dim in DIMENSIONS
    }
    missing_rsids: List[str] = []

    # Count favorable alleles per dimension for scoring
    favorable_counts: Dict[str, int] = {dim: 0 for dim in DIMENSIONS}
    max_possible: Dict[str, int] = {dim: 0 for dim in DIMENSIONS}

    found = 0
    not_found = 0

    for rsid, snp_info in FITNESS_VARIANTS.items():
        rsid_lower = rsid.lower()
        max_possible[snp_info.dimension] += 2  # max 2 favorable alleles per SNP

        if rsid_lower in variants:
            found += 1
            _chrom, _pos, genotype = variants[rsid_lower]
            genotype = genotype.upper().replace("-", "")

            fav_count = count_allele(genotype, snp_info.favorable_allele)
            favorable_counts[snp_info.dimension] += fav_count

            interpretation = _build_interpretation(genotype, snp_info, fav_count)

            results_by_dimension[snp_info.dimension].append(
                FitnessVariantResult(
                    rsid=rsid,
                    gene=snp_info.gene,
                    dimension=snp_info.dimension,
                    trait=snp_info.trait,
                    genotype=genotype,
                    favorable_allele_count=fav_count,
                    effect=snp_info.effect if fav_count > 0 else (
                        f"No {snp_info.favorable_allele} allele detected; "
                        f"typical or less favorable variant for this trait"
                    ),
                    interpretation=interpretation,
                    evidence=snp_info.evidence,
                )
            )
        else:
            not_found += 1
            missing_rsids.append(f"{rsid} ({snp_info.gene})")

    # Calculate dimension scores (0-100)
    dimension_scores: Dict[str, float] = {}
    dimension_max: Dict[str, float] = {}
    for dim in DIMENSIONS:
        mp = max_possible[dim]
        if mp > 0:
            dimension_scores[dim] = round((favorable_counts[dim] / mp) * 100, 1)
        else:
            dimension_scores[dim] = 0.0
        dimension_max[dim] = float(mp)

    return FitnessAnalysisResult(
        total_checked=len(FITNESS_VARIANTS),
        found=found,
        not_found=not_found,
        dimension_scores=dimension_scores,
        dimension_max_scores=dimension_max,
        results_by_dimension=results_by_dimension,
        missing_rsids=missing_rsids,
    )


def _build_interpretation(genotype: str, snp: FitnessSNP, fav_count: int) -> str:
    """Build a human-readable interpretation for a single variant result."""
    fav = snp.favorable_allele.upper()
    risk = snp.risk_allele.upper()

    if fav_count == 2:
        return (
            f"You carry two copies of the favorable {fav} allele (homozygous). "
            f"This is the most advantageous genotype for {snp.trait.lower()}."
        )
    elif fav_count == 1:
        return (
            f"You carry one copy of the favorable {fav} allele (heterozygous). "
            f"You have an intermediate genetic predisposition for {snp.trait.lower()}."
        )
    else:
        return (
            f"You carry two copies of the {risk} allele. "
            f"You may have a less favorable genetic predisposition for {snp.trait.lower()}, "
            f"though training and environment play major roles."
        )


# ─────────────────────────────────────────────────────────────────────────────
# JSON OUTPUT
# ─────────────────────────────────────────────────────────────────────────────

def generate_fitness_json(result: FitnessAnalysisResult) -> dict:
    """
    Generate a JSON-serializable dict for the frontend fitness report.

    Args:
        result: FitnessAnalysisResult from analyze_fitness()

    Returns:
        Dict matching the fitness report JSON schema
    """
    dimensions_list = []

    for dim in DIMENSIONS:
        findings = result.results_by_dimension.get(dim, [])
        variants_json = []

        for f in findings:
            variants_json.append({
                "rsid": f.rsid,
                "gene": f.gene,
                "trait": f.trait,
                "genotype": f.genotype,
                "favorableAlleleCount": f.favorable_allele_count,
                "effect": f.effect,
                "interpretation": f.interpretation,
                "evidence": f.evidence,
            })

        dimensions_list.append({
            "name": dim,
            "emoji": DIMENSION_EMOJI.get(dim, ""),
            "score": result.dimension_scores.get(dim, 0.0),
            "maxScore": result.dimension_max_scores.get(dim, 0.0),
            "variants": variants_json,
        })

    return {
        "summary": {
            "totalChecked": result.total_checked,
            "found": result.found,
            "notFound": result.not_found,
            "powerScore": result.dimension_scores.get("Power/Sprint", 0.0),
            "enduranceScore": result.dimension_scores.get("Endurance", 0.0),
            "recoveryScore": result.dimension_scores.get("Recovery", 0.0),
            "injuryRiskScore": result.dimension_scores.get("Injury Risk", 0.0),
        },
        "dimensions": dimensions_list,
    }


# ─────────────────────────────────────────────────────────────────────────────
# MARKDOWN REPORT (optional, for compatibility)
# ─────────────────────────────────────────────────────────────────────────────

def generate_fitness_report(result: FitnessAnalysisResult, subject_name: str = "Subject") -> str:
    """Generate a Markdown fitness report."""
    from datetime import datetime

    report = f"""# Fitness & Athletic Performance Report

**Subject:** {subject_name}
**Generated:** {datetime.now().strftime("%B %d, %Y")}

---

## Overview

This report analyzes your genetic variants related to athletic performance, exercise recovery, and injury susceptibility. Genetics contributes to athletic potential, but training, nutrition, and lifestyle are equally important.

| Metric | Value |
|--------|-------|
| Total Fitness SNPs Checked | {result.total_checked} |
| SNPs Found in Your Data | {result.found} |
| SNPs Not Available | {result.not_found} |

### Dimension Scores

| Dimension | Score |
|-----------|-------|
"""

    for dim in DIMENSIONS:
        score = result.dimension_scores.get(dim, 0.0)
        emoji = DIMENSION_EMOJI.get(dim, "")
        report += f"| {emoji} {dim} | {score}/100 |\n"

    report += "\n---\n\n"

    for dim in DIMENSIONS:
        findings = result.results_by_dimension.get(dim, [])
        if not findings:
            continue

        emoji = DIMENSION_EMOJI.get(dim, "")
        report += f"## {emoji} {dim}\n\n"

        for f in findings:
            fav_icon = {0: "\u26aa", 1: "\U0001f535", 2: "\U0001f7e2"}.get(
                f.favorable_allele_count, "\u26aa"
            )
            report += f"{fav_icon} **{f.gene}** ({f.rsid})\n\n"
            report += f"- **Genotype:** {f.genotype}\n"
            report += f"- **Trait:** {f.trait}\n"
            report += f"- **Interpretation:** {f.interpretation}\n"
            report += f"- **Evidence:** {f.evidence}\n\n"

        report += "---\n\n"

    report += """## Important Notes

1. **Genetics is one factor** - Athletic performance depends heavily on training, nutrition, sleep, and mental resilience. Genetic predisposition does not determine outcomes.

2. **Scores are relative** - Dimension scores reflect the proportion of favorable alleles detected, not absolute athletic ability.

3. **Injury risk is probabilistic** - A higher injury risk score suggests genetic susceptibility but does not predict injuries. Proper warm-up, technique, and progressive loading remain the best prevention.

4. **Consult professionals** - For personalized training programs, work with qualified coaches and sports medicine professionals.

---
*Generated by GeneHealth Analysis Platform*
"""

    return report
