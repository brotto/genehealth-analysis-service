"""
Big Five Personality Analyzer
Analyzes genetic variants related to the Big Five (OCEAN) personality dimensions:
Openness, Conscientiousness, Extraversion, Agreeableness, and Neuroticism.

Unlike health-focused reports, each SNP here is scored DIMENSIONALLY: alleles push
a trait HIGHER or LOWER on its axis rather than being "favorable" or "risk".

References:
- Munafo et al. (2008) Biol Psychiatry - DRD4 and novelty seeking
- Ebstein et al. (1996) Nat Genet - DRD4 7R and novelty seeking
- Boutwell et al. (2017) Personal Individ Differ - CADM2 and risk-taking
- Day et al. (2016) Nat Commun - CADM2 GWAS risk behaviors
- Smeland et al. (2017) Mol Psychiatry - DRD2 rs1076560 and cognition
- Cubells & Zabetian (2004) Psychopharmacology - DBH and extraversion
- Blum et al. (1996) Pharmacogenetics - DRD2/ANKK1 Taq1A
- Kosfeld et al. (2005) Nature - OXTR and trust
- Rodrigues et al. (2009) PNAS - OXTR rs53576 empathy
- Walum et al. (2008) PNAS - AVPR1A and pair-bonding
- Caspi et al. (2003) Science - 5-HTTLPR x stress
- Sen et al. (2004) Am J Med Genet - BDNF Val66Met and neuroticism
- Stein et al. (2005) Psychol Med - COMT Val158Met and neuroticism
- Terracciano et al. (2010) Mol Psychiatry - GWAS of Big Five
- de Moor et al. (2012) Mol Psychiatry - GWAS of personality
- Hashimoto et al. (2006) Hum Mol Genet - DISC1 and traits
- Enoch et al. (2006) Am J Med Genet - GABRA2 and anxiety/neuroticism
- Zhou et al. (2008) Nature - NPY and emotional processing
- Juhasz et al. (2009) Neuropsychopharmacology - CNR1 and personality
- Gatt et al. (2009) Mol Psychiatry - BDNF/5-HTT neuroticism interactions
- Perlis et al. (2007) Mol Psychiatry - CREB1 and emotional traits
"""

from typing import Dict, List, Tuple, Any
from dataclasses import dataclass, field


# -----------------------------------------------------------------------------
# COMPLEMENT MAP
# -----------------------------------------------------------------------------

COMPLEMENT = {"A": "T", "T": "A", "C": "G", "G": "C"}


def complement_allele(allele: str) -> str:
    """Return the complement of a single nucleotide."""
    return COMPLEMENT.get(allele.upper(), allele.upper())


def count_allele(genotype: str, allele: str) -> int:
    """Count occurrences of an allele, handling complement strands."""
    genotype = genotype.upper().replace("-", "")
    allele = allele.upper()
    comp = complement_allele(allele)

    count = genotype.count(allele)
    if count == 0 and comp != allele:
        count = genotype.count(comp)

    return min(count, 2)


# -----------------------------------------------------------------------------
# DIMENSION DEFINITIONS
# -----------------------------------------------------------------------------

DIMENSIONS = [
    "Openness",
    "Conscientiousness",
    "Extraversion",
    "Agreeableness",
    "Neuroticism",
]

DIMENSION_EMOJI = {
    "Openness": "\U0001f3a8",           # artist palette
    "Conscientiousness": "\U0001f4cb",  # clipboard
    "Extraversion": "\U0001f389",       # party popper
    "Agreeableness": "\U0001f91d",      # handshake
    "Neuroticism": "\U0001f30a",        # water wave
}


# -----------------------------------------------------------------------------
# BIG FIVE VARIANT DATABASE
# -----------------------------------------------------------------------------

@dataclass
class BigFiveSNP:
    rsid: str
    gene: str
    dimension: str
    trait: str
    higher_allele: str      # allele that pushes the trait HIGHER on the dimension
    lower_allele: str       # allele that pushes the trait LOWER on the dimension
    effect_higher: str      # description when higher allele is present
    effect_lower: str       # description when lower allele is present
    evidence: str           # "Strong", "Moderate", "Preliminary"


BIG_FIVE_VARIANTS: Dict[str, BigFiveSNP] = {
    # ---- OPENNESS ----
    "rs1800955": BigFiveSNP(
        rsid="rs1800955", gene="DRD4", dimension="Openness",
        trait="Dopamine D4 receptor promoter - novelty seeking",
        higher_allele="T", lower_allele="C",
        effect_higher="T allele associated with higher novelty seeking, "
                      "curiosity, and openness to new experiences",
        effect_lower="C allele associated with more conventional preferences "
                     "and lower novelty-seeking tendencies",
        evidence="Moderate",
    ),
    "rs3758653": BigFiveSNP(
        rsid="rs3758653", gene="DRD4", dimension="Openness",
        trait="DRD4 receptor variant - exploration behavior",
        higher_allele="C", lower_allele="T",
        effect_higher="C allele linked to greater openness and exploratory "
                      "behavior through dopaminergic pathways",
        effect_lower="T allele associated with more cautious cognitive "
                     "exploration patterns",
        evidence="Preliminary",
    ),
    "rs17518584": BigFiveSNP(
        rsid="rs17518584", gene="CADM2", dimension="Openness",
        trait="Cell adhesion molecule 2 - risk-taking and exploration",
        higher_allele="T", lower_allele="C",
        effect_higher="T allele associated with higher risk-taking propensity "
                      "and openness to new experiences (GWAS-significant)",
        effect_lower="C allele linked to lower risk-taking and more cautious "
                     "openness profile",
        evidence="Strong",
    ),
    "rs4680_opn": BigFiveSNP(
        rsid="rs4818", gene="COMT", dimension="Openness",
        trait="COMT haplotype - prefrontal cognitive flexibility",
        higher_allele="C", lower_allele="G",
        effect_higher="C allele associated with higher cognitive flexibility "
                      "and openness to abstract ideas",
        effect_lower="G allele associated with more focused, less exploratory "
                     "cognitive style",
        evidence="Moderate",
    ),
    "rs6265_opn": BigFiveSNP(
        rsid="rs11030104", gene="BDNF", dimension="Openness",
        trait="BDNF expression proxy - neuroplasticity and learning openness",
        higher_allele="A", lower_allele="G",
        effect_higher="A allele linked to higher trait openness and intellectual "
                      "curiosity in GWAS of personality",
        effect_lower="G allele associated with lower trait openness scores",
        evidence="Preliminary",
    ),

    # ---- CONSCIENTIOUSNESS ----
    "rs1076560": BigFiveSNP(
        rsid="rs1076560", gene="DRD2", dimension="Conscientiousness",
        trait="DRD2 splicing - self-regulation and impulse control",
        higher_allele="G", lower_allele="T",
        effect_higher="G allele associated with balanced striatal D2 signaling "
                      "supporting self-control and conscientious behavior",
        effect_lower="T allele linked to altered D2 signaling and lower "
                     "impulse-control tendencies",
        evidence="Moderate",
    ),
    "rs1800497_c": BigFiveSNP(
        rsid="rs1800497", gene="DRD2/ANKK1", dimension="Conscientiousness",
        trait="DRD2 receptor density - reward sensitivity and planning",
        higher_allele="G", lower_allele="A",
        effect_higher="G allele (A2) associated with normal D2 density and "
                      "better long-term planning / lower impulsivity",
        effect_lower="A allele (A1) linked to reduced D2 density and higher "
                     "reward-driven impulsivity",
        evidence="Moderate",
    ),
    "rs806377": BigFiveSNP(
        rsid="rs806377", gene="CNR1", dimension="Conscientiousness",
        trait="Cannabinoid receptor 1 - self-regulation",
        higher_allele="C", lower_allele="T",
        effect_higher="C allele associated with steadier endocannabinoid tone "
                      "supporting self-regulation and persistence",
        effect_lower="T allele linked to more variable self-regulation",
        evidence="Preliminary",
    ),
    "rs2760118": BigFiveSNP(
        rsid="rs2760118", gene="CREB1", dimension="Conscientiousness",
        trait="CREB1 - learning, memory, goal persistence",
        higher_allele="C", lower_allele="T",
        effect_higher="C allele associated with robust CREB-mediated learning "
                      "and goal-directed persistence",
        effect_lower="T allele linked to lower persistence scores in some "
                     "personality studies",
        evidence="Preliminary",
    ),
    "rs4938723": BigFiveSNP(
        rsid="rs4938723", gene="NCAM1", dimension="Conscientiousness",
        trait="Neural cell adhesion molecule - executive function",
        higher_allele="C", lower_allele="T",
        effect_higher="C allele linked to higher conscientiousness in GWAS "
                      "of personality traits",
        effect_lower="T allele linked to lower conscientiousness scores",
        evidence="Preliminary",
    ),

    # ---- EXTRAVERSION ----
    "rs1611115": BigFiveSNP(
        rsid="rs1611115", gene="DBH", dimension="Extraversion",
        trait="Dopamine beta-hydroxylase - norepinephrine synthesis",
        higher_allele="C", lower_allele="T",
        effect_higher="C allele associated with higher DBH activity, greater "
                      "arousal, sociability and extraversion",
        effect_lower="T allele linked to lower DBH activity and more "
                     "introverted profile",
        evidence="Moderate",
    ),
    "rs1800497_e": BigFiveSNP(
        rsid="rs2283265", gene="DRD2", dimension="Extraversion",
        trait="DRD2 intronic variant - reward sensitivity",
        higher_allele="G", lower_allele="T",
        effect_higher="G allele associated with higher reward sensitivity and "
                      "social approach behaviors typical of extraversion",
        effect_lower="T allele linked to lower social reward sensitivity",
        evidence="Moderate",
    ),
    "rs7412746": BigFiveSNP(
        rsid="rs7412746", gene="MAOA", dimension="Extraversion",
        trait="MAOA - monoamine turnover and assertiveness",
        higher_allele="C", lower_allele="T",
        effect_higher="C allele linked to higher assertiveness and social "
                      "dominance aspects of extraversion",
        effect_lower="T allele associated with less assertive social style",
        evidence="Preliminary",
    ),
    "rs4532": BigFiveSNP(
        rsid="rs4532", gene="DRD1", dimension="Extraversion",
        trait="Dopamine D1 receptor - positive affect",
        higher_allele="C", lower_allele="T",
        effect_higher="C allele associated with elevated positive affect and "
                      "sociability",
        effect_lower="T allele linked to lower trait positive affect",
        evidence="Preliminary",
    ),
    "rs6277": BigFiveSNP(
        rsid="rs6277", gene="DRD2", dimension="Extraversion",
        trait="DRD2 C957T - striatal dopamine signaling",
        higher_allele="T", lower_allele="C",
        effect_higher="T allele linked to higher striatal D2 availability and "
                      "stronger reward/approach tendencies",
        effect_lower="C allele associated with lower reward approach",
        evidence="Moderate",
    ),

    # ---- AGREEABLENESS ----
    "rs53576": BigFiveSNP(
        rsid="rs53576", gene="OXTR", dimension="Agreeableness",
        trait="Oxytocin receptor - empathy and prosociality",
        higher_allele="G", lower_allele="A",
        effect_higher="G allele associated with higher empathy, trust, and "
                      "prosocial agreeable behavior",
        effect_lower="A allele associated with lower social sensitivity and "
                     "more socially distant style",
        evidence="Strong",
    ),
    "rs2254298": BigFiveSNP(
        rsid="rs2254298", gene="OXTR", dimension="Agreeableness",
        trait="OXTR intronic variant - social bonding",
        higher_allele="G", lower_allele="A",
        effect_higher="G allele linked to stronger attachment behavior and "
                      "agreeable interpersonal style",
        effect_lower="A allele associated with lower attachment security",
        evidence="Moderate",
    ),
    "rs3": BigFiveSNP(
        rsid="rs11174811", gene="AVPR1A", dimension="Agreeableness",
        trait="Arginine vasopressin receptor 1A - pair bonding",
        higher_allele="A", lower_allele="C",
        effect_higher="A allele associated with higher pair-bonding, altruism "
                      "and prosocial behavior",
        effect_lower="C allele linked to lower altruistic tendencies",
        evidence="Moderate",
    ),
    "rs237887": BigFiveSNP(
        rsid="rs237887", gene="OXTR", dimension="Agreeableness",
        trait="OXTR variant - face recognition and trust",
        higher_allele="A", lower_allele="G",
        effect_higher="A allele linked to better face recognition and higher "
                      "interpersonal trust",
        effect_lower="G allele associated with lower social recognition",
        evidence="Preliminary",
    ),
    "rs6323_a": BigFiveSNP(
        rsid="rs6323", gene="MAOA", dimension="Agreeableness",
        trait="MAOA activity - aggression regulation",
        higher_allele="G", lower_allele="T",
        effect_higher="G allele (high-activity MAOA) associated with calmer "
                      "temperament and higher agreeableness",
        effect_lower="T allele (low-activity MAOA) linked to higher reactive "
                     "aggression and lower agreeableness",
        evidence="Moderate",
    ),

    # ---- NEUROTICISM ----
    "rs25531": BigFiveSNP(
        rsid="rs25531", gene="SLC6A4", dimension="Neuroticism",
        trait="5-HTTLPR modifier - emotional reactivity",
        higher_allele="G", lower_allele="A",
        effect_higher="G allele (in LG haplotype) associated with lower "
                      "serotonin transporter expression and higher neuroticism",
        effect_lower="A allele associated with higher 5-HTT expression and "
                     "lower trait neuroticism",
        evidence="Strong",
    ),
    "rs4795541": BigFiveSNP(
        rsid="rs4795541", gene="SLC6A4", dimension="Neuroticism",
        trait="5-HTTLPR short/long - stress sensitivity",
        higher_allele="C", lower_allele="T",
        effect_higher="C (short) allele associated with higher amygdala "
                      "reactivity and trait neuroticism",
        effect_lower="T (long) allele linked to lower neuroticism scores",
        evidence="Strong",
    ),
    "rs6265": BigFiveSNP(
        rsid="rs6265", gene="BDNF", dimension="Neuroticism",
        trait="BDNF Val66Met - stress vulnerability",
        higher_allele="A", lower_allele="G",
        effect_higher="A (Met) allele associated with altered BDNF secretion "
                      "and higher trait neuroticism",
        effect_lower="G (Val) allele linked to lower neuroticism and greater "
                     "emotional stability",
        evidence="Strong",
    ),
    "rs4680": BigFiveSNP(
        rsid="rs4680", gene="COMT", dimension="Neuroticism",
        trait="COMT Val158Met - emotional pain sensitivity",
        higher_allele="A", lower_allele="G",
        effect_higher="A (Met) allele ('worrier' phenotype) associated with "
                      "higher emotional reactivity and neuroticism",
        effect_lower="G (Val) allele ('warrior' phenotype) linked to lower "
                     "neuroticism and better stress coping",
        evidence="Strong",
    ),
    "rs1360780": BigFiveSNP(
        rsid="rs1360780", gene="FKBP5", dimension="Neuroticism",
        trait="HPA axis feedback - stress sensitivity",
        higher_allele="T", lower_allele="C",
        effect_higher="T allele associated with prolonged cortisol response "
                      "and higher trait neuroticism",
        effect_lower="C allele linked to faster stress recovery and lower "
                     "neuroticism",
        evidence="Strong",
    ),
}


# -----------------------------------------------------------------------------
# ANALYSIS RESULT DATACLASSES
# -----------------------------------------------------------------------------

@dataclass
class BigFiveVariantResult:
    rsid: str
    gene: str
    dimension: str
    trait: str
    genotype: str
    direction: str          # "higher", "lower", or "neutral"
    higher_allele_count: int
    effect: str
    interpretation: str
    evidence: str


@dataclass
class BigFiveAnalysisResult:
    total_checked: int
    found: int
    not_found: int
    dimension_scores: Dict[str, float]          # 0-100 per OCEAN dimension
    dimension_max_scores: Dict[str, float]
    results_by_dimension: Dict[str, List[BigFiveVariantResult]]
    missing_rsids: List[str]


# -----------------------------------------------------------------------------
# CORE ANALYSIS
# -----------------------------------------------------------------------------

def analyze_big_five(
    variants: Dict[str, Tuple[str, str, str]]
) -> BigFiveAnalysisResult:
    """
    Analyze genetic variants for Big Five personality dimensions.

    Each SNP contributes 0, 1, or 2 "higher" alleles. Dimension scores are
    normalized to 0-100 where 50 represents a neutral/balanced profile and
    higher values indicate a genetic push toward the HIGHER end of that trait.
    """
    results_by_dimension: Dict[str, List[BigFiveVariantResult]] = {
        dim: [] for dim in DIMENSIONS
    }
    missing_rsids: List[str] = []

    higher_counts: Dict[str, int] = {dim: 0 for dim in DIMENSIONS}
    max_possible: Dict[str, int] = {dim: 0 for dim in DIMENSIONS}

    found = 0
    not_found = 0

    for key, snp_info in BIG_FIVE_VARIANTS.items():
        # The dict key is a unique identifier; the real rsID lives on the object
        rsid = snp_info.rsid
        rsid_lower = rsid.lower()
        max_possible[snp_info.dimension] += 2

        if rsid_lower in variants:
            found += 1
            _chrom, _pos, genotype = variants[rsid_lower]
            genotype = genotype.upper().replace("-", "")

            higher_count = count_allele(genotype, snp_info.higher_allele)
            higher_counts[snp_info.dimension] += higher_count

            # direction classification
            if higher_count == 2:
                direction = "higher"
                effect_text = snp_info.effect_higher
            elif higher_count == 0:
                direction = "lower"
                effect_text = snp_info.effect_lower
            else:
                direction = "neutral"
                effect_text = (
                    f"Heterozygous: carries one {snp_info.higher_allele} and "
                    f"one {snp_info.lower_allele} allele, producing an "
                    f"intermediate effect on {snp_info.trait.lower()}."
                )

            interpretation = _build_interpretation(
                genotype, snp_info, higher_count, direction
            )

            results_by_dimension[snp_info.dimension].append(
                BigFiveVariantResult(
                    rsid=rsid,
                    gene=snp_info.gene,
                    dimension=snp_info.dimension,
                    trait=snp_info.trait,
                    genotype=genotype,
                    direction=direction,
                    higher_allele_count=higher_count,
                    effect=effect_text,
                    interpretation=interpretation,
                    evidence=snp_info.evidence,
                )
            )
        else:
            not_found += 1
            missing_rsids.append(f"{rsid} ({snp_info.gene})")

    # Dimension scores: 0-100, where 50 = neutral baseline
    dimension_scores: Dict[str, float] = {}
    dimension_max: Dict[str, float] = {}
    for dim in DIMENSIONS:
        mp = max_possible[dim]
        if mp > 0:
            dimension_scores[dim] = round((higher_counts[dim] / mp) * 100, 1)
        else:
            dimension_scores[dim] = 50.0  # neutral default if no data
        dimension_max[dim] = float(mp)

    return BigFiveAnalysisResult(
        total_checked=len(BIG_FIVE_VARIANTS),
        found=found,
        not_found=not_found,
        dimension_scores=dimension_scores,
        dimension_max_scores=dimension_max,
        results_by_dimension=results_by_dimension,
        missing_rsids=missing_rsids,
    )


def _build_interpretation(
    genotype: str, snp: BigFiveSNP, higher_count: int, direction: str
) -> str:
    """Build a human-readable interpretation for a single variant."""
    high = snp.higher_allele.upper()
    low = snp.lower_allele.upper()
    trait_lower = snp.trait.lower()

    if direction == "higher":
        return (
            f"You carry two copies of the {high} allele (homozygous). "
            f"This genotype is associated with a genetic push TOWARD the higher "
            f"end of {snp.dimension.lower()} via {trait_lower}."
        )
    elif direction == "lower":
        return (
            f"You carry two copies of the {low} allele (homozygous). "
            f"This genotype is associated with a genetic push TOWARD the lower "
            f"end of {snp.dimension.lower()} via {trait_lower}."
        )
    else:
        return (
            f"You are heterozygous ({high}/{low}) at this locus. "
            f"The effect on {snp.dimension.lower()} is intermediate - neither "
            f"a strong push higher nor lower for {trait_lower}."
        )


def _dimension_interpretation(dim: str, score: float) -> str:
    """High-level interpretation text for a dimension score."""
    if score >= 70:
        level = "high"
    elif score >= 55:
        level = "moderately high"
    elif score >= 45:
        level = "balanced"
    elif score >= 30:
        level = "moderately low"
    else:
        level = "low"

    descriptions = {
        "Openness": {
            "high": "Your variants suggest strong curiosity, creativity, and "
                    "openness to new ideas and experiences.",
            "moderately high": "Your variants lean toward curiosity and "
                               "willingness to explore novel experiences.",
            "balanced": "Your variants suggest a balanced mix of openness and "
                        "preference for the familiar.",
            "moderately low": "Your variants lean toward practicality and "
                              "preference for familiar routines.",
            "low": "Your variants suggest a preference for tradition, routine "
                   "and concrete thinking over novelty.",
        },
        "Conscientiousness": {
            "high": "Your variants suggest strong self-discipline, organization, "
                    "and goal-directed persistence.",
            "moderately high": "Your variants lean toward structure, planning, "
                               "and reliable follow-through.",
            "balanced": "Your variants suggest a balanced profile between "
                        "structure and flexibility.",
            "moderately low": "Your variants lean toward spontaneity and "
                              "flexible, less-planned approaches.",
            "low": "Your variants suggest a more spontaneous, less structured "
                   "and more impulsive approach to tasks.",
        },
        "Extraversion": {
            "high": "Your variants suggest high sociability, assertiveness, "
                    "and positive affect in social settings.",
            "moderately high": "Your variants lean toward sociability and "
                               "enjoyment of social activity.",
            "balanced": "Your variants suggest an ambivert profile with "
                        "flexibility between social and solitary modes.",
            "moderately low": "Your variants lean toward introversion and "
                              "preference for calm, smaller-group settings.",
            "low": "Your variants suggest a strongly introverted profile with "
                   "preference for solitude and low-stimulation environments.",
        },
        "Agreeableness": {
            "high": "Your variants suggest strong empathy, trust, and "
                    "cooperative prosocial behavior.",
            "moderately high": "Your variants lean toward warmth, cooperation, "
                               "and interpersonal sensitivity.",
            "balanced": "Your variants suggest a balanced mix of cooperation "
                        "and self-interest.",
            "moderately low": "Your variants lean toward directness and "
                              "skepticism in interpersonal relationships.",
            "low": "Your variants suggest a more competitive, skeptical, or "
                   "detached interpersonal style.",
        },
        "Neuroticism": {
            "high": "Your variants suggest higher emotional reactivity and "
                    "sensitivity to stress and negative emotions.",
            "moderately high": "Your variants lean toward greater emotional "
                               "sensitivity and stress reactivity.",
            "balanced": "Your variants suggest a balanced emotional reactivity "
                        "profile.",
            "moderately low": "Your variants lean toward emotional stability "
                              "and even-tempered responses.",
            "low": "Your variants suggest strong emotional stability and "
                   "resilience to stress.",
        },
    }

    return descriptions.get(dim, {}).get(level, "Balanced profile.")


# -----------------------------------------------------------------------------
# JSON OUTPUT
# -----------------------------------------------------------------------------

def generate_big_five_json(result: BigFiveAnalysisResult) -> dict:
    """Generate JSON-serializable dict for the frontend Big Five report."""
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
                "direction": f.direction,
                "higherAlleleCount": f.higher_allele_count,
                "effect": f.effect,
                "interpretation": f.interpretation,
                "evidence": f.evidence,
            })

        score = result.dimension_scores.get(dim, 50.0)
        dimensions_list.append({
            "name": dim,
            "emoji": DIMENSION_EMOJI.get(dim, ""),
            "score": score,
            "maxScore": result.dimension_max_scores.get(dim, 0.0),
            "interpretation": _dimension_interpretation(dim, score),
            "variants": variants_json,
        })

    return {
        "summary": {
            "totalChecked": result.total_checked,
            "found": result.found,
            "notFound": result.not_found,
            "opennessScore": result.dimension_scores.get("Openness", 50.0),
            "conscientiousnessScore": result.dimension_scores.get("Conscientiousness", 50.0),
            "extraversionScore": result.dimension_scores.get("Extraversion", 50.0),
            "agreeablenessScore": result.dimension_scores.get("Agreeableness", 50.0),
            "neuroticismScore": result.dimension_scores.get("Neuroticism", 50.0),
        },
        "dimensions": dimensions_list,
    }


# -----------------------------------------------------------------------------
# MARKDOWN REPORT (optional)
# -----------------------------------------------------------------------------

def generate_big_five_report(
    result: BigFiveAnalysisResult, subject_name: str = "Subject"
) -> str:
    """Generate a Markdown Big Five personality report."""
    from datetime import datetime

    report = f"""# Big Five Personality Genetic Report

**Subject:** {subject_name}
**Generated:** {datetime.now().strftime("%B %d, %Y")}

---

## Overview

This report analyzes your genetic variants associated with the Big Five (OCEAN) personality dimensions. Personality is shaped by a complex interaction of genetics, upbringing, culture, and life experience. Genetic variants contribute only a modest portion of individual differences in personality.

| Metric | Value |
|--------|-------|
| Total SNPs Checked | {result.total_checked} |
| SNPs Found in Your Data | {result.found} |
| SNPs Not Available | {result.not_found} |

### OCEAN Dimension Scores

| Dimension | Score |
|-----------|-------|
"""

    for dim in DIMENSIONS:
        score = result.dimension_scores.get(dim, 50.0)
        emoji = DIMENSION_EMOJI.get(dim, "")
        report += f"| {emoji} {dim} | {score}/100 |\n"

    report += "\n---\n\n"

    for dim in DIMENSIONS:
        findings = result.results_by_dimension.get(dim, [])
        if not findings:
            continue

        emoji = DIMENSION_EMOJI.get(dim, "")
        score = result.dimension_scores.get(dim, 50.0)
        interp = _dimension_interpretation(dim, score)
        report += f"## {emoji} {dim} ({score}/100)\n\n"
        report += f"{interp}\n\n"

        for f in findings:
            arrow = {"higher": "\u2b06\ufe0f", "lower": "\u2b07\ufe0f",
                     "neutral": "\u27a1\ufe0f"}.get(f.direction, "\u27a1\ufe0f")
            report += f"{arrow} **{f.gene}** ({f.rsid})\n\n"
            report += f"- **Genotype:** {f.genotype}\n"
            report += f"- **Direction:** {f.direction}\n"
            report += f"- **Trait:** {f.trait}\n"
            report += f"- **Interpretation:** {f.interpretation}\n"
            report += f"- **Evidence:** {f.evidence}\n\n"

        report += "---\n\n"

    report += """## Important Notes

1. **Genetics explains only a small portion of personality** - Twin studies estimate heritability of Big Five traits around 40-50%. Most of this is polygenic, involving thousands of variants each with tiny effects. This report samples only well-studied candidate genes.

2. **Not a personality test** - This report does not replace validated personality instruments like the NEO-PI-R or BFI-2. For an accurate personality assessment, take a validated self-report inventory.

3. **Environment and experience shape personality** - Parenting, culture, relationships, and life events play a major role in who you become.

4. **Personality is not destiny** - Traits can shift across the lifespan. Self-awareness and intentional effort can meaningfully change behavior.

---
*Generated by GeneHealth Analysis Platform*
"""

    return report
