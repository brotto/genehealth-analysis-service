"""
Mental Wellbeing Analyzer
Analyzes genetic variants related to mental health and emotional wellbeing across
four dimensions: Depression Risk, Anxiety, Stress Resilience, and Emotional Regulation.

References:
- Caspi et al. (2003) Science - 5-HTTLPR x stress interaction
- Hek et al. (2013) Biol Psychiatry - 5-HTTLPR meta-analysis
- Egan et al. (2003) Cell - BDNF Val66Met
- Chen et al. (2006) Science - BDNF Met allele and hippocampal function
- Stein et al. (2008) Biol Psychiatry - COMT Val158Met and anxiety
- Binder et al. (2008) JAMA - FKBP5 and PTSD
- Klengel et al. (2013) Nat Neurosci - FKBP5 epigenetic mechanism
- Bradley et al. (2008) Arch Gen Psychiatry - CRHR1 and depression
- Sabol et al. (1998) Hum Genet - MAOA uVNTR
- van Rossum et al. (2006) Biol Psychiatry - NR3C1 BclI
- Rodrigues et al. (2009) PNAS - OXTR rs53576 and empathy
- Gao et al. (2012) Mol Psychiatry - TPH2 meta-analysis
- Blum et al. (1996) Pharmacogenetics - DRD2 Taq1A
- Enoch et al. (2006) Am J Med Genet - GABRA2 and anxiety
- Ferreira et al. (2008) Nat Genet - CACNA1C and mood disorders
- Hashimoto et al. (2006) Hum Mol Genet - DISC1 and depression
- Zhou et al. (2008) Nature - NPY haplotype and stress response
- Hariri et al. (2009) Arch Gen Psychiatry - CNR1 and depression
- Schumacher et al. (2009) Mol Psychiatry - ANK3 and bipolar
- Domschke et al. (2011) Int J Neuropsychopharmacol - NPSR1 and panic
"""

from typing import Dict, List, Tuple, Any
from dataclasses import dataclass, field


# -----------------------------------------------------------------------------
# COMPLEMENT MAP (for strand-flip handling)
# -----------------------------------------------------------------------------

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

    if count == 0 and comp != allele:
        count = genotype.count(comp)

    return min(count, 2)


# -----------------------------------------------------------------------------
# CATEGORY DEFINITIONS
# -----------------------------------------------------------------------------

CATEGORIES = ["Depression Risk", "Anxiety", "Stress Resilience", "Emotional Regulation"]

CATEGORY_EMOJI = {
    "Depression Risk": "\U0001f327\ufe0f",      # cloud with rain
    "Anxiety": "\U0001f630",                     # anxious face with sweat
    "Stress Resilience": "\U0001f6e1\ufe0f",     # shield
    "Emotional Regulation": "\U0001f49a",        # green heart
}


# -----------------------------------------------------------------------------
# MENTAL WELLBEING VARIANT DATABASE
# -----------------------------------------------------------------------------

@dataclass
class MentalWellbeingSNP:
    rsid: str
    gene: str
    category: str
    trait: str
    risk_allele: str        # allele associated with vulnerability
    favorable_allele: str   # allele associated with resilience / protection
    effect: str             # description when favorable allele is present
    evidence: str           # "Strong", "Moderate", or "Preliminary"


MENTAL_WELLBEING_VARIANTS: Dict[str, MentalWellbeingSNP] = {
    # -- DEPRESSION RISK --
    "rs25531": MentalWellbeingSNP(
        rsid="rs25531", gene="SLC6A4", category="Depression Risk",
        trait="Serotonin transporter expression (5-HTTLPR modifier)",
        risk_allele="G", favorable_allele="A",
        effect="A allele associated with higher serotonin transporter expression "
               "and reduced susceptibility to stress-induced depression",
        evidence="Strong",
    ),
    "rs4795541": MentalWellbeingSNP(
        rsid="rs4795541", gene="SLC6A4", category="Depression Risk",
        trait="5-HTTLPR long/short variant proxy",
        risk_allele="C", favorable_allele="T",
        effect="T (long) allele associated with higher serotonin reuptake capacity "
               "and lower risk of depression following life stressors",
        evidence="Strong",
    ),
    "rs6265": MentalWellbeingSNP(
        rsid="rs6265", gene="BDNF", category="Depression Risk",
        trait="BDNF Val66Met - neuroplasticity and mood",
        risk_allele="A", favorable_allele="G",
        effect="Val/Val (G allele) genotype supports normal BDNF secretion, "
               "hippocampal function, and resilience against depression",
        evidence="Strong",
    ),
    "rs1800497": MentalWellbeingSNP(
        rsid="rs1800497", gene="DRD2/ANKK1", category="Depression Risk",
        trait="Dopamine D2 receptor density (Taq1A)",
        risk_allele="A", favorable_allele="G",
        effect="A1 (A) allele reduces D2 receptor density; G allele supports "
               "normal dopaminergic reward signaling and lower anhedonia risk",
        evidence="Moderate",
    ),
    "rs821616": MentalWellbeingSNP(
        rsid="rs821616", gene="DISC1", category="Depression Risk",
        trait="DISC1 Ser704Cys - neurodevelopment and mood",
        risk_allele="A", favorable_allele="T",
        effect="Ser (T) allele associated with normal DISC1 function and "
               "lower risk of major depression and psychiatric disorders",
        evidence="Moderate",
    ),

    # -- ANXIETY --
    "rs4680": MentalWellbeingSNP(
        rsid="rs4680", gene="COMT", category="Anxiety",
        trait="Catechol-O-methyltransferase Val158Met",
        risk_allele="A", favorable_allele="G",
        effect="Val (G) allele provides faster dopamine clearance in prefrontal "
               "cortex, associated with lower anxiety sensitivity ('warrior' phenotype)",
        evidence="Strong",
    ),
    "rs6323": MentalWellbeingSNP(
        rsid="rs6323", gene="MAOA", category="Anxiety",
        trait="Monoamine oxidase A activity",
        risk_allele="T", favorable_allele="G",
        effect="G allele associated with higher MAOA activity, more efficient "
               "neurotransmitter metabolism and reduced trait anxiety",
        evidence="Moderate",
    ),
    "rs279858": MentalWellbeingSNP(
        rsid="rs279858", gene="GABRA2", category="Anxiety",
        trait="GABA-A receptor alpha-2 subunit",
        risk_allele="C", favorable_allele="T",
        effect="T allele associated with normal GABAergic inhibitory tone and "
               "lower generalized anxiety and alcohol-use risk",
        evidence="Moderate",
    ),
    "rs324981": MentalWellbeingSNP(
        rsid="rs324981", gene="NPSR1", category="Anxiety",
        trait="Neuropeptide S receptor 1 - arousal regulation",
        risk_allele="A", favorable_allele="T",
        effect="T allele linked to balanced arousal response and reduced "
               "risk of panic disorder and anxiety sensitivity",
        evidence="Preliminary",
    ),
    "rs1006737": MentalWellbeingSNP(
        rsid="rs1006737", gene="CACNA1C", category="Anxiety",
        trait="L-type calcium channel - mood circuit excitability",
        risk_allele="A", favorable_allele="G",
        effect="G allele associated with typical neuronal calcium signaling "
               "and lower risk of anxiety and bipolar spectrum disorders",
        evidence="Strong",
    ),

    # -- STRESS RESILIENCE --
    "rs1360780": MentalWellbeingSNP(
        rsid="rs1360780", gene="FKBP5", category="Stress Resilience",
        trait="Glucocorticoid receptor sensitivity - HPA axis",
        risk_allele="T", favorable_allele="C",
        effect="C allele supports efficient cortisol feedback and faster HPA "
               "axis recovery; lower PTSD risk after trauma",
        evidence="Strong",
    ),
    "rs3800373": MentalWellbeingSNP(
        rsid="rs3800373", gene="FKBP5", category="Stress Resilience",
        trait="FKBP5 co-chaperone expression",
        risk_allele="C", favorable_allele="A",
        effect="A allele (major) associated with typical FKBP5 expression "
               "and better stress recovery and PTSD resilience",
        evidence="Strong",
    ),
    "rs110402": MentalWellbeingSNP(
        rsid="rs110402", gene="CRHR1", category="Stress Resilience",
        trait="CRH receptor 1 - stress hormone signaling",
        risk_allele="G", favorable_allele="A",
        effect="A allele associated with protective effect against depression "
               "in individuals with childhood adversity",
        evidence="Moderate",
    ),
    "rs41423247": MentalWellbeingSNP(
        rsid="rs41423247", gene="NR3C1", category="Stress Resilience",
        trait="Glucocorticoid receptor BclI variant",
        risk_allele="C", favorable_allele="G",
        effect="G allele associated with modulated cortisol sensitivity and "
               "resilience to chronic stress effects",
        evidence="Moderate",
    ),
    "rs16147": MentalWellbeingSNP(
        rsid="rs16147", gene="NPY", category="Stress Resilience",
        trait="Neuropeptide Y expression - stress buffering",
        risk_allele="T", favorable_allele="C",
        effect="C allele associated with higher NPY expression, greater "
               "stress resilience and lower amygdala reactivity",
        evidence="Moderate",
    ),

    # -- EMOTIONAL REGULATION --
    "rs53576": MentalWellbeingSNP(
        rsid="rs53576", gene="OXTR", category="Emotional Regulation",
        trait="Oxytocin receptor - social bonding and empathy",
        risk_allele="A", favorable_allele="G",
        effect="GG genotype associated with higher empathy, better social "
               "sensitivity, and stronger emotional support seeking",
        evidence="Strong",
    ),
    "rs4570625": MentalWellbeingSNP(
        rsid="rs4570625", gene="TPH2", category="Emotional Regulation",
        trait="Tryptophan hydroxylase 2 - brain serotonin synthesis",
        risk_allele="T", favorable_allele="G",
        effect="G allele associated with normal brain serotonin synthesis "
               "and balanced amygdala response to emotional stimuli",
        evidence="Moderate",
    ),
    "rs1049353": MentalWellbeingSNP(
        rsid="rs1049353", gene="CNR1", category="Emotional Regulation",
        trait="Cannabinoid receptor 1 - emotional processing",
        risk_allele="A", favorable_allele="G",
        effect="G allele associated with typical endocannabinoid signaling "
               "and more stable emotional processing",
        evidence="Moderate",
    ),
    "rs10994336": MentalWellbeingSNP(
        rsid="rs10994336", gene="ANK3", category="Emotional Regulation",
        trait="Ankyrin 3 - neuronal excitability and mood stability",
        risk_allele="T", favorable_allele="C",
        effect="C allele associated with normal sodium channel clustering "
               "and reduced risk of mood instability",
        evidence="Moderate",
    ),
    "rs6314": MentalWellbeingSNP(
        rsid="rs6314", gene="HTR2A", category="Emotional Regulation",
        trait="Serotonin 2A receptor - emotional reactivity",
        risk_allele="A", favorable_allele="G",
        effect="G allele (His452) associated with normal 5-HT2A signaling "
               "and balanced emotional regulation",
        evidence="Preliminary",
    ),
}


# -----------------------------------------------------------------------------
# ANALYSIS RESULT DATACLASSES
# -----------------------------------------------------------------------------

@dataclass
class MentalWellbeingVariantResult:
    rsid: str
    gene: str
    category: str
    trait: str
    genotype: str
    favorable_allele_count: int  # 0, 1, or 2
    effect: str
    interpretation: str
    evidence: str


@dataclass
class MentalWellbeingAnalysisResult:
    total_checked: int
    found: int
    not_found: int
    category_scores: Dict[str, float]
    category_max_scores: Dict[str, float]
    overall_score: float
    results_by_category: Dict[str, List[MentalWellbeingVariantResult]]
    missing_rsids: List[str]


# -----------------------------------------------------------------------------
# CORE ANALYSIS
# -----------------------------------------------------------------------------

def analyze_mental_wellbeing(
    variants: Dict[str, Tuple[str, str, str]]
) -> MentalWellbeingAnalysisResult:
    """
    Analyze genetic variants for mental wellbeing traits.

    Args:
        variants: Dictionary mapping rsID (lowercase) to (chromosome, position, genotype)

    Returns:
        MentalWellbeingAnalysisResult with category scores and per-variant results
    """
    results_by_category: Dict[str, List[MentalWellbeingVariantResult]] = {
        cat: [] for cat in CATEGORIES
    }
    missing_rsids: List[str] = []

    favorable_counts: Dict[str, int] = {cat: 0 for cat in CATEGORIES}
    max_possible: Dict[str, int] = {cat: 0 for cat in CATEGORIES}

    found = 0
    not_found = 0

    for rsid, snp_info in MENTAL_WELLBEING_VARIANTS.items():
        rsid_lower = rsid.lower()
        max_possible[snp_info.category] += 2

        if rsid_lower in variants:
            found += 1
            _chrom, _pos, genotype = variants[rsid_lower]
            genotype = genotype.upper().replace("-", "")

            fav_count = count_allele(genotype, snp_info.favorable_allele)
            favorable_counts[snp_info.category] += fav_count

            interpretation = _build_interpretation(genotype, snp_info, fav_count)

            results_by_category[snp_info.category].append(
                MentalWellbeingVariantResult(
                    rsid=rsid,
                    gene=snp_info.gene,
                    category=snp_info.category,
                    trait=snp_info.trait,
                    genotype=genotype,
                    favorable_allele_count=fav_count,
                    effect=snp_info.effect if fav_count > 0 else (
                        f"No {snp_info.favorable_allele} allele detected; "
                        f"genotype may confer greater vulnerability for this trait"
                    ),
                    interpretation=interpretation,
                    evidence=snp_info.evidence,
                )
            )
        else:
            not_found += 1
            missing_rsids.append(f"{rsid} ({snp_info.gene})")

    # Category scores (0-100)
    category_scores: Dict[str, float] = {}
    category_max: Dict[str, float] = {}
    for cat in CATEGORIES:
        mp = max_possible[cat]
        if mp > 0:
            category_scores[cat] = round((favorable_counts[cat] / mp) * 100, 1)
        else:
            category_scores[cat] = 0.0
        category_max[cat] = float(mp)

    # Overall wellbeing score: average across categories that have data
    scored_cats = [category_scores[c] for c in CATEGORIES if max_possible[c] > 0]
    overall = round(sum(scored_cats) / len(scored_cats), 1) if scored_cats else 0.0

    return MentalWellbeingAnalysisResult(
        total_checked=len(MENTAL_WELLBEING_VARIANTS),
        found=found,
        not_found=not_found,
        category_scores=category_scores,
        category_max_scores=category_max,
        overall_score=overall,
        results_by_category=results_by_category,
        missing_rsids=missing_rsids,
    )


def _build_interpretation(
    genotype: str, snp: MentalWellbeingSNP, fav_count: int
) -> str:
    """Build a human-readable interpretation for a single variant result."""
    fav = snp.favorable_allele.upper()
    risk = snp.risk_allele.upper()

    if fav_count == 2:
        return (
            f"You carry two copies of the protective {fav} allele (homozygous). "
            f"This is the most favorable genotype for {snp.trait.lower()}, "
            f"suggesting greater genetic resilience in this dimension."
        )
    elif fav_count == 1:
        return (
            f"You carry one copy of the protective {fav} allele (heterozygous). "
            f"You have an intermediate genetic predisposition regarding "
            f"{snp.trait.lower()}."
        )
    else:
        return (
            f"You carry two copies of the {risk} allele. "
            f"You may have greater genetic susceptibility related to "
            f"{snp.trait.lower()}, though environment, lifestyle, and support "
            f"systems strongly modulate actual outcomes."
        )


# -----------------------------------------------------------------------------
# JSON OUTPUT
# -----------------------------------------------------------------------------

def generate_mental_wellbeing_json(result: MentalWellbeingAnalysisResult) -> dict:
    """
    Generate a JSON-serializable dict for the frontend mental wellbeing report.
    """
    categories_list = []

    for cat in CATEGORIES:
        findings = result.results_by_category.get(cat, [])
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

        categories_list.append({
            "name": cat,
            "emoji": CATEGORY_EMOJI.get(cat, ""),
            "score": result.category_scores.get(cat, 0.0),
            "maxScore": result.category_max_scores.get(cat, 0.0),
            "variants": variants_json,
        })

    return {
        "summary": {
            "totalChecked": result.total_checked,
            "found": result.found,
            "notFound": result.not_found,
            "overallScore": result.overall_score,
            "depressionRiskScore": result.category_scores.get("Depression Risk", 0.0),
            "anxietyScore": result.category_scores.get("Anxiety", 0.0),
            "stressResilienceScore": result.category_scores.get("Stress Resilience", 0.0),
            "emotionalRegulationScore": result.category_scores.get("Emotional Regulation", 0.0),
        },
        "categories": categories_list,
    }


# -----------------------------------------------------------------------------
# MARKDOWN REPORT (optional)
# -----------------------------------------------------------------------------

def generate_mental_wellbeing_report(
    result: MentalWellbeingAnalysisResult, subject_name: str = "Subject"
) -> str:
    """Generate a Markdown mental wellbeing report."""
    from datetime import datetime

    report = f"""# Mental Wellbeing Report

**Subject:** {subject_name}
**Generated:** {datetime.now().strftime("%B %d, %Y")}

---

## Overview

This report analyzes your genetic variants related to mental health and emotional wellbeing. Mental health is shaped by a complex interplay of genetics, environment, life experiences, and coping resources. Genetic predisposition is only one factor.

| Metric | Value |
|--------|-------|
| Total SNPs Checked | {result.total_checked} |
| SNPs Found in Your Data | {result.found} |
| SNPs Not Available | {result.not_found} |
| Overall Wellbeing Score | {result.overall_score}/100 |

### Category Scores

| Category | Score |
|----------|-------|
"""

    for cat in CATEGORIES:
        score = result.category_scores.get(cat, 0.0)
        emoji = CATEGORY_EMOJI.get(cat, "")
        report += f"| {emoji} {cat} | {score}/100 |\n"

    report += "\n---\n\n"

    for cat in CATEGORIES:
        findings = result.results_by_category.get(cat, [])
        if not findings:
            continue

        emoji = CATEGORY_EMOJI.get(cat, "")
        report += f"## {emoji} {cat}\n\n"

        for f in findings:
            icon = {0: "\u26aa", 1: "\U0001f535", 2: "\U0001f7e2"}.get(
                f.favorable_allele_count, "\u26aa"
            )
            report += f"{icon} **{f.gene}** ({f.rsid})\n\n"
            report += f"- **Genotype:** {f.genotype}\n"
            report += f"- **Trait:** {f.trait}\n"
            report += f"- **Interpretation:** {f.interpretation}\n"
            report += f"- **Evidence:** {f.evidence}\n\n"

        report += "---\n\n"

    report += """## Important Notes

1. **Genetics is one of many factors** - Mental health is profoundly influenced by environment, relationships, life events, and personal choices. Genetic predisposition does not determine mental health outcomes.

2. **This is not a diagnosis** - These results do not diagnose or predict mental illness. A higher susceptibility score does not mean you will develop a condition.

3. **Protective factors matter** - Social support, therapy, exercise, sleep, and stress management can significantly offset genetic vulnerabilities.

4. **Seek professional help** - If you are experiencing mental health difficulties, please consult a qualified mental health professional.

---
*Generated by GeneHealth Analysis Platform*
"""

    return report
