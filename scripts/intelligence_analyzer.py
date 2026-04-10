"""
Intelligence Analyzer
Analyzes genetic variants associated with cognitive abilities across five
dimensions: General Intelligence, Memory, Processing Speed, Verbal Ability,
and Learning & Plasticity.

IMPORTANT DISCLAIMER
====================
Intelligence is overwhelmingly shaped by environmental, educational, social,
nutritional, and health factors. The genetic contribution to any single
cognitive outcome from common SNPs is very small (typically each variant
explains <0.1% of variance). These results describe statistical TENDENCIES
observed in research populations, they do NOT determine or predict an
individual's intelligence, capabilities, or potential. Treat this report as
an educational curiosity, not a measurement of ability.

References:
- Dick et al. (2007) Behav Genet - CHRM2 and cognition
- Gosso et al. (2006) Genes Brain Behav - SNAP25
- Deary et al. (2002) Neurosci Lett - APOE and cognitive aging
- Papassotiropoulos et al. (2006) Science - KIBRA/WWC1 and memory
- Egan et al. (2001) Cell - BDNF Val66Met and memory
- Egan et al. (2003) PNAS - COMT Val158Met and working memory
- Lai et al. (2001) Nature - FOXP2 and language
- Vernes et al. (2008) NEJM - CNTNAP2 and language development
- Hill et al. (2018) Mol Psychiatry - GWAS of intelligence
- Savage et al. (2018) Nat Genet - GWAS of intelligence (248k)
- Lee et al. (2018) Nat Genet - GWAS of educational attainment
- Ferreira et al. (2008) Nat Genet - ANK3, CACNA1C
- Williams et al. (2011) Mol Psychiatry - ZNF804A
"""

from typing import Dict, List, Tuple, Any
from dataclasses import dataclass, field


# ─────────────────────────────────────────────────────────────────────────────
# IMPORTANT DISCLAIMER (included in JSON + markdown output)
# ─────────────────────────────────────────────────────────────────────────────

INTELLIGENCE_DISCLAIMER = (
    "Intelligence is MOSTLY shaped by environment, education, nutrition, sleep, "
    "health, social context, and effort. The genetic contribution from common "
    "SNPs is very small — each variant explains well under 1% of variance in "
    "cognitive outcomes. This report describes statistical TENDENCIES observed "
    "in research populations, it does NOT measure, determine, or predict your "
    "intelligence, learning capacity, or potential. Treat it as educational "
    "curiosity, not an assessment of ability."
)


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

    if count == 0 and comp != allele:
        count = genotype.count(comp)

    return min(count, 2)


# ─────────────────────────────────────────────────────────────────────────────
# CATEGORY DEFINITIONS
# ─────────────────────────────────────────────────────────────────────────────

CATEGORIES = [
    "General Intelligence",
    "Memory",
    "Processing Speed",
    "Verbal Ability",
    "Learning & Plasticity",
]

CATEGORY_EMOJI = {
    "General Intelligence": "\U0001f9e0",   # 🧠
    "Memory": "\U0001f4ad",                  # 💭
    "Processing Speed": "\u26a1",            # ⚡
    "Verbal Ability": "\U0001f4ac",          # 💬
    "Learning & Plasticity": "\U0001f4da",   # 📚
}


# ─────────────────────────────────────────────────────────────────────────────
# INTELLIGENCE VARIANT DATABASE
# ─────────────────────────────────────────────────────────────────────────────

@dataclass
class IntelligenceSNP:
    rsid: str
    gene: str
    category: str
    trait: str
    risk_allele: str        # allele associated with LESS favorable tendency
    favorable_allele: str   # allele associated with MORE favorable tendency
    effect: str             # description when favorable allele is present
    evidence: str           # "Strong", "Moderate", or "Preliminary"


INTELLIGENCE_VARIANTS: Dict[str, IntelligenceSNP] = {
    # ── GENERAL INTELLIGENCE ──
    "rs8191992": IntelligenceSNP(
        rsid="rs8191992", gene="CHRM2", category="General Intelligence",
        trait="Muscarinic acetylcholine receptor (cognitive performance)",
        risk_allele="A", favorable_allele="T",
        effect="T allele linked in some studies to slightly higher scores "
               "on general cognitive ability tasks",
        evidence="Preliminary",
    ),
    "rs362549": IntelligenceSNP(
        rsid="rs362549", gene="SNAP25", category="General Intelligence",
        trait="Synaptic vesicle release (SNAP-25)",
        risk_allele="A", favorable_allele="C",
        effect="C allele associated with modestly higher IQ in pediatric "
               "cohorts; role in synaptic transmission",
        evidence="Preliminary",
    ),
    "rs429358": IntelligenceSNP(
        rsid="rs429358", gene="APOE", category="General Intelligence",
        trait="APOE E4 variant & cognitive aging",
        risk_allele="C", favorable_allele="T",
        effect="T allele (non-E4) linked to better preservation of "
               "cognitive function in later life",
        evidence="Moderate",
    ),
    "rs2251499": IntelligenceSNP(
        rsid="rs2251499", gene="GWAS-IQ locus", category="General Intelligence",
        trait="Intelligence GWAS locus (Savage et al. 2018)",
        risk_allele="T", favorable_allele="C",
        effect="C allele tagged by large GWAS of intelligence (~248k samples); "
               "very small per-variant contribution",
        evidence="Moderate",
    ),
    "rs4851266": IntelligenceSNP(
        rsid="rs4851266", gene="LINC01104", category="General Intelligence",
        trait="Educational attainment / cognition GWAS hit",
        risk_allele="C", favorable_allele="T",
        effect="T allele associated with slightly higher educational "
               "attainment in large GWAS; tiny effect size",
        evidence="Preliminary",
    ),

    # ── MEMORY ──
    "rs17070145": IntelligenceSNP(
        rsid="rs17070145", gene="KIBRA/WWC1", category="Memory",
        trait="Episodic memory performance",
        risk_allele="C", favorable_allele="T",
        effect="T allele associated with better episodic memory recall "
               "across multiple independent cohorts",
        evidence="Moderate",
    ),
    "rs4680": IntelligenceSNP(
        rsid="rs4680", gene="COMT", category="Memory",
        trait="Catechol-O-methyltransferase Val158Met (prefrontal dopamine)",
        risk_allele="G", favorable_allele="A",
        effect="A allele (Met158) associated with lower COMT activity, higher "
               "prefrontal dopamine, and better working memory performance",
        evidence="Moderate",
    ),
    "rs6265": IntelligenceSNP(
        rsid="rs6265", gene="BDNF", category="Memory",
        trait="BDNF Val66Met (activity-dependent secretion)",
        risk_allele="A", favorable_allele="G",
        effect="G allele (Val66) supports normal activity-dependent BDNF "
               "secretion, linked to better episodic memory",
        evidence="Moderate",
    ),
    "rs6439886": IntelligenceSNP(
        rsid="rs6439886", gene="CLSTN2", category="Memory",
        trait="Calsyntenin-2 & episodic memory",
        risk_allele="C", favorable_allele="T",
        effect="T allele linked in candidate-gene studies to better "
               "episodic memory task performance",
        evidence="Preliminary",
    ),
    "rs4343": IntelligenceSNP(
        rsid="rs4343", gene="ACE", category="Memory",
        trait="ACE activity and memory tasks",
        risk_allele="G", favorable_allele="A",
        effect="A allele associated with more favorable memory performance "
               "in some elderly cohorts; mechanism via brain RAS",
        evidence="Preliminary",
    ),

    # ── PROCESSING SPEED ──
    "rs1800497": IntelligenceSNP(
        rsid="rs1800497", gene="DRD2/ANKK1", category="Processing Speed",
        trait="Dopamine D2 receptor density (Taq1A)",
        risk_allele="A", favorable_allele="G",
        effect="G allele (A2) associated with higher striatal D2 receptor "
               "availability and faster processing on some tasks",
        evidence="Moderate",
    ),
    "rs1018381": IntelligenceSNP(
        rsid="rs1018381", gene="DTNBP1", category="Processing Speed",
        trait="Dysbindin & cognitive processing",
        risk_allele="T", favorable_allele="C",
        effect="C allele linked to more typical dysbindin expression and "
               "slightly better processing speed",
        evidence="Preliminary",
    ),
    "rs2760118": IntelligenceSNP(
        rsid="rs2760118", gene="ALDH5A1", category="Processing Speed",
        trait="Succinic semialdehyde dehydrogenase (SSADH)",
        risk_allele="T", favorable_allele="C",
        effect="C allele associated with normal GABA degradation and "
               "slightly better processing speed in small studies",
        evidence="Preliminary",
    ),
    "rs4938021": IntelligenceSNP(
        rsid="rs4938021", gene="ANKK1", category="Processing Speed",
        trait="Dopaminergic signaling & reaction time",
        risk_allele="T", favorable_allele="C",
        effect="C allele linked to faster reaction times on attention "
               "tasks in candidate-gene studies",
        evidence="Preliminary",
    ),

    # ── VERBAL ABILITY ──
    "rs2253478": IntelligenceSNP(
        rsid="rs2253478", gene="FOXP2", category="Verbal Ability",
        trait="Transcription factor for speech & language",
        risk_allele="A", favorable_allele="G",
        effect="G allele tags typical FOXP2 regulation, a gene critical "
               "for speech development and language processing",
        evidence="Moderate",
    ),
    "rs7794745": IntelligenceSNP(
        rsid="rs7794745", gene="CNTNAP2", category="Verbal Ability",
        trait="Language development & neural connectivity",
        risk_allele="T", favorable_allele="A",
        effect="A allele associated with typical language acquisition "
               "trajectories; variant implicated in specific language impairment",
        evidence="Moderate",
    ),
    "rs3743205": IntelligenceSNP(
        rsid="rs3743205", gene="DYX1C1", category="Verbal Ability",
        trait="Dyslexia susceptibility & reading ability",
        risk_allele="A", favorable_allele="G",
        effect="G allele is the more common reference variant, linked to "
               "typical reading and verbal processing development",
        evidence="Preliminary",
    ),
    "rs6803202": IntelligenceSNP(
        rsid="rs6803202", gene="ROBO1", category="Verbal Ability",
        trait="Axon guidance & phonological buffer",
        risk_allele="A", favorable_allele="G",
        effect="G allele linked to typical ROBO1 function supporting "
               "phonological short-term memory",
        evidence="Preliminary",
    ),

    # ── LEARNING & PLASTICITY ──
    "rs3924999": IntelligenceSNP(
        rsid="rs3924999", gene="NRG1", category="Learning & Plasticity",
        trait="Neuregulin-1 & synaptic plasticity",
        risk_allele="A", favorable_allele="G",
        effect="G allele associated with typical NRG1 signaling supporting "
               "myelination and synaptic plasticity",
        evidence="Preliminary",
    ),
    "rs2032582": IntelligenceSNP(
        rsid="rs2032582", gene="NCAM1", category="Learning & Plasticity",
        trait="Neural cell adhesion molecule & learning",
        risk_allele="T", favorable_allele="G",
        effect="G allele linked to typical NCAM1 expression, supporting "
               "synaptic remodeling during learning",
        evidence="Preliminary",
    ),
    "rs1006737": IntelligenceSNP(
        rsid="rs1006737", gene="CACNA1C", category="Learning & Plasticity",
        trait="L-type calcium channel & neural plasticity",
        risk_allele="A", favorable_allele="G",
        effect="G allele associated with typical calcium signaling supporting "
               "long-term potentiation and learning",
        evidence="Moderate",
    ),
    "rs10994336": IntelligenceSNP(
        rsid="rs10994336", gene="ANK3", category="Learning & Plasticity",
        trait="Ankyrin-3 & neuronal excitability",
        risk_allele="T", favorable_allele="C",
        effect="C allele linked to typical ANK3 function supporting axonal "
               "firing and neural plasticity",
        evidence="Preliminary",
    ),
    "rs1344706": IntelligenceSNP(
        rsid="rs1344706", gene="ZNF804A", category="Learning & Plasticity",
        trait="Neural connectivity & cognitive flexibility",
        risk_allele="A", favorable_allele="C",
        effect="C allele associated with typical connectivity patterns "
               "supporting cognitive flexibility",
        evidence="Preliminary",
    ),
    "rs1800532": IntelligenceSNP(
        rsid="rs1800532", gene="TPH1", category="Learning & Plasticity",
        trait="Tryptophan hydroxylase & serotonergic plasticity",
        risk_allele="A", favorable_allele="C",
        effect="C allele linked to typical serotonin synthesis supporting "
               "mood-dependent learning processes",
        evidence="Preliminary",
    ),
    "rs1800955": IntelligenceSNP(
        rsid="rs1800955", gene="DRD4", category="Learning & Plasticity",
        trait="DRD4 promoter & novelty-driven learning",
        risk_allele="T", favorable_allele="C",
        effect="C allele associated with typical D4 receptor expression "
               "involved in exploratory learning and attention",
        evidence="Preliminary",
    ),
}


# ─────────────────────────────────────────────────────────────────────────────
# ANALYSIS RESULT DATACLASSES
# ─────────────────────────────────────────────────────────────────────────────

@dataclass
class IntelligenceVariantResult:
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
class IntelligenceAnalysisResult:
    total_checked: int
    found: int
    not_found: int
    category_scores: Dict[str, float]            # category -> 0-100 score
    category_max_scores: Dict[str, float]        # category -> max possible
    cognitive_score: float                       # aggregate cognitive score 0-100
    results_by_category: Dict[str, List[IntelligenceVariantResult]]
    missing_rsids: List[str]


# ─────────────────────────────────────────────────────────────────────────────
# CORE ANALYSIS
# ─────────────────────────────────────────────────────────────────────────────

def analyze_intelligence(
    variants: Dict[str, Tuple[str, str, str]]
) -> IntelligenceAnalysisResult:
    """
    Analyze genetic variants for cognitive tendencies.

    IMPORTANT: This describes statistical tendencies observed in research
    populations, NOT an individual's intelligence or potential. Environment,
    education, nutrition, sleep, and effort dominate real-world outcomes.

    Args:
        variants: Dictionary mapping rsID (lowercase) to (chromosome, position, genotype)

    Returns:
        IntelligenceAnalysisResult with category scores and per-variant results
    """
    results_by_category: Dict[str, List[IntelligenceVariantResult]] = {
        cat: [] for cat in CATEGORIES
    }
    missing_rsids: List[str] = []

    favorable_counts: Dict[str, int] = {cat: 0 for cat in CATEGORIES}
    max_possible: Dict[str, int] = {cat: 0 for cat in CATEGORIES}

    found = 0
    not_found = 0

    for rsid, snp_info in INTELLIGENCE_VARIANTS.items():
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
                IntelligenceVariantResult(
                    rsid=rsid,
                    gene=snp_info.gene,
                    category=snp_info.category,
                    trait=snp_info.trait,
                    genotype=genotype,
                    favorable_allele_count=fav_count,
                    effect=snp_info.effect if fav_count > 0 else (
                        f"No {snp_info.favorable_allele} allele detected; "
                        f"tendency described does not apply to this genotype"
                    ),
                    interpretation=interpretation,
                    evidence=snp_info.evidence,
                )
            )
        else:
            not_found += 1
            missing_rsids.append(f"{rsid} ({snp_info.gene})")

    # Calculate category scores (0-100)
    category_scores: Dict[str, float] = {}
    category_max: Dict[str, float] = {}
    for cat in CATEGORIES:
        mp = max_possible[cat]
        if mp > 0:
            category_scores[cat] = round((favorable_counts[cat] / mp) * 100, 1)
        else:
            category_scores[cat] = 0.0
        category_max[cat] = float(mp)

    # Aggregate cognitive score = average across categories with data
    non_zero_scores = [s for s in category_scores.values() if s > 0]
    cognitive_score = (
        round(sum(non_zero_scores) / len(non_zero_scores), 1)
        if non_zero_scores
        else 0.0
    )

    return IntelligenceAnalysisResult(
        total_checked=len(INTELLIGENCE_VARIANTS),
        found=found,
        not_found=not_found,
        category_scores=category_scores,
        category_max_scores=category_max,
        cognitive_score=cognitive_score,
        results_by_category=results_by_category,
        missing_rsids=missing_rsids,
    )


def _build_interpretation(
    genotype: str, snp: IntelligenceSNP, fav_count: int
) -> str:
    """Build a human-readable interpretation for a single variant result.

    Uses 'tendency' language to avoid deterministic framing.
    """
    fav = snp.favorable_allele.upper()
    risk = snp.risk_allele.upper()

    if fav_count == 2:
        return (
            f"You carry two copies of the {fav} allele (homozygous). "
            f"Studies describe a mild statistical tendency toward "
            f"{snp.trait.lower()} with this genotype. Real-world impact is "
            f"very small compared to environment, education, and effort."
        )
    elif fav_count == 1:
        return (
            f"You carry one copy of the {fav} allele (heterozygous). "
            f"This represents an intermediate genetic tendency related to "
            f"{snp.trait.lower()}. Environmental factors remain dominant."
        )
    else:
        return (
            f"You carry two copies of the {risk} allele. "
            f"The statistical tendency described for this variant does not "
            f"apply strongly to you; this says nothing about your actual "
            f"cognitive abilities, which depend overwhelmingly on environment, "
            f"education, and effort."
        )


# ─────────────────────────────────────────────────────────────────────────────
# JSON OUTPUT
# ─────────────────────────────────────────────────────────────────────────────

def generate_intelligence_json(result: IntelligenceAnalysisResult) -> dict:
    """
    Generate a JSON-serializable dict for the frontend intelligence report.
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
        "disclaimer": INTELLIGENCE_DISCLAIMER,
        "summary": {
            "totalChecked": result.total_checked,
            "found": result.found,
            "notFound": result.not_found,
            "cognitiveScore": result.cognitive_score,
            "generalScore": result.category_scores.get("General Intelligence", 0.0),
            "memoryScore": result.category_scores.get("Memory", 0.0),
            "processingScore": result.category_scores.get("Processing Speed", 0.0),
            "verbalScore": result.category_scores.get("Verbal Ability", 0.0),
            "learningScore": result.category_scores.get("Learning & Plasticity", 0.0),
        },
        "categories": categories_list,
    }


# ─────────────────────────────────────────────────────────────────────────────
# MARKDOWN REPORT (optional, for compatibility)
# ─────────────────────────────────────────────────────────────────────────────

def generate_intelligence_report(
    result: IntelligenceAnalysisResult,
    subject_name: str = "Subject",
) -> str:
    """Generate a Markdown intelligence report."""
    from datetime import datetime

    report = f"""# Cognitive Tendencies Report

**Subject:** {subject_name}
**Generated:** {datetime.now().strftime("%B %d, %Y")}

---

## IMPORTANT DISCLAIMER

> {INTELLIGENCE_DISCLAIMER}

---

## Overview

This report explores genetic variants that have been associated, in research
populations, with small statistical tendencies in cognitive tasks. It is
emphatically NOT a measurement of your intelligence. Environment, education,
nutrition, sleep, social context, and effort have vastly larger effects on
cognition than any combination of common SNPs.

| Metric | Value |
|--------|-------|
| Total SNPs Checked | {result.total_checked} |
| SNPs Found in Your Data | {result.found} |
| SNPs Not Available | {result.not_found} |
| Aggregate Cognitive Score (curiosity only) | {result.cognitive_score}/100 |

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

1. **Genes do not determine intelligence** - Twin and adoption studies show
   heritability for cognitive traits, but heritability is a population statistic,
   not a prediction for individuals, and most heritability is not captured by
   common SNPs.

2. **Effect sizes are tiny** - Even the largest GWAS of intelligence find that
   each individual SNP explains well under 0.1% of variance in cognitive scores.

3. **Environment dominates** - Schooling, nutrition, sleep, stress, social
   support, and engagement with challenging tasks have vastly larger effects on
   cognitive development and performance than any genetic variant.

4. **Neuroplasticity is lifelong** - The brain adapts throughout life. Learning,
   practice, physical exercise, and social engagement continue to shape
   cognitive ability at every age.

5. **Do not use this report to self-label** - These variants describe
   tendencies in statistical aggregates, not traits of individuals. Treat this
   report as educational curiosity only.

---
*Generated by GeneHealth Analysis Platform*
"""

    return report
