"""
Sports Analyzer
Analyzes genetic variants related to specific SPORT SUITABILITY across five
categories: Sprint & Explosive, Endurance, Strength, Mixed/Team, and Precision.

Unlike fitness_analyzer.py which produces general athletic dimension scores,
this module translates genotype into SPECIFIC SPORT RECOMMENDATIONS, helping
the user understand which disciplines their genetic profile best matches.

References:
- Yang et al. (2003) AJHG - ACTN3 R577X
- Ahmetov & Fedotovskaya (2015) Biol Sport - Genes and athletic performance
- Schuelke et al. (2004) NEJM - MSTN (myostatin) and muscle mass
- Williams et al. (2000) Nature - ACE I/D and endurance
- Eynon et al. (2011) - PPARGC1A and endurance
- Voisin et al. (2014) - Review of endurance polymorphisms
- Hagberg et al. (2011) - AGT and strength training response
- Kostek et al. (2005) - IGF1 and muscle strength response
- Montgomery et al. (1999) - ACE and altitude performance
- Egan et al. (2015) - EPAS1 and high-altitude adaptation
- Bouchard et al. (2011) - HERITAGE study on training response
- Chen et al. (2020) BDNF and motor learning
- Stein et al. (2006) COMT and cognitive/motor performance
- Tsai et al. (2010) DRD2 and reaction time
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

    if count == 0 and comp != allele:
        count = genotype.count(comp)

    return min(count, 2)


# ─────────────────────────────────────────────────────────────────────────────
# CATEGORY DEFINITIONS
# ─────────────────────────────────────────────────────────────────────────────

CATEGORIES = [
    "Sprint & Explosive Sports",
    "Endurance Sports (Long Distance)",
    "Strength Sports",
    "Mixed/Team Sports",
    "Precision Sports",
]

CATEGORY_EMOJI = {
    "Sprint & Explosive Sports": "\U0001f4a5",       # 💥
    "Endurance Sports (Long Distance)": "\U0001f3c3", # 🏃
    "Strength Sports": "\U0001f4aa",                  # 💪
    "Mixed/Team Sports": "\u26bd",                    # ⚽
    "Precision Sports": "\U0001f3af",                 # 🎯
}

CATEGORY_RECOMMENDED_SPORTS = {
    "Sprint & Explosive Sports": [
        "100m / 200m sprint",
        "Olympic weightlifting",
        "High jump & long jump",
        "Volleyball (spiker)",
        "Short-track speed skating",
    ],
    "Endurance Sports (Long Distance)": [
        "Marathon & ultra-running",
        "Road cycling & time trial",
        "Triathlon (Ironman)",
        "Cross-country skiing",
        "Open-water swimming",
    ],
    "Strength Sports": [
        "Powerlifting",
        "Strongman competitions",
        "Rugby (forward)",
        "Shot put & discus throw",
        "Bodybuilding",
    ],
    "Mixed/Team Sports": [
        "Soccer / football (midfielder)",
        "Basketball",
        "Tennis",
        "CrossFit",
        "Mixed martial arts (MMA)",
    ],
    "Precision Sports": [
        "Archery",
        "Shooting sports",
        "Golf",
        "Table tennis",
        "Fencing",
    ],
}


# ─────────────────────────────────────────────────────────────────────────────
# SPORTS VARIANT DATABASE
# ─────────────────────────────────────────────────────────────────────────────

@dataclass
class SportsSNP:
    rsid: str
    gene: str
    category: str
    trait: str
    risk_allele: str        # less favorable for this SPORT category
    favorable_allele: str   # more favorable for this SPORT category
    effect: str
    evidence: str           # "Strong", "Moderate", or "Preliminary"


SPORTS_VARIANTS: Dict[str, SportsSNP] = {
    # ── SPRINT & EXPLOSIVE SPORTS ──
    "rs1815739": SportsSNP(
        rsid="rs1815739", gene="ACTN3", category="Sprint & Explosive Sports",
        trait="Alpha-actinin-3 in fast-twitch fibers (R577X)",
        risk_allele="T", favorable_allele="C",
        effect="CC (RR) genotype is overrepresented in elite sprinters and "
               "power athletes; enables explosive force production",
        evidence="Strong",
    ),
    "rs699": SportsSNP(
        rsid="rs699", gene="AGT", category="Sprint & Explosive Sports",
        trait="Angiotensinogen M235T (explosive power)",
        risk_allele="G", favorable_allele="A",
        effect="A allele (235T) is associated with greater strength and "
               "power gains in response to explosive training",
        evidence="Moderate",
    ),
    "rs17602729": SportsSNP(
        rsid="rs17602729", gene="AMPD1", category="Sprint & Explosive Sports",
        trait="AMP deaminase activity for anaerobic bursts",
        risk_allele="A", favorable_allele="G",
        effect="G allele preserves AMPD1 activity, essential for rapid "
               "ATP regeneration during short, maximal efforts",
        evidence="Moderate",
    ),
    "rs1805086": SportsSNP(
        rsid="rs1805086", gene="MSTN", category="Sprint & Explosive Sports",
        trait="Myostatin K153R (muscle mass ceiling)",
        risk_allele="G", favorable_allele="A",
        effect="A (K allele) associated with typical myostatin function; "
               "variation linked with muscle hypertrophy potential for power",
        evidence="Preliminary",
    ),

    # ── ENDURANCE SPORTS (LONG DISTANCE) ──
    "rs4343": SportsSNP(
        rsid="rs4343", gene="ACE", category="Endurance Sports (Long Distance)",
        trait="ACE I/D tag SNP (cardiovascular efficiency)",
        risk_allele="A", favorable_allele="G",
        effect="G allele tagging the I insertion is associated with "
               "superior endurance performance and aerobic efficiency",
        evidence="Strong",
    ),
    "rs8192678": SportsSNP(
        rsid="rs8192678", gene="PPARGC1A",
        category="Endurance Sports (Long Distance)",
        trait="PGC-1alpha Gly482Ser (mitochondrial biogenesis)",
        risk_allele="A", favorable_allele="G",
        effect="Gly482 (G) drives mitochondrial biogenesis, higher VO2max "
               "and superior long-distance performance",
        evidence="Strong",
    ),
    "rs11549465": SportsSNP(
        rsid="rs11549465", gene="HIF1A",
        category="Endurance Sports (Long Distance)",
        trait="Hypoxia response Pro582Ser",
        risk_allele="C", favorable_allele="T",
        effect="T allele (Ser582) improves oxygen sensing and endurance "
               "capacity, particularly at altitude",
        evidence="Moderate",
    ),
    "rs1867785": SportsSNP(
        rsid="rs1867785", gene="EPAS1",
        category="Endurance Sports (Long Distance)",
        trait="HIF-2alpha altitude/endurance adaptation",
        risk_allele="A", favorable_allele="G",
        effect="G allele linked to favorable EPAS1 signaling for oxygen "
               "delivery in endurance athletes",
        evidence="Preliminary",
    ),
    "rs1042713": SportsSNP(
        rsid="rs1042713", gene="ADRB2",
        category="Endurance Sports (Long Distance)",
        trait="Beta-2 adrenergic receptor Arg16Gly",
        risk_allele="A", favorable_allele="G",
        effect="Gly16 (G) promotes bronchodilation and improved aerobic "
               "capacity crucial for sustained endurance events",
        evidence="Moderate",
    ),

    # ── STRENGTH SPORTS ──
    "rs35767": SportsSNP(
        rsid="rs35767", gene="IGF1", category="Strength Sports",
        trait="IGF-1 promoter variant (growth/strength response)",
        risk_allele="A", favorable_allele="G",
        effect="G allele associated with higher IGF-1 levels and greater "
               "muscle strength gains from resistance training",
        evidence="Moderate",
    ),
    "rs2296135": SportsSNP(
        rsid="rs2296135", gene="IL15RA", category="Strength Sports",
        trait="IL-15 receptor alpha (muscle hypertrophy)",
        risk_allele="C", favorable_allele="A",
        effect="A allele linked to greater muscle hypertrophy response "
               "and improved isometric strength in resistance trainers",
        evidence="Moderate",
    ),
    "rs10783485": SportsSNP(
        rsid="rs10783485", gene="ACVR1B", category="Strength Sports",
        trait="Activin receptor 1B (muscle strength)",
        risk_allele="G", favorable_allele="A",
        effect="A allele associated with higher baseline muscular "
               "strength and greater gains after strength training",
        evidence="Preliminary",
    ),
    "rs1799983": SportsSNP(
        rsid="rs1799983", gene="NOS3", category="Strength Sports",
        trait="Endothelial NOS Glu298Asp (muscle blood flow)",
        risk_allele="T", favorable_allele="G",
        effect="G allele (Glu298) supports nitric oxide production, blood "
               "flow to working muscles and strength expression",
        evidence="Moderate",
    ),

    # ── MIXED / TEAM SPORTS ──
    "rs4253778": SportsSNP(
        rsid="rs4253778", gene="PPARA", category="Mixed/Team Sports",
        trait="PPAR-alpha fat/carb substrate utilization",
        risk_allele="C", favorable_allele="G",
        effect="G allele balances fat and carbohydrate oxidation, ideal "
               "for repeated high-intensity efforts in team sports",
        evidence="Moderate",
    ),
    "rs12722": SportsSNP(
        rsid="rs12722", gene="COL5A1", category="Mixed/Team Sports",
        trait="Collagen V (joint flexibility for cutting sports)",
        risk_allele="T", favorable_allele="C",
        effect="C allele associated with better tendon flexibility and "
               "reduced ACL/Achilles injury risk in cutting & pivoting sports",
        evidence="Strong",
    ),
    "rs1042714": SportsSNP(
        rsid="rs1042714", gene="ADRB2", category="Mixed/Team Sports",
        trait="Beta-2 adrenergic Gln27Glu (cardiovascular response)",
        risk_allele="C", favorable_allele="G",
        effect="G allele (Glu27) linked to improved cardiovascular "
               "response during intermittent high-intensity efforts",
        evidence="Preliminary",
    ),
    "rs2070744": SportsSNP(
        rsid="rs2070744", gene="NOS3", category="Mixed/Team Sports",
        trait="NOS3 T-786C promoter (endothelial function)",
        risk_allele="C", favorable_allele="T",
        effect="T allele promotes higher eNOS expression, supporting "
               "blood flow during mixed aerobic-anaerobic efforts",
        evidence="Moderate",
    ),

    # ── PRECISION SPORTS ──
    "rs6265": SportsSNP(
        rsid="rs6265", gene="BDNF", category="Precision Sports",
        trait="BDNF Val66Met (motor learning & skill acquisition)",
        risk_allele="A", favorable_allele="G",
        effect="G allele (Val66) associated with enhanced motor learning, "
               "procedural memory, and fine skill acquisition",
        evidence="Strong",
    ),
    "rs4680": SportsSNP(
        rsid="rs4680", gene="COMT", category="Precision Sports",
        trait="COMT Val158Met (dopamine & focus under pressure)",
        risk_allele="A", favorable_allele="G",
        effect="G allele (Val158) linked to stable dopamine under stress, "
               "supporting focus and precision in pressured situations",
        evidence="Moderate",
    ),
    "rs1076560": SportsSNP(
        rsid="rs1076560", gene="DRD2", category="Precision Sports",
        trait="Dopamine receptor D2 (reaction time & accuracy)",
        risk_allele="A", favorable_allele="C",
        effect="C allele associated with faster reaction time and "
               "improved motor precision in target-based sports",
        evidence="Moderate",
    ),
    "rs25531": SportsSNP(
        rsid="rs25531", gene="SLC6A4", category="Precision Sports",
        trait="Serotonin transporter (stress regulation)",
        risk_allele="G", favorable_allele="A",
        effect="A allele linked to better emotional regulation under "
               "performance pressure, aiding precision and focus",
        evidence="Preliminary",
    ),
    "rs8191992": SportsSNP(
        rsid="rs8191992", gene="CHRM2", category="Precision Sports",
        trait="Muscarinic ACh receptor M2 (attention & concentration)",
        risk_allele="A", favorable_allele="T",
        effect="T allele associated with improved sustained attention "
               "and concentration, beneficial for precision disciplines",
        evidence="Preliminary",
    ),
}


# ─────────────────────────────────────────────────────────────────────────────
# ANALYSIS RESULT DATACLASSES
# ─────────────────────────────────────────────────────────────────────────────

@dataclass
class SportsVariantResult:
    rsid: str
    gene: str
    category: str
    trait: str
    genotype: str
    favorable_allele_count: int
    effect: str
    interpretation: str
    evidence: str


@dataclass
class SportsAnalysisResult:
    total_checked: int
    found: int
    not_found: int
    category_scores: Dict[str, float]
    category_max_scores: Dict[str, float]
    results_by_category: Dict[str, List[SportsVariantResult]]
    top_category: str
    recommendation_text: str
    missing_rsids: List[str]


# ─────────────────────────────────────────────────────────────────────────────
# CORE ANALYSIS
# ─────────────────────────────────────────────────────────────────────────────

def analyze_sports(
    variants: Dict[str, Tuple[str, str, str]]
) -> SportsAnalysisResult:
    """
    Analyze genetic variants for sport suitability.

    Args:
        variants: Dictionary mapping rsID (lowercase) to (chromosome, position, genotype)

    Returns:
        SportsAnalysisResult with per-category scores and sport recommendations
    """
    results_by_category: Dict[str, List[SportsVariantResult]] = {
        cat: [] for cat in CATEGORIES
    }
    missing_rsids: List[str] = []

    favorable_counts: Dict[str, int] = {cat: 0 for cat in CATEGORIES}
    max_possible: Dict[str, int] = {cat: 0 for cat in CATEGORIES}

    found = 0
    not_found = 0

    for rsid, snp_info in SPORTS_VARIANTS.items():
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
                SportsVariantResult(
                    rsid=rsid,
                    gene=snp_info.gene,
                    category=snp_info.category,
                    trait=snp_info.trait,
                    genotype=genotype,
                    favorable_allele_count=fav_count,
                    effect=snp_info.effect if fav_count > 0 else (
                        f"No {snp_info.favorable_allele} allele detected; "
                        f"typical or less favorable variant for this sport category"
                    ),
                    interpretation=interpretation,
                    evidence=snp_info.evidence,
                )
            )
        else:
            not_found += 1
            missing_rsids.append(f"{rsid} ({snp_info.gene})")

    # Calculate category scores
    category_scores: Dict[str, float] = {}
    category_max: Dict[str, float] = {}
    for cat in CATEGORIES:
        mp = max_possible[cat]
        if mp > 0:
            category_scores[cat] = round((favorable_counts[cat] / mp) * 100, 1)
        else:
            category_scores[cat] = 0.0
        category_max[cat] = float(mp)

    # Determine top category
    top_category = max(category_scores, key=lambda c: category_scores[c]) \
        if category_scores else CATEGORIES[0]

    recommendation_text = _build_recommendation_text(top_category, category_scores)

    return SportsAnalysisResult(
        total_checked=len(SPORTS_VARIANTS),
        found=found,
        not_found=not_found,
        category_scores=category_scores,
        category_max_scores=category_max,
        results_by_category=results_by_category,
        top_category=top_category,
        recommendation_text=recommendation_text,
        missing_rsids=missing_rsids,
    )


def _build_interpretation(genotype: str, snp: SportsSNP, fav_count: int) -> str:
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
            f"Intermediate predisposition for {snp.trait.lower()}."
        )
    else:
        return (
            f"You carry two copies of the {risk} allele. "
            f"Less favorable predisposition for {snp.trait.lower()}, though "
            f"training remains the dominant factor."
        )


def _build_recommendation_text(top_category: str, scores: Dict[str, float]) -> str:
    """Build a narrative recommendation based on the top sport category."""
    top_sports = CATEGORY_RECOMMENDED_SPORTS.get(top_category, [])
    sports_list = ", ".join(top_sports[:3]) if top_sports else "various sports"
    top_score = scores.get(top_category, 0.0)

    return (
        f"Your genetic profile shows the strongest signal for "
        f"{top_category} (score {top_score}/100). Sports like {sports_list} "
        f"may suit your genetic strengths particularly well. Remember that "
        f"training, coaching, and personal preference remain the most "
        f"important determinants of athletic success."
    )


# ─────────────────────────────────────────────────────────────────────────────
# JSON OUTPUT
# ─────────────────────────────────────────────────────────────────────────────

def generate_sports_json(result: SportsAnalysisResult) -> dict:
    """
    Generate a JSON-serializable dict for the frontend sports report.

    Args:
        result: SportsAnalysisResult from analyze_sports()

    Returns:
        Dict matching the sports report JSON schema
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
            "recommendedSports": CATEGORY_RECOMMENDED_SPORTS.get(cat, []),
        })

    return {
        "summary": {
            "totalChecked": result.total_checked,
            "found": result.found,
            "notFound": result.not_found,
            "topCategory": result.top_category,
            "topCategoryEmoji": CATEGORY_EMOJI.get(result.top_category, ""),
            "topCategoryScore": result.category_scores.get(result.top_category, 0.0),
            "recommendation": result.recommendation_text,
        },
        "categories": categories_list,
    }


# ─────────────────────────────────────────────────────────────────────────────
# MARKDOWN REPORT (optional, for compatibility)
# ─────────────────────────────────────────────────────────────────────────────

def generate_sports_report(result: SportsAnalysisResult, subject_name: str = "Subject") -> str:
    """Generate a Markdown sports suitability report."""
    from datetime import datetime

    report = f"""# Sports Suitability Report

**Subject:** {subject_name}
**Generated:** {datetime.now().strftime("%B %d, %Y")}

---

## Overview

This report translates your genetic profile into SPORT RECOMMENDATIONS.
Unlike a general fitness profile, it suggests which disciplines best match
your genetic strengths across five categories.

| Metric | Value |
|--------|-------|
| Total Sports SNPs Checked | {result.total_checked} |
| SNPs Found in Your Data | {result.found} |
| SNPs Not Available | {result.not_found} |
| Top Category | {result.top_category} |

### {result.recommendation_text}

---

"""

    for cat in CATEGORIES:
        findings = result.results_by_category.get(cat, [])
        emoji = CATEGORY_EMOJI.get(cat, "")
        score = result.category_scores.get(cat, 0.0)
        sports = CATEGORY_RECOMMENDED_SPORTS.get(cat, [])

        report += f"## {emoji} {cat}  —  {score}/100\n\n"
        report += "**Recommended sports:**\n"
        for s in sports:
            report += f"- {s}\n"
        report += "\n"

        if findings:
            for f in findings:
                fav_icon = {0: "\u26aa", 1: "\U0001f535", 2: "\U0001f7e2"}.get(
                    f.favorable_allele_count, "\u26aa"
                )
                report += f"{fav_icon} **{f.gene}** ({f.rsid}) — {f.genotype}\n"
                report += f"- {f.trait}\n"
                report += f"- {f.interpretation}\n\n"

        report += "---\n\n"

    report += """## Important Notes

1. **Genetics is one factor** — Sport suitability depends heavily on training,
   coaching, personal preference, and opportunity. Genetic profile is a hint,
   not a verdict.

2. **Scores are relative** — Category scores reflect the proportion of
   favorable alleles detected in a limited SNP panel.

3. **Enjoyment matters most** — The best sport for you is often the one you
   love and will stick with consistently.

---
*Generated by GeneHealth Analysis Platform*
"""

    return report
