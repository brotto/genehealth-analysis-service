"""
Organ & Blood Wellness Analyzer
Analyzes genetic variants related to organ health and blood markers across five
categories: Kidney Health, Liver Function, Cardiovascular, Blood & Iron, and Cholesterol.

References:
- Kottgen et al. (2010) Nat Genet - UMOD, SHROOM3 and eGFR
- Genovese et al. (2010) Science - APOL1 and kidney disease
- Romeo et al. (2008) Nat Genet - PNPLA3 and NAFLD
- Kozlitina et al. (2014) Nat Genet - TM6SF2 and liver disease
- Feder et al. (1996) Nat Genet - HFE and hereditary hemochromatosis
- Bennet et al. (2007) JAMA - APOE and cardiovascular risk
- Frosst et al. (1995) Nat Genet - MTHFR C677T and homocysteine
- Ridker et al. (1995) NEJM - Factor V Leiden (F5)
- Poort et al. (1996) Blood - Prothrombin G20210A (F2)
- Benjafield et al. (2003) Am J Hypertens - ACE
- Hingorani et al. (1999) Hypertension - NOS3 Glu298Asp
- Benyamin et al. (2014) Nat Commun - TMPRSS6 and iron status
- Cohen et al. (2006) NEJM - PCSK9 and LDL cholesterol
- Thompson et al. (2008) NEJM - CETP and HDL
- Musunuru et al. (2010) Nature - SORT1/1p13 and LDL
- Jeunemaitre et al. (1992) Cell - AGT M235T
- Bonnardeaux et al. (1994) Hypertension - AGTR1 A1166C
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
    "Kidney Health",
    "Liver Function",
    "Cardiovascular",
    "Blood & Iron",
    "Cholesterol",
]

CATEGORY_EMOJI = {
    "Kidney Health": "\U0001fad8",   # 🫘
    "Liver Function": "\U0001f33f",  # 🌿
    "Cardiovascular": "\u2764\ufe0f", # ❤️
    "Blood & Iron": "\U0001fa78",    # 🩸
    "Cholesterol": "\U0001f9c8",     # 🧈
}


# ─────────────────────────────────────────────────────────────────────────────
# ORGAN & BLOOD WELLNESS VARIANT DATABASE
# ─────────────────────────────────────────────────────────────────────────────

@dataclass
class OrganBloodSNP:
    rsid: str
    gene: str
    category: str
    trait: str
    risk_allele: str        # allele associated with LESS favorable outcome
    favorable_allele: str   # allele associated with MORE favorable outcome
    effect: str             # description when favorable allele is present
    evidence: str           # "Strong", "Moderate", or "Preliminary"


ORGAN_BLOOD_VARIANTS: Dict[str, OrganBloodSNP] = {
    # ── KIDNEY HEALTH ──
    "rs12917707": OrganBloodSNP(
        rsid="rs12917707", gene="UMOD", category="Kidney Health",
        trait="Uromodulin expression & kidney filtration (eGFR)",
        risk_allele="G", favorable_allele="T",
        effect="T allele associated with lower UMOD expression, higher eGFR, "
               "and reduced risk of chronic kidney disease",
        evidence="Strong",
    ),
    "rs17319721": OrganBloodSNP(
        rsid="rs17319721", gene="SHROOM3", category="Kidney Health",
        trait="Glomerular filtration rate (eGFR)",
        risk_allele="A", favorable_allele="G",
        effect="G allele associated with higher eGFR and reduced risk of "
               "chronic kidney disease in GWAS meta-analyses",
        evidence="Strong",
    ),
    "rs73885319": OrganBloodSNP(
        rsid="rs73885319", gene="APOL1", category="Kidney Health",
        trait="APOL1 G1 risk variant (focal segmental glomerulosclerosis)",
        risk_allele="G", favorable_allele="A",
        effect="A allele (reference) protects from APOL1-associated nephropathy "
               "seen primarily in individuals of West African ancestry",
        evidence="Strong",
    ),
    "rs4293393": OrganBloodSNP(
        rsid="rs4293393", gene="UMOD", category="Kidney Health",
        trait="Uromodulin promoter & kidney function",
        risk_allele="A", favorable_allele="G",
        effect="G allele linked to better kidney filtration and reduced "
               "hypertension-related kidney disease risk",
        evidence="Moderate",
    ),

    # ── LIVER FUNCTION ──
    "rs738409": OrganBloodSNP(
        rsid="rs738409", gene="PNPLA3", category="Liver Function",
        trait="Hepatic triglyceride content & NAFLD risk",
        risk_allele="G", favorable_allele="C",
        effect="C allele (Ile148) associated with lower liver fat content "
               "and reduced risk of non-alcoholic fatty liver disease",
        evidence="Strong",
    ),
    "rs58542926": OrganBloodSNP(
        rsid="rs58542926", gene="TM6SF2", category="Liver Function",
        trait="Hepatic lipid export (E167K)",
        risk_allele="T", favorable_allele="C",
        effect="C allele (Glu167) supports normal VLDL secretion and "
               "reduced risk of steatohepatitis and liver fibrosis",
        evidence="Strong",
    ),
    "rs1800562": OrganBloodSNP(
        rsid="rs1800562", gene="HFE", category="Liver Function",
        trait="HFE C282Y (hereditary hemochromatosis)",
        risk_allele="A", favorable_allele="G",
        effect="G allele (Cys282) supports normal iron regulation and "
               "protects from iron overload and hepatic injury",
        evidence="Strong",
    ),
    "rs2642438": OrganBloodSNP(
        rsid="rs2642438", gene="MARC1", category="Liver Function",
        trait="Mitochondrial amidoxime-reducing component 1",
        risk_allele="G", favorable_allele="A",
        effect="A allele (Thr165) associated with lower risk of cirrhosis "
               "and all-cause liver mortality",
        evidence="Moderate",
    ),

    # ── CARDIOVASCULAR ──
    "rs429358": OrganBloodSNP(
        rsid="rs429358", gene="APOE", category="Cardiovascular",
        trait="APOE isoform (E4 determinant)",
        risk_allele="C", favorable_allele="T",
        effect="T allele avoids the E4 isoform, associated with lower LDL "
               "cholesterol and reduced cardiovascular and Alzheimer risk",
        evidence="Strong",
    ),
    "rs7412": OrganBloodSNP(
        rsid="rs7412", gene="APOE", category="Cardiovascular",
        trait="APOE isoform (E2 determinant)",
        risk_allele="C", favorable_allele="T",
        effect="T allele defines APOE2, associated with lower LDL cholesterol "
               "and reduced coronary artery disease risk in most contexts",
        evidence="Strong",
    ),
    "rs1801133": OrganBloodSNP(
        rsid="rs1801133", gene="MTHFR", category="Cardiovascular",
        trait="C677T variant & homocysteine metabolism",
        risk_allele="T", favorable_allele="C",
        effect="C allele (Ala222) supports normal MTHFR enzymatic activity "
               "and efficient homocysteine remethylation",
        evidence="Strong",
    ),
    "rs1801131": OrganBloodSNP(
        rsid="rs1801131", gene="MTHFR", category="Cardiovascular",
        trait="A1298C variant & folate metabolism",
        risk_allele="G", favorable_allele="T",
        effect="T allele (Glu429) maintains normal MTHFR activity supporting "
               "folate and homocysteine balance",
        evidence="Moderate",
    ),
    "rs6025": OrganBloodSNP(
        rsid="rs6025", gene="F5", category="Cardiovascular",
        trait="Factor V Leiden (venous thrombosis risk)",
        risk_allele="T", favorable_allele="C",
        effect="C allele (Arg506) retains normal Factor V inactivation by "
               "activated protein C, lowering venous thromboembolism risk",
        evidence="Strong",
    ),
    "rs1799963": OrganBloodSNP(
        rsid="rs1799963", gene="F2", category="Cardiovascular",
        trait="Prothrombin G20210A (thrombosis risk)",
        risk_allele="A", favorable_allele="G",
        effect="G allele maintains normal prothrombin expression and "
               "reduced risk of venous thrombosis compared to A carriers",
        evidence="Strong",
    ),
    "rs4343": OrganBloodSNP(
        rsid="rs4343", gene="ACE", category="Cardiovascular",
        trait="Angiotensin-converting enzyme activity (I/D proxy)",
        risk_allele="G", favorable_allele="A",
        effect="A allele (proxy for I allele) associated with lower ACE "
               "activity and reduced hypertension and LVH risk",
        evidence="Moderate",
    ),
    "rs1799983": OrganBloodSNP(
        rsid="rs1799983", gene="NOS3", category="Cardiovascular",
        trait="Endothelial NO synthase (Glu298Asp)",
        risk_allele="T", favorable_allele="G",
        effect="G allele (Glu298) supports normal nitric oxide production, "
               "vascular tone regulation, and endothelial function",
        evidence="Moderate",
    ),

    # ── BLOOD & IRON ──
    "rs855791": OrganBloodSNP(
        rsid="rs855791", gene="TMPRSS6", category="Blood & Iron",
        trait="Iron homeostasis (matriptase-2, V736A)",
        risk_allele="A", favorable_allele="G",
        effect="G allele (Val736) associated with higher hemoglobin and "
               "serum iron levels, reducing risk of iron-deficiency anemia",
        evidence="Strong",
    ),
    "rs1799945": OrganBloodSNP(
        rsid="rs1799945", gene="HFE", category="Blood & Iron",
        trait="HFE H63D variant & iron absorption",
        risk_allele="G", favorable_allele="C",
        effect="C allele (His63) preserves normal HFE function and balanced "
               "intestinal iron absorption",
        evidence="Strong",
    ),
    "rs7385804": OrganBloodSNP(
        rsid="rs7385804", gene="TFR2", category="Blood & Iron",
        trait="Transferrin receptor 2 & erythrocyte indices",
        risk_allele="C", favorable_allele="A",
        effect="A allele associated with favorable erythrocyte parameters "
               "and efficient iron delivery to developing red cells",
        evidence="Moderate",
    ),

    # ── CHOLESTEROL ──
    "rs11591147": OrganBloodSNP(
        rsid="rs11591147", gene="PCSK9", category="Cholesterol",
        trait="PCSK9 R46L loss-of-function",
        risk_allele="G", favorable_allele="T",
        effect="T allele (Leu46) reduces PCSK9 activity, lowering LDL "
               "cholesterol ~15% and markedly reducing CAD risk",
        evidence="Strong",
    ),
    "rs670": OrganBloodSNP(
        rsid="rs670", gene="APOA1", category="Cholesterol",
        trait="APOA1 promoter & HDL levels",
        risk_allele="G", favorable_allele="A",
        effect="A allele (-75A) associated with higher HDL cholesterol "
               "and favorable lipid profile",
        evidence="Moderate",
    ),
    "rs708272": OrganBloodSNP(
        rsid="rs708272", gene="CETP", category="Cholesterol",
        trait="Cholesteryl ester transfer protein (TaqIB)",
        risk_allele="G", favorable_allele="A",
        effect="A allele (B2) associated with lower CETP activity, higher "
               "HDL cholesterol and reduced cardiovascular risk",
        evidence="Strong",
    ),
    "rs1800588": OrganBloodSNP(
        rsid="rs1800588", gene="LIPC", category="Cholesterol",
        trait="Hepatic lipase promoter (-514 C/T)",
        risk_allele="C", favorable_allele="T",
        effect="T allele associated with lower hepatic lipase activity "
               "and higher HDL cholesterol",
        evidence="Moderate",
    ),
    "rs12740374": OrganBloodSNP(
        rsid="rs12740374", gene="SORT1", category="Cholesterol",
        trait="1p13 locus & LDL cholesterol (SORT1/sortilin)",
        risk_allele="G", favorable_allele="T",
        effect="T allele increases SORT1 expression, reduces LDL "
               "cholesterol and lowers myocardial infarction risk",
        evidence="Strong",
    ),

    # ── ADDITIONAL CARDIOVASCULAR / BLOOD PRESSURE ──
    "rs699": OrganBloodSNP(
        rsid="rs699", gene="AGT", category="Cardiovascular",
        trait="Angiotensinogen M235T & blood pressure",
        risk_allele="A", favorable_allele="G",
        effect="G allele (Met235) associated with lower plasma angiotensinogen "
               "levels and reduced risk of essential hypertension",
        evidence="Moderate",
    ),
    "rs5186": OrganBloodSNP(
        rsid="rs5186", gene="AGTR1", category="Cardiovascular",
        trait="AT1 receptor A1166C & vascular reactivity",
        risk_allele="C", favorable_allele="A",
        effect="A allele associated with lower AT1 receptor activity, "
               "reduced vasoconstriction and lower hypertension risk",
        evidence="Moderate",
    ),
    "rs1799998": OrganBloodSNP(
        rsid="rs1799998", gene="CYP11B2", category="Cardiovascular",
        trait="Aldosterone synthase (-344 T/C)",
        risk_allele="C", favorable_allele="T",
        effect="T allele associated with balanced aldosterone production "
               "and favorable blood pressure regulation",
        evidence="Preliminary",
    ),
}


# ─────────────────────────────────────────────────────────────────────────────
# ANALYSIS RESULT DATACLASSES
# ─────────────────────────────────────────────────────────────────────────────

@dataclass
class OrganBloodVariantResult:
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
class OrganBloodAnalysisResult:
    total_checked: int
    found: int
    not_found: int
    category_scores: Dict[str, float]            # category -> 0-100 score
    category_max_scores: Dict[str, float]        # category -> max possible
    overall_score: float                          # aggregate wellness score 0-100
    results_by_category: Dict[str, List[OrganBloodVariantResult]]
    missing_rsids: List[str]


# ─────────────────────────────────────────────────────────────────────────────
# CORE ANALYSIS
# ─────────────────────────────────────────────────────────────────────────────

def analyze_organ_blood_wellness(
    variants: Dict[str, Tuple[str, str, str]]
) -> OrganBloodAnalysisResult:
    """
    Analyze genetic variants for organ and blood wellness.

    Args:
        variants: Dictionary mapping rsID (lowercase) to (chromosome, position, genotype)

    Returns:
        OrganBloodAnalysisResult with category scores and per-variant results
    """
    results_by_category: Dict[str, List[OrganBloodVariantResult]] = {
        cat: [] for cat in CATEGORIES
    }
    missing_rsids: List[str] = []

    favorable_counts: Dict[str, int] = {cat: 0 for cat in CATEGORIES}
    max_possible: Dict[str, int] = {cat: 0 for cat in CATEGORIES}

    found = 0
    not_found = 0

    for rsid, snp_info in ORGAN_BLOOD_VARIANTS.items():
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
                OrganBloodVariantResult(
                    rsid=rsid,
                    gene=snp_info.gene,
                    category=snp_info.category,
                    trait=snp_info.trait,
                    genotype=genotype,
                    favorable_allele_count=fav_count,
                    effect=snp_info.effect if fav_count > 0 else (
                        f"No {snp_info.favorable_allele} allele detected; "
                        f"less favorable variant for this trait"
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

    # Overall wellness score = average of category scores (weighted equally)
    non_zero_scores = [s for s in category_scores.values() if s > 0]
    overall_score = (
        round(sum(non_zero_scores) / len(non_zero_scores), 1)
        if non_zero_scores
        else 0.0
    )

    return OrganBloodAnalysisResult(
        total_checked=len(ORGAN_BLOOD_VARIANTS),
        found=found,
        not_found=not_found,
        category_scores=category_scores,
        category_max_scores=category_max,
        overall_score=overall_score,
        results_by_category=results_by_category,
        missing_rsids=missing_rsids,
    )


def _build_interpretation(genotype: str, snp: OrganBloodSNP, fav_count: int) -> str:
    """Build a human-readable interpretation for a single variant result."""
    fav = snp.favorable_allele.upper()
    risk = snp.risk_allele.upper()

    if fav_count == 2:
        return (
            f"You carry two copies of the protective {fav} allele (homozygous). "
            f"This is the most favorable genotype for {snp.trait.lower()}."
        )
    elif fav_count == 1:
        return (
            f"You carry one copy of the protective {fav} allele (heterozygous). "
            f"You have an intermediate genetic predisposition for {snp.trait.lower()}."
        )
    else:
        return (
            f"You carry two copies of the {risk} allele. "
            f"You may have a higher genetic predisposition related to {snp.trait.lower()}, "
            f"though lifestyle, diet, and regular medical monitoring play major roles."
        )


# ─────────────────────────────────────────────────────────────────────────────
# JSON OUTPUT
# ─────────────────────────────────────────────────────────────────────────────

def generate_organ_blood_json(result: OrganBloodAnalysisResult) -> dict:
    """
    Generate a JSON-serializable dict for the frontend organ & blood wellness report.
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
            "kidneyScore": result.category_scores.get("Kidney Health", 0.0),
            "liverScore": result.category_scores.get("Liver Function", 0.0),
            "cardiovascularScore": result.category_scores.get("Cardiovascular", 0.0),
            "bloodIronScore": result.category_scores.get("Blood & Iron", 0.0),
            "cholesterolScore": result.category_scores.get("Cholesterol", 0.0),
        },
        "categories": categories_list,
    }


# ─────────────────────────────────────────────────────────────────────────────
# MARKDOWN REPORT (optional, for compatibility)
# ─────────────────────────────────────────────────────────────────────────────

def generate_organ_blood_report(
    result: OrganBloodAnalysisResult,
    subject_name: str = "Subject",
) -> str:
    """Generate a Markdown organ & blood wellness report."""
    from datetime import datetime

    report = f"""# Organ & Blood Wellness Report

**Subject:** {subject_name}
**Generated:** {datetime.now().strftime("%B %d, %Y")}

---

## Overview

This report analyzes your genetic variants related to kidney and liver function,
cardiovascular health, blood and iron metabolism, and cholesterol regulation.
Genetics influences baseline predisposition; lifestyle, diet, and medical care
remain the dominant factors shaping long-term organ health.

| Metric | Value |
|--------|-------|
| Total SNPs Checked | {result.total_checked} |
| SNPs Found in Your Data | {result.found} |
| SNPs Not Available | {result.not_found} |
| Overall Wellness Score | {result.overall_score}/100 |

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

1. **This is not a diagnosis** - Genetic predisposition is only one component of
   organ health. Many environmental and lifestyle factors matter more.

2. **Consult medical professionals** - Any concerning findings should be discussed
   with your physician. Blood tests, imaging, and clinical evaluation provide
   the ground truth for organ function.

3. **Actionable lifestyle factors** - Diet, exercise, sleep, hydration, not smoking,
   moderate alcohol intake, and routine health screenings are the most effective
   levers for maintaining organ wellness.

---
*Generated by GeneHealth Analysis Platform*
"""

    return report
