"""
Longevity & Aging Pathways Analyzer
Analyzes genetic variants associated with aging pathways, longevity, and
healthspan based on validated GWAS findings and biological mechanisms.

EIGHT SECTIONS:
  1. APOE & Brain/Heart Aging       — rs429358, rs7412 (epsilon haplotype)
  2. FOXO3 & Stress Resistance      — rs2802292
  3. Telomere Maintenance           — TERT rs2736100, TERC rs10936599
  4. Cardiovascular Longevity       — 9p21 rs1333049, LDLR rs6511720, SORT1 rs646776
  5. Metabolic Health               — TCF7L2 rs7903146
  6. Inflammation & Immune Aging    — IL6 rs1800795, CRP rs1205
  7. Oxidative Stress Defense       — SOD2 rs4880, GPX1 rs1050450, PON1 rs662
  8. Cellular Repair & Maintenance  — SIRT1 rs7895833, MTHFR rs1801133, CETP rs5882,
                                      TP53 rs1042522, KLOTHO rs9536314

SCIENTIFIC BASIS:
  - Deelen et al. (2019) Nature Communications 10:3669 — GWAS meta-analysis
    confirming APOE and FOXO3 as genome-wide significant for longevity
  - Timmers et al. (2019) eLife 8:e39856 — Genomics of 1M parent lifespans
  - Timmers et al. (2020) Nature Communications 11:3570 — Multivariate
    genomic scan for healthspan/lifespan/longevity
  - Willcox et al. (2008) PNAS — FOXO3 longevity association
  - Kaplanis et al. (2018) Science — Heritability of lifespan ~7%
  - Ruby et al. (2018) Genetics — Heritability inflated by assortative mating

ETHICAL REQUIREMENTS:
  - APOE e4: Carries Alzheimer's risk implications. STRONG disclaimers required.
  - NEVER predict lifespan or time of death.
  - Frame as "aging pathways" not "longevity prediction."
  - Genetics accounts for ~10% of lifespan variation.

CAVEATS:
  - Effect sizes are small for most variants (OR 1.05-1.2), except APOE
  - Environment, lifestyle, healthcare access dominate longevity outcomes
  - Not a medical diagnosis — consult healthcare providers
"""

from typing import Dict, List, Tuple, Any, Optional
import json
from datetime import datetime, timezone


# ---------------------------------------------------------------------------
# UTILITIES
# ---------------------------------------------------------------------------

COMPLEMENT: Dict[str, str] = {"A": "T", "T": "A", "C": "G", "G": "C"}


def _complement(allele: str) -> str:
    return COMPLEMENT.get(allele.upper(), allele.upper())


def _count_allele(genotype: str, allele: str) -> int:
    genotype = genotype.upper().replace("-", "")
    allele = allele.upper()
    comp = _complement(allele)
    count = genotype.count(allele)
    if count == 0 and comp != allele:
        count = genotype.count(comp)
    return min(count, 2)


def _analyze_snp_list(
    variants: Dict[str, Tuple[str, str, str]],
    snp_list: List[Dict[str, Any]],
) -> List[Dict[str, Any]]:
    results = []
    for snp in snp_list:
        rsid = snp["rsid"].lower()
        raw = variants.get(rsid)
        if raw is None:
            chrom, pos, genotype = None, None, None
        elif isinstance(raw, str):
            chrom, pos, genotype = None, None, raw
        else:
            chrom, pos, genotype = raw

        if genotype:
            genotype = genotype.upper().replace("-", "")
            geno_data = snp["genotypes"].get(genotype)

            if not geno_data:
                comp_geno = "".join(_complement(a) for a in genotype)
                geno_data = snp["genotypes"].get(comp_geno)
                if not geno_data and len(genotype) == 2:
                    geno_data = snp["genotypes"].get(genotype[::-1])
                if not geno_data:
                    rev_comp = comp_geno[::-1]
                    geno_data = snp["genotypes"].get(rev_comp)

            if not geno_data:
                geno_data = {
                    "label": genotype,
                    "nickname": "Uncatalogued variant",
                    "interpretation": (
                        f"Your genotype {genotype} was detected but is not "
                        "catalogued for this SNP. It may be a rare allele."
                    ),
                    "score": 0.3,
                }

            results.append({
                **snp,
                "genotype": genotype,
                "genotypeData": geno_data,
                "found": True,
            })
        else:
            results.append({
                **snp,
                "genotype": "Not available",
                "genotypeData": {
                    "label": "N/A",
                    "nickname": "Not available",
                    "interpretation": (
                        "This SNP was not found in your genomic file. "
                        "Kits from different companies cover distinct sets of variants."
                    ),
                    "score": None,
                },
                "found": False,
            })

        results[-1]["category"] = snp.get("category", "")

    return results


# ---------------------------------------------------------------------------
# APOE EPSILON HAPLOTYPE DETERMINATION
# ---------------------------------------------------------------------------

def _determine_apoe_epsilon(variants: Dict) -> Dict[str, Any]:
    """
    Determine APOE epsilon genotype from rs429358 (C/T) and rs7412 (C/T).

    Haplotypes:
      epsilon2 = T at rs429358 + T at rs7412
      epsilon3 = T at rs429358 + C at rs7412  (most common)
      epsilon4 = C at rs429358 + C at rs7412

    Returns dict with genotype string, risk level, and interpretation.
    """
    raw_429 = variants.get("rs429358")
    raw_7412 = variants.get("rs7412")

    def _get_geno(raw):
        if raw is None:
            return None
        if isinstance(raw, str):
            return raw.upper().replace("-", "")
        return raw[2].upper().replace("-", "") if raw[2] else None

    g429 = _get_geno(raw_429)
    g7412 = _get_geno(raw_7412)

    if not g429 or not g7412:
        return {
            "genotype": "Not determined",
            "alleles": [],
            "risk_level": "unknown",
            "interpretation": (
                "APOE epsilon genotype could not be determined. One or both required "
                "SNPs (rs429358, rs7412) were not found in your DNA file. This is common "
                "with some genotyping providers."
            ),
            "warnings": [],
            "found": False,
            "score": 0.5,
        }

    # Determine alleles on each chromosome
    # rs429358: T=e2/e3, C=e4
    # rs7412:   T=e2,    C=e3/e4
    def _alleles_at(g429_allele: str, g7412_allele: str) -> str:
        # Handle complement (some arrays report on opposite strand)
        t429 = g429_allele.upper()
        t7412 = g7412_allele.upper()
        # Normalize to forward strand
        if t429 in ("A", "G"):
            t429 = _complement(t429)
        if t7412 in ("A", "G"):
            t7412 = _complement(t7412)

        if t429 == "T" and t7412 == "T":
            return "e2"
        elif t429 == "T" and t7412 == "C":
            return "e3"
        elif t429 == "C" and t7412 == "C":
            return "e4"
        else:
            return "e?"  # C/T at 429, T at 7412 doesn't map to standard

    # Get individual alleles
    a1_429, a2_429 = g429[0], g429[1] if len(g429) >= 2 else g429[0]
    a1_7412, a2_7412 = g7412[0], g7412[1] if len(g7412) >= 2 else g7412[0]

    allele1 = _alleles_at(a1_429, a1_7412)
    allele2 = _alleles_at(a2_429, a2_7412)

    # Sort for consistent display
    alleles = sorted([allele1, allele2])
    genotype_str = f"{alleles[0]}/{alleles[1]}"

    # Score and interpretation
    e4_count = alleles.count("e4")
    e2_count = alleles.count("e2")

    warnings = []

    if e4_count == 2:
        score = 0.1
        risk_level = "high_risk"
        interpretation = (
            "You carry two copies of the APOE epsilon-4 allele (e4/e4). This is the genotype most "
            "strongly associated with reduced longevity in genetic studies. APOE e4 homozygotes have "
            "approximately 50% reduced odds of reaching age 100 compared to e3/e3 carriers. "
            "This genotype is also associated with significantly increased risk for late-onset "
            "Alzheimer's disease and cardiovascular disease. However, many e4/e4 carriers live "
            "long, healthy lives — genetics is one factor among many."
        )
        warnings = [
            "IMPORTANT: APOE e4/e4 is associated with significantly increased risk for "
            "late-onset Alzheimer's disease (10-15x higher than e3/e3). This is NOT a diagnosis. "
            "Many e4/e4 carriers never develop Alzheimer's. Consult a genetic counselor or "
            "healthcare provider to discuss these findings.",
            "Lifestyle factors such as regular exercise, Mediterranean diet, cognitive engagement, "
            "quality sleep, and cardiovascular health management can significantly modify risk "
            "regardless of APOE genotype."
        ]
    elif e4_count == 1:
        score = 0.3
        risk_level = "elevated"
        interpretation = (
            "You carry one copy of the APOE epsilon-4 allele. Heterozygous e4 carriers have "
            "moderately reduced odds of extreme longevity compared to e3/e3. This genotype is "
            "associated with increased (but not deterministic) risk for late-onset Alzheimer's "
            "disease and cardiovascular disease. Many e3/e4 carriers live well into their 80s and 90s."
        )
        warnings = [
            "APOE e3/e4 is associated with moderately increased risk for late-onset Alzheimer's "
            "disease (2-3x higher than e3/e3). This is NOT a diagnosis. Consult a healthcare "
            "provider if you have concerns.",
        ]
    elif e2_count == 2:
        score = 0.9
        risk_level = "protective"
        interpretation = (
            "You carry two copies of the APOE epsilon-2 allele (e2/e2). This is the genotype most "
            "strongly associated with longevity. APOE e2 carriers have approximately 1.5-2x "
            "increased odds of reaching extreme old age. e2 is associated with lower LDL cholesterol "
            "and reduced cardiovascular and Alzheimer's risk."
        )
    elif e2_count == 1 and e4_count == 0:
        score = 0.75
        risk_level = "favorable"
        interpretation = (
            "You carry one APOE epsilon-2 allele and one epsilon-3 allele. The e2 allele is "
            "associated with increased longevity, lower LDL cholesterol, and reduced risk for "
            "cardiovascular and neurodegenerative disease. This is a favorable genotype for aging."
        )
    elif "e?" in alleles:
        score = 0.5
        risk_level = "unknown"
        interpretation = (
            "Your APOE genotype could not be fully resolved. This may be due to strand orientation "
            "differences in your genotyping data. The raw alleles were detected but do not map to "
            "standard epsilon haplotypes."
        )
    else:  # e3/e3
        score = 0.55
        risk_level = "average"
        interpretation = (
            "You carry the most common APOE genotype (e3/e3), found in approximately 60% of the "
            "population. This genotype confers average longevity odds — neither the protective "
            "effect of e2 nor the risk associated with e4. Your APOE genotype is neutral for aging."
        )

    return {
        "genotype": genotype_str,
        "alleles": alleles,
        "risk_level": risk_level,
        "interpretation": interpretation,
        "warnings": warnings,
        "found": True,
        "score": score,
        "rs429358": g429,
        "rs7412": g7412,
    }


# ---------------------------------------------------------------------------
# SNP DATABASE: FOXO3 — Stress Resistance
# ---------------------------------------------------------------------------

FOXO3_SNPS: List[Dict[str, Any]] = [
    {
        "rsid": "rs2802292",
        "gene": "FOXO3",
        "name": "Longevity-Associated Stress Resistance",
        "chromosome": "6",
        "position": 108985784,
        "effectSize": "moderate",
        "category": "foxo3StressResistance",
        "genotypes": {
            "TT": {
                "label": "T/T",
                "nickname": "Longevity Variant (Homozygous)",
                "interpretation": (
                    "Two copies of the T allele at FOXO3 rs2802292. This is the genotype most "
                    "consistently associated with longevity across multiple populations (European, "
                    "Japanese, Chinese, Hawaiian). The T allele enhances FOXO3 transcription factor "
                    "activity, which promotes cellular stress resistance, DNA repair, apoptosis of "
                    "damaged cells, and antioxidant defense. OR ~1.17-1.26 for longevity (Deelen 2019)."
                ),
                "score": 0.9,
            },
            "GT": {
                "label": "G/T",
                "nickname": "Longevity Variant (Heterozygous)",
                "interpretation": (
                    "One copy of the longevity-associated T allele at FOXO3. Intermediate "
                    "enhancement of FOXO3 activity. This genotype still confers a modest longevity "
                    "advantage compared to G/G homozygotes."
                ),
                "score": 0.65,
            },
            "GG": {
                "label": "G/G",
                "nickname": "Common Variant",
                "interpretation": (
                    "The most common genotype at FOXO3 rs2802292. Standard FOXO3 transcription "
                    "factor activity. No enhanced longevity association from this variant, but FOXO3 "
                    "pathway activity is influenced by many factors including caloric restriction, "
                    "exercise, and fasting."
                ),
                "score": 0.4,
            },
        },
    },
]


# ---------------------------------------------------------------------------
# SNP DATABASE: Telomere Maintenance
# ---------------------------------------------------------------------------

TELOMERE_SNPS: List[Dict[str, Any]] = [
    {
        "rsid": "rs2736100",
        "gene": "TERT",
        "name": "Telomerase Reverse Transcriptase",
        "chromosome": "5",
        "position": 1286516,
        "effectSize": "moderate",
        "category": "telomereMaintenance",
        "genotypes": {
            "CC": {
                "label": "C/C",
                "nickname": "Enhanced Telomere Maintenance",
                "interpretation": (
                    "Two copies of the C allele at TERT, associated with longer leukocyte telomere "
                    "length. TERT encodes the catalytic subunit of telomerase, the enzyme that "
                    "maintains telomere length. Longer telomeres are generally associated with "
                    "slower cellular aging, though the relationship is complex."
                ),
                "score": 0.85,
            },
            "AC": {
                "label": "A/C",
                "nickname": "Intermediate Telomere Maintenance",
                "interpretation": (
                    "One copy of the C allele. Intermediate telomere length association. "
                    "Telomere length is also influenced by lifestyle factors such as stress, "
                    "exercise, sleep quality, and diet."
                ),
                "score": 0.6,
            },
            "AA": {
                "label": "A/A",
                "nickname": "Standard Telomere Maintenance",
                "interpretation": (
                    "Common genotype associated with average telomere length. Telomere shortening "
                    "rate is strongly modifiable through lifestyle — regular exercise, stress "
                    "management, and adequate sleep can preserve telomere length."
                ),
                "score": 0.35,
            },
        },
    },
    {
        "rsid": "rs10936599",
        "gene": "TERC",
        "name": "Telomerase RNA Component",
        "chromosome": "3",
        "position": 169774313,
        "effectSize": "moderate",
        "category": "telomereMaintenance",
        "genotypes": {
            "CC": {
                "label": "C/C",
                "nickname": "Longer Telomeres",
                "interpretation": (
                    "Associated with longer telomere length. TERC provides the RNA template "
                    "for telomere extension. This genotype supports more efficient telomere "
                    "maintenance during cell division."
                ),
                "score": 0.8,
            },
            "CT": {
                "label": "C/T",
                "nickname": "Average Telomere Length",
                "interpretation": (
                    "Intermediate effect on telomere length. One copy of each allele provides "
                    "moderate telomerase RNA template efficiency."
                ),
                "score": 0.55,
            },
            "TT": {
                "label": "T/T",
                "nickname": "Shorter Telomeres",
                "interpretation": (
                    "Associated with shorter leukocyte telomere length (~0.08 SD per allele). "
                    "Shorter telomeres are associated with accelerated cellular aging, but "
                    "this is modifiable through healthy lifestyle choices."
                ),
                "score": 0.3,
            },
        },
    },
]


# ---------------------------------------------------------------------------
# SNP DATABASE: Cardiovascular Longevity
# ---------------------------------------------------------------------------

CARDIOVASCULAR_SNPS: List[Dict[str, Any]] = [
    {
        "rsid": "rs1333049",
        "gene": "CDKN2B-AS1 (9p21)",
        "name": "Cellular Senescence & Cardiovascular Risk",
        "chromosome": "9",
        "position": 22125503,
        "effectSize": "moderate",
        "category": "cardiovascularLongevity",
        "genotypes": {
            "GG": {
                "label": "G/G",
                "nickname": "Favorable Cardiovascular Aging",
                "interpretation": (
                    "The 9p21 locus is one of the most replicated cardiovascular risk regions "
                    "in the genome. This genotype is associated with lower risk for coronary "
                    "artery disease and healthier vascular aging. The region regulates CDKN2B, "
                    "which controls cellular senescence."
                ),
                "score": 0.8,
            },
            "GC": {
                "label": "G/C",
                "nickname": "Average Cardiovascular Risk",
                "interpretation": (
                    "Heterozygous at the 9p21 locus. Intermediate cardiovascular risk. "
                    "Standard lifestyle modifications (exercise, diet, not smoking) are "
                    "effective at managing cardiovascular health."
                ),
                "score": 0.5,
            },
            "CC": {
                "label": "C/C",
                "nickname": "Elevated Cardiovascular Risk",
                "interpretation": (
                    "Two copies of the risk allele at 9p21. Associated with ~25-30% increased "
                    "risk for coronary artery disease. This is the most replicated genetic "
                    "cardiovascular risk factor. Proactive cardiovascular health management "
                    "can significantly modify this risk."
                ),
                "score": 0.25,
            },
        },
    },
    {
        "rsid": "rs6511720",
        "gene": "LDLR",
        "name": "LDL Receptor — Cholesterol Clearance",
        "chromosome": "19",
        "position": 11202306,
        "effectSize": "moderate",
        "category": "cardiovascularLongevity",
        "genotypes": {
            "TT": {
                "label": "T/T",
                "nickname": "Enhanced LDL Clearance",
                "interpretation": (
                    "Two copies of the T allele at LDLR, associated with lower LDL cholesterol "
                    "levels and improved cardiovascular longevity. Identified in the Timmers 2019 "
                    "healthspan GWAS."
                ),
                "score": 0.85,
            },
            "GT": {
                "label": "G/T",
                "nickname": "Moderate LDL Clearance",
                "interpretation": (
                    "One copy of each allele. Intermediate LDL clearance efficiency. "
                    "Diet and exercise significantly influence cholesterol levels."
                ),
                "score": 0.55,
            },
            "GG": {
                "label": "G/G",
                "nickname": "Standard LDL Clearance",
                "interpretation": (
                    "Common genotype with standard LDL receptor activity. Cholesterol management "
                    "through diet, exercise, and medical guidance as needed."
                ),
                "score": 0.35,
            },
        },
    },
    {
        "rsid": "rs646776",
        "gene": "CELSR2/SORT1",
        "name": "Lipid Metabolism & Sorting",
        "chromosome": "1",
        "position": 109818530,
        "effectSize": "moderate",
        "category": "cardiovascularLongevity",
        "genotypes": {
            "TT": {
                "label": "T/T",
                "nickname": "Lower LDL Lipid Profile",
                "interpretation": (
                    "Associated with lower LDL cholesterol through enhanced hepatic lipid sorting. "
                    "The SORT1 locus is one of the strongest GWAS signals for LDL levels and "
                    "is associated with improved healthspan."
                ),
                "score": 0.8,
            },
            "CT": {
                "label": "C/T",
                "nickname": "Moderate Lipid Profile",
                "interpretation": (
                    "Intermediate effect on LDL cholesterol. One copy of each allele."
                ),
                "score": 0.55,
            },
            "CC": {
                "label": "C/C",
                "nickname": "Standard Lipid Metabolism",
                "interpretation": (
                    "Common genotype. Standard hepatic lipid sorting. Lifestyle and dietary "
                    "choices have significant impact on lipid levels."
                ),
                "score": 0.35,
            },
        },
    },
]


# ---------------------------------------------------------------------------
# SNP DATABASE: Metabolic Health
# ---------------------------------------------------------------------------

METABOLIC_SNPS: List[Dict[str, Any]] = [
    {
        "rsid": "rs7903146",
        "gene": "TCF7L2",
        "name": "Metabolic Health & Insulin Signaling",
        "chromosome": "10",
        "position": 114758349,
        "effectSize": "large",
        "category": "metabolicHealth",
        "genotypes": {
            "CC": {
                "label": "C/C",
                "nickname": "Favorable Metabolic Aging",
                "interpretation": (
                    "No risk alleles at TCF7L2, the strongest known genetic risk factor for "
                    "type 2 diabetes. This genotype supports healthy insulin signaling and "
                    "glucose metabolism, contributing to metabolic healthspan."
                ),
                "score": 0.85,
            },
            "CT": {
                "label": "C/T",
                "nickname": "Moderate Metabolic Risk",
                "interpretation": (
                    "One copy of the T risk allele. Modestly increased risk for type 2 diabetes "
                    "(OR ~1.4 per allele). Maintaining healthy weight, regular exercise, and "
                    "balanced diet are particularly important for metabolic health."
                ),
                "score": 0.45,
            },
            "TT": {
                "label": "T/T",
                "nickname": "Elevated Metabolic Risk",
                "interpretation": (
                    "Two copies of the T risk allele at TCF7L2. This genotype is associated with "
                    "~1.8x increased risk for type 2 diabetes through impaired beta-cell function "
                    "and insulin secretion. Proactive metabolic health management (weight control, "
                    "exercise, blood sugar monitoring) can significantly reduce this risk."
                ),
                "score": 0.2,
            },
        },
    },
]


# ---------------------------------------------------------------------------
# SNP DATABASE: Inflammation & Immune Aging
# ---------------------------------------------------------------------------

INFLAMMATION_SNPS: List[Dict[str, Any]] = [
    {
        "rsid": "rs1800795",
        "gene": "IL6",
        "name": "Interleukin-6 — Inflammatory Signaling",
        "chromosome": "7",
        "position": 22766246,
        "effectSize": "moderate",
        "category": "inflammationImmuneAging",
        "genotypes": {
            "CC": {
                "label": "C/C (-174C)",
                "nickname": "Lower IL-6 Production",
                "interpretation": (
                    "Two copies of the C allele at the IL6 promoter. Associated with lower "
                    "baseline interleukin-6 production. Lower chronic inflammation ('inflammaging') "
                    "is consistently associated with healthier aging and reduced risk for "
                    "age-related diseases."
                ),
                "score": 0.8,
            },
            "GC": {
                "label": "G/C",
                "nickname": "Intermediate IL-6 Levels",
                "interpretation": (
                    "Heterozygous genotype. Intermediate IL-6 production. Anti-inflammatory "
                    "lifestyle factors (omega-3 fatty acids, exercise, stress reduction) can "
                    "modulate inflammatory signaling."
                ),
                "score": 0.55,
            },
            "GG": {
                "label": "G/G (-174G)",
                "nickname": "Higher IL-6 Production",
                "interpretation": (
                    "Two copies of the G allele. Associated with higher baseline IL-6 production "
                    "and potentially greater chronic inflammation. Anti-inflammatory diet and "
                    "regular exercise are particularly beneficial."
                ),
                "score": 0.3,
            },
        },
    },
    {
        "rsid": "rs1205",
        "gene": "CRP",
        "name": "C-Reactive Protein — Inflammation Biomarker",
        "chromosome": "1",
        "position": 159682233,
        "effectSize": "moderate",
        "category": "inflammationImmuneAging",
        "genotypes": {
            "TT": {
                "label": "T/T",
                "nickname": "Lower CRP Levels",
                "interpretation": (
                    "Associated with lower baseline C-reactive protein, a key marker of systemic "
                    "inflammation. Lower CRP is linked to reduced cardiovascular risk and healthier "
                    "aging across multiple GWAS studies."
                ),
                "score": 0.8,
            },
            "CT": {
                "label": "C/T",
                "nickname": "Intermediate CRP Levels",
                "interpretation": (
                    "One copy of each allele. Intermediate CRP levels. Lifestyle factors "
                    "strongly influence CRP — weight management and exercise are effective."
                ),
                "score": 0.55,
            },
            "CC": {
                "label": "C/C",
                "nickname": "Higher CRP Levels",
                "interpretation": (
                    "Associated with higher baseline CRP levels. This does not mean you have "
                    "active inflammation — it reflects a genetic tendency. Regular blood tests "
                    "can monitor actual CRP levels, and lifestyle modifications are effective."
                ),
                "score": 0.3,
            },
        },
    },
]


# ---------------------------------------------------------------------------
# SNP DATABASE: Oxidative Stress Defense
# ---------------------------------------------------------------------------

OXIDATIVE_SNPS: List[Dict[str, Any]] = [
    {
        "rsid": "rs4880",
        "gene": "SOD2",
        "name": "Superoxide Dismutase 2 — Mitochondrial Antioxidant",
        "chromosome": "6",
        "position": 160113872,
        "effectSize": "moderate",
        "category": "oxidativeStressDefense",
        "genotypes": {
            "TT": {
                "label": "T/T (Ala/Ala)",
                "nickname": "Efficient Mitochondrial Defense",
                "interpretation": (
                    "Ala/Ala genotype — the alanine form is efficiently transported into "
                    "the mitochondrial matrix where SOD2 neutralizes superoxide radicals. "
                    "This supports more effective mitochondrial antioxidant defense."
                ),
                "score": 0.8,
            },
            "CT": {
                "label": "C/T (Ala/Val)",
                "nickname": "Moderate Mitochondrial Defense",
                "interpretation": (
                    "One copy of each allele. Intermediate SOD2 mitochondrial transport "
                    "efficiency. Adequate dietary antioxidants support this genotype well."
                ),
                "score": 0.55,
            },
            "CC": {
                "label": "C/C (Val/Val)",
                "nickname": "Reduced Mitochondrial Import",
                "interpretation": (
                    "Val/Val genotype — the valine form has reduced import into the "
                    "mitochondrial matrix, potentially leading to less efficient neutralization "
                    "of superoxide radicals. Antioxidant-rich diet (fruits, vegetables) is "
                    "particularly beneficial."
                ),
                "score": 0.3,
            },
        },
    },
    {
        "rsid": "rs1050450",
        "gene": "GPX1",
        "name": "Glutathione Peroxidase 1 — Selenium Antioxidant",
        "chromosome": "3",
        "position": 49394834,
        "effectSize": "small",
        "category": "oxidativeStressDefense",
        "genotypes": {
            "CC": {
                "label": "C/C (Pro/Pro)",
                "nickname": "Standard GPX1 Activity",
                "interpretation": (
                    "Pro/Pro genotype at GPX1. Standard glutathione peroxidase activity. "
                    "This selenium-dependent enzyme neutralizes hydrogen peroxide and lipid "
                    "hydroperoxides."
                ),
                "score": 0.65,
            },
            "CT": {
                "label": "C/T (Pro/Leu)",
                "nickname": "Reduced GPX1 Activity",
                "interpretation": (
                    "One copy of the Leu variant. Slightly reduced GPX1 enzyme activity. "
                    "Adequate selenium intake supports this genotype."
                ),
                "score": 0.45,
            },
            "TT": {
                "label": "T/T (Leu/Leu)",
                "nickname": "Lower GPX1 Activity",
                "interpretation": (
                    "Leu/Leu genotype with more significantly reduced GPX1 activity. "
                    "Ensuring adequate selenium intake (Brazil nuts, seafood, whole grains) "
                    "is particularly important."
                ),
                "score": 0.3,
            },
        },
    },
    {
        "rsid": "rs662",
        "gene": "PON1",
        "name": "Paraoxonase 1 — Lipid Antioxidant",
        "chromosome": "7",
        "position": 95308134,
        "effectSize": "moderate",
        "category": "oxidativeStressDefense",
        "genotypes": {
            "AA": {
                "label": "A/A (Gln/Gln)",
                "nickname": "Higher PON1 Activity",
                "interpretation": (
                    "Gln/Gln genotype. Higher paraoxonase activity for certain substrates. "
                    "PON1 protects LDL from oxidation and is associated with cardiovascular "
                    "protection."
                ),
                "score": 0.75,
            },
            "AG": {
                "label": "A/G (Gln/Arg)",
                "nickname": "Intermediate PON1 Activity",
                "interpretation": (
                    "Heterozygous genotype. Intermediate PON1 enzyme activity. This enzyme "
                    "is also modifiable by diet — olive oil and polyphenols enhance PON1."
                ),
                "score": 0.55,
            },
            "GG": {
                "label": "G/G (Arg/Arg)",
                "nickname": "Altered PON1 Specificity",
                "interpretation": (
                    "Arg/Arg genotype — different substrate specificity for PON1. "
                    "A Mediterranean diet rich in olive oil and polyphenols can help "
                    "support PON1 function."
                ),
                "score": 0.35,
            },
        },
    },
]


# ---------------------------------------------------------------------------
# SNP DATABASE: Cellular Repair & Maintenance
# ---------------------------------------------------------------------------

CELLULAR_REPAIR_SNPS: List[Dict[str, Any]] = [
    {
        "rsid": "rs7895833",
        "gene": "SIRT1",
        "name": "Sirtuin 1 — Cellular Aging Regulator",
        "chromosome": "10",
        "position": 69644218,
        "effectSize": "small",
        "category": "cellularRepairMaintenance",
        "genotypes": {
            "AA": {
                "label": "A/A",
                "nickname": "Enhanced SIRT1 Expression",
                "interpretation": (
                    "Associated with higher SIRT1 expression. SIRT1 is a NAD+-dependent "
                    "deacetylase that promotes DNA repair, reduces inflammation, and enhances "
                    "metabolic efficiency. It is activated by caloric restriction and exercise."
                ),
                "score": 0.8,
            },
            "AG": {
                "label": "A/G",
                "nickname": "Moderate SIRT1 Expression",
                "interpretation": (
                    "Heterozygous genotype. Intermediate SIRT1 activity. Fasting, exercise, "
                    "and resveratrol-containing foods may support SIRT1 pathway activation."
                ),
                "score": 0.55,
            },
            "GG": {
                "label": "G/G",
                "nickname": "Standard SIRT1 Expression",
                "interpretation": (
                    "Common genotype. Standard SIRT1 expression. The SIRT1 pathway is highly "
                    "responsive to lifestyle — caloric restriction and intermittent fasting "
                    "are powerful activators regardless of genotype."
                ),
                "score": 0.35,
            },
        },
    },
    {
        "rsid": "rs1801133",
        "gene": "MTHFR",
        "name": "Methylation & Folate Metabolism",
        "chromosome": "1",
        "position": 11856378,
        "effectSize": "moderate",
        "category": "cellularRepairMaintenance",
        "genotypes": {
            "CC": {
                "label": "C/C (677C — wild type)",
                "nickname": "Full MTHFR Activity",
                "interpretation": (
                    "Normal MTHFR enzyme activity (~100%). Efficient conversion of folate to "
                    "methylfolate, supporting DNA methylation, homocysteine metabolism, and "
                    "epigenetic maintenance — all important for healthy cellular aging."
                ),
                "score": 0.8,
            },
            "CT": {
                "label": "C/T (677CT)",
                "nickname": "Reduced MTHFR Activity (~65%)",
                "interpretation": (
                    "One copy of the T variant. MTHFR enzyme activity reduced to ~65%. "
                    "Adequate folate intake (leafy greens or methylfolate supplementation) "
                    "can compensate. Monitor homocysteine levels."
                ),
                "score": 0.5,
            },
            "TT": {
                "label": "T/T (677TT)",
                "nickname": "Significantly Reduced MTHFR (~30%)",
                "interpretation": (
                    "Two copies of the T variant. MTHFR activity reduced to ~30%. May lead "
                    "to elevated homocysteine (a cardiovascular risk factor) if folate intake "
                    "is inadequate. Methylfolate supplementation and B-vitamin intake are "
                    "particularly important for this genotype."
                ),
                "score": 0.25,
            },
        },
    },
    {
        "rsid": "rs5882",
        "gene": "CETP",
        "name": "HDL Cholesterol & Longevity",
        "chromosome": "16",
        "position": 57017292,
        "effectSize": "moderate",
        "category": "cellularRepairMaintenance",
        "genotypes": {
            "GG": {
                "label": "G/G (Ile/Ile)",
                "nickname": "Standard CETP Activity",
                "interpretation": (
                    "Standard cholesteryl ester transfer protein activity. Normal HDL/LDL balance."
                ),
                "score": 0.45,
            },
            "AG": {
                "label": "A/G (Ile/Val)",
                "nickname": "Moderately Favorable HDL Profile",
                "interpretation": (
                    "One copy of the Val variant. Reduced CETP activity leads to larger, "
                    "more protective HDL particles. Found at higher frequency in some "
                    "centenarian populations (Barzilai, Ashkenazi studies)."
                ),
                "score": 0.65,
            },
            "AA": {
                "label": "A/A (Val/Val)",
                "nickname": "Longevity-Associated HDL Profile",
                "interpretation": (
                    "Two copies of the Val variant (Ile405Val). Significantly reduced CETP "
                    "activity, resulting in larger HDL particles and favorable lipid profile. "
                    "This genotype was enriched in Ashkenazi Jewish centenarians (Barzilai 2003)."
                ),
                "score": 0.85,
            },
        },
    },
    {
        "rsid": "rs1042522",
        "gene": "TP53",
        "name": "Tumor Protein p53 — DNA Damage Response",
        "chromosome": "17",
        "position": 7579472,
        "effectSize": "small",
        "category": "cellularRepairMaintenance",
        "genotypes": {
            "GG": {
                "label": "G/G (Arg/Arg)",
                "nickname": "Stronger Apoptotic Response",
                "interpretation": (
                    "Arg/Arg at codon 72 of TP53. The Arg form induces apoptosis more "
                    "efficiently, which is protective against cancer but may accelerate "
                    "tissue aging through increased cell death. A balanced trade-off."
                ),
                "score": 0.55,
            },
            "GC": {
                "label": "G/C (Arg/Pro)",
                "nickname": "Balanced p53 Response",
                "interpretation": (
                    "Heterozygous genotype. Balance between apoptotic efficiency (Arg) "
                    "and DNA repair/cell cycle arrest (Pro). Some studies suggest this "
                    "heterozygous balance may be favorable for longevity."
                ),
                "score": 0.65,
            },
            "CC": {
                "label": "C/C (Pro/Pro)",
                "nickname": "Stronger DNA Repair Response",
                "interpretation": (
                    "Pro/Pro at codon 72. The Pro form is better at inducing cell cycle "
                    "arrest and DNA repair rather than apoptosis. May preserve tissue "
                    "integrity but with potentially reduced cancer surveillance."
                ),
                "score": 0.5,
            },
        },
    },
    {
        "rsid": "rs9536314",
        "gene": "KLOTHO",
        "name": "Klotho — Anti-Aging Hormone",
        "chromosome": "13",
        "position": 33628138,
        "effectSize": "moderate",
        "category": "cellularRepairMaintenance",
        "genotypes": {
            "TT": {
                "label": "T/T",
                "nickname": "Standard Klotho Levels",
                "interpretation": (
                    "Common genotype. Standard circulating Klotho levels. Klotho is a "
                    "transmembrane protein that functions as an anti-aging hormone, "
                    "suppressing oxidative stress and regulating calcium/phosphate metabolism."
                ),
                "score": 0.45,
            },
            "GT": {
                "label": "G/T (KL-VS heterozygous)",
                "nickname": "Enhanced Klotho — Longevity Variant",
                "interpretation": (
                    "One copy of the KL-VS variant. Heterozygous carriers have higher "
                    "circulating Klotho and are associated with better cognitive function "
                    "and increased longevity. The KL-VS variant enhances Klotho secretion."
                ),
                "score": 0.8,
            },
            "GG": {
                "label": "G/G (KL-VS homozygous)",
                "nickname": "Paradoxical — Reduced Benefit",
                "interpretation": (
                    "Two copies of KL-VS. Paradoxically, homozygous carriers do NOT show "
                    "the same longevity benefit as heterozygous carriers. The mechanism "
                    "is not fully understood, but Klotho levels may not increase proportionally "
                    "with two copies."
                ),
                "score": 0.4,
            },
        },
    },
]


# ---------------------------------------------------------------------------
# SCORING ENGINE
# ---------------------------------------------------------------------------

SECTION_WEIGHTS = {
    "apoeBrainHeart": 0.25,
    "foxo3StressResistance": 0.15,
    "telomereMaintenance": 0.12,
    "cardiovascularLongevity": 0.15,
    "metabolicHealth": 0.10,
    "inflammationImmuneAging": 0.08,
    "oxidativeStressDefense": 0.08,
    "cellularRepairMaintenance": 0.07,
}


def _section_score(results: List[Dict]) -> Tuple[float, float]:
    """Return (sum_of_scores, max_possible_score) from analyzed SNP results."""
    scores = []
    for r in results:
        s = r.get("genotypeData", {}).get("score")
        if s is not None:
            scores.append(s)
    if not scores:
        return 0.5, 1.0  # Default to average if no data
    return sum(scores) / len(scores), 1.0


def _score_label(score: float) -> str:
    if score >= 0.75:
        return "Favorable"
    elif score >= 0.55:
        return "Above Average"
    elif score >= 0.40:
        return "Average"
    elif score >= 0.25:
        return "Below Average"
    else:
        return "Unfavorable"


def _overall_score_label(score: int) -> str:
    if score >= 75:
        return "Favorable Genetic Aging Profile"
    elif score >= 60:
        return "Above Average Aging Profile"
    elif score >= 45:
        return "Average Aging Profile"
    elif score >= 30:
        return "Below Average Aging Profile"
    else:
        return "Genetic Risk Factors Identified"


# ---------------------------------------------------------------------------
# ACTIONABLE INSIGHTS GENERATOR
# ---------------------------------------------------------------------------

def _generate_insights(
    apoe: Dict, sections: Dict[str, List[Dict]]
) -> List[Dict[str, Any]]:
    insights = []

    # APOE-based insights
    if apoe.get("risk_level") in ("high_risk", "elevated"):
        insights.append({
            "category": "Brain Health",
            "recommendation": (
                "Prioritize cardiovascular exercise (150+ min/week), Mediterranean diet, "
                "quality sleep (7-9 hours), and cognitive engagement. These interventions "
                "have the strongest evidence for modifying APOE-related risk."
            ),
            "basedOn": f"APOE {apoe['genotype']}",
            "evidenceLevel": "Strong (multiple RCTs and cohort studies)",
            "priority": "high",
        })

    # Telomere-based insights
    tel_scores = [
        r["genotypeData"]["score"] for r in sections.get("telomere", [])
        if r["found"] and r["genotypeData"].get("score") is not None
    ]
    if tel_scores and sum(tel_scores) / len(tel_scores) < 0.5:
        insights.append({
            "category": "Telomere Protection",
            "recommendation": (
                "Your telomere maintenance variants suggest focusing on stress management, "
                "regular aerobic exercise, adequate sleep, and an antioxidant-rich diet. "
                "These have been shown to slow telomere shortening."
            ),
            "basedOn": "TERT and TERC variants",
            "evidenceLevel": "Moderate (observational studies)",
            "priority": "moderate",
        })

    # Inflammation insights
    inf_scores = [
        r["genotypeData"]["score"] for r in sections.get("inflammation", [])
        if r["found"] and r["genotypeData"].get("score") is not None
    ]
    if inf_scores and sum(inf_scores) / len(inf_scores) < 0.5:
        insights.append({
            "category": "Anti-Inflammatory Lifestyle",
            "recommendation": (
                "Your genetic variants suggest a tendency toward higher inflammatory markers. "
                "An anti-inflammatory diet (rich in omega-3, turmeric, leafy greens, berries), "
                "regular exercise, and stress reduction can significantly lower inflammation."
            ),
            "basedOn": "IL6 and CRP variants",
            "evidenceLevel": "Strong (interventional evidence)",
            "priority": "moderate",
        })

    # Metabolic insights
    met_scores = [
        r["genotypeData"]["score"] for r in sections.get("metabolic", [])
        if r["found"] and r["genotypeData"].get("score") is not None
    ]
    if met_scores and sum(met_scores) / len(met_scores) < 0.5:
        insights.append({
            "category": "Metabolic Health",
            "recommendation": (
                "Your TCF7L2 variant is associated with increased metabolic risk. Maintaining "
                "a healthy weight, regular physical activity, limiting refined carbohydrates, "
                "and monitoring blood glucose are evidence-based protective strategies."
            ),
            "basedOn": "TCF7L2 rs7903146",
            "evidenceLevel": "Strong (GWAS replicated, interventional evidence)",
            "priority": "high",
        })

    # Oxidative stress insights
    ox_scores = [
        r["genotypeData"]["score"] for r in sections.get("oxidative", [])
        if r["found"] and r["genotypeData"].get("score") is not None
    ]
    if ox_scores and sum(ox_scores) / len(ox_scores) < 0.5:
        insights.append({
            "category": "Antioxidant Support",
            "recommendation": (
                "Your oxidative stress defense variants suggest supporting your antioxidant "
                "systems through diet: colorful fruits and vegetables, Brazil nuts (selenium), "
                "olive oil, and green tea."
            ),
            "basedOn": "SOD2, GPX1, and PON1 variants",
            "evidenceLevel": "Moderate (mechanistic + observational)",
            "priority": "low",
        })

    # MTHFR insight
    for r in sections.get("cellular", []):
        if r["rsid"] == "rs1801133" and r["found"]:
            gt = r["genotype"]
            if gt in ("TT", "AA"):  # TT on forward strand, AA on complement
                insights.append({
                    "category": "Methylation Support",
                    "recommendation": (
                        "Your MTHFR 677TT genotype significantly reduces folate metabolism. "
                        "Consider methylfolate (not folic acid) supplementation, increase leafy "
                        "green intake, and monitor homocysteine levels with your doctor."
                    ),
                    "basedOn": "MTHFR rs1801133 TT",
                    "evidenceLevel": "Strong (well-established biochemistry)",
                    "priority": "high",
                })
                break

    # General longevity insight (always include)
    insights.append({
        "category": "General Longevity",
        "recommendation": (
            "The strongest evidence-based longevity strategies apply regardless of genotype: "
            "regular physical activity (both aerobic and strength training), a whole-food diet "
            "rich in plants, maintaining social connections, quality sleep, moderate alcohol or "
            "none, not smoking, and regular preventive healthcare."
        ),
        "basedOn": "Overall aging pathway analysis",
        "evidenceLevel": "Very Strong (meta-analyses of lifestyle interventions)",
        "priority": "high",
    })

    return insights


# ---------------------------------------------------------------------------
# MAIN ANALYZER
# ---------------------------------------------------------------------------

def analyze_longevity_aging(variants: Dict) -> Dict[str, Any]:
    """
    Main entry point. Analyzes genetic variants for longevity-associated pathways.

    Args:
        variants: Dict mapping rsid -> (chrom, pos, genotype) or rsid -> genotype_str

    Returns:
        Dict with full longevity report data matching TypeScript LongevityAgingReport interface
    """

    # 1. APOE epsilon determination (special logic)
    apoe = _determine_apoe_epsilon(variants)

    # 2. Analyze all other sections
    foxo3_results = _analyze_snp_list(variants, FOXO3_SNPS)
    telomere_results = _analyze_snp_list(variants, TELOMERE_SNPS)
    cardio_results = _analyze_snp_list(variants, CARDIOVASCULAR_SNPS)
    metabolic_results = _analyze_snp_list(variants, METABOLIC_SNPS)
    inflammation_results = _analyze_snp_list(variants, INFLAMMATION_SNPS)
    oxidative_results = _analyze_snp_list(variants, OXIDATIVE_SNPS)
    cellular_results = _analyze_snp_list(variants, CELLULAR_REPAIR_SNPS)

    # 3. Section scores
    sections_raw = {
        "telomere": telomere_results,
        "inflammation": inflammation_results,
        "metabolic": metabolic_results,
        "oxidative": oxidative_results,
        "cellular": cellular_results,
    }

    apoe_score = apoe["score"]
    foxo3_score, foxo3_max = _section_score(foxo3_results)
    telomere_score, telomere_max = _section_score(telomere_results)
    cardio_score, cardio_max = _section_score(cardio_results)
    metabolic_score, metabolic_max = _section_score(metabolic_results)
    inflammation_score, inflammation_max = _section_score(inflammation_results)
    oxidative_score, oxidative_max = _section_score(oxidative_results)
    cellular_score, cellular_max = _section_score(cellular_results)

    # 4. Weighted overall longevity score (0-100)
    weighted = (
        apoe_score * SECTION_WEIGHTS["apoeBrainHeart"]
        + foxo3_score * SECTION_WEIGHTS["foxo3StressResistance"]
        + telomere_score * SECTION_WEIGHTS["telomereMaintenance"]
        + cardio_score * SECTION_WEIGHTS["cardiovascularLongevity"]
        + metabolic_score * SECTION_WEIGHTS["metabolicHealth"]
        + inflammation_score * SECTION_WEIGHTS["inflammationImmuneAging"]
        + oxidative_score * SECTION_WEIGHTS["oxidativeStressDefense"]
        + cellular_score * SECTION_WEIGHTS["cellularRepairMaintenance"]
    )
    longevity_score = round(weighted * 100)

    # 5. Count variants
    all_results = (
        foxo3_results + telomere_results + cardio_results
        + metabolic_results + inflammation_results
        + oxidative_results + cellular_results
    )
    total_analyzed = len(all_results) + 2  # +2 for APOE SNPs
    total_found = sum(1 for r in all_results if r["found"]) + (2 if apoe["found"] else 0)

    # 6. Key findings
    key_findings = []
    key_findings.append(f"APOE genotype: {apoe['genotype']} ({apoe['risk_level'].replace('_', ' ')})")

    for r in foxo3_results:
        if r["found"]:
            key_findings.append(
                f"FOXO3 longevity variant: {r['genotype']} — {r['genotypeData']['nickname']}"
            )

    if longevity_score >= 65:
        key_findings.append("Your overall genetic aging profile is above average")
    elif longevity_score <= 35:
        key_findings.append("Several aging-related risk factors were identified — see actionable insights")

    # 7. Build section objects
    def _build_section(
        title: str, description: str, results: List[Dict],
        score: float, max_score: float, warnings: List[str] = None
    ) -> Dict:
        return {
            "title": title,
            "description": description,
            "score": round(score * 100),
            "maxScore": round(max_score * 100),
            "label": _score_label(score),
            "interpretation": _build_section_interpretation(title, score, results),
            "variants": [
                {
                    "rsid": r["rsid"],
                    "gene": r["gene"],
                    "genotype": r["genotype"],
                    "trait": r["name"],
                    "effect": r.get("effectSize", "moderate"),
                    "interpretation": r["genotypeData"]["interpretation"],
                    "evidence": _evidence_tier(r["rsid"]),
                    "score": round((r["genotypeData"].get("score") or 0.5) * 100),
                    "found": r["found"],
                }
                for r in results
            ],
            **({"warnings": warnings} if warnings else {}),
        }

    apoe_section = {
        "title": "APOE & Brain/Heart Aging",
        "description": (
            "APOE is the strongest known genetic determinant of human longevity. "
            "Your epsilon genotype influences brain aging, cardiovascular health, "
            "and lipid metabolism."
        ),
        "score": round(apoe_score * 100),
        "maxScore": 100,
        "label": _score_label(apoe_score),
        "interpretation": apoe["interpretation"],
        "variants": [
            {
                "rsid": "rs429358 + rs7412",
                "gene": "APOE",
                "genotype": apoe["genotype"],
                "trait": "APOE Epsilon Haplotype",
                "effect": "large",
                "interpretation": apoe["interpretation"],
                "evidence": "Tier 1",
                "score": round(apoe_score * 100),
                "found": apoe["found"],
            }
        ],
        **({"warnings": apoe["warnings"]} if apoe.get("warnings") else {}),
    }

    # 8. Actionable insights
    insights = _generate_insights(apoe, sections_raw)

    # 9. Disclaimers
    disclaimers = [
        "This report analyzes genetic variants associated with aging pathways. Genetics accounts for approximately 10% of lifespan variation. Lifestyle, environment, healthcare access, and chance are far more influential.",
        "This is NOT a prediction of lifespan or health outcomes. It reflects genetic tendencies in biological pathways associated with aging research.",
        "APOE genotype information carries implications for Alzheimer's disease risk. If you are concerned about your APOE results, please consult a genetic counselor or healthcare provider.",
        "All effect sizes referenced are from population-level studies. Individual outcomes vary enormously and are not determined by genetics alone.",
        "This report is for educational and informational purposes only. It is not a medical diagnosis. Consult healthcare professionals for medical decisions.",
    ]

    return {
        "report_type": "longevity_aging",
        "version": "1.0",
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "summary": {
            "totalVariantsAnalyzed": total_analyzed,
            "totalVariantsFound": total_found,
            "longevityScore": longevity_score,
            "longevityScoreLabel": _overall_score_label(longevity_score),
            "keyFindings": key_findings,
            "apoeGenotype": apoe["genotype"],
        },
        "apoeBrainHeart": apoe_section,
        "foxo3StressResistance": _build_section(
            "FOXO3 & Stress Resistance",
            "FOXO3 is the second most replicated longevity gene. It enhances cellular stress resistance, DNA repair, and antioxidant defense.",
            foxo3_results, foxo3_score, foxo3_max,
        ),
        "telomereMaintenance": _build_section(
            "Telomere Maintenance",
            "Telomeres are protective caps on chromosome ends that shorten with each cell division. TERT and TERC maintain telomere length.",
            telomere_results, telomere_score, telomere_max,
        ),
        "cardiovascularLongevity": _build_section(
            "Cardiovascular Longevity",
            "Heart and vascular health is the strongest modifiable factor for lifespan. These variants influence cholesterol, vascular aging, and cardiac risk.",
            cardio_results, cardio_score, cardio_max,
        ),
        "metabolicHealth": _build_section(
            "Metabolic Health & Aging",
            "Metabolic health — insulin sensitivity, glucose regulation, and energy metabolism — is central to healthy aging.",
            metabolic_results, metabolic_score, metabolic_max,
        ),
        "inflammationImmuneAging": _build_section(
            "Inflammation & Immune Aging",
            "Chronic low-grade inflammation ('inflammaging') accelerates biological aging. These variants influence baseline inflammatory signaling.",
            inflammation_results, inflammation_score, inflammation_max,
        ),
        "oxidativeStressDefense": _build_section(
            "Oxidative Stress Defense",
            "Oxidative damage from free radicals contributes to cellular aging. These enzymes form your innate antioxidant defense system.",
            oxidative_results, oxidative_score, oxidative_max,
        ),
        "cellularRepairMaintenance": _build_section(
            "Cellular Repair & Maintenance",
            "DNA repair, epigenetic maintenance, and cellular quality control systems are essential for preventing age-related decline.",
            cellular_results, cellular_score, cellular_max,
        ),
        "actionableInsights": insights,
        "disclaimers": disclaimers,
    }


def _build_section_interpretation(title: str, score: float, results: List[Dict]) -> str:
    found = [r for r in results if r["found"]]
    if not found:
        return f"No variants were found for the {title} section in your DNA file."
    label = _score_label(score)
    return (
        f"Based on {len(found)} variant{'s' if len(found) > 1 else ''} analyzed, "
        f"your {title.lower()} profile is rated as '{label}' "
        f"(score: {round(score * 100)}/100)."
    )


# Evidence tier lookup
_TIER_MAP = {
    "rs429358": "Tier 1", "rs7412": "Tier 1",
    "rs2802292": "Tier 1",
    "rs2736100": "Tier 1", "rs10936599": "Tier 1",
    "rs1333049": "Tier 2", "rs6511720": "Tier 2", "rs646776": "Tier 2",
    "rs7903146": "Tier 2",
    "rs1800795": "Tier 2", "rs1205": "Tier 2",
    "rs4880": "Tier 3", "rs1050450": "Tier 3", "rs662": "Tier 3",
    "rs7895833": "Tier 3", "rs1801133": "Tier 2", "rs5882": "Tier 3",
    "rs1042522": "Tier 3", "rs9536314": "Tier 2",
}


def _evidence_tier(rsid: str) -> str:
    return _TIER_MAP.get(rsid.lower(), "Tier 3")


def generate_longevity_aging_json(result: Dict) -> str:
    """Convert the analysis result to JSON string."""
    return json.dumps(result, indent=2, ensure_ascii=False)
