"""
Comprehensive SNP Database
Curated database of clinically significant genetic variants
"""

from typing import Dict, List, Optional
from dataclasses import dataclass


@dataclass
class SNPInfo:
    rsid: str
    gene: str
    category: str
    risk_allele: str
    normal_allele: str
    significance: str  # high, moderate, low, beneficial
    condition: str
    description: str
    recommendations: List[str]


# Curated database of significant SNPs
SNP_DATABASE: Dict[str, SNPInfo] = {
    # Cardiovascular Health
    "rs1801133": SNPInfo(
        rsid="rs1801133",
        gene="MTHFR",
        category="Cardiovascular",
        risk_allele="T",
        normal_allele="C",
        significance="moderate",
        condition="MTHFR C677T Variant",
        description="Affects folate metabolism and homocysteine levels. Associated with increased cardiovascular disease risk.",
        recommendations=[
            "Consider methylfolate supplementation",
            "Monitor homocysteine levels",
            "Ensure adequate B12 intake"
        ]
    ),
    "rs429358": SNPInfo(
        rsid="rs429358",
        gene="APOE",
        category="Cardiovascular/Neurological",
        risk_allele="C",
        normal_allele="T",
        significance="high",
        condition="APOE4 Variant",
        description="Associated with increased risk of Alzheimer's disease and cardiovascular disease.",
        recommendations=[
            "Regular cardiovascular monitoring",
            "Consider cognitive assessments",
            "Heart-healthy diet",
            "Regular exercise"
        ]
    ),
    "rs7412": SNPInfo(
        rsid="rs7412",
        gene="APOE",
        category="Cardiovascular/Neurological",
        risk_allele="T",
        normal_allele="C",
        significance="beneficial",
        condition="APOE2 Variant",
        description="Associated with reduced cardiovascular risk and possible longevity benefits.",
        recommendations=[]
    ),
    "rs1800566": SNPInfo(
        rsid="rs1800566",
        gene="NQO1",
        category="Detoxification",
        risk_allele="T",
        normal_allele="C",
        significance="moderate",
        condition="NQO1*2 Variant",
        description="Reduced ability to detoxify certain compounds. May affect drug metabolism.",
        recommendations=[
            "Avoid excessive oxidative stress",
            "Consider antioxidant supplementation",
            "Discuss medication metabolism with healthcare provider"
        ]
    ),

    # Pharmacogenomics
    "rs1799853": SNPInfo(
        rsid="rs1799853",
        gene="CYP2C9",
        category="Drug Metabolism",
        risk_allele="T",
        normal_allele="C",
        significance="moderate",
        condition="CYP2C9*2 Variant",
        description="Reduced metabolism of warfarin, NSAIDs, and other medications.",
        recommendations=[
            "May need lower doses of warfarin",
            "Inform healthcare providers before surgery",
            "Caution with certain pain medications"
        ]
    ),
    "rs1057910": SNPInfo(
        rsid="rs1057910",
        gene="CYP2C9",
        category="Drug Metabolism",
        risk_allele="C",
        normal_allele="A",
        significance="moderate",
        condition="CYP2C9*3 Variant",
        description="Significantly reduced metabolism of many medications including warfarin.",
        recommendations=[
            "May need significantly lower doses of warfarin",
            "Inform all healthcare providers",
            "Consider pharmacogenomic testing before new medications"
        ]
    ),
    "rs4244285": SNPInfo(
        rsid="rs4244285",
        gene="CYP2C19",
        category="Drug Metabolism",
        risk_allele="A",
        normal_allele="G",
        significance="high",
        condition="CYP2C19*2 Variant",
        description="Poor metabolizer of clopidogrel (Plavix), PPIs, and some antidepressants.",
        recommendations=[
            "Alternative antiplatelet therapy may be needed",
            "Inform cardiologist before procedures",
            "May need alternative medications"
        ]
    ),
    "rs12248560": SNPInfo(
        rsid="rs12248560",
        gene="CYP2C19",
        category="Drug Metabolism",
        risk_allele="T",
        normal_allele="C",
        significance="moderate",
        condition="CYP2C19*17 Variant",
        description="Ultra-rapid metabolizer. May need higher doses of some medications.",
        recommendations=[
            "Standard doses may be less effective",
            "Discuss with healthcare provider",
            "Monitor medication effectiveness"
        ]
    ),
    "rs4986893": SNPInfo(
        rsid="rs4986893",
        gene="CYP2C19",
        category="Drug Metabolism",
        risk_allele="A",
        normal_allele="G",
        significance="high",
        condition="CYP2C19*3 Variant",
        description="Non-functional enzyme. Poor metabolizer of many medications.",
        recommendations=[
            "Significant drug metabolism impairment",
            "Alternative medications may be needed",
            "Comprehensive pharmacogenomic review recommended"
        ]
    ),

    # Metabolic Health
    "rs1801282": SNPInfo(
        rsid="rs1801282",
        gene="PPARG",
        category="Metabolic",
        risk_allele="G",
        normal_allele="C",
        significance="moderate",
        condition="PPARG Pro12Ala Variant",
        description="Associated with insulin sensitivity and type 2 diabetes risk.",
        recommendations=[
            "Monitor blood glucose regularly",
            "Maintain healthy weight",
            "Regular exercise",
            "Low glycemic diet"
        ]
    ),
    "rs7903146": SNPInfo(
        rsid="rs7903146",
        gene="TCF7L2",
        category="Metabolic",
        risk_allele="T",
        normal_allele="C",
        significance="high",
        condition="TCF7L2 Risk Variant",
        description="One of the strongest genetic risk factors for type 2 diabetes.",
        recommendations=[
            "Regular diabetes screening",
            "Strict glycemic control",
            "Lifestyle modifications critical",
            "Consider early intervention"
        ]
    ),
    "rs1800562": SNPInfo(
        rsid="rs1800562",
        gene="HFE",
        category="Metabolic",
        risk_allele="A",
        normal_allele="G",
        significance="high",
        condition="Hereditary Hemochromatosis (C282Y)",
        description="Risk for iron overload disorder. Two copies significantly increase risk.",
        recommendations=[
            "Regular iron and ferritin testing",
            "Avoid iron supplements unless deficient",
            "Limit vitamin C with meals",
            "Consider therapeutic phlebotomy if elevated"
        ]
    ),

    # Nutrition & Vitamins
    "rs12934922": SNPInfo(
        rsid="rs12934922",
        gene="BCMO1",
        category="Nutrition",
        risk_allele="T",
        normal_allele="A",
        significance="moderate",
        condition="Beta-Carotene Conversion Variant",
        description="Reduced ability to convert beta-carotene to vitamin A.",
        recommendations=[
            "May need preformed vitamin A (retinol)",
            "Include animal sources of vitamin A",
            "Consider retinol supplementation"
        ]
    ),
    "rs7041": SNPInfo(
        rsid="rs7041",
        gene="GC",
        category="Nutrition",
        risk_allele="T",
        normal_allele="G",
        significance="moderate",
        condition="Vitamin D Binding Protein Variant",
        description="May affect vitamin D transport and bioavailability.",
        recommendations=[
            "Regular vitamin D level testing",
            "May need higher supplementation",
            "Sun exposure for natural synthesis"
        ]
    ),
    "rs602662": SNPInfo(
        rsid="rs602662",
        gene="FUT2",
        category="Nutrition",
        risk_allele="A",
        normal_allele="G",
        significance="moderate",
        condition="Vitamin B12 Absorption Variant",
        description="May affect B12 absorption from food sources.",
        recommendations=[
            "Regular B12 level monitoring",
            "May benefit from sublingual B12",
            "Consider methylcobalamin form"
        ]
    ),

    # Inflammation & Immunity
    "rs1800795": SNPInfo(
        rsid="rs1800795",
        gene="IL6",
        category="Inflammation",
        risk_allele="C",
        normal_allele="G",
        significance="moderate",
        condition="IL-6 Promoter Variant",
        description="Associated with increased inflammatory response.",
        recommendations=[
            "Anti-inflammatory diet",
            "Omega-3 fatty acids",
            "Regular exercise",
            "Stress management"
        ]
    ),
    "rs1800896": SNPInfo(
        rsid="rs1800896",
        gene="IL10",
        category="Inflammation",
        risk_allele="A",
        normal_allele="G",
        significance="moderate",
        condition="IL-10 Variant",
        description="Affects anti-inflammatory cytokine production.",
        recommendations=[
            "Support immune balance",
            "Anti-inflammatory lifestyle",
            "Consider probiotics"
        ]
    ),

    # Detoxification
    "rs1695": SNPInfo(
        rsid="rs1695",
        gene="GSTP1",
        category="Detoxification",
        risk_allele="G",
        normal_allele="A",
        significance="moderate",
        condition="GSTP1 Variant",
        description="Reduced glutathione S-transferase activity. Affects detoxification.",
        recommendations=[
            "Support glutathione production",
            "Cruciferous vegetables",
            "NAC or glutathione supplementation",
            "Limit toxin exposure"
        ]
    ),
    "rs4680": SNPInfo(
        rsid="rs4680",
        gene="COMT",
        category="Neurological",
        risk_allele="A",
        normal_allele="G",
        significance="moderate",
        condition="COMT Val158Met Variant",
        description="Affects dopamine and estrogen metabolism. The 'worrier' vs 'warrior' gene.",
        recommendations=[
            "SAMe may help if slow COMT",
            "Stress management techniques",
            "Avoid excessive caffeine if slow COMT"
        ]
    ),

    # Sleep & Circadian
    "rs73598374": SNPInfo(
        rsid="rs73598374",
        gene="ADA",
        category="Sleep",
        risk_allele="T",
        normal_allele="C",
        significance="low",
        condition="Adenosine Deaminase Variant",
        description="Associated with deep sleep architecture.",
        recommendations=[
            "Maintain consistent sleep schedule",
            "Optimize sleep environment"
        ]
    ),

    # Caffeine Metabolism
    "rs762551": SNPInfo(
        rsid="rs762551",
        gene="CYP1A2",
        category="Drug Metabolism",
        risk_allele="C",
        normal_allele="A",
        significance="low",
        condition="Caffeine Metabolism Variant",
        description="Slow caffeine metabolizer. Caffeine has longer effects.",
        recommendations=[
            "Limit afternoon caffeine",
            "May be more sensitive to caffeine effects",
            "Consider cutting caffeine earlier in day"
        ]
    ),

    # Celiac & Gluten
    "rs2187668": SNPInfo(
        rsid="rs2187668",
        gene="HLA-DQ2",
        category="Autoimmune",
        risk_allele="T",
        normal_allele="C",
        significance="moderate",
        condition="Celiac Disease Risk (HLA-DQ2)",
        description="Genetic predisposition to celiac disease.",
        recommendations=[
            "Consider celiac antibody testing if symptoms present",
            "Monitor for gluten sensitivity symptoms",
            "Genetic risk does not mean disease will develop"
        ]
    ),

    # Cancer Risk
    "rs1042522": SNPInfo(
        rsid="rs1042522",
        gene="TP53",
        category="Cancer",
        risk_allele="C",
        normal_allele="G",
        significance="moderate",
        condition="p53 Codon 72 Variant",
        description="May affect p53 tumor suppressor function.",
        recommendations=[
            "Regular cancer screenings",
            "Healthy lifestyle",
            "Avoid known carcinogens"
        ]
    ),

    # Cardiovascular - Additional
    "rs1333049": SNPInfo(
        rsid="rs1333049",
        gene="9p21",
        category="Cardiovascular",
        risk_allele="C",
        normal_allele="G",
        significance="high",
        condition="9p21 Coronary Artery Disease Risk",
        description="One of the strongest genetic markers for coronary artery disease.",
        recommendations=[
            "Aggressive cardiovascular risk management",
            "Regular cardiac assessments",
            "Strict blood pressure control",
            "Lipid management"
        ]
    ),
    "rs10757274": SNPInfo(
        rsid="rs10757274",
        gene="9p21",
        category="Cardiovascular",
        risk_allele="G",
        normal_allele="A",
        significance="high",
        condition="9p21 Myocardial Infarction Risk",
        description="Associated with increased risk of heart attack.",
        recommendations=[
            "Heart-healthy lifestyle critical",
            "Regular cardiac monitoring",
            "Know heart attack warning signs"
        ]
    ),

    # Clotting Factors
    "rs6025": SNPInfo(
        rsid="rs6025",
        gene="F5",
        category="Blood Clotting",
        risk_allele="T",
        normal_allele="C",
        significance="high",
        condition="Factor V Leiden",
        description="Increased risk of blood clots (deep vein thrombosis, pulmonary embolism).",
        recommendations=[
            "Avoid prolonged immobility",
            "Discuss with doctor before surgery or travel",
            "Caution with estrogen-containing medications",
            "Stay hydrated"
        ]
    ),
    "rs1799963": SNPInfo(
        rsid="rs1799963",
        gene="F2",
        category="Blood Clotting",
        risk_allele="A",
        normal_allele="G",
        significance="high",
        condition="Prothrombin G20210A",
        description="Increased risk of venous thromboembolism.",
        recommendations=[
            "Avoid prolonged immobility",
            "Inform healthcare providers",
            "Caution with hormonal therapies",
            "Compression stockings for long travel"
        ]
    ),
}


def get_snp_info(rsid: str) -> Optional[SNPInfo]:
    """Get SNP information by rsID"""
    return SNP_DATABASE.get(rsid.lower())


def get_snps_by_category(category: str) -> List[SNPInfo]:
    """Get all SNPs in a category"""
    return [snp for snp in SNP_DATABASE.values() if snp.category == category]


def get_all_categories() -> List[str]:
    """Get list of all categories"""
    return list(set(snp.category for snp in SNP_DATABASE.values()))


def get_high_significance_snps() -> List[SNPInfo]:
    """Get all high significance SNPs"""
    return [snp for snp in SNP_DATABASE.values() if snp.significance == "high"]
