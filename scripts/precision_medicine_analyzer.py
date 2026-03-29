"""
Precision Medicine / Pharmacogenomics Analyzer
Analyzes genetic variants related to drug metabolism and response.

Covers 30+ pharmacogenomic SNPs across 15 genes and assesses 50+ medications
organized by therapeutic category (cardiovascular, pain, psychiatry, oncology,
infectious disease, gastroenterology, hormonal, transplant, other).

Evidence levels follow CPIC (Clinical Pharmacogenetics Implementation
Consortium) and PharmGKB classification:
  1A = Strong evidence, CPIC guideline available
  1B = Strong evidence, moderate clinical context
  2A = Moderate evidence, known pharmacogene
  2B = Moderate evidence, limited clinical context

References:
- Caudle et al. (2014) Clin Pharmacol Ther - CPIC guidelines
- Relling & Klein (2011) Clin Pharmacol Ther - CPIC overview
- Whirl-Carrillo et al. (2012) Clin Pharmacol Ther - PharmGKB
- Scott et al. (2013) Clin Pharmacol Ther - CYP2C19 & clopidogrel
- Crews et al. (2014) Clin Pharmacol Ther - CYP2D6 & codeine
- Johnson et al. (2017) Clin Pharmacol Ther - CYP2C9/VKORC1 & warfarin
- Amstutz et al. (2018) Clin Pharmacol Ther - DPYD & fluoropyrimidines
- Relling et al. (2019) Clin Pharmacol Ther - TPMT/NUDT15 & thiopurines
"""

import json
from datetime import datetime, timezone
from typing import Any, Dict, List, Optional, Tuple


# ─────────────────────────────────────────────────────────────────────────────
# COMPLEMENT MAP (for strand-flip handling)
# ─────────────────────────────────────────────────────────────────────────────

COMPLEMENT = {"A": "T", "T": "A", "C": "G", "G": "C"}


def _complement_allele(allele: str) -> str:
    """Return the complement of a single nucleotide."""
    return COMPLEMENT.get(allele.upper(), allele.upper())


def _count_allele(genotype: str, allele: str) -> int:
    """
    Count occurrences of an allele in a genotype, handling complement strands.
    For ambiguous SNPs (A/T or C/G), checks both the allele and its complement.
    """
    genotype = genotype.upper().replace("-", "")
    allele = allele.upper()
    comp = _complement_allele(allele)

    count = genotype.count(allele)

    # Only count complement if it differs from the original allele and the
    # genotype doesn't already contain the original.
    if count == 0 and comp != allele:
        count = genotype.count(comp)

    return min(count, 2)


def _lookup_snp(
    variants: Dict[str, Tuple[str, str, str]], rsid: str
) -> Optional[str]:
    """
    Look up a genotype from the variants dict.
    Handles both tuple format (chrom, pos, genotype) and plain string genotype.

    Returns the genotype string (uppercase) or None if not found.
    """
    rsid_lower = rsid.lower()
    if rsid_lower not in variants:
        return None

    val = variants[rsid_lower]
    if isinstance(val, (list, tuple)) and len(val) >= 3:
        genotype = str(val[2])
    elif isinstance(val, str):
        genotype = val
    else:
        return None

    genotype = genotype.upper().replace("-", "").strip()
    if not genotype or genotype in ("--", "NA", "N/A", "00", "DD", "II"):
        return None
    return genotype


# ─────────────────────────────────────────────────────────────────────────────
# PHARMACOGENOMIC SNP DATABASE
# ─────────────────────────────────────────────────────────────────────────────

PHARMA_DATABASE: Dict[str, Dict[str, str]] = {
    # CYP2C19 - Clopidogrel, SSRIs, PPIs
    "rs4244285": {
        "gene": "CYP2C19", "star": "*2", "name": "CYP2C19*2",
        "function": "Loss of function", "risk_allele": "A",
        "evidence": "1A", "source": "CPIC",
    },
    "rs4986893": {
        "gene": "CYP2C19", "star": "*3", "name": "CYP2C19*3",
        "function": "Loss of function", "risk_allele": "A",
        "evidence": "1A", "source": "CPIC",
    },
    "rs12248560": {
        "gene": "CYP2C19", "star": "*17", "name": "CYP2C19*17",
        "function": "Increased function", "risk_allele": "T",
        "evidence": "1A", "source": "CPIC",
    },
    # CYP2D6 - Codeine, Tamoxifen, TCAs, SSRIs
    "rs3892097": {
        "gene": "CYP2D6", "star": "*4", "name": "CYP2D6*4",
        "function": "Loss of function", "risk_allele": "A",
        "evidence": "1A", "source": "CPIC",
    },
    "rs5030655": {
        "gene": "CYP2D6", "star": "*6", "name": "CYP2D6*6",
        "function": "Loss of function", "risk_allele": "DEL",
        "evidence": "1A", "source": "CPIC",
    },
    "rs16947": {
        "gene": "CYP2D6", "star": "*2", "name": "CYP2D6*2",
        "function": "Normal function", "risk_allele": "A",
        "evidence": "1A", "source": "CPIC",
    },
    "rs1065852": {
        "gene": "CYP2D6", "star": "*10", "name": "CYP2D6*10",
        "function": "Decreased function", "risk_allele": "T",
        "evidence": "1A", "source": "CPIC",
    },
    # CYP2C9 - Warfarin, NSAIDs
    "rs1799853": {
        "gene": "CYP2C9", "star": "*2", "name": "CYP2C9*2",
        "function": "Decreased function", "risk_allele": "T",
        "evidence": "1A", "source": "CPIC",
    },
    "rs1057910": {
        "gene": "CYP2C9", "star": "*3", "name": "CYP2C9*3",
        "function": "Loss of function", "risk_allele": "C",
        "evidence": "1A", "source": "CPIC",
    },
    # CYP1A2 - Caffeine, Clozapine
    "rs762551": {
        "gene": "CYP1A2", "star": "*1F", "name": "CYP1A2*1F",
        "function": "Inducible / ultra-rapid", "risk_allele": "A",
        "evidence": "2A", "source": "CPIC",
    },
    # CYP3A5 - Tacrolimus, many drugs
    "rs776746": {
        "gene": "CYP3A5", "star": "*3", "name": "CYP3A5*3",
        "function": "Loss of function (non-expressor)", "risk_allele": "G",
        "evidence": "1A", "source": "CPIC",
    },
    # VKORC1 - Warfarin
    "rs9923231": {
        "gene": "VKORC1", "star": "-1639G>A", "name": "VKORC1 -1639G>A",
        "function": "Increased warfarin sensitivity", "risk_allele": "A",
        "evidence": "1A", "source": "CPIC",
    },
    # DPYD - 5-FU, Capecitabine
    "rs3918290": {
        "gene": "DPYD", "star": "*2A", "name": "DPYD*2A",
        "function": "Loss of function", "risk_allele": "A",
        "evidence": "1A", "source": "CPIC",
    },
    "rs55886062": {
        "gene": "DPYD", "star": "*13", "name": "DPYD*13",
        "function": "Loss of function", "risk_allele": "A",
        "evidence": "1A", "source": "CPIC",
    },
    # TPMT - Azathioprine, 6-MP
    "rs1800462": {
        "gene": "TPMT", "star": "*2", "name": "TPMT*2",
        "function": "Loss of function", "risk_allele": "G",
        "evidence": "1A", "source": "CPIC",
    },
    "rs1800460": {
        "gene": "TPMT", "star": "*3B", "name": "TPMT*3B",
        "function": "Loss of function", "risk_allele": "G",
        "evidence": "1A", "source": "CPIC",
    },
    "rs1142345": {
        "gene": "TPMT", "star": "*3C", "name": "TPMT*3C",
        "function": "Loss of function", "risk_allele": "G",
        "evidence": "1A", "source": "CPIC",
    },
    # UGT1A1 - Irinotecan, Atazanavir
    "rs8175347": {
        "gene": "UGT1A1", "star": "*28", "name": "UGT1A1*28",
        "function": "Decreased function (7 TA repeats)", "risk_allele": "7",
        "evidence": "1A", "source": "CPIC",
    },
    # SLCO1B1 - Statins
    "rs4149056": {
        "gene": "SLCO1B1", "star": "*5", "name": "SLCO1B1*5",
        "function": "Decreased transporter function", "risk_allele": "C",
        "evidence": "1A", "source": "CPIC",
    },
    # ABCB1 (P-glycoprotein) - many drugs
    "rs1045642": {
        "gene": "ABCB1", "star": "3435C>T", "name": "ABCB1 3435C>T",
        "function": "Altered P-glycoprotein expression", "risk_allele": "T",
        "evidence": "2A", "source": "PharmGKB",
    },
    # HLA-B*57:01 - Abacavir
    "rs2395029": {
        "gene": "HLA-B", "star": "*57:01 tag", "name": "HLA-B*57:01",
        "function": "Abacavir hypersensitivity", "risk_allele": "G",
        "evidence": "1A", "source": "CPIC",
    },
    # OPRM1 - Opioid response
    "rs1799971": {
        "gene": "OPRM1", "star": "A118G", "name": "OPRM1 A118G",
        "function": "Altered opioid receptor binding", "risk_allele": "G",
        "evidence": "2A", "source": "PharmGKB",
    },
    # COMT - Pain sensitivity, catecholamine metabolism
    "rs4680": {
        "gene": "COMT", "star": "Val158Met", "name": "COMT Val158Met",
        "function": "Altered catechol-O-methyltransferase activity",
        "risk_allele": "A",
        "evidence": "2A", "source": "PharmGKB",
    },
    # MTHFR - Methotrexate, folate
    "rs1801133": {
        "gene": "MTHFR", "star": "C677T", "name": "MTHFR C677T",
        "function": "Decreased enzyme activity", "risk_allele": "T",
        "evidence": "2A", "source": "PharmGKB",
    },
    "rs1801131": {
        "gene": "MTHFR", "star": "A1298C", "name": "MTHFR A1298C",
        "function": "Decreased enzyme activity", "risk_allele": "G",
        "evidence": "2B", "source": "PharmGKB",
    },
    # Factor V Leiden - Thrombosis risk with estrogen
    "rs6025": {
        "gene": "F5", "star": "Leiden", "name": "Factor V Leiden",
        "function": "Activated protein C resistance", "risk_allele": "T",
        "evidence": "1A", "source": "CPIC",
    },
    # IFNL3/IL28B - Hepatitis C treatment
    "rs12979860": {
        "gene": "IFNL3", "star": "IL28B", "name": "IFNL3/IL28B",
        "function": "Hepatitis C treatment response", "risk_allele": "T",
        "evidence": "1A", "source": "CPIC",
    },
    # NAT2 - Isoniazid
    "rs1801280": {
        "gene": "NAT2", "star": "*5", "name": "NAT2*5",
        "function": "Slow acetylator", "risk_allele": "C",
        "evidence": "2A", "source": "PharmGKB",
    },
    # IFNL4 - Hepatitis C
    "rs368234815": {
        "gene": "IFNL4", "star": "ss469415590", "name": "IFNL4 \u0394G",
        "function": "Interferon lambda 4 production", "risk_allele": "TG",
        "evidence": "2A", "source": "PharmGKB",
    },
    # CYP2B6 - Efavirenz
    "rs3745274": {
        "gene": "CYP2B6", "star": "*6", "name": "CYP2B6*6",
        "function": "Decreased function", "risk_allele": "T",
        "evidence": "1A", "source": "CPIC",
    },
    # NUDT15 - Thiopurines
    "rs116855232": {
        "gene": "NUDT15", "star": "*3", "name": "NUDT15*3",
        "function": "Loss of function", "risk_allele": "T",
        "evidence": "1A", "source": "CPIC",
    },
}


# ─────────────────────────────────────────────────────────────────────────────
# GENE METADATA
# ─────────────────────────────────────────────────────────────────────────────

GENE_FULL_NAMES: Dict[str, str] = {
    "CYP2C19": "Cytochrome P450 2C19",
    "CYP2D6": "Cytochrome P450 2D6",
    "CYP2C9": "Cytochrome P450 2C9",
    "CYP1A2": "Cytochrome P450 1A2",
    "CYP3A5": "Cytochrome P450 3A5",
    "CYP2B6": "Cytochrome P450 2B6",
    "VKORC1": "Vitamin K Epoxide Reductase Complex Subunit 1",
    "DPYD": "Dihydropyrimidine Dehydrogenase",
    "TPMT": "Thiopurine S-Methyltransferase",
    "NUDT15": "Nudix Hydrolase 15",
    "UGT1A1": "UDP Glucuronosyltransferase Family 1 Member A1",
    "SLCO1B1": "Solute Carrier Organic Anion Transporter 1B1",
    "ABCB1": "ATP Binding Cassette Subfamily B Member 1 (P-glycoprotein)",
    "HLA-B": "Human Leukocyte Antigen B",
    "OPRM1": "Opioid Receptor Mu 1",
    "COMT": "Catechol-O-Methyltransferase",
    "MTHFR": "Methylenetetrahydrofolate Reductase",
    "F5": "Coagulation Factor V",
    "IFNL3": "Interferon Lambda 3",
    "IFNL4": "Interferon Lambda 4",
    "NAT2": "N-Acetyltransferase 2",
}

# SNPs belonging to each gene for classification
GENE_SNPS: Dict[str, List[str]] = {
    "CYP2C19": ["rs4244285", "rs4986893", "rs12248560"],
    "CYP2D6": ["rs3892097", "rs5030655", "rs16947", "rs1065852"],
    "CYP2C9": ["rs1799853", "rs1057910"],
    "CYP1A2": ["rs762551"],
    "CYP3A5": ["rs776746"],
    "CYP2B6": ["rs3745274"],
    "VKORC1": ["rs9923231"],
    "DPYD": ["rs3918290", "rs55886062"],
    "TPMT": ["rs1800462", "rs1800460", "rs1142345"],
    "NUDT15": ["rs116855232"],
    "UGT1A1": ["rs8175347"],
    "SLCO1B1": ["rs4149056"],
    "ABCB1": ["rs1045642"],
    "HLA-B": ["rs2395029"],
    "OPRM1": ["rs1799971"],
    "COMT": ["rs4680"],
    "MTHFR": ["rs1801133", "rs1801131"],
    "F5": ["rs6025"],
    "IFNL3": ["rs12979860"],
    "IFNL4": ["rs368234815"],
    "NAT2": ["rs1801280"],
}

# Population frequency notes per gene/status
POPULATION_FREQUENCY: Dict[str, str] = {
    "CYP2C19": "About 2-5% of Caucasians and 12-23% of Asians are CYP2C19 poor metabolizers",
    "CYP2D6": "About 5-10% of Caucasians are CYP2D6 poor metabolizers; 1-2% are ultra-rapid",
    "CYP2C9": "About 1-3% of Caucasians carry CYP2C9*3/*3 (poor metabolizer) genotype",
    "CYP1A2": "About 45-55% of people carry the CYP1A2*1F allele (rapid inducible metabolism)",
    "CYP3A5": "About 80-90% of Caucasians are CYP3A5 non-expressors (*3/*3)",
    "CYP2B6": "About 3-8% of people are CYP2B6 poor metabolizers",
    "VKORC1": "About 37% of Caucasians carry at least one A allele at VKORC1 -1639",
    "DPYD": "About 3-5% of people carry at least one DPYD reduced-function allele",
    "TPMT": "About 10% of people are TPMT intermediate metabolizers; 0.3% are poor",
    "NUDT15": "About 2% of Asians carry NUDT15*3; rare in Caucasians",
    "UGT1A1": "About 10-16% of Caucasians are homozygous for UGT1A1*28",
    "SLCO1B1": "About 15-20% of people carry at least one SLCO1B1*5 allele",
    "ABCB1": "The ABCB1 3435T allele frequency is approximately 50% in Caucasians",
    "HLA-B": "HLA-B*57:01 is found in about 5-8% of Caucasians",
    "OPRM1": "About 10-15% of Caucasians carry the OPRM1 118G allele",
    "COMT": "COMT Val158Met heterozygosity is found in about 50% of people",
    "MTHFR": "About 10-15% of Caucasians are homozygous for MTHFR C677T",
    "F5": "Factor V Leiden is found in about 3-8% of Caucasians",
    "IFNL3": "IFNL3 CC genotype (favorable for HCV treatment) is found in about 40% of Caucasians",
    "IFNL4": "IFNL4 TT/TT genotype is found in about 70% of Caucasians",
    "NAT2": "About 40-60% of Caucasians are slow acetylators (NAT2)",
}


# ─────────────────────────────────────────────────────────────────────────────
# DRUG DATABASE (50+ medications)
# ─────────────────────────────────────────────────────────────────────────────

DRUG_DATABASE: Dict[str, Dict[str, Any]] = {
    # ── CARDIOVASCULAR ──
    "warfarin": {
        "name": "Warfarin",
        "brand_names": ["Coumadin", "Jantoven"],
        "category": "cardiovascular",
        "genes": ["CYP2C9", "VKORC1"],
        "snps": ["rs1799853", "rs1057910", "rs9923231"],
        "indication": "Blood thinner for stroke/DVT prevention",
        "risk_level": "high",
        "cpic_level": "1A",
        "fda_label": True,
        "pubmed_ids": ["28198005", "25974703"],
        "interpretations": {
            "normal": "Standard dosing expected to be appropriate. Initiate per clinical protocol.",
            "intermediate": "Consider 25-50% dose reduction. Monitor INR closely during initiation.",
            "poor": "Consider 50-75% dose reduction. High bleeding risk. Monitor INR very closely.",
            "rapid": "May need higher doses to achieve therapeutic INR. Monitor closely.",
        },
        "alternatives": ["Apixaban", "Rivaroxaban", "Dabigatran (DOACs not affected by CYP2C9/VKORC1)"],
    },
    "clopidogrel": {
        "name": "Clopidogrel",
        "brand_names": ["Plavix"],
        "category": "cardiovascular",
        "genes": ["CYP2C19"],
        "snps": ["rs4244285", "rs4986893", "rs12248560"],
        "indication": "Antiplatelet for acute coronary syndrome and stent placement",
        "risk_level": "high",
        "cpic_level": "1A",
        "fda_label": True,
        "pubmed_ids": ["23698643", "22992668"],
        "interpretations": {
            "normal": "Standard clopidogrel dosing expected to be effective.",
            "intermediate": "Reduced platelet inhibition. Consider alternative antiplatelet (prasugrel or ticagrelor) if undergoing PCI.",
            "poor": "Significantly reduced activation of clopidogrel. Use alternative antiplatelet therapy (prasugrel or ticagrelor). FDA boxed warning.",
            "rapid": "Enhanced platelet inhibition. Standard dosing appropriate; may have increased bleeding risk.",
        },
        "alternatives": ["Prasugrel", "Ticagrelor"],
    },
    "simvastatin": {
        "name": "Simvastatin",
        "brand_names": ["Zocor"],
        "category": "cardiovascular",
        "genes": ["SLCO1B1"],
        "snps": ["rs4149056"],
        "indication": "Cholesterol-lowering statin for cardiovascular risk reduction",
        "risk_level": "high",
        "cpic_level": "1A",
        "fda_label": True,
        "pubmed_ids": ["22617227", "24918167"],
        "interpretations": {
            "normal": "Standard dosing appropriate. Low risk of statin-induced myopathy.",
            "intermediate": "Increased myopathy risk. Avoid simvastatin >20 mg/day. Consider alternative statin.",
            "poor": "High myopathy risk. Avoid simvastatin or use lowest dose. Strongly consider alternative statin (rosuvastatin, pravastatin).",
            "rapid": "Standard dosing appropriate.",
        },
        "alternatives": ["Rosuvastatin", "Pravastatin", "Fluvastatin"],
    },
    "atorvastatin": {
        "name": "Atorvastatin",
        "brand_names": ["Lipitor"],
        "category": "cardiovascular",
        "genes": ["SLCO1B1"],
        "snps": ["rs4149056"],
        "indication": "Cholesterol-lowering statin for cardiovascular risk reduction",
        "risk_level": "moderate",
        "cpic_level": "1A",
        "fda_label": False,
        "pubmed_ids": ["22617227"],
        "interpretations": {
            "normal": "Standard dosing appropriate.",
            "intermediate": "Moderate myopathy risk. Monitor for muscle symptoms. Consider dose adjustment.",
            "poor": "Elevated myopathy risk. Use lower dose or consider rosuvastatin/pravastatin.",
            "rapid": "Standard dosing appropriate.",
        },
        "alternatives": ["Rosuvastatin", "Pravastatin"],
    },
    "metoprolol": {
        "name": "Metoprolol",
        "brand_names": ["Lopressor", "Toprol-XL"],
        "category": "cardiovascular",
        "genes": ["CYP2D6"],
        "snps": ["rs3892097", "rs5030655", "rs16947", "rs1065852"],
        "indication": "Beta-blocker for hypertension, heart failure, and arrhythmias",
        "risk_level": "moderate",
        "cpic_level": "1A",
        "fda_label": False,
        "pubmed_ids": ["30801677"],
        "interpretations": {
            "normal": "Standard dosing expected to be appropriate.",
            "intermediate": "Slightly increased drug exposure. Consider standard dose with monitoring.",
            "poor": "Significantly increased drug levels. Consider 50% dose reduction or use alternative beta-blocker (atenolol, bisoprolol).",
            "rapid": "Reduced drug exposure. May need higher doses or alternative agent.",
        },
        "alternatives": ["Atenolol", "Bisoprolol", "Carvedilol"],
    },
    # ── PAIN MANAGEMENT ──
    "codeine": {
        "name": "Codeine",
        "brand_names": ["Tylenol #3", "Various"],
        "category": "pain",
        "genes": ["CYP2D6"],
        "snps": ["rs3892097", "rs5030655", "rs16947", "rs1065852"],
        "indication": "Mild to moderate pain relief",
        "risk_level": "high",
        "cpic_level": "1A",
        "fda_label": True,
        "pubmed_ids": ["25974703", "31006110"],
        "interpretations": {
            "normal": "Standard codeine metabolism. Expected normal analgesic effect.",
            "intermediate": "Reduced conversion to morphine. May experience decreased pain relief. Consider alternative analgesics.",
            "poor": "Minimal conversion to morphine. Codeine will be ineffective for pain. Use alternative analgesics.",
            "rapid": "CAUTION: Excessive conversion to morphine. Risk of respiratory depression and toxicity, especially in children. AVOID codeine.",
        },
        "alternatives": ["Acetaminophen", "NSAIDs (check CYP2C9)", "Non-opioid analgesics"],
    },
    "tramadol": {
        "name": "Tramadol",
        "brand_names": ["Ultram", "ConZip"],
        "category": "pain",
        "genes": ["CYP2D6"],
        "snps": ["rs3892097", "rs5030655", "rs16947", "rs1065852"],
        "indication": "Moderate to moderately severe pain",
        "risk_level": "high",
        "cpic_level": "1A",
        "fda_label": True,
        "pubmed_ids": ["31006110"],
        "interpretations": {
            "normal": "Standard tramadol metabolism expected.",
            "intermediate": "Reduced formation of active metabolite (O-desmethyltramadol). Decreased efficacy possible.",
            "poor": "Greatly reduced formation of active metabolite. Tramadol likely ineffective. Use alternative analgesic.",
            "rapid": "CAUTION: Increased formation of active metabolite. Risk of toxicity. Consider alternative or reduced dose.",
        },
        "alternatives": ["Acetaminophen", "Non-opioid analgesics"],
    },
    "ibuprofen": {
        "name": "Ibuprofen",
        "brand_names": ["Advil", "Motrin"],
        "category": "pain",
        "genes": ["CYP2C9"],
        "snps": ["rs1799853", "rs1057910"],
        "indication": "Pain, inflammation, and fever reduction",
        "risk_level": "moderate",
        "cpic_level": "2A",
        "fda_label": False,
        "pubmed_ids": ["32189324"],
        "interpretations": {
            "normal": "Standard NSAID metabolism. Standard dosing appropriate.",
            "intermediate": "Slightly increased drug exposure. Use lowest effective dose. Monitor for GI bleeding.",
            "poor": "Significantly increased drug exposure and half-life. Increased risk of GI bleeding. Use lowest dose for shortest duration or consider alternative.",
            "rapid": "Standard dosing appropriate.",
        },
        "alternatives": ["Acetaminophen", "Celecoxib (lower GI risk)"],
    },
    "morphine": {
        "name": "Morphine / Opioids (general)",
        "brand_names": ["MS Contin", "Various"],
        "category": "pain",
        "genes": ["OPRM1", "COMT"],
        "snps": ["rs1799971", "rs4680"],
        "indication": "Moderate to severe pain management",
        "risk_level": "moderate",
        "cpic_level": "2A",
        "fda_label": False,
        "pubmed_ids": ["24458010"],
        "interpretations": {
            "normal": "Standard opioid receptor function and pain sensitivity.",
            "intermediate": "May have altered opioid response. Dose titration may need adjustment based on clinical response.",
            "poor": "Reduced opioid receptor binding affinity. May require higher doses for pain relief. Monitor closely.",
            "rapid": "Standard dosing appropriate with monitoring.",
        },
        "alternatives": ["Non-opioid analgesics", "Multimodal pain management"],
    },
    # ── PSYCHIATRY / MENTAL HEALTH ──
    "citalopram": {
        "name": "Citalopram",
        "brand_names": ["Celexa"],
        "category": "psychiatry",
        "genes": ["CYP2C19"],
        "snps": ["rs4244285", "rs4986893", "rs12248560"],
        "indication": "Depression (SSRI antidepressant)",
        "risk_level": "high",
        "cpic_level": "1A",
        "fda_label": True,
        "pubmed_ids": ["25974703"],
        "interpretations": {
            "normal": "Standard dosing expected to be effective.",
            "intermediate": "Increased drug levels possible. Standard starting dose, but monitor for side effects.",
            "poor": "Significantly increased drug levels. FDA recommends maximum dose of 20 mg/day. Monitor QTc interval.",
            "rapid": "Reduced drug exposure. May need dose increase if inadequate response. Consider alternative SSRI.",
        },
        "alternatives": ["Sertraline", "Fluoxetine", "Bupropion"],
    },
    "escitalopram": {
        "name": "Escitalopram",
        "brand_names": ["Lexapro"],
        "category": "psychiatry",
        "genes": ["CYP2C19"],
        "snps": ["rs4244285", "rs4986893", "rs12248560"],
        "indication": "Depression and generalized anxiety (SSRI)",
        "risk_level": "high",
        "cpic_level": "1A",
        "fda_label": True,
        "pubmed_ids": ["25974703"],
        "interpretations": {
            "normal": "Standard dosing expected to be effective.",
            "intermediate": "Increased drug levels possible. Standard starting dose with monitoring.",
            "poor": "Significantly increased drug levels. Consider 50% dose reduction. Monitor QTc interval.",
            "rapid": "Reduced drug exposure. May need dose increase. Consider alternative SSRI.",
        },
        "alternatives": ["Sertraline", "Fluoxetine", "Bupropion"],
    },
    "sertraline": {
        "name": "Sertraline",
        "brand_names": ["Zoloft"],
        "category": "psychiatry",
        "genes": ["CYP2C19"],
        "snps": ["rs4244285", "rs4986893", "rs12248560"],
        "indication": "Depression, anxiety, OCD, PTSD (SSRI)",
        "risk_level": "moderate",
        "cpic_level": "1A",
        "fda_label": False,
        "pubmed_ids": ["25974703"],
        "interpretations": {
            "normal": "Standard dosing expected to be effective.",
            "intermediate": "Minor increase in drug levels. Standard dosing with monitoring.",
            "poor": "Increased drug levels. Consider lower starting dose or alternative SSRI.",
            "rapid": "Reduced drug exposure. Consider dose titration upward if inadequate response.",
        },
        "alternatives": ["Fluoxetine", "Bupropion", "Venlafaxine"],
    },
    "fluoxetine": {
        "name": "Fluoxetine",
        "brand_names": ["Prozac"],
        "category": "psychiatry",
        "genes": ["CYP2D6"],
        "snps": ["rs3892097", "rs5030655", "rs16947", "rs1065852"],
        "indication": "Depression, OCD, bulimia, panic disorder (SSRI)",
        "risk_level": "moderate",
        "cpic_level": "1A",
        "fda_label": False,
        "pubmed_ids": ["25974703"],
        "interpretations": {
            "normal": "Standard dosing expected to be appropriate.",
            "intermediate": "Slightly increased drug levels. Standard dosing with monitoring for side effects.",
            "poor": "Increased drug levels and prolonged half-life. Consider lower starting dose or alternative.",
            "rapid": "Reduced drug levels. May need dose increase if inadequate response.",
        },
        "alternatives": ["Sertraline", "Citalopram", "Bupropion"],
    },
    "paroxetine": {
        "name": "Paroxetine",
        "brand_names": ["Paxil", "Seroxat"],
        "category": "psychiatry",
        "genes": ["CYP2D6"],
        "snps": ["rs3892097", "rs5030655", "rs16947", "rs1065852"],
        "indication": "Depression, anxiety, OCD, panic disorder (SSRI)",
        "risk_level": "moderate",
        "cpic_level": "1A",
        "fda_label": False,
        "pubmed_ids": ["25974703"],
        "interpretations": {
            "normal": "Standard dosing expected to be appropriate.",
            "intermediate": "Slightly increased drug levels. Standard starting dose appropriate.",
            "poor": "Increased drug levels. Consider 50% dose reduction or select alternative SSRI not primarily metabolized by CYP2D6.",
            "rapid": "Reduced drug levels. May need dose increase. Monitor for efficacy.",
        },
        "alternatives": ["Sertraline", "Citalopram", "Escitalopram"],
    },
    "amitriptyline": {
        "name": "Amitriptyline",
        "brand_names": ["Elavil"],
        "category": "psychiatry",
        "genes": ["CYP2D6", "CYP2C19"],
        "snps": ["rs3892097", "rs5030655", "rs16947", "rs1065852",
                 "rs4244285", "rs4986893", "rs12248560"],
        "indication": "Depression, neuropathic pain, migraine prevention (TCA)",
        "risk_level": "high",
        "cpic_level": "1A",
        "fda_label": False,
        "pubmed_ids": ["27997040"],
        "interpretations": {
            "normal": "Standard dosing expected to be appropriate.",
            "intermediate": "Increased drug levels. Consider 25% dose reduction. Monitor for side effects (sedation, anticholinergic effects).",
            "poor": "Significantly increased drug levels and toxicity risk. Avoid amitriptyline or reduce dose by 50%. Consider alternative (SSRI, SNRI).",
            "rapid": "Reduced drug levels. May need dose increase. Consider TDM (therapeutic drug monitoring).",
        },
        "alternatives": ["Nortriptyline", "Desipramine", "SSRI/SNRI alternatives"],
    },
    "nortriptyline": {
        "name": "Nortriptyline",
        "brand_names": ["Pamelor"],
        "category": "psychiatry",
        "genes": ["CYP2D6"],
        "snps": ["rs3892097", "rs5030655", "rs16947", "rs1065852"],
        "indication": "Depression, neuropathic pain (TCA)",
        "risk_level": "high",
        "cpic_level": "1A",
        "fda_label": False,
        "pubmed_ids": ["27997040"],
        "interpretations": {
            "normal": "Standard dosing expected to be appropriate. Target therapeutic window 50-150 ng/mL.",
            "intermediate": "Consider 25% dose reduction. Use TDM to guide dosing.",
            "poor": "Significantly increased drug levels. Reduce dose by 50% or avoid. Use TDM. Consider alternative.",
            "rapid": "Reduced drug levels. May need dose increase. Use TDM to verify therapeutic levels.",
        },
        "alternatives": ["Desipramine", "SSRI/SNRI alternatives"],
    },
    "venlafaxine": {
        "name": "Venlafaxine",
        "brand_names": ["Effexor", "Effexor XR"],
        "category": "psychiatry",
        "genes": ["CYP2D6"],
        "snps": ["rs3892097", "rs5030655", "rs16947", "rs1065852"],
        "indication": "Depression, anxiety, panic disorder (SNRI)",
        "risk_level": "moderate",
        "cpic_level": "1A",
        "fda_label": False,
        "pubmed_ids": ["25974703"],
        "interpretations": {
            "normal": "Standard dosing expected to be appropriate.",
            "intermediate": "Minor alteration in parent/metabolite ratio. Standard dosing with monitoring.",
            "poor": "Increased venlafaxine levels, decreased active metabolite. Monitor for side effects (nausea, hypertension). Consider alternative SNRI (desvenlafaxine).",
            "rapid": "Increased conversion to active metabolite. Standard dosing likely appropriate.",
        },
        "alternatives": ["Desvenlafaxine (not CYP2D6 dependent)", "Duloxetine", "SSRI alternatives"],
    },
    "clozapine": {
        "name": "Clozapine",
        "brand_names": ["Clozaril"],
        "category": "psychiatry",
        "genes": ["CYP1A2"],
        "snps": ["rs762551"],
        "indication": "Treatment-resistant schizophrenia",
        "risk_level": "high",
        "cpic_level": "2A",
        "fda_label": False,
        "pubmed_ids": ["30801677"],
        "interpretations": {
            "normal": "Standard dosing with TDM. Smoking induces CYP1A2 and reduces clozapine levels.",
            "intermediate": "Standard dosing with TDM.",
            "poor": "Increased clozapine levels. Higher risk of dose-dependent side effects (seizures, sedation). Consider lower dose.",
            "rapid": "Significantly reduced clozapine levels, especially if smoking. May need higher doses. Mandatory TDM.",
        },
        "alternatives": [],
    },
    # ── ONCOLOGY ──
    "fluorouracil": {
        "name": "5-Fluorouracil / Capecitabine",
        "brand_names": ["Adrucil", "Efudex", "Xeloda"],
        "category": "oncology",
        "genes": ["DPYD"],
        "snps": ["rs3918290", "rs55886062"],
        "indication": "Colorectal, breast, head/neck, and other cancers",
        "risk_level": "high",
        "cpic_level": "1A",
        "fda_label": True,
        "pubmed_ids": ["29152729", "30348537"],
        "interpretations": {
            "normal": "Normal DPD enzyme activity. Standard dosing expected to be tolerated.",
            "intermediate": "Reduced DPD activity. REDUCE DOSE by 50%. Increased risk of severe toxicity (myelosuppression, mucositis, diarrhea, hand-foot syndrome).",
            "poor": "LIFE-THREATENING: Absent or near-absent DPD activity. AVOID 5-FU and capecitabine entirely. Use alternative chemotherapy regimen.",
            "rapid": "Standard dosing appropriate.",
        },
        "alternatives": ["Alternative chemotherapy regimens (consult oncologist)"],
    },
    "azathioprine": {
        "name": "Azathioprine / 6-Mercaptopurine",
        "brand_names": ["Imuran", "Azasan", "Purinethol"],
        "category": "oncology",
        "genes": ["TPMT", "NUDT15"],
        "snps": ["rs1800462", "rs1800460", "rs1142345", "rs116855232"],
        "indication": "Autoimmune diseases, organ transplant, acute lymphoblastic leukemia",
        "risk_level": "high",
        "cpic_level": "1A",
        "fda_label": True,
        "pubmed_ids": ["30447069", "23422873"],
        "interpretations": {
            "normal": "Normal TPMT/NUDT15 activity. Standard dosing appropriate.",
            "intermediate": "Reduced thiopurine metabolism. REDUCE DOSE by 30-50%. Monitor blood counts closely.",
            "poor": "LIFE-THREATENING: Greatly reduced thiopurine metabolism. REDUCE DOSE by 90% or AVOID entirely. Severe myelosuppression risk.",
            "rapid": "Standard dosing appropriate.",
        },
        "alternatives": ["Mycophenolate mofetil (for immunosuppression)"],
    },
    "irinotecan": {
        "name": "Irinotecan",
        "brand_names": ["Camptosar"],
        "category": "oncology",
        "genes": ["UGT1A1"],
        "snps": ["rs8175347"],
        "indication": "Colorectal and other cancers",
        "risk_level": "high",
        "cpic_level": "1A",
        "fda_label": True,
        "pubmed_ids": ["29385237"],
        "interpretations": {
            "normal": "Standard dosing appropriate.",
            "intermediate": "Increased risk of neutropenia. Consider dose reduction if experiencing toxicity.",
            "poor": "High risk of severe neutropenia and diarrhea. REDUCE initial dose by 30%. Monitor closely.",
            "rapid": "Standard dosing appropriate.",
        },
        "alternatives": ["Alternative chemotherapy regimens (consult oncologist)"],
    },
    "tamoxifen": {
        "name": "Tamoxifen",
        "brand_names": ["Nolvadex", "Soltamox"],
        "category": "oncology",
        "genes": ["CYP2D6"],
        "snps": ["rs3892097", "rs5030655", "rs16947", "rs1065852"],
        "indication": "Hormone receptor-positive breast cancer treatment and prevention",
        "risk_level": "high",
        "cpic_level": "1A",
        "fda_label": False,
        "pubmed_ids": ["29385237", "29385238"],
        "interpretations": {
            "normal": "Normal conversion to endoxifen (active metabolite). Standard dosing appropriate.",
            "intermediate": "Reduced endoxifen levels. Consider increased dose (40 mg/day) or alternative (aromatase inhibitor if postmenopausal).",
            "poor": "Significantly reduced endoxifen levels. Tamoxifen may be ineffective. Use alternative hormonal therapy (aromatase inhibitor). Avoid concomitant CYP2D6 inhibitors.",
            "rapid": "Enhanced endoxifen production. Standard dosing appropriate.",
        },
        "alternatives": ["Aromatase inhibitors (anastrozole, letrozole) for postmenopausal patients"],
    },
    # ── INFECTIOUS DISEASE ──
    "abacavir": {
        "name": "Abacavir",
        "brand_names": ["Ziagen", "component of Epzicom/Triumeq"],
        "category": "infectious_disease",
        "genes": ["HLA-B"],
        "snps": ["rs2395029"],
        "indication": "HIV antiretroviral therapy",
        "risk_level": "high",
        "cpic_level": "1A",
        "fda_label": True,
        "pubmed_ids": ["22378157", "24561393"],
        "interpretations": {
            "normal": "HLA-B*57:01 not detected (based on tag SNP). Low risk of hypersensitivity. Standard use appropriate.",
            "intermediate": "N/A",
            "poor": "LIFE-THREATENING: HLA-B*57:01 likely present. AVOID abacavir entirely. Hypersensitivity reaction risk (fever, rash, GI symptoms, potentially fatal). FDA boxed warning.",
            "rapid": "N/A",
        },
        "alternatives": ["Tenofovir-based regimens (TAF or TDF)"],
    },
    "efavirenz": {
        "name": "Efavirenz",
        "brand_names": ["Sustiva", "Stocrin"],
        "category": "infectious_disease",
        "genes": ["CYP2B6"],
        "snps": ["rs3745274"],
        "indication": "HIV antiretroviral therapy (NNRTI)",
        "risk_level": "moderate",
        "cpic_level": "1A",
        "fda_label": False,
        "pubmed_ids": ["31006110"],
        "interpretations": {
            "normal": "Standard dosing (600 mg/day) expected to be tolerated.",
            "intermediate": "Increased drug levels. Higher risk of CNS side effects (vivid dreams, dizziness). Consider reduced dose (400 mg/day).",
            "poor": "Significantly increased drug levels. High risk of CNS toxicity. Reduce to 400 mg/day or consider alternative NNRTI.",
            "rapid": "Standard dosing appropriate.",
        },
        "alternatives": ["Doravirine", "Rilpivirine", "Integrase inhibitor-based regimens"],
    },
    "isoniazid": {
        "name": "Isoniazid",
        "brand_names": ["Nydrazid"],
        "category": "infectious_disease",
        "genes": ["NAT2"],
        "snps": ["rs1801280"],
        "indication": "Tuberculosis treatment and prophylaxis",
        "risk_level": "moderate",
        "cpic_level": "2A",
        "fda_label": False,
        "pubmed_ids": ["30801677"],
        "interpretations": {
            "normal": "Rapid acetylator. Standard dosing appropriate. Lower risk of hepatotoxicity.",
            "intermediate": "Intermediate acetylator. Standard dosing with routine liver function monitoring.",
            "poor": "Slow acetylator. Increased risk of hepatotoxicity and peripheral neuropathy. Monitor liver function. Ensure adequate pyridoxine (vitamin B6) supplementation.",
            "rapid": "Rapid acetylator. Standard or slightly increased dose may be needed.",
        },
        "alternatives": ["Rifampin-based regimens (for latent TB)"],
    },
    # ── GASTROENTEROLOGY ──
    "omeprazole": {
        "name": "Omeprazole / PPIs",
        "brand_names": ["Prilosec", "Nexium (esomeprazole)", "Prevacid (lansoprazole)"],
        "category": "gastroenterology",
        "genes": ["CYP2C19"],
        "snps": ["rs4244285", "rs4986893", "rs12248560"],
        "indication": "Gastric acid suppression for GERD, ulcers, H. pylori eradication",
        "risk_level": "moderate",
        "cpic_level": "1A",
        "fda_label": False,
        "pubmed_ids": ["25974703"],
        "interpretations": {
            "normal": "Standard dosing appropriate. Good acid suppression expected.",
            "intermediate": "Enhanced acid suppression. Standard dosing appropriate.",
            "poor": "Very effective acid suppression at standard doses. May achieve higher H. pylori eradication rates.",
            "rapid": "Reduced acid suppression. May need increased dose (double standard dose) or switch to rabeprazole (less CYP2C19-dependent). Higher H. pylori treatment failure risk.",
        },
        "alternatives": ["Rabeprazole (less CYP2C19 dependent)", "H2 blockers (famotidine)"],
    },
    "methotrexate": {
        "name": "Methotrexate",
        "brand_names": ["Trexall", "Rheumatrex"],
        "category": "gastroenterology",
        "genes": ["MTHFR"],
        "snps": ["rs1801133", "rs1801131"],
        "indication": "Rheumatoid arthritis, psoriasis, some cancers, ectopic pregnancy",
        "risk_level": "moderate",
        "cpic_level": "2A",
        "fda_label": False,
        "pubmed_ids": ["24458010"],
        "interpretations": {
            "normal": "Standard dosing expected to be tolerated. Ensure adequate folate supplementation.",
            "intermediate": "Slightly reduced folate metabolism. Ensure adequate folate/leucovorin supplementation. Monitor for toxicity.",
            "poor": "Significantly reduced folate metabolism. Increased toxicity risk (mucositis, myelosuppression). Ensure high-dose folate supplementation and close monitoring.",
            "rapid": "Standard dosing with routine folate supplementation.",
        },
        "alternatives": ["Leflunomide", "Biologics (for RA)"],
    },
    # ── HORMONAL ──
    "oral_contraceptives": {
        "name": "Oral Contraceptives / Estrogen HRT",
        "brand_names": ["Various combined oral contraceptives", "Estrogen HRT"],
        "category": "hormonal",
        "genes": ["F5"],
        "snps": ["rs6025"],
        "indication": "Contraception, hormone replacement therapy, menstrual regulation",
        "risk_level": "high",
        "cpic_level": "1A",
        "fda_label": True,
        "pubmed_ids": ["24458010"],
        "interpretations": {
            "normal": "No Factor V Leiden detected. Standard thrombosis risk with estrogen use.",
            "intermediate": "Heterozygous Factor V Leiden. 3-8x increased VTE risk with estrogen. Discuss with provider. Consider progestin-only or non-hormonal alternatives.",
            "poor": "Homozygous Factor V Leiden. 50-80x increased VTE risk with estrogen. AVOID estrogen-containing products. Use progestin-only or non-hormonal contraception.",
            "rapid": "N/A",
        },
        "alternatives": ["Progestin-only methods (IUD, minipill)", "Non-hormonal contraception", "Copper IUD"],
    },
    # ── TRANSPLANT ──
    "tacrolimus": {
        "name": "Tacrolimus",
        "brand_names": ["Prograf", "Envarsus XR"],
        "category": "transplant",
        "genes": ["CYP3A5"],
        "snps": ["rs776746"],
        "indication": "Organ transplant immunosuppression",
        "risk_level": "high",
        "cpic_level": "1A",
        "fda_label": False,
        "pubmed_ids": ["25801146"],
        "interpretations": {
            "normal": "CYP3A5 expressor. Requires 1.5-2x higher tacrolimus dose to achieve target trough levels. Start at 0.3 mg/kg/day.",
            "intermediate": "CYP3A5 intermediate expressor. Requires standard to slightly higher doses. Start at 0.25 mg/kg/day. Monitor trough levels.",
            "poor": "CYP3A5 non-expressor (most common in Caucasians). Standard dosing appropriate. Start at 0.15-0.2 mg/kg/day. Monitor trough levels.",
            "rapid": "N/A",
        },
        "alternatives": ["Cyclosporine (monitor CYP3A5 as well)"],
    },
    # ── OTHER ──
    "caffeine": {
        "name": "Caffeine",
        "brand_names": ["Coffee", "Tea", "Energy drinks", "Supplements"],
        "category": "other",
        "genes": ["CYP1A2"],
        "snps": ["rs762551"],
        "indication": "Stimulant / daily consumption",
        "risk_level": "low",
        "cpic_level": "2A",
        "fda_label": False,
        "pubmed_ids": ["16522833"],
        "interpretations": {
            "normal": "Moderate caffeine metabolism. Standard consumption likely well tolerated.",
            "intermediate": "Standard caffeine metabolism.",
            "poor": "Slow caffeine metabolizer. Caffeine stays in your system longer. Higher risk of insomnia, anxiety, and cardiovascular effects from excessive intake. Limit to 1-2 cups/day.",
            "rapid": "Fast caffeine metabolizer. You clear caffeine quickly. May tolerate higher intake without sleep disruption. Coffee may have cardiovascular protective benefits in fast metabolizers.",
        },
        "alternatives": [],
    },
    "atomoxetine": {
        "name": "Atomoxetine",
        "brand_names": ["Strattera"],
        "category": "other",
        "genes": ["CYP2D6"],
        "snps": ["rs3892097", "rs5030655", "rs16947", "rs1065852"],
        "indication": "ADHD (non-stimulant)",
        "risk_level": "moderate",
        "cpic_level": "1A",
        "fda_label": True,
        "pubmed_ids": ["30801677"],
        "interpretations": {
            "normal": "Standard dosing expected to be appropriate.",
            "intermediate": "Slightly increased drug levels. Standard dosing with monitoring.",
            "poor": "Significantly increased drug levels and longer half-life. Start at lowest dose. FDA label recommends dose adjustments for known CYP2D6 PMs.",
            "rapid": "Reduced drug levels. May need higher doses for adequate effect.",
        },
        "alternatives": ["Stimulant medications (methylphenidate, amphetamines)"],
    },
    # Additional drugs to reach 50+
    "celecoxib": {
        "name": "Celecoxib",
        "brand_names": ["Celebrex"],
        "category": "pain",
        "genes": ["CYP2C9"],
        "snps": ["rs1799853", "rs1057910"],
        "indication": "Pain and inflammation (COX-2 inhibitor)",
        "risk_level": "moderate",
        "cpic_level": "1A",
        "fda_label": True,
        "pubmed_ids": ["32189324"],
        "interpretations": {
            "normal": "Standard dosing appropriate.",
            "intermediate": "Increased drug exposure. Consider starting at lowest dose (100 mg/day).",
            "poor": "Significantly increased drug exposure. Start at 25-50% of normal dose. Consider alternative.",
            "rapid": "Standard dosing appropriate.",
        },
        "alternatives": ["Acetaminophen", "Non-CYP2C9 analgesics"],
    },
    "flurbiprofen": {
        "name": "Flurbiprofen",
        "brand_names": ["Ansaid"],
        "category": "pain",
        "genes": ["CYP2C9"],
        "snps": ["rs1799853", "rs1057910"],
        "indication": "Pain and inflammation (NSAID)",
        "risk_level": "moderate",
        "cpic_level": "2A",
        "fda_label": False,
        "pubmed_ids": ["32189324"],
        "interpretations": {
            "normal": "Standard dosing appropriate.",
            "intermediate": "Increased drug exposure. Use lowest effective dose.",
            "poor": "Significantly increased drug exposure. Avoid or use lowest dose for shortest duration.",
            "rapid": "Standard dosing appropriate.",
        },
        "alternatives": ["Acetaminophen"],
    },
    "piroxicam": {
        "name": "Piroxicam",
        "brand_names": ["Feldene"],
        "category": "pain",
        "genes": ["CYP2C9"],
        "snps": ["rs1799853", "rs1057910"],
        "indication": "Arthritis pain (NSAID)",
        "risk_level": "moderate",
        "cpic_level": "2A",
        "fda_label": False,
        "pubmed_ids": ["32189324"],
        "interpretations": {
            "normal": "Standard dosing appropriate.",
            "intermediate": "Increased drug exposure. Use lowest effective dose.",
            "poor": "Significantly increased drug exposure. Avoid or reduce dose by 50%.",
            "rapid": "Standard dosing appropriate.",
        },
        "alternatives": ["Acetaminophen", "Non-CYP2C9 analgesics"],
    },
    "doxepin": {
        "name": "Doxepin",
        "brand_names": ["Sinequan", "Silenor"],
        "category": "psychiatry",
        "genes": ["CYP2D6", "CYP2C19"],
        "snps": ["rs3892097", "rs5030655", "rs16947", "rs1065852",
                 "rs4244285", "rs4986893", "rs12248560"],
        "indication": "Depression, insomnia, anxiety (TCA)",
        "risk_level": "high",
        "cpic_level": "1A",
        "fda_label": False,
        "pubmed_ids": ["27997040"],
        "interpretations": {
            "normal": "Standard dosing expected to be appropriate.",
            "intermediate": "Consider 25% dose reduction. Monitor for anticholinergic effects.",
            "poor": "Significantly increased drug levels. Avoid or reduce dose by 50%. Consider alternative.",
            "rapid": "Reduced drug levels. May need dose increase with TDM.",
        },
        "alternatives": ["SSRI/SNRI alternatives", "Trazodone (for insomnia)"],
    },
    "imipramine": {
        "name": "Imipramine",
        "brand_names": ["Tofranil"],
        "category": "psychiatry",
        "genes": ["CYP2D6", "CYP2C19"],
        "snps": ["rs3892097", "rs5030655", "rs16947", "rs1065852",
                 "rs4244285", "rs4986893", "rs12248560"],
        "indication": "Depression, enuresis (TCA)",
        "risk_level": "high",
        "cpic_level": "1A",
        "fda_label": False,
        "pubmed_ids": ["27997040"],
        "interpretations": {
            "normal": "Standard dosing expected to be appropriate.",
            "intermediate": "Consider 25% dose reduction. Use TDM to guide dosing.",
            "poor": "Significantly increased drug levels. Avoid or reduce dose by 50%. Consider alternative.",
            "rapid": "Reduced drug levels. May need dose increase. Use TDM.",
        },
        "alternatives": ["Nortriptyline", "Desipramine", "SSRI/SNRI alternatives"],
    },
    "trimipramine": {
        "name": "Trimipramine",
        "brand_names": ["Surmontil"],
        "category": "psychiatry",
        "genes": ["CYP2D6"],
        "snps": ["rs3892097", "rs5030655", "rs16947", "rs1065852"],
        "indication": "Depression (TCA)",
        "risk_level": "moderate",
        "cpic_level": "1A",
        "fda_label": False,
        "pubmed_ids": ["27997040"],
        "interpretations": {
            "normal": "Standard dosing expected to be appropriate.",
            "intermediate": "Consider dose reduction with TDM.",
            "poor": "Avoid or reduce dose by 50%. Consider SSRI alternative.",
            "rapid": "May need dose increase. Use TDM.",
        },
        "alternatives": ["SSRI/SNRI alternatives"],
    },
    "fluvoxamine": {
        "name": "Fluvoxamine",
        "brand_names": ["Luvox"],
        "category": "psychiatry",
        "genes": ["CYP2D6"],
        "snps": ["rs3892097", "rs5030655", "rs16947", "rs1065852"],
        "indication": "OCD, social anxiety (SSRI)",
        "risk_level": "moderate",
        "cpic_level": "2A",
        "fda_label": False,
        "pubmed_ids": ["25974703"],
        "interpretations": {
            "normal": "Standard dosing appropriate.",
            "intermediate": "Minor increase in drug levels. Standard dosing with monitoring.",
            "poor": "Increased drug levels. Consider lower starting dose.",
            "rapid": "Reduced drug levels. Monitor for efficacy.",
        },
        "alternatives": ["Sertraline", "Escitalopram"],
    },
    "aripiprazole": {
        "name": "Aripiprazole",
        "brand_names": ["Abilify"],
        "category": "psychiatry",
        "genes": ["CYP2D6"],
        "snps": ["rs3892097", "rs5030655", "rs16947", "rs1065852"],
        "indication": "Schizophrenia, bipolar disorder, augmentation for depression",
        "risk_level": "moderate",
        "cpic_level": "1A",
        "fda_label": True,
        "pubmed_ids": ["30801677"],
        "interpretations": {
            "normal": "Standard dosing appropriate.",
            "intermediate": "Standard dosing appropriate. Monitor for side effects.",
            "poor": "Increased drug levels. FDA label recommends 50% dose reduction (e.g., from 10 mg to 5 mg).",
            "rapid": "Reduced drug levels. Standard dosing; monitor efficacy.",
        },
        "alternatives": ["Quetiapine", "Olanzapine"],
    },
    "ondansetron": {
        "name": "Ondansetron",
        "brand_names": ["Zofran"],
        "category": "other",
        "genes": ["CYP2D6"],
        "snps": ["rs3892097", "rs5030655", "rs16947", "rs1065852"],
        "indication": "Nausea and vomiting prevention",
        "risk_level": "moderate",
        "cpic_level": "1A",
        "fda_label": False,
        "pubmed_ids": ["30801677"],
        "interpretations": {
            "normal": "Standard dosing appropriate.",
            "intermediate": "Standard dosing appropriate.",
            "poor": "Increased drug levels but generally well tolerated at standard doses.",
            "rapid": "Reduced drug exposure. May have decreased antiemetic effect. Consider granisetron as alternative.",
        },
        "alternatives": ["Granisetron", "Palonosetron"],
    },
    "lansoprazole": {
        "name": "Lansoprazole",
        "brand_names": ["Prevacid"],
        "category": "gastroenterology",
        "genes": ["CYP2C19"],
        "snps": ["rs4244285", "rs4986893", "rs12248560"],
        "indication": "Gastric acid suppression (PPI)",
        "risk_level": "low",
        "cpic_level": "1A",
        "fda_label": False,
        "pubmed_ids": ["25974703"],
        "interpretations": {
            "normal": "Standard dosing appropriate.",
            "intermediate": "Enhanced acid suppression at standard doses.",
            "poor": "Excellent acid suppression at standard doses.",
            "rapid": "Reduced acid suppression. Consider double dose or rabeprazole.",
        },
        "alternatives": ["Rabeprazole", "Famotidine"],
    },
    "pantoprazole": {
        "name": "Pantoprazole",
        "brand_names": ["Protonix"],
        "category": "gastroenterology",
        "genes": ["CYP2C19"],
        "snps": ["rs4244285", "rs4986893", "rs12248560"],
        "indication": "Gastric acid suppression (PPI)",
        "risk_level": "low",
        "cpic_level": "1A",
        "fda_label": False,
        "pubmed_ids": ["25974703"],
        "interpretations": {
            "normal": "Standard dosing appropriate.",
            "intermediate": "Enhanced acid suppression at standard doses.",
            "poor": "Excellent acid suppression at standard doses.",
            "rapid": "Reduced acid suppression. May need dose increase.",
        },
        "alternatives": ["Rabeprazole", "Famotidine"],
    },
    "hepatitis_c_treatment": {
        "name": "PEG-IFN / Ribavirin (Hepatitis C)",
        "brand_names": ["Pegasys + Copegus", "PegIntron + Rebetol"],
        "category": "infectious_disease",
        "genes": ["IFNL3"],
        "snps": ["rs12979860"],
        "indication": "Chronic Hepatitis C treatment (historical - DAAs now preferred)",
        "risk_level": "moderate",
        "cpic_level": "1A",
        "fda_label": False,
        "pubmed_ids": ["24458010"],
        "interpretations": {
            "normal": "IFNL3 CC genotype. Higher probability of SVR with interferon-based therapy. DAAs now standard of care.",
            "intermediate": "IFNL3 CT genotype. Intermediate SVR probability with interferon-based therapy. DAAs recommended.",
            "poor": "IFNL3 TT genotype. Lower probability of SVR with interferon-based therapy. DAAs strongly recommended.",
            "rapid": "N/A",
        },
        "alternatives": ["Direct-acting antivirals (sofosbuvir, ledipasvir, glecaprevir/pibrentasvir)"],
    },
    # Additional drugs to reach 50+
    "diclofenac": {
        "name": "Diclofenac",
        "brand_names": ["Voltaren", "Cataflam"],
        "category": "pain",
        "genes": ["CYP2C9"],
        "snps": ["rs1799853", "rs1057910"],
        "indication": "Pain and inflammation (NSAID)",
        "risk_level": "moderate",
        "cpic_level": "2A",
        "fda_label": False,
        "pubmed_ids": ["32189324"],
        "interpretations": {
            "normal": "Standard dosing appropriate.",
            "intermediate": "Increased drug exposure. Use lowest effective dose. Monitor for GI bleeding.",
            "poor": "Significantly increased drug exposure. Avoid or use lowest dose for shortest duration.",
            "rapid": "Standard dosing appropriate.",
        },
        "alternatives": ["Acetaminophen", "Non-CYP2C9 analgesics"],
    },
    "meloxicam": {
        "name": "Meloxicam",
        "brand_names": ["Mobic"],
        "category": "pain",
        "genes": ["CYP2C9"],
        "snps": ["rs1799853", "rs1057910"],
        "indication": "Arthritis pain (NSAID, COX-2 preferential)",
        "risk_level": "moderate",
        "cpic_level": "2A",
        "fda_label": False,
        "pubmed_ids": ["32189324"],
        "interpretations": {
            "normal": "Standard dosing appropriate.",
            "intermediate": "Increased drug exposure. Use lowest effective dose.",
            "poor": "Significantly increased drug exposure. Consider dose reduction by 50%.",
            "rapid": "Standard dosing appropriate.",
        },
        "alternatives": ["Acetaminophen"],
    },
    "desipramine": {
        "name": "Desipramine",
        "brand_names": ["Norpramin"],
        "category": "psychiatry",
        "genes": ["CYP2D6"],
        "snps": ["rs3892097", "rs5030655", "rs16947", "rs1065852"],
        "indication": "Depression (TCA)",
        "risk_level": "high",
        "cpic_level": "1A",
        "fda_label": False,
        "pubmed_ids": ["27997040"],
        "interpretations": {
            "normal": "Standard dosing expected to be appropriate.",
            "intermediate": "Consider 25% dose reduction. Use TDM.",
            "poor": "Significantly increased drug levels. Avoid or reduce dose by 50%. Consider SSRI/SNRI.",
            "rapid": "Reduced drug levels. May need dose increase. Use TDM.",
        },
        "alternatives": ["SSRI/SNRI alternatives"],
    },
    "clomipramine": {
        "name": "Clomipramine",
        "brand_names": ["Anafranil"],
        "category": "psychiatry",
        "genes": ["CYP2D6", "CYP2C19"],
        "snps": ["rs3892097", "rs5030655", "rs16947", "rs1065852",
                 "rs4244285", "rs4986893", "rs12248560"],
        "indication": "OCD, depression (TCA)",
        "risk_level": "high",
        "cpic_level": "1A",
        "fda_label": False,
        "pubmed_ids": ["27997040"],
        "interpretations": {
            "normal": "Standard dosing expected to be appropriate.",
            "intermediate": "Consider 25% dose reduction. Monitor for side effects.",
            "poor": "Significantly increased drug levels. Avoid or reduce dose by 50%. Consider SSRI for OCD.",
            "rapid": "Reduced drug levels. May need dose increase with TDM.",
        },
        "alternatives": ["Fluvoxamine (for OCD)", "SSRI alternatives"],
    },
    "risperidone": {
        "name": "Risperidone",
        "brand_names": ["Risperdal"],
        "category": "psychiatry",
        "genes": ["CYP2D6"],
        "snps": ["rs3892097", "rs5030655", "rs16947", "rs1065852"],
        "indication": "Schizophrenia, bipolar disorder, irritability in autism",
        "risk_level": "moderate",
        "cpic_level": "2A",
        "fda_label": False,
        "pubmed_ids": ["30801677"],
        "interpretations": {
            "normal": "Standard dosing appropriate.",
            "intermediate": "Slightly increased active moiety levels. Standard dosing with monitoring.",
            "poor": "Increased active moiety levels. Consider lower starting dose or use paliperidone.",
            "rapid": "Reduced active moiety. Monitor for efficacy.",
        },
        "alternatives": ["Paliperidone (active metabolite, not CYP2D6 dependent)", "Aripiprazole"],
    },
    "haloperidol": {
        "name": "Haloperidol",
        "brand_names": ["Haldol"],
        "category": "psychiatry",
        "genes": ["CYP2D6"],
        "snps": ["rs3892097", "rs5030655", "rs16947", "rs1065852"],
        "indication": "Schizophrenia, acute psychosis, Tourette syndrome",
        "risk_level": "moderate",
        "cpic_level": "2A",
        "fda_label": False,
        "pubmed_ids": ["30801677"],
        "interpretations": {
            "normal": "Standard dosing appropriate.",
            "intermediate": "Slightly increased drug levels. Standard dosing with monitoring.",
            "poor": "Increased drug levels. Consider lower dose. Monitor for EPS.",
            "rapid": "Reduced drug levels. May need standard to higher doses.",
        },
        "alternatives": ["Aripiprazole", "Quetiapine"],
    },
    "hydrocodone": {
        "name": "Hydrocodone",
        "brand_names": ["Vicodin", "Norco", "Lortab"],
        "category": "pain",
        "genes": ["CYP2D6"],
        "snps": ["rs3892097", "rs5030655", "rs16947", "rs1065852"],
        "indication": "Moderate to severe pain",
        "risk_level": "high",
        "cpic_level": "1A",
        "fda_label": False,
        "pubmed_ids": ["31006110"],
        "interpretations": {
            "normal": "Standard conversion to hydromorphone. Normal analgesic effect expected.",
            "intermediate": "Reduced conversion to active metabolite. May need dose adjustment.",
            "poor": "Greatly reduced activation. Consider alternative opioid or non-opioid analgesic.",
            "rapid": "Increased conversion to active metabolite. Risk of respiratory depression. Use with caution.",
        },
        "alternatives": ["Morphine", "Non-opioid analgesics"],
    },
    "oxycodone": {
        "name": "Oxycodone",
        "brand_names": ["OxyContin", "Percocet"],
        "category": "pain",
        "genes": ["CYP2D6"],
        "snps": ["rs3892097", "rs5030655", "rs16947", "rs1065852"],
        "indication": "Moderate to severe pain",
        "risk_level": "moderate",
        "cpic_level": "2A",
        "fda_label": False,
        "pubmed_ids": ["31006110"],
        "interpretations": {
            "normal": "Standard dosing appropriate. Oxycodone less dependent on CYP2D6 than codeine.",
            "intermediate": "Minor effect on oxymorphone formation. Standard dosing with monitoring.",
            "poor": "Reduced oxymorphone formation but oxycodone retains direct analgesic activity. Monitor pain control.",
            "rapid": "Increased oxymorphone formation. Monitor for excess sedation.",
        },
        "alternatives": ["Morphine", "Non-opioid analgesics"],
    },
    "rabeprazole": {
        "name": "Rabeprazole",
        "brand_names": ["Aciphex"],
        "category": "gastroenterology",
        "genes": ["CYP2C19"],
        "snps": ["rs4244285", "rs4986893", "rs12248560"],
        "indication": "Gastric acid suppression (PPI, less CYP2C19-dependent)",
        "risk_level": "low",
        "cpic_level": "1A",
        "fda_label": False,
        "pubmed_ids": ["25974703"],
        "interpretations": {
            "normal": "Standard dosing appropriate. Good choice regardless of CYP2C19 status.",
            "intermediate": "Standard dosing appropriate. Less affected by CYP2C19 variation.",
            "poor": "Standard dosing appropriate. Rabeprazole is preferred PPI for CYP2C19 PMs.",
            "rapid": "Better acid suppression than other PPIs in rapid metabolizers. Preferred PPI choice.",
        },
        "alternatives": [],
    },
}


# ─────────────────────────────────────────────────────────────────────────────
# DRUG INTERACTION DATABASE
# ─────────────────────────────────────────────────────────────────────────────

DRUG_INTERACTIONS: List[Dict[str, Any]] = [
    {
        "drugs": ["Codeine", "Paroxetine"],
        "gene": "CYP2D6",
        "risk": "Paroxetine is a strong CYP2D6 inhibitor, further reducing codeine-to-morphine "
                "conversion. Combined with a reduced metabolizer genotype, codeine will be "
                "effectively ineffective for pain relief.",
        "recommendation": "Avoid this combination. Use alternative analgesic or alternative SSRI.",
    },
    {
        "drugs": ["Codeine", "Fluoxetine"],
        "gene": "CYP2D6",
        "risk": "Fluoxetine is a strong CYP2D6 inhibitor, further reducing codeine activation. "
                "Pain relief will be significantly diminished.",
        "recommendation": "Avoid this combination. Use alternative analgesic.",
    },
    {
        "drugs": ["Tamoxifen", "Paroxetine"],
        "gene": "CYP2D6",
        "risk": "Paroxetine inhibits CYP2D6, reducing tamoxifen conversion to its active "
                "metabolite endoxifen. This may reduce tamoxifen's cancer-preventive efficacy.",
        "recommendation": "Avoid this combination. Use SSRI that does not inhibit CYP2D6 "
                         "(e.g., citalopram, venlafaxine).",
    },
    {
        "drugs": ["Tamoxifen", "Fluoxetine"],
        "gene": "CYP2D6",
        "risk": "Fluoxetine inhibits CYP2D6, reducing endoxifen levels and potentially "
                "decreasing tamoxifen efficacy for breast cancer prevention.",
        "recommendation": "Avoid this combination. Use alternative antidepressant.",
    },
    {
        "drugs": ["Warfarin", "Ibuprofen"],
        "gene": "CYP2C9",
        "risk": "Both metabolized by CYP2C9. Combined use increases bleeding risk, "
                "especially in CYP2C9 intermediate/poor metabolizers.",
        "recommendation": "Use acetaminophen for pain when on warfarin. If NSAID needed, "
                         "use lowest dose for shortest duration with close INR monitoring.",
    },
    {
        "drugs": ["Clopidogrel", "Omeprazole / PPIs"],
        "gene": "CYP2C19",
        "risk": "Omeprazole inhibits CYP2C19, potentially reducing clopidogrel activation. "
                "This is especially relevant for CYP2C19 intermediate metabolizers.",
        "recommendation": "If PPI needed with clopidogrel, prefer pantoprazole (weaker CYP2C19 inhibitor).",
    },
    {
        "drugs": ["Metoprolol", "Paroxetine"],
        "gene": "CYP2D6",
        "risk": "Paroxetine inhibits CYP2D6, increasing metoprolol levels. Risk of "
                "bradycardia and hypotension, especially in CYP2D6 poor metabolizers.",
        "recommendation": "Monitor heart rate and blood pressure. Consider atenolol or "
                         "bisoprolol (not CYP2D6 dependent).",
    },
    {
        "drugs": ["Tramadol", "Paroxetine"],
        "gene": "CYP2D6",
        "risk": "Paroxetine inhibits CYP2D6, reducing tramadol activation. Also increases "
                "serotonin syndrome risk when combined with SSRIs.",
        "recommendation": "Avoid this combination due to reduced efficacy and serotonin syndrome risk.",
    },
]


# ─────────────────────────────────────────────────────────────────────────────
# BIOMARKER RECOMMENDATIONS
# ─────────────────────────────────────────────────────────────────────────────

BIOMARKER_RECOMMENDATIONS: List[Dict[str, Any]] = [
    {
        "test": "INR / Prothrombin Time (PT)",
        "genes": ["CYP2C9", "VKORC1"],
        "reason": "Your CYP2C9 and VKORC1 status affects warfarin metabolism and sensitivity",
        "frequency": "If prescribed warfarin: weekly during initiation, then monthly once stable",
        "condition": "non_normal",  # only show if CYP2C9 or VKORC1 is non-normal
    },
    {
        "test": "Homocysteine level",
        "genes": ["MTHFR"],
        "reason": "MTHFR variant may affect folate metabolism and homocysteine levels",
        "frequency": "Annually or as directed",
        "condition": "non_normal",
    },
    {
        "test": "Complete Blood Count (CBC)",
        "genes": ["DPYD", "TPMT", "NUDT15", "UGT1A1"],
        "reason": "Your pharmacogenomic profile warrants close blood count monitoring if chemotherapy is prescribed",
        "frequency": "Per oncology protocol; more frequently if carrying variant alleles",
        "condition": "non_normal",
    },
    {
        "test": "Liver Function Tests (ALT, AST)",
        "genes": ["NAT2"],
        "reason": "NAT2 slow acetylator status increases hepatotoxicity risk with isoniazid",
        "frequency": "Monthly during isoniazid treatment",
        "condition": "non_normal",
    },
    {
        "test": "Therapeutic Drug Monitoring (TDM) for TCAs",
        "genes": ["CYP2D6"],
        "reason": "CYP2D6 status affects tricyclic antidepressant levels",
        "frequency": "At steady state (5-7 days) and after dose changes",
        "condition": "non_normal",
    },
    {
        "test": "Tacrolimus trough levels",
        "genes": ["CYP3A5"],
        "reason": "CYP3A5 expressor status significantly affects tacrolimus dosing requirements",
        "frequency": "Per transplant protocol",
        "condition": "always",
    },
    {
        "test": "Lipid panel + CK (Creatine Kinase)",
        "genes": ["SLCO1B1"],
        "reason": "SLCO1B1 variant increases statin-induced myopathy risk",
        "frequency": "Baseline and 4-12 weeks after starting statin; report muscle symptoms immediately",
        "condition": "non_normal",
    },
    {
        "test": "D-dimer / Thrombophilia panel",
        "genes": ["F5"],
        "reason": "Factor V Leiden increases thrombosis risk, especially with estrogen use",
        "frequency": "Before starting estrogen therapy; as clinically indicated",
        "condition": "non_normal",
    },
]


# ─────────────────────────────────────────────────────────────────────────────
# METABOLIZER STATUS CLASSIFICATION
# ─────────────────────────────────────────────────────────────────────────────

def _classify_metabolizer(
    gene: str, found_variants: Dict[str, Dict[str, Any]]
) -> Dict[str, Any]:
    """
    Classify the metabolizer status for a given gene based on detected variants.

    Args:
        gene: Gene symbol (e.g., "CYP2D6")
        found_variants: Dict of rsid -> {"genotype": str, "risk_allele_count": int,
                                          "snp_info": dict}

    Returns:
        Dict with "status", "label", and "clinical_significance" keys.
    """
    # Count loss-of-function (LOF) and gain-of-function (GOF) alleles
    lof_alleles = 0
    gof_alleles = 0
    decreased_alleles = 0

    for rsid, vdata in found_variants.items():
        snp_info = vdata["snp_info"]
        risk_count = vdata["risk_allele_count"]
        func = snp_info.get("function", "").lower()

        if "loss of function" in func:
            lof_alleles += risk_count
        elif "decreased function" in func or "decreased" in func:
            decreased_alleles += risk_count
        elif "increased function" in func or "ultra-rapid" in func or "inducible" in func:
            gof_alleles += risk_count

    # Special handling per gene
    if gene == "CYP2C19":
        return _classify_cyp2c19(lof_alleles, gof_alleles, decreased_alleles)
    elif gene == "CYP2D6":
        return _classify_cyp2d6(lof_alleles, gof_alleles, decreased_alleles)
    elif gene == "CYP1A2":
        return _classify_cyp1a2(gof_alleles)
    elif gene == "CYP3A5":
        return _classify_cyp3a5(lof_alleles)
    elif gene == "VKORC1":
        return _classify_vkorc1(found_variants)
    elif gene == "F5":
        return _classify_factor_v(found_variants)
    elif gene == "HLA-B":
        return _classify_hla_b(found_variants)
    elif gene in ("OPRM1", "COMT"):
        return _classify_receptor_gene(gene, found_variants)
    elif gene in ("IFNL3", "IFNL4"):
        return _classify_ifnl(gene, found_variants)
    elif gene == "NAT2":
        return _classify_nat2(found_variants)
    else:
        # Generic classification for DPYD, TPMT, NUDT15, CYP2C9, UGT1A1,
        # SLCO1B1, ABCB1, MTHFR, CYP2B6
        return _classify_generic(gene, lof_alleles, gof_alleles, decreased_alleles)


def _classify_cyp2c19(lof: int, gof: int, dec: int) -> Dict[str, Any]:
    """
    CYP2C19 classification per CPIC:
    - *1/*1 = Normal (EM)
    - *1/*17 or *17/*17 = Rapid/Ultra-Rapid
    - *1/*2 or *1/*3 = Intermediate
    - *2/*2, *2/*3, *3/*3 = Poor
    - *2/*17 = Intermediate (conflicting, but CPIC classifies as IM)
    """
    total_lof = lof + dec
    if total_lof >= 2:
        return {
            "status": "poor", "label": "Poor Metabolizer",
            "clinical_significance": (
                "You have two loss-of-function alleles in CYP2C19. You metabolize "
                "CYP2C19 substrates much slower than normal. This significantly affects "
                "clopidogrel (reduced activation), SSRIs (increased levels), and PPIs (enhanced effect)."
            ),
        }
    elif total_lof == 1 and gof == 0:
        return {
            "status": "intermediate", "label": "Intermediate Metabolizer",
            "clinical_significance": (
                "You carry one loss-of-function CYP2C19 allele. You metabolize CYP2C19 "
                "substrates slower than normal. Dose adjustments may be needed for "
                "clopidogrel and some SSRIs."
            ),
        }
    elif total_lof == 1 and gof >= 1:
        # *2/*17 scenario - CPIC classifies as intermediate
        return {
            "status": "intermediate", "label": "Intermediate Metabolizer",
            "clinical_significance": (
                "You carry both a loss-of-function and a gain-of-function CYP2C19 allele. "
                "The net effect is intermediate metabolism. Clinical guidance follows "
                "intermediate metabolizer recommendations."
            ),
        }
    elif gof >= 2:
        return {
            "status": "rapid", "label": "Ultra-Rapid Metabolizer",
            "clinical_significance": (
                "You carry two increased-function CYP2C19*17 alleles. You metabolize "
                "CYP2C19 substrates much faster than normal. PPIs may be less effective "
                "at standard doses. Clopidogrel activation is enhanced."
            ),
        }
    elif gof == 1:
        return {
            "status": "rapid", "label": "Rapid Metabolizer",
            "clinical_significance": (
                "You carry one CYP2C19*17 increased-function allele. You metabolize "
                "CYP2C19 substrates faster than normal. PPIs may need higher doses."
            ),
        }
    else:
        return {
            "status": "normal", "label": "Normal Metabolizer",
            "clinical_significance": (
                "No CYP2C19 loss-of-function or gain-of-function variants detected. "
                "Standard metabolism of CYP2C19 substrates expected."
            ),
        }


def _classify_cyp2d6(lof: int, gof: int, dec: int) -> Dict[str, Any]:
    """
    CYP2D6 classification (approximate - full star allele calling requires
    copy number and structural variant data not available from consumer tests).
    Based on available SNPs: *4 (LOF), *6 (LOF), *10 (decreased), *2 (normal).
    """
    total_reduced = lof + dec
    if lof >= 2:
        return {
            "status": "poor", "label": "Poor Metabolizer",
            "clinical_significance": (
                "You likely have significantly reduced CYP2D6 enzyme activity. This affects "
                "many medications including codeine (ineffective), tamoxifen (reduced efficacy), "
                "TCAs (increased toxicity risk), and SSRIs. Note: CYP2D6 cannot be fully "
                "characterized from SNP data alone; clinical genotyping may detect additional variants."
            ),
        }
    elif total_reduced >= 2 and lof >= 1:
        return {
            "status": "intermediate", "label": "Intermediate Metabolizer",
            "clinical_significance": (
                "You carry CYP2D6 variants suggesting reduced enzyme activity. You likely "
                "metabolize CYP2D6 substrates slower than normal. Dose adjustments may be "
                "needed for codeine, some antidepressants, and tamoxifen."
            ),
        }
    elif total_reduced == 1:
        return {
            "status": "intermediate", "label": "Intermediate Metabolizer",
            "clinical_significance": (
                "You carry one CYP2D6 reduced-function allele. You have intermediate "
                "enzyme activity. Some medications may need dose adjustment."
            ),
        }
    elif dec >= 2:
        return {
            "status": "intermediate", "label": "Intermediate Metabolizer",
            "clinical_significance": (
                "You carry CYP2D6 decreased-function variants. Your enzyme activity is "
                "reduced but not absent. Consider dose adjustments for sensitive CYP2D6 substrates."
            ),
        }
    else:
        return {
            "status": "normal", "label": "Normal Metabolizer",
            "clinical_significance": (
                "No CYP2D6 loss-of-function variants detected from available SNPs. "
                "Standard CYP2D6 metabolism expected. Note: Consumer DNA tests cannot detect "
                "gene duplications (ultra-rapid) or all known variants; clinical genotyping "
                "provides more comprehensive results."
            ),
        }


def _classify_cyp1a2(gof: int) -> Dict[str, Any]:
    """CYP1A2*1F classification: AA = ultra-rapid inducible, AC = intermediate, CC = normal."""
    if gof >= 2:
        return {
            "status": "rapid", "label": "Ultra-Rapid Metabolizer (Inducible)",
            "clinical_significance": (
                "You are homozygous for CYP1A2*1F (AA). Your CYP1A2 enzyme is highly "
                "inducible by smoking, cruciferous vegetables, and charbroiled foods. "
                "You metabolize caffeine and clozapine faster, especially if you smoke."
            ),
        }
    elif gof == 1:
        return {
            "status": "intermediate", "label": "Intermediate Metabolizer",
            "clinical_significance": (
                "You carry one CYP1A2*1F allele (AC). Intermediate inducibility of "
                "CYP1A2 enzyme. Moderate caffeine metabolism."
            ),
        }
    else:
        return {
            "status": "normal", "label": "Normal Metabolizer",
            "clinical_significance": (
                "You are CYP1A2*1/*1 (CC). Normal, non-inducible CYP1A2 activity. "
                "Slower caffeine metabolism compared to *1F carriers."
            ),
        }


def _classify_cyp3a5(lof: int) -> Dict[str, Any]:
    """CYP3A5*3 classification: *3/*3 = non-expressor, *1/*3 = intermediate, *1/*1 = expressor."""
    if lof >= 2:
        return {
            "status": "poor", "label": "Non-Expressor",
            "clinical_significance": (
                "You are CYP3A5*3/*3 (non-expressor). You do not produce functional "
                "CYP3A5 enzyme. This is the most common phenotype in Caucasians (80-90%). "
                "Standard tacrolimus dosing applies."
            ),
        }
    elif lof == 1:
        return {
            "status": "intermediate", "label": "Intermediate Expressor",
            "clinical_significance": (
                "You are CYP3A5*1/*3 (heterozygous). You produce some CYP3A5 enzyme. "
                "May require moderately higher tacrolimus doses."
            ),
        }
    else:
        return {
            "status": "normal", "label": "Expressor",
            "clinical_significance": (
                "You are CYP3A5*1/*1 (expressor). You produce full CYP3A5 enzyme. "
                "You will likely need significantly higher tacrolimus doses (1.5-2x) "
                "to achieve therapeutic levels."
            ),
        }


def _classify_vkorc1(found_variants: Dict[str, Dict[str, Any]]) -> Dict[str, Any]:
    """VKORC1 -1639G>A: AA = high sensitivity, GA = intermediate, GG = normal."""
    risk_count = 0
    for rsid, vdata in found_variants.items():
        risk_count += vdata["risk_allele_count"]

    if risk_count >= 2:
        return {
            "status": "poor", "label": "High Warfarin Sensitivity",
            "clinical_significance": (
                "You are homozygous for VKORC1 -1639A (AA). You have high sensitivity "
                "to warfarin and typically require significantly lower doses (about 50% "
                "of the standard dose). Close INR monitoring is essential."
            ),
        }
    elif risk_count == 1:
        return {
            "status": "intermediate", "label": "Intermediate Warfarin Sensitivity",
            "clinical_significance": (
                "You are heterozygous for VKORC1 -1639G>A. You have intermediate warfarin "
                "sensitivity and may need a moderate dose reduction (about 25%)."
            ),
        }
    else:
        return {
            "status": "normal", "label": "Normal Warfarin Sensitivity",
            "clinical_significance": (
                "You are VKORC1 -1639GG. Standard warfarin sensitivity. Dosing based "
                "on clinical algorithms (which include CYP2C9 status as well)."
            ),
        }


def _classify_factor_v(found_variants: Dict[str, Dict[str, Any]]) -> Dict[str, Any]:
    """Factor V Leiden: TT = homozygous, GT = heterozygous, GG = normal."""
    risk_count = 0
    for rsid, vdata in found_variants.items():
        risk_count += vdata["risk_allele_count"]

    if risk_count >= 2:
        return {
            "status": "poor", "label": "Homozygous Factor V Leiden",
            "clinical_significance": (
                "You are homozygous for Factor V Leiden (TT). You have a significantly "
                "elevated risk of venous thromboembolism (50-80x baseline). Estrogen-containing "
                "contraceptives and HRT should be AVOIDED. Discuss anticoagulation strategy "
                "with your healthcare provider."
            ),
        }
    elif risk_count == 1:
        return {
            "status": "intermediate", "label": "Heterozygous Factor V Leiden",
            "clinical_significance": (
                "You are heterozygous for Factor V Leiden (GT). You have a moderately "
                "elevated risk of venous thromboembolism (3-8x baseline). Estrogen use "
                "further increases this risk. Discuss alternatives with your provider."
            ),
        }
    else:
        return {
            "status": "normal", "label": "No Factor V Leiden",
            "clinical_significance": (
                "No Factor V Leiden variant detected. Standard thrombosis risk from "
                "a Factor V perspective."
            ),
        }


def _classify_hla_b(found_variants: Dict[str, Dict[str, Any]]) -> Dict[str, Any]:
    """HLA-B*57:01 tag SNP: presence of risk allele suggests HLA-B*57:01 carrier status."""
    risk_count = 0
    for rsid, vdata in found_variants.items():
        risk_count += vdata["risk_allele_count"]

    if risk_count >= 1:
        return {
            "status": "poor", "label": "HLA-B*57:01 Likely Positive",
            "clinical_significance": (
                "Tag SNP rs2395029 suggests you likely carry HLA-B*57:01. "
                "AVOID abacavir (HIV medication) due to LIFE-THREATENING hypersensitivity "
                "reaction risk. Confirmatory HLA typing recommended if abacavir is considered."
            ),
        }
    else:
        return {
            "status": "normal", "label": "HLA-B*57:01 Likely Negative",
            "clinical_significance": (
                "Tag SNP does not indicate HLA-B*57:01 carrier status. Low risk of "
                "abacavir hypersensitivity. Note: confirmatory HLA typing is recommended "
                "before abacavir initiation per clinical guidelines."
            ),
        }


def _classify_receptor_gene(
    gene: str, found_variants: Dict[str, Dict[str, Any]]
) -> Dict[str, Any]:
    """Classify OPRM1 and COMT (not traditional metabolizers)."""
    risk_count = 0
    for rsid, vdata in found_variants.items():
        risk_count += vdata["risk_allele_count"]

    if gene == "OPRM1":
        if risk_count >= 2:
            return {
                "status": "poor", "label": "Altered Opioid Response (Homozygous G)",
                "clinical_significance": (
                    "You are homozygous for OPRM1 118G. You may have reduced opioid "
                    "receptor binding affinity, potentially requiring higher opioid doses "
                    "for pain relief. Response to naltrexone for alcohol dependence may be enhanced."
                ),
            }
        elif risk_count == 1:
            return {
                "status": "intermediate", "label": "Altered Opioid Response (Heterozygous)",
                "clinical_significance": (
                    "You are heterozygous for OPRM1 A118G. You may have mildly altered "
                    "opioid receptor function. Clinical significance is moderate."
                ),
            }
        else:
            return {
                "status": "normal", "label": "Typical Opioid Response",
                "clinical_significance": (
                    "Standard OPRM1 genotype (AA). Typical opioid receptor binding expected."
                ),
            }
    else:  # COMT
        if risk_count >= 2:
            return {
                "status": "poor", "label": "Low COMT Activity (Met/Met)",
                "clinical_significance": (
                    "You are homozygous for COMT Val158Met (Met/Met). You have lower "
                    "COMT enzyme activity, resulting in higher dopamine and catecholamine "
                    "levels. May have increased pain sensitivity but better cognitive performance "
                    "under low-stress conditions. May respond differently to pain medications."
                ),
            }
        elif risk_count == 1:
            return {
                "status": "intermediate", "label": "Intermediate COMT Activity (Val/Met)",
                "clinical_significance": (
                    "You are heterozygous for COMT Val158Met. You have intermediate "
                    "COMT enzyme activity. Balanced dopamine metabolism."
                ),
            }
        else:
            return {
                "status": "normal", "label": "Normal COMT Activity (Val/Val)",
                "clinical_significance": (
                    "You are COMT Val/Val. Normal COMT enzyme activity with efficient "
                    "catecholamine clearance. May have higher stress resilience but slightly "
                    "lower baseline dopamine levels."
                ),
            }


def _classify_ifnl(
    gene: str, found_variants: Dict[str, Dict[str, Any]]
) -> Dict[str, Any]:
    """Classify IFNL3/IFNL4 for HCV treatment response."""
    risk_count = 0
    for rsid, vdata in found_variants.items():
        risk_count += vdata["risk_allele_count"]

    if risk_count >= 2:
        return {
            "status": "poor", "label": "Unfavorable HCV Treatment Response",
            "clinical_significance": (
                f"Homozygous for {gene} unfavorable allele. Lower probability of "
                "sustained virological response (SVR) with interferon-based therapy. "
                "Direct-acting antivirals (DAAs) are now standard of care and highly "
                "effective regardless of IFNL genotype."
            ),
        }
    elif risk_count == 1:
        return {
            "status": "intermediate", "label": "Intermediate HCV Treatment Response",
            "clinical_significance": (
                f"Heterozygous for {gene} variant. Intermediate SVR probability with "
                "interferon-based therapy. DAAs recommended as first-line treatment."
            ),
        }
    else:
        return {
            "status": "normal", "label": "Favorable HCV Treatment Response",
            "clinical_significance": (
                f"Favorable {gene} genotype. Higher probability of SVR with "
                "interferon-based therapy. DAAs remain preferred treatment."
            ),
        }


def _classify_nat2(found_variants: Dict[str, Dict[str, Any]]) -> Dict[str, Any]:
    """NAT2 slow vs rapid acetylator classification."""
    risk_count = 0
    for rsid, vdata in found_variants.items():
        risk_count += vdata["risk_allele_count"]

    if risk_count >= 2:
        return {
            "status": "poor", "label": "Slow Acetylator",
            "clinical_significance": (
                "You are likely a slow NAT2 acetylator. Isoniazid and other NAT2 "
                "substrates are cleared more slowly, increasing the risk of "
                "hepatotoxicity. Monitor liver function during isoniazid treatment."
            ),
        }
    elif risk_count == 1:
        return {
            "status": "intermediate", "label": "Intermediate Acetylator",
            "clinical_significance": (
                "You carry one NAT2 slow-acetylator allele. Intermediate acetylation "
                "rate. Standard dosing with routine monitoring."
            ),
        }
    else:
        return {
            "status": "normal", "label": "Rapid Acetylator",
            "clinical_significance": (
                "No NAT2 slow-acetylator variants detected. Standard isoniazid "
                "dosing appropriate."
            ),
        }


def _classify_generic(gene: str, lof: int, gof: int, dec: int) -> Dict[str, Any]:
    """Generic metabolizer classification for genes not requiring special logic."""
    total_reduced = lof + dec

    if lof >= 2:
        return {
            "status": "poor", "label": "Poor Metabolizer",
            "clinical_significance": (
                f"Two loss-of-function alleles detected in {gene}. Significantly "
                f"reduced enzyme activity. Dose adjustments needed for {gene} substrates."
            ),
        }
    elif total_reduced >= 2:
        return {
            "status": "poor", "label": "Poor Metabolizer",
            "clinical_significance": (
                f"Multiple reduced-function alleles detected in {gene}. Significantly "
                f"reduced enzyme activity."
            ),
        }
    elif total_reduced == 1:
        return {
            "status": "intermediate", "label": "Intermediate Metabolizer",
            "clinical_significance": (
                f"One reduced-function allele detected in {gene}. Moderately reduced "
                f"enzyme activity. Some medications may need dose adjustment."
            ),
        }
    elif gof >= 2:
        return {
            "status": "rapid", "label": "Ultra-Rapid Metabolizer",
            "clinical_significance": (
                f"Two increased-function alleles detected in {gene}. Increased "
                f"enzyme activity. Standard doses may be less effective."
            ),
        }
    elif gof == 1:
        return {
            "status": "rapid", "label": "Rapid Metabolizer",
            "clinical_significance": (
                f"One increased-function allele detected in {gene}. Slightly "
                f"increased enzyme activity."
            ),
        }
    else:
        return {
            "status": "normal", "label": "Normal Metabolizer",
            "clinical_significance": (
                f"No actionable variants detected in {gene}. Standard enzyme "
                f"activity expected."
            ),
        }


# ─────────────────────────────────────────────────────────────────────────────
# DRUG ASSESSMENT
# ─────────────────────────────────────────────────────────────────────────────

def _assess_drug(
    drug_key: str,
    drug_info: Dict[str, Any],
    enzyme_panel: Dict[str, Dict[str, Any]],
) -> Dict[str, Any]:
    """
    Generate a drug assessment card based on the user's enzyme panel.

    For drugs that depend on multiple genes (e.g., warfarin: CYP2C9 + VKORC1,
    amitriptyline: CYP2D6 + CYP2C19), we use the worst-case metabolizer status.
    """
    relevant_genes = drug_info["genes"]
    statuses = []

    for gene in relevant_genes:
        if gene in enzyme_panel:
            statuses.append(enzyme_panel[gene]["status"])

    # Determine overall metabolizer status (worst-case for multi-gene drugs)
    status_priority = {"poor": 0, "intermediate": 1, "rapid": 2, "normal": 3}
    if statuses:
        worst_status = min(statuses, key=lambda s: status_priority.get(s, 3))
    else:
        worst_status = "normal"

    # Get interpretation
    interpretations = drug_info.get("interpretations", {})
    recommendation = interpretations.get(worst_status, interpretations.get("normal", ""))

    # Determine traffic light
    if worst_status == "normal":
        traffic_light = "green"
    elif worst_status in ("intermediate", "rapid"):
        traffic_light = "yellow"
    else:
        traffic_light = "red"

    # Determine risk level based on drug risk + metabolizer status
    base_risk = drug_info.get("risk_level", "low")
    if worst_status == "poor" and base_risk == "high":
        effective_risk = "high"
    elif worst_status == "poor" and base_risk == "moderate":
        effective_risk = "high"
    elif worst_status in ("intermediate", "rapid") and base_risk == "high":
        effective_risk = "moderate"
    elif worst_status == "normal":
        effective_risk = "low"
    else:
        effective_risk = base_risk

    # Check for critical alert
    critical_alert = False
    if worst_status == "poor" and base_risk == "high":
        critical_alert = True
    if "LIFE-THREATENING" in recommendation or "AVOID" in recommendation:
        critical_alert = True

    # Dose guidance
    if worst_status == "poor":
        dose_guidance = "Dose reduction likely needed or avoid drug entirely. Consult pharmacogenomic guidance."
    elif worst_status == "intermediate":
        dose_guidance = "Consider dose adjustment. Monitor for efficacy and side effects."
    elif worst_status == "rapid":
        dose_guidance = "Standard or increased dose may be needed. Monitor for efficacy."
    else:
        dose_guidance = "Standard dosing protocol expected to be appropriate."

    return {
        "drug_name": drug_info["name"],
        "brand_names": drug_info.get("brand_names", []),
        "category": drug_info["category"],
        "indication": drug_info.get("indication", ""),
        "genes": relevant_genes,
        "metabolizer_status": worst_status,
        "risk_level": effective_risk,
        "traffic_light": traffic_light,
        "recommendation": recommendation,
        "dose_guidance": dose_guidance,
        "alternatives": drug_info.get("alternatives", []),
        "evidence_level": drug_info.get("cpic_level", ""),
        "source": "CPIC Guideline" if drug_info.get("cpic_level", "").startswith("1") else "PharmGKB",
        "fda_label": drug_info.get("fda_label", False),
        "critical_alert": critical_alert,
        "pubmed_ids": drug_info.get("pubmed_ids", []),
    }


# ─────────────────────────────────────────────────────────────────────────────
# CRITICAL ALERTS
# ─────────────────────────────────────────────────────────────────────────────

def _generate_critical_alerts(
    drug_cards: List[Dict[str, Any]],
    enzyme_panel: Dict[str, Dict[str, Any]],
) -> List[Dict[str, Any]]:
    """
    Generate critical alerts for high-risk drug/gene combinations.
    """
    alerts: List[Dict[str, Any]] = []

    for card in drug_cards:
        if not card.get("critical_alert"):
            continue

        recommendation = card.get("recommendation", "")

        if "LIFE-THREATENING" in recommendation:
            severity = "life_threatening"
        elif "AVOID" in recommendation:
            severity = "contraindicated"
        else:
            severity = "high_risk"

        # Build action text
        if "AVOID" in recommendation:
            action = f"AVOID {card['drug_name']} or use with extreme caution under specialist supervision"
        elif "REDUCE DOSE" in recommendation:
            action = f"REDUCE DOSE of {card['drug_name']} per pharmacogenomic guidelines"
        else:
            action = f"Consult pharmacogenomics specialist before using {card['drug_name']}"

        gene_str = ", ".join(card.get("genes", []))

        alerts.append({
            "severity": severity,
            "drug": card["drug_name"],
            "gene": gene_str,
            "metabolizer_status": card["metabolizer_status"],
            "message": recommendation,
            "action": action,
        })

    # Sort by severity
    severity_order = {"life_threatening": 0, "contraindicated": 1, "high_risk": 2}
    alerts.sort(key=lambda a: severity_order.get(a["severity"], 3))

    return alerts


# ─────────────────────────────────────────────────────────────────────────────
# DRUG INTERACTION RISKS
# ─────────────────────────────────────────────────────────────────────────────

def _generate_interaction_risks(
    drug_cards: List[Dict[str, Any]],
    enzyme_panel: Dict[str, Dict[str, Any]],
) -> List[Dict[str, Any]]:
    """
    Generate drug-drug interaction warnings relevant to the user's genotype.
    Only include interactions where the user's metabolizer status makes the
    interaction clinically more significant.
    """
    results: List[Dict[str, Any]] = []

    for interaction in DRUG_INTERACTIONS:
        gene = interaction["gene"]
        if gene not in enzyme_panel:
            continue

        status = enzyme_panel[gene]["status"]
        # Interactions are more significant for non-normal metabolizers
        if status == "normal":
            # Still include but note it's less concerning
            relevance = "Note: Your normal metabolizer status reduces but does not eliminate this interaction risk."
        else:
            relevance = (
                f"Your {enzyme_panel[gene]['label']} status for {gene} makes this "
                f"interaction particularly significant."
            )

        results.append({
            "drugs": interaction["drugs"],
            "gene": gene,
            "metabolizer_status": status,
            "risk": interaction["risk"],
            "relevance": relevance,
            "recommendation": interaction["recommendation"],
        })

    return results


# ─────────────────────────────────────────────────────────────────────────────
# BIOMARKER PLAN
# ─────────────────────────────────────────────────────────────────────────────

def _generate_biomarker_plan(
    enzyme_panel: Dict[str, Dict[str, Any]],
) -> List[Dict[str, Any]]:
    """
    Generate recommended monitoring tests based on the user's enzyme panel.
    """
    plan: List[Dict[str, Any]] = []

    for rec in BIOMARKER_RECOMMENDATIONS:
        relevant_genes = rec["genes"]
        condition = rec.get("condition", "non_normal")

        # Check if any relevant gene has non-normal status
        include = False
        gene_statuses = []
        for gene in relevant_genes:
            if gene in enzyme_panel:
                status = enzyme_panel[gene]["status"]
                gene_statuses.append(f"{gene}: {enzyme_panel[gene]['label']}")
                if condition == "always":
                    include = True
                elif condition == "non_normal" and status != "normal":
                    include = True

        if include:
            plan.append({
                "test": rec["test"],
                "reason": rec["reason"],
                "frequency": rec["frequency"],
                "gene": ", ".join(relevant_genes),
                "your_status": "; ".join(gene_statuses) if gene_statuses else "Not tested",
            })

    return plan


# ─────────────────────────────────────────────────────────────────────────────
# CATEGORY SUMMARIES
# ─────────────────────────────────────────────────────────────────────────────

CATEGORY_LABELS = {
    "cardiovascular": "Cardiovascular",
    "pain": "Pain Management",
    "psychiatry": "Mental Health / Psychiatry",
    "oncology": "Oncology",
    "infectious_disease": "Infectious Disease",
    "gastroenterology": "Gastroenterology",
    "hormonal": "Hormonal",
    "transplant": "Transplant",
    "other": "Other",
}


def _generate_category_summaries(
    drug_cards: List[Dict[str, Any]],
) -> Dict[str, Dict[str, Any]]:
    """Generate per-category drug assessment summaries."""
    summaries: Dict[str, Dict[str, Any]] = {}

    for cat_key, cat_label in CATEGORY_LABELS.items():
        cat_cards = [c for c in drug_cards if c["category"] == cat_key]
        if not cat_cards:
            continue

        alerts = sum(1 for c in cat_cards if c.get("critical_alert"))
        yellow_count = sum(1 for c in cat_cards if c["traffic_light"] == "yellow")
        red_count = sum(1 for c in cat_cards if c["traffic_light"] == "red")
        green_count = sum(1 for c in cat_cards if c["traffic_light"] == "green")

        if red_count > 0:
            summary = (
                f"{red_count} medication(s) in this category have significant "
                f"pharmacogenomic implications. Review the drug cards below for details."
            )
        elif yellow_count > 0:
            summary = (
                f"{yellow_count} medication(s) in this category may need dose "
                f"adjustments based on your genetic profile."
            )
        else:
            summary = (
                f"All {green_count} medication(s) in this category are expected "
                f"to work normally with standard dosing based on your profile."
            )

        summaries[cat_key] = {
            "label": cat_label,
            "drugs_analyzed": len(cat_cards),
            "alerts": alerts,
            "green": green_count,
            "yellow": yellow_count,
            "red": red_count,
            "summary": summary,
        }

    return summaries


# ─────────────────────────────────────────────────────────────────────────────
# PHARMACOGENOMIC PASSPORT
# ─────────────────────────────────────────────────────────────────────────────

# Key drugs per gene for the passport summary
GENE_KEY_DRUGS: Dict[str, str] = {
    "CYP2C19": "Clopidogrel, SSRIs, PPIs",
    "CYP2D6": "Codeine, Tamoxifen, SSRIs, TCAs",
    "CYP2C9": "Warfarin, NSAIDs",
    "CYP1A2": "Caffeine, Clozapine",
    "CYP3A5": "Tacrolimus",
    "CYP2B6": "Efavirenz",
    "VKORC1": "Warfarin",
    "DPYD": "5-FU, Capecitabine",
    "TPMT": "Azathioprine, 6-MP",
    "NUDT15": "Azathioprine, 6-MP",
    "UGT1A1": "Irinotecan, Atazanavir",
    "SLCO1B1": "Simvastatin, Statins",
    "ABCB1": "Multiple drugs (P-glycoprotein)",
    "HLA-B": "Abacavir",
    "OPRM1": "Opioids",
    "COMT": "Pain medications, Catecholamine drugs",
    "MTHFR": "Methotrexate, Folate",
    "F5": "Oral contraceptives, Estrogen",
    "IFNL3": "HCV treatment",
    "IFNL4": "HCV treatment",
    "NAT2": "Isoniazid",
}


def _generate_passport(
    enzyme_panel: Dict[str, Dict[str, Any]],
) -> Dict[str, Any]:
    """Generate a pharmacogenomic passport card for healthcare providers."""
    genes_list = []
    for gene in sorted(enzyme_panel.keys()):
        panel = enzyme_panel[gene]
        genes_list.append({
            "gene": gene,
            "full_name": GENE_FULL_NAMES.get(gene, gene),
            "status": panel["label"],
            "key_drugs": GENE_KEY_DRUGS.get(gene, ""),
        })

    return {
        "description": (
            "A summary card you can share with your healthcare provider. "
            "Print or save this section for your medical records."
        ),
        "genes": genes_list,
    }


# ─────────────────────────────────────────────────────────────────────────────
# MAIN ANALYSIS FUNCTION
# ─────────────────────────────────────────────────────────────────────────────

def analyze_precision_medicine(
    variants: Dict[str, Tuple[str, str, str]],
) -> Dict[str, Any]:
    """
    Analyze genetic variants for pharmacogenomic implications.

    Args:
        variants: Dictionary mapping rsID (lowercase) to
                  (chromosome, position, genotype) tuple.

    Returns:
        Complete precision medicine report as a JSON-serializable dict.
    """
    generated_at = datetime.now(timezone.utc).isoformat()

    # ── Step 1: Scan all pharmacogenomic SNPs ──
    all_found_variants: Dict[str, Dict[str, Any]] = {}
    all_missing_rsids: List[str] = []
    total_variants_found = 0

    for rsid, snp_info in PHARMA_DATABASE.items():
        genotype = _lookup_snp(variants, rsid)
        if genotype is not None:
            total_variants_found += 1
            risk_allele = snp_info["risk_allele"]

            # Handle special alleles (DEL, repeat counts, multi-char)
            if len(risk_allele) <= 1:
                risk_count = _count_allele(genotype, risk_allele)
            elif risk_allele == "DEL":
                # Deletion allele - check for missing allele indicator
                risk_count = genotype.count("-") if "-" in genotype else 0
            else:
                # Multi-character alleles (e.g., "TG" for IFNL4)
                risk_count = genotype.count(risk_allele)

            # Determine zygosity
            if risk_count >= 2:
                zygosity = "homozygous"
            elif risk_count == 1:
                zygosity = "heterozygous"
            else:
                zygosity = "wild_type"

            all_found_variants[rsid] = {
                "genotype": genotype,
                "risk_allele_count": risk_count,
                "zygosity": zygosity,
                "snp_info": snp_info,
            }
        else:
            all_missing_rsids.append(f"{rsid} ({snp_info['gene']} {snp_info['star']})")

    # ── Step 2: Build enzyme panel (per-gene classification) ──
    enzyme_panel: Dict[str, Dict[str, Any]] = {}

    for gene, snp_list in GENE_SNPS.items():
        gene_found_variants: Dict[str, Dict[str, Any]] = {}
        gene_variants_detail: List[Dict[str, Any]] = []
        gene_not_found: List[str] = []

        for rsid in snp_list:
            if rsid in all_found_variants:
                gene_found_variants[rsid] = all_found_variants[rsid]
                vdata = all_found_variants[rsid]
                gene_variants_detail.append({
                    "rsid": rsid,
                    "genotype": vdata["genotype"],
                    "star": vdata["snp_info"]["star"],
                    "function": vdata["snp_info"]["function"],
                    "zygosity": vdata["zygosity"],
                    "risk_allele_count": vdata["risk_allele_count"],
                })
            else:
                gene_not_found.append(rsid)

        # Classify even if no variants found (defaults to normal)
        classification = _classify_metabolizer(gene, gene_found_variants)

        # Count drugs affected
        drugs_affected = sum(
            1 for d in DRUG_DATABASE.values() if gene in d["genes"]
        )

        enzyme_panel[gene] = {
            "gene": gene,
            "full_name": GENE_FULL_NAMES.get(gene, gene),
            "status": classification["status"],
            "label": classification["label"],
            "clinical_significance": classification["clinical_significance"],
            "variants_found": gene_variants_detail,
            "variants_not_found": gene_not_found,
            "drugs_affected_count": drugs_affected,
            "population_frequency": POPULATION_FREQUENCY.get(gene, ""),
            "has_data": len(gene_found_variants) > 0,
        }

    # ── Step 3: Assess all drugs ──
    drug_cards: List[Dict[str, Any]] = []
    for drug_key, drug_info in DRUG_DATABASE.items():
        card = _assess_drug(drug_key, drug_info, enzyme_panel)
        drug_cards.append(card)

    # Sort drug cards: red first, then yellow, then green
    light_order = {"red": 0, "yellow": 1, "green": 2}
    drug_cards.sort(key=lambda c: light_order.get(c["traffic_light"], 3))

    # ── Step 4: Generate critical alerts ──
    critical_alerts = _generate_critical_alerts(drug_cards, enzyme_panel)

    # ── Step 5: Generate interaction risks ──
    interaction_risks = _generate_interaction_risks(drug_cards, enzyme_panel)

    # ── Step 6: Generate biomarker plan ──
    biomarker_plan = _generate_biomarker_plan(enzyme_panel)

    # ── Step 7: Category summaries ──
    category_summaries = _generate_category_summaries(drug_cards)

    # ── Step 8: Pharmacogenomic passport ──
    passport = _generate_passport(enzyme_panel)

    # ── Step 9: Summary ──
    genes_analyzed = len(enzyme_panel)
    non_normal_genes = [
        g for g, p in enzyme_panel.items() if p["status"] != "normal" and p["has_data"]
    ]
    drugs_with_alerts = sum(1 for c in drug_cards if c["traffic_light"] != "green")

    if non_normal_genes:
        affected_drug_count = sum(
            enzyme_panel[g]["drugs_affected_count"] for g in non_normal_genes
        )
        # Deduplicate (some drugs affected by multiple genes)
        overview = (
            f"You have {len(non_normal_genes)} gene(s) with non-standard metabolism "
            f"({', '.join(non_normal_genes)}) that may affect your response to "
            f"multiple medications. Review the drug cards below for specific guidance."
        )
    else:
        overview = (
            "Based on the variants detected, your pharmacogenomic profile appears "
            "largely normal. Standard drug dosing is expected to be appropriate for "
            "most medications. Always consult your healthcare provider for personalized advice."
        )

    summary = {
        "total_genes_analyzed": genes_analyzed,
        "total_variants_found": total_variants_found,
        "total_variants_checked": len(PHARMA_DATABASE),
        "total_drugs_assessed": len(drug_cards),
        "high_risk_alerts": len(critical_alerts),
        "drugs_needing_attention": drugs_with_alerts,
        "non_normal_genes": len(non_normal_genes),
        "metabolizer_overview": overview,
    }

    # ── Build final report ──
    report = {
        "report_type": "precision_medicine",
        "version": "1.0",
        "generated_at": generated_at,
        "disclaimer": (
            "IMPORTANT: This pharmacogenomics report is for INFORMATIONAL and EDUCATIONAL "
            "purposes only. It is NOT a substitute for professional medical advice, diagnosis, "
            "or treatment. Genetic information from consumer DNA tests covers only a subset of "
            "known pharmacogenomic variants and CANNOT replace clinical-grade pharmacogenomic "
            "testing. Never change, start, or stop any medication based solely on this report. "
            "Always consult your healthcare provider and a qualified pharmacist or clinical "
            "pharmacogeneticist before making any medication decisions. Drug response is "
            "influenced by many factors beyond genetics, including age, weight, organ function, "
            "other medications, and overall health status."
        ),
        "summary": summary,
        "enzyme_panel": [
            enzyme_panel[gene] for gene in sorted(enzyme_panel.keys())
        ],
        "drug_cards": drug_cards,
        "critical_alerts": critical_alerts,
        "drug_interaction_risks": interaction_risks,
        "biomarker_plan": biomarker_plan,
        "category_summaries": category_summaries,
        "pharmacogenomic_passport": passport,
    }

    return report


# ─────────────────────────────────────────────────────────────────────────────
# JSON WRAPPER
# ─────────────────────────────────────────────────────────────────────────────

def generate_precision_medicine_json(
    variants: Dict[str, Tuple[str, str, str]],
) -> str:
    """
    Analyze pharmacogenomic variants and return a JSON string.

    This is the main entry point for the analysis pipeline.

    Args:
        variants: Dictionary mapping rsID (lowercase) to
                  (chromosome, position, genotype) tuple.

    Returns:
        JSON string of the complete precision medicine report.
    """
    report = analyze_precision_medicine(variants)
    return json.dumps(report, indent=2, ensure_ascii=False)
