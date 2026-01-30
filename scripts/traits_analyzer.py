"""
Traits Analyzer
Analyzes genetic variants for traits, wellness, and personal characteristics
"""

from typing import Dict, List, Tuple, Optional
from dataclasses import dataclass

from .traits_snp_database import TRAITS_DATABASE, TraitSNP, TRAIT_CATEGORIES


@dataclass
class TraitResult:
    rsid: str
    gene: str
    category: str
    trait: str
    genotype: str
    has_risk_allele: bool
    risk_allele_count: int  # 0, 1, or 2
    effect: str
    interpretation: str
    description: str


@dataclass
class TraitsAnalysisResult:
    total_traits_checked: int
    traits_found: int
    traits_not_found: int
    results_by_category: Dict[str, List[TraitResult]]
    missing_by_category: Dict[str, List[str]]  # SNPs not in genome


def analyze_traits(
    variants: Dict[str, Tuple[str, str, str]]  # rsid -> (chromosome, position, genotype)
) -> TraitsAnalysisResult:
    """
    Analyze genetic variants against the traits database

    Args:
        variants: Dictionary mapping rsID to (chromosome, position, genotype)

    Returns:
        TraitsAnalysisResult with categorized findings
    """
    results_by_category: Dict[str, List[TraitResult]] = {cat: [] for cat in TRAIT_CATEGORIES}
    missing_by_category: Dict[str, List[str]] = {cat: [] for cat in TRAIT_CATEGORIES}

    traits_found = 0
    traits_not_found = 0

    # Check each trait SNP
    for rsid, trait_info in TRAITS_DATABASE.items():
        rsid_lower = rsid.lower()

        if rsid_lower in variants:
            traits_found += 1
            chrom, pos, genotype = variants[rsid_lower]
            genotype = genotype.upper().replace('-', '')

            # Analyze the genotype
            result = analyze_single_trait(rsid_lower, genotype, trait_info)
            results_by_category[trait_info.category].append(result)
        else:
            traits_not_found += 1
            missing_by_category[trait_info.category].append(
                f"{rsid} ({trait_info.gene}): {trait_info.trait}"
            )

    return TraitsAnalysisResult(
        total_traits_checked=len(TRAITS_DATABASE),
        traits_found=traits_found,
        traits_not_found=traits_not_found,
        results_by_category=results_by_category,
        missing_by_category=missing_by_category
    )


def analyze_single_trait(rsid: str, genotype: str, trait_info: TraitSNP) -> TraitResult:
    """Analyze a single trait SNP"""
    risk_allele = trait_info.risk_allele.upper()
    risk_count = genotype.count(risk_allele)
    has_risk = risk_count > 0

    # Generate interpretation based on genotype
    if risk_count == 0:
        interpretation = f"You do not carry the {risk_allele} allele associated with this trait."
    elif risk_count == 1:
        interpretation = f"You carry one copy of the {risk_allele} allele (heterozygous). {trait_info.effect}"
    else:
        interpretation = f"You carry two copies of the {risk_allele} allele (homozygous). {trait_info.effect}"

    return TraitResult(
        rsid=rsid.upper() if not rsid.startswith('rs') else rsid,
        gene=trait_info.gene,
        category=trait_info.category,
        trait=trait_info.trait,
        genotype=genotype,
        has_risk_allele=has_risk,
        risk_allele_count=risk_count,
        effect=trait_info.effect if has_risk else "Typical/common variant",
        interpretation=interpretation,
        description=trait_info.description
    )


def generate_traits_report(result: TraitsAnalysisResult, subject_name: str = "Subject") -> str:
    """Generate a comprehensive traits report"""
    from datetime import datetime

    report = f"""# Personal Traits & Wellness Report

**Subject:** {subject_name}
**Generated:** {datetime.now().strftime("%B %d, %Y")}

---

## Overview

This report analyzes your genetic variants related to personal traits, wellness, and characteristics beyond disease risk. These insights are based on peer-reviewed research but should be considered informational - genetics is just one factor influencing these traits.

| Metric | Value |
|--------|-------|
| Total Trait SNPs Checked | {result.total_traits_checked} |
| SNPs Found in Your Data | {result.traits_found} |
| SNPs Not Available | {result.traits_not_found} |

---

"""

    # Generate sections for each category with findings
    for category in TRAIT_CATEGORIES:
        findings = result.results_by_category.get(category, [])
        if not findings:
            continue

        report += f"## {get_category_emoji(category)} {category}\n\n"

        # Group by trait
        traits_dict: Dict[str, List[TraitResult]] = {}
        for finding in findings:
            if finding.trait not in traits_dict:
                traits_dict[finding.trait] = []
            traits_dict[finding.trait].append(finding)

        for trait_name, trait_findings in traits_dict.items():
            report += f"### {trait_name}\n\n"

            for finding in trait_findings:
                status_icon = "🔵" if finding.has_risk_allele else "⚪"
                report += f"{status_icon} **{finding.gene}** ({finding.rsid})\n\n"
                report += f"- **Genotype:** {finding.genotype}\n"
                report += f"- **Interpretation:** {finding.interpretation}\n"
                report += f"- **Background:** {finding.description}\n\n"

        report += "---\n\n"

    # Section for missing data
    has_missing = any(missing for missing in result.missing_by_category.values())
    if has_missing:
        report += """## Data Not Available

The following traits could not be analyzed because the SNPs were not found in your genome file. This is normal - different DNA testing services test different SNPs.

"""
        for category, missing in result.missing_by_category.items():
            if missing:
                report += f"### {category}\n"
                for item in missing[:5]:  # Limit to 5 per category
                    report += f"- {item}\n"
                if len(missing) > 5:
                    report += f"- *...and {len(missing) - 5} more*\n"
                report += "\n"

    report += """---

## Important Notes

1. **Genetics is not destiny** - These are tendencies and probabilities, not certainties. Environment, lifestyle, and other genetic factors also play major roles.

2. **Research is ongoing** - Scientific understanding of these associations continues to evolve. Some findings may be refined or revised.

3. **Individual variation** - Even with "risk" variants, many people do not exhibit the associated traits, and vice versa.

4. **Consult professionals** - For mental health, cognitive, or significant wellness concerns, consult qualified healthcare providers.

---
*Generated by GeneHealth Analysis Platform*
"""

    return report


def get_category_emoji(category: str) -> str:
    """Return an emoji for each category"""
    emojis = {
        "Cognitive": "🧠",
        "Metabolism": "⚡",
        "Sleep": "😴",
        "Physical": "👤",
        "Athletic": "🏃",
        "Mental Health": "💭",
        "Longevity": "⏳",
        "Sensitivity": "🎯"
    }
    return emojis.get(category, "📊")
