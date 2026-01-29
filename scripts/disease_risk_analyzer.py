"""
Disease Risk Analyzer
Analyzes genetic variants against ClinVar and curated databases
"""

import os
import csv
from typing import Dict, List, Optional, Tuple
from dataclasses import dataclass

from .comprehensive_snp_database import SNP_DATABASE, SNPInfo


@dataclass
class VariantMatch:
    rsid: str
    genotype: str
    gene: str
    condition: str
    clinical_significance: str
    risk_level: str
    description: str
    recommendations: List[str]
    category: str


@dataclass
class DiseaseRiskResult:
    total_variants_analyzed: int
    clinvar_matches: int
    high_risk_variants: List[VariantMatch]
    moderate_risk_variants: List[VariantMatch]
    low_risk_variants: List[VariantMatch]
    beneficial_variants: List[VariantMatch]
    pharmacogenomic_variants: List[VariantMatch]
    categories: Dict[str, int]


class DiseaseRiskAnalyzer:
    """Analyzes genetic variants for disease risk"""

    def __init__(self, clinvar_path: Optional[str] = None):
        self.clinvar_data: Dict[str, dict] = {}
        self.snp_database = SNP_DATABASE

        if clinvar_path and os.path.exists(clinvar_path):
            self._load_clinvar(clinvar_path)

    def _load_clinvar(self, path: str):
        """Load ClinVar data from TSV file"""
        try:
            with open(path, 'r', encoding='utf-8') as f:
                reader = csv.DictReader(f, delimiter='\t')
                for row in reader:
                    rsid = row.get('RS# (dbSNP)', '').strip()
                    if rsid:
                        self.clinvar_data[f"rs{rsid}"] = {
                            'gene': row.get('GeneSymbol', ''),
                            'condition': row.get('PhenotypeList', ''),
                            'clinical_significance': row.get('ClinicalSignificance', ''),
                            'review_status': row.get('ReviewStatus', ''),
                            'chromosome': row.get('Chromosome', ''),
                            'position': row.get('PositionVCF', ''),
                        }
        except Exception as e:
            print(f"Warning: Could not load ClinVar data: {e}")

    def analyze_variants(
        self,
        variants: Dict[str, Tuple[str, str, str]]  # rsid -> (chromosome, position, genotype)
    ) -> DiseaseRiskResult:
        """
        Analyze a set of genetic variants

        Args:
            variants: Dictionary mapping rsID to (chromosome, position, genotype)

        Returns:
            DiseaseRiskResult with categorized findings
        """
        high_risk = []
        moderate_risk = []
        low_risk = []
        beneficial = []
        pharmacogenomic = []
        categories: Dict[str, int] = {}
        clinvar_matches = 0

        for rsid, (chrom, pos, genotype) in variants.items():
            rsid_lower = rsid.lower()

            # Check curated SNP database first
            if rsid_lower in self.snp_database:
                snp_info = self.snp_database[rsid_lower]
                match = self._evaluate_snp_match(rsid_lower, genotype, snp_info)

                if match:
                    # Categorize by risk level
                    if match.risk_level == "high":
                        high_risk.append(match)
                    elif match.risk_level == "moderate":
                        moderate_risk.append(match)
                    elif match.risk_level == "low":
                        low_risk.append(match)
                    elif match.risk_level == "beneficial":
                        beneficial.append(match)

                    # Track pharmacogenomic variants separately
                    if match.category == "Drug Metabolism":
                        pharmacogenomic.append(match)

                    # Count by category
                    categories[match.category] = categories.get(match.category, 0) + 1

            # Check ClinVar data
            if rsid_lower in self.clinvar_data:
                clinvar_matches += 1
                clinvar_info = self.clinvar_data[rsid_lower]

                # Only add if not already found in curated database
                if rsid_lower not in self.snp_database:
                    match = self._create_clinvar_match(rsid_lower, genotype, clinvar_info)
                    if match:
                        if "pathogenic" in match.clinical_significance.lower():
                            high_risk.append(match)
                        elif "likely pathogenic" in match.clinical_significance.lower():
                            moderate_risk.append(match)
                        categories["ClinVar"] = categories.get("ClinVar", 0) + 1

        return DiseaseRiskResult(
            total_variants_analyzed=len(variants),
            clinvar_matches=clinvar_matches,
            high_risk_variants=high_risk,
            moderate_risk_variants=moderate_risk,
            low_risk_variants=low_risk,
            beneficial_variants=beneficial,
            pharmacogenomic_variants=pharmacogenomic,
            categories=categories
        )

    def _evaluate_snp_match(
        self,
        rsid: str,
        genotype: str,
        snp_info: SNPInfo
    ) -> Optional[VariantMatch]:
        """Evaluate if genotype matches risk allele"""
        genotype = genotype.upper().replace('-', '')

        # Count risk alleles
        risk_count = genotype.count(snp_info.risk_allele.upper())

        if risk_count == 0:
            return None  # No risk allele present

        # Determine risk level based on zygosity
        if risk_count == 2:  # Homozygous risk
            risk_level = snp_info.significance
            description = f"Homozygous for risk variant. {snp_info.description}"
        else:  # Heterozygous
            # Generally lower risk for heterozygous
            risk_level = "moderate" if snp_info.significance == "high" else "low"
            description = f"Heterozygous carrier. {snp_info.description}"

        return VariantMatch(
            rsid=rsid,
            genotype=genotype,
            gene=snp_info.gene,
            condition=snp_info.condition,
            clinical_significance=snp_info.significance,
            risk_level=risk_level,
            description=description,
            recommendations=snp_info.recommendations,
            category=snp_info.category
        )

    def _create_clinvar_match(
        self,
        rsid: str,
        genotype: str,
        clinvar_info: dict
    ) -> Optional[VariantMatch]:
        """Create a match from ClinVar data"""
        clinical_sig = clinvar_info.get('clinical_significance', '')

        # Only include pathogenic or likely pathogenic
        if not any(term in clinical_sig.lower() for term in ['pathogenic', 'risk factor']):
            return None

        return VariantMatch(
            rsid=rsid,
            genotype=genotype,
            gene=clinvar_info.get('gene', 'Unknown'),
            condition=clinvar_info.get('condition', 'Unknown condition'),
            clinical_significance=clinical_sig,
            risk_level="high" if "pathogenic" in clinical_sig.lower() else "moderate",
            description=f"ClinVar classification: {clinical_sig}",
            recommendations=["Consult with a genetic counselor for interpretation"],
            category="ClinVar"
        )


def analyze_disease_risks(
    variants: Dict[str, Tuple[str, str, str]],
    clinvar_path: Optional[str] = None
) -> DiseaseRiskResult:
    """
    Convenience function to analyze disease risks

    Args:
        variants: Dictionary mapping rsID to (chromosome, position, genotype)
        clinvar_path: Optional path to ClinVar data file

    Returns:
        DiseaseRiskResult with categorized findings
    """
    analyzer = DiseaseRiskAnalyzer(clinvar_path)
    return analyzer.analyze_variants(variants)
