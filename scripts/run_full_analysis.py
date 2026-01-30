"""
Full Analysis Pipeline
Main orchestrator for genetic analysis
"""

import os
import re
from typing import Dict, Tuple, Any

from .disease_risk_analyzer import analyze_disease_risks
from .generate_exhaustive_report import generate_reports
from .traits_analyzer import analyze_traits, generate_traits_report


def parse_genome_file(content: str, source_format: str) -> Dict[str, Tuple[str, str, str]]:
    """
    Parse genome file content into variant dictionary

    Args:
        content: Raw genome file content
        source_format: Format identifier (23andme, ancestry, etc.)

    Returns:
        Dictionary mapping rsID to (chromosome, position, genotype)
    """
    variants: Dict[str, Tuple[str, str, str]] = {}
    lines = content.strip().split('\n')

    for line in lines:
        # Skip comments and empty lines
        if line.startswith('#') or not line.strip():
            continue

        # Skip header lines
        if line.lower().startswith('rsid') or line.lower().startswith('"rsid'):
            continue

        try:
            parts = parse_line(line, source_format)
            if parts:
                rsid, chrom, pos, genotype = parts
                # Normalize rsID
                if rsid.lower().startswith('rs'):
                    variants[rsid.lower()] = (chrom, pos, genotype)
        except Exception:
            continue  # Skip malformed lines

    return variants


def parse_line(line: str, source_format: str) -> Tuple[str, str, str, str] | None:
    """Parse a single line based on format"""

    if source_format == '23andme':
        parts = line.split('\t')
        if len(parts) >= 4:
            return (parts[0], parts[1], parts[2], parts[3])

    elif source_format == 'ancestry':
        parts = line.split('\t')
        if len(parts) >= 5:
            # Ancestry has split alleles
            genotype = parts[3] + parts[4]
            genotype = genotype.replace('0', '-')
            return (parts[0], parts[1], parts[2], genotype)

    elif source_format == 'myheritage':
        # Remove quotes and split by comma
        line = line.replace('"', '')
        parts = line.split(',')
        if len(parts) >= 4:
            return (parts[0], parts[1], parts[2], parts[3])

    elif source_format == 'ftdna':
        parts = line.split(',')
        if len(parts) >= 4:
            return (parts[0], parts[1], parts[2], parts[3])

    elif source_format == 'genera':
        # Can be tab or comma separated
        if '\t' in line:
            parts = line.split('\t')
        else:
            parts = line.split(',')
        if len(parts) >= 4:
            return (parts[0], parts[1], parts[2], parts[3])

    elif source_format == 'nebula':
        # VCF format
        parts = line.split('\t')
        if len(parts) >= 10 and not line.startswith('#'):
            chrom = parts[0].replace('chr', '')
            pos = parts[1]
            rsid = parts[2] if parts[2] != '.' else f"chr{chrom}:{pos}"
            ref = parts[3]
            alt = parts[4]

            # Parse genotype from sample column
            format_fields = parts[8].split(':')
            sample_fields = parts[9].split(':')

            gt_index = format_fields.index('GT') if 'GT' in format_fields else 0
            gt = sample_fields[gt_index]

            # Convert genotype
            if gt in ('0/0', '0|0'):
                genotype = ref + ref
            elif gt in ('0/1', '0|1', '1/0', '1|0'):
                genotype = ref + alt
            elif gt in ('1/1', '1|1'):
                genotype = alt + alt
            else:
                genotype = '--'

            return (rsid, chrom, pos, genotype)

    else:
        # Generic format - try common patterns
        # Tab-separated
        if '\t' in line:
            parts = line.split('\t')
            if len(parts) >= 4:
                return (parts[0], parts[1], parts[2], parts[3])

        # Comma-separated
        parts = line.split(',')
        if len(parts) >= 4:
            return (parts[0], parts[1], parts[2], parts[3])

    return None


def run_analysis(genome_content: str, source_format: str) -> Dict[str, Any]:
    """
    Run the complete analysis pipeline

    Args:
        genome_content: Raw genome file content
        source_format: Format identifier

    Returns:
        Dictionary with analysis results
    """
    # Parse genome file
    variants = parse_genome_file(genome_content, source_format)

    if not variants:
        raise ValueError("No valid variants found in genome file")

    # Run disease risk analysis
    clinvar_path = os.environ.get('CLINVAR_PATH', '/app/data/clinvar_alleles.tsv')
    if not os.path.exists(clinvar_path):
        clinvar_path = None

    risk_result = analyze_disease_risks(variants, clinvar_path)

    # Run traits analysis
    traits_result = analyze_traits(variants)
    print(f"Traits analysis: {traits_result.traits_found} found, {traits_result.traits_not_found} not available")

    # Generate findings summary
    findings_summary = {
        "totalSnpsAnalyzed": risk_result.total_variants_analyzed,
        "clinvarMatches": risk_result.clinvar_matches,
        "highRiskVariants": len(risk_result.high_risk_variants),
        "moderateRiskVariants": len(risk_result.moderate_risk_variants),
        "pharmacogenomicVariants": len(risk_result.pharmacogenomic_variants),
        "traitsAnalyzed": traits_result.traits_found,
        "categories": [
            {"category": cat, "count": count, "highRisk": 0}
            for cat, count in risk_result.categories.items()
        ]
    }

    # Generate reports (disease-focused)
    reports = generate_reports(risk_result)

    # Add traits report
    reports["traits"] = generate_traits_report(traits_result)

    return {
        "snp_count": len(variants),
        "findings_summary": findings_summary,
        "reports": reports,
        "risk_result": risk_result,
        "traits_result": traits_result
    }


if __name__ == "__main__":
    # Test with sample data
    sample_data = """# Sample 23andMe data
rsid	chromosome	position	genotype
rs1801133	1	11856378	CT
rs429358	19	44908684	TC
rs4244285	10	94781859	GA
"""
    result = run_analysis(sample_data, "23andme")
    print(f"Analyzed {result['snp_count']} variants")
    print(f"High risk: {result['findings_summary']['highRiskVariants']}")
