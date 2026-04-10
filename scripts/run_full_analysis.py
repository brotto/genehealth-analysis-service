"""
Full Analysis Pipeline
Main orchestrator for genetic analysis
"""

import os
import re
import json
from typing import Dict, Tuple, Any

from .disease_risk_analyzer import analyze_disease_risks
from .generate_exhaustive_report import generate_reports
from .traits_analyzer import analyze_traits, generate_traits_report, generate_traits_json
from .neanderthal_analyzer import analyze_neanderthal
from .ancestry_analyzer import analyze_ancestry
from .fitness_analyzer import analyze_fitness, generate_fitness_json
from .nutrition_analyzer import analyze_nutrition, generate_nutrition_json
from .skincare_analyzer import analyze_skincare, generate_skincare_json
from .ancient_analyzer import analyze_ancient, generate_ancient_json
from .viking_analyzer import analyze_viking, generate_viking_json
from .italian_analyzer import analyze_italian, generate_italian_json
from .evolution_analyzer import analyze_evolution, generate_evolution_json
from .haplogroup_analyzer import analyze_haplogroup, generate_haplogroup_json
from .roh_analyzer import analyze_roh, generate_roh_json
from .immune_microbiome_analyzer import analyze_immune_microbiome, generate_immune_microbiome_json
from .historical_connections_analyzer import analyze_historical_connections, generate_historical_connections_json
from .mind_spirit_analyzer import analyze_mind_spirit, generate_mind_spirit_json
from .dream_sleep_analyzer import analyze_dream_sleep, generate_dream_sleep_json
from .precision_medicine_analyzer import analyze_precision_medicine, generate_precision_medicine_json
from .longevity_aging_analyzer import analyze_longevity_aging, generate_longevity_aging_json
from .steroid_pharmacogenomics_analyzer import analyze_steroid_pharmacogenomics, generate_steroid_pharmacogenomics_json
from .mental_wellbeing_analyzer import analyze_mental_wellbeing, generate_mental_wellbeing_json
from .big_five_personality_analyzer import analyze_big_five, generate_big_five_json
from .organ_blood_wellness_analyzer import analyze_organ_blood_wellness, generate_organ_blood_json
from .intelligence_analyzer import analyze_intelligence, generate_intelligence_json
from .sports_analyzer import analyze_sports, generate_sports_json
from .fun_traits_analyzer import analyze_fun_traits, generate_fun_traits_json
from .celtic_analyzer import analyze_celtic, generate_celtic_json
from .british_isles_analyzer import analyze_british_isles, generate_british_isles_json
from .native_american_analyzer import analyze_native_american, generate_native_american_json
from .latino_breakdown_analyzer import analyze_latino_breakdown, generate_latino_breakdown_json
from .african_breakdown_analyzer import analyze_african_breakdown, generate_african_breakdown_json
from .european_breakdown_analyzer import analyze_european_breakdown, generate_european_breakdown_json
from .hunter_gatherer_farmer_analyzer import analyze_hunter_gatherer_farmer, generate_hunter_gatherer_farmer_json
from .deep_ancestry_analyzer import analyze_deep_ancestry, generate_deep_ancestry_json


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

    # Add traits report (JSON for visual display)
    reports["traits"] = json.dumps(generate_traits_json(traits_result))

    # Run neanderthal analysis
    try:
        neanderthal_result = analyze_neanderthal(variants)
        reports["neanderthal"] = json.dumps(neanderthal_result)
        print(f"Neanderthal analysis: {neanderthal_result['summary']['variantsCarried']} variants carried")
        findings_summary["neanderthalVariantsCarried"] = neanderthal_result["summary"]["variantsCarried"]
        findings_summary["neanderthalPercentage"] = neanderthal_result["summary"]["estimatedPercentage"]
    except Exception as e:
        print(f"Neanderthal analysis failed: {e}")

    # Run ancestry analysis
    try:
        ancestry_result = analyze_ancestry(variants)
        reports["ancestry"] = json.dumps(ancestry_result)
        print(f"Ancestry analysis: {ancestry_result['summary']['aimsGenotyped']} AIMs genotyped")
        findings_summary["ancestryReliability"] = ancestry_result["summary"]["reliability"]
    except Exception as e:
        print(f"Ancestry analysis failed: {e}")

    # Run fitness analysis
    try:
        fitness_result = analyze_fitness(variants)
        reports["fitness"] = json.dumps(generate_fitness_json(fitness_result))
        print(f"Fitness analysis: {fitness_result.found} variants found")
    except Exception as e:
        print(f"Fitness analysis failed: {e}")

    # Run nutrition analysis
    try:
        nutrition_result = analyze_nutrition(variants)
        reports["nutrition"] = json.dumps(generate_nutrition_json(nutrition_result))
        print(f"Nutrition analysis: {nutrition_result['found']} variants found")
    except Exception as e:
        print(f"Nutrition analysis failed: {e}")

    # Run skincare analysis
    try:
        skincare_result = analyze_skincare(variants)
        reports["skincare"] = json.dumps(generate_skincare_json(skincare_result))
        print(f"Skincare analysis: {skincare_result['found']} variants found")
    except Exception as e:
        print(f"Skincare analysis failed: {e}")

    # Run ancient bloodlines analysis (includes Viking + Italian sub-analyses)
    try:
        ancient_result = analyze_ancient(variants)
        ancient_json = generate_ancient_json(ancient_result)
        print(f"Ancient bloodlines analysis: {len(ancient_result['used_snps'])} AIMs genotyped")

        # Add Viking heritage sub-analysis
        try:
            viking_result = analyze_viking(variants)
            ancient_json["vikingHeritage"] = generate_viking_json(viking_result)
            print(f"Viking heritage analysis: {len(viking_result['used_snps'])} SNPs used")
        except Exception as e:
            print(f"Viking heritage sub-analysis failed: {e}")

        # Add Italian ancestry sub-analysis
        try:
            italian_result = analyze_italian(variants)
            ancient_json["italianAncestry"] = generate_italian_json(italian_result)
            print(f"Italian ancestry analysis: {len(italian_result['used_snps'])} SNPs used")
        except Exception as e:
            print(f"Italian ancestry sub-analysis failed: {e}")

        reports["ancient"] = json.dumps(ancient_json)
    except Exception as e:
        print(f"Ancient bloodlines analysis failed: {e}")

    # Run evolution analysis (adaptive sweeps + archaic introgression)
    try:
        evolution_result = analyze_evolution(variants)
        reports["evolution"] = json.dumps(generate_evolution_json(evolution_result))
        print(f"Evolution analysis: {evolution_result['snps_used']} SNPs used")
    except Exception as e:
        print(f"Evolution analysis failed: {e}")

    # Run haplogroup analysis (mitochondrial + Y chromosome)
    try:
        haplogroup_result = analyze_haplogroup(variants)
        reports["haplogroup"] = json.dumps(generate_haplogroup_json(haplogroup_result))
        mt = haplogroup_result.get('mitochondrial', {}).get('haplogroup', 'Unknown')
        y = haplogroup_result.get('y_chromosome', {}).get('haplogroup', 'Unknown')
        print(f"Haplogroup analysis: mt={mt}, Y={y}")
    except Exception as e:
        print(f"Haplogroup analysis failed: {e}")

    # Run ROH analysis (runs of homozygosity)
    try:
        roh_result = analyze_roh(variants)
        reports["roh"] = json.dumps(generate_roh_json(roh_result))
        print(f"ROH analysis: FROH={roh_result['froh']:.4f}, {roh_result['total_segments']} segments")
    except Exception as e:
        print(f"ROH analysis failed: {e}")

    # Run immune & microbiome analysis
    try:
        immune_result = analyze_immune_microbiome(variants)
        reports["immune_microbiome"] = json.dumps(generate_immune_microbiome_json(immune_result))
        print(f"Immune/Microbiome analysis: {immune_result['immune_snps_typed']} immune SNPs, {immune_result['microbiome_snps_typed']} microbiome SNPs")
    except Exception as e:
        print(f"Immune/Microbiome analysis failed: {e}")

    # Run historical connections analysis (haplogroup matching + Otzi comparison)
    try:
        hist_result = analyze_historical_connections(variants)
        reports["historical_connections"] = json.dumps(generate_historical_connections_json(hist_result))
        y_hg = hist_result['user_y_hg'] or 'N/A'
        mt_hg = hist_result['user_mt_hg'] or 'N/A'
        y_matches = len(hist_result['y_connections'])
        mt_matches = len(hist_result['mt_connections'])
        print(f"Historical connections: Y={y_hg}, mt={mt_hg}, {y_matches} Y-DNA matches, {mt_matches} mtDNA matches")
    except Exception as e:
        print(f"Historical connections analysis failed: {e}")

    # Run mind & spirit analysis (personality, mental health, spiritual sensitivity)
    try:
        ms_result = analyze_mind_spirit(variants)
        reports["mind_spirit"] = json.dumps(generate_mind_spirit_json(ms_result))
        p_count = len([t for t in ms_result['personality_traits'] if t['genotype'] != 'Not available'])
        m_count = len([t for t in ms_result['mental_health_traits'] if t['genotype'] != 'Not available'])
        s_count = len([t for t in ms_result['spiritual_traits'] if t['genotype'] != 'Not available'])
        print(f"Mind & Spirit: {p_count} personality, {m_count} mental health, {s_count} spiritual SNPs found")
    except Exception as e:
        print(f"Mind & Spirit analysis failed: {e}")

    # Run dream & sleep architecture analysis
    print("Starting dream & sleep analysis...")
    try:
        ds_result = analyze_dream_sleep(variants)
        reports["dream_sleep"] = json.dumps(generate_dream_sleep_json(ds_result))
        dv_count = len([t for t in ds_result.get('dream_vividness', []) if t.get('found')])
        sa_count = len([t for t in ds_result.get('sleep_architecture', []) if t.get('found')])
        ch_count = len([t for t in ds_result.get('chronotype', []) if t.get('found')])
        np_count = len([t for t in ds_result.get('nocturnal_phenomena', []) if t.get('found')])
        print(f"Dream & Sleep: {dv_count} dream vividness, {sa_count} sleep architecture, {ch_count} chronotype, {np_count} nocturnal phenomena SNPs found")
    except Exception as e:
        print(f"Dream & Sleep analysis failed: {e}")

    # Run precision medicine / pharmacogenomics analysis
    print("Starting precision medicine analysis...")
    try:
        pm_result = analyze_precision_medicine(variants)
        reports["precision_medicine"] = json.dumps(pm_result)
        genes_found = len([e for e in pm_result.get('enzyme_panel', []) if e.get('variants_found')])
        alerts = pm_result.get('summary', {}).get('high_risk_alerts', 0)
        drugs = pm_result.get('summary', {}).get('total_drugs_assessed', 0)
        print(f"Precision Medicine: {genes_found} genes with variants, {drugs} drugs assessed, {alerts} critical alerts")
    except Exception as e:
        print(f"Precision Medicine analysis failed: {e}")

    # Run longevity & aging pathways analysis
    print("Starting longevity & aging analysis...")
    try:
        longevity_result = analyze_longevity_aging(variants)
        reports["longevity_aging"] = generate_longevity_aging_json(longevity_result)
        score = longevity_result.get('summary', {}).get('longevityScore', 0)
        found = longevity_result.get('summary', {}).get('totalVariantsFound', 0)
        apoe = longevity_result.get('summary', {}).get('apoeGenotype', 'Unknown')
        print(f"Longevity & Aging: score={score}, {found} variants found, APOE={apoe}")
    except Exception as e:
        print(f"Longevity & Aging analysis failed: {e}")

    # Run steroid pharmacogenomics analysis (hormonal response profile)
    print("Starting steroid pharmacogenomics analysis...")
    try:
        spg_result = analyze_steroid_pharmacogenomics(variants)
        reports["steroid_pharmacogenomics"] = json.dumps(generate_steroid_pharmacogenomics_json(spg_result))
        found = len(spg_result['markers_found'])
        score = spg_result['composite_scores']['overall_androgen_sensitivity']
        label = spg_result['composite_scores']['overall_label']
        print(f"Steroid Pharmacogenomics: {found} markers found, overall score={score} ({label})")
    except Exception as e:
        print(f"Steroid Pharmacogenomics analysis failed: {e}")

    # Run mental wellbeing analysis
    try:
        mw_result = analyze_mental_wellbeing(variants)
        reports["mental_wellbeing"] = json.dumps(generate_mental_wellbeing_json(mw_result))
        print(f"Mental Wellbeing: {mw_result.found} variants found")
    except Exception as e:
        print(f"Mental Wellbeing analysis failed: {e}")

    # Run Big Five personality analysis
    try:
        bf_result = analyze_big_five(variants)
        reports["big_five_personality"] = json.dumps(generate_big_five_json(bf_result))
        print(f"Big Five Personality: {bf_result.found} variants found")
    except Exception as e:
        print(f"Big Five Personality analysis failed: {e}")

    # Run Organ & Blood Wellness analysis
    try:
        ob_result = analyze_organ_blood_wellness(variants)
        reports["organ_blood_wellness"] = json.dumps(generate_organ_blood_json(ob_result))
        print(f"Organ & Blood Wellness: {ob_result.found} variants found")
    except Exception as e:
        print(f"Organ & Blood Wellness analysis failed: {e}")

    # Run Intelligence analysis
    try:
        iq_result = analyze_intelligence(variants)
        reports["intelligence"] = json.dumps(generate_intelligence_json(iq_result))
        print(f"Intelligence: {iq_result.found} variants found")
    except Exception as e:
        print(f"Intelligence analysis failed: {e}")

    # Run Sports analysis
    try:
        sp_result = analyze_sports(variants)
        reports["sports"] = json.dumps(generate_sports_json(sp_result))
        print(f"Sports: {sp_result.found} variants found")
    except Exception as e:
        print(f"Sports analysis failed: {e}")

    # Run Fun Traits analysis
    try:
        ft_result = analyze_fun_traits(variants)
        reports["fun_traits"] = json.dumps(generate_fun_traits_json(ft_result))
        print(f"Fun Traits: {ft_result.found} variants found")
    except Exception as e:
        print(f"Fun Traits analysis failed: {e}")

    # ── PHASE 3: ANCESTRY REPORTS ────────────────────────────────────────────

    # Promote Viking heritage to standalone report
    try:
        viking_result_standalone = analyze_viking(variants)
        reports["viking"] = json.dumps(generate_viking_json(viking_result_standalone))
        print(f"Viking (standalone): {len(viking_result_standalone['used_snps'])} SNPs used")
    except Exception as e:
        print(f"Viking standalone analysis failed: {e}")

    # Promote Italian ancestry to standalone report
    try:
        italian_result_standalone = analyze_italian(variants)
        reports["italian_ancestry"] = json.dumps(generate_italian_json(italian_result_standalone))
        print(f"Italian (standalone): {len(italian_result_standalone['used_snps'])} SNPs used")
    except Exception as e:
        print(f"Italian standalone analysis failed: {e}")

    # Celtic ancestry
    try:
        celtic_result = analyze_celtic(variants)
        reports["celtic"] = json.dumps(generate_celtic_json(celtic_result))
        print(f"Celtic ancestry: {len(celtic_result['used_snps'])} SNPs used")
    except Exception as e:
        print(f"Celtic analysis failed: {e}")

    # British Isles ancestry
    try:
        british_result = analyze_british_isles(variants)
        reports["british_isles"] = json.dumps(generate_british_isles_json(british_result))
        print(f"British Isles: {len(british_result['used_snps'])} SNPs used")
    except Exception as e:
        print(f"British Isles analysis failed: {e}")

    # Native American ancestry
    try:
        na_result = analyze_native_american(variants)
        reports["native_american"] = json.dumps(generate_native_american_json(na_result))
        print(f"Native American: {len(na_result['used_snps'])} SNPs used")
    except Exception as e:
        print(f"Native American analysis failed: {e}")

    # Latino breakdown
    try:
        latino_result = analyze_latino_breakdown(variants)
        reports["latino_breakdown"] = json.dumps(generate_latino_breakdown_json(latino_result))
        print(f"Latino breakdown: {len(latino_result['used_snps'])} SNPs used")
    except Exception as e:
        print(f"Latino breakdown analysis failed: {e}")

    # African breakdown
    try:
        african_result = analyze_african_breakdown(variants)
        reports["african_breakdown"] = json.dumps(generate_african_breakdown_json(african_result))
        print(f"African breakdown: {len(african_result['used_snps'])} SNPs used")
    except Exception as e:
        print(f"African breakdown analysis failed: {e}")

    # European breakdown
    try:
        european_result = analyze_european_breakdown(variants)
        reports["european_breakdown"] = json.dumps(generate_european_breakdown_json(european_result))
        print(f"European breakdown: {len(european_result['used_snps'])} SNPs used")
    except Exception as e:
        print(f"European breakdown analysis failed: {e}")

    # Hunter-Gatherer vs Farmer (3-wave model refinement)
    try:
        hgf_result = analyze_hunter_gatherer_farmer(variants)
        reports["hunter_gatherer_farmer"] = json.dumps(generate_hunter_gatherer_farmer_json(hgf_result))
        print(f"Hunter-Gatherer vs Farmer: {len(hgf_result['used_snps'])} SNPs used")
    except Exception as e:
        print(f"Hunter-Gatherer vs Farmer analysis failed: {e}")

    # Deep ancestry (out-of-Africa + archaic admixture)
    try:
        deep_result = analyze_deep_ancestry(variants)
        reports["deep_ancestry"] = json.dumps(generate_deep_ancestry_json(deep_result))
        print(f"Deep ancestry: {len(deep_result['used_snps'])} SNPs used")
    except Exception as e:
        print(f"Deep ancestry analysis failed: {e}")

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
