"""
Neanderthal Introgression Analyzer
Analyzes DNA data for Neanderthal-derived variants.
Returns structured JSON for frontend visualization.

References:
- Vernot & Akey (2014) Science
- Prufer et al. (2014) Nature
- Sankararaman et al. (2014) Nature
- Dannemann & Kelso (2017) AJHG
- Simonti et al. (2016) Science
- Zeberg & Paabo (2020, 2021) Nature
"""

from typing import Dict, Tuple, Any, List
from dataclasses import dataclass, field, asdict

# ─────────────────────────────────────────────────────────────────────────────
# NEANDERTHAL VARIANT CATALOG
# ─────────────────────────────────────────────────────────────────────────────

NEANDERTHAL_VARIANTS = {
    # ── SKIN / HAIR / PIGMENTATION ──
    "rs10756819":  {"gene": "BNC2",    "chr": "9",  "nea_allele": "T", "category": "Pigmentation",
                    "phenotype": "Skin color (lighter skin in Europeans)", "confidence": "High",
                    "reference": "Vernot & Akey 2014"},
    "rs2153271":   {"gene": "BNC2",    "chr": "9",  "nea_allele": "T", "category": "Pigmentation",
                    "phenotype": "Skin color (lighter skin in Europeans)", "confidence": "High",
                    "reference": "Vernot & Akey 2014"},
    "rs10962897":  {"gene": "BNC2",    "chr": "9",  "nea_allele": "C", "category": "Pigmentation",
                    "phenotype": "Skin color / freckling", "confidence": "High",
                    "reference": "Vernot & Akey 2014"},
    "rs2048683":   {"gene": "BNC2",    "chr": "9",  "nea_allele": "T", "category": "Pigmentation",
                    "phenotype": "Skin pigmentation", "confidence": "High",
                    "reference": "Vernot & Akey 2014"},
    "rs10814359":  {"gene": "HYAL2",   "chr": "3",  "nea_allele": "A", "category": "Pigmentation",
                    "phenotype": "Skin adaptation / UV response", "confidence": "High",
                    "reference": "Vernot & Akey 2014"},
    "rs1800407":   {"gene": "OCA2",    "chr": "15", "nea_allele": "T", "category": "Pigmentation",
                    "phenotype": "Blue/green eye color", "confidence": "Moderate",
                    "reference": "Dannemann & Kelso 2017"},
    "rs12913832":  {"gene": "HERC2",   "chr": "15", "nea_allele": "A", "category": "Pigmentation",
                    "phenotype": "Blue eye color (primary determinant)", "confidence": "Moderate",
                    "reference": "Dannemann & Kelso 2017"},
    "rs1805007":   {"gene": "MC1R",    "chr": "16", "nea_allele": "T", "category": "Pigmentation",
                    "phenotype": "Red hair / fair skin / sun sensitivity", "confidence": "Moderate",
                    "reference": "Prufer et al. 2014"},
    "rs1805008":   {"gene": "MC1R",    "chr": "16", "nea_allele": "T", "category": "Pigmentation",
                    "phenotype": "Red hair color (R151C variant)", "confidence": "Moderate",
                    "reference": "Prufer et al. 2014"},
    "rs1426654":   {"gene": "SLC24A5", "chr": "15", "nea_allele": "A", "category": "Pigmentation",
                    "phenotype": "Skin lightening (major European pigmentation allele)", "confidence": "High",
                    "reference": "Sankararaman et al. 2014"},
    "rs2733832":   {"gene": "TYRP1",   "chr": "9",  "nea_allele": "T", "category": "Pigmentation",
                    "phenotype": "Hair and skin color variation", "confidence": "Moderate",
                    "reference": "Vernot & Akey 2014"},

    # ── IMMUNE SYSTEM ──
    "rs1131454":   {"gene": "STAT2",   "chr": "12", "nea_allele": "A", "category": "Immune",
                    "phenotype": "Antiviral innate immunity / interferon signaling", "confidence": "High",
                    "reference": "Sankararaman et al. 2014"},
    "rs2066807":   {"gene": "STAT2",   "chr": "12", "nea_allele": "G", "category": "Immune",
                    "phenotype": "Antiviral innate immunity", "confidence": "High",
                    "reference": "Mendez et al. 2012"},
    "rs10774671":  {"gene": "OAS1",    "chr": "12", "nea_allele": "G", "category": "Immune",
                    "phenotype": "Antiviral defense (COVID-19 protection)", "confidence": "High",
                    "reference": "Sams et al. 2016; Zeberg et al. 2021"},
    "rs2660":      {"gene": "OAS1",    "chr": "12", "nea_allele": "G", "category": "Immune",
                    "phenotype": "OAS1 antiviral activity", "confidence": "High",
                    "reference": "Zeberg et al. 2021"},
    "rs1293740":   {"gene": "OAS3",    "chr": "12", "nea_allele": "G", "category": "Immune",
                    "phenotype": "Antiviral innate immunity", "confidence": "High",
                    "reference": "Sams et al. 2016"},
    "rs2072136":   {"gene": "OAS2",    "chr": "12", "nea_allele": "A", "category": "Immune",
                    "phenotype": "Antiviral defense", "confidence": "Moderate",
                    "reference": "Sams et al. 2016"},
    "rs2248374":   {"gene": "HLA-C",   "chr": "6",  "nea_allele": "G", "category": "Immune",
                    "phenotype": "Immune cell recognition (NK cell regulation)", "confidence": "High",
                    "reference": "Abi-Rached et al. 2011"},
    "rs9264942":   {"gene": "HLA-C",   "chr": "6",  "nea_allele": "T", "category": "Immune",
                    "phenotype": "HLA-C expression level / HIV control", "confidence": "High",
                    "reference": "Abi-Rached et al. 2011"},
    "rs2523393":   {"gene": "HLA-B",   "chr": "6",  "nea_allele": "C", "category": "Immune",
                    "phenotype": "Immune response / pathogen resistance", "confidence": "Moderate",
                    "reference": "Abi-Rached et al. 2011"},
    "rs3094626":   {"gene": "HLA-A",   "chr": "6",  "nea_allele": "G", "category": "Immune",
                    "phenotype": "Immune response (MHC class I)", "confidence": "Moderate",
                    "reference": "Abi-Rached et al. 2011"},
    "rs20575":     {"gene": "TNFRSF10A","chr": "8", "nea_allele": "C", "category": "Immune",
                    "phenotype": "Apoptosis regulation / tumor necrosis factor signaling", "confidence": "Moderate",
                    "reference": "Simonti et al. 2016"},
    "rs1141747":   {"gene": "IFIT1",   "chr": "10", "nea_allele": "T", "category": "Immune",
                    "phenotype": "Interferon-induced antiviral response", "confidence": "Moderate",
                    "reference": "Sams et al. 2016"},
    "rs73885319":  {"gene": "APOL1",   "chr": "22", "nea_allele": "G", "category": "Immune",
                    "phenotype": "Kidney function / sleeping sickness resistance", "confidence": "Moderate",
                    "reference": "Prufer et al. 2014"},

    # ── METABOLISM ──
    "rs13342232":  {"gene": "SLC16A11","chr": "17", "nea_allele": "C", "category": "Metabolism",
                    "phenotype": "Type 2 diabetes risk (Neanderthal haplotype)", "confidence": "High",
                    "reference": "SIGMA Consortium 2014"},
    "rs75418188":  {"gene": "SLC16A11","chr": "17", "nea_allele": "A", "category": "Metabolism",
                    "phenotype": "Type 2 diabetes risk / fat metabolism", "confidence": "High",
                    "reference": "SIGMA Consortium 2014"},
    "rs117767867": {"gene": "SLC16A11","chr": "17", "nea_allele": "T", "category": "Metabolism",
                    "phenotype": "Type 2 diabetes / liver lipid metabolism", "confidence": "High",
                    "reference": "SIGMA Consortium 2014"},
    "rs4253778":   {"gene": "PPARA",   "chr": "22", "nea_allele": "C", "category": "Metabolism",
                    "phenotype": "Fatty acid oxidation / lipid metabolism", "confidence": "Moderate",
                    "reference": "Simonti et al. 2016"},
    "rs135539":    {"gene": "PPARA",   "chr": "22", "nea_allele": "G", "category": "Metabolism",
                    "phenotype": "Lipid metabolism", "confidence": "Moderate",
                    "reference": "Simonti et al. 2016"},
    "rs1137101":   {"gene": "LEPR",    "chr": "1",  "nea_allele": "G", "category": "Metabolism",
                    "phenotype": "Leptin signaling / body weight regulation", "confidence": "Moderate",
                    "reference": "Simonti et al. 2016"},
    "rs4988235":   {"gene": "LCT",     "chr": "2",  "nea_allele": "A", "category": "Metabolism",
                    "phenotype": "Lactase persistence / dairy digestion", "confidence": "Moderate",
                    "reference": "Simonti et al. 2016"},
    "rs2271385":   {"gene": "SPTLC1",  "chr": "9",  "nea_allele": "G", "category": "Metabolism",
                    "phenotype": "Sphingolipid metabolism", "confidence": "Moderate",
                    "reference": "Vernot & Akey 2014"},
    "rs10830963":  {"gene": "MTNR1B",  "chr": "11", "nea_allele": "G", "category": "Metabolism",
                    "phenotype": "Melatonin receptor / glucose metabolism", "confidence": "Moderate",
                    "reference": "Dannemann & Kelso 2017"},
    "rs3088291":   {"gene": "WARS2",   "chr": "1",  "nea_allele": "T", "category": "Metabolism",
                    "phenotype": "Mitochondrial tRNA synthesis", "confidence": "Moderate",
                    "reference": "Racimo et al. 2017"},

    # ── NEUROLOGICAL ──
    "rs11932595":  {"gene": "CLOCK",   "chr": "4",  "nea_allele": "C", "category": "Neurological",
                    "phenotype": "Circadian rhythm / sleep timing", "confidence": "High",
                    "reference": "Simonti et al. 2016"},
    "rs1801260":   {"gene": "CLOCK",   "chr": "4",  "nea_allele": "C", "category": "Neurological",
                    "phenotype": "Circadian rhythm / morningness-eveningness", "confidence": "High",
                    "reference": "Simonti et al. 2016"},
    "rs3749474":   {"gene": "CLOCK",   "chr": "4",  "nea_allele": "T", "category": "Neurological",
                    "phenotype": "Circadian rhythm regulation", "confidence": "Moderate",
                    "reference": "Simonti et al. 2016"},
    "rs1020117":   {"gene": "POU2F3",  "chr": "11", "nea_allele": "A", "category": "Neurological",
                    "phenotype": "Neural development", "confidence": "Moderate",
                    "reference": "Vernot & Akey 2014"},
    "rs2305480":   {"gene": "GRIN2A",  "chr": "16", "nea_allele": "G", "category": "Neurological",
                    "phenotype": "Glutamate receptor (epilepsy / learning)", "confidence": "Moderate",
                    "reference": "Simonti et al. 2016"},
    "rs10500418":  {"gene": "RBFOX1",  "chr": "16", "nea_allele": "C", "category": "Neurological",
                    "phenotype": "RNA splicing regulation (neurological functions)", "confidence": "Moderate",
                    "reference": "Vernot & Akey 2014"},
    "rs6803513":   {"gene": "ROBO2",   "chr": "3",  "nea_allele": "G", "category": "Neurological",
                    "phenotype": "Axon guidance / brain development", "confidence": "Moderate",
                    "reference": "Vernot & Akey 2014"},
    "rs2301893":   {"gene": "CNTN5",   "chr": "11", "nea_allele": "A", "category": "Neurological",
                    "phenotype": "Neuronal contact / synaptic function", "confidence": "Moderate",
                    "reference": "Vernot & Akey 2014"},
    "rs7801835":   {"gene": "AUTS2",   "chr": "7",  "nea_allele": "T", "category": "Neurological",
                    "phenotype": "Cognitive development", "confidence": "Moderate",
                    "reference": "Simonti et al. 2016"},
    "rs7782257":   {"gene": "FOXP2",   "chr": "7",  "nea_allele": "G", "category": "Neurological",
                    "phenotype": "Speech and language development", "confidence": "Moderate",
                    "reference": "Sankararaman et al. 2014"},
    "rs10748842":  {"gene": "NRG3",    "chr": "10", "nea_allele": "C", "category": "Neurological",
                    "phenotype": "Neuregulin signaling", "confidence": "Moderate",
                    "reference": "Simonti et al. 2016"},
    "rs2835740":   {"gene": "DYRK1A",  "chr": "21", "nea_allele": "A", "category": "Neurological",
                    "phenotype": "Cognitive development (DYRK1A kinase)", "confidence": "Moderate",
                    "reference": "Simonti et al. 2016"},

    # ── RESPIRATORY / COVID ──
    "rs35044562":  {"gene": "LZTFL1",  "chr": "3",  "nea_allele": "A", "category": "Respiratory",
                    "phenotype": "Severe COVID-19 risk (Neanderthal haplotype on chr3)", "confidence": "High",
                    "reference": "Zeberg & Paabo 2020"},
    "rs17713054":  {"gene": "LZTFL1",  "chr": "3",  "nea_allele": "G", "category": "Respiratory",
                    "phenotype": "Severe COVID-19 risk", "confidence": "High",
                    "reference": "Zeberg & Paabo 2020"},
    "rs11385942":  {"gene": "SLC6A20", "chr": "3",  "nea_allele": "GA","category": "Respiratory",
                    "phenotype": "Severe COVID-19 / respiratory failure risk", "confidence": "High",
                    "reference": "Zeberg & Paabo 2020"},

    # ── COAGULATION ──
    "rs2289252":   {"gene": "F11",     "chr": "4",  "nea_allele": "T", "category": "Coagulation",
                    "phenotype": "Blood clotting (coagulation factor XI level)", "confidence": "High",
                    "reference": "Simonti et al. 2016"},
    "rs4253417":   {"gene": "F11",     "chr": "4",  "nea_allele": "C", "category": "Coagulation",
                    "phenotype": "Blood clotting factor XI", "confidence": "High",
                    "reference": "Zeberg et al. 2020"},

    # ── PAIN ──
    "rs6746030":   {"gene": "SCN9A",   "chr": "2",  "nea_allele": "A", "category": "Pain",
                    "phenotype": "Pain sensitivity (Nav1.7 sodium channel)", "confidence": "High",
                    "reference": "Zeberg et al. 2020"},
    "rs8065080":   {"gene": "TRPV1",   "chr": "17", "nea_allele": "C", "category": "Pain",
                    "phenotype": "Pain / temperature sensing (capsaicin receptor)", "confidence": "Moderate",
                    "reference": "Zeberg et al. 2020"},

    # ── MORPHOLOGY ──
    "rs984222":    {"gene": "TBX15",   "chr": "1",  "nea_allele": "G", "category": "Morphology",
                    "phenotype": "Fat distribution / body shape", "confidence": "High",
                    "reference": "Racimo et al. 2017"},
    "rs4926488":   {"gene": "TBX15",   "chr": "1",  "nea_allele": "G", "category": "Morphology",
                    "phenotype": "Body fat distribution", "confidence": "High",
                    "reference": "Racimo et al. 2017"},
    "rs3827760":   {"gene": "EDAR",    "chr": "2",  "nea_allele": "A", "category": "Morphology",
                    "phenotype": "Hair thickness / sweat gland density", "confidence": "Moderate",
                    "reference": "Gittelman et al. 2016"},

    # ── SKIN STRUCTURE ──
    "rs3829241":   {"gene": "KRT71",   "chr": "12", "nea_allele": "A", "category": "Skin structure",
                    "phenotype": "Keratin (hair shaft structure / curliness)", "confidence": "Moderate",
                    "reference": "Vernot & Akey 2014"},
    "rs34400929":  {"gene": "KRT84",   "chr": "12", "nea_allele": "T", "category": "Skin structure",
                    "phenotype": "Keratin structure", "confidence": "Moderate",
                    "reference": "Vernot & Akey 2014"},

    # ── CARDIOVASCULAR ──
    "rs216009":    {"gene": "CACNA1C", "chr": "12", "nea_allele": "G", "category": "Cardiovascular",
                    "phenotype": "Cardiac / neurological calcium channel", "confidence": "Moderate",
                    "reference": "Simonti et al. 2016"},
    "rs6746206":   {"gene": "PLN",     "chr": "6",  "nea_allele": "A", "category": "Cardiovascular",
                    "phenotype": "Cardiac muscle regulation", "confidence": "Moderate",
                    "reference": "Simonti et al. 2016"},

    # ── OTHER ──
    "rs3748960":   {"gene": "ADAMTS3", "chr": "4",  "nea_allele": "C", "category": "Connective tissue",
                    "phenotype": "Extracellular matrix / joint development", "confidence": "Moderate",
                    "reference": "Vernot & Akey 2014"},
    "rs2078742":   {"gene": "LARGE1",  "chr": "22", "nea_allele": "T", "category": "Cellular",
                    "phenotype": "Protein glycosylation / muscular dystrophy pathway", "confidence": "Moderate",
                    "reference": "Vernot & Akey 2014"},
    "rs28399433":  {"gene": "CYP2A6",  "chr": "19", "nea_allele": "A", "category": "Pharmacogenomics",
                    "phenotype": "Nicotine metabolism / drug metabolism", "confidence": "Moderate",
                    "reference": "Simonti et al. 2016"},
    "rs13419896":  {"gene": "EPAS1",   "chr": "2",  "nea_allele": "C", "category": "Altitude adaptation",
                    "phenotype": "Hypoxia response / high-altitude adaptation", "confidence": "Moderate",
                    "reference": "Huerta-Sanchez et al. 2014"},
    "rs1805495":   {"gene": "PDE6A",   "chr": "5",  "nea_allele": "T", "category": "Sensory",
                    "phenotype": "Photoreceptor function", "confidence": "Moderate",
                    "reference": "Vernot & Akey 2014"},
    "rs187238":    {"gene": "IL18",    "chr": "11", "nea_allele": "C", "category": "Immune",
                    "phenotype": "Interleukin-18 / inflammatory response", "confidence": "Moderate",
                    "reference": "Simonti et al. 2016"},
    "rs1800629":   {"gene": "TNF",     "chr": "6",  "nea_allele": "A", "category": "Immune",
                    "phenotype": "TNF-alpha / inflammatory cytokine", "confidence": "Moderate",
                    "reference": "Simonti et al. 2016"},
    "rs7217270":   {"gene": "TRPV3",   "chr": "17", "nea_allele": "A", "category": "Sensory",
                    "phenotype": "Temperature sensation (thermosensory TRP channel)", "confidence": "Moderate",
                    "reference": "Zeberg et al. 2020"},
}


def _determine_zygosity(user_alleles: str, nea_allele: str) -> str:
    """Determine if the user carries the Neanderthal allele."""
    if not user_alleles or len(user_alleles) < 2 or user_alleles in ('--', '00', 'II', 'DD'):
        return 'unknown'
    nea = nea_allele.upper()
    a1, a2 = user_alleles[0].upper(), user_alleles[1].upper()
    count = sum(1 for a in [a1, a2] if a == nea)
    if count == 2:
        return 'homozygous'
    elif count == 1:
        return 'heterozygous'
    return 'not_carried'


def _estimate_neanderthal_percentage(carried: int, total_checked: int) -> float:
    """Rough estimate of Neanderthal ancestry fraction (1-4% range)."""
    if total_checked == 0:
        return 0.0
    carry_rate = carried / total_checked
    estimated_pct = carry_rate * (2.0 / 0.45)
    return min(round(estimated_pct, 2), 4.0)


def analyze_neanderthal(variants: Dict[str, Tuple[str, str, str]]) -> Dict[str, Any]:
    """
    Analyze variants for Neanderthal introgression.

    Args:
        variants: Dict mapping rsID -> (chromosome, position, genotype)
                  as parsed by run_full_analysis.parse_genome_file()

    Returns:
        Structured JSON dict for the neanderthal report
    """
    carried_variants: List[Dict[str, Any]] = []
    not_carried_variants: List[Dict[str, Any]] = []
    not_found_rsids: List[str] = []
    category_stats: Dict[str, Dict[str, int]] = {}

    for rsid, info in NEANDERTHAL_VARIANTS.items():
        rsid_lower = rsid.lower()
        variant_data = variants.get(rsid_lower)

        if variant_data is None:
            not_found_rsids.append(rsid)
            continue

        _chrom, _pos, genotype = variant_data
        zygosity = _determine_zygosity(genotype, info['nea_allele'])

        entry = {
            'rsid': rsid,
            'gene': info['gene'],
            'chromosome': info['chr'],
            'neanderthalAllele': info['nea_allele'],
            'userGenotype': genotype,
            'zygosity': zygosity,
            'category': info['category'],
            'phenotype': info['phenotype'],
            'confidence': info['confidence'],
            'reference': info['reference'],
        }

        if zygosity in ('homozygous', 'heterozygous'):
            carried_variants.append(entry)
        elif zygosity == 'not_carried':
            not_carried_variants.append(entry)

        # Track category stats
        cat = info['category']
        if cat not in category_stats:
            category_stats[cat] = {'carried': 0, 'total': 0}
        category_stats[cat]['total'] += 1
        if zygosity in ('homozygous', 'heterozygous'):
            category_stats[cat]['carried'] += 1

    total_checked = len(carried_variants) + len(not_carried_variants)
    n_carried = len(carried_variants)
    n_homo = sum(1 for v in carried_variants if v['zygosity'] == 'homozygous')
    n_het = sum(1 for v in carried_variants if v['zygosity'] == 'heterozygous')

    carry_rate = round(n_carried / total_checked, 3) if total_checked > 0 else 0
    estimated_pct = _estimate_neanderthal_percentage(n_carried, total_checked)

    # Build categories array sorted by carried count desc
    categories = []
    for cat, stats in sorted(category_stats.items(), key=lambda x: x[1]['carried'], reverse=True):
        categories.append({
            'name': cat,
            'carried': stats['carried'],
            'total': stats['total'],
        })

    # Sort carried variants by category
    carried_variants.sort(key=lambda v: (v['category'], v['gene']))

    # High-confidence highlights
    highlights = [v for v in carried_variants if v['confidence'] == 'High']

    return {
        'summary': {
            'catalogSize': len(NEANDERTHAL_VARIANTS),
            'variantsChecked': total_checked,
            'variantsNotOnChip': len(not_found_rsids),
            'variantsCarried': n_carried,
            'homozygous': n_homo,
            'heterozygous': n_het,
            'notCarried': len(not_carried_variants),
            'carryRate': carry_rate,
            'estimatedPercentage': estimated_pct,
        },
        'categories': categories,
        'carriedVariants': carried_variants,
        'highlights': highlights,
    }
