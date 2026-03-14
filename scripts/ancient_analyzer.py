"""
Ancient Bloodlines Analyzer
Estimates genetic affinity to ancient human populations using Ancestry-Informative
Markers (AIMs) with allele frequencies derived from published archaeogenomics literature.

ANCIENT POPULATIONS:
  WHG          Western Hunter-Gatherer (~8000 BCE)
  EEF          Early European Farmer / Anatolian Neolithic (~6000-4000 BCE)
  Steppe       Bronze Age Steppe / Yamnaya (~3000-2500 BCE)
  CHG          Caucasus Hunter-Gatherer (~10000-6000 BCE)
  EHG          Eastern Hunter-Gatherer (~8000-5000 BCE)
  Natufian     Levantine Hunter-Gatherer (~12000-9000 BCE)
  Ancient_EAS  Ancient East Asian proxy (~8000-3000 BCE)
  Ancient_AFR  Ancient Sub-Saharan African proxy (modern YRI as reference)

Based on:
  - Lazaridis et al. 2014 (Nature): Three ancestral populations for Europeans
  - Haak et al. 2015 (Nature): Massive migration from the steppe
  - Mathieson et al. 2015 (Nature): Selection in 230 ancient Eurasians
  - Olalde et al. 2014 (Nature): WHG pigmentation (Cheddar Man)
  - Lazaridis et al. 2016 (Nature): Natufian and Basal Eurasian ancestry

Returns structured JSON for frontend visualization.
"""

from typing import Dict, Tuple, Any, Optional
import numpy as np
from scipy.optimize import nnls

# ---------------------------------------------------------------------------
# POPULATION DEFINITIONS
# ---------------------------------------------------------------------------

ANCIENT_POPS = [
    "WHG", "EEF", "Steppe", "CHG", "EHG", "Natufian", "Ancient_EAS", "Ancient_AFR"
]

ANCIENT_POP_LABELS = {
    "WHG":         "Western Hunter-Gatherer (~8000 BCE)",
    "EEF":         "Early European Farmer / Anatolian Neolithic (~6000-4000 BCE)",
    "Steppe":      "Bronze Age Steppe / Yamnaya Indo-Europeans (~3000-2500 BCE)",
    "CHG":         "Caucasus Hunter-Gatherer (~10000-6000 BCE)",
    "EHG":         "Eastern Hunter-Gatherer (~8000-5000 BCE)",
    "Natufian":    "Natufian / Levantine Hunter-Gatherer (~12000-9000 BCE)",
    "Ancient_EAS": "Ancient East Asian (~8000-3000 BCE)",
    "Ancient_AFR": "Ancient Sub-Saharan African (proxy: modern Yoruba/YRI)",
}

ANCIENT_POP_PERIODS = {
    "WHG":         "~8000 BCE",
    "EEF":         "~6000-4000 BCE",
    "Steppe":      "~3000-2500 BCE",
    "CHG":         "~10000-6000 BCE",
    "EHG":         "~8000-5000 BCE",
    "Natufian":    "~12000-9000 BCE",
    "Ancient_EAS": "~8000-3000 BCE",
    "Ancient_AFR": "proxy: modern Yoruba",
}

ANCIENT_POP_COLORS = {
    "WHG":         "#3B82F6",  # blue
    "EEF":         "#F59E0B",  # amber
    "Steppe":      "#EF4444",  # red
    "CHG":         "#8B5CF6",  # violet
    "EHG":         "#06B6D4",  # cyan
    "Natufian":    "#10B981",  # emerald
    "Ancient_EAS": "#F97316",  # orange
    "Ancient_AFR": "#14B8A6",  # teal
}

ANCIENT_POP_DESCRIPTIONS = {
    "WHG": (
        "Western Hunter-Gatherers were the pre-farming inhabitants of Europe before "
        "agriculture arrived from Anatolia. They lived as foragers in western and "
        "central Europe from ~40,000 BCE until ~6,000 BCE when farmers displaced "
        "and admixed with them. Key genetic findings: they had DARK SKIN despite "
        "living in Europe (SLC24A5 ancestral allele dominant), but BLUE EYES "
        "(HERC2/OCA2 derived allele already present at ~70%). This means the "
        "evolution of blue eyes in Europeans PRECEDED the evolution of light skin. "
        "They had NO lactase persistence — dairy was not part of their diet. "
        "Famous individuals: Cheddar Man (England, ~7100 BCE), La Braña 1 "
        "(Spain, ~7000 BCE), Loschbour (Luxembourg, ~8000 BCE)."
    ),
    "EEF": (
        "Early European Farmers descended from Anatolian Neolithic populations who "
        "crossed the Aegean and Bosphorus from what is now Turkey, bringing "
        "agriculture to Europe ~8,000-9,000 years ago. They replaced and partially "
        "admixed with the WHG population. Key genetic findings: they had LIGHT SKIN "
        "(SLC24A5 near-fixed, ~95%) but mostly DARK/BROWN EYES — a combination that "
        "is uncommon today. They had NO lactase persistence (dairy fermentation was "
        "practiced but digestion of raw milk was poor). EEF ancestry is the primary "
        "component of modern Southern Europeans (Italians, Iberians, Greeks). "
        "Famous ancient individuals: Stuttgart woman (Germany, ~5500 BCE), "
        "Ötzi the Iceman (Alps, ~3300 BCE — essentially pure EEF)."
    ),
    "Steppe": (
        "Bronze Age Steppe people (Yamnaya culture and related groups) were pastoral "
        "nomads from the Pontic-Caspian steppe (modern Ukraine/Russia/Kazakhstan) "
        "who swept across Europe and Central Asia ~5,000 years ago. They are the "
        "ancestors of most Indo-European language speakers today. Genetically, they "
        "themselves were a mix of CHG and EHG ancestry. Key finding: they are "
        "associated with the RISE OF LACTASE PERSISTENCE in Europe — their "
        "cattle-herding lifestyle created strong selection for milk digestion. They "
        "brought lactase persistence alleles, Indo-European languages, wheel "
        "technology, and horse domestication to Europe. Yamnaya ancestry is highest "
        "in modern Northern and Eastern Europeans (Scandinavians, Slavs) and lowest "
        "in Sardinians and Basques."
    ),
    "CHG": (
        "Caucasus Hunter-Gatherers inhabited the South Caucasus (modern Georgia, "
        "Armenia, Azerbaijan) from at least ~30,000 BCE. They contributed heavily to "
        "the Yamnaya Bronze Age Steppe people (~50% CHG + ~50% EHG = Yamnaya). CHG "
        "are genetically intermediate between Middle Eastern and European populations, "
        "reflecting their geographic position as a bridge between the Near East and "
        "Europe. They are also the major ancestor of many modern South Caucasian "
        "populations (Georgians, Armenians). Key ancient individuals: Kotias Klde "
        "(Georgia, ~9,700 BCE), Satsurblia (Georgia, ~13,300 BCE)."
    ),
    "EHG": (
        "Eastern Hunter-Gatherers inhabited the eastern European plain and western "
        "Siberia, roughly present-day Russia. They were closely related to WHG but "
        "had a distinct component from Ancient North Eurasians (ANE) — the population "
        "that also contributed to Native Americans. EHG + CHG admixed to form the "
        "Yamnaya Steppe people. Key ancient individuals: Karelia hunter-gatherer "
        "(Russia, ~7,500 BCE), Samara hunter-gatherer (Russia, ~5,500 BCE)."
    ),
    "Natufian": (
        "Natufians were the Levantine hunter-gatherers of the ancient Near East "
        "(modern Israel, Palestine, Jordan, Lebanon) who lived immediately before "
        "and during the invention of agriculture (~15,000-9,000 BCE). They are the "
        "primary ancestor of the Early European Farmers (EEF). Natufians had dark "
        "skin (ancestral SLC24A5 allele), dark eyes, and a mixed diet of hunted "
        "game and wild cereals. Their descendants became the first farmers, and "
        "then spread that farming culture to Anatolia and ultimately Europe. "
        "Genetically, they cluster with modern Middle Eastern populations more "
        "than with Europeans — they represent the 'Basal Eurasian' lineage that "
        "separated from other Eurasians before ~50,000 BCE."
    ),
    "Ancient_EAS": (
        "Ancient East Asian populations from the Holocene period (Yellow River basin, "
        "Korean Peninsula, Japanese archipelago ancestors). The EDAR 370A variant "
        "(rs3827760) was already at high frequency in ancient East Asians, indicating "
        "this selection occurred deep in East Asian prehistory. The ALDH2*2 allele "
        "(Asian alcohol flush reaction) and dry earwax allele were also present at "
        "moderate frequencies. These ancient East Asians are the ancestors of modern "
        "Han Chinese, Japanese, Korean, and related populations."
    ),
    "Ancient_AFR": (
        "Modern West African populations (Yoruba) serve as the best available proxy "
        "for ancient sub-Saharan African ancestry in this analysis. The Duffy-null "
        "allele (rs2814778*C) — conferring resistance to Plasmodium vivax malaria — "
        "is present at ~97% in this reference group, making it the most powerful "
        "single-marker identifier of African ancestry. Note: Africa contains the "
        "most genetic diversity of any continent; this reference represents West "
        "African ancestry specifically, not all African lineages."
    ),
}

# ---------------------------------------------------------------------------
# AIM ALLELE FREQUENCIES FOR ANCIENT POPULATIONS
# Derived from Mathieson 2015, Haak 2015, Lazaridis 2014/2016, gnomAD v3
# For each rsid: (ref, alt, weight, {pop: alt_freq})
# ---------------------------------------------------------------------------

AIMS_ANCIENT = {
    # Group 1: Ultra-high-FST continental markers
    "rs2814778": ("T", "C", 4.0, {
        "WHG": 0.00, "EEF": 0.00, "Steppe": 0.00, "CHG": 0.00,
        "EHG": 0.00, "Natufian": 0.01, "Ancient_EAS": 0.00, "Ancient_AFR": 0.95,
    }),
    "rs3827760": ("A", "G", 4.0, {
        "WHG": 0.00, "EEF": 0.00, "Steppe": 0.02, "CHG": 0.01,
        "EHG": 0.01, "Natufian": 0.01, "Ancient_EAS": 0.85, "Ancient_AFR": 0.00,
    }),
    "rs671": ("G", "A", 4.0, {
        "WHG": 0.00, "EEF": 0.00, "Steppe": 0.00, "CHG": 0.00,
        "EHG": 0.00, "Natufian": 0.00, "Ancient_EAS": 0.22, "Ancient_AFR": 0.00,
    }),
    "rs1229984": ("G", "A", 3.5, {
        "WHG": 0.02, "EEF": 0.04, "Steppe": 0.05, "CHG": 0.08,
        "EHG": 0.03, "Natufian": 0.12, "Ancient_EAS": 0.60, "Ancient_AFR": 0.02,
    }),
    "rs17822931": ("C", "T", 3.5, {
        "WHG": 0.01, "EEF": 0.02, "Steppe": 0.02, "CHG": 0.04,
        "EHG": 0.02, "Natufian": 0.04, "Ancient_EAS": 0.65, "Ancient_AFR": 0.01,
    }),
    "rs1426654": ("A", "G", 2.0, {
        "WHG": 0.12, "EEF": 0.96, "Steppe": 0.88, "CHG": 0.75,
        "EHG": 0.45, "Natufian": 0.18, "Ancient_EAS": 0.01, "Ancient_AFR": 0.02,
    }),

    # Group 2: European pigmentation / lactase markers
    "rs16891982": ("C", "G", 2.5, {
        "WHG": 0.30, "EEF": 0.55, "Steppe": 0.72, "CHG": 0.45,
        "EHG": 0.50, "Natufian": 0.30, "Ancient_EAS": 0.02, "Ancient_AFR": 0.01,
    }),
    "rs12913832": ("A", "G", 2.5, {
        "WHG": 0.72, "EEF": 0.26, "Steppe": 0.45, "CHG": 0.35,
        "EHG": 0.60, "Natufian": 0.14, "Ancient_EAS": 0.06, "Ancient_AFR": 0.04,
    }),
    "rs4988235": ("A", "G", 2.5, {
        "WHG": 0.02, "EEF": 0.02, "Steppe": 0.22, "CHG": 0.08,
        "EHG": 0.05, "Natufian": 0.03, "Ancient_EAS": 0.02, "Ancient_AFR": 0.11,
    }),
    "rs12203592": ("C", "T", 2.0, {
        "WHG": 0.22, "EEF": 0.20, "Steppe": 0.20, "CHG": 0.18,
        "EHG": 0.20, "Natufian": 0.14, "Ancient_EAS": 0.00, "Ancient_AFR": 0.00,
    }),
    "rs1805007": ("C", "T", 2.0, {
        "WHG": 0.07, "EEF": 0.05, "Steppe": 0.08, "CHG": 0.04,
        "EHG": 0.07, "Natufian": 0.03, "Ancient_EAS": 0.00, "Ancient_AFR": 0.00,
    }),
    "rs1800562": ("G", "A", 2.5, {
        "WHG": 0.01, "EEF": 0.01, "Steppe": 0.02, "CHG": 0.01,
        "EHG": 0.01, "Natufian": 0.01, "Ancient_EAS": 0.00, "Ancient_AFR": 0.00,
    }),
    "rs1799945": ("C", "G", 2.0, {
        "WHG": 0.06, "EEF": 0.06, "Steppe": 0.08, "CHG": 0.05,
        "EHG": 0.06, "Natufian": 0.04, "Ancient_EAS": 0.01, "Ancient_AFR": 0.01,
    }),
    "rs1042602": ("C", "A", 2.0, {
        "WHG": 0.22, "EEF": 0.32, "Steppe": 0.35, "CHG": 0.26,
        "EHG": 0.26, "Natufian": 0.22, "Ancient_EAS": 0.02, "Ancient_AFR": 0.03,
    }),
    "rs1800414": ("C", "T", 0.4, {
        "WHG": 0.30, "EEF": 0.22, "Steppe": 0.35, "CHG": 0.20,
        "EHG": 0.28, "Natufian": 0.16, "Ancient_EAS": 0.65, "Ancient_AFR": 0.16,
    }),

    # Group 3: EAS/AMR markers
    "rs4918664": ("A", "G", 2.0, {
        "WHG": 0.11, "EEF": 0.12, "Steppe": 0.14, "CHG": 0.20,
        "EHG": 0.15, "Natufian": 0.30, "Ancient_EAS": 0.84, "Ancient_AFR": 0.02,
    }),
    "rs1800497": ("G", "A", 1.5, {
        "WHG": 0.20, "EEF": 0.19, "Steppe": 0.22, "CHG": 0.24,
        "EHG": 0.22, "Natufian": 0.26, "Ancient_EAS": 0.38, "Ancient_AFR": 0.39,
    }),
    "rs1799971": ("A", "G", 1.5, {
        "WHG": 0.10, "EEF": 0.11, "Steppe": 0.15, "CHG": 0.18,
        "EHG": 0.14, "Natufian": 0.18, "Ancient_EAS": 0.38, "Ancient_AFR": 0.04,
    }),
    "rs7554936": ("C", "T", 1.0, {
        "WHG": 0.65, "EEF": 0.66, "Steppe": 0.68, "CHG": 0.66,
        "EHG": 0.67, "Natufian": 0.62, "Ancient_EAS": 0.85, "Ancient_AFR": 0.03,
    }),

    # Group 4: Population-differentiated markers
    "rs8050136": ("C", "A", 1.5, {
        "WHG": 0.38, "EEF": 0.40, "Steppe": 0.40, "CHG": 0.36,
        "EHG": 0.38, "Natufian": 0.35, "Ancient_EAS": 0.12, "Ancient_AFR": 0.46,
    }),
    "rs4680": ("G", "A", 1.0, {
        "WHG": 0.44, "EEF": 0.43, "Steppe": 0.44, "CHG": 0.40,
        "EHG": 0.43, "Natufian": 0.38, "Ancient_EAS": 0.28, "Ancient_AFR": 0.23,
    }),
    "rs1726866": ("G", "A", 1.0, {
        "WHG": 0.42, "EEF": 0.42, "Steppe": 0.43, "CHG": 0.40,
        "EHG": 0.41, "Natufian": 0.40, "Ancient_EAS": 0.28, "Ancient_AFR": 0.24,
    }),
    "rs713598": ("G", "C", 1.0, {
        "WHG": 0.40, "EEF": 0.40, "Steppe": 0.41, "CHG": 0.38,
        "EHG": 0.39, "Natufian": 0.37, "Ancient_EAS": 0.26, "Ancient_AFR": 0.21,
    }),
    "rs1801133": ("G", "A", 1.0, {
        "WHG": 0.26, "EEF": 0.38, "Steppe": 0.32, "CHG": 0.32,
        "EHG": 0.28, "Natufian": 0.32, "Ancient_EAS": 0.23, "Ancient_AFR": 0.10,
    }),
    "rs1801131": ("A", "C", 1.0, {
        "WHG": 0.22, "EEF": 0.24, "Steppe": 0.26, "CHG": 0.22,
        "EHG": 0.23, "Natufian": 0.21, "Ancient_EAS": 0.16, "Ancient_AFR": 0.10,
    }),
    "rs1800629": ("G", "A", 1.0, {
        "WHG": 0.12, "EEF": 0.17, "Steppe": 0.13, "CHG": 0.20,
        "EHG": 0.12, "Natufian": 0.24, "Ancient_EAS": 0.07, "Ancient_AFR": 0.20,
    }),
    "rs4129267": ("T", "C", 1.0, {
        "WHG": 0.32, "EEF": 0.42, "Steppe": 0.40, "CHG": 0.48,
        "EHG": 0.36, "Natufian": 0.52, "Ancient_EAS": 0.33, "Ancient_AFR": 0.18,
    }),
    "rs2435357": ("G", "C", 1.0, {
        "WHG": 0.20, "EEF": 0.28, "Steppe": 0.32, "CHG": 0.40,
        "EHG": 0.22, "Natufian": 0.46, "Ancient_EAS": 0.30, "Ancient_AFR": 0.09,
    }),
    "rs2228479": ("A", "G", 1.0, {
        "WHG": 0.05, "EEF": 0.06, "Steppe": 0.09, "CHG": 0.15,
        "EHG": 0.07, "Natufian": 0.19, "Ancient_EAS": 0.16, "Ancient_AFR": 0.05,
    }),

    # Group 5: HLA/population-differentiated markers
    "rs2076533": ("C", "A", 1.5, {
        "WHG": 0.12, "EEF": 0.18, "Steppe": 0.15, "CHG": 0.30,
        "EHG": 0.12, "Natufian": 0.42, "Ancient_EAS": 0.09, "Ancient_AFR": 0.07,
    }),
    "rs9267531": ("G", "A", 1.5, {
        "WHG": 0.15, "EEF": 0.22, "Steppe": 0.19, "CHG": 0.35,
        "EHG": 0.16, "Natufian": 0.44, "Ancient_EAS": 0.11, "Ancient_AFR": 0.09,
    }),
    "rs6090989": ("G", "A", 1.2, {
        "WHG": 0.10, "EEF": 0.18, "Steppe": 0.14, "CHG": 0.28,
        "EHG": 0.11, "Natufian": 0.38, "Ancient_EAS": 0.07, "Ancient_AFR": 0.04,
    }),
    "rs3135027": ("G", "A", 1.0, {
        "WHG": 0.16, "EEF": 0.25, "Steppe": 0.21, "CHG": 0.32,
        "EHG": 0.17, "Natufian": 0.38, "Ancient_EAS": 0.14, "Ancient_AFR": 0.11,
    }),

    # Group 6: Additional markers
    "rs10497520": ("C", "T", 0.8, {
        "WHG": 0.35, "EEF": 0.26, "Steppe": 0.38, "CHG": 0.28,
        "EHG": 0.36, "Natufian": 0.20, "Ancient_EAS": 0.08, "Ancient_AFR": 0.05,
    }),
    "rs1540771": ("A", "G", 0.8, {
        "WHG": 0.58, "EEF": 0.55, "Steppe": 0.62, "CHG": 0.50,
        "EHG": 0.60, "Natufian": 0.44, "Ancient_EAS": 0.20, "Ancient_AFR": 0.12,
    }),
    "rs1937845": ("A", "G", 0.8, {
        "WHG": 0.10, "EEF": 0.20, "Steppe": 0.25, "CHG": 0.38,
        "EHG": 0.15, "Natufian": 0.44, "Ancient_EAS": 0.22, "Ancient_AFR": 0.18,
    }),
    "rs2073711": ("C", "T", 0.8, {
        "WHG": 0.20, "EEF": 0.35, "Steppe": 0.30, "CHG": 0.40,
        "EHG": 0.24, "Natufian": 0.46, "Ancient_EAS": 0.31, "Ancient_AFR": 0.27,
    }),
    "rs7149477": ("A", "G", 0.8, {
        "WHG": 0.15, "EEF": 0.30, "Steppe": 0.22, "CHG": 0.38,
        "EHG": 0.18, "Natufian": 0.44, "Ancient_EAS": 0.19, "Ancient_AFR": 0.11,
    }),
}


# ---------------------------------------------------------------------------
# HELPERS
# ---------------------------------------------------------------------------

def _get_dosage(genotype: str, alt: str) -> Optional[int]:
    """Convert genotype string to alt allele dosage (0, 1, 2)."""
    if not genotype or genotype in ("--", "00", "NN", "NC", "DI", "II", "DD"):
        return None
    g = genotype.upper().strip()
    a = alt.upper()
    cnt = g.count(a)
    if cnt > 2 or len(g) < 2:
        return None
    return cnt


# ---------------------------------------------------------------------------
# ANCIENT ANCESTRY ESTIMATION (Weighted NNLS with Tikhonov regularization)
# ---------------------------------------------------------------------------

def _estimate_proportions(
    variants: Dict[str, Tuple[str, str, str]],
) -> Tuple[Dict[str, float], list, float]:
    """
    Run weighted NNLS decomposition into ancient population components.

    Returns (proportions_dict, used_snps_list, residual).
    """
    n_pops = len(ANCIENT_POPS)
    rows_A = []
    rows_f = []
    rows_w = []
    used_snps = []

    for rsid, (ref, alt, weight, pop_freqs) in AIMS_ANCIENT.items():
        variant_data = variants.get(rsid.lower())
        if variant_data is None:
            continue
        _chrom, _pos, genotype = variant_data
        dosage = _get_dosage(genotype, alt)
        if dosage is None:
            continue
        obs_freq = dosage / 2.0
        row = [pop_freqs.get(p, 0.0) for p in ANCIENT_POPS]
        rows_A.append(row)
        rows_f.append(obs_freq)
        rows_w.append(weight)
        used_snps.append(rsid)

    if len(used_snps) < 5:
        return {p: 1.0 / n_pops for p in ANCIENT_POPS}, used_snps, 0.0

    A = np.array(rows_A)
    f = np.array(rows_f)
    w = np.array(rows_w)

    # Weighted NNLS with Tikhonov regularization
    reg = 0.002
    A_w = A * w[:, np.newaxis]
    f_w = f * w
    A_aug = np.vstack([A_w, np.eye(n_pops) * reg])
    f_aug = np.concatenate([f_w, np.zeros(n_pops)])

    x, resid = nnls(A_aug, f_aug)
    total = x.sum()
    if total > 0:
        x /= total
    else:
        x = np.ones(n_pops) / n_pops

    props = {p: float(x[i]) for i, p in enumerate(ANCIENT_POPS)}
    return props, used_snps, float(resid)


def _bootstrap_ci(
    variants: Dict[str, Tuple[str, str, str]],
    n: int = 300,
    seed: int = 42,
) -> Tuple[Optional[np.ndarray], Optional[np.ndarray]]:
    """Bootstrap 95% CI for ancient ancestry estimates."""
    n_pops = len(ANCIENT_POPS)
    valid = []

    for rsid, (ref, alt, weight, pop_freqs) in AIMS_ANCIENT.items():
        variant_data = variants.get(rsid.lower())
        if variant_data is None:
            continue
        _chrom, _pos, genotype = variant_data
        dosage = _get_dosage(genotype, alt)
        if dosage is not None:
            valid.append((
                rsid,
                dosage / 2.0,
                weight,
                [pop_freqs.get(p, 0.0) for p in ANCIENT_POPS],
            ))

    if len(valid) < 8:
        return None, None

    rng = np.random.default_rng(seed)
    all_x = []

    for _ in range(n):
        boot = rng.choice(len(valid), size=len(valid), replace=True)
        A_b = np.array([valid[i][3] for i in boot])
        f_b = np.array([valid[i][1] for i in boot])
        w_b = np.array([valid[i][2] for i in boot])
        A_bw = A_b * w_b[:, np.newaxis]
        f_bw = f_b * w_b
        x, _ = nnls(A_bw, f_bw)
        s = x.sum()
        all_x.append(x / s if s > 0 else np.ones(n_pops) / n_pops)

    all_x = np.array(all_x)
    return np.percentile(all_x, 2.5, axis=0), np.percentile(all_x, 97.5, axis=0)


# ---------------------------------------------------------------------------
# PIGMENTATION INTERPRETATION
# ---------------------------------------------------------------------------

def _interpret_pigmentation(
    variants: Dict[str, Tuple[str, str, str]],
) -> Dict[str, Dict[str, str]]:
    """Interpret key pigmentation SNPs that define ancient population appearance."""
    results: Dict[str, Dict[str, str]] = {}

    # SLC24A5 skin color (rs1426654): G=light, A=dark/ancestral
    variant_data = variants.get("rs1426654")
    if variant_data is not None:
        _c, _p, geno = variant_data
        if geno and geno not in ("--", "00", "NN"):
            g_count = geno.upper().count("G")
            if g_count == 2:
                results["skin"] = {
                    "status": "light",
                    "detail": "GG -- homozygous light skin (European/Asian allele)",
                }
            elif g_count == 1:
                results["skin"] = {
                    "status": "mixed",
                    "detail": "AG -- heterozygous; one light, one ancestral allele",
                }
            else:
                results["skin"] = {
                    "status": "dark",
                    "detail": "AA -- ancestral dark skin allele (rare in modern Europeans)",
                }

    # HERC2/OCA2 eye color (rs12913832): G=blue, A=brown
    variant_data = variants.get("rs12913832")
    if variant_data is not None:
        _c, _p, geno = variant_data
        if geno and geno not in ("--", "00", "NN"):
            g_count = geno.upper().count("G")
            if g_count == 2:
                results["eyes"] = {
                    "status": "blue",
                    "detail": "GG -- strongly blue/light eyes",
                }
            elif g_count == 1:
                results["eyes"] = {
                    "status": "mixed",
                    "detail": "AG -- likely hazel/green/mixed",
                }
            else:
                results["eyes"] = {
                    "status": "brown",
                    "detail": "AA -- brown eyes",
                }

    # LCT lactase (rs4988235): G=persistent, A=non-persistent
    variant_data = variants.get("rs4988235")
    if variant_data is not None:
        _c, _p, geno = variant_data
        if geno and geno not in ("--", "00", "NN"):
            g_count = geno.upper().count("G")
            if g_count >= 1:
                geno_label = "GG" if g_count == 2 else "AG"
                results["lactase"] = {
                    "status": "persistent",
                    "detail": f"{geno_label} -- can digest milk as adult",
                }
            else:
                results["lactase"] = {
                    "status": "non-persistent",
                    "detail": "AA -- ancestral; likely lactose intolerance in adults",
                }

    return results


# ---------------------------------------------------------------------------
# PUBLIC API
# ---------------------------------------------------------------------------

def analyze_ancient(
    variants: Dict[str, Tuple[str, str, str]],
) -> Dict[str, Any]:
    """
    Analyze variants for ancient population affinity.

    Args:
        variants: Dict mapping rsID (lowercase) -> (chromosome, position, genotype)

    Returns:
        Dict with proportions, CI, pigmentation, used_snps, residual.
    """
    props, used_snps, resid = _estimate_proportions(variants)
    ci_low, ci_high = _bootstrap_ci(variants)
    pigmentation = _interpret_pigmentation(variants)

    return {
        "proportions": props,
        "ci_low": ci_low,
        "ci_high": ci_high,
        "used_snps": used_snps,
        "residual": resid,
        "pigmentation": pigmentation,
    }


def generate_ancient_json(result: Dict[str, Any]) -> Dict[str, Any]:
    """
    Convert raw analysis result into the structured JSON for the frontend.

    Args:
        result: Output from analyze_ancient()

    Returns:
        Structured JSON dict matching the report schema.
    """
    props = result["proportions"]
    ci_low = result["ci_low"]
    ci_high = result["ci_high"]
    used_snps = result["used_snps"]
    resid = result["residual"]
    pigmentation = result["pigmentation"]

    panel_size = len(AIMS_ANCIENT)
    aims_genotyped = len(used_snps)

    # Build sorted population list
    sorted_pops = sorted(ANCIENT_POPS, key=lambda p: -props[p])
    populations = []
    for pop in sorted_pops:
        idx = ANCIENT_POPS.index(pop)
        proportion = round(props[pop], 4)
        entry = {
            "code": pop,
            "label": ANCIENT_POP_LABELS[pop],
            "period": ANCIENT_POP_PERIODS[pop],
            "proportion": proportion,
            "percentage": round(proportion * 100, 1),
            "ciLow": round(float(ci_low[idx]) * 100, 1) if ci_low is not None else None,
            "ciHigh": round(float(ci_high[idx]) * 100, 1) if ci_high is not None else None,
            "description": ANCIENT_POP_DESCRIPTIONS[pop],
            "color": ANCIENT_POP_COLORS[pop],
        }
        populations.append(entry)

    # Three-wave European model
    whg_p = props.get("WHG", 0)
    eef_p = props.get("EEF", 0)
    steppe_p = props.get("Steppe", 0)
    chg_p = props.get("CHG", 0)
    ehg_p = props.get("EHG", 0)
    natufian_p = props.get("Natufian", 0)

    three_wave = {
        "hunterGatherer": round(whg_p + ehg_p, 4),
        "anatolianFarmer": round(eef_p + natufian_p, 4),
        "steppePastoralist": round(steppe_p + chg_p, 4),
    }

    return {
        "summary": {
            "panelSize": panel_size,
            "aimsGenotyped": aims_genotyped,
            "coverage": round(aims_genotyped / panel_size * 100, 1),
            "residual": round(resid, 4),
        },
        "populations": populations,
        "pigmentation": pigmentation,
        "threeWaveModel": three_wave,
    }
