"""
Italian Ancient Ancestry Analyzer
Estimates genetic affinity to ancient Italian populations across archaeological
periods, using Ancestry-Informative Markers (AIMs) calibrated to published
Italian archaeogenomics data.

KEY REFERENCE:
  Antonio et al. 2019 (Science): "Ancient Rome: A genetic crossroads of Europe
  and the Mediterranean"
    - 127 ancient individuals from Italy spanning Mesolithic to Medieval period

SUPPORTING REFERENCES:
  - Mathieson et al. 2015 (Nature): European genetic transformations
  - Haak et al. 2015 (Nature): Bronze Age Steppe expansion
  - Olalde et al. 2019 (Science): Iberian and European Bronze Age
  - Lazaridis et al. 2017 (Nature): Ancient Anatolian/EEF ancestry
  - Posth et al. 2018 (Nature Communications): Southern European prehistory

Returns structured JSON for frontend visualization.
"""

from typing import Dict, Tuple, Any, Optional
import numpy as np
from scipy.optimize import nnls

# ---------------------------------------------------------------------------
# POPULATION DEFINITIONS
# ---------------------------------------------------------------------------

ITALIAN_POPS = [
    "Mesolithic_Italy",    # ~10000-6000 BCE: WHG-dominant
    "Neolithic_Italy",     # ~5500-3000 BCE: EEF-dominant (Anatolian farmers)
    "BronzeAge_Italy",     # ~1700-700 BCE: EEF + Steppe admixture arriving
    "Iron_Rome",           # ~700-0 BCE: Etruscan/Latin period; EEF+Steppe base
    "Imperial_Rome",       # ~0-400 CE: Peak cosmopolitan diversity (MID, EEF, Steppe)
    "Medieval_Italy",      # ~600-1200 CE: More local/Northern EUR profile
    "Modern_Italian",      # Contemporary Italian (TSI/Southern Italian reference)
]

ITALIAN_POP_LABELS = {
    "Mesolithic_Italy":  "Mesolithic Italy (~10000-6000 BCE): Hunter-Gatherers",
    "Neolithic_Italy":   "Neolithic Italy (~5500-3000 BCE): First Farmers",
    "BronzeAge_Italy":   "Bronze Age Italy (~1700-700 BCE): Indo-European Arrivals",
    "Iron_Rome":         "Iron Age / Early Rome (~700-0 BCE): Etruscans & Latins",
    "Imperial_Rome":     "Imperial Rome (~0-400 CE): Crossroads of the Ancient World",
    "Medieval_Italy":    "Medieval Italy (~600-1200 CE): Post-Roman Transformations",
    "Modern_Italian":    "Modern Italian (Tuscany/Southern Italy reference)",
}

ITALIAN_POP_COLORS = {
    "Mesolithic_Italy":  "#3B82F6",  # blue
    "Neolithic_Italy":   "#F59E0B",  # amber
    "BronzeAge_Italy":   "#F97316",  # orange
    "Iron_Rome":         "#EF4444",  # red
    "Imperial_Rome":     "#8B5CF6",  # violet
    "Medieval_Italy":    "#10B981",  # emerald
    "Modern_Italian":    "#06B6D4",  # cyan
}

ITALIAN_POP_DESCRIPTIONS = {
    "Mesolithic_Italy": (
        "Mesolithic Italians (10,000-6,000 BCE) were hunter-gatherers genetically "
        "similar to Western Hunter-Gatherers (WHG) across Europe. They inhabited "
        "the Italian peninsula after the Last Glacial Maximum, exploiting coastal "
        "and forest resources. Like other WHG, they had DARK SKIN (SLC24A5 "
        "ancestral allele) but likely blue eyes (HERC2 derived allele present) — "
        "an unusual combination from our modern perspective. They had NO lactase "
        "persistence. Their descendants largely disappeared genetically with the "
        "Neolithic farming revolution, contributing only ~20-30% to later Italian "
        "populations. Key sites: Grotta del Romito (Calabria), Arene Candide "
        "(Liguria), Riparo Villabruna (Veneto)."
    ),
    "Neolithic_Italy": (
        "Neolithic Italians (5,500-3,000 BCE) were the Early European Farmers "
        "(EEF) who arrived from Anatolia by sea, replacing most of the WHG "
        "population. They introduced agriculture, pottery (Impressed Ware culture), "
        "and permanent villages to Italy. Genetically, they were predominantly "
        "Anatolian Neolithic (~80-90%) with minor WHG admixture. They had LIGHT "
        "SKIN (SLC24A5 near-fixed from their Anatolian origin) but DARK/BROWN "
        "EYES — and NO lactase persistence. Otzi the Iceman (Alps, 3,300 BCE) "
        "is the most famous ancient Italian: essentially pure Neolithic EEF. "
        "Key sites: Grotta dell'Uzzo (Sicily), Arene Candide Neolithic layer."
    ),
    "BronzeAge_Italy": (
        "Bronze Age Italy (1,700-700 BCE) saw the arrival of Steppe/Indo-European "
        "ancestry, likely via the Bell Beaker and later Urnfield cultures. Antonio "
        "et al. 2019 found ~15-20% Steppe ancestry in Bronze Age Italian samples, "
        "rising from essentially zero in Neolithic. This Steppe wave brought "
        "Indo-European languages (the ancestor of Latin, Oscan, Umbrian), "
        "metallurgy, horse use, and rising lactase persistence. The Bronze Age "
        "cultures of Italy include the Terramare (Po Valley), Nuragic (Sardinia), "
        "and Proto-Villanovan — leading to the Iron Age cultures of Etruscans, "
        "Latins, and Italic peoples. Sardinians today remain closest to Neolithic "
        "EEF with minimal Bronze Age Steppe admixture, making them a unique "
        "genetic window into early European farmers."
    ),
    "Iron_Rome": (
        "Iron Age and Early Roman period (700-0 BCE) encompasses the Etruscan "
        "civilization (900-300 BCE) and the rise of Rome. Antonio et al. 2019 "
        "found that pre-Imperial Romans had an ancestry profile dominated by "
        "EEF + Bronze Age Steppe, with some Middle Eastern admixture — consistent "
        "with Etruscan origins being partly Anatolian (debated) and Rome's early "
        "trade with Greece and the Eastern Mediterranean. The Etruscans, long "
        "mysterious in origin, appear to have local Italian ancestry from Neolithic "
        "farmers with additional Eastern Mediterranean contacts. Latin and Italic "
        "tribes are genomically similar to modern Southern Italians with slightly "
        "more Steppe ancestry."
    ),
    "Imperial_Rome": (
        "Imperial Rome (0-400 CE) was, genetically, one of the most cosmopolitan "
        "cities in ancient history. Antonio et al. 2019's key finding: at the "
        "height of the Empire, the population of Rome itself showed MASSIVE "
        "Middle Eastern and Eastern Mediterranean ancestry — shifting from a "
        "typical South European profile to one resembling modern Syrians and "
        "Lebanese. This reflects historical reality: Rome was a city of migrants, "
        "with enslaved peoples, soldiers, traders, and free citizens arriving "
        "from across the Empire — North Africa, the Levant, Anatolia, Persia, "
        "and beyond. The catacombs of Rome have yielded DNA from people with "
        "ancestry spanning the entire ancient world. Many of Rome's most famous "
        "emperors had provincial origins (Trajan: Spanish, Hadrian: Spanish, "
        "Septimius Severus: North African, Philip the Arab: Syrian). "
        "Interpreting 'Roman' as an ethnic identity for this period is "
        "historically and genetically incorrect — it was a political identity."
    ),
    "Medieval_Italy": (
        "Medieval Italy (600-1200 CE) shows a genetic shift back toward a more "
        "local European profile, following the collapse of the Western Roman "
        "Empire. This reflects the demographic decline of Rome's cosmopolitan "
        "immigrant population, replaced by local Italian and Germanic (Lombard, "
        "Gothic) settlers. The Lombards who invaded northern Italy in 568 CE "
        "brought a Northern European (Scandinavian-adjacent) genetic component "
        "that persists slightly in Po Valley populations today. Byzantine "
        "control of southern Italy and Sicily maintained some Eastern Mediterranean "
        "connections. Norman conquest of Sicily (1071-1194 CE) brought Northern "
        "European ancestry to the south. Medieval Italian genetics thus represent "
        "a complex blend of surviving Imperial-era populations plus Germanic and "
        "Northern European admixture."
    ),
    "Modern_Italian": (
        "Modern Italians (Tuscany/Southern Italy reference: TSI from 1000 Genomes) "
        "represent the contemporary genetic profile of Italy. They are "
        "predominantly EEF (Anatolian Neolithic farmer) with Bronze Age Steppe "
        "admixture (~15-20%), residual WHG ancestry (~15%), and variable Middle "
        "Eastern/North African admixture. Southern Italians and Sicilians have "
        "more Middle Eastern ancestry (reflecting ancient Greek colonization, Arab "
        "rule of Sicily 831-1072 CE, and Norman-era Saracen contacts). Northern "
        "Italians have more Steppe and Central European ancestry. Sardinians remain "
        "the most EEF-like population in Europe, essentially a genetic snapshot "
        "of the Neolithic farmers, and are often used as the EEF reference "
        "population in archaeogenomics."
    ),
}

# ---------------------------------------------------------------------------
# AIM ALLELE FREQUENCIES FOR ANCIENT ITALIAN PERIODS
# Frequencies derived from:
#   - Antonio et al. 2019 (Science) mean frequencies across samples per period
#   - Supplemented with Mathieson 2015, Haak 2015 for early periods
#   - gnomAD TSI for modern Italian reference
# For each rsid: (ref, alt, weight, {pop: alt_freq})
# ---------------------------------------------------------------------------

AIMS_ITALIAN = {
    # rsid: (ref, alt, weight, {pop: alt_freq})
    "rs2814778": ("T", "C", 4.0, {
        # Duffy-null (AFR): absent in ancient Italians; slight signal in Imperial Rome
        "Mesolithic_Italy": 0.00, "Neolithic_Italy": 0.00, "BronzeAge_Italy": 0.00,
        "Iron_Rome": 0.00, "Imperial_Rome": 0.02, "Medieval_Italy": 0.01,
        "Modern_Italian": 0.01,
    }),
    "rs3827760": ("A", "G", 4.0, {
        # EDAR (EAS): absent in Italy across all periods
        "Mesolithic_Italy": 0.00, "Neolithic_Italy": 0.00, "BronzeAge_Italy": 0.00,
        "Iron_Rome": 0.00, "Imperial_Rome": 0.01, "Medieval_Italy": 0.00,
        "Modern_Italian": 0.01,
    }),
    "rs671": ("G", "A", 4.0, {
        # ALDH2 (EAS): absent
        "Mesolithic_Italy": 0.00, "Neolithic_Italy": 0.00, "BronzeAge_Italy": 0.00,
        "Iron_Rome": 0.00, "Imperial_Rome": 0.00, "Medieval_Italy": 0.00,
        "Modern_Italian": 0.00,
    }),
    "rs1229984": ("G", "A", 3.5, {
        # ADH1B: low EUR ~4-5%
        "Mesolithic_Italy": 0.02, "Neolithic_Italy": 0.04, "BronzeAge_Italy": 0.05,
        "Iron_Rome": 0.05, "Imperial_Rome": 0.08, "Medieval_Italy": 0.05,
        "Modern_Italian": 0.05,
    }),
    "rs17822931": ("C", "T", 3.5, {
        # ABCC11 dry earwax: very low EUR
        "Mesolithic_Italy": 0.01, "Neolithic_Italy": 0.02, "BronzeAge_Italy": 0.02,
        "Iron_Rome": 0.02, "Imperial_Rome": 0.04, "Medieval_Italy": 0.02,
        "Modern_Italian": 0.02,
    }),
    "rs1426654": ("A", "G", 2.0, {
        # SLC24A5 skin lightening
        "Mesolithic_Italy": 0.12, "Neolithic_Italy": 0.96, "BronzeAge_Italy": 0.97,
        "Iron_Rome": 0.98, "Imperial_Rome": 0.96, "Medieval_Italy": 0.98,
        "Modern_Italian": 0.99,
    }),
    "rs16891982": ("C", "G", 2.5, {
        # SLC45A2: moderate in S.EUR; lower than N.EUR
        "Mesolithic_Italy": 0.30, "Neolithic_Italy": 0.50, "BronzeAge_Italy": 0.60,
        "Iron_Rome": 0.65, "Imperial_Rome": 0.55, "Medieval_Italy": 0.68,
        "Modern_Italian": 0.71,
    }),
    "rs12913832": ("A", "G", 2.5, {
        # HERC2 blue eyes: LOW in S.EUR/Mediterranean populations
        "Mesolithic_Italy": 0.68, "Neolithic_Italy": 0.22, "BronzeAge_Italy": 0.30,
        "Iron_Rome": 0.32, "Imperial_Rome": 0.20, "Medieval_Italy": 0.35,
        "Modern_Italian": 0.29,
    }),
    "rs4988235": ("A", "G", 2.5, {
        # LCT lactase persistence: LOW in ancient Italy, rising slowly
        "Mesolithic_Italy": 0.01, "Neolithic_Italy": 0.02, "BronzeAge_Italy": 0.08,
        "Iron_Rome": 0.14, "Imperial_Rome": 0.12, "Medieval_Italy": 0.22,
        "Modern_Italian": 0.34,
    }),
    "rs12203592": ("C", "T", 2.0, {
        # IRF4 freckles/light hair: ~20-23% EUR
        "Mesolithic_Italy": 0.22, "Neolithic_Italy": 0.20, "BronzeAge_Italy": 0.22,
        "Iron_Rome": 0.22, "Imperial_Rome": 0.18, "Medieval_Italy": 0.23,
        "Modern_Italian": 0.23,
    }),
    "rs1805007": ("C", "T", 2.0, {
        # MC1R R151C red hair: low in S.EUR (~6%)
        "Mesolithic_Italy": 0.07, "Neolithic_Italy": 0.05, "BronzeAge_Italy": 0.06,
        "Iron_Rome": 0.06, "Imperial_Rome": 0.05, "Medieval_Italy": 0.07,
        "Modern_Italian": 0.06,
    }),
    "rs1800562": ("G", "A", 2.5, {
        # HFE C282Y: very low in S.EUR (~1-2%); rises northward
        "Mesolithic_Italy": 0.01, "Neolithic_Italy": 0.01, "BronzeAge_Italy": 0.01,
        "Iron_Rome": 0.02, "Imperial_Rome": 0.01, "Medieval_Italy": 0.02,
        "Modern_Italian": 0.01,
    }),
    "rs1799945": ("C", "G", 2.0, {
        # HFE H63D: ~10-11% S.EUR
        "Mesolithic_Italy": 0.06, "Neolithic_Italy": 0.08, "BronzeAge_Italy": 0.10,
        "Iron_Rome": 0.10, "Imperial_Rome": 0.09, "Medieval_Italy": 0.11,
        "Modern_Italian": 0.11,
    }),
    "rs1042602": ("C", "A", 2.0, {
        # TYR pigmentation: ~38% TSI
        "Mesolithic_Italy": 0.22, "Neolithic_Italy": 0.32, "BronzeAge_Italy": 0.36,
        "Iron_Rome": 0.38, "Imperial_Rome": 0.32, "Medieval_Italy": 0.39,
        "Modern_Italian": 0.38,
    }),
    "rs1800414": ("C", "T", 0.4, {
        # OCA2 downweighted
        "Mesolithic_Italy": 0.28, "Neolithic_Italy": 0.20, "BronzeAge_Italy": 0.26,
        "Iron_Rome": 0.28, "Imperial_Rome": 0.24, "Medieval_Italy": 0.28,
        "Modern_Italian": 0.22,
    }),
    "rs4918664": ("A", "G", 2.0, {
        # CYP2C19: low EUR ~10-12%; higher in MID (~35%)
        "Mesolithic_Italy": 0.11, "Neolithic_Italy": 0.12, "BronzeAge_Italy": 0.12,
        "Iron_Rome": 0.14, "Imperial_Rome": 0.22, "Medieval_Italy": 0.14,
        "Modern_Italian": 0.10,
    }),
    "rs1800497": ("G", "A", 1.5, {
        # DRD2 Taq1A: ~17% S.EUR
        "Mesolithic_Italy": 0.20, "Neolithic_Italy": 0.18, "BronzeAge_Italy": 0.18,
        "Iron_Rome": 0.18, "Imperial_Rome": 0.22, "Medieval_Italy": 0.18,
        "Modern_Italian": 0.17,
    }),
    "rs1799971": ("A", "G", 1.5, {
        # OPRM1: ~9% S.EUR
        "Mesolithic_Italy": 0.10, "Neolithic_Italy": 0.10, "BronzeAge_Italy": 0.11,
        "Iron_Rome": 0.11, "Imperial_Rome": 0.16, "Medieval_Italy": 0.11,
        "Modern_Italian": 0.09,
    }),
    "rs7554936": ("C", "T", 1.0, {
        # not-AFR marker: high in all
        "Mesolithic_Italy": 0.65, "Neolithic_Italy": 0.66, "BronzeAge_Italy": 0.67,
        "Iron_Rome": 0.67, "Imperial_Rome": 0.64, "Medieval_Italy": 0.67,
        "Modern_Italian": 0.66,
    }),
    "rs8050136": ("C", "A", 1.5, {
        # FTO: ~38% S.EUR
        "Mesolithic_Italy": 0.38, "Neolithic_Italy": 0.40, "BronzeAge_Italy": 0.40,
        "Iron_Rome": 0.40, "Imperial_Rome": 0.38, "Medieval_Italy": 0.40,
        "Modern_Italian": 0.38,
    }),
    "rs4680": ("G", "A", 1.0, {
        # COMT: ~44% S.EUR
        "Mesolithic_Italy": 0.44, "Neolithic_Italy": 0.43, "BronzeAge_Italy": 0.44,
        "Iron_Rome": 0.44, "Imperial_Rome": 0.40, "Medieval_Italy": 0.44,
        "Modern_Italian": 0.44,
    }),
    "rs1726866": ("G", "A", 1.0, {
        "Mesolithic_Italy": 0.42, "Neolithic_Italy": 0.42, "BronzeAge_Italy": 0.42,
        "Iron_Rome": 0.43, "Imperial_Rome": 0.40, "Medieval_Italy": 0.43,
        "Modern_Italian": 0.41,
    }),
    "rs713598": ("G", "C", 1.0, {
        "Mesolithic_Italy": 0.40, "Neolithic_Italy": 0.40, "BronzeAge_Italy": 0.40,
        "Iron_Rome": 0.41, "Imperial_Rome": 0.38, "Medieval_Italy": 0.41,
        "Modern_Italian": 0.39,
    }),
    "rs1801133": ("G", "A", 1.0, {
        # MTHFR C677T: higher in S.EUR (~42% TSI)
        "Mesolithic_Italy": 0.28, "Neolithic_Italy": 0.38, "BronzeAge_Italy": 0.40,
        "Iron_Rome": 0.41, "Imperial_Rome": 0.36, "Medieval_Italy": 0.40,
        "Modern_Italian": 0.42,
    }),
    "rs1801131": ("A", "C", 1.0, {
        # MTHFR A1298C
        "Mesolithic_Italy": 0.22, "Neolithic_Italy": 0.24, "BronzeAge_Italy": 0.26,
        "Iron_Rome": 0.27, "Imperial_Rome": 0.25, "Medieval_Italy": 0.27,
        "Modern_Italian": 0.26,
    }),
    "rs1800629": ("G", "A", 1.0, {
        # TNF: ~16% S.EUR; higher in MID (~24%)
        "Mesolithic_Italy": 0.12, "Neolithic_Italy": 0.17, "BronzeAge_Italy": 0.17,
        "Iron_Rome": 0.18, "Imperial_Rome": 0.22, "Medieval_Italy": 0.17,
        "Modern_Italian": 0.16,
    }),
    "rs4129267": ("T", "C", 1.0, {
        # IL6R: ~40% S.EUR; higher in MID (~52%)
        "Mesolithic_Italy": 0.32, "Neolithic_Italy": 0.42, "BronzeAge_Italy": 0.42,
        "Iron_Rome": 0.43, "Imperial_Rome": 0.50, "Medieval_Italy": 0.42,
        "Modern_Italian": 0.40,
    }),
    "rs2435357": ("G", "C", 1.0, {
        # RET: ~24% EUR; higher MID ~46%
        "Mesolithic_Italy": 0.20, "Neolithic_Italy": 0.28, "BronzeAge_Italy": 0.28,
        "Iron_Rome": 0.30, "Imperial_Rome": 0.40, "Medieval_Italy": 0.28,
        "Modern_Italian": 0.24,
    }),
    "rs2228479": ("A", "G", 1.0, {
        # MC1R Val60Leu: ~5% S.EUR; higher SAS/MID ~19-24%
        "Mesolithic_Italy": 0.05, "Neolithic_Italy": 0.07, "BronzeAge_Italy": 0.08,
        "Iron_Rome": 0.09, "Imperial_Rome": 0.14, "Medieval_Italy": 0.08,
        "Modern_Italian": 0.05,
    }),
    "rs2076533": ("C", "A", 1.5, {
        # HLA-B region: MID+ASH enriched; slightly elevated in Imperial Rome
        "Mesolithic_Italy": 0.12, "Neolithic_Italy": 0.20, "BronzeAge_Italy": 0.21,
        "Iron_Rome": 0.24, "Imperial_Rome": 0.38, "Medieval_Italy": 0.22,
        "Modern_Italian": 0.21,
    }),
    "rs9267531": ("G", "A", 1.5, {
        # HLA region
        "Mesolithic_Italy": 0.15, "Neolithic_Italy": 0.24, "BronzeAge_Italy": 0.25,
        "Iron_Rome": 0.28, "Imperial_Rome": 0.42, "Medieval_Italy": 0.25,
        "Modern_Italian": 0.24,
    }),
    "rs6090989": ("G", "A", 1.2, {
        # MID+ASH enriched: diagnostic for Imperial Rome
        "Mesolithic_Italy": 0.10, "Neolithic_Italy": 0.18, "BronzeAge_Italy": 0.20,
        "Iron_Rome": 0.22, "Imperial_Rome": 0.38, "Medieval_Italy": 0.20,
        "Modern_Italian": 0.22,
    }),
    "rs3135027": ("G", "A", 1.0, {
        # HLA region
        "Mesolithic_Italy": 0.16, "Neolithic_Italy": 0.25, "BronzeAge_Italy": 0.27,
        "Iron_Rome": 0.30, "Imperial_Rome": 0.40, "Medieval_Italy": 0.27,
        "Modern_Italian": 0.28,
    }),
    "rs10497520": ("C", "T", 0.8, {
        # N-S EUR gradient: lower in S.EUR
        "Mesolithic_Italy": 0.35, "Neolithic_Italy": 0.26, "BronzeAge_Italy": 0.28,
        "Iron_Rome": 0.28, "Imperial_Rome": 0.24, "Medieval_Italy": 0.30,
        "Modern_Italian": 0.26,
    }),
    "rs1540771": ("A", "G", 0.8, {
        # EUR-enriched: ~55% S.EUR
        "Mesolithic_Italy": 0.58, "Neolithic_Italy": 0.55, "BronzeAge_Italy": 0.56,
        "Iron_Rome": 0.57, "Imperial_Rome": 0.52, "Medieval_Italy": 0.57,
        "Modern_Italian": 0.55,
    }),
    "rs1937845": ("A", "G", 0.8, {
        # FGFR3: ~17% S.EUR; higher MID ~40%
        "Mesolithic_Italy": 0.10, "Neolithic_Italy": 0.20, "BronzeAge_Italy": 0.22,
        "Iron_Rome": 0.24, "Imperial_Rome": 0.35, "Medieval_Italy": 0.22,
        "Modern_Italian": 0.17,
    }),
    "rs2073711": ("C", "T", 0.8, {
        # SAS+MID-enriched
        "Mesolithic_Italy": 0.20, "Neolithic_Italy": 0.35, "BronzeAge_Italy": 0.36,
        "Iron_Rome": 0.38, "Imperial_Rome": 0.46, "Medieval_Italy": 0.36,
        "Modern_Italian": 0.24,
    }),
    "rs7149477": ("A", "G", 0.8, {
        # S.EUR + MID-enriched; higher in Italy than N.EUR
        "Mesolithic_Italy": 0.15, "Neolithic_Italy": 0.30, "BronzeAge_Italy": 0.32,
        "Iron_Rome": 0.35, "Imperial_Rome": 0.44, "Medieval_Italy": 0.33,
        "Modern_Italian": 0.31,
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
# NNLS DECOMPOSITION
# ---------------------------------------------------------------------------

def _estimate_proportions(
    variants: Dict[str, Tuple[str, str, str]],
) -> Tuple[Dict[str, float], list, float]:
    """
    Run weighted NNLS decomposition into ancient Italian period components.

    Returns (proportions_dict, used_snps_list, residual).
    """
    n_pops = len(ITALIAN_POPS)
    rows_A = []
    rows_f = []
    rows_w = []
    used_snps = []

    for rsid, (ref, alt, weight, pop_freqs) in AIMS_ITALIAN.items():
        variant_data = variants.get(rsid.lower())
        if variant_data is None:
            continue
        _chrom, _pos, genotype = variant_data
        dosage = _get_dosage(genotype, alt)
        if dosage is None:
            continue
        obs_freq = dosage / 2.0
        row = [pop_freqs.get(p, 0.0) for p in ITALIAN_POPS]
        rows_A.append(row)
        rows_f.append(obs_freq)
        rows_w.append(weight)
        used_snps.append(rsid)

    if len(used_snps) < 5:
        return {p: 1.0 / n_pops for p in ITALIAN_POPS}, used_snps, 0.0

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

    props = {p: float(x[i]) for i, p in enumerate(ITALIAN_POPS)}
    return props, used_snps, float(resid)


# ---------------------------------------------------------------------------
# AFFINITY SCORES (min-max normalized distances)
# ---------------------------------------------------------------------------

def _compute_affinity_scores(
    variants: Dict[str, Tuple[str, str, str]],
) -> Dict[str, float]:
    """
    Compute relative affinity to each ancient Italian period using
    min-max normalization of weighted Euclidean distances.
    """
    distances = {}
    for pop in ITALIAN_POPS:
        wsd = 0.0
        wsum = 0.0
        for rsid, (ref, alt, weight, pop_freqs) in AIMS_ITALIAN.items():
            variant_data = variants.get(rsid.lower())
            if variant_data is None:
                continue
            _chrom, _pos, genotype = variant_data
            dosage = _get_dosage(genotype, alt)
            if dosage is None:
                continue
            obs = dosage / 2.0
            wsd += weight * (obs - pop_freqs[pop]) ** 2
            wsum += weight
        if wsum > 0:
            distances[pop] = np.sqrt(wsd / wsum)

    if not distances:
        return {p: 0.0 for p in ITALIAN_POPS}

    min_d = min(distances.values())
    max_d = max(distances.values())
    rng = max_d - min_d

    if rng > 0.001:
        return {p: (1.0 - (distances[p] - min_d) / rng) * 100 for p in distances}
    return {p: 100.0 / len(distances) for p in distances}


# ---------------------------------------------------------------------------
# KEY MARKERS
# ---------------------------------------------------------------------------

def _interpret_key_markers(
    variants: Dict[str, Tuple[str, str, str]],
) -> Dict[str, Dict[str, str]]:
    """Interpret key markers that distinguish ancient Italian periods."""
    results: Dict[str, Dict[str, str]] = {}

    # SLC24A5 skin color (rs1426654): G=light (EEF), A=dark (WHG)
    variant_data = variants.get("rs1426654")
    if variant_data is not None:
        _c, _p, geno = variant_data
        if geno and geno not in ("--", "00", "NN"):
            g_count = geno.upper().count("G")
            if g_count == 2:
                results["skin"] = {
                    "status": "light",
                    "detail": "GG -- characteristic of EEF Neolithic farmers; WHG was dark",
                }
            elif g_count == 1:
                results["skin"] = {
                    "status": "mixed",
                    "detail": "AG -- one ancestral/one derived allele",
                }
            else:
                results["skin"] = {
                    "status": "dark",
                    "detail": "AA -- ancestral/dark; extremely rare in modern Italians",
                }

    # HERC2 eye color (rs12913832): G=blue (WHG-like), A=brown (EEF/Mediterranean)
    variant_data = variants.get("rs12913832")
    if variant_data is not None:
        _c, _p, geno = variant_data
        if geno and geno not in ("--", "00", "NN"):
            g_count = geno.upper().count("G")
            if g_count == 2:
                results["eyes"] = {
                    "status": "blue",
                    "detail": "GG -- blue/light; WHG heritage or Northern European component",
                }
            elif g_count == 1:
                results["eyes"] = {
                    "status": "mixed",
                    "detail": "AG -- hazel/green/mixed (common in Italy)",
                }
            else:
                results["eyes"] = {
                    "status": "brown",
                    "detail": "AA -- brown; characteristic Mediterranean/EEF profile",
                }

    # LCT lactase (rs4988235): G=persistent (Steppe/Bronze Age), A=ancestral
    variant_data = variants.get("rs4988235")
    if variant_data is not None:
        _c, _p, geno = variant_data
        if geno and geno not in ("--", "00", "NN"):
            g_count = geno.upper().count("G")
            if g_count == 2:
                results["lactase"] = {
                    "status": "persistent",
                    "detail": "GG -- fully persistent; Northern/Steppe ancestry",
                }
            elif g_count == 1:
                results["lactase"] = {
                    "status": "partial",
                    "detail": "AG -- one LP allele; partial persistence",
                }
            else:
                results["lactase"] = {
                    "status": "non-persistent",
                    "detail": "AA -- ancestral non-persistent; typical ancient Mediterranean",
                }

    # MTHFR C677T (rs1801133): A allele elevated in Southern European/EEF
    variant_data = variants.get("rs1801133")
    if variant_data is not None:
        _c, _p, geno = variant_data
        if geno and geno not in ("--", "00", "NN"):
            a_count = geno.upper().count("A")
            if a_count > 0:
                geno_label = "AA" if a_count == 2 else "GA"
                results["mthfr"] = {
                    "status": "present",
                    "detail": f"{geno_label} -- MTHFR 677T allele present; elevated in Southern Europe/EEF",
                }
            else:
                results["mthfr"] = {
                    "status": "absent",
                    "detail": "GG -- ancestral; MTHFR C677T absent",
                }

    return results


# ---------------------------------------------------------------------------
# PUBLIC API
# ---------------------------------------------------------------------------

def analyze_italian(
    variants: Dict[str, Tuple[str, str, str]],
) -> Dict[str, Any]:
    """
    Analyze variants for ancient Italian population affinity.

    Args:
        variants: Dict mapping rsID (lowercase) -> (chromosome, position, genotype)

    Returns:
        Dict with proportions, affinity_scores, used_snps, residual, key_markers.
    """
    props, used_snps, resid = _estimate_proportions(variants)
    affinity_scores = _compute_affinity_scores(variants)
    key_markers = _interpret_key_markers(variants)

    return {
        "proportions": props,
        "affinity_scores": affinity_scores,
        "used_snps": used_snps,
        "residual": resid,
        "key_markers": key_markers,
    }


def generate_italian_json(result: Dict[str, Any]) -> Dict[str, Any]:
    """
    Convert raw analysis result into structured JSON for the frontend.

    Args:
        result: Output from analyze_italian()

    Returns:
        Structured JSON dict matching the Italian ancestry report schema.
    """
    props = result["proportions"]
    affinity_scores = result["affinity_scores"]
    used_snps = result["used_snps"]
    resid = result["residual"]
    key_markers = result["key_markers"]

    panel_size = len(AIMS_ITALIAN)
    snps_used = len(used_snps)

    # Build period list sorted by affinity score descending
    sorted_pops = sorted(ITALIAN_POPS, key=lambda p: -affinity_scores.get(p, 0))
    periods = []
    for pop in sorted_pops:
        proportion = round(props.get(pop, 0.0), 4)
        entry = {
            "code": pop,
            "label": ITALIAN_POP_LABELS[pop],
            "affinityScore": round(affinity_scores.get(pop, 0.0), 1),
            "proportion": proportion,
            "percentage": round(proportion * 100, 1),
            "description": ITALIAN_POP_DESCRIPTIONS[pop],
            "color": ITALIAN_POP_COLORS[pop],
        }
        periods.append(entry)

    # Top match is first in sorted list
    top_pop = sorted_pops[0]
    top_match = {
        "code": top_pop,
        "label": ITALIAN_POP_LABELS[top_pop],
        "description": ITALIAN_POP_DESCRIPTIONS[top_pop],
    }

    return {
        "summary": {
            "panelSize": panel_size,
            "snpsUsed": snps_used,
            "coverage": round(snps_used / panel_size * 100, 1),
        },
        "periods": periods,
        "topMatch": top_match,
        "keyMarkers": key_markers,
    }
