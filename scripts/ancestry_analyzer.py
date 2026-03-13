"""
Ancestry Admixture Analyzer
Estimates admixture proportions using Ancestry-Informative Markers (AIMs)
with Weighted NNLS (Non-Negative Least Squares).
Returns structured JSON for frontend visualization.

Based on:
- gnomAD v3.1.2 verified population frequencies
- 1000 Genomes Phase 3 superpopulation data
- Behar et al. 2010, Price et al. 2008 (Ashkenazi)
- HGDP (Rosenberg 2002, Li 2008) for Middle Eastern
"""

from typing import Dict, Tuple, Any, List, Optional
import numpy as np
from scipy.optimize import nnls

# ─────────────────────────────────────────────────────────────────────────────
# POPULATION DEFINITIONS
# ─────────────────────────────────────────────────────────────────────────────

POPS = ["EUR_N", "EUR_S", "EUR_E", "SAS", "MID", "EAS", "AMR", "AFR", "ASH"]

POP_LABELS = {
    "EUR_N": "Northern European",
    "EUR_S": "Mediterranean / Southern European",
    "EUR_E": "Eastern European",
    "SAS":   "South Asian",
    "MID":   "Middle Eastern",
    "EAS":   "East Asian",
    "AMR":   "Native American",
    "AFR":   "Sub-Saharan African",
    "ASH":   "Ashkenazi Jewish",
}

POP_COLORS = {
    "EUR_N": "#4A90D9",
    "EUR_S": "#E8A838",
    "EUR_E": "#7B68EE",
    "SAS":   "#2ECC71",
    "MID":   "#E74C3C",
    "EAS":   "#F39C12",
    "AMR":   "#8B4513",
    "AFR":   "#1ABC9C",
    "ASH":   "#9B59B6",
}

POP_NARRATIVES = {
    "EUR_N": (
        "Northern European ancestry traces to ancient Scandinavian and northern "
        "Germanic peoples. Defined by high frequencies of lactase persistence, "
        "blue-eye alleles, and Northern European skin-lightening variants. "
        "In Latin America, arrived through 19th-20th century immigration: "
        "Germans to southern Brazil, British to Argentina, Scandinavians to Uruguay."
    ),
    "EUR_S": (
        "Mediterranean/Southern European ancestry descends from Early European "
        "Farmers who migrated from Anatolia ~7,000-9,000 years ago. The dominant "
        "European component in Latin America, arriving through Spanish and Portuguese "
        "colonization from 1500 onward, and later through mass Italian immigration."
    ),
    "EUR_E": (
        "Eastern European ancestry reflects populations of Central and Eastern "
        "Europe — a blend of Early European Farmers, Western Hunter-Gatherers, "
        "and significant Bronze Age Yamnaya steppe ancestry. In Latin America, "
        "significant immigration of Poles and Ukrainians to southern Brazil and Argentina."
    ),
    "SAS": (
        "South Asian ancestry represents the genetic legacy of the ancient Indus "
        "Valley Civilization and Indo-Aryan migrations from the Central Asian Steppe "
        "around 2000-1500 BCE. Present in Latin America primarily through 19th century "
        "Indian indentured labor in Trinidad, Guyana, and Suriname."
    ),
    "MID": (
        "Middle Eastern ancestry traces to the ancient populations of the Fertile "
        "Crescent — the cradle of civilization where agriculture first developed "
        "~10,000 BCE. In Latin America, substantial Lebanese, Syrian, and Palestinian "
        "immigration from ~1880-1950."
    ),
    "EAS": (
        "East Asian ancestry traces to ancient populations of China, Japan, and Korea. "
        "Characterized by the EDAR hair thickness variant and ALDH2 alcohol flush variant. "
        "In Latin America, arrived through Chinese workers to Cuba/Peru and Japanese "
        "immigration to Brazil (the largest Japanese diaspora outside Japan)."
    ),
    "AMR": (
        "Native American ancestry traces to the first peoples of the Americas, who "
        "crossed the Bering Land Bridge from Siberia ~15,000-20,000 years ago. "
        "They gave rise to the Maya, Aztec, and Inca civilizations. In modern "
        "Latin Americans, this reflects indigenous heritage."
    ),
    "AFR": (
        "Sub-Saharan African ancestry represents humanity's oldest genetic lineages. "
        "In Latin America, comes primarily from the transatlantic slave trade "
        "(1500-1850 CE). Brazil received the largest number of enslaved Africans "
        "in the Americas (~4-5 million)."
    ),
    "ASH": (
        "Ashkenazi Jewish ancestry is genetically distinctive, shaped by ancient "
        "Middle Eastern origin, substantial Southern European admixture, and "
        "~600-800 years of endogamy in Central/Eastern Europe. In Latin America, "
        "significant immigration to Argentina, Brazil, Uruguay, and Mexico from ~1900-1940."
    ),
}

# ─────────────────────────────────────────────────────────────────────────────
# VERIFIED AIM REFERENCE PANEL (gnomAD v3.1.2 / 1000 Genomes / literature)
# ─────────────────────────────────────────────────────────────────────────────

AIMS = {
    "rs2814778": {
        "ref": "T", "alt": "C", "gene": "DARC/ACKR1", "weight": 4.0,
        "EUR_N": 0.00, "EUR_S": 0.01, "EUR_E": 0.01, "SAS": 0.01,
        "MID":  0.06,  "EAS": 0.00,  "AMR":  0.05,  "AFR": 0.97, "ASH": 0.01,
    },
    "rs3827760": {
        "ref": "A", "alt": "G", "gene": "EDAR", "weight": 4.0,
        "EUR_N": 0.01, "EUR_S": 0.01, "EUR_E": 0.01, "SAS": 0.07,
        "MID":  0.02,  "EAS": 0.93,  "AMR":  0.70,  "AFR": 0.00, "ASH": 0.02,
    },
    "rs671": {
        "ref": "G", "alt": "A", "gene": "ALDH2", "weight": 4.0,
        "EUR_N": 0.00, "EUR_S": 0.00, "EUR_E": 0.00, "SAS": 0.00,
        "MID":  0.00,  "EAS": 0.28,  "AMR":  0.05,  "AFR": 0.00, "ASH": 0.00,
    },
    "rs1229984": {
        "ref": "G", "alt": "A", "gene": "ADH1B", "weight": 3.5,
        "EUR_N": 0.04, "EUR_S": 0.05, "EUR_E": 0.04, "SAS": 0.12,
        "MID":  0.16,  "EAS": 0.72,  "AMR":  0.28,  "AFR": 0.02, "ASH": 0.14,
    },
    "rs17822931": {
        "ref": "C", "alt": "T", "gene": "ABCC11", "weight": 3.5,
        "EUR_N": 0.02, "EUR_S": 0.02, "EUR_E": 0.02, "SAS": 0.10,
        "MID":  0.05,  "EAS": 0.75,  "AMR":  0.50,  "AFR": 0.01, "ASH": 0.03,
    },
    "rs1426654": {
        "ref": "A", "alt": "G", "gene": "SLC24A5", "weight": 2.0,
        "EUR_N": 0.99, "EUR_S": 0.99, "EUR_E": 0.99, "SAS": 0.95,
        "MID":  0.97,  "EAS": 0.01,  "AMR":  0.62,  "AFR": 0.02, "ASH": 0.99,
    },
    "rs16891982": {
        "ref": "C", "alt": "G", "gene": "SLC45A2", "weight": 2.5,
        "EUR_N": 0.97, "EUR_S": 0.71, "EUR_E": 0.87, "SAS": 0.24,
        "MID":  0.30,  "EAS": 0.02,  "AMR":  0.41,  "AFR": 0.01, "ASH": 0.80,
    },
    "rs12913832": {
        "ref": "A", "alt": "G", "gene": "HERC2/OCA2", "weight": 2.5,
        "EUR_N": 0.78, "EUR_S": 0.29, "EUR_E": 0.52, "SAS": 0.13,
        "MID":  0.10,  "EAS": 0.06,  "AMR":  0.12,  "AFR": 0.04, "ASH": 0.33,
    },
    "rs4988235": {
        "ref": "A", "alt": "G", "gene": "MCM6/LCT", "weight": 2.5,
        "EUR_N": 0.74, "EUR_S": 0.34, "EUR_E": 0.42, "SAS": 0.27,
        "MID":  0.20,  "EAS": 0.02,  "AMR":  0.18,  "AFR": 0.11, "ASH": 0.52,
    },
    "rs12203592": {
        "ref": "C", "alt": "T", "gene": "IRF4", "weight": 2.0,
        "EUR_N": 0.21, "EUR_S": 0.23, "EUR_E": 0.21, "SAS": 0.05,
        "MID":  0.07,  "EAS": 0.00,  "AMR":  0.07,  "AFR": 0.00, "ASH": 0.21,
    },
    "rs1805007": {
        "ref": "C", "alt": "T", "gene": "MC1R", "weight": 2.0,
        "EUR_N": 0.11, "EUR_S": 0.06, "EUR_E": 0.08, "SAS": 0.01,
        "MID":  0.01,  "EAS": 0.00,  "AMR":  0.02,  "AFR": 0.00, "ASH": 0.07,
    },
    "rs1800562": {
        "ref": "G", "alt": "A", "gene": "HFE", "weight": 2.5,
        "EUR_N": 0.06, "EUR_S": 0.02, "EUR_E": 0.03, "SAS": 0.00,
        "MID":  0.01,  "EAS": 0.00,  "AMR":  0.01,  "AFR": 0.00, "ASH": 0.02,
    },
    "rs1799945": {
        "ref": "C", "alt": "G", "gene": "HFE", "weight": 2.0,
        "EUR_N": 0.14, "EUR_S": 0.11, "EUR_E": 0.13, "SAS": 0.03,
        "MID":  0.05,  "EAS": 0.01,  "AMR":  0.04,  "AFR": 0.01, "ASH": 0.09,
    },
    "rs1042602": {
        "ref": "C", "alt": "A", "gene": "TYR", "weight": 2.0,
        "EUR_N": 0.44, "EUR_S": 0.38, "EUR_E": 0.41, "SAS": 0.20,
        "MID":  0.22,  "EAS": 0.02,  "AMR":  0.14,  "AFR": 0.03, "ASH": 0.38,
    },
    "rs1800414": {
        "ref": "C", "alt": "T", "gene": "OCA2", "weight": 0.4,
        "EUR_N": 0.42, "EUR_S": 0.20, "EUR_E": 0.33, "SAS": 0.09,
        "MID":  0.14,  "EAS": 0.70,  "AMR":  0.58,  "AFR": 0.16, "ASH": 0.23,
    },
    "rs4918664": {
        "ref": "A", "alt": "G", "gene": "CYP2C19", "weight": 2.0,
        "EUR_N": 0.12, "EUR_S": 0.10, "EUR_E": 0.11, "SAS": 0.39,
        "MID":  0.35,  "EAS": 0.87,  "AMR":  0.40,  "AFR": 0.02, "ASH": 0.20,
    },
    "rs1800497": {
        "ref": "G", "alt": "A", "gene": "ANKK1/DRD2", "weight": 1.5,
        "EUR_N": 0.22, "EUR_S": 0.17, "EUR_E": 0.19, "SAS": 0.31,
        "MID":  0.28,  "EAS": 0.41,  "AMR":  0.31,  "AFR": 0.39, "ASH": 0.19,
    },
    "rs1799971": {
        "ref": "A", "alt": "G", "gene": "OPRM1", "weight": 1.5,
        "EUR_N": 0.13, "EUR_S": 0.09, "EUR_E": 0.11, "SAS": 0.36,
        "MID":  0.18,  "EAS": 0.40,  "AMR":  0.20,  "AFR": 0.04, "ASH": 0.11,
    },
    "rs7554936": {
        "ref": "C", "alt": "T", "gene": "Unknown", "weight": 1.0,
        "EUR_N": 0.67, "EUR_S": 0.66, "EUR_E": 0.67, "SAS": 0.78,
        "MID":  0.65,  "EAS": 0.86,  "AMR":  0.68,  "AFR": 0.03, "ASH": 0.67,
    },
    "rs8050136": {
        "ref": "C", "alt": "A", "gene": "FTO", "weight": 1.5,
        "EUR_N": 0.43, "EUR_S": 0.38, "EUR_E": 0.40, "SAS": 0.32,
        "MID":  0.36,  "EAS": 0.12,  "AMR":  0.23,  "AFR": 0.46, "ASH": 0.40,
    },
    "rs4680": {
        "ref": "G", "alt": "A", "gene": "COMT", "weight": 1.0,
        "EUR_N": 0.48, "EUR_S": 0.44, "EUR_E": 0.46, "SAS": 0.36,
        "MID":  0.38,  "EAS": 0.28,  "AMR":  0.33,  "AFR": 0.23, "ASH": 0.46,
    },
    "rs1726866": {
        "ref": "G", "alt": "A", "gene": "TAS2R38", "weight": 1.0,
        "EUR_N": 0.46, "EUR_S": 0.41, "EUR_E": 0.44, "SAS": 0.36,
        "MID":  0.40,  "EAS": 0.28,  "AMR":  0.33,  "AFR": 0.24, "ASH": 0.43,
    },
    "rs713598": {
        "ref": "G", "alt": "C", "gene": "TAS2R38", "weight": 1.0,
        "EUR_N": 0.44, "EUR_S": 0.39, "EUR_E": 0.42, "SAS": 0.33,
        "MID":  0.37,  "EAS": 0.26,  "AMR":  0.31,  "AFR": 0.21, "ASH": 0.41,
    },
    "rs1801133": {
        "ref": "G", "alt": "A", "gene": "MTHFR", "weight": 1.0,
        "EUR_N": 0.33, "EUR_S": 0.42, "EUR_E": 0.37, "SAS": 0.26,
        "MID":  0.32,  "EAS": 0.23,  "AMR":  0.28,  "AFR": 0.10, "ASH": 0.37,
    },
    "rs1801131": {
        "ref": "A", "alt": "C", "gene": "MTHFR", "weight": 1.0,
        "EUR_N": 0.30, "EUR_S": 0.26, "EUR_E": 0.28, "SAS": 0.23,
        "MID":  0.21,  "EAS": 0.16,  "AMR":  0.20,  "AFR": 0.10, "ASH": 0.28,
    },
    "rs1800629": {
        "ref": "G", "alt": "A", "gene": "TNF", "weight": 1.0,
        "EUR_N": 0.13, "EUR_S": 0.16, "EUR_E": 0.14, "SAS": 0.12,
        "MID":  0.24,  "EAS": 0.07,  "AMR":  0.11,  "AFR": 0.20, "ASH": 0.22,
    },
    "rs4129267": {
        "ref": "T", "alt": "C", "gene": "IL6R", "weight": 1.0,
        "EUR_N": 0.36, "EUR_S": 0.40, "EUR_E": 0.38, "SAS": 0.56,
        "MID":  0.52,  "EAS": 0.33,  "AMR":  0.38,  "AFR": 0.18, "ASH": 0.50,
    },
    "rs2435357": {
        "ref": "G", "alt": "C", "gene": "RET", "weight": 1.0,
        "EUR_N": 0.26, "EUR_S": 0.24, "EUR_E": 0.25, "SAS": 0.56,
        "MID":  0.46,  "EAS": 0.30,  "AMR":  0.29,  "AFR": 0.09, "ASH": 0.36,
    },
    "rs2228479": {
        "ref": "A", "alt": "G", "gene": "MC1R", "weight": 1.0,
        "EUR_N": 0.08, "EUR_S": 0.05, "EUR_E": 0.06, "SAS": 0.24,
        "MID":  0.19,  "EAS": 0.16,  "AMR":  0.11,  "AFR": 0.05, "ASH": 0.11,
    },
    "rs2076533": {
        "ref": "C", "alt": "A", "gene": "HLA-B region", "weight": 1.5,
        "EUR_N": 0.14, "EUR_S": 0.21, "EUR_E": 0.17, "SAS": 0.24,
        "MID":  0.44,  "EAS": 0.09,  "AMR":  0.14,  "AFR": 0.07, "ASH": 0.50,
    },
    "rs9267531": {
        "ref": "G", "alt": "A", "gene": "HLA region", "weight": 1.5,
        "EUR_N": 0.19, "EUR_S": 0.24, "EUR_E": 0.21, "SAS": 0.29,
        "MID":  0.47,  "EAS": 0.11,  "AMR":  0.17,  "AFR": 0.09, "ASH": 0.45,
    },
    "rs6090989": {
        "ref": "G", "alt": "A", "gene": "Unknown", "weight": 1.2,
        "EUR_N": 0.17, "EUR_S": 0.22, "EUR_E": 0.19, "SAS": 0.14,
        "MID":  0.40,  "EAS": 0.07,  "AMR":  0.11,  "AFR": 0.04, "ASH": 0.55,
    },
    "rs3135027": {
        "ref": "G", "alt": "A", "gene": "HLA region", "weight": 1.0,
        "EUR_N": 0.21, "EUR_S": 0.28, "EUR_E": 0.24, "SAS": 0.33,
        "MID":  0.40,  "EAS": 0.14,  "AMR":  0.19,  "AFR": 0.11, "ASH": 0.42,
    },
    "rs10497520": {
        "ref": "C", "alt": "T", "gene": "Unknown", "weight": 0.8,
        "EUR_N": 0.40, "EUR_S": 0.26, "EUR_E": 0.33, "SAS": 0.17,
        "MID":  0.19,  "EAS": 0.08,  "AMR":  0.14,  "AFR": 0.05, "ASH": 0.30,
    },
    "rs1540771": {
        "ref": "A", "alt": "G", "gene": "Unknown", "weight": 0.8,
        "EUR_N": 0.63, "EUR_S": 0.55, "EUR_E": 0.59, "SAS": 0.42,
        "MID":  0.45,  "EAS": 0.20,  "AMR":  0.31,  "AFR": 0.12, "ASH": 0.56,
    },
    "rs1937845": {
        "ref": "A", "alt": "G", "gene": "FGFR3", "weight": 0.8,
        "EUR_N": 0.14, "EUR_S": 0.17, "EUR_E": 0.16, "SAS": 0.47,
        "MID":  0.40,  "EAS": 0.22,  "AMR":  0.26,  "AFR": 0.18, "ASH": 0.32,
    },
    "rs2073711": {
        "ref": "C", "alt": "T", "gene": "Unknown", "weight": 0.8,
        "EUR_N": 0.19, "EUR_S": 0.24, "EUR_E": 0.21, "SAS": 0.57,
        "MID":  0.47,  "EAS": 0.31,  "AMR":  0.36,  "AFR": 0.27, "ASH": 0.39,
    },
    "rs7149477": {
        "ref": "A", "alt": "G", "gene": "Unknown", "weight": 0.8,
        "EUR_N": 0.17, "EUR_S": 0.31, "EUR_E": 0.23, "SAS": 0.37,
        "MID":  0.44,  "EAS": 0.19,  "AMR":  0.24,  "AFR": 0.11, "ASH": 0.40,
    },
}


def _get_dosage(genotype: str, alt: str) -> Optional[int]:
    """Convert genotype string to alt allele dosage (0, 1, 2)."""
    if not genotype or genotype in ('--', '00', 'NN', 'NC', 'DI', 'II', 'DD'):
        return None
    g = genotype.upper().strip()
    a = alt.upper()
    cnt = g.count(a)
    if cnt > 2 or len(g) < 2:
        return None
    return cnt


def _detect_discordant(dosages: Dict[str, Optional[int]]) -> set:
    """Detect markers strongly inconsistent with majority signal."""
    flagged = set()
    slc24a5 = dosages.get("rs1426654")
    if slc24a5 == 0:
        european_evidence = []
        for em in ["rs12913832", "rs4988235", "rs16891982", "rs1042602", "rs12203592"]:
            d = dosages.get(em)
            if d is not None:
                european_evidence.append(d / 2.0)
        if european_evidence and np.mean(european_evidence) > 0.35:
            flagged.add("rs1426654")
    return flagged


def analyze_ancestry(variants: Dict[str, Tuple[str, str, str]]) -> Dict[str, Any]:
    """
    Analyze variants for ancestry admixture proportions.

    Args:
        variants: Dict mapping rsID -> (chromosome, position, genotype)

    Returns:
        Structured JSON dict for the ancestry report
    """
    # Extract dosages from parsed variants
    dosages: Dict[str, Optional[int]] = {}
    for rsid, info in AIMS.items():
        rsid_lower = rsid.lower()
        variant_data = variants.get(rsid_lower)
        if variant_data is None:
            dosages[rsid] = None
        else:
            _chrom, _pos, genotype = variant_data
            dosages[rsid] = _get_dosage(genotype, info["alt"])

    # Detect discordant markers
    discordant = _detect_discordant(dosages)

    # Build reference matrix
    rsids = list(AIMS.keys())
    weights = np.array([AIMS[r].get("weight", 1.0) for r in rsids])
    A = np.array([[AIMS[r].get(p, 0.0) for p in POPS] for r in rsids])

    # Find valid markers
    valid = [(i, rsids[i]) for i in range(len(rsids)) if dosages.get(rsids[i]) is not None]
    used_rsids = [r for _, r in valid]
    missing_rsids = [r for r in rsids if dosages.get(r) is None]

    if len(valid) < 5:
        # Not enough data for meaningful estimation
        equal_prop = round(1.0 / len(POPS), 4)
        proportions = {p: equal_prop for p in POPS}
        return _build_result(proportions, None, None, used_rsids, missing_rsids, dosages, discordant, reliable=False)

    # Build weighted system
    idx = [i for i, _ in valid]
    f = np.array([dosages[r] / 2.0 for _, r in valid])
    w = weights[idx].copy()

    # Down-weight discordant markers
    for k, rsid in enumerate(used_rsids):
        if rsid in discordant:
            w[k] *= 0.02

    # Weighted NNLS with Tikhonov regularization
    reg = 0.002
    A_w = A[idx, :] * w[:, np.newaxis]
    f_w = f * w
    A_aug = np.vstack([A_w, np.eye(len(POPS)) * reg])
    f_aug = np.concatenate([f_w, np.zeros(len(POPS))])

    x, resid = nnls(A_aug, f_aug)
    total = x.sum()
    if total > 0:
        x /= total
    else:
        x = np.ones(len(POPS)) / len(POPS)

    proportions = {p: round(float(x[j]), 4) for j, p in enumerate(POPS)}

    # Bootstrap confidence intervals
    ci_low, ci_high = _bootstrap_ci(dosages, A, rsids, weights, discordant)

    return _build_result(proportions, ci_low, ci_high, used_rsids, missing_rsids, dosages, discordant, reliable=True)


def _bootstrap_ci(dosages, A, rsids, weights, discordant, n=300, seed=42):
    """Bootstrap resampling for 95% confidence intervals."""
    valid = [i for i, r in enumerate(rsids) if dosages.get(r) is not None]
    if len(valid) < 8:
        return None, None

    eff_w = weights.copy()
    for i, r in enumerate(rsids):
        if r in discordant:
            eff_w[i] *= 0.02

    rng = np.random.default_rng(seed)
    all_x = []
    for _ in range(n):
        boot = rng.choice(valid, size=len(valid), replace=True)
        w = eff_w[boot]
        A_b = A[boot, :] * w[:, np.newaxis]
        f_b = np.array([dosages[rsids[i]] / 2.0 for i in boot]) * w
        x, _ = nnls(A_b, f_b)
        s = x.sum()
        all_x.append(x / s if s > 0 else np.ones(len(POPS)) / len(POPS))

    all_x = np.array(all_x)
    return np.percentile(all_x, 2.5, axis=0), np.percentile(all_x, 97.5, axis=0)


def _build_result(proportions, ci_low, ci_high, used, missing, dosages, discordant, reliable):
    """Build the structured JSON result."""
    # Sort populations by proportion descending
    sorted_pops = sorted(proportions.items(), key=lambda x: x[1], reverse=True)

    populations = []
    for pop, prop in sorted_pops:
        ji = POPS.index(pop)
        entry = {
            "code": pop,
            "label": POP_LABELS[pop],
            "proportion": prop,
            "percentage": round(prop * 100, 1),
            "color": POP_COLORS[pop],
        }
        if ci_low is not None and ci_high is not None:
            entry["ciLow"] = round(float(ci_low[ji]) * 100, 1)
            entry["ciHigh"] = round(float(ci_high[ji]) * 100, 1)
        if prop >= 0.005:
            entry["narrative"] = POP_NARRATIVES[pop]
        populations.append(entry)

    # Determine reliability
    n_used = len(used)
    if n_used >= 30:
        reliability = "high"
        reliability_label = "Good AIM coverage; continental estimates reliable"
    elif n_used >= 15:
        reliability = "moderate"
        reliability_label = "Continental ancestry reliable; sub-continental is approximate"
    elif n_used >= 8:
        reliability = "low"
        reliability_label = "Broad continental direction only"
    else:
        reliability = "very_low"
        reliability_label = "Insufficient AIMs; results unreliable"

    return {
        "summary": {
            "panelSize": len(AIMS),
            "aimsGenotyped": n_used,
            "aimsMissing": len(missing),
            "coverage": round(n_used / len(AIMS) * 100, 1),
            "discordantMarkers": len(discordant),
            "reliability": reliability,
            "reliabilityLabel": reliability_label,
        },
        "populations": populations,
    }
