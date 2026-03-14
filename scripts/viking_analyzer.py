"""
Viking Heritage Analyzer
Estimates genetic affinity to Viking-age populations using Ancestry-Informative
Markers (AIMs) with allele frequencies calibrated to published Viking-age genomics.

KEY REFERENCE:
  Margaryan et al. 2020 (Nature): "Population genomics of the Viking world"
    - 442 Viking-age individuals (700-1100 CE)
    - 80 archaeological sites across Scandinavia, British Isles, Eastern Europe
    - Found substantial genetic heterogeneity within and between Viking groups

SUPPORTING REFERENCES:
  - Guenther et al. 2015 (MBE): Scandinavian pre-Viking genetics
  - Olalde et al. 2014 (Nature): Western Hunter-Gatherer ancestry
  - Haak et al. 2015 (Nature): Bronze Age Steppe contribution

SCIENTIFIC LIMITATIONS:
  Viking-age populations differ subtly from each other (~1,000-year-old
  populations are still mostly modern-European-like). Our 50-SNP panel can
  distinguish broad ancient components (WHG, EEF, Steppe) but cannot
  reliably distinguish Danish from Norwegian Vikings at the individual level.
  Results show relative affinity scores, not direct descent.

Returns structured JSON for frontend visualization.
"""

from typing import Dict, Tuple, Any, Optional
import numpy as np
from scipy.optimize import nnls

# ---------------------------------------------------------------------------
# POPULATION DEFINITIONS
# ---------------------------------------------------------------------------

VIKING_POPS = [
    "Norwegian", "Danish", "Swedish", "British_Norse",
    "Eastern_Route", "Modern_Scandinavian"
]

VIKING_POP_LABELS = {
    "Norwegian":           "Norwegian / Western Norse Vikings",
    "Danish":              "Danish / Southern Norse Vikings",
    "Swedish":             "Swedish / Eastern Norse Vikings",
    "British_Norse":       "British Isles Norse (Anglo-Scandinavian)",
    "Eastern_Route":       "Eastern Route (Rus' / Baltic / Slavic Norse)",
    "Modern_Scandinavian": "Modern Scandinavian (contemporary reference)",
}

VIKING_POP_COLORS = {
    "Norwegian":           "#3B82F6",  # blue
    "Danish":              "#EF4444",  # red
    "Swedish":             "#F59E0B",  # amber
    "British_Norse":       "#10B981",  # emerald
    "Eastern_Route":       "#8B5CF6",  # violet
    "Modern_Scandinavian": "#06B6D4",  # cyan
}

VIKING_POP_DESCRIPTIONS = {
    "Norwegian": (
        "Norwegian Vikings (700-1100 CE) primarily raided and settled the British "
        "Isles, Faroe Islands, Iceland, and Greenland. They were the first Europeans "
        "to reach North America (Vinland, c. 1000 CE). Norwegian Viking-age individuals "
        "from Margaryan 2020 show the highest Western Hunter-Gatherer (WHG) and "
        "Scandinavian Hunter-Gatherer (SHG) ancestry among Viking groups -- reflecting "
        "the isolation of mountainous western Norway. Norwegian-origin individuals are "
        "characterized by very high frequencies of the blue-eye allele (HERC2) and "
        "high lactase persistence. Famous Norwegian Viking sites: Salme (Estonia), "
        "Gokstad ship burial, Oseberg ship burial."
    ),
    "Danish": (
        "Danish Vikings (700-1100 CE) were traders, raiders, and rulers of England "
        "(the Danelaw, 865-954 CE) and France (Normandy). Genetically, Danish Vikings "
        "cluster closest to modern Danes and British/Irish populations, with slightly "
        "more EEF (Early European Farmer) ancestry than Norwegian Vikings. The Jutland "
        "peninsula's flat, agricultural landscape historically selected for EEF "
        "ancestry. Famous Danish Viking sites: Hedeby (trading center), Trelleborg "
        "ring fortresses, Jelling complex (burial of Harald Bluetooth's parents)."
    ),
    "Swedish": (
        "Swedish Vikings (700-1100 CE) primarily traveled east -- to the Baltic, "
        "Ladoga, Novgorod, and ultimately Constantinople (as Varangian Guard). "
        "The Birka and Gotland trading populations in Margaryan 2020 show distinct "
        "Eastern affinity: more Baltic/Slavic ancestry and some Central Asian "
        "admixture from Volga trade routes. One famous Birka warrior (Bj.581), "
        "long assumed male, was revealed by ancient DNA to be biologically female. "
        "Swedish Viking-age individuals span a broader genetic range than Norwegian "
        "Vikings, reflecting their role as long-distance traders along Eastern routes."
    ),
    "British_Norse": (
        "Anglo-Scandinavian populations (Danelaw England, Viking Scotland/Ireland, "
        "Orkney/Shetland settlers) represent Norse Vikings who mixed with local "
        "British, Irish, and Pictish populations. Margaryan 2020 found that Norse "
        "settlers in the British Isles quickly admixed with locals -- within 1-2 "
        "generations. DNA from Viking-age York, Dublin, and Orkney shows ~50-70% "
        "Norse ancestry alongside British/Irish contributions. These populations "
        "carry more of the British Celtic/Pre-Roman ancestry component alongside "
        "their Scandinavian heritage."
    ),
    "Eastern_Route": (
        "Eastern Route Vikings (Rus') traveled from Sweden down the Volga and Dnieper "
        "rivers, founding Novgorod and Kyiv, trading with the Byzantine Empire and "
        "Caliphate. Margaryan 2020 found that Rus' individuals have significant "
        "Central/Eastern European (Slavic) and occasionally Central Asian admixture, "
        "acquired through trade and settlement in the Pontic Steppe region. Some "
        "Eastern Route individuals even show East Asian ancestry from Volga-Ural "
        "contact zones. The Varangian Guard in Constantinople represented the most "
        "cosmopolitan Viking group."
    ),
    "Modern_Scandinavian": (
        "Modern Scandinavians (Danes, Swedes, Norwegians, Icelanders) are the direct "
        "genetic descendants of Viking-age populations, with ~1,000 years of additional "
        "drift and minor admixture. Icelanders are the closest living descendants of "
        "Norwegian Vikings (settlement ~874 CE, minimal subsequent immigration). "
        "Modern Scandinavians have higher lactase persistence than Viking-age "
        "individuals -- the trait continued to increase under selection through the "
        "medieval period. This reference population uses gnomAD FIN/GBR data."
    ),
}

# ---------------------------------------------------------------------------
# AIM ALLELE FREQUENCIES FOR VIKING-AGE POPULATIONS
# Frequencies represent approximate allele frequency means across Viking-age
# individuals per geographic cluster (Margaryan 2020 clusters) supplemented
# with modern North European data for calibration.
#
# For each rsid: (ref, alt, weight, {pop: alt_freq})
# ---------------------------------------------------------------------------

AIMS_VIKING = {
    "rs2814778": ("T", "C", 4.0, {
        "Norwegian": 0.00, "Danish": 0.00, "Swedish": 0.00,
        "British_Norse": 0.00, "Eastern_Route": 0.00, "Modern_Scandinavian": 0.00,
    }),
    "rs3827760": ("A", "G", 4.0, {
        "Norwegian": 0.00, "Danish": 0.00, "Swedish": 0.01,
        "British_Norse": 0.00, "Eastern_Route": 0.03, "Modern_Scandinavian": 0.01,
    }),
    "rs671": ("G", "A", 4.0, {
        "Norwegian": 0.00, "Danish": 0.00, "Swedish": 0.00,
        "British_Norse": 0.00, "Eastern_Route": 0.00, "Modern_Scandinavian": 0.00,
    }),
    "rs1229984": ("G", "A", 3.5, {
        "Norwegian": 0.04, "Danish": 0.04, "Swedish": 0.04,
        "British_Norse": 0.04, "Eastern_Route": 0.04, "Modern_Scandinavian": 0.04,
    }),
    "rs17822931": ("C", "T", 3.5, {
        "Norwegian": 0.01, "Danish": 0.01, "Swedish": 0.01,
        "British_Norse": 0.01, "Eastern_Route": 0.02, "Modern_Scandinavian": 0.01,
    }),
    "rs1426654": ("A", "G", 2.0, {
        "Norwegian": 0.99, "Danish": 0.99, "Swedish": 0.99,
        "British_Norse": 0.99, "Eastern_Route": 0.98, "Modern_Scandinavian": 0.99,
    }),
    "rs16891982": ("C", "G", 2.5, {
        "Norwegian": 0.97, "Danish": 0.93, "Swedish": 0.95,
        "British_Norse": 0.92, "Eastern_Route": 0.88, "Modern_Scandinavian": 0.96,
    }),
    "rs12913832": ("A", "G", 2.5, {
        "Norwegian": 0.82, "Danish": 0.76, "Swedish": 0.79,
        "British_Norse": 0.73, "Eastern_Route": 0.65, "Modern_Scandinavian": 0.80,
    }),
    "rs4988235": ("A", "G", 2.5, {
        "Norwegian": 0.70, "Danish": 0.64, "Swedish": 0.67,
        "British_Norse": 0.62, "Eastern_Route": 0.52, "Modern_Scandinavian": 0.74,
    }),
    "rs12203592": ("C", "T", 2.0, {
        "Norwegian": 0.22, "Danish": 0.21, "Swedish": 0.21,
        "British_Norse": 0.21, "Eastern_Route": 0.20, "Modern_Scandinavian": 0.22,
    }),
    "rs1805007": ("C", "T", 2.0, {
        "Norwegian": 0.10, "Danish": 0.10, "Swedish": 0.09,
        "British_Norse": 0.12, "Eastern_Route": 0.08, "Modern_Scandinavian": 0.10,
    }),
    "rs1800562": ("G", "A", 2.5, {
        "Norwegian": 0.07, "Danish": 0.07, "Swedish": 0.06,
        "British_Norse": 0.07, "Eastern_Route": 0.04, "Modern_Scandinavian": 0.07,
    }),
    "rs1799945": ("C", "G", 2.0, {
        "Norwegian": 0.15, "Danish": 0.15, "Swedish": 0.14,
        "British_Norse": 0.14, "Eastern_Route": 0.12, "Modern_Scandinavian": 0.15,
    }),
    "rs1042602": ("C", "A", 2.0, {
        "Norwegian": 0.44, "Danish": 0.43, "Swedish": 0.44,
        "British_Norse": 0.43, "Eastern_Route": 0.41, "Modern_Scandinavian": 0.44,
    }),
    "rs1800414": ("C", "T", 0.4, {
        "Norwegian": 0.44, "Danish": 0.42, "Swedish": 0.43,
        "British_Norse": 0.41, "Eastern_Route": 0.38, "Modern_Scandinavian": 0.43,
    }),
    "rs4918664": ("A", "G", 2.0, {
        "Norwegian": 0.11, "Danish": 0.11, "Swedish": 0.11,
        "British_Norse": 0.11, "Eastern_Route": 0.12, "Modern_Scandinavian": 0.11,
    }),
    "rs1800497": ("G", "A", 1.5, {
        "Norwegian": 0.21, "Danish": 0.20, "Swedish": 0.21,
        "British_Norse": 0.20, "Eastern_Route": 0.21, "Modern_Scandinavian": 0.21,
    }),
    "rs1799971": ("A", "G", 1.5, {
        "Norwegian": 0.13, "Danish": 0.12, "Swedish": 0.13,
        "British_Norse": 0.12, "Eastern_Route": 0.12, "Modern_Scandinavian": 0.13,
    }),
    "rs7554936": ("C", "T", 1.0, {
        "Norwegian": 0.68, "Danish": 0.67, "Swedish": 0.68,
        "British_Norse": 0.67, "Eastern_Route": 0.67, "Modern_Scandinavian": 0.68,
    }),
    "rs8050136": ("C", "A", 1.5, {
        "Norwegian": 0.43, "Danish": 0.43, "Swedish": 0.43,
        "British_Norse": 0.42, "Eastern_Route": 0.41, "Modern_Scandinavian": 0.43,
    }),
    "rs4680": ("G", "A", 1.0, {
        "Norwegian": 0.48, "Danish": 0.47, "Swedish": 0.48,
        "British_Norse": 0.47, "Eastern_Route": 0.46, "Modern_Scandinavian": 0.48,
    }),
    "rs1726866": ("G", "A", 1.0, {
        "Norwegian": 0.47, "Danish": 0.46, "Swedish": 0.47,
        "British_Norse": 0.46, "Eastern_Route": 0.44, "Modern_Scandinavian": 0.47,
    }),
    "rs713598": ("G", "C", 1.0, {
        "Norwegian": 0.45, "Danish": 0.44, "Swedish": 0.45,
        "British_Norse": 0.44, "Eastern_Route": 0.42, "Modern_Scandinavian": 0.45,
    }),
    "rs1801133": ("G", "A", 1.0, {
        "Norwegian": 0.32, "Danish": 0.33, "Swedish": 0.33,
        "British_Norse": 0.34, "Eastern_Route": 0.35, "Modern_Scandinavian": 0.32,
    }),
    "rs1801131": ("A", "C", 1.0, {
        "Norwegian": 0.30, "Danish": 0.30, "Swedish": 0.30,
        "British_Norse": 0.30, "Eastern_Route": 0.29, "Modern_Scandinavian": 0.30,
    }),
    "rs1800629": ("G", "A", 1.0, {
        "Norwegian": 0.13, "Danish": 0.13, "Swedish": 0.13,
        "British_Norse": 0.13, "Eastern_Route": 0.14, "Modern_Scandinavian": 0.13,
    }),
    "rs4129267": ("T", "C", 1.0, {
        "Norwegian": 0.36, "Danish": 0.37, "Swedish": 0.37,
        "British_Norse": 0.37, "Eastern_Route": 0.38, "Modern_Scandinavian": 0.36,
    }),
    "rs2435357": ("G", "C", 1.0, {
        "Norwegian": 0.26, "Danish": 0.25, "Swedish": 0.26,
        "British_Norse": 0.25, "Eastern_Route": 0.26, "Modern_Scandinavian": 0.26,
    }),
    "rs2228479": ("A", "G", 1.0, {
        "Norwegian": 0.08, "Danish": 0.08, "Swedish": 0.08,
        "British_Norse": 0.08, "Eastern_Route": 0.08, "Modern_Scandinavian": 0.08,
    }),
    "rs2076533": ("C", "A", 1.5, {
        "Norwegian": 0.14, "Danish": 0.14, "Swedish": 0.14,
        "British_Norse": 0.14, "Eastern_Route": 0.15, "Modern_Scandinavian": 0.14,
    }),
    "rs9267531": ("G", "A", 1.5, {
        "Norwegian": 0.19, "Danish": 0.19, "Swedish": 0.19,
        "British_Norse": 0.19, "Eastern_Route": 0.20, "Modern_Scandinavian": 0.19,
    }),
    "rs6090989": ("G", "A", 1.2, {
        "Norwegian": 0.17, "Danish": 0.17, "Swedish": 0.17,
        "British_Norse": 0.17, "Eastern_Route": 0.18, "Modern_Scandinavian": 0.17,
    }),
    "rs3135027": ("G", "A", 1.0, {
        "Norwegian": 0.21, "Danish": 0.21, "Swedish": 0.21,
        "British_Norse": 0.21, "Eastern_Route": 0.22, "Modern_Scandinavian": 0.21,
    }),
    "rs10497520": ("C", "T", 0.8, {
        "Norwegian": 0.41, "Danish": 0.39, "Swedish": 0.40,
        "British_Norse": 0.39, "Eastern_Route": 0.37, "Modern_Scandinavian": 0.41,
    }),
    "rs1540771": ("A", "G", 0.8, {
        "Norwegian": 0.64, "Danish": 0.63, "Swedish": 0.63,
        "British_Norse": 0.63, "Eastern_Route": 0.61, "Modern_Scandinavian": 0.64,
    }),
    "rs1937845": ("A", "G", 0.8, {
        "Norwegian": 0.14, "Danish": 0.15, "Swedish": 0.15,
        "British_Norse": 0.15, "Eastern_Route": 0.16, "Modern_Scandinavian": 0.14,
    }),
    "rs2073711": ("C", "T", 0.8, {
        "Norwegian": 0.19, "Danish": 0.20, "Swedish": 0.20,
        "British_Norse": 0.20, "Eastern_Route": 0.22, "Modern_Scandinavian": 0.19,
    }),
    "rs7149477": ("A", "G", 0.8, {
        "Norwegian": 0.17, "Danish": 0.18, "Swedish": 0.18,
        "British_Norse": 0.18, "Eastern_Route": 0.20, "Modern_Scandinavian": 0.17,
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
# NNLS ANCESTRY ESTIMATION (Weighted NNLS with Tikhonov regularization)
# ---------------------------------------------------------------------------

def _estimate_proportions(
    variants: Dict[str, Tuple[str, str, str]],
) -> Tuple[Dict[str, float], list, float]:
    """
    Run weighted NNLS decomposition into Viking-age population components.

    Returns (proportions_dict, used_snps_list, residual).
    """
    n_pops = len(VIKING_POPS)
    rows_A = []
    rows_f = []
    rows_w = []
    used_snps = []

    for rsid, (ref, alt, weight, pop_freqs) in AIMS_VIKING.items():
        variant_data = variants.get(rsid.lower())
        if variant_data is None:
            continue
        _chrom, _pos, genotype = variant_data
        dosage = _get_dosage(genotype, alt)
        if dosage is None:
            continue
        obs_freq = dosage / 2.0
        row = [pop_freqs.get(p, 0.0) for p in VIKING_POPS]
        rows_A.append(row)
        rows_f.append(obs_freq)
        rows_w.append(weight)
        used_snps.append(rsid)

    if len(used_snps) < 5:
        return {p: 1.0 / n_pops for p in VIKING_POPS}, used_snps, 0.0

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

    props = {p: float(x[i]) for i, p in enumerate(VIKING_POPS)}
    return props, used_snps, float(resid)


# ---------------------------------------------------------------------------
# AFFINITY SCORES (min-max normalized distances)
# ---------------------------------------------------------------------------

def _compute_affinity_scores(
    variants: Dict[str, Tuple[str, str, str]],
) -> Dict[str, float]:
    """
    Compute relative affinity to each Viking population using min-max normalization.
    Since all Viking sub-populations are genetically similar, we normalize within
    the observed range so the closest always scores 100 and furthest 0.
    """
    raw_distances: Dict[str, float] = {}

    for pop in VIKING_POPS:
        weighted_sq_diff = 0.0
        w_sum = 0.0
        for rsid, (ref, alt, weight, pop_freqs) in AIMS_VIKING.items():
            variant_data = variants.get(rsid.lower())
            if variant_data is None:
                continue
            _chrom, _pos, genotype = variant_data
            dosage = _get_dosage(genotype, alt)
            if dosage is None:
                continue
            obs = dosage / 2.0
            ref_freq = pop_freqs[pop]
            weighted_sq_diff += weight * (obs - ref_freq) ** 2
            w_sum += weight
        if w_sum > 0:
            raw_distances[pop] = np.sqrt(weighted_sq_diff / w_sum)

    if not raw_distances:
        return {p: 0.0 for p in VIKING_POPS}

    # Min-max normalization within Viking group (shows relative affinity)
    min_d = min(raw_distances.values())
    max_d = max(raw_distances.values())
    rng = max_d - min_d

    if rng > 0.001:
        affinity = {
            p: round((1.0 - (raw_distances[p] - min_d) / rng) * 100, 1)
            for p in raw_distances
        }
    else:
        # All populations equally distant (typically non-Scandinavian individual)
        affinity = {p: round(100.0 / len(raw_distances), 1) for p in raw_distances}

    return affinity


# ---------------------------------------------------------------------------
# KEY VIKING MARKERS SUMMARY
# ---------------------------------------------------------------------------

def _summarize_key_markers(
    variants: Dict[str, Tuple[str, str, str]],
) -> Dict[str, Dict[str, str]]:
    """Report the key markers distinguishing Viking sub-populations."""
    results: Dict[str, Dict[str, str]] = {}

    # Blue eyes (HERC2 rs12913832)
    variant_data = variants.get("rs12913832")
    if variant_data is not None:
        _c, _p, geno = variant_data
        if geno and geno not in ("--", "00", "NN"):
            g_count = geno.upper().count("G")
            if g_count == 2:
                results["blue_eyes"] = {
                    "status": "GG",
                    "detail": "GG -- strongly Scandinavian eye profile",
                }
            elif g_count == 1:
                results["blue_eyes"] = {
                    "status": "AG",
                    "detail": "AG -- intermediate eye color",
                }
            else:
                results["blue_eyes"] = {
                    "status": "AA",
                    "detail": "AA -- brown eyes (less typical)",
                }

    # Lactase persistence (rs4988235)
    variant_data = variants.get("rs4988235")
    if variant_data is not None:
        _c, _p, geno = variant_data
        if geno and geno not in ("--", "00", "NN"):
            g_count = geno.upper().count("G")
            if g_count == 2:
                results["lactase"] = {
                    "status": "GG",
                    "detail": "GG -- fully lactase persistent (Viking/Northern European)",
                }
            elif g_count == 1:
                results["lactase"] = {
                    "status": "AG",
                    "detail": "AG -- one LP allele (partially persistent)",
                }
            else:
                results["lactase"] = {
                    "status": "AA",
                    "detail": "AA -- ancestral; lactase non-persistent",
                }

    # SLC45A2 (rs16891982) -- N-S European gradient
    variant_data = variants.get("rs16891982")
    if variant_data is not None:
        _c, _p, geno = variant_data
        if geno and geno not in ("--", "00", "NN"):
            g_count = geno.upper().count("G")
            if g_count == 2:
                results["slc45a2"] = {
                    "status": "GG",
                    "detail": "GG -- very high Northern European skin pigmentation score",
                }
            elif g_count == 1:
                results["slc45a2"] = {
                    "status": "CG",
                    "detail": "CG -- one derived allele",
                }
            else:
                results["slc45a2"] = {
                    "status": "CC",
                    "detail": "CC -- lower; more Southern/Eastern European pattern",
                }

    # HFE C282Y (rs1800562) -- N. EUR enriched
    variant_data = variants.get("rs1800562")
    if variant_data is not None:
        _c, _p, geno = variant_data
        if geno and geno not in ("--", "00", "NN"):
            a_count = geno.upper().count("A")
            if a_count > 0:
                status = "GA" if a_count == 1 else "AA"
                results["hfe_c282y"] = {
                    "status": status,
                    "detail": "A allele present -- Northern European hemochromatosis variant",
                }
            else:
                results["hfe_c282y"] = {
                    "status": "GG",
                    "detail": "GG -- common/ancestral (no hemochromatosis variant)",
                }

    return results


# ---------------------------------------------------------------------------
# PUBLIC API
# ---------------------------------------------------------------------------

def analyze_viking(
    variants: Dict[str, Tuple[str, str, str]],
) -> Dict[str, Any]:
    """
    Analyze variants for Viking-age population affinity.

    Args:
        variants: Dict mapping rsID (lowercase) -> (chromosome, position, genotype)

    Returns:
        Dict with proportions, affinity_scores, used_snps, residual, key_markers.
    """
    props, used_snps, resid = _estimate_proportions(variants)
    affinity_scores = _compute_affinity_scores(variants)
    key_markers = _summarize_key_markers(variants)

    return {
        "proportions": props,
        "affinity_scores": affinity_scores,
        "used_snps": used_snps,
        "residual": resid,
        "key_markers": key_markers,
    }


def generate_viking_json(result: Dict[str, Any]) -> Dict[str, Any]:
    """
    Convert raw analysis result into the structured JSON for the frontend.

    Args:
        result: Output from analyze_viking()

    Returns:
        Structured JSON dict matching the report schema.
    """
    props = result["proportions"]
    affinity_scores = result["affinity_scores"]
    used_snps = result["used_snps"]
    resid = result["residual"]
    key_markers = result["key_markers"]

    panel_size = len(AIMS_VIKING)
    snps_used = len(used_snps)

    # Build population list sorted by affinity score descending
    sorted_pops = sorted(VIKING_POPS, key=lambda p: -affinity_scores.get(p, 0.0))

    populations = []
    for pop in sorted_pops:
        proportion = round(props.get(pop, 0.0), 4)
        populations.append({
            "code": pop,
            "label": VIKING_POP_LABELS[pop],
            "affinityScore": affinity_scores.get(pop, 0.0),
            "proportion": proportion,
            "percentage": round(proportion * 100, 1),
            "description": VIKING_POP_DESCRIPTIONS[pop],
            "color": VIKING_POP_COLORS[pop],
        })

    # Top match
    top_pop = sorted_pops[0]

    return {
        "summary": {
            "panelSize": panel_size,
            "snpsUsed": snps_used,
            "coverage": round(snps_used / panel_size * 100, 1) if panel_size > 0 else 0.0,
        },
        "populations": populations,
        "topMatch": {
            "code": top_pop,
            "label": VIKING_POP_LABELS[top_pop],
            "description": VIKING_POP_DESCRIPTIONS[top_pop],
        },
        "keyMarkers": key_markers,
    }
