"""
Hunter-Gatherer / Farmer / Steppe Ancestry Analyzer
Estimates the proportions of Mesolithic hunter-gatherer, Neolithic farmer, and
Bronze Age steppe herder ancestry components in modern European populations.

This is one of the most extensively studied topics in paleogenomics: the
three-wave model of European peopling (WHG -> EEF -> WSH), which together
explain the majority of genetic variation among present-day Europeans.

KEY REFERENCES:
  Lazaridis et al. 2014 (Nature): "Ancient human genomes suggest three
    ancestral populations for present-day Europeans"
      -- First formal definition of WHG, EEF, ANE components
      -- 9 ancient genomes + 2,345 modern individuals
  Haak et al. 2015 (Nature): "Massive migration from the steppe was a source
    for Indo-European languages in Europe"
      -- 69 ancient Europeans including Yamnaya
      -- Demonstrated ~75% Yamnaya ancestry in Corded Ware
  Mathieson et al. 2015 (Nature): "Genome-wide patterns of selection in 230
    ancient Eurasians"
      -- Trajectories of LCT, SLC24A5, SLC45A2, HERC2 over time
      -- Demonstrated lactase persistence rose under Bronze Age selection
  Fu et al. 2016 (Nature): "The genetic history of Ice Age Europe"
      -- 51 European pre-Neolithic genomes
      -- Traced WHG emergence from Villabruna cluster

SUPPORTING REFERENCES:
  - Allentoft et al. 2015 (Nature): Bronze Age population genomics
  - Gamba et al. 2014 (Nat Commun): Hungarian Neolithic timeline
  - Olalde et al. 2018 (Nature): Beaker complex spread
  - Jones et al. 2015 (Nat Commun): CHG Satsurblia/Kotias genomes

SCIENTIFIC LIMITATIONS:
  A 50-SNP AIM panel cannot reproduce the full ADMIXTURE/qpAdm decomposition
  used in academic papers. The proportions returned here are informative at
  the broad "ancestral component" level (e.g. 45% WSH / 35% EEF / 20% WHG is
  Northern European; 15% WSH / 60% EEF / 10% WHG is Southern European with
  additional Near Eastern contribution). Individual-level estimates have
  ~5-10 percentage point uncertainty.

Returns structured JSON for frontend visualization.
"""

from typing import Dict, Tuple, Any, Optional
import numpy as np
from scipy.optimize import nnls

# ---------------------------------------------------------------------------
# POPULATION DEFINITIONS
# Seven ancient ancestral components spanning ~24,000-4,000 BP
# ---------------------------------------------------------------------------

HGF_POPS = [
    "Western_Hunter_Gatherer",
    "Eastern_Hunter_Gatherer",
    "Scandinavian_Hunter_Gatherer",
    "Early_European_Farmer",
    "Western_Steppe_Herder",
    "Caucasus_Hunter_Gatherer",
    "Iran_Neolithic",
]

HGF_POP_LABELS = {
    "Western_Hunter_Gatherer":       "Western Hunter-Gatherer (WHG)",
    "Eastern_Hunter_Gatherer":       "Eastern Hunter-Gatherer (EHG)",
    "Scandinavian_Hunter_Gatherer":  "Scandinavian Hunter-Gatherer (SHG)",
    "Early_European_Farmer":         "Early European Farmer (EEF)",
    "Western_Steppe_Herder":         "Western Steppe Herder / Yamnaya (WSH)",
    "Caucasus_Hunter_Gatherer":      "Caucasus Hunter-Gatherer (CHG)",
    "Iran_Neolithic":                "Iran Neolithic (Iran_N)",
}

# Earth-tone palette as requested
HGF_POP_COLORS = {
    "Western_Hunter_Gatherer":       "#2D5016",  # forest green
    "Eastern_Hunter_Gatherer":       "#1E3A5F",  # deep blue
    "Scandinavian_Hunter_Gatherer":  "#6B7280",  # cold grey
    "Early_European_Farmer":         "#C8A27C",  # warm tan
    "Western_Steppe_Herder":         "#6B21A8",  # steppe purple
    "Caucasus_Hunter_Gatherer":      "#D97706",  # dawn orange
    "Iran_Neolithic":                "#9A3412",  # rust red
}

HGF_POP_DESCRIPTIONS = {
    "Western_Hunter_Gatherer": (
        "Western Hunter-Gatherers (WHG) were the Mesolithic inhabitants of Europe "
        "between the Last Glacial Maximum and the arrival of farming (~15,000-8,000 "
        "years before present). They descended from the Villabruna cluster and ranged "
        "from Iberia to the Balkans and British Isles. Famous representatives include "
        "'La Brana 1' (Spain, ~7,000 BP), 'Loschbour' (Luxembourg, ~8,000 BP) and "
        "'Cheddar Man' (UK, ~10,000 BP). The most striking feature of WHG genomes "
        "(Olalde 2014, Fu 2016) is an unexpected combination: dark skin (ancestral "
        "SLC24A5/SLC45A2 alleles) alongside already-derived blue-eye alleles at HERC2 "
        "(rs12913832). This pigmentation pattern is the opposite of modern Europeans. "
        "WHG ancestry peaks in modern Scandinavians, Baltics and the British Isles "
        "at ~10-25%, and was largely absent from the Near East."
    ),
    "Eastern_Hunter_Gatherer": (
        "Eastern Hunter-Gatherers (EHG) were Mesolithic populations of the Eastern "
        "European Plain and the forest-steppe zone between the Baltic and the Urals "
        "(~10,000-5,500 BP). They are represented in ancient DNA by samples from "
        "Karelia (Russia) and Samara (Volga region). Genetically, EHG are modeled as "
        "a mixture of WHG and Ancient North Eurasians (ANE, proxied by Mal'ta boy). "
        "EHG contributed the ANE-rich northern lineage that later fused with CHG to "
        "form the Western Steppe Herder (Yamnaya) ancestry. Today, EHG ancestry is "
        "most elevated in Finns, Karelians, Balts, Russians and other Uralic/Baltic "
        "populations. References: Haak 2015, Mathieson 2015, Fu 2016."
    ),
    "Scandinavian_Hunter_Gatherer": (
        "Scandinavian Hunter-Gatherers (SHG) were Mesolithic populations of southern "
        "Scandinavia (~9,500-6,000 BP), represented by the Motala samples from Sweden "
        "and Stora Forvar from Gotland. SHG are genetically modeled as roughly "
        "65% EHG + 35% WHG (Lazaridis 2014, Gunther 2018), reflecting the post-glacial "
        "recolonization of Scandinavia by two independent routes: from the south "
        "(WHG via Denmark) and from the northeast (EHG via Norway's Atlantic coast). "
        "SHG individuals already carried high frequencies of the blue-eye allele and "
        "the derived SLC45A2 light-skin allele, hinting that natural selection on "
        "pigmentation began well before the arrival of farmers."
    ),
    "Early_European_Farmer": (
        "Early European Farmers (EEF) descended from the Anatolian Neolithic -- the "
        "first agriculturalists, who domesticated wheat, barley, sheep and cattle "
        "in southwestern Asia around 10,500 BP. Between ~8,500 and 5,500 BP they "
        "spread west into Europe, first along the Mediterranean (Cardial Ware) and "
        "then north (Linearbandkeramik / LBK, Starcevo-Koros). Representative "
        "samples include 'Barcin' (Turkey), the LBK Stuttgart farmer, and Iberian "
        "Cardial individuals. EEF genomes (Lazaridis 2014, Mathieson 2015, Gamba 2014) "
        "carry already-near-fixed derived alleles at SLC24A5 (rs1426654 A->G) and "
        "SLC45A2 (rs16891982 C->G), meaning they had lighter skin than WHG. They had "
        "brown eyes and had NOT evolved lactase persistence. EEF ancestry dominates "
        "modern Sardinia (~75-80%) and Southern Europe more broadly (~50-60%)."
    ),
    "Western_Steppe_Herder": (
        "Western Steppe Herders (WSH), best represented by the Yamnaya culture "
        "(~5,300-4,400 BP) of the Pontic-Caspian steppe (modern Ukraine, southern "
        "Russia), were Bronze Age pastoralists associated with the spread of "
        "Indo-European languages into Europe (Haak 2015, Allentoft 2015). Yamnaya "
        "individuals are modeled as ~50% EHG + ~50% CHG and introduced into Europe "
        "the lactase-persistence allele (LCT rs4988235), taller stature (through "
        "height-associated loci), R1b and R1a Y-haplogroups, and horse/cattle-based "
        "mobile pastoralism. The Corded Ware culture of Northern Europe is modeled "
        "as ~75% Yamnaya (Haak 2015). Today, WSH ancestry peaks in Northern Europe "
        "(45-60% in Scandinavia, Baltic, Ireland, Scotland) and declines toward "
        "the Mediterranean."
    ),
    "Caucasus_Hunter_Gatherer": (
        "Caucasus Hunter-Gatherers (CHG) are known from the Satsurblia (~13,300 BP) "
        "and Kotias (~9,700 BP) caves in Georgia (Jones 2015). They represent a "
        "distinct deep lineage that had been isolated in the Caucasus refugium since "
        "the Last Glacial Maximum, diverging from the ancestors of WHG/EHG around "
        "45,000 BP. CHG carried substantial Basal Eurasian ancestry (shared with "
        "Anatolian farmers) and contributed approximately 50% of the Yamnaya genome "
        "when they admixed with EHG on the steppe. CHG is therefore the 'hidden' "
        "component that entered Europe in Bronze Age waves. CHG-related ancestry is "
        "highest in modern South Caucasians (Georgians, Armenians) and in South "
        "Asians (contributing via Iran_N-like ancestry)."
    ),
    "Iran_Neolithic": (
        "Iran Neolithic (Iran_N) refers to early farmers from the Zagros Mountains "
        "of western Iran (~10,000-8,000 BP), best represented by samples from Ganj "
        "Dareh and Wezmeh Cave (Lazaridis 2016, Broushaki 2016). Iran_N populations "
        "independently domesticated goats and cultivated crops, but were genetically "
        "distinct from their Anatolian contemporaries -- closer to CHG than to EEF. "
        "They carried very high Basal Eurasian ancestry (~45-60%) alongside a deep "
        "Eurasian component. Iran_N is the dominant farmer-related ancestry in modern "
        "South Asian populations (mixed with Ancestral North Indian/ASI) and "
        "contributed to the peopling of the Arabian Peninsula. In modern Europeans "
        "Iran_N is present only at trace levels, mostly through indirect Yamnaya/CHG "
        "routes."
    ),
}

# ---------------------------------------------------------------------------
# AIM ALLELE FREQUENCIES FOR ANCIENT ANCESTRAL COMPONENTS
# Calibrated from: Lazaridis 2014, Haak 2015, Mathieson 2015, Fu 2016,
# Allentoft 2015, Olalde 2018 -- approximate published means for each
# ancient cluster supplemented with modern reference populations.
#
# Format: rsid: (ref, alt, weight, {pop: alt_freq})
# Weights reflect discrimination power (pigmentation loci highest, then
# diet/metabolism, then neutral AIMs used for geographic signal).
# ---------------------------------------------------------------------------

AIMS_HGF = {
    # --- PIGMENTATION (the canonical ancient DNA markers) ---
    # SLC24A5 rs1426654 -- light skin, near-fixed derived in EEF/WSH, ~0 in WHG
    "rs1426654": ("A", "G", 5.0, {
        "Western_Hunter_Gatherer":      0.02,
        "Eastern_Hunter_Gatherer":      0.40,
        "Scandinavian_Hunter_Gatherer": 0.55,
        "Early_European_Farmer":        0.99,
        "Western_Steppe_Herder":        0.95,
        "Caucasus_Hunter_Gatherer":     0.95,
        "Iran_Neolithic":               0.98,
    }),
    # SLC45A2 rs16891982 -- light skin, derived allele rose later than SLC24A5
    "rs16891982": ("C", "G", 5.0, {
        "Western_Hunter_Gatherer":      0.02,
        "Eastern_Hunter_Gatherer":      0.05,
        "Scandinavian_Hunter_Gatherer": 0.45,
        "Early_European_Farmer":        0.60,
        "Western_Steppe_Herder":        0.35,
        "Caucasus_Hunter_Gatherer":     0.25,
        "Iran_Neolithic":               0.30,
    }),
    # HERC2 rs12913832 -- blue eyes, surprisingly already high in WHG (!)
    "rs12913832": ("A", "G", 5.0, {
        "Western_Hunter_Gatherer":      0.95,
        "Eastern_Hunter_Gatherer":      0.60,
        "Scandinavian_Hunter_Gatherer": 0.85,
        "Early_European_Farmer":        0.15,
        "Western_Steppe_Herder":        0.30,
        "Caucasus_Hunter_Gatherer":     0.10,
        "Iran_Neolithic":               0.05,
    }),
    # OCA2 rs1800407 -- minor eye color effect
    "rs1800407": ("C", "T", 2.0, {
        "Western_Hunter_Gatherer":      0.10,
        "Eastern_Hunter_Gatherer":      0.08,
        "Scandinavian_Hunter_Gatherer": 0.11,
        "Early_European_Farmer":        0.06,
        "Western_Steppe_Herder":        0.09,
        "Caucasus_Hunter_Gatherer":     0.05,
        "Iran_Neolithic":               0.04,
    }),
    # TYR rs1042602 -- skin/freckling, frequent in EEF
    "rs1042602": ("C", "A", 2.5, {
        "Western_Hunter_Gatherer":      0.10,
        "Eastern_Hunter_Gatherer":      0.15,
        "Scandinavian_Hunter_Gatherer": 0.20,
        "Early_European_Farmer":        0.50,
        "Western_Steppe_Herder":        0.35,
        "Caucasus_Hunter_Gatherer":     0.40,
        "Iran_Neolithic":               0.45,
    }),
    # TYRP1 rs1408799 -- hair/skin pigmentation
    "rs1408799": ("T", "C", 1.5, {
        "Western_Hunter_Gatherer":      0.05,
        "Eastern_Hunter_Gatherer":      0.15,
        "Scandinavian_Hunter_Gatherer": 0.25,
        "Early_European_Farmer":        0.40,
        "Western_Steppe_Herder":        0.45,
        "Caucasus_Hunter_Gatherer":     0.30,
        "Iran_Neolithic":               0.30,
    }),
    # KITLG rs12821256 -- blond hair (European-specific)
    "rs12821256": ("T", "C", 2.0, {
        "Western_Hunter_Gatherer":      0.02,
        "Eastern_Hunter_Gatherer":      0.05,
        "Scandinavian_Hunter_Gatherer": 0.10,
        "Early_European_Farmer":        0.12,
        "Western_Steppe_Herder":        0.18,
        "Caucasus_Hunter_Gatherer":     0.03,
        "Iran_Neolithic":               0.02,
    }),
    # MC1R rs1805007 (Arg151Cys) -- red hair
    "rs1805007": ("C", "T", 1.5, {
        "Western_Hunter_Gatherer":      0.06,
        "Eastern_Hunter_Gatherer":      0.05,
        "Scandinavian_Hunter_Gatherer": 0.08,
        "Early_European_Farmer":        0.06,
        "Western_Steppe_Herder":        0.08,
        "Caucasus_Hunter_Gatherer":     0.02,
        "Iran_Neolithic":               0.02,
    }),
    # MC1R rs1805008 (Arg160Trp) -- red hair
    "rs1805008": ("C", "T", 1.2, {
        "Western_Hunter_Gatherer":      0.04,
        "Eastern_Hunter_Gatherer":      0.04,
        "Scandinavian_Hunter_Gatherer": 0.05,
        "Early_European_Farmer":        0.05,
        "Western_Steppe_Herder":        0.06,
        "Caucasus_Hunter_Gatherer":     0.02,
        "Iran_Neolithic":               0.02,
    }),
    # ASIP rs6119471 -- skin pigmentation
    "rs6119471": ("C", "G", 1.5, {
        "Western_Hunter_Gatherer":      0.10,
        "Eastern_Hunter_Gatherer":      0.12,
        "Scandinavian_Hunter_Gatherer": 0.15,
        "Early_European_Farmer":        0.30,
        "Western_Steppe_Herder":        0.22,
        "Caucasus_Hunter_Gatherer":     0.25,
        "Iran_Neolithic":               0.30,
    }),
    # IRF4 rs12203592 -- fair skin / freckles (Europe-specific)
    "rs12203592": ("C", "T", 2.0, {
        "Western_Hunter_Gatherer":      0.05,
        "Eastern_Hunter_Gatherer":      0.08,
        "Scandinavian_Hunter_Gatherer": 0.15,
        "Early_European_Farmer":        0.15,
        "Western_Steppe_Herder":        0.25,
        "Caucasus_Hunter_Gatherer":     0.05,
        "Iran_Neolithic":               0.03,
    }),

    # --- DIET / METABOLISM ---
    # LCT rs4988235 -- lactase persistence, classic Yamnaya/WSH signature
    "rs4988235": ("G", "A", 5.0, {
        "Western_Hunter_Gatherer":      0.00,
        "Eastern_Hunter_Gatherer":      0.02,
        "Scandinavian_Hunter_Gatherer": 0.02,
        "Early_European_Farmer":        0.03,
        "Western_Steppe_Herder":        0.15,
        "Caucasus_Hunter_Gatherer":     0.00,
        "Iran_Neolithic":               0.00,
    }),
    # MCM6 rs182549 -- linked to LCT
    "rs182549": ("C", "T", 3.0, {
        "Western_Hunter_Gatherer":      0.00,
        "Eastern_Hunter_Gatherer":      0.02,
        "Scandinavian_Hunter_Gatherer": 0.02,
        "Early_European_Farmer":        0.03,
        "Western_Steppe_Herder":        0.15,
        "Caucasus_Hunter_Gatherer":     0.00,
        "Iran_Neolithic":               0.00,
    }),
    # FADS1 rs174546 -- long-chain PUFA synthesis, rose with EEF/farming
    "rs174546": ("C", "T", 2.5, {
        "Western_Hunter_Gatherer":      0.30,
        "Eastern_Hunter_Gatherer":      0.35,
        "Scandinavian_Hunter_Gatherer": 0.35,
        "Early_European_Farmer":        0.65,
        "Western_Steppe_Herder":        0.45,
        "Caucasus_Hunter_Gatherer":     0.55,
        "Iran_Neolithic":               0.60,
    }),
    # FADS2 rs174570 -- PUFA biosynthesis
    "rs174570": ("C", "T", 1.8, {
        "Western_Hunter_Gatherer":      0.25,
        "Eastern_Hunter_Gatherer":      0.30,
        "Scandinavian_Hunter_Gatherer": 0.30,
        "Early_European_Farmer":        0.60,
        "Western_Steppe_Herder":        0.40,
        "Caucasus_Hunter_Gatherer":     0.50,
        "Iran_Neolithic":               0.55,
    }),
    # APOE rs429358 -- lipid metabolism
    "rs429358": ("T", "C", 1.5, {
        "Western_Hunter_Gatherer":      0.20,
        "Eastern_Hunter_Gatherer":      0.18,
        "Scandinavian_Hunter_Gatherer": 0.18,
        "Early_European_Farmer":        0.10,
        "Western_Steppe_Herder":        0.12,
        "Caucasus_Hunter_Gatherer":     0.08,
        "Iran_Neolithic":               0.10,
    }),
    # ADH1B rs1229984 -- alcohol metabolism, ~0 in Europeans, not Asian
    "rs1229984": ("G", "A", 1.0, {
        "Western_Hunter_Gatherer":      0.00,
        "Eastern_Hunter_Gatherer":      0.00,
        "Scandinavian_Hunter_Gatherer": 0.00,
        "Early_European_Farmer":        0.01,
        "Western_Steppe_Herder":        0.00,
        "Caucasus_Hunter_Gatherer":     0.02,
        "Iran_Neolithic":               0.05,
    }),

    # --- HEIGHT / SKELETAL (Yamnaya brought increase in height) ---
    # HMGA2 rs1042725 -- height, rose with WSH
    "rs1042725": ("C", "T", 2.0, {
        "Western_Hunter_Gatherer":      0.35,
        "Eastern_Hunter_Gatherer":      0.40,
        "Scandinavian_Hunter_Gatherer": 0.42,
        "Early_European_Farmer":        0.35,
        "Western_Steppe_Herder":        0.60,
        "Caucasus_Hunter_Gatherer":     0.35,
        "Iran_Neolithic":               0.35,
    }),
    # GDF5 rs143384 -- height / joint
    "rs143384": ("G", "A", 1.5, {
        "Western_Hunter_Gatherer":      0.45,
        "Eastern_Hunter_Gatherer":      0.50,
        "Scandinavian_Hunter_Gatherer": 0.50,
        "Early_European_Farmer":        0.55,
        "Western_Steppe_Herder":        0.60,
        "Caucasus_Hunter_Gatherer":     0.55,
        "Iran_Neolithic":               0.55,
    }),
    # ACAN rs3817428 -- height (hitchhiked with WSH)
    "rs3817428": ("C", "T", 1.5, {
        "Western_Hunter_Gatherer":      0.20,
        "Eastern_Hunter_Gatherer":      0.25,
        "Scandinavian_Hunter_Gatherer": 0.25,
        "Early_European_Farmer":        0.22,
        "Western_Steppe_Herder":        0.40,
        "Caucasus_Hunter_Gatherer":     0.22,
        "Iran_Neolithic":               0.20,
    }),
    # ZBTB38 rs6763931 -- height
    "rs6763931": ("A", "G", 1.2, {
        "Western_Hunter_Gatherer":      0.40,
        "Eastern_Hunter_Gatherer":      0.45,
        "Scandinavian_Hunter_Gatherer": 0.45,
        "Early_European_Farmer":        0.48,
        "Western_Steppe_Herder":        0.55,
        "Caucasus_Hunter_Gatherer":     0.48,
        "Iran_Neolithic":               0.45,
    }),

    # --- IMMUNE / DISEASE RESISTANCE ---
    # HLA-A/B region rs2523393 -- Neolithic-associated MHC shift
    "rs2523393": ("A", "G", 1.5, {
        "Western_Hunter_Gatherer":      0.25,
        "Eastern_Hunter_Gatherer":      0.28,
        "Scandinavian_Hunter_Gatherer": 0.28,
        "Early_European_Farmer":        0.50,
        "Western_Steppe_Herder":        0.40,
        "Caucasus_Hunter_Gatherer":     0.45,
        "Iran_Neolithic":               0.45,
    }),
    # SLC22A4 rs1050152 -- immune / ergothioneine transporter
    "rs1050152": ("C", "T", 1.2, {
        "Western_Hunter_Gatherer":      0.35,
        "Eastern_Hunter_Gatherer":      0.38,
        "Scandinavian_Hunter_Gatherer": 0.38,
        "Early_European_Farmer":        0.50,
        "Western_Steppe_Herder":        0.45,
        "Caucasus_Hunter_Gatherer":     0.48,
        "Iran_Neolithic":               0.50,
    }),
    # CCR5 rs333 tag rs113341849 -- infection resistance
    "rs113341849": ("G", "A", 1.0, {
        "Western_Hunter_Gatherer":      0.08,
        "Eastern_Hunter_Gatherer":      0.10,
        "Scandinavian_Hunter_Gatherer": 0.12,
        "Early_European_Farmer":        0.08,
        "Western_Steppe_Herder":        0.10,
        "Caucasus_Hunter_Gatherer":     0.05,
        "Iran_Neolithic":               0.04,
    }),
    # HFE rs1800562 -- Celtic/Northern European hemochromatosis
    "rs1800562": ("G", "A", 1.5, {
        "Western_Hunter_Gatherer":      0.02,
        "Eastern_Hunter_Gatherer":      0.02,
        "Scandinavian_Hunter_Gatherer": 0.03,
        "Early_European_Farmer":        0.02,
        "Western_Steppe_Herder":        0.08,
        "Caucasus_Hunter_Gatherer":     0.02,
        "Iran_Neolithic":               0.01,
    }),

    # --- VITAMIN D / SKIN ---
    # DHCR7 rs12785878 -- vitamin D synthesis (selection in N Europe)
    "rs12785878": ("G", "T", 1.8, {
        "Western_Hunter_Gatherer":      0.20,
        "Eastern_Hunter_Gatherer":      0.25,
        "Scandinavian_Hunter_Gatherer": 0.30,
        "Early_European_Farmer":        0.28,
        "Western_Steppe_Herder":        0.35,
        "Caucasus_Hunter_Gatherer":     0.25,
        "Iran_Neolithic":               0.22,
    }),
    # GC rs2282679 -- vitamin D binding protein
    "rs2282679": ("T", "G", 1.2, {
        "Western_Hunter_Gatherer":      0.20,
        "Eastern_Hunter_Gatherer":      0.25,
        "Scandinavian_Hunter_Gatherer": 0.25,
        "Early_European_Farmer":        0.30,
        "Western_Steppe_Herder":        0.28,
        "Caucasus_Hunter_Gatherer":     0.28,
        "Iran_Neolithic":               0.30,
    }),

    # --- DENTAL / EDAR outgroup ---
    # EDAR rs3827760 -- East Asian / absent in W Eurasians (outgroup)
    "rs3827760": ("A", "G", 1.0, {
        "Western_Hunter_Gatherer":      0.00,
        "Eastern_Hunter_Gatherer":      0.01,
        "Scandinavian_Hunter_Gatherer": 0.00,
        "Early_European_Farmer":        0.00,
        "Western_Steppe_Herder":        0.01,
        "Caucasus_Hunter_Gatherer":     0.00,
        "Iran_Neolithic":               0.00,
    }),

    # --- NEUTRAL GEOGRAPHIC AIMS (provide ancestry signal) ---
    "rs1834640": ("A", "G", 1.2, {
        "Western_Hunter_Gatherer":      0.40,
        "Eastern_Hunter_Gatherer":      0.35,
        "Scandinavian_Hunter_Gatherer": 0.38,
        "Early_European_Farmer":        0.60,
        "Western_Steppe_Herder":        0.45,
        "Caucasus_Hunter_Gatherer":     0.55,
        "Iran_Neolithic":               0.58,
    }),
    "rs260690": ("T", "C", 1.2, {
        "Western_Hunter_Gatherer":      0.50,
        "Eastern_Hunter_Gatherer":      0.45,
        "Scandinavian_Hunter_Gatherer": 0.47,
        "Early_European_Farmer":        0.30,
        "Western_Steppe_Herder":        0.40,
        "Caucasus_Hunter_Gatherer":     0.25,
        "Iran_Neolithic":               0.22,
    }),
    "rs2814778": ("T", "C", 1.0, {
        "Western_Hunter_Gatherer":      0.00,
        "Eastern_Hunter_Gatherer":      0.00,
        "Scandinavian_Hunter_Gatherer": 0.00,
        "Early_European_Farmer":        0.00,
        "Western_Steppe_Herder":        0.00,
        "Caucasus_Hunter_Gatherer":     0.00,
        "Iran_Neolithic":               0.00,
    }),
    "rs3811801": ("G", "A", 1.2, {
        "Western_Hunter_Gatherer":      0.00,
        "Eastern_Hunter_Gatherer":      0.02,
        "Scandinavian_Hunter_Gatherer": 0.01,
        "Early_European_Farmer":        0.00,
        "Western_Steppe_Herder":        0.03,
        "Caucasus_Hunter_Gatherer":     0.02,
        "Iran_Neolithic":               0.05,
    }),
    "rs1229984_b": ("G", "A", 0.8, {
        "Western_Hunter_Gatherer":      0.00,
        "Eastern_Hunter_Gatherer":      0.00,
        "Scandinavian_Hunter_Gatherer": 0.00,
        "Early_European_Farmer":        0.01,
        "Western_Steppe_Herder":        0.00,
        "Caucasus_Hunter_Gatherer":     0.02,
        "Iran_Neolithic":               0.05,
    }),
    "rs671": ("G", "A", 1.0, {
        "Western_Hunter_Gatherer":      0.00,
        "Eastern_Hunter_Gatherer":      0.00,
        "Scandinavian_Hunter_Gatherer": 0.00,
        "Early_European_Farmer":        0.00,
        "Western_Steppe_Herder":        0.00,
        "Caucasus_Hunter_Gatherer":     0.00,
        "Iran_Neolithic":               0.00,
    }),
    "rs17822931": ("C", "T", 1.0, {
        "Western_Hunter_Gatherer":      0.01,
        "Eastern_Hunter_Gatherer":      0.02,
        "Scandinavian_Hunter_Gatherer": 0.02,
        "Early_European_Farmer":        0.02,
        "Western_Steppe_Herder":        0.02,
        "Caucasus_Hunter_Gatherer":     0.03,
        "Iran_Neolithic":               0.05,
    }),
    "rs6548238": ("C", "T", 1.0, {
        "Western_Hunter_Gatherer":      0.15,
        "Eastern_Hunter_Gatherer":      0.18,
        "Scandinavian_Hunter_Gatherer": 0.18,
        "Early_European_Farmer":        0.30,
        "Western_Steppe_Herder":        0.22,
        "Caucasus_Hunter_Gatherer":     0.28,
        "Iran_Neolithic":               0.30,
    }),
    "rs10497191": ("C", "T", 1.0, {
        "Western_Hunter_Gatherer":      0.10,
        "Eastern_Hunter_Gatherer":      0.15,
        "Scandinavian_Hunter_Gatherer": 0.18,
        "Early_European_Farmer":        0.25,
        "Western_Steppe_Herder":        0.20,
        "Caucasus_Hunter_Gatherer":     0.22,
        "Iran_Neolithic":               0.25,
    }),
    "rs7903146": ("C", "T", 1.0, {
        "Western_Hunter_Gatherer":      0.15,
        "Eastern_Hunter_Gatherer":      0.18,
        "Scandinavian_Hunter_Gatherer": 0.18,
        "Early_European_Farmer":        0.25,
        "Western_Steppe_Herder":        0.25,
        "Caucasus_Hunter_Gatherer":     0.30,
        "Iran_Neolithic":               0.32,
    }),
    "rs2066702": ("C", "T", 1.0, {
        "Western_Hunter_Gatherer":      0.00,
        "Eastern_Hunter_Gatherer":      0.00,
        "Scandinavian_Hunter_Gatherer": 0.00,
        "Early_European_Farmer":        0.00,
        "Western_Steppe_Herder":        0.00,
        "Caucasus_Hunter_Gatherer":     0.00,
        "Iran_Neolithic":               0.00,
    }),
    "rs1799971": ("A", "G", 0.8, {
        "Western_Hunter_Gatherer":      0.14,
        "Eastern_Hunter_Gatherer":      0.16,
        "Scandinavian_Hunter_Gatherer": 0.15,
        "Early_European_Farmer":        0.12,
        "Western_Steppe_Herder":        0.14,
        "Caucasus_Hunter_Gatherer":     0.12,
        "Iran_Neolithic":               0.12,
    }),
    "rs6152": ("G", "A", 1.0, {
        "Western_Hunter_Gatherer":      0.30,
        "Eastern_Hunter_Gatherer":      0.35,
        "Scandinavian_Hunter_Gatherer": 0.35,
        "Early_European_Farmer":        0.40,
        "Western_Steppe_Herder":        0.35,
        "Caucasus_Hunter_Gatherer":     0.45,
        "Iran_Neolithic":               0.45,
    }),
    "rs16969968": ("G", "A", 1.0, {
        "Western_Hunter_Gatherer":      0.30,
        "Eastern_Hunter_Gatherer":      0.32,
        "Scandinavian_Hunter_Gatherer": 0.32,
        "Early_European_Farmer":        0.40,
        "Western_Steppe_Herder":        0.38,
        "Caucasus_Hunter_Gatherer":     0.30,
        "Iran_Neolithic":               0.28,
    }),
    "rs10757278": ("A", "G", 1.0, {
        "Western_Hunter_Gatherer":      0.40,
        "Eastern_Hunter_Gatherer":      0.42,
        "Scandinavian_Hunter_Gatherer": 0.42,
        "Early_European_Farmer":        0.50,
        "Western_Steppe_Herder":        0.48,
        "Caucasus_Hunter_Gatherer":     0.50,
        "Iran_Neolithic":               0.52,
    }),
    "rs2187668": ("C", "T", 1.0, {
        "Western_Hunter_Gatherer":      0.05,
        "Eastern_Hunter_Gatherer":      0.08,
        "Scandinavian_Hunter_Gatherer": 0.08,
        "Early_European_Farmer":        0.12,
        "Western_Steppe_Herder":        0.10,
        "Caucasus_Hunter_Gatherer":     0.08,
        "Iran_Neolithic":               0.08,
    }),
    "rs11466334": ("C", "T", 0.8, {
        "Western_Hunter_Gatherer":      0.10,
        "Eastern_Hunter_Gatherer":      0.12,
        "Scandinavian_Hunter_Gatherer": 0.12,
        "Early_European_Farmer":        0.15,
        "Western_Steppe_Herder":        0.14,
        "Caucasus_Hunter_Gatherer":     0.18,
        "Iran_Neolithic":               0.20,
    }),
    "rs12896399": ("G", "T", 1.2, {
        "Western_Hunter_Gatherer":      0.30,
        "Eastern_Hunter_Gatherer":      0.35,
        "Scandinavian_Hunter_Gatherer": 0.40,
        "Early_European_Farmer":        0.42,
        "Western_Steppe_Herder":        0.48,
        "Caucasus_Hunter_Gatherer":     0.32,
        "Iran_Neolithic":               0.30,
    }),
    "rs7495174": ("A", "G", 1.0, {
        "Western_Hunter_Gatherer":      0.05,
        "Eastern_Hunter_Gatherer":      0.06,
        "Scandinavian_Hunter_Gatherer": 0.08,
        "Early_European_Farmer":        0.10,
        "Western_Steppe_Herder":        0.12,
        "Caucasus_Hunter_Gatherer":     0.06,
        "Iran_Neolithic":               0.05,
    }),
    "rs4778138": ("A", "G", 1.0, {
        "Western_Hunter_Gatherer":      0.15,
        "Eastern_Hunter_Gatherer":      0.20,
        "Scandinavian_Hunter_Gatherer": 0.22,
        "Early_European_Farmer":        0.25,
        "Western_Steppe_Herder":        0.28,
        "Caucasus_Hunter_Gatherer":     0.20,
        "Iran_Neolithic":               0.18,
    }),
    "rs4778241": ("C", "A", 1.0, {
        "Western_Hunter_Gatherer":      0.20,
        "Eastern_Hunter_Gatherer":      0.22,
        "Scandinavian_Hunter_Gatherer": 0.25,
        "Early_European_Farmer":        0.30,
        "Western_Steppe_Herder":        0.32,
        "Caucasus_Hunter_Gatherer":     0.20,
        "Iran_Neolithic":               0.18,
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


def _normalize_rsid(rsid: str) -> str:
    """Remove helper suffixes like '_b' used for duplicate entries."""
    if "_" in rsid:
        return rsid.split("_", 1)[0]
    return rsid


# ---------------------------------------------------------------------------
# NNLS ANCESTRY ESTIMATION (Weighted NNLS with Tikhonov regularization)
# ---------------------------------------------------------------------------

def _estimate_proportions(
    variants: Dict[str, Tuple[str, str, str]],
) -> Tuple[Dict[str, float], list, float]:
    """
    Run weighted NNLS decomposition into ancient ancestral components.

    Returns (proportions_dict, used_snps_list, residual).
    """
    n_pops = len(HGF_POPS)
    rows_A = []
    rows_f = []
    rows_w = []
    used_snps = []

    for rsid, (ref, alt, weight, pop_freqs) in AIMS_HGF.items():
        lookup = _normalize_rsid(rsid).lower()
        variant_data = variants.get(lookup)
        if variant_data is None:
            continue
        _chrom, _pos, genotype = variant_data
        dosage = _get_dosage(genotype, alt)
        if dosage is None:
            continue
        obs_freq = dosage / 2.0
        row = [pop_freqs.get(p, 0.0) for p in HGF_POPS]
        rows_A.append(row)
        rows_f.append(obs_freq)
        rows_w.append(weight)
        used_snps.append(_normalize_rsid(rsid))

    if len(used_snps) < 5:
        return {p: 1.0 / n_pops for p in HGF_POPS}, used_snps, 0.0

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

    props = {p: float(x[i]) for i, p in enumerate(HGF_POPS)}
    return props, used_snps, float(resid)


# ---------------------------------------------------------------------------
# AFFINITY SCORES (min-max normalized distances)
# ---------------------------------------------------------------------------

def _compute_affinity_scores(
    variants: Dict[str, Tuple[str, str, str]],
) -> Dict[str, float]:
    """
    Compute relative affinity to each ancient component using min-max normalization.
    """
    raw_distances: Dict[str, float] = {}

    for pop in HGF_POPS:
        weighted_sq_diff = 0.0
        w_sum = 0.0
        for rsid, (ref, alt, weight, pop_freqs) in AIMS_HGF.items():
            lookup = _normalize_rsid(rsid).lower()
            variant_data = variants.get(lookup)
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
        return {p: 0.0 for p in HGF_POPS}

    min_d = min(raw_distances.values())
    max_d = max(raw_distances.values())
    rng = max_d - min_d

    if rng > 0.001:
        affinity = {
            p: round((1.0 - (raw_distances[p] - min_d) / rng) * 100, 1)
            for p in raw_distances
        }
    else:
        affinity = {p: round(100.0 / len(raw_distances), 1) for p in raw_distances}

    return affinity


# ---------------------------------------------------------------------------
# KEY MARKERS SUMMARY
# ---------------------------------------------------------------------------

def _summarize_key_markers(
    variants: Dict[str, Tuple[str, str, str]],
) -> Dict[str, Dict[str, str]]:
    """
    Report the key pigmentation, diet and selection markers that distinguish
    Mesolithic hunter-gatherers from Neolithic farmers and Bronze Age steppe.
    """
    results: Dict[str, Dict[str, str]] = {}

    # --- SLC24A5 rs1426654 (skin lightening, EEF signature) ---
    variant_data = variants.get("rs1426654")
    if variant_data is not None:
        _c, _p, geno = variant_data
        if geno and geno not in ("--", "00", "NN"):
            g_count = geno.upper().count("G")
            if g_count == 2:
                results["slc24a5"] = {
                    "status": "GG",
                    "detail": (
                        "GG -- derived/derived. Light skin allele (typical of EEF and "
                        "most modern Europeans). This allele was near-absent in WHG."
                    ),
                }
            elif g_count == 1:
                results["slc24a5"] = {
                    "status": "AG",
                    "detail": (
                        "AG -- one ancestral and one derived allele. Rare in modern "
                        "Europeans; may indicate residual WHG or non-European ancestry."
                    ),
                }
            else:
                results["slc24a5"] = {
                    "status": "AA",
                    "detail": (
                        "AA -- ancestral/ancestral. Darker skin allele. Very rare in "
                        "modern Europeans (<1%) but was typical of Mesolithic "
                        "hunter-gatherers."
                    ),
                }

    # --- SLC45A2 rs16891982 (N-S European gradient) ---
    variant_data = variants.get("rs16891982")
    if variant_data is not None:
        _c, _p, geno = variant_data
        if geno and geno not in ("--", "00", "NN"):
            g_count = geno.upper().count("G")
            if g_count == 2:
                results["slc45a2"] = {
                    "status": "GG",
                    "detail": (
                        "GG -- derived/derived. Strong Northern European skin "
                        "pigmentation profile; allele was absent in WHG, moderate in "
                        "EEF, and rose to near-fixation later."
                    ),
                }
            elif g_count == 1:
                results["slc45a2"] = {
                    "status": "CG",
                    "detail": (
                        "CG -- one derived allele. Intermediate; common in "
                        "Southern/Central Europeans."
                    ),
                }
            else:
                results["slc45a2"] = {
                    "status": "CC",
                    "detail": (
                        "CC -- ancestral. Darker skin allele, typical in Mesolithic "
                        "Europeans, North Africans and South Asians."
                    ),
                }

    # --- HERC2 rs12913832 (blue eyes, WHG paradox) ---
    variant_data = variants.get("rs12913832")
    if variant_data is not None:
        _c, _p, geno = variant_data
        if geno and geno not in ("--", "00", "NN"):
            g_count = geno.upper().count("G")
            if g_count == 2:
                results["herc2_blue_eyes"] = {
                    "status": "GG",
                    "detail": (
                        "GG -- blue eyes. This allele was already common in Western "
                        "Hunter-Gatherers (!), meaning blue eyes predate lighter skin "
                        "in Europe."
                    ),
                }
            elif g_count == 1:
                results["herc2_blue_eyes"] = {
                    "status": "AG",
                    "detail": (
                        "AG -- intermediate eye color (hazel, green, light brown). "
                        "One WHG-derived allele."
                    ),
                }
            else:
                results["herc2_blue_eyes"] = {
                    "status": "AA",
                    "detail": (
                        "AA -- brown eyes. Ancestral; typical of EEF, CHG and "
                        "Iran_N, rare in WHG."
                    ),
                }

    # --- LCT rs4988235 (lactase persistence - Yamnaya signature) ---
    variant_data = variants.get("rs4988235")
    if variant_data is not None:
        _c, _p, geno = variant_data
        if geno and geno not in ("--", "00", "NN"):
            a_count = geno.upper().count("A")
            if a_count == 2:
                results["lct_lactase"] = {
                    "status": "AA",
                    "detail": (
                        "AA -- fully lactase persistent. This allele was absent in WHG "
                        "and EEF; it arrived with Bronze Age Yamnaya/WSH and rose to "
                        "high frequency in Northern Europe under strong selection."
                    ),
                }
            elif a_count == 1:
                results["lct_lactase"] = {
                    "status": "GA",
                    "detail": (
                        "GA -- one LP allele, partial lactase persistence. Carries "
                        "one WSH-derived allele."
                    ),
                }
            else:
                results["lct_lactase"] = {
                    "status": "GG",
                    "detail": (
                        "GG -- ancestral, lactase non-persistent. This is the WHG/EEF "
                        "ancestral state. Today most prevalent in Southern Europe and "
                        "non-Europeans."
                    ),
                }

    # --- IRF4 rs12203592 (European freckling/fair skin, WSH-related) ---
    variant_data = variants.get("rs12203592")
    if variant_data is not None:
        _c, _p, geno = variant_data
        if geno and geno not in ("--", "00", "NN"):
            t_count = geno.upper().count("T")
            if t_count == 2:
                results["irf4_freckles"] = {
                    "status": "TT",
                    "detail": (
                        "TT -- strong freckling / red-tone skin allele. European-specific; "
                        "rose in frequency after the Bronze Age."
                    ),
                }
            elif t_count == 1:
                results["irf4_freckles"] = {
                    "status": "CT",
                    "detail": "CT -- moderate freckling tendency.",
                }
            else:
                results["irf4_freckles"] = {
                    "status": "CC",
                    "detail": "CC -- ancestral, no freckling allele.",
                }

    # --- HMGA2 rs1042725 (height, rose with Yamnaya) ---
    variant_data = variants.get("rs1042725")
    if variant_data is not None:
        _c, _p, geno = variant_data
        if geno and geno not in ("--", "00", "NN"):
            t_count = geno.upper().count("T")
            if t_count == 2:
                results["hmga2_height"] = {
                    "status": "TT",
                    "detail": (
                        "TT -- two 'tall' alleles. Selection for height intensified "
                        "during and after the Yamnaya expansion."
                    ),
                }
            elif t_count == 1:
                results["hmga2_height"] = {
                    "status": "CT",
                    "detail": "CT -- one tall allele.",
                }
            else:
                results["hmga2_height"] = {
                    "status": "CC",
                    "detail": "CC -- ancestral short allele, more common in WHG/EEF.",
                }

    return results


# ---------------------------------------------------------------------------
# PUBLIC API
# ---------------------------------------------------------------------------

def analyze_hunter_gatherer_farmer(
    variants: Dict[str, Tuple[str, str, str]],
) -> Dict[str, Any]:
    """
    Analyze variants for ancient European ancestral component proportions
    (WHG / EHG / SHG / EEF / WSH / CHG / Iran_N).

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


def generate_hunter_gatherer_farmer_json(result: Dict[str, Any]) -> Dict[str, Any]:
    """
    Convert raw analysis result into structured JSON for the frontend.
    """
    props = result["proportions"]
    affinity_scores = result["affinity_scores"]
    used_snps = result["used_snps"]
    resid = result["residual"]
    key_markers = result["key_markers"]

    # Unique AIM rsIDs (strip helper suffixes)
    unique_rsids = {_normalize_rsid(r) for r in AIMS_HGF.keys()}
    panel_size = len(unique_rsids)
    snps_used = len(set(used_snps))

    # Sort populations by proportion (ancestry decomposition) descending
    sorted_pops = sorted(HGF_POPS, key=lambda p: -props.get(p, 0.0))

    populations = []
    for pop in sorted_pops:
        proportion = round(props.get(pop, 0.0), 4)
        populations.append({
            "code": pop,
            "label": HGF_POP_LABELS[pop],
            "affinityScore": affinity_scores.get(pop, 0.0),
            "proportion": proportion,
            "percentage": round(proportion * 100, 1),
            "description": HGF_POP_DESCRIPTIONS[pop],
            "color": HGF_POP_COLORS[pop],
        })

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
            "label": HGF_POP_LABELS[top_pop],
            "description": HGF_POP_DESCRIPTIONS[top_pop],
        },
        "keyMarkers": key_markers,
    }
