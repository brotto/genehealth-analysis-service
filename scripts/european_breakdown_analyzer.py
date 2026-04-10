"""
European Breakdown Analyzer
Estimates fine-grained European ancestry (Northern / Southern / Eastern /
Baltic / Mediterranean / Basque / Ashkenazi / Caucasian) using Ancestry-
Informative Markers (AIMs) calibrated to published European genomics data.

KEY REFERENCES:
  Novembre et al. 2008 (Nature): "Genes mirror geography within Europe"
    - 3,192 individuals across Europe; PCA of 500k SNPs
    - Demonstrated that genetic structure in Europe follows geographic
      distance with striking precision (the famous PC1-PC2 map of Europe)

  Lazaridis et al. 2014 (Nature): "Ancient human genomes suggest three
    ancestral populations for present-day Europeans"
    - Identified WHG, EEF, ANE/Steppe ancestral components

  Lazaridis et al. 2016 (Nature): "Genomic insights into the origin of
    farming in the ancient Near East"

  Leslie et al. 2015 (Nature): "The fine-scale genetic structure of the
    British population" - The 'People of the British Isles' study

  Lao et al. 2008 (Current Biology): "Correlation between genetic and
    geographic structure in Europe" - 299k SNPs, 2,457 individuals

SUPPORTING REFERENCES:
  - Haak et al. 2015 (Nature): Bronze Age Steppe migration into Europe
  - Olalde et al. 2018 (Nature): Beaker complex and British genome replacement
  - Behar et al. 2010 (Nature): Jewish people's diaspora genetics
  - Bray et al. 2010 (AJHG): Ashkenazi founder effect
  - Bergstrom et al. 2020 (Science): HGDP high-coverage panel
  - Eiberg et al. 2008 (HG): HERC2/OCA2 blue-eye ancestry founder (~10,000 yrs)

SCIENTIFIC LIMITATIONS:
  European populations are genetically similar compared to continental
  contrasts (most Europeans are ~99.9% identical at the SNP level). Our
  ~60-SNP panel can distinguish broad N-S-E-W regional trends but will
  struggle to resolve e.g. Swedish vs. Danish or Northern vs. Southern
  Italian at the individual level. Results are RELATIVE affinity scores,
  not direct descent. Basque and Ashkenazi Jewish clusters are the most
  genetically distinctive European populations and are the easiest to
  distinguish with modest SNP panels. 'Caucasian West Asian' (Armenian,
  Georgian) is included as a border reference population, not implying
  these are European per se.

Returns structured JSON for frontend visualization.
"""

from typing import Dict, Tuple, Any, Optional
import numpy as np
from scipy.optimize import nnls

# ---------------------------------------------------------------------------
# POPULATION DEFINITIONS
# ---------------------------------------------------------------------------

EUROPEAN_POPS = [
    "Northwest_European",
    "Scandinavian",
    "Iberian",
    "Basque",
    "Italian_Mediterranean",
    "Balkan_Southeast",
    "Eastern_European_Slavic",
    "Baltic_Finnic",
    "Ashkenazi_Jewish",
    "Caucasian_West_Asian",
]

EUROPEAN_POP_LABELS = {
    "Northwest_European":       "Northwest European (British, Irish, Dutch, German, French)",
    "Scandinavian":             "Scandinavian (Norwegian, Swedish, Danish, Icelandic)",
    "Iberian":                  "Iberian (Spanish, Portuguese)",
    "Basque":                   "Basque (Euskaldunak)",
    "Italian_Mediterranean":    "Italian / Mediterranean (Southern Italian, Sicilian, Sardinian)",
    "Balkan_Southeast":         "Balkan / Southeast European (Greek, Albanian, South Slav)",
    "Eastern_European_Slavic":  "Eastern European Slavic (Polish, Ukrainian, Russian, Belarusian)",
    "Baltic_Finnic":            "Baltic / Finnic (Finnish, Estonian, Latvian)",
    "Ashkenazi_Jewish":         "Ashkenazi Jewish",
    "Caucasian_West_Asian":     "Caucasian / West Asian (Armenian, Georgian)",
}

# Geographic-themed palette: N Europe cool blues/whites, S Europe warm
# reds/oranges, E Europe purples, Mediterranean golds
EUROPEAN_POP_COLORS = {
    "Northwest_European":       "#3B82F6",  # blue
    "Scandinavian":             "#60A5FA",  # light blue / frost
    "Iberian":                  "#F97316",  # orange
    "Basque":                   "#DC2626",  # crimson red
    "Italian_Mediterranean":    "#EAB308",  # Mediterranean gold
    "Balkan_Southeast":         "#EA580C",  # burnt orange
    "Eastern_European_Slavic":  "#8B5CF6",  # violet
    "Baltic_Finnic":            "#06B6D4",  # cyan
    "Ashkenazi_Jewish":         "#A855F7",  # purple
    "Caucasian_West_Asian":     "#EC4899",  # pink-rose
}

EUROPEAN_POP_DESCRIPTIONS = {
    "Northwest_European": (
        "Northwest Europeans (British, Irish, Dutch, German, Northern French) "
        "carry the highest proportion of Bronze Age Steppe (Yamnaya-related) "
        "ancestry in Europe, alongside substantial Early European Farmer (EEF) "
        "and Western Hunter-Gatherer (WHG) components. Leslie 2015 (the People "
        "of the British Isles study) identified 17 fine-scale clusters within "
        "Britain alone, reflecting post-Roman Anglo-Saxon and Danish migrations. "
        "They have very high frequencies of lactase persistence (rs4988235 T "
        "allele), the blue-eye HERC2 rs12913832 G allele, and European skin-"
        "pigmentation alleles (SLC45A2, SLC24A5). Olalde 2018 showed the "
        "Beaker complex (~4,500 years ago) replaced ~90% of British Neolithic "
        "farmer DNA, making modern British primarily descended from Bell Beaker "
        "peoples of continental origin."
    ),
    "Scandinavian": (
        "Scandinavians (Norwegian, Swedish, Danish, Icelandic, Faroese) represent "
        "a Northern European cluster genetically adjacent to Northwest Europeans "
        "but with even higher Western Hunter-Gatherer (WHG) and Scandinavian "
        "Hunter-Gatherer (SHG) ancestry. Margaryan 2020 (Nature) showed that "
        "Viking-age Scandinavians were already genetically heterogeneous, with "
        "Swedish Vikings carrying more Baltic/Slavic ancestry and Norwegian "
        "Vikings carrying more WHG. Icelanders are the direct founder population "
        "of Norwegian Vikings (~874 CE) with minimal subsequent admixture, "
        "making them among the most genetically homogeneous European groups. "
        "Scandinavians have the world's highest frequencies of blue eyes and "
        "red hair (MC1R variants)."
    ),
    "Iberian": (
        "Iberians (Spanish, Portuguese, excluding Basques) reflect one of Europe's "
        "most complex genetic histories. Olalde 2019 (Science) showed that Bronze "
        "Age Steppe migration replaced ~40% of Iberian Y-chromosomes ~4,500 years "
        "ago but left female mitochondrial lineages largely intact. Modern Iberians "
        "carry additional North African (~1-10%, higher in Southern Spain/Portugal) "
        "ancestry from the Umayyad and Almoravid periods (711-1492 CE), and Roman-"
        "era Italian admixture. They have lower lactase persistence than Northern "
        "Europeans (~30-40% LP) and higher frequencies of the Mediterranean HLA "
        "A*02-B*50-DRB1*07 haplotype. Iberian populations also preserve rare "
        "pre-Indo-European Y-chromosome lineages (especially Portugal's I2a1a)."
    ),
    "Basque": (
        "The Basque (Euskaldunak) are a genetic and linguistic isolate in "
        "southwestern France and northern Spain. Their language, Euskara, has "
        "no known linguistic relatives and may represent a pre-Indo-European "
        "language of Western Europe. Gunther 2015 and Olalde 2019 showed that "
        "Basques carry the highest proportion of Early European Farmer (EEF) "
        "ancestry in modern Europe combined with WHG, and the LOWEST proportion "
        "of Bronze Age Steppe ancestry (~10-15% vs 25-50% in other Europeans). "
        "This makes Basques the closest living genetic proxies to pre-Indo-"
        "European Neolithic farmers of Western Europe. They have distinctive "
        "blood-group frequencies (highest Rh-negative ~25% in Europe), unique "
        "HLA haplotypes, and a 40% frequency of the M173 Y-chromosome lineage."
    ),
    "Italian_Mediterranean": (
        "Southern Italians and Sicilians cluster closer to Greeks and Levantines "
        "than to Northern Europeans, reflecting the Magna Graecia colonization "
        "(~8th c. BCE), Byzantine rule, and deeper Neolithic farmer ancestry. "
        "Sardinians are a special case: they are the closest living population "
        "to Otzi the Iceman (~5,300 years ago) and early Neolithic European "
        "farmers generally, due to Sardinia's geographic isolation and minimal "
        "Bronze Age Steppe admixture. Sicilians carry 1-5% North African ancestry "
        "from Aghlabid / Kalbid periods (~827-1091 CE). Mediterranean populations "
        "have significantly higher FADS1/2 omega-3 metabolism alleles (selection "
        "for olive-oil based diets) and lower lactase persistence than Northern "
        "Europeans."
    ),
    "Balkan_Southeast": (
        "Balkan and Southeast European populations (Greek, Albanian, Bulgarian, "
        "Serbian, Macedonian, Romanian) represent a genetic bridge between "
        "Northern Europe, Mediterranean, and the Near East. They carry elevated "
        "Neolithic Anatolian farmer ancestry, Bronze Age Steppe migration, and "
        "some Byzantine-era Near Eastern admixture. Mathieson 2018 (Nature) "
        "identified Balkan Neolithic individuals as carrying one of the earliest "
        "forms of Neolithic farmer ancestry that later spread west. Modern Greeks "
        "are remarkably similar to Bronze Age Mycenaeans (Lazaridis 2017, Nature). "
        "South Slavs (Serbian, Bulgarian) received Slavic migration admixture "
        "during the 6-7th c. CE. The region shows one of Europe's highest HLA "
        "diversities."
    ),
    "Eastern_European_Slavic": (
        "Eastern European Slavs (Polish, Ukrainian, Russian, Belarusian) share a "
        "common genetic signature distinguished by higher Bronze Age Steppe (Yamnaya "
        "/ Corded Ware) ancestry than most other Europeans and by a detectable "
        "component of ancient East European Hunter-Gatherer (EHG) ancestry. "
        "Modern Russians, especially from northern Russia, also carry ~5-15% "
        "Finno-Ugric / Uralic ancestry from pre-Slavic inhabitants. Slavic "
        "expansion from the 5-7th c. CE carried this signature across Eastern "
        "Europe. Polish populations cluster closest to Lithuanians and Germans; "
        "Ukrainians bridge Russian and Balkan; Belarusians cluster closest to "
        "Lithuanians and Poles. Kushniarevich 2015 and Behar 2013 provide "
        "detailed Slavic genetic structure."
    ),
    "Baltic_Finnic": (
        "Baltic Finnic speakers (Finnish, Estonian, Karelian, Veps) and Baltic "
        "speakers (Latvian, Lithuanian) form a distinctive Northeast European "
        "cluster. Finns are the most genetically isolated and homogeneous "
        "European population, with a pronounced founder effect and bottleneck "
        "(~2,000 years ago) -- the 'Finnish disease heritage' includes dozens "
        "of rare Mendelian disorders enriched by drift. Lamnidis 2018 (Nat Comm) "
        "and Saag 2019 showed that Finns carry distinctive East Asian (Siberian) "
        "ancestry (~5-10%) absent in most other Europeans, reflecting ancient "
        "Uralic-speaking migration from the Volga-Ural region. Finns also have "
        "the highest lactase persistence frequency in the world (~80%), reflecting "
        "strong selection under dairy pastoralism in high latitudes. Estonians "
        "and Latvians are intermediate between Finns and Lithuanians."
    ),
    "Ashkenazi_Jewish": (
        "Ashkenazi Jews (Central/Eastern European Jewish diaspora) form a "
        "distinctive genetic cluster with ~45-55% Middle Eastern (Levantine) "
        "ancestry, ~45-55% European (mostly Southern/Mediterranean + Northern) "
        "ancestry, and a well-documented founder effect/bottleneck (~350 founders "
        "~800 years ago). Behar 2010 and Atzmon 2010 showed that Ashkenazi and "
        "Sephardi Jews cluster together and are closest to Levantine Arabs. The "
        "population has characteristic founder mutations including BRCA1 "
        "185delAG (1% carrier frequency), Tay-Sachs HEXA, Niemann-Pick, Gaucher, "
        "Bloom syndrome, Canavan, and Familial Mediterranean Fever. Carmi 2014 "
        "sequenced 128 Ashkenazi genomes to characterize the bottleneck. "
        "Ashkenazi Jews have distinctive HLA class II haplotypes and a "
        "characteristic Y-chromosome signature (J1/J2 and R1a1a 'Levite' cluster)."
    ),
    "Caucasian_West_Asian": (
        "Caucasian / West Asian populations (Armenian, Georgian, Abkhaz, Chechen, "
        "Azerbaijani) are included as a border reference. They are technically "
        "not European but cluster adjacent to Europeans and share substantial "
        "ancestry with Neolithic Anatolian farmers, who are the main source of "
        "European farmer (EEF) ancestry. Armenians show remarkable genetic "
        "continuity with Bronze Age populations of the region. Georgians have "
        "been proposed as the closest modern proxies to the Caucasus "
        "Hunter-Gatherers (CHG) identified by Jones 2015 (Nat Comm) -- a key "
        "component contributing to European Bronze Age Steppe ancestry via "
        "Yamnaya admixture. People with apparent 'Mediterranean' ancestry often "
        "score here due to the shared Neolithic Anatolian component, plus "
        "West Asian admixture from the Byzantine and early Islamic periods."
    ),
}

# ---------------------------------------------------------------------------
# AIM ALLELE FREQUENCIES FOR EUROPEAN POPULATIONS
# Frequencies compiled from gnomAD (non-Finnish European, Finnish), 1000
# Genomes (CEU, GBR, FIN, TSI, IBS), HGDP (Basque, Sardinian, Russian,
# Adygei, Orcadian), Novembre 2008 AIM panel, Lazaridis 2014/2016 supplementary,
# Lao 2008 pan-European AIMs, and Behar 2010 Jewish populations.
#
# For each rsid: (ref, alt, weight, {pop: alt_freq})
# ---------------------------------------------------------------------------

AIMS_EUROPEAN = {
    # --- LACTASE PERSISTENCE (N-S gradient) ---
    # rs4988235: LCT -13910 C>T. Strong N-S gradient in Europe.
    # FIN ~80%, GBR ~70%, CEU ~72%, TSI ~30%, IBS ~35%, Sardinian ~6%.
    "rs4988235": ("G", "A", 4.0, {
        "Northwest_European": 0.72, "Scandinavian": 0.76,
        "Iberian": 0.42, "Basque": 0.48,
        "Italian_Mediterranean": 0.30, "Balkan_Southeast": 0.35,
        "Eastern_European_Slavic": 0.55, "Baltic_Finnic": 0.80,
        "Ashkenazi_Jewish": 0.35, "Caucasian_West_Asian": 0.25,
    }),
    # --- EYE COLOR ---
    # rs12913832: HERC2 blue-eye (A=ancestral dark, G=derived blue).
    # Scandinavia ~85%, UK ~78%, Iberia ~45%, Italy ~25%, Ashkenazi ~40%.
    "rs12913832": ("A", "G", 3.5, {
        "Northwest_European": 0.78, "Scandinavian": 0.85,
        "Iberian": 0.45, "Basque": 0.52,
        "Italian_Mediterranean": 0.28, "Balkan_Southeast": 0.35,
        "Eastern_European_Slavic": 0.62, "Baltic_Finnic": 0.82,
        "Ashkenazi_Jewish": 0.40, "Caucasian_West_Asian": 0.20,
    }),
    # --- SKIN PIGMENTATION ---
    # rs1426654: SLC24A5 A111T. ~99% in Europeans overall -- minimal variation.
    "rs1426654": ("A", "G", 1.2, {
        "Northwest_European": 0.995, "Scandinavian": 0.998,
        "Iberian": 0.995, "Basque": 0.995,
        "Italian_Mediterranean": 0.99, "Balkan_Southeast": 0.98,
        "Eastern_European_Slavic": 0.99, "Baltic_Finnic": 0.998,
        "Ashkenazi_Jewish": 0.97, "Caucasian_West_Asian": 0.90,
    }),
    # rs16891982: SLC45A2 F374L (derived G = light skin). N-S gradient.
    "rs16891982": ("C", "G", 3.5, {
        "Northwest_European": 0.94, "Scandinavian": 0.96,
        "Iberian": 0.85, "Basque": 0.88,
        "Italian_Mediterranean": 0.70, "Balkan_Southeast": 0.78,
        "Eastern_European_Slavic": 0.92, "Baltic_Finnic": 0.95,
        "Ashkenazi_Jewish": 0.82, "Caucasian_West_Asian": 0.55,
    }),
    # rs1042602: TYR S192Y.
    "rs1042602": ("C", "A", 2.5, {
        "Northwest_European": 0.42, "Scandinavian": 0.45,
        "Iberian": 0.35, "Basque": 0.40,
        "Italian_Mediterranean": 0.30, "Balkan_Southeast": 0.33,
        "Eastern_European_Slavic": 0.40, "Baltic_Finnic": 0.48,
        "Ashkenazi_Jewish": 0.32, "Caucasian_West_Asian": 0.22,
    }),
    # rs1800407: OCA2 R305W (light iris color).
    "rs1800407": ("C", "T", 1.8, {
        "Northwest_European": 0.10, "Scandinavian": 0.12,
        "Iberian": 0.05, "Basque": 0.07,
        "Italian_Mediterranean": 0.04, "Balkan_Southeast": 0.05,
        "Eastern_European_Slavic": 0.09, "Baltic_Finnic": 0.13,
        "Ashkenazi_Jewish": 0.06, "Caucasian_West_Asian": 0.03,
    }),
    # rs12203592: IRF4 freckles and sun sensitivity.
    "rs12203592": ("C", "T", 2.5, {
        "Northwest_European": 0.20, "Scandinavian": 0.24,
        "Iberian": 0.08, "Basque": 0.10,
        "Italian_Mediterranean": 0.05, "Balkan_Southeast": 0.07,
        "Eastern_European_Slavic": 0.15, "Baltic_Finnic": 0.22,
        "Ashkenazi_Jewish": 0.08, "Caucasian_West_Asian": 0.04,
    }),
    # rs1800414: OCA2 H615R (East Asian-enriched, very low in Europeans).
    "rs1800414": ("T", "C", 0.8, {
        "Northwest_European": 0.00, "Scandinavian": 0.00,
        "Iberian": 0.00, "Basque": 0.00,
        "Italian_Mediterranean": 0.00, "Balkan_Southeast": 0.00,
        "Eastern_European_Slavic": 0.00, "Baltic_Finnic": 0.01,
        "Ashkenazi_Jewish": 0.00, "Caucasian_West_Asian": 0.00,
    }),
    # rs885479: MC1R red hair variant.
    "rs885479": ("G", "A", 1.5, {
        "Northwest_European": 0.15, "Scandinavian": 0.14,
        "Iberian": 0.10, "Basque": 0.12,
        "Italian_Mediterranean": 0.08, "Balkan_Southeast": 0.09,
        "Eastern_European_Slavic": 0.12, "Baltic_Finnic": 0.13,
        "Ashkenazi_Jewish": 0.10, "Caucasian_West_Asian": 0.07,
    }),
    # rs1805007: MC1R R151C (red hair/pale skin, Celtic-enriched).
    "rs1805007": ("C", "T", 2.5, {
        "Northwest_European": 0.12, "Scandinavian": 0.10,
        "Iberian": 0.05, "Basque": 0.07,
        "Italian_Mediterranean": 0.03, "Balkan_Southeast": 0.04,
        "Eastern_European_Slavic": 0.08, "Baltic_Finnic": 0.08,
        "Ashkenazi_Jewish": 0.05, "Caucasian_West_Asian": 0.02,
    }),
    # rs1805008: MC1R R160W.
    "rs1805008": ("C", "T", 2.0, {
        "Northwest_European": 0.08, "Scandinavian": 0.07,
        "Iberian": 0.03, "Basque": 0.04,
        "Italian_Mediterranean": 0.02, "Balkan_Southeast": 0.03,
        "Eastern_European_Slavic": 0.06, "Baltic_Finnic": 0.06,
        "Ashkenazi_Jewish": 0.04, "Caucasian_West_Asian": 0.01,
    }),
    # rs1805009: MC1R D294H.
    "rs1805009": ("G", "C", 1.5, {
        "Northwest_European": 0.02, "Scandinavian": 0.02,
        "Iberian": 0.01, "Basque": 0.02,
        "Italian_Mediterranean": 0.01, "Balkan_Southeast": 0.01,
        "Eastern_European_Slavic": 0.02, "Baltic_Finnic": 0.02,
        "Ashkenazi_Jewish": 0.01, "Caucasian_West_Asian": 0.005,
    }),
    # --- VITAMIN D SYNTHESIS (adapted in N Europeans) ---
    # rs12785878: NADSYN1/DHCR7 vit D. Higher derived allele in N Europe.
    "rs12785878": ("T", "G", 2.5, {
        "Northwest_European": 0.72, "Scandinavian": 0.76,
        "Iberian": 0.55, "Basque": 0.58,
        "Italian_Mediterranean": 0.50, "Balkan_Southeast": 0.52,
        "Eastern_European_Slavic": 0.68, "Baltic_Finnic": 0.75,
        "Ashkenazi_Jewish": 0.55, "Caucasian_West_Asian": 0.42,
    }),
    # rs10741657: CYP2R1 vit D.
    "rs10741657": ("G", "A", 1.8, {
        "Northwest_European": 0.38, "Scandinavian": 0.40,
        "Iberian": 0.30, "Basque": 0.32,
        "Italian_Mediterranean": 0.28, "Balkan_Southeast": 0.30,
        "Eastern_European_Slavic": 0.35, "Baltic_Finnic": 0.42,
        "Ashkenazi_Jewish": 0.30, "Caucasian_West_Asian": 0.22,
    }),
    # --- OMEGA-3 METABOLISM (adapted in Mediterraneans) ---
    # rs174537: FADS1 cluster. Higher derived allele in populations with
    # plant-based diets (Mediterranean, South Asian).
    "rs174537": ("G", "T", 2.5, {
        "Northwest_European": 0.35, "Scandinavian": 0.32,
        "Iberian": 0.48, "Basque": 0.45,
        "Italian_Mediterranean": 0.55, "Balkan_Southeast": 0.50,
        "Eastern_European_Slavic": 0.40, "Baltic_Finnic": 0.30,
        "Ashkenazi_Jewish": 0.48, "Caucasian_West_Asian": 0.52,
    }),
    # rs1535: FADS2.
    "rs1535": ("A", "G", 1.8, {
        "Northwest_European": 0.32, "Scandinavian": 0.30,
        "Iberian": 0.42, "Basque": 0.40,
        "Italian_Mediterranean": 0.50, "Balkan_Southeast": 0.45,
        "Eastern_European_Slavic": 0.36, "Baltic_Finnic": 0.28,
        "Ashkenazi_Jewish": 0.45, "Caucasian_West_Asian": 0.48,
    }),
    # --- HFE IRON (Northern European enriched) ---
    # rs1800562: HFE C282Y. Northern European hemochromatosis variant.
    "rs1800562": ("G", "A", 3.0, {
        "Northwest_European": 0.08, "Scandinavian": 0.09,
        "Iberian": 0.03, "Basque": 0.07,
        "Italian_Mediterranean": 0.02, "Balkan_Southeast": 0.02,
        "Eastern_European_Slavic": 0.04, "Baltic_Finnic": 0.06,
        "Ashkenazi_Jewish": 0.02, "Caucasian_West_Asian": 0.01,
    }),
    # rs1799945: HFE H63D.
    "rs1799945": ("C", "G", 1.5, {
        "Northwest_European": 0.15, "Scandinavian": 0.14,
        "Iberian": 0.21, "Basque": 0.30,  # Basque particularly enriched
        "Italian_Mediterranean": 0.12, "Balkan_Southeast": 0.14,
        "Eastern_European_Slavic": 0.14, "Baltic_Finnic": 0.12,
        "Ashkenazi_Jewish": 0.14, "Caucasian_West_Asian": 0.08,
    }),
    # --- HAIR TEXTURE / EDAR ---
    # rs3827760: EDAR V370A. East Asian; detectable only in Finnish/Russian.
    "rs3827760": ("A", "G", 2.0, {
        "Northwest_European": 0.00, "Scandinavian": 0.00,
        "Iberian": 0.00, "Basque": 0.00,
        "Italian_Mediterranean": 0.00, "Balkan_Southeast": 0.00,
        "Eastern_European_Slavic": 0.01, "Baltic_Finnic": 0.04,
        "Ashkenazi_Jewish": 0.00, "Caucasian_West_Asian": 0.00,
    }),
    # --- EARWAX ---
    # rs17822931: ABCC11 G/A wet/dry earwax. Very low in Europeans, but
    # detectable in Finnic populations due to East Asian admixture.
    "rs17822931": ("C", "T", 2.0, {
        "Northwest_European": 0.01, "Scandinavian": 0.01,
        "Iberian": 0.01, "Basque": 0.01,
        "Italian_Mediterranean": 0.01, "Balkan_Southeast": 0.01,
        "Eastern_European_Slavic": 0.02, "Baltic_Finnic": 0.05,
        "Ashkenazi_Jewish": 0.02, "Caucasian_West_Asian": 0.05,
    }),
    # --- TLR SELECTION ---
    # rs5743611: TLR1. Subject to selective sweep in Europeans.
    "rs5743611": ("C", "G", 1.5, {
        "Northwest_European": 0.28, "Scandinavian": 0.30,
        "Iberian": 0.22, "Basque": 0.25,
        "Italian_Mediterranean": 0.20, "Balkan_Southeast": 0.22,
        "Eastern_European_Slavic": 0.25, "Baltic_Finnic": 0.32,
        "Ashkenazi_Jewish": 0.20, "Caucasian_West_Asian": 0.18,
    }),
    # rs5743810: TLR6 P249S.
    "rs5743810": ("A", "G", 1.5, {
        "Northwest_European": 0.40, "Scandinavian": 0.42,
        "Iberian": 0.35, "Basque": 0.38,
        "Italian_Mediterranean": 0.30, "Balkan_Southeast": 0.33,
        "Eastern_European_Slavic": 0.38, "Baltic_Finnic": 0.42,
        "Ashkenazi_Jewish": 0.32, "Caucasian_West_Asian": 0.25,
    }),
    # rs5744168: TLR5 stop-gain (F616L). European selective sweep.
    "rs5744168": ("C", "T", 1.5, {
        "Northwest_European": 0.10, "Scandinavian": 0.11,
        "Iberian": 0.08, "Basque": 0.09,
        "Italian_Mediterranean": 0.07, "Balkan_Southeast": 0.08,
        "Eastern_European_Slavic": 0.10, "Baltic_Finnic": 0.12,
        "Ashkenazi_Jewish": 0.08, "Caucasian_West_Asian": 0.06,
    }),
    # --- APOE (alzheimer risk, cosmopolitan pattern) ---
    "rs429358": ("T", "C", 1.5, {
        "Northwest_European": 0.15, "Scandinavian": 0.19,
        "Iberian": 0.11, "Basque": 0.13,
        "Italian_Mediterranean": 0.08, "Balkan_Southeast": 0.10,
        "Eastern_European_Slavic": 0.14, "Baltic_Finnic": 0.22,
        "Ashkenazi_Jewish": 0.09, "Caucasian_West_Asian": 0.10,
    }),
    "rs7412": ("C", "T", 1.0, {
        "Northwest_European": 0.08, "Scandinavian": 0.09,
        "Iberian": 0.07, "Basque": 0.07,
        "Italian_Mediterranean": 0.06, "Balkan_Southeast": 0.07,
        "Eastern_European_Slavic": 0.07, "Baltic_Finnic": 0.06,
        "Ashkenazi_Jewish": 0.06, "Caucasian_West_Asian": 0.07,
    }),
    # --- ADH METABOLISM ---
    # rs1229984: ADH1B His48Arg. Mostly East Asian; notable in Ashkenazi
    # and West Asian populations (Near Eastern introgression).
    "rs1229984": ("G", "A", 2.5, {
        "Northwest_European": 0.04, "Scandinavian": 0.04,
        "Iberian": 0.06, "Basque": 0.05,
        "Italian_Mediterranean": 0.08, "Balkan_Southeast": 0.06,
        "Eastern_European_Slavic": 0.05, "Baltic_Finnic": 0.05,
        "Ashkenazi_Jewish": 0.18, "Caucasian_West_Asian": 0.15,
    }),
    # rs671: ALDH2*2. East Asian only, near-zero in Europeans.
    "rs671": ("G", "A", 1.0, {
        "Northwest_European": 0.00, "Scandinavian": 0.00,
        "Iberian": 0.00, "Basque": 0.00,
        "Italian_Mediterranean": 0.00, "Balkan_Southeast": 0.00,
        "Eastern_European_Slavic": 0.00, "Baltic_Finnic": 0.00,
        "Ashkenazi_Jewish": 0.00, "Caucasian_West_Asian": 0.01,
    }),
    # --- CCR5 delta32 (Northern/Baltic enriched) ---
    "rs333": ("I", "D", 3.0, {
        "Northwest_European": 0.10, "Scandinavian": 0.12,
        "Iberian": 0.06, "Basque": 0.07,
        "Italian_Mediterranean": 0.05, "Balkan_Southeast": 0.05,
        "Eastern_European_Slavic": 0.12, "Baltic_Finnic": 0.17,
        "Ashkenazi_Jewish": 0.10, "Caucasian_West_Asian": 0.03,
    }),
    # --- FOLATE METABOLISM ---
    # rs1801133: MTHFR C677T. Mediterranean enriched.
    "rs1801133": ("G", "A", 2.0, {
        "Northwest_European": 0.32, "Scandinavian": 0.30,
        "Iberian": 0.42, "Basque": 0.38,
        "Italian_Mediterranean": 0.48, "Balkan_Southeast": 0.42,
        "Eastern_European_Slavic": 0.35, "Baltic_Finnic": 0.32,
        "Ashkenazi_Jewish": 0.42, "Caucasian_West_Asian": 0.45,
    }),
    # rs1801131: MTHFR A1298C.
    "rs1801131": ("A", "C", 1.2, {
        "Northwest_European": 0.30, "Scandinavian": 0.29,
        "Iberian": 0.32, "Basque": 0.30,
        "Italian_Mediterranean": 0.28, "Balkan_Southeast": 0.30,
        "Eastern_European_Slavic": 0.30, "Baltic_Finnic": 0.28,
        "Ashkenazi_Jewish": 0.32, "Caucasian_West_Asian": 0.32,
    }),
    # --- BITTER TASTE TAS2R38 ---
    "rs713598": ("G", "C", 1.5, {
        "Northwest_European": 0.47, "Scandinavian": 0.45,
        "Iberian": 0.50, "Basque": 0.48,
        "Italian_Mediterranean": 0.52, "Balkan_Southeast": 0.50,
        "Eastern_European_Slavic": 0.46, "Baltic_Finnic": 0.44,
        "Ashkenazi_Jewish": 0.50, "Caucasian_West_Asian": 0.52,
    }),
    "rs1726866": ("G", "A", 1.5, {
        "Northwest_European": 0.46, "Scandinavian": 0.44,
        "Iberian": 0.48, "Basque": 0.47,
        "Italian_Mediterranean": 0.50, "Balkan_Southeast": 0.49,
        "Eastern_European_Slavic": 0.45, "Baltic_Finnic": 0.43,
        "Ashkenazi_Jewish": 0.49, "Caucasian_West_Asian": 0.50,
    }),
    # --- HLA (diversity marker -- proxy for immune gradient) ---
    # rs9264942: HLA-C haplotype tag.
    "rs9264942": ("T", "C", 1.2, {
        "Northwest_European": 0.28, "Scandinavian": 0.30,
        "Iberian": 0.32, "Basque": 0.30,
        "Italian_Mediterranean": 0.35, "Balkan_Southeast": 0.33,
        "Eastern_European_Slavic": 0.30, "Baltic_Finnic": 0.28,
        "Ashkenazi_Jewish": 0.35, "Caucasian_West_Asian": 0.38,
    }),
    # rs3135391: HLA-DRB1 tag.
    "rs3135391": ("G", "A", 1.2, {
        "Northwest_European": 0.22, "Scandinavian": 0.24,
        "Iberian": 0.18, "Basque": 0.20,
        "Italian_Mediterranean": 0.17, "Balkan_Southeast": 0.19,
        "Eastern_European_Slavic": 0.23, "Baltic_Finnic": 0.28,
        "Ashkenazi_Jewish": 0.25, "Caucasian_West_Asian": 0.15,
    }),
    # --- ASHKENAZI-ENRICHED (founder mutation proxies) ---
    # rs80357906: BRCA1 185delAG surrogate SNP (linked).
    "rs80357906": ("G", "A", 3.0, {
        "Northwest_European": 0.001, "Scandinavian": 0.001,
        "Iberian": 0.001, "Basque": 0.001,
        "Italian_Mediterranean": 0.001, "Balkan_Southeast": 0.001,
        "Eastern_European_Slavic": 0.002, "Baltic_Finnic": 0.001,
        "Ashkenazi_Jewish": 0.010, "Caucasian_West_Asian": 0.001,
    }),
    # rs121908001: HEXA Tay-Sachs proxy SNP.
    "rs121908001": ("G", "A", 2.0, {
        "Northwest_European": 0.001, "Scandinavian": 0.001,
        "Iberian": 0.001, "Basque": 0.001,
        "Italian_Mediterranean": 0.002, "Balkan_Southeast": 0.002,
        "Eastern_European_Slavic": 0.003, "Baltic_Finnic": 0.001,
        "Ashkenazi_Jewish": 0.013, "Caucasian_West_Asian": 0.002,
    }),
    # --- BASQUE/IBERIAN ENRICHED (Rh-negative proxy, pre-IE markers) ---
    # rs590394: Basque-enriched intergenic marker (Lao 2008 Basque component).
    "rs590394": ("T", "C", 2.5, {
        "Northwest_European": 0.30, "Scandinavian": 0.28,
        "Iberian": 0.38, "Basque": 0.55,
        "Italian_Mediterranean": 0.32, "Balkan_Southeast": 0.30,
        "Eastern_European_Slavic": 0.28, "Baltic_Finnic": 0.25,
        "Ashkenazi_Jewish": 0.30, "Caucasian_West_Asian": 0.25,
    }),
    # rs3118914: Basque-enriched marker.
    "rs3118914": ("A", "G", 2.0, {
        "Northwest_European": 0.22, "Scandinavian": 0.20,
        "Iberian": 0.28, "Basque": 0.45,
        "Italian_Mediterranean": 0.22, "Balkan_Southeast": 0.20,
        "Eastern_European_Slavic": 0.20, "Baltic_Finnic": 0.18,
        "Ashkenazi_Jewish": 0.22, "Caucasian_West_Asian": 0.18,
    }),
    # --- SARDINIAN / NEOLITHIC ANATOLIAN PROXIES ---
    # rs11549407: diagnostic for high EEF ancestry (Sardinian/Basque enriched).
    "rs11549407": ("A", "G", 1.5, {
        "Northwest_European": 0.08, "Scandinavian": 0.07,
        "Iberian": 0.12, "Basque": 0.15,
        "Italian_Mediterranean": 0.18, "Balkan_Southeast": 0.13,
        "Eastern_European_Slavic": 0.07, "Baltic_Finnic": 0.05,
        "Ashkenazi_Jewish": 0.10, "Caucasian_West_Asian": 0.14,
    }),
    # --- BALTIC/FINNIC (SIBERIAN ANCESTRY PROXIES) ---
    # rs1024611: CCL2 variant enriched in Finnic populations.
    "rs1024611": ("A", "G", 2.0, {
        "Northwest_European": 0.18, "Scandinavian": 0.20,
        "Iberian": 0.18, "Basque": 0.20,
        "Italian_Mediterranean": 0.22, "Balkan_Southeast": 0.20,
        "Eastern_European_Slavic": 0.22, "Baltic_Finnic": 0.35,
        "Ashkenazi_Jewish": 0.20, "Caucasian_West_Asian": 0.18,
    }),
    # rs2814778 Duffy (control - all European populations are Duffy positive)
    "rs2814778": ("T", "C", 2.5, {
        "Northwest_European": 0.00, "Scandinavian": 0.00,
        "Iberian": 0.01, "Basque": 0.00,
        "Italian_Mediterranean": 0.01, "Balkan_Southeast": 0.01,
        "Eastern_European_Slavic": 0.00, "Baltic_Finnic": 0.00,
        "Ashkenazi_Jewish": 0.01, "Caucasian_West_Asian": 0.02,
    }),
    # --- CYP2D6 drug metabolism ---
    "rs1065852": ("G", "A", 1.5, {
        "Northwest_European": 0.22, "Scandinavian": 0.20,
        "Iberian": 0.18, "Basque": 0.18,
        "Italian_Mediterranean": 0.18, "Balkan_Southeast": 0.20,
        "Eastern_European_Slavic": 0.23, "Baltic_Finnic": 0.25,
        "Ashkenazi_Jewish": 0.20, "Caucasian_West_Asian": 0.20,
    }),
    # --- PRNP prion ---
    "rs1799990": ("A", "G", 1.2, {
        "Northwest_European": 0.37, "Scandinavian": 0.35,
        "Iberian": 0.35, "Basque": 0.33,
        "Italian_Mediterranean": 0.38, "Balkan_Southeast": 0.36,
        "Eastern_European_Slavic": 0.38, "Baltic_Finnic": 0.42,
        "Ashkenazi_Jewish": 0.35, "Caucasian_West_Asian": 0.38,
    }),
    # --- LIPID METABOLISM ---
    # rs662799: APOA5 (Mediterranean protective).
    "rs662799": ("A", "G", 1.3, {
        "Northwest_European": 0.06, "Scandinavian": 0.05,
        "Iberian": 0.08, "Basque": 0.07,
        "Italian_Mediterranean": 0.10, "Balkan_Southeast": 0.09,
        "Eastern_European_Slavic": 0.07, "Baltic_Finnic": 0.06,
        "Ashkenazi_Jewish": 0.08, "Caucasian_West_Asian": 0.11,
    }),
    # rs4880: SOD2.
    "rs4880": ("C", "T", 1.0, {
        "Northwest_European": 0.50, "Scandinavian": 0.50,
        "Iberian": 0.48, "Basque": 0.48,
        "Italian_Mediterranean": 0.46, "Balkan_Southeast": 0.48,
        "Eastern_European_Slavic": 0.50, "Baltic_Finnic": 0.52,
        "Ashkenazi_Jewish": 0.48, "Caucasian_West_Asian": 0.45,
    }),
    # --- G6PD (Mediterranean variant - rs5030868 enriched in southern Eur) ---
    "rs5030868": ("G", "A", 2.0, {
        "Northwest_European": 0.001, "Scandinavian": 0.001,
        "Iberian": 0.005, "Basque": 0.003,
        "Italian_Mediterranean": 0.015, "Balkan_Southeast": 0.018,
        "Eastern_European_Slavic": 0.002, "Baltic_Finnic": 0.001,
        "Ashkenazi_Jewish": 0.010, "Caucasian_West_Asian": 0.012,
    }),
    # --- HBB beta-thalassemia variant (Mediterranean enriched) ---
    "rs33931746": ("C", "T", 2.0, {
        "Northwest_European": 0.001, "Scandinavian": 0.001,
        "Iberian": 0.008, "Basque": 0.003,
        "Italian_Mediterranean": 0.025, "Balkan_Southeast": 0.020,
        "Eastern_European_Slavic": 0.003, "Baltic_Finnic": 0.001,
        "Ashkenazi_Jewish": 0.004, "Caucasian_West_Asian": 0.008,
    }),
    # --- CFTR F508del proxy (Northern European enriched) ---
    "rs113993960": ("A", "T", 2.5, {
        "Northwest_European": 0.025, "Scandinavian": 0.023,
        "Iberian": 0.015, "Basque": 0.018,
        "Italian_Mediterranean": 0.012, "Balkan_Southeast": 0.012,
        "Eastern_European_Slavic": 0.020, "Baltic_Finnic": 0.022,
        "Ashkenazi_Jewish": 0.015, "Caucasian_West_Asian": 0.008,
    }),
    # --- MCM6 context SNP ---
    "rs182549": ("C", "T", 1.8, {
        "Northwest_European": 0.72, "Scandinavian": 0.76,
        "Iberian": 0.42, "Basque": 0.48,
        "Italian_Mediterranean": 0.30, "Balkan_Southeast": 0.35,
        "Eastern_European_Slavic": 0.55, "Baltic_Finnic": 0.80,
        "Ashkenazi_Jewish": 0.35, "Caucasian_West_Asian": 0.25,
    }),
    # --- ABO blood group (weak but informative) ---
    "rs8176719": ("I", "D", 1.0, {
        "Northwest_European": 0.43, "Scandinavian": 0.42,
        "Iberian": 0.45, "Basque": 0.55,  # Basque unique O freq
        "Italian_Mediterranean": 0.45, "Balkan_Southeast": 0.44,
        "Eastern_European_Slavic": 0.42, "Baltic_Finnic": 0.40,
        "Ashkenazi_Jewish": 0.43, "Caucasian_West_Asian": 0.40,
    }),
    # --- RH blood group (Basque rare RH-neg) ---
    "rs590787": ("C", "T", 1.2, {
        "Northwest_European": 0.39, "Scandinavian": 0.38,
        "Iberian": 0.36, "Basque": 0.52,  # Basque >50% Rh-neg
        "Italian_Mediterranean": 0.32, "Balkan_Southeast": 0.32,
        "Eastern_European_Slavic": 0.35, "Baltic_Finnic": 0.36,
        "Ashkenazi_Jewish": 0.30, "Caucasian_West_Asian": 0.25,
    }),
    # --- KITLG pigmentation ---
    "rs12821256": ("T", "C", 1.5, {
        "Northwest_European": 0.14, "Scandinavian": 0.15,
        "Iberian": 0.10, "Basque": 0.12,
        "Italian_Mediterranean": 0.08, "Balkan_Southeast": 0.09,
        "Eastern_European_Slavic": 0.13, "Baltic_Finnic": 0.16,
        "Ashkenazi_Jewish": 0.10, "Caucasian_West_Asian": 0.06,
    }),
    # --- SLC16A12 ---
    "rs1800497": ("G", "A", 1.0, {
        "Northwest_European": 0.20, "Scandinavian": 0.21,
        "Iberian": 0.22, "Basque": 0.22,
        "Italian_Mediterranean": 0.22, "Balkan_Southeast": 0.22,
        "Eastern_European_Slavic": 0.22, "Baltic_Finnic": 0.22,
        "Ashkenazi_Jewish": 0.22, "Caucasian_West_Asian": 0.24,
    }),
    # --- CYP3A5*3 (European high-expresser) ---
    "rs776746": ("C", "T", 1.5, {
        "Northwest_European": 0.08, "Scandinavian": 0.08,
        "Iberian": 0.12, "Basque": 0.10,
        "Italian_Mediterranean": 0.12, "Balkan_Southeast": 0.12,
        "Eastern_European_Slavic": 0.09, "Baltic_Finnic": 0.08,
        "Ashkenazi_Jewish": 0.15, "Caucasian_West_Asian": 0.18,
    }),
    # --- HLA-B27 (ankylosing spondylitis, Northern European enriched) ---
    "rs13202464": ("A", "G", 1.2, {
        "Northwest_European": 0.09, "Scandinavian": 0.12,
        "Iberian": 0.06, "Basque": 0.07,
        "Italian_Mediterranean": 0.05, "Balkan_Southeast": 0.06,
        "Eastern_European_Slavic": 0.08, "Baltic_Finnic": 0.14,
        "Ashkenazi_Jewish": 0.03, "Caucasian_West_Asian": 0.02,
    }),
    # --- BRCA2 European-specific founder proxy ---
    "rs80359550": ("G", "A", 1.8, {
        "Northwest_European": 0.002, "Scandinavian": 0.002,
        "Iberian": 0.001, "Basque": 0.001,
        "Italian_Mediterranean": 0.001, "Balkan_Southeast": 0.001,
        "Eastern_European_Slavic": 0.002, "Baltic_Finnic": 0.002,
        "Ashkenazi_Jewish": 0.012, "Caucasian_West_Asian": 0.001,
    }),
    # --- SEPT9 for Ashkenazi distinction ---
    "rs2070744": ("T", "C", 1.0, {
        "Northwest_European": 0.40, "Scandinavian": 0.42,
        "Iberian": 0.35, "Basque": 0.38,
        "Italian_Mediterranean": 0.33, "Balkan_Southeast": 0.35,
        "Eastern_European_Slavic": 0.38, "Baltic_Finnic": 0.42,
        "Ashkenazi_Jewish": 0.36, "Caucasian_West_Asian": 0.32,
    }),
    # --- TYR R402Q ---
    "rs1126809": ("G", "A", 1.5, {
        "Northwest_European": 0.28, "Scandinavian": 0.30,
        "Iberian": 0.22, "Basque": 0.25,
        "Italian_Mediterranean": 0.18, "Balkan_Southeast": 0.22,
        "Eastern_European_Slavic": 0.28, "Baltic_Finnic": 0.30,
        "Ashkenazi_Jewish": 0.24, "Caucasian_West_Asian": 0.15,
    }),
    # --- OCA2 European eye-color tag ---
    "rs7495174": ("G", "A", 1.2, {
        "Northwest_European": 0.82, "Scandinavian": 0.85,
        "Iberian": 0.70, "Basque": 0.75,
        "Italian_Mediterranean": 0.62, "Balkan_Southeast": 0.68,
        "Eastern_European_Slavic": 0.78, "Baltic_Finnic": 0.84,
        "Ashkenazi_Jewish": 0.72, "Caucasian_West_Asian": 0.55,
    }),
    # --- ANK3 (N European enriched) ---
    "rs10994336": ("C", "T", 1.0, {
        "Northwest_European": 0.08, "Scandinavian": 0.09,
        "Iberian": 0.05, "Basque": 0.06,
        "Italian_Mediterranean": 0.04, "Balkan_Southeast": 0.05,
        "Eastern_European_Slavic": 0.07, "Baltic_Finnic": 0.09,
        "Ashkenazi_Jewish": 0.06, "Caucasian_West_Asian": 0.03,
    }),
    # --- IL7R ---
    "rs6897932": ("C", "T", 1.0, {
        "Northwest_European": 0.25, "Scandinavian": 0.26,
        "Iberian": 0.22, "Basque": 0.24,
        "Italian_Mediterranean": 0.20, "Balkan_Southeast": 0.22,
        "Eastern_European_Slavic": 0.24, "Baltic_Finnic": 0.28,
        "Ashkenazi_Jewish": 0.22, "Caucasian_West_Asian": 0.18,
    }),
    # --- Y-chromosome R1a / R1b proxy (Eastern vs Western Europe) ---
    "rs17250845": ("C", "T", 1.5, {
        "Northwest_European": 0.18, "Scandinavian": 0.22,
        "Iberian": 0.22, "Basque": 0.28,  # R1b-heavy
        "Italian_Mediterranean": 0.17, "Balkan_Southeast": 0.14,
        "Eastern_European_Slavic": 0.38, "Baltic_Finnic": 0.30,  # R1a-heavy
        "Ashkenazi_Jewish": 0.12, "Caucasian_West_Asian": 0.08,
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
    Run weighted NNLS decomposition into European sub-population components.

    Returns (proportions_dict, used_snps_list, residual).
    """
    n_pops = len(EUROPEAN_POPS)
    rows_A = []
    rows_f = []
    rows_w = []
    used_snps = []

    for rsid, (ref, alt, weight, pop_freqs) in AIMS_EUROPEAN.items():
        variant_data = variants.get(rsid.lower())
        if variant_data is None:
            continue
        _chrom, _pos, genotype = variant_data
        dosage = _get_dosage(genotype, alt)
        if dosage is None:
            continue
        obs_freq = dosage / 2.0
        row = [pop_freqs.get(p, 0.0) for p in EUROPEAN_POPS]
        rows_A.append(row)
        rows_f.append(obs_freq)
        rows_w.append(weight)
        used_snps.append(rsid)

    if len(used_snps) < 5:
        return {p: 1.0 / n_pops for p in EUROPEAN_POPS}, used_snps, 0.0

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

    props = {p: float(x[i]) for i, p in enumerate(EUROPEAN_POPS)}
    return props, used_snps, float(resid)


# ---------------------------------------------------------------------------
# AFFINITY SCORES (min-max normalized distances)
# ---------------------------------------------------------------------------

def _compute_affinity_scores(
    variants: Dict[str, Tuple[str, str, str]],
) -> Dict[str, float]:
    """
    Compute relative affinity to each European population using min-max
    normalization. Closest population scores 100, furthest scores 0.
    """
    raw_distances: Dict[str, float] = {}

    for pop in EUROPEAN_POPS:
        weighted_sq_diff = 0.0
        w_sum = 0.0
        for rsid, (ref, alt, weight, pop_freqs) in AIMS_EUROPEAN.items():
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
        return {p: 0.0 for p in EUROPEAN_POPS}

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
# KEY EUROPEAN MARKERS SUMMARY
# ---------------------------------------------------------------------------

def _summarize_key_markers(
    variants: Dict[str, Tuple[str, str, str]],
) -> Dict[str, Dict[str, str]]:
    """Report the key markers distinguishing European sub-populations."""
    results: Dict[str, Dict[str, str]] = {}

    # Lactase persistence rs4988235 (European variant) - N-S gradient
    variant_data = variants.get("rs4988235")
    if variant_data is not None:
        _c, _p, geno = variant_data
        if geno and geno not in ("--", "00", "NN"):
            a_count = geno.upper().count("A")
            if a_count == 2:
                results["lactase"] = {
                    "status": "AA",
                    "detail": "AA (T/T) -- fully lactase persistent. Strongest in "
                              "Northern/Baltic Europeans; lowest in Mediterranean.",
                }
            elif a_count == 1:
                results["lactase"] = {
                    "status": "AG",
                    "detail": "AG (C/T) -- one LP allele, still lactase persistent. "
                              "Intermediate N-S position.",
                }
            else:
                results["lactase"] = {
                    "status": "GG",
                    "detail": "GG (C/C) -- ancestral; lactase non-persistent. "
                              "Southern European, Sardinian, or non-European pattern.",
                }

    # Blue eyes rs12913832 (HERC2) - N-S gradient with Scandinavian peak
    variant_data = variants.get("rs12913832")
    if variant_data is not None:
        _c, _p, geno = variant_data
        if geno and geno not in ("--", "00", "NN"):
            g_count = geno.upper().count("G")
            if g_count == 2:
                results["blue_eyes"] = {
                    "status": "GG",
                    "detail": "GG -- homozygous blue-eye allele. Strongest signal in "
                              "Scandinavian and Baltic-Finnic populations.",
                }
            elif g_count == 1:
                results["blue_eyes"] = {
                    "status": "AG",
                    "detail": "AG -- intermediate eye color; green/hazel/variable. "
                              "Typical Central/Eastern European pattern.",
                }
            else:
                results["blue_eyes"] = {
                    "status": "AA",
                    "detail": "AA -- brown eyes. Typical Mediterranean/Ashkenazi/"
                              "Caucasian pattern.",
                }

    # SLC45A2 rs16891982 - European skin pigmentation (N-S gradient)
    variant_data = variants.get("rs16891982")
    if variant_data is not None:
        _c, _p, geno = variant_data
        if geno and geno not in ("--", "00", "NN"):
            g_count = geno.upper().count("G")
            if g_count == 2:
                results["slc45a2"] = {
                    "status": "GG",
                    "detail": "GG -- homozygous derived (light skin). Very high in "
                              "Northern Europeans.",
                }
            elif g_count == 1:
                results["slc45a2"] = {
                    "status": "CG",
                    "detail": "CG -- heterozygous; Mediterranean or Caucasian pattern.",
                }
            else:
                results["slc45a2"] = {
                    "status": "CC",
                    "detail": "CC -- ancestral; Mediterranean, Middle Eastern, or "
                              "non-European pattern.",
                }

    # HFE C282Y (rs1800562) - Northern European hemochromatosis
    variant_data = variants.get("rs1800562")
    if variant_data is not None:
        _c, _p, geno = variant_data
        if geno and geno not in ("--", "00", "NN"):
            a_count = geno.upper().count("A")
            if a_count == 2:
                results["hfe_c282y"] = {
                    "status": "AA",
                    "detail": "AA -- hemochromatosis homozygote. Very strong "
                              "Northwest European (especially Irish/British) signal. "
                              "Clinically significant iron overload risk.",
                }
            elif a_count == 1:
                results["hfe_c282y"] = {
                    "status": "GA",
                    "detail": "GA -- hemochromatosis carrier. Northern European "
                              "marker.",
                }
            else:
                results["hfe_c282y"] = {
                    "status": "GG",
                    "detail": "GG -- ancestral; no hemochromatosis variant.",
                }

    # CCR5 delta32 (rs333) - Baltic/Northern European HIV resistance
    variant_data = variants.get("rs333")
    if variant_data is not None:
        _c, _p, geno = variant_data
        if geno and geno not in ("--", "00", "NN", "II"):
            d_count = geno.upper().count("D")
            if d_count == 2:
                results["ccr5_delta32"] = {
                    "status": "DD",
                    "detail": "DD -- homozygous CCR5-delta32. HIV resistant. "
                              "Peaks in Baltic/Finnic populations (17%).",
                }
            elif d_count == 1:
                results["ccr5_delta32"] = {
                    "status": "ID",
                    "detail": "ID -- heterozygous CCR5-delta32. Slower HIV progression. "
                              "Northern/Eastern European marker.",
                }
            else:
                results["ccr5_delta32"] = {
                    "status": "II",
                    "detail": "II -- intact CCR5; no delta-32 variant.",
                }

    # Basque ABO/RH (rs590787) - distinctive Basque blood group
    variant_data = variants.get("rs590787")
    if variant_data is not None:
        _c, _p, geno = variant_data
        if geno and geno not in ("--", "00", "NN"):
            t_count = geno.upper().count("T")
            if t_count >= 1:
                results["basque_rh"] = {
                    "status": "T+",
                    "detail": f"T allele present ({'homozygous' if t_count == 2 else 'heterozygous'}) "
                              "-- associated with elevated Rh-negative frequency; "
                              "Basque-enriched signal.",
                }

    # MCM6 context (rs182549) - confirms LP signal
    variant_data = variants.get("rs182549")
    if variant_data is not None:
        _c, _p, geno = variant_data
        if geno and geno not in ("--", "00", "NN"):
            t_count = geno.upper().count("T")
            if t_count == 2:
                results["mcm6_context"] = {
                    "status": "TT",
                    "detail": "TT -- MCM6 context allele supporting LP signal; "
                              "typical Northern European.",
                }

    return results


# ---------------------------------------------------------------------------
# PUBLIC API
# ---------------------------------------------------------------------------

def analyze_european_breakdown(
    variants: Dict[str, Tuple[str, str, str]],
) -> Dict[str, Any]:
    """
    Analyze variants for fine-grained European ancestry breakdown.

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


def generate_european_breakdown_json(result: Dict[str, Any]) -> Dict[str, Any]:
    """
    Convert raw analysis result into the structured JSON for the frontend.

    Args:
        result: Output from analyze_european_breakdown()

    Returns:
        Structured JSON dict matching the report schema.
    """
    props = result["proportions"]
    affinity_scores = result["affinity_scores"]
    used_snps = result["used_snps"]
    resid = result["residual"]
    key_markers = result["key_markers"]

    panel_size = len(AIMS_EUROPEAN)
    snps_used = len(used_snps)

    # Build population list sorted by affinity score descending
    sorted_pops = sorted(EUROPEAN_POPS, key=lambda p: -affinity_scores.get(p, 0.0))

    populations = []
    for pop in sorted_pops:
        proportion = round(props.get(pop, 0.0), 4)
        populations.append({
            "code": pop,
            "label": EUROPEAN_POP_LABELS[pop],
            "affinityScore": affinity_scores.get(pop, 0.0),
            "proportion": proportion,
            "percentage": round(proportion * 100, 1),
            "description": EUROPEAN_POP_DESCRIPTIONS[pop],
            "color": EUROPEAN_POP_COLORS[pop],
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
            "label": EUROPEAN_POP_LABELS[top_pop],
            "description": EUROPEAN_POP_DESCRIPTIONS[top_pop],
        },
        "keyMarkers": key_markers,
    }
