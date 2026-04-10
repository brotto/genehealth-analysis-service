"""
Deep Ancestry / Out-of-Africa Analyzer
Estimates the oldest layers of human ancestry: the first Out-of-Africa
migration (~60,000-100,000 years ago), archaic admixture (Neanderthal and
Denisovan), and the deep Pleistocene lineages that seeded today's
continental populations.

KEY REFERENCES:
  Prufer et al. 2014 (Nature): "The complete genome sequence of a Neanderthal
    from the Altai Mountains"
      -- High-coverage (52x) Altai Neanderthal genome
      -- Estimated ~1.5-2.1% Neanderthal ancestry in non-Africans
  Meyer et al. 2012 (Science): "A high-coverage genome sequence from an
    archaic Denisovan individual"
      -- 31x Denisova cave finger bone
      -- Identified Denisovan introgression in Oceanian populations
  Mallick et al. 2016 (Nature): "The Simons Genome Diversity Project: 300
    genomes from 142 diverse populations"
      -- Global catalog of human genetic variation
      -- Refined estimates of archaic admixture and deep splits
  Sankararaman et al. 2014 (Nature): "The genomic landscape of Neanderthal
    ancestry in present-day humans"
      -- Mapped Neanderthal haplotypes across the human genome
      -- Tracked selection on introgressed tracts

SUPPORTING REFERENCES:
  - Huerta-Sanchez et al. 2014 (Nature): EPAS1 Denisovan introgression in Tibetans
  - Fu et al. 2014 (Nature): Ust-Ishim Siberian (45,000 BP)
  - Raghavan et al. 2014 (Nature): Mal'ta boy / ANE lineage (24,000 BP)
  - Yang et al. 2017 (Curr Biol): Tianyuan East Asian (40,000 BP)
  - Malaspinas et al. 2016 (Nature): Aboriginal Australian / Papuan Denisovan
  - Mondal et al. 2019 (Nat Commun): South Asian Denisovan signal
  - Lipson et al. 2020 (Nature): Shum Laka (W Africa, 8,000 BP)

SCIENTIFIC LIMITATIONS:
  A 40-SNP AIM panel cannot match whole-genome tract-based methods
  (Sankararaman et al.) for precise archaic ancestry estimation. Our
  Neanderthal and Denisovan percentages are coarse proxies -- they should
  be interpreted as low/average/elevated bands rather than exact fractions.
  Similarly, "Basal Eurasian" ancestry (Lazaridis 2014) cannot be cleanly
  isolated from modern reference populations at the individual level.

Returns structured JSON for frontend visualization.
"""

from typing import Dict, Tuple, Any, Optional
import numpy as np
from scipy.optimize import nnls

# ---------------------------------------------------------------------------
# POPULATION DEFINITIONS
# Eight deep-time ancestral lineages spanning ~200,000-15,000 BP
# ---------------------------------------------------------------------------

DEEP_POPS = [
    "Sub_Saharan_African",
    "Khoisan_Deep",
    "Out_of_Africa_Basal",
    "Ancestral_North_Eurasian",
    "East_Eurasian_Core",
    "Papuan_Australasian",
    "Native_American_Deep",
    "Central_Asian_Steppe_Pleistocene",
]

DEEP_POP_LABELS = {
    "Sub_Saharan_African":              "Sub-Saharan African (modern/Bantu)",
    "Khoisan_Deep":                     "Khoisan Deep (oldest modern human lineage)",
    "Out_of_Africa_Basal":              "Basal Eurasian (first Out-of-Africa)",
    "Ancestral_North_Eurasian":         "Ancestral North Eurasian (ANE)",
    "East_Eurasian_Core":               "East Eurasian Core",
    "Papuan_Australasian":              "Papuan / Australasian (Denisovan carriers)",
    "Native_American_Deep":             "Native American Deep (ANE + East Eurasian)",
    "Central_Asian_Steppe_Pleistocene": "Central Asian Steppe Pleistocene",
}

# Colors as requested: deep violet for archaic, cyan for Denisovan,
# dawn red for ANE, jade for E Eurasian, sand for Basal Eurasian, etc.
DEEP_POP_COLORS = {
    "Sub_Saharan_African":              "#3E2723",  # deep brown
    "Khoisan_Deep":                     "#5D4037",  # bark brown
    "Out_of_Africa_Basal":              "#D4A574",  # sand
    "Ancestral_North_Eurasian":         "#B71C1C",  # dawn red
    "East_Eurasian_Core":               "#00695C",  # jade
    "Papuan_Australasian":              "#00BCD4",  # cyan (Denisovan-rich)
    "Native_American_Deep":             "#BF360C",  # Amazonian red-brown
    "Central_Asian_Steppe_Pleistocene": "#4A148C",  # deep violet
}

DEEP_POP_DESCRIPTIONS = {
    "Sub_Saharan_African": (
        "Sub-Saharan African populations (Yoruba, Mandinka, Luhya, Mbuti, and others) "
        "represent the deepest-rooted human genomes and carry the greatest genetic "
        "diversity on Earth. They lack the Neanderthal admixture found in all "
        "non-Africans (though trace amounts have been detected in some East Africans "
        "due to recent Eurasian back-migration). The Bantu expansion (~3,000 BP) "
        "spread the West African gene pool across much of central, eastern and "
        "southern Africa. Reference datasets: 1000 Genomes YRI (Yoruba Ibadan), MSL "
        "(Mende Sierra Leone), GWD (Gambian), and the Simons Genome Diversity Project "
        "African panel (Mallick 2016)."
    ),
    "Khoisan_Deep": (
        "Khoisan populations (Ju/'hoansi, !Xun, Hadza, Sandawe) represent the oldest "
        "surviving branch of modern humans, having diverged from the ancestors of "
        "all other present-day populations around 200,000 years ago (Schlebusch 2017; "
        "Mallick 2016). Their languages contain click consonants and their genomes "
        "carry ancient segments that are unique to this lineage. Khoisan populations "
        "have essentially no Neanderthal ancestry (their ancestors never left Africa) "
        "and preserve very high allelic diversity, including some alleles that have "
        "gone extinct in every other population. A trace Khoisan signal in a genome "
        "typically indicates Southern African ancestry."
    ),
    "Out_of_Africa_Basal": (
        "Basal Eurasian is a cryptic deep lineage first identified by Lazaridis 2014 "
        "in Anatolian and Levantine early farmers. It is a population that split off "
        "from the main Out-of-Africa trunk BEFORE the separation of East Eurasians "
        "and West Eurasians (and crucially, before the main pulse of Neanderthal "
        "admixture). Basal Eurasians therefore have LOWER Neanderthal ancestry than "
        "other non-Africans -- a distinctive signature. The component contributes "
        "substantially to modern Near Easterners (~40-50%), Anatolian farmers "
        "(~45%) and, via EEF, to present-day Europeans. It is hypothesized to have "
        "lived in a Pleistocene refugium in Arabia or the Near East."
    ),
    "Ancestral_North_Eurasian": (
        "Ancestral North Eurasians (ANE) are a late Upper Paleolithic population "
        "best represented by the 'Mal'ta boy' (MA-1) from Siberia, who lived ~24,000 "
        "years ago (Raghavan 2014). ANE-related populations roamed the Siberian and "
        "Central Asian steppes during the Last Glacial Maximum. They contributed: "
        "~40% ancestry to Native Americans (Raghavan 2014); ~50% to Eastern "
        "Hunter-Gatherers, and via EHG to the Yamnaya/Western Steppe Herders "
        "(Haak 2015); and trace amounts (~5-15%) to most modern Europeans. ANE "
        "ancestry is also elevated in South Asians (via the ANI component) and "
        "Caucasians. It is the 'hidden thread' linking the First Americans, "
        "modern Europeans, and Central Asians."
    ),
    "East_Eurasian_Core": (
        "The East Eurasian Core represents the deep ancestral lineage of East and "
        "Southeast Asian populations, best represented by the ~40,000-year-old "
        "Tianyuan individual from China (Yang 2017) and modern Han Chinese, Japanese, "
        "Korean, and Southeast Asian populations. This lineage split from Western "
        "Eurasians early in the Out-of-Africa journey (~45-50kya) and subsequently "
        "differentiated into northern (Amur basin / Baikal) and southern (Yangtze / "
        "Hoabinhian) subclades. East Eurasian ancestry dominates modern East and "
        "Southeast Asia, and -- via mixing with ANE -- contributed roughly half of "
        "the Native American deep component. It is characterized by the EDAR 370A "
        "allele (rs3827760) that affects hair thickness, sweat glands and dentition."
    ),
    "Papuan_Australasian": (
        "Papuan, Aboriginal Australian, Melanesian and some Southeast Asian populations "
        "descend from a very early Out-of-Africa wave that reached Sahul (Australia + "
        "New Guinea) by ~50,000 BP (Malaspinas 2016). These populations carry the "
        "highest known Denisovan ancestry -- ~3-6% -- acquired through admixture with "
        "Denisovans somewhere in Island Southeast Asia (Reich 2010; Meyer 2012). They "
        "represent one of the oldest continuous human occupations outside Africa and "
        "are genetically distinct from mainland East Asians at the deepest level. A "
        "trace Papuan/Australasian signal in a non-Oceanian individual typically "
        "indicates remote Denisovan-enriched ancestry and is often found in Filipinos, "
        "Negritos, Philippine Ayta, and some South Asian groups (Mondal 2019)."
    ),
    "Native_American_Deep": (
        "Native American Deep ancestry formed in Beringia around 24,000-16,000 years "
        "ago from the admixture of Ancestral North Eurasian (~40%) and an East "
        "Eurasian-related lineage (~60%). This composite population crossed into the "
        "Americas between ~16,000 and 13,000 BP and diversified into the two main "
        "subgroups (Northern and Southern Native Americans). Modern Indigenous peoples "
        "of the Americas preserve this unique ANE-East Eurasian signature, with some "
        "populations (e.g. Surui, Karitiana) also showing a cryptic 'Y-signal' shared "
        "with Papuans/Australasians, suggesting an even older wave (Skoglund 2015). "
        "Small fractions of this ancestry appear in modern Siberians and, through "
        "recent admixture, in Mestizo/Hispanic populations."
    ),
    "Central_Asian_Steppe_Pleistocene": (
        "Deep Pleistocene Central Asian Steppe populations, best represented by the "
        "Ust-Ishim individual (Western Siberia, ~45,000 BP; Fu 2014) and the Kostenki "
        "14 individual (Russia, ~37,000 BP; Seguin-Orlando 2014). These are the "
        "oldest modern human genomes outside Africa and the Near East. They lived "
        "before the major splits between Europeans and East Asians had fully formed "
        "and carried slightly longer Neanderthal introgression tracts than modern "
        "humans -- indicating the Neanderthal admixture event happened only a few "
        "thousand years before. Trace signals of this deep Pleistocene component in "
        "modern genomes indicate unusually deep Eurasian ancestry, sometimes seen "
        "in Siberians, northern Europeans and some East Asian groups."
    ),
}

# ---------------------------------------------------------------------------
# AIM ALLELE FREQUENCIES FOR DEEP ANCESTRAL LINEAGES
# Calibrated from: Prufer 2014, Meyer 2012, Mallick 2016, Sankararaman 2014,
# 1000 Genomes Phase 3, SGDP, HGDP, Raghavan 2014, Fu 2014/2016.
#
# Format: rsid: (ref, alt, weight, {pop: alt_freq})
# Weights reflect discrimination power: Neanderthal-tagging and Denisovan-
# introgressed alleles get the highest weights.
# ---------------------------------------------------------------------------

AIMS_DEEP = {
    # --- NEANDERTHAL-INTROGRESSED / TAG VARIANTS ---
    # BNC2 rs10946808 -- Neanderthal-derived allele (skin color gene)
    "rs10946808": ("G", "A", 4.0, {
        "Sub_Saharan_African":              0.00,
        "Khoisan_Deep":                     0.00,
        "Out_of_Africa_Basal":              0.05,
        "Ancestral_North_Eurasian":         0.40,
        "East_Eurasian_Core":               0.35,
        "Papuan_Australasian":              0.40,
        "Native_American_Deep":             0.38,
        "Central_Asian_Steppe_Pleistocene": 0.45,
    }),
    # POU2F3 rs10831496 -- Neanderthal-introgressed
    "rs10831496": ("C", "T", 3.5, {
        "Sub_Saharan_African":              0.00,
        "Khoisan_Deep":                     0.00,
        "Out_of_Africa_Basal":              0.10,
        "Ancestral_North_Eurasian":         0.35,
        "East_Eurasian_Core":               0.55,
        "Papuan_Australasian":              0.35,
        "Native_American_Deep":             0.40,
        "Central_Asian_Steppe_Pleistocene": 0.40,
    }),
    # HYAL2 rs12488302 -- Neanderthal-derived, enriched in Eurasians
    "rs12488302": ("G", "A", 3.5, {
        "Sub_Saharan_African":              0.00,
        "Khoisan_Deep":                     0.00,
        "Out_of_Africa_Basal":              0.08,
        "Ancestral_North_Eurasian":         0.35,
        "East_Eurasian_Core":               0.40,
        "Papuan_Australasian":              0.30,
        "Native_American_Deep":             0.35,
        "Central_Asian_Steppe_Pleistocene": 0.38,
    }),
    # OAS1 rs1131454 -- Neanderthal-introgressed immune gene
    "rs1131454": ("A", "G", 3.0, {
        "Sub_Saharan_African":              0.00,
        "Khoisan_Deep":                     0.00,
        "Out_of_Africa_Basal":              0.20,
        "Ancestral_North_Eurasian":         0.42,
        "East_Eurasian_Core":               0.45,
        "Papuan_Australasian":              0.45,
        "Native_American_Deep":             0.45,
        "Central_Asian_Steppe_Pleistocene": 0.48,
    }),
    # SLC16A11 rs13342232 -- Neanderthal-introgressed, Native American enriched
    "rs13342232": ("C", "T", 3.5, {
        "Sub_Saharan_African":              0.00,
        "Khoisan_Deep":                     0.00,
        "Out_of_Africa_Basal":              0.02,
        "Ancestral_North_Eurasian":         0.10,
        "East_Eurasian_Core":               0.08,
        "Papuan_Australasian":              0.05,
        "Native_American_Deep":             0.45,
        "Central_Asian_Steppe_Pleistocene": 0.10,
    }),
    # STAT2 rs2066807 -- Neanderthal-introgressed immune
    "rs2066807": ("C", "G", 3.0, {
        "Sub_Saharan_African":              0.00,
        "Khoisan_Deep":                     0.00,
        "Out_of_Africa_Basal":              0.02,
        "Ancestral_North_Eurasian":         0.05,
        "East_Eurasian_Core":               0.08,
        "Papuan_Australasian":              0.55,
        "Native_American_Deep":             0.05,
        "Central_Asian_Steppe_Pleistocene": 0.05,
    }),
    # SPPL2C rs4792897 -- Neanderthal allele near MAPT
    "rs4792897": ("C", "T", 2.5, {
        "Sub_Saharan_African":              0.00,
        "Khoisan_Deep":                     0.00,
        "Out_of_Africa_Basal":              0.05,
        "Ancestral_North_Eurasian":         0.15,
        "East_Eurasian_Core":               0.18,
        "Papuan_Australasian":              0.20,
        "Native_American_Deep":             0.18,
        "Central_Asian_Steppe_Pleistocene": 0.20,
    }),
    # TLR1/6/10 rs5743618 -- Neanderthal-introgressed immune cluster
    "rs5743618": ("T", "G", 2.5, {
        "Sub_Saharan_African":              0.00,
        "Khoisan_Deep":                     0.00,
        "Out_of_Africa_Basal":              0.15,
        "Ancestral_North_Eurasian":         0.45,
        "East_Eurasian_Core":               0.40,
        "Papuan_Australasian":              0.45,
        "Native_American_Deep":             0.50,
        "Central_Asian_Steppe_Pleistocene": 0.50,
    }),

    # --- DENISOVAN-INTROGRESSED ---
    # EPAS1 rs13419896 -- famous Denisovan high-altitude allele in Tibetans
    "rs13419896": ("A", "G", 5.0, {
        "Sub_Saharan_African":              0.00,
        "Khoisan_Deep":                     0.00,
        "Out_of_Africa_Basal":              0.00,
        "Ancestral_North_Eurasian":         0.02,
        "East_Eurasian_Core":               0.15,  # elevated in Tibetans/Sherpas
        "Papuan_Australasian":              0.10,
        "Native_American_Deep":             0.03,
        "Central_Asian_Steppe_Pleistocene": 0.02,
    }),
    # TBX15/WARS2 rs2298080 -- Denisovan body-shape allele (high in Inuit)
    "rs2298080": ("G", "A", 3.0, {
        "Sub_Saharan_African":              0.00,
        "Khoisan_Deep":                     0.00,
        "Out_of_Africa_Basal":              0.05,
        "Ancestral_North_Eurasian":         0.25,
        "East_Eurasian_Core":               0.30,
        "Papuan_Australasian":              0.35,
        "Native_American_Deep":             0.45,
        "Central_Asian_Steppe_Pleistocene": 0.30,
    }),
    # WARS2 rs2001207 -- Denisovan-tagged
    "rs2001207": ("T", "C", 2.5, {
        "Sub_Saharan_African":              0.00,
        "Khoisan_Deep":                     0.00,
        "Out_of_Africa_Basal":              0.03,
        "Ancestral_North_Eurasian":         0.15,
        "East_Eurasian_Core":               0.20,
        "Papuan_Australasian":              0.40,
        "Native_American_Deep":             0.20,
        "Central_Asian_Steppe_Pleistocene": 0.15,
    }),

    # --- AFRICAN vs NON-AFRICAN DISCRIMINATION ---
    # DARC rs2814778 -- Duffy null; near-fixed in Sub-Saharan Africa
    "rs2814778": ("T", "C", 4.5, {
        "Sub_Saharan_African":              0.95,
        "Khoisan_Deep":                     0.40,
        "Out_of_Africa_Basal":              0.05,
        "Ancestral_North_Eurasian":         0.00,
        "East_Eurasian_Core":               0.00,
        "Papuan_Australasian":              0.00,
        "Native_American_Deep":             0.00,
        "Central_Asian_Steppe_Pleistocene": 0.00,
    }),
    # HBB rs334 -- sickle cell (African-enriched)
    "rs334": ("A", "T", 3.0, {
        "Sub_Saharan_African":              0.12,
        "Khoisan_Deep":                     0.00,
        "Out_of_Africa_Basal":              0.02,
        "Ancestral_North_Eurasian":         0.00,
        "East_Eurasian_Core":               0.00,
        "Papuan_Australasian":              0.00,
        "Native_American_Deep":             0.00,
        "Central_Asian_Steppe_Pleistocene": 0.00,
    }),
    # APOL1 rs73885319 -- African kidney variant
    "rs73885319": ("A", "G", 3.0, {
        "Sub_Saharan_African":              0.20,
        "Khoisan_Deep":                     0.01,
        "Out_of_Africa_Basal":              0.00,
        "Ancestral_North_Eurasian":         0.00,
        "East_Eurasian_Core":               0.00,
        "Papuan_Australasian":              0.00,
        "Native_American_Deep":             0.00,
        "Central_Asian_Steppe_Pleistocene": 0.00,
    }),
    # LCT rs145946881 -- African-specific lactase variant
    "rs145946881": ("C", "G", 2.5, {
        "Sub_Saharan_African":              0.18,
        "Khoisan_Deep":                     0.00,
        "Out_of_Africa_Basal":              0.00,
        "Ancestral_North_Eurasian":         0.00,
        "East_Eurasian_Core":               0.00,
        "Papuan_Australasian":              0.00,
        "Native_American_Deep":             0.00,
        "Central_Asian_Steppe_Pleistocene": 0.00,
    }),

    # --- EAST EURASIAN MARKERS ---
    # EDAR rs3827760 -- classic East Asian / absent elsewhere
    "rs3827760": ("A", "G", 4.5, {
        "Sub_Saharan_African":              0.00,
        "Khoisan_Deep":                     0.00,
        "Out_of_Africa_Basal":              0.00,
        "Ancestral_North_Eurasian":         0.05,
        "East_Eurasian_Core":               0.85,
        "Papuan_Australasian":              0.10,
        "Native_American_Deep":             0.70,
        "Central_Asian_Steppe_Pleistocene": 0.35,
    }),
    # ABCC11 rs17822931 -- dry earwax / East Asian
    "rs17822931": ("C", "T", 3.5, {
        "Sub_Saharan_African":              0.00,
        "Khoisan_Deep":                     0.00,
        "Out_of_Africa_Basal":              0.05,
        "Ancestral_North_Eurasian":         0.15,
        "East_Eurasian_Core":               0.80,
        "Papuan_Australasian":              0.10,
        "Native_American_Deep":             0.15,
        "Central_Asian_Steppe_Pleistocene": 0.25,
    }),
    # ALDH2 rs671 -- Asian alcohol flush (East Asian exclusive)
    "rs671": ("G", "A", 3.5, {
        "Sub_Saharan_African":              0.00,
        "Khoisan_Deep":                     0.00,
        "Out_of_Africa_Basal":              0.00,
        "Ancestral_North_Eurasian":         0.00,
        "East_Eurasian_Core":               0.25,
        "Papuan_Australasian":              0.00,
        "Native_American_Deep":             0.00,
        "Central_Asian_Steppe_Pleistocene": 0.05,
    }),
    # ADH1B rs1229984 -- East Asian alcohol dehydrogenase
    "rs1229984": ("G", "A", 3.0, {
        "Sub_Saharan_African":              0.00,
        "Khoisan_Deep":                     0.00,
        "Out_of_Africa_Basal":              0.03,
        "Ancestral_North_Eurasian":         0.02,
        "East_Eurasian_Core":               0.65,
        "Papuan_Australasian":              0.05,
        "Native_American_Deep":             0.05,
        "Central_Asian_Steppe_Pleistocene": 0.10,
    }),

    # --- ANE / BASAL EURASIAN PROXIES ---
    # KITLG rs12821256 -- blond hair, ANE/European
    "rs12821256": ("T", "C", 2.5, {
        "Sub_Saharan_African":              0.00,
        "Khoisan_Deep":                     0.00,
        "Out_of_Africa_Basal":              0.02,
        "Ancestral_North_Eurasian":         0.15,
        "East_Eurasian_Core":               0.00,
        "Papuan_Australasian":              0.00,
        "Native_American_Deep":             0.00,
        "Central_Asian_Steppe_Pleistocene": 0.10,
    }),
    # TLR5 rs5744174 -- differentiates Basal Eurasian / non-African
    "rs5744174": ("T", "C", 1.5, {
        "Sub_Saharan_African":              0.25,
        "Khoisan_Deep":                     0.20,
        "Out_of_Africa_Basal":              0.45,
        "Ancestral_North_Eurasian":         0.40,
        "East_Eurasian_Core":               0.50,
        "Papuan_Australasian":              0.40,
        "Native_American_Deep":             0.45,
        "Central_Asian_Steppe_Pleistocene": 0.40,
    }),
    # HERC2 rs12913832 -- blue eyes, ANE/West Eurasian
    "rs12913832": ("A", "G", 3.0, {
        "Sub_Saharan_African":              0.00,
        "Khoisan_Deep":                     0.00,
        "Out_of_Africa_Basal":              0.02,
        "Ancestral_North_Eurasian":         0.30,
        "East_Eurasian_Core":               0.00,
        "Papuan_Australasian":              0.00,
        "Native_American_Deep":             0.00,
        "Central_Asian_Steppe_Pleistocene": 0.15,
    }),

    # --- PIGMENTATION GRADIENTS ---
    # SLC24A5 rs1426654 -- lighter-skin Out-of-Africa variant
    "rs1426654": ("A", "G", 3.5, {
        "Sub_Saharan_African":              0.02,
        "Khoisan_Deep":                     0.05,
        "Out_of_Africa_Basal":              0.90,
        "Ancestral_North_Eurasian":         0.60,
        "East_Eurasian_Core":               0.02,
        "Papuan_Australasian":              0.02,
        "Native_American_Deep":             0.10,
        "Central_Asian_Steppe_Pleistocene": 0.45,
    }),
    # SLC45A2 rs16891982 -- West Eurasian lighter skin
    "rs16891982": ("C", "G", 2.5, {
        "Sub_Saharan_African":              0.00,
        "Khoisan_Deep":                     0.00,
        "Out_of_Africa_Basal":              0.25,
        "Ancestral_North_Eurasian":         0.15,
        "East_Eurasian_Core":               0.00,
        "Papuan_Australasian":              0.00,
        "Native_American_Deep":             0.00,
        "Central_Asian_Steppe_Pleistocene": 0.10,
    }),
    # MC1R rs885479 -- East Asian skin color variant
    "rs885479": ("G", "A", 2.0, {
        "Sub_Saharan_African":              0.00,
        "Khoisan_Deep":                     0.00,
        "Out_of_Africa_Basal":              0.00,
        "Ancestral_North_Eurasian":         0.02,
        "East_Eurasian_Core":               0.55,
        "Papuan_Australasian":              0.05,
        "Native_American_Deep":             0.08,
        "Central_Asian_Steppe_Pleistocene": 0.15,
    }),

    # --- IMMUNE / DISEASE ---
    # G6PD rs1050828 -- African G6PD deficiency
    "rs1050828": ("C", "T", 2.0, {
        "Sub_Saharan_African":              0.15,
        "Khoisan_Deep":                     0.00,
        "Out_of_Africa_Basal":              0.02,
        "Ancestral_North_Eurasian":         0.00,
        "East_Eurasian_Core":               0.00,
        "Papuan_Australasian":              0.00,
        "Native_American_Deep":             0.00,
        "Central_Asian_Steppe_Pleistocene": 0.00,
    }),
    # CCR5 rs1800023 -- Eurasian enrichment
    "rs1800023": ("G", "A", 1.5, {
        "Sub_Saharan_African":              0.02,
        "Khoisan_Deep":                     0.00,
        "Out_of_Africa_Basal":              0.10,
        "Ancestral_North_Eurasian":         0.12,
        "East_Eurasian_Core":               0.08,
        "Papuan_Australasian":              0.08,
        "Native_American_Deep":             0.10,
        "Central_Asian_Steppe_Pleistocene": 0.12,
    }),
    # FUT2 rs601338 -- secretor status
    "rs601338": ("G", "A", 1.5, {
        "Sub_Saharan_African":              0.40,
        "Khoisan_Deep":                     0.35,
        "Out_of_Africa_Basal":              0.40,
        "Ancestral_North_Eurasian":         0.50,
        "East_Eurasian_Core":               0.05,
        "Papuan_Australasian":              0.30,
        "Native_American_Deep":             0.25,
        "Central_Asian_Steppe_Pleistocene": 0.40,
    }),

    # --- Y-HAPLOGROUP DEEP LINEAGE PROXIES (autosomal tags) ---
    # rs3911 -- linked to Y-haplogroup C-M130 clade (Australasian/Oceanian)
    "rs3911": ("A", "G", 1.2, {
        "Sub_Saharan_African":              0.05,
        "Khoisan_Deep":                     0.05,
        "Out_of_Africa_Basal":              0.15,
        "Ancestral_North_Eurasian":         0.20,
        "East_Eurasian_Core":               0.30,
        "Papuan_Australasian":              0.50,
        "Native_American_Deep":             0.20,
        "Central_Asian_Steppe_Pleistocene": 0.25,
    }),
    # rs2032599 -- African deep Y-haplogroup tag
    "rs2032599": ("C", "T", 1.0, {
        "Sub_Saharan_African":              0.40,
        "Khoisan_Deep":                     0.55,
        "Out_of_Africa_Basal":              0.15,
        "Ancestral_North_Eurasian":         0.10,
        "East_Eurasian_Core":               0.10,
        "Papuan_Australasian":              0.10,
        "Native_American_Deep":             0.10,
        "Central_Asian_Steppe_Pleistocene": 0.10,
    }),

    # --- mtDNA deep haplogroup proxies (autosomal linked) ---
    # rs28359168 -- informative for L0/L1 (African) vs M/N (non-African) proxy
    "rs28359168": ("A", "G", 1.2, {
        "Sub_Saharan_African":              0.65,
        "Khoisan_Deep":                     0.80,
        "Out_of_Africa_Basal":              0.15,
        "Ancestral_North_Eurasian":         0.10,
        "East_Eurasian_Core":               0.10,
        "Papuan_Australasian":              0.12,
        "Native_American_Deep":             0.10,
        "Central_Asian_Steppe_Pleistocene": 0.10,
    }),

    # --- NEUTRAL AIMs WITH STRONG CONTINENTAL SIGNAL ---
    "rs260690": ("T", "C", 1.5, {
        "Sub_Saharan_African":              0.60,
        "Khoisan_Deep":                     0.55,
        "Out_of_Africa_Basal":              0.30,
        "Ancestral_North_Eurasian":         0.40,
        "East_Eurasian_Core":               0.10,
        "Papuan_Australasian":              0.20,
        "Native_American_Deep":             0.15,
        "Central_Asian_Steppe_Pleistocene": 0.30,
    }),
    "rs1834640": ("A", "G", 1.5, {
        "Sub_Saharan_African":              0.70,
        "Khoisan_Deep":                     0.65,
        "Out_of_Africa_Basal":              0.55,
        "Ancestral_North_Eurasian":         0.35,
        "East_Eurasian_Core":               0.50,
        "Papuan_Australasian":              0.40,
        "Native_American_Deep":             0.25,
        "Central_Asian_Steppe_Pleistocene": 0.40,
    }),
    "rs3811801": ("G", "A", 2.0, {
        "Sub_Saharan_African":              0.00,
        "Khoisan_Deep":                     0.00,
        "Out_of_Africa_Basal":              0.00,
        "Ancestral_North_Eurasian":         0.02,
        "East_Eurasian_Core":               0.70,
        "Papuan_Australasian":              0.00,
        "Native_American_Deep":             0.05,
        "Central_Asian_Steppe_Pleistocene": 0.15,
    }),
    "rs3916235": ("A", "G", 1.0, {
        "Sub_Saharan_African":              0.30,
        "Khoisan_Deep":                     0.30,
        "Out_of_Africa_Basal":              0.25,
        "Ancestral_North_Eurasian":         0.35,
        "East_Eurasian_Core":               0.45,
        "Papuan_Australasian":              0.30,
        "Native_American_Deep":             0.40,
        "Central_Asian_Steppe_Pleistocene": 0.35,
    }),
    "rs6548238": ("C", "T", 1.0, {
        "Sub_Saharan_African":              0.05,
        "Khoisan_Deep":                     0.02,
        "Out_of_Africa_Basal":              0.30,
        "Ancestral_North_Eurasian":         0.22,
        "East_Eurasian_Core":               0.15,
        "Papuan_Australasian":              0.10,
        "Native_American_Deep":             0.15,
        "Central_Asian_Steppe_Pleistocene": 0.20,
    }),
    "rs7349332": ("G", "A", 1.0, {
        "Sub_Saharan_African":              0.40,
        "Khoisan_Deep":                     0.35,
        "Out_of_Africa_Basal":              0.45,
        "Ancestral_North_Eurasian":         0.50,
        "East_Eurasian_Core":               0.35,
        "Papuan_Australasian":              0.40,
        "Native_American_Deep":             0.30,
        "Central_Asian_Steppe_Pleistocene": 0.40,
    }),
    "rs2430561": ("T", "A", 1.0, {
        "Sub_Saharan_African":              0.45,
        "Khoisan_Deep":                     0.40,
        "Out_of_Africa_Basal":              0.30,
        "Ancestral_North_Eurasian":         0.35,
        "East_Eurasian_Core":               0.25,
        "Papuan_Australasian":              0.20,
        "Native_American_Deep":             0.25,
        "Central_Asian_Steppe_Pleistocene": 0.30,
    }),
    "rs1408801": ("G", "A", 1.0, {
        "Sub_Saharan_African":              0.25,
        "Khoisan_Deep":                     0.25,
        "Out_of_Africa_Basal":              0.40,
        "Ancestral_North_Eurasian":         0.30,
        "East_Eurasian_Core":               0.45,
        "Papuan_Australasian":              0.35,
        "Native_American_Deep":             0.40,
        "Central_Asian_Steppe_Pleistocene": 0.35,
    }),
}

# Total Neanderthal-tagging SNPs and Denisovan-tagging SNPs used for
# dedicated percentage estimates (separate from the NNLS decomposition).
NEANDERTHAL_TAG_SNPS = [
    "rs10946808", "rs10831496", "rs12488302", "rs1131454",
    "rs13342232", "rs2066807", "rs4792897", "rs5743618",
]

DENISOVAN_TAG_SNPS = [
    "rs13419896", "rs2298080", "rs2001207",
]

# Reference frequency assumed for a "typical non-African" individual
# (calibrated to the ~1.8% published mean for Eurasians; used for scaling)
NEANDERTHAL_EURASIAN_REF = 1.8   # percent
DENISOVAN_OCEANIAN_REF = 4.5     # percent, Papuan/Melanesian reference


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
# NNLS ANCESTRY ESTIMATION
# ---------------------------------------------------------------------------

def _estimate_proportions(
    variants: Dict[str, Tuple[str, str, str]],
) -> Tuple[Dict[str, float], list, float]:
    """
    Run weighted NNLS decomposition into deep ancestral components.
    Returns (proportions_dict, used_snps_list, residual).
    """
    n_pops = len(DEEP_POPS)
    rows_A = []
    rows_f = []
    rows_w = []
    used_snps = []

    for rsid, (ref, alt, weight, pop_freqs) in AIMS_DEEP.items():
        variant_data = variants.get(rsid.lower())
        if variant_data is None:
            continue
        _chrom, _pos, genotype = variant_data
        dosage = _get_dosage(genotype, alt)
        if dosage is None:
            continue
        obs_freq = dosage / 2.0
        row = [pop_freqs.get(p, 0.0) for p in DEEP_POPS]
        rows_A.append(row)
        rows_f.append(obs_freq)
        rows_w.append(weight)
        used_snps.append(rsid)

    if len(used_snps) < 5:
        return {p: 1.0 / n_pops for p in DEEP_POPS}, used_snps, 0.0

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

    props = {p: float(x[i]) for i, p in enumerate(DEEP_POPS)}
    return props, used_snps, float(resid)


# ---------------------------------------------------------------------------
# AFFINITY SCORES
# ---------------------------------------------------------------------------

def _compute_affinity_scores(
    variants: Dict[str, Tuple[str, str, str]],
) -> Dict[str, float]:
    """Min-max normalized affinity to each deep lineage."""
    raw_distances: Dict[str, float] = {}

    for pop in DEEP_POPS:
        weighted_sq_diff = 0.0
        w_sum = 0.0
        for rsid, (ref, alt, weight, pop_freqs) in AIMS_DEEP.items():
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
        return {p: 0.0 for p in DEEP_POPS}

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
# ARCHAIC ADMIXTURE ESTIMATES (Neanderthal / Denisovan)
# ---------------------------------------------------------------------------

def _estimate_neanderthal_percentage(
    variants: Dict[str, Tuple[str, str, str]],
) -> float:
    """
    Estimate Neanderthal ancestry percentage from introgressed tag SNPs.
    Calibrated so that an 'average non-African' individual (dosage ratio
    ~1.0 vs the Eurasian reference) returns ~1.8%.
    """
    observed_sum = 0.0
    reference_sum = 0.0

    for rsid in NEANDERTHAL_TAG_SNPS:
        if rsid not in AIMS_DEEP:
            continue
        ref, alt, weight, pop_freqs = AIMS_DEEP[rsid]
        variant_data = variants.get(rsid.lower())
        if variant_data is None:
            continue
        _c, _p, genotype = variant_data
        dosage = _get_dosage(genotype, alt)
        if dosage is None:
            continue

        # Use an average of non-African populations as the Eurasian reference
        ref_freq = np.mean([
            pop_freqs["Out_of_Africa_Basal"],
            pop_freqs["Ancestral_North_Eurasian"],
            pop_freqs["East_Eurasian_Core"],
            pop_freqs["Papuan_Australasian"],
            pop_freqs["Native_American_Deep"],
            pop_freqs["Central_Asian_Steppe_Pleistocene"],
        ])
        observed_sum += dosage / 2.0
        reference_sum += ref_freq

    if reference_sum <= 0:
        return 0.0

    ratio = observed_sum / reference_sum
    percentage = round(ratio * NEANDERTHAL_EURASIAN_REF, 2)
    return max(0.0, min(percentage, 5.0))


def _estimate_denisovan_percentage(
    variants: Dict[str, Tuple[str, str, str]],
) -> float:
    """
    Estimate Denisovan ancestry percentage from introgressed tag SNPs.
    Calibrated so an average Papuan individual returns ~4.5% and most
    non-Oceanians return <0.2%.
    """
    observed_sum = 0.0
    reference_sum = 0.0

    for rsid in DENISOVAN_TAG_SNPS:
        if rsid not in AIMS_DEEP:
            continue
        ref, alt, weight, pop_freqs = AIMS_DEEP[rsid]
        variant_data = variants.get(rsid.lower())
        if variant_data is None:
            continue
        _c, _p, genotype = variant_data
        dosage = _get_dosage(genotype, alt)
        if dosage is None:
            continue

        ref_freq = pop_freqs["Papuan_Australasian"]
        observed_sum += dosage / 2.0
        reference_sum += ref_freq

    if reference_sum <= 0:
        return 0.0

    ratio = observed_sum / reference_sum
    percentage = round(ratio * DENISOVAN_OCEANIAN_REF, 2)
    return max(0.0, min(percentage, 6.5))


# ---------------------------------------------------------------------------
# KEY MARKERS SUMMARY
# ---------------------------------------------------------------------------

def _summarize_key_markers(
    variants: Dict[str, Tuple[str, str, str]],
) -> Dict[str, Dict[str, str]]:
    """Report key deep-ancestry markers: Neanderthal, Denisovan, ANE, Basal."""
    results: Dict[str, Dict[str, str]] = {}

    # --- Neanderthal-introgressed BNC2 ---
    variant_data = variants.get("rs10946808")
    if variant_data is not None:
        _c, _p, geno = variant_data
        if geno and geno not in ("--", "00", "NN"):
            a_count = geno.upper().count("A")
            if a_count == 2:
                results["neanderthal_bnc2"] = {
                    "status": "AA",
                    "detail": (
                        "AA -- carrying two Neanderthal-introgressed alleles in BNC2 "
                        "(a skin color gene). Typical of individuals with elevated "
                        "Neanderthal tract density."
                    ),
                }
            elif a_count == 1:
                results["neanderthal_bnc2"] = {
                    "status": "GA",
                    "detail": (
                        "GA -- one Neanderthal-introgressed BNC2 allele. Common in "
                        "Eurasians (the derived allele is absent in Sub-Saharan Africans)."
                    ),
                }
            else:
                results["neanderthal_bnc2"] = {
                    "status": "GG",
                    "detail": (
                        "GG -- ancestral. No Neanderthal-derived BNC2 allele. Common "
                        "in Africans and in Eurasians who did not inherit this tract."
                    ),
                }

    # --- Denisovan EPAS1 high-altitude allele ---
    variant_data = variants.get("rs13419896")
    if variant_data is not None:
        _c, _p, geno = variant_data
        if geno and geno not in ("--", "00", "NN"):
            g_count = geno.upper().count("G")
            if g_count >= 1:
                results["denisovan_epas1"] = {
                    "status": "AG" if g_count == 1 else "GG",
                    "detail": (
                        "Denisovan-derived EPAS1 allele present. This famous allele "
                        "(Huerta-Sanchez 2014) is nearly exclusive to Tibetan and "
                        "Sherpa populations and confers improved high-altitude "
                        "adaptation. Very rare in other groups."
                    ),
                }
            else:
                results["denisovan_epas1"] = {
                    "status": "AA",
                    "detail": (
                        "AA -- ancestral. Does not carry the Denisovan-introgressed "
                        "EPAS1 high-altitude allele (the standard state outside Tibet)."
                    ),
                }

    # --- DARC rs2814778 (Duffy null, African marker) ---
    variant_data = variants.get("rs2814778")
    if variant_data is not None:
        _c, _p, geno = variant_data
        if geno and geno not in ("--", "00", "NN"):
            c_count = geno.upper().count("C")
            if c_count == 2:
                results["darc_duffy"] = {
                    "status": "CC",
                    "detail": (
                        "CC -- Duffy-null. This allele is near-fixed in Sub-Saharan "
                        "Africans (malaria resistance via P. vivax) and extremely "
                        "rare elsewhere. A strong African-ancestry signature."
                    ),
                }
            elif c_count == 1:
                results["darc_duffy"] = {
                    "status": "TC",
                    "detail": (
                        "TC -- one Duffy-null allele. Indicates partial Sub-Saharan "
                        "African ancestry."
                    ),
                }
            else:
                results["darc_duffy"] = {
                    "status": "TT",
                    "detail": (
                        "TT -- Duffy-positive (ancestral). Standard state outside "
                        "Sub-Saharan Africa."
                    ),
                }

    # --- EDAR rs3827760 (East Asian marker, absent elsewhere) ---
    variant_data = variants.get("rs3827760")
    if variant_data is not None:
        _c, _p, geno = variant_data
        if geno and geno not in ("--", "00", "NN"):
            g_count = geno.upper().count("G")
            if g_count == 2:
                results["edar_east_asian"] = {
                    "status": "GG",
                    "detail": (
                        "GG -- East Asian EDAR 370A allele homozygous. Affects hair "
                        "thickness, sweat glands, and shovel-shaped incisors. "
                        "Strong East Eurasian / Native American signature."
                    ),
                }
            elif g_count == 1:
                results["edar_east_asian"] = {
                    "status": "AG",
                    "detail": (
                        "AG -- one East Asian EDAR allele. Suggests partial East "
                        "Eurasian or Native American ancestry."
                    ),
                }
            else:
                results["edar_east_asian"] = {
                    "status": "AA",
                    "detail": (
                        "AA -- ancestral. Standard state in Africans, Europeans, "
                        "Middle Easterners, and South Asians."
                    ),
                }

    # --- SLC16A11 rs13342232 (Native American Neanderthal signal) ---
    variant_data = variants.get("rs13342232")
    if variant_data is not None:
        _c, _p, geno = variant_data
        if geno and geno not in ("--", "00", "NN"):
            t_count = geno.upper().count("T")
            if t_count >= 1:
                status = "CT" if t_count == 1 else "TT"
                results["slc16a11_amerindian"] = {
                    "status": status,
                    "detail": (
                        "Carries the Neanderthal-introgressed SLC16A11 risk allele, "
                        "present at ~45-50% frequency in Indigenous Americans and "
                        "Mestizos but rare elsewhere. A characteristic Native "
                        "American deep signature."
                    ),
                }
            else:
                results["slc16a11_amerindian"] = {
                    "status": "CC",
                    "detail": (
                        "CC -- ancestral. Standard in non-Amerindian populations."
                    ),
                }

    # --- HERC2 rs12913832 (ANE/European blue eyes) ---
    variant_data = variants.get("rs12913832")
    if variant_data is not None:
        _c, _p, geno = variant_data
        if geno and geno not in ("--", "00", "NN"):
            g_count = geno.upper().count("G")
            if g_count >= 1:
                results["ane_european"] = {
                    "status": "AG" if g_count == 1 else "GG",
                    "detail": (
                        "Carries the European/ANE blue-eye allele at HERC2. "
                        "This allele traces back to Western Hunter-Gatherers and "
                        "is absent in Africans and East Asians."
                    ),
                }
            else:
                results["ane_european"] = {
                    "status": "AA",
                    "detail": (
                        "AA -- ancestral. Brown eyes; standard in non-European "
                        "populations worldwide."
                    ),
                }

    return results


# ---------------------------------------------------------------------------
# PUBLIC API
# ---------------------------------------------------------------------------

def analyze_deep_ancestry(
    variants: Dict[str, Tuple[str, str, str]],
) -> Dict[str, Any]:
    """
    Analyze variants for deep Out-of-Africa ancestry and archaic admixture.

    Args:
        variants: Dict mapping rsID (lowercase) -> (chromosome, position, genotype)

    Returns:
        Dict with proportions, affinity_scores, used_snps, residual, key_markers,
        neanderthal_percentage, denisovan_percentage.
    """
    props, used_snps, resid = _estimate_proportions(variants)
    affinity_scores = _compute_affinity_scores(variants)
    key_markers = _summarize_key_markers(variants)
    neanderthal_pct = _estimate_neanderthal_percentage(variants)
    denisovan_pct = _estimate_denisovan_percentage(variants)

    return {
        "proportions": props,
        "affinity_scores": affinity_scores,
        "used_snps": used_snps,
        "residual": resid,
        "key_markers": key_markers,
        "neanderthal_percentage": neanderthal_pct,
        "denisovan_percentage": denisovan_pct,
    }


def generate_deep_ancestry_json(result: Dict[str, Any]) -> Dict[str, Any]:
    """Convert raw analysis result into structured JSON for the frontend."""
    props = result["proportions"]
    affinity_scores = result["affinity_scores"]
    used_snps = result["used_snps"]
    resid = result["residual"]
    key_markers = result["key_markers"]
    neanderthal_pct = result.get("neanderthal_percentage", 0.0)
    denisovan_pct = result.get("denisovan_percentage", 0.0)

    panel_size = len(AIMS_DEEP)
    snps_used = len(used_snps)

    # Sort by proportion (deep-time decomposition) descending
    sorted_pops = sorted(DEEP_POPS, key=lambda p: -props.get(p, 0.0))

    populations = []
    for pop in sorted_pops:
        proportion = round(props.get(pop, 0.0), 4)
        populations.append({
            "code": pop,
            "label": DEEP_POP_LABELS[pop],
            "affinityScore": affinity_scores.get(pop, 0.0),
            "proportion": proportion,
            "percentage": round(proportion * 100, 1),
            "description": DEEP_POP_DESCRIPTIONS[pop],
            "color": DEEP_POP_COLORS[pop],
        })

    top_pop = sorted_pops[0]

    return {
        "summary": {
            "panelSize": panel_size,
            "snpsUsed": snps_used,
            "coverage": round(snps_used / panel_size * 100, 1) if panel_size > 0 else 0.0,
            "neanderthalPercentage": neanderthal_pct,
            "denisovanPercentage": denisovan_pct,
        },
        "populations": populations,
        "topMatch": {
            "code": top_pop,
            "label": DEEP_POP_LABELS[top_pop],
            "description": DEEP_POP_DESCRIPTIONS[top_pop],
        },
        "keyMarkers": key_markers,
    }
