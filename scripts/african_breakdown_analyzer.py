"""
African Breakdown Analyzer
Estimates fine-grained African ancestry components using Ancestry-Informative
Markers (AIMs) calibrated to published continental and sub-Saharan genomics.

KEY REFERENCES:
  Tishkoff et al. 2009 (Science): "The Genetic Structure and History of Africans
    and African Americans"
    - 121 African populations, 4 African American populations, 60 non-African
    - 1,327 polymorphic markers; 14 ancestral clusters identified
    - Established that the highest genetic diversity on Earth resides in Africa

  Schlebusch et al. 2020 (Current Biology): "Genomic variation in seven Khoe-San
    groups reveals adaptation and complex African history"
    - Southern African Khoe-San branched from other modern humans ~200,000 years ago

  Prado-Martinez et al. 2013 (Nature): Great ape genome comparisons providing
    deep-time calibration of African divergence events

SUPPORTING REFERENCES:
  - Henn et al. 2011 (PNAS): Southern African Khoisan basal divergence
  - Pagani et al. 2015 (AJHG): East African (Horn of Africa) population structure
  - Pickrell et al. 2014 (Nat Comm): Ancient admixture in southern Africa
  - Patin et al. 2014 (Nat Comm): Pygmy rainforest hunter-gatherer ancestry
  - Busby et al. 2016 (eLife): Admixture in sub-Saharan Africa
  - Lazaridis et al. 2014 (Nature): Basal Eurasian backflow into East Africa
  - Fan et al. 2023 (Cell): Whole-genome sequencing of 180 individuals from 12
    indigenous African populations (updated Tishkoff framework)

SCIENTIFIC LIMITATIONS:
  Africa contains the deepest genetic splits in Homo sapiens. Our 50-AIM panel
  can distinguish broad components (West African Niger-Congo, East African
  Nilotic/Cushitic, Southern African Khoisan, Central African Pygmy, North
  African Berber) but cannot resolve finer ethnolinguistic groups at the
  individual level. Results show relative affinity scores, not direct descent.
  Malaria-adaptive variants (HbS, Duffy-null, G6PD) are informative BUT are
  under strong selection and must be interpreted alongside neutral markers.

Returns structured JSON for frontend visualization.
"""

from typing import Dict, Tuple, Any, Optional
import numpy as np
from scipy.optimize import nnls

# ---------------------------------------------------------------------------
# POPULATION DEFINITIONS
# ---------------------------------------------------------------------------

AFRICAN_POPS = [
    "West_African_Niger_Congo",
    "Central_African_Bantu",
    "East_African_Cushitic",
    "East_African_Nilotic",
    "Southern_African_Khoisan",
    "North_African_Berber",
    "Pygmy_Rainforest",
    "Horn_of_Africa_Semitic",
]

AFRICAN_POP_LABELS = {
    "West_African_Niger_Congo": "West African / Niger-Congo (Yoruba, Igbo, Akan, Fulani)",
    "Central_African_Bantu":    "Central African Bantu (Luba, Kongo, Kikuyu)",
    "East_African_Cushitic":    "East African Cushitic (Amhara, Oromo, Somali)",
    "East_African_Nilotic":     "East African Nilotic (Maasai, Nuer, Dinka)",
    "Southern_African_Khoisan": "Southern African Khoe-San (Ju|'hoan, Nama, Karretjie)",
    "North_African_Berber":     "North African Berber / Amazigh",
    "Pygmy_Rainforest":         "Central African Rainforest Pygmy (Mbuti, Biaka, Baka)",
    "Horn_of_Africa_Semitic":   "Horn of Africa Ethio-Semitic (Tigrinya, Amharic)",
}

# Warm earth palette: ochre, terracotta, sienna, umber, gold, olive
AFRICAN_POP_COLORS = {
    "West_African_Niger_Congo": "#B8450F",  # sienna red
    "Central_African_Bantu":    "#8B4513",  # saddle brown / umber
    "East_African_Cushitic":    "#D2691E",  # chocolate / terracotta
    "East_African_Nilotic":     "#CD853F",  # peru / tan
    "Southern_African_Khoisan": "#DAA520",  # goldenrod
    "North_African_Berber":     "#BDB76B",  # dark khaki / olive-gold
    "Pygmy_Rainforest":         "#556B2F",  # dark olive green
    "Horn_of_Africa_Semitic":   "#A0522D",  # sienna
}

AFRICAN_POP_DESCRIPTIONS = {
    "West_African_Niger_Congo": (
        "The West African Niger-Congo cluster (Yoruba, Igbo, Akan, Fulani, Mandinka) "
        "is the dominant ancestry component in African-American populations due to the "
        "transatlantic slave trade (1525-1866). West Africans have the highest frequency "
        "of the Duffy-null allele (rs2814778 CC approaches 100%) -- a malaria-resistance "
        "variant that blocks Plasmodium vivax invasion of red blood cells. They also "
        "carry the highest frequency of the sickle-cell allele (HbS, rs334 T), which "
        "in heterozygous form protects against P. falciparum malaria. The APOL1 G1 "
        "and G2 risk alleles (rs73885319, rs71785313), which confer kidney disease "
        "risk but also protect against trypanosomiasis, are essentially West African "
        "specific. Tishkoff 2009 found that Yoruba represent one of the most "
        "comprehensively sampled African reference populations in human genomics."
    ),
    "Central_African_Bantu": (
        "The Bantu expansion (~3,000-5,000 years ago) carried Niger-Congo languages "
        "and agriculture from the Cameroon-Nigeria border across equatorial and "
        "southern Africa. Central African Bantu speakers (Luba, Kongo, Kikuyu) "
        "share substantial ancestry with West Africans but also carry admixture "
        "from pre-Bantu rainforest hunter-gatherers (Pygmies). Tishkoff 2009 and "
        "Patin 2014 showed that modern Bantu-speaking farmers in Central Africa "
        "typically carry 5-25% Pygmy ancestry, reflecting 2,000+ years of contact. "
        "The Bantu expansion is one of the most significant demographic events of "
        "the Holocene and is why nearly 350 million people today speak a Bantu "
        "language across sub-Equatorial Africa."
    ),
    "East_African_Cushitic": (
        "Cushitic speakers (Afro-Asiatic family: Amhara, Oromo, Somali, Beja) occupy "
        "the Horn of Africa and show a distinctive genetic signature: ~40-60% "
        "sub-Saharan African ancestry combined with ~35-50% ancient West Eurasian "
        "ancestry, likely introduced from the Middle East via Bronze Age backflow "
        "(~3,000 BCE). This 'Basal Eurasian'-related backflow was first documented "
        "by Lazaridis 2014 and Pagani 2015. Cushitic populations independently "
        "evolved lactase persistence at the rs145946881 G variant (distinct from "
        "the European rs4988235 A variant), reflecting convergent evolution for "
        "dairy pastoralism. They have high frequencies of SLC24A5 (rs1426654) "
        "lighter-skin alleles introduced via the Eurasian backflow."
    ),
    "East_African_Nilotic": (
        "Nilotic speakers (Nilo-Saharan family: Maasai, Nuer, Dinka, Luo, Samburu) "
        "are tall, slender, pastoralist peoples of the East African savannas. They "
        "are genetically distinct from both Niger-Congo Bantu and Cushitic Afro-"
        "Asiatic speakers, carrying ancestry from an ancient East African component "
        "plus some Cushitic admixture. The Dinka of South Sudan are among the "
        "tallest populations on Earth. Nilotic populations independently evolved "
        "lactase persistence through the same rs145946881 variant as Cushitic "
        "pastoralists. They carry low frequencies of Duffy-null compared to West "
        "Africans but retain the HbS sickle-cell variant where malaria is endemic. "
        "Tishkoff 2007 identified convergent LP evolution in East African pastoralists."
    ),
    "Southern_African_Khoisan": (
        "The Khoe-San (Ju|'hoan, !Kung, Nama, Karretjie, Hadza-adjacent) represent "
        "the most genetically diverse human population on Earth. Schlebusch 2020 "
        "and Skoglund 2017 estimate that the Khoe-San lineage diverged from other "
        "modern humans as early as ~200,000-350,000 years ago -- pre-dating the "
        "split between all other African and non-African populations. Every other "
        "living human is more closely related to each other than any of us is to "
        "Khoe-San. They carry deeply divergent Y-chromosome (A00, A0, B) and "
        "mtDNA (L0) lineages. Morphologically they have distinctive features "
        "including reduced pigmentation relative to other Africans (lower HbS, "
        "lower Duffy-null because they historically inhabited non-malarial "
        "highlands). The Khoisan click languages are thought to preserve "
        "phonological features of unprecedented antiquity."
    ),
    "North_African_Berber": (
        "Berbers / Amazigh (Kabyle, Chaouis, Rif, Tuareg) are the indigenous people "
        "of North Africa (Morocco, Algeria, Tunisia, Libya) pre-dating the Arab "
        "conquest. Genetically they are a mix of ancient North African (with deep "
        "ties to Paleolithic Iberomaurusian ~15,000 BCE), Near Eastern Neolithic "
        "farmer ancestry, and sub-Saharan African admixture from trans-Saharan "
        "contact. Van de Loosdrecht 2018 (Science) showed that Late Pleistocene "
        "North Africans (Taforalt, Morocco) were a mix of Levantine (Natufian-like) "
        "and sub-Saharan components -- making Berbers one of the oldest continuous "
        "populations on Earth. They carry intermediate frequencies of European "
        "skin-pigmentation alleles (SLC24A5, SLC45A2) and typically lack the "
        "Duffy-null mutation characteristic of sub-Saharan Africans."
    ),
    "Pygmy_Rainforest": (
        "Central African rainforest foragers (Mbuti of Ituri, Biaka/Aka of "
        "Central African Republic, Baka of Cameroon) are among the most "
        "genetically and culturally distinctive peoples in Africa. Patin 2014 "
        "and Hsieh 2016 estimate that Western (Biaka/Baka) and Eastern (Mbuti) "
        "Pygmy lineages diverged from agricultural Bantu ancestors ~60,000-90,000 "
        "years ago and from each other ~20,000 years ago. Their short stature "
        "(average 1.4-1.5m) is a convergently evolved adaptation to the "
        "rainforest environment, driven by independent genetic architectures in "
        "Western vs Eastern Pygmies (Jarvis 2012). They carry unique MHC/HLA "
        "diversity and distinctive patterns of archaic introgression from now-"
        "extinct ghost African lineages (Hammer 2011; Lorente-Galdos 2019)."
    ),
    "Horn_of_Africa_Semitic": (
        "Ethio-Semitic speakers (Amharic, Tigrinya, Gurage, Tigre) occupy the "
        "Ethiopian and Eritrean highlands. Their linguistic ancestors arrived "
        "from Southern Arabia ~3,000 years ago, but genetically they are "
        "predominantly indigenous East African with ~40-50% West Eurasian "
        "admixture layered on top -- making them genetically similar to "
        "neighboring Cushitic speakers but with additional Arabian Peninsula "
        "influence. Pagani 2012 found that Ethiopian populations carry a "
        "distinctive combination of sub-Saharan, Near Eastern, and Arabian "
        "Peninsula ancestries. Ethiopian Highland Amhara and Tigrinya are "
        "well-known for high-altitude adaptation and remain a key population "
        "for understanding early human out-of-Africa migrations and later "
        "Eurasian backflow."
    ),
}

# ---------------------------------------------------------------------------
# AIM ALLELE FREQUENCIES FOR AFRICAN POPULATIONS
# Frequencies represent approximate means across published African reference
# populations (Tishkoff 2009, Fan 2023, 1000 Genomes YRI/LWK/MSL/ESN/GWD,
# HGDP Biaka/Mbuti/Mandenka/San/Yoruba, Pagani 2015, Schlebusch 2020).
#
# For each rsid: (ref, alt, weight, {pop: alt_freq})
# ---------------------------------------------------------------------------

AIMS_AFRICAN = {
    # --- CORE MALARIA / AFRICAN-DIAGNOSTIC MARKERS ---
    # rs2814778: Duffy-null (DARC/ACKR1 -67T>C). ~100% in West/Central Africans,
    # very low elsewhere. THE single most discriminating African AIM.
    "rs2814778": ("T", "C", 5.0, {
        "West_African_Niger_Congo": 0.99, "Central_African_Bantu": 0.95,
        "East_African_Cushitic": 0.50, "East_African_Nilotic": 0.80,
        "Southern_African_Khoisan": 0.20, "North_African_Berber": 0.15,
        "Pygmy_Rainforest": 0.97, "Horn_of_Africa_Semitic": 0.45,
    }),
    # rs334: HbS sickle cell (HBB Glu6Val). Heterozygote advantage against
    # P. falciparum; tracks malaria endemicity in Africa.
    "rs334": ("A", "T", 4.5, {
        "West_African_Niger_Congo": 0.15, "Central_African_Bantu": 0.12,
        "East_African_Cushitic": 0.05, "East_African_Nilotic": 0.08,
        "Southern_African_Khoisan": 0.00, "North_African_Berber": 0.02,
        "Pygmy_Rainforest": 0.08, "Horn_of_Africa_Semitic": 0.03,
    }),
    # rs73885319: APOL1 G1 allele (S342G). Confers CKD risk + trypanosomiasis
    # resistance. Nearly West African specific.
    "rs73885319": ("A", "G", 4.5, {
        "West_African_Niger_Congo": 0.22, "Central_African_Bantu": 0.10,
        "East_African_Cushitic": 0.01, "East_African_Nilotic": 0.02,
        "Southern_African_Khoisan": 0.00, "North_African_Berber": 0.01,
        "Pygmy_Rainforest": 0.04, "Horn_of_Africa_Semitic": 0.01,
    }),
    # rs71785313: APOL1 G2 allele (N388del/Y389del 6bp indel). West African.
    "rs71785313": ("I", "D", 4.0, {
        "West_African_Niger_Congo": 0.13, "Central_African_Bantu": 0.06,
        "East_African_Cushitic": 0.01, "East_African_Nilotic": 0.02,
        "Southern_African_Khoisan": 0.00, "North_African_Berber": 0.00,
        "Pygmy_Rainforest": 0.03, "Horn_of_Africa_Semitic": 0.01,
    }),
    # --- LACTASE PERSISTENCE (AFRICAN-SPECIFIC VARIANT) ---
    # rs145946881: C-14010G African-specific LP variant (Tishkoff 2007).
    # Dominant in East African pastoralists; absent in West/Central Africa.
    "rs145946881": ("C", "G", 4.0, {
        "West_African_Niger_Congo": 0.00, "Central_African_Bantu": 0.01,
        "East_African_Cushitic": 0.38, "East_African_Nilotic": 0.45,
        "Southern_African_Khoisan": 0.02, "North_African_Berber": 0.05,
        "Pygmy_Rainforest": 0.00, "Horn_of_Africa_Semitic": 0.30,
    }),
    # rs4988235: European LP variant (-13910 T). Present at non-trivial
    # frequency in North African Berbers and Horn of Africa (Eurasian backflow).
    "rs4988235": ("G", "A", 3.0, {
        "West_African_Niger_Congo": 0.00, "Central_African_Bantu": 0.00,
        "East_African_Cushitic": 0.08, "East_African_Nilotic": 0.02,
        "Southern_African_Khoisan": 0.00, "North_African_Berber": 0.18,
        "Pygmy_Rainforest": 0.00, "Horn_of_Africa_Semitic": 0.10,
    }),
    # --- SKIN PIGMENTATION ---
    # rs1426654: SLC24A5 A111T (derived G light-skin allele). ~99% in Europeans.
    # Gradient in Africa: low in sub-Saharan, high in Berbers/Ethio-Semitic.
    "rs1426654": ("A", "G", 3.5, {
        "West_African_Niger_Congo": 0.02, "Central_African_Bantu": 0.03,
        "East_African_Cushitic": 0.48, "East_African_Nilotic": 0.08,
        "Southern_African_Khoisan": 0.05, "North_African_Berber": 0.82,
        "Pygmy_Rainforest": 0.02, "Horn_of_Africa_Semitic": 0.45,
    }),
    # rs16891982: SLC45A2 F374L (derived G light-skin allele). Very low in
    # sub-Saharan Africans; intermediate in Berbers.
    "rs16891982": ("C", "G", 3.0, {
        "West_African_Niger_Congo": 0.01, "Central_African_Bantu": 0.01,
        "East_African_Cushitic": 0.05, "East_African_Nilotic": 0.02,
        "Southern_African_Khoisan": 0.02, "North_African_Berber": 0.55,
        "Pygmy_Rainforest": 0.01, "Horn_of_Africa_Semitic": 0.05,
    }),
    # rs1042602: TYR S192Y. Low frequency in sub-Saharan Africans.
    "rs1042602": ("C", "A", 2.5, {
        "West_African_Niger_Congo": 0.01, "Central_African_Bantu": 0.01,
        "East_African_Cushitic": 0.10, "East_African_Nilotic": 0.02,
        "Southern_African_Khoisan": 0.01, "North_African_Berber": 0.32,
        "Pygmy_Rainforest": 0.01, "Horn_of_Africa_Semitic": 0.08,
    }),
    # rs1800404: OCA2. Pigmentation variant with sub-Saharan patterning.
    "rs1800404": ("C", "T", 2.0, {
        "West_African_Niger_Congo": 0.85, "Central_African_Bantu": 0.83,
        "East_African_Cushitic": 0.60, "East_African_Nilotic": 0.78,
        "Southern_African_Khoisan": 0.62, "North_African_Berber": 0.50,
        "Pygmy_Rainforest": 0.88, "Horn_of_Africa_Semitic": 0.58,
    }),
    # rs1448484: Diagnostic skin pigmentation SNP in African populations.
    "rs1448484": ("C", "T", 2.0, {
        "West_African_Niger_Congo": 0.88, "Central_African_Bantu": 0.85,
        "East_African_Cushitic": 0.50, "East_African_Nilotic": 0.72,
        "Southern_African_Khoisan": 0.55, "North_African_Berber": 0.30,
        "Pygmy_Rainforest": 0.90, "Horn_of_Africa_Semitic": 0.48,
    }),
    # rs2733832: TYRP1 skin color variant.
    "rs2733832": ("T", "C", 2.0, {
        "West_African_Niger_Congo": 0.85, "Central_African_Bantu": 0.82,
        "East_African_Cushitic": 0.55, "East_African_Nilotic": 0.72,
        "Southern_African_Khoisan": 0.60, "North_African_Berber": 0.32,
        "Pygmy_Rainforest": 0.88, "Horn_of_Africa_Semitic": 0.50,
    }),
    # --- HAIR / EDAR / EVS ---
    # rs3827760: EDAR V370A. ~0% in Africans, distinguishes from East Asians/NE.
    "rs3827760": ("A", "G", 2.0, {
        "West_African_Niger_Congo": 0.00, "Central_African_Bantu": 0.00,
        "East_African_Cushitic": 0.00, "East_African_Nilotic": 0.00,
        "Southern_African_Khoisan": 0.00, "North_African_Berber": 0.00,
        "Pygmy_Rainforest": 0.00, "Horn_of_Africa_Semitic": 0.00,
    }),
    # --- G6PD DEFICIENCY MALARIA ADAPTIVE ---
    # rs1050828: G6PD A- (V68M) African variant - common in sub-Saharan males.
    "rs1050828": ("C", "T", 3.5, {
        "West_African_Niger_Congo": 0.15, "Central_African_Bantu": 0.14,
        "East_African_Cushitic": 0.05, "East_African_Nilotic": 0.12,
        "Southern_African_Khoisan": 0.01, "North_African_Berber": 0.02,
        "Pygmy_Rainforest": 0.10, "Horn_of_Africa_Semitic": 0.04,
    }),
    # rs1050829: G6PD N126D (A allele marker - African variant).
    "rs1050829": ("T", "C", 3.0, {
        "West_African_Niger_Congo": 0.30, "Central_African_Bantu": 0.28,
        "East_African_Cushitic": 0.08, "East_African_Nilotic": 0.20,
        "Southern_African_Khoisan": 0.02, "North_African_Berber": 0.03,
        "Pygmy_Rainforest": 0.22, "Horn_of_Africa_Semitic": 0.06,
    }),
    # --- HEMOGLOBIN / THALASSEMIA ---
    # rs33930165: HbC (beta6 Glu->Lys). West African, esp. Burkina Faso.
    "rs33930165": ("C", "T", 3.0, {
        "West_African_Niger_Congo": 0.08, "Central_African_Bantu": 0.02,
        "East_African_Cushitic": 0.01, "East_African_Nilotic": 0.01,
        "Southern_African_Khoisan": 0.00, "North_African_Berber": 0.00,
        "Pygmy_Rainforest": 0.01, "Horn_of_Africa_Semitic": 0.00,
    }),
    # --- IMMUNE / HLA REGION ---
    # rs9264942: HLA-C HIV control variant. Gradient across African populations.
    "rs9264942": ("T", "C", 1.5, {
        "West_African_Niger_Congo": 0.08, "Central_African_Bantu": 0.10,
        "East_African_Cushitic": 0.25, "East_African_Nilotic": 0.15,
        "Southern_African_Khoisan": 0.20, "North_African_Berber": 0.30,
        "Pygmy_Rainforest": 0.08, "Horn_of_Africa_Semitic": 0.28,
    }),
    # --- BITTER TASTE / TAS2R ---
    # rs1726866: TAS2R38 V296A. African sub-population variation.
    "rs1726866": ("G", "A", 1.0, {
        "West_African_Niger_Congo": 0.50, "Central_African_Bantu": 0.48,
        "East_African_Cushitic": 0.42, "East_African_Nilotic": 0.45,
        "Southern_African_Khoisan": 0.40, "North_African_Berber": 0.45,
        "Pygmy_Rainforest": 0.50, "Horn_of_Africa_Semitic": 0.44,
    }),
    # rs713598: TAS2R38 P49A.
    "rs713598": ("G", "C", 1.0, {
        "West_African_Niger_Congo": 0.55, "Central_African_Bantu": 0.52,
        "East_African_Cushitic": 0.48, "East_African_Nilotic": 0.50,
        "Southern_African_Khoisan": 0.45, "North_African_Berber": 0.48,
        "Pygmy_Rainforest": 0.55, "Horn_of_Africa_Semitic": 0.48,
    }),
    # --- EYE COLOR ---
    # rs12913832: HERC2 blue-eye. Essentially zero in sub-Saharan Africans,
    # nonzero in North African Berbers and Horn of Africa.
    "rs12913832": ("A", "G", 2.5, {
        "West_African_Niger_Congo": 0.00, "Central_African_Bantu": 0.00,
        "East_African_Cushitic": 0.05, "East_African_Nilotic": 0.00,
        "Southern_African_Khoisan": 0.00, "North_African_Berber": 0.15,
        "Pygmy_Rainforest": 0.00, "Horn_of_Africa_Semitic": 0.04,
    }),
    # --- ADH METABOLISM ---
    # rs1229984: ADH1B His48Arg. Asian-specific mostly; very low in Africa.
    "rs1229984": ("G", "A", 1.5, {
        "West_African_Niger_Congo": 0.00, "Central_African_Bantu": 0.01,
        "East_African_Cushitic": 0.02, "East_African_Nilotic": 0.00,
        "Southern_African_Khoisan": 0.01, "North_African_Berber": 0.03,
        "Pygmy_Rainforest": 0.00, "Horn_of_Africa_Semitic": 0.02,
    }),
    # --- EARWAX / APOCRINE ---
    # rs17822931: ABCC11 G->A wet/dry earwax. Nearly 100% wet in Africa.
    "rs17822931": ("C", "T", 1.5, {
        "West_African_Niger_Congo": 0.00, "Central_African_Bantu": 0.00,
        "East_African_Cushitic": 0.00, "East_African_Nilotic": 0.00,
        "Southern_African_Khoisan": 0.00, "North_African_Berber": 0.01,
        "Pygmy_Rainforest": 0.00, "Horn_of_Africa_Semitic": 0.00,
    }),
    # --- APOE ---
    # rs429358: APOE C/T (epsilon 4 tag). African ancestral freq is higher.
    "rs429358": ("T", "C", 1.5, {
        "West_African_Niger_Congo": 0.29, "Central_African_Bantu": 0.28,
        "East_African_Cushitic": 0.18, "East_African_Nilotic": 0.25,
        "Southern_African_Khoisan": 0.20, "North_African_Berber": 0.12,
        "Pygmy_Rainforest": 0.30, "Horn_of_Africa_Semitic": 0.16,
    }),
    # --- Y-CHROMOSOME PROXIES (E1b1a family vs A/B basal) ---
    # rs9786283: marker tagging Y E-M2 (E1b1a) haplogroup family (Bantu expansion).
    "rs9786283": ("C", "T", 1.5, {
        "West_African_Niger_Congo": 0.70, "Central_African_Bantu": 0.72,
        "East_African_Cushitic": 0.15, "East_African_Nilotic": 0.35,
        "Southern_African_Khoisan": 0.08, "North_African_Berber": 0.12,
        "Pygmy_Rainforest": 0.45, "Horn_of_Africa_Semitic": 0.20,
    }),
    # --- CHROMOSOME 7 (SLC6A4) ---
    "rs25531": ("A", "G", 1.0, {
        "West_African_Niger_Congo": 0.30, "Central_African_Bantu": 0.28,
        "East_African_Cushitic": 0.25, "East_African_Nilotic": 0.28,
        "Southern_African_Khoisan": 0.20, "North_African_Berber": 0.10,
        "Pygmy_Rainforest": 0.30, "Horn_of_Africa_Semitic": 0.22,
    }),
    # --- KHOISAN-SPECIFIC DEEPLY DIVERGED PROXIES ---
    # rs4540055: proxy for Khoisan-enriched variation near VKORC1 region.
    "rs4540055": ("T", "C", 2.5, {
        "West_African_Niger_Congo": 0.05, "Central_African_Bantu": 0.08,
        "East_African_Cushitic": 0.12, "East_African_Nilotic": 0.10,
        "Southern_African_Khoisan": 0.45, "North_African_Berber": 0.15,
        "Pygmy_Rainforest": 0.08, "Horn_of_Africa_Semitic": 0.13,
    }),
    # rs2250072: Khoisan-enriched SNP (proxy from Schlebusch 2020).
    "rs2250072": ("G", "A", 2.0, {
        "West_African_Niger_Congo": 0.10, "Central_African_Bantu": 0.12,
        "East_African_Cushitic": 0.15, "East_African_Nilotic": 0.14,
        "Southern_African_Khoisan": 0.55, "North_African_Berber": 0.20,
        "Pygmy_Rainforest": 0.13, "Horn_of_Africa_Semitic": 0.16,
    }),
    # --- PYGMY RAINFOREST HEIGHT / STATURE ---
    # rs2871865: IGF1R region (short-stature associated in Pygmies).
    "rs2871865": ("A", "G", 2.0, {
        "West_African_Niger_Congo": 0.22, "Central_African_Bantu": 0.25,
        "East_African_Cushitic": 0.20, "East_African_Nilotic": 0.18,
        "Southern_African_Khoisan": 0.22, "North_African_Berber": 0.20,
        "Pygmy_Rainforest": 0.58, "Horn_of_Africa_Semitic": 0.20,
    }),
    # rs314280: LIN28B region (Pygmy-enriched height locus).
    "rs314280": ("T", "C", 1.5, {
        "West_African_Niger_Congo": 0.30, "Central_African_Bantu": 0.32,
        "East_African_Cushitic": 0.28, "East_African_Nilotic": 0.28,
        "Southern_African_Khoisan": 0.30, "North_African_Berber": 0.28,
        "Pygmy_Rainforest": 0.60, "Horn_of_Africa_Semitic": 0.28,
    }),
    # --- HFE IRON OVERLOAD (LOW IN AFRICANS) ---
    "rs1800562": ("G", "A", 1.5, {
        "West_African_Niger_Congo": 0.00, "Central_African_Bantu": 0.00,
        "East_African_Cushitic": 0.01, "East_African_Nilotic": 0.00,
        "Southern_African_Khoisan": 0.00, "North_African_Berber": 0.02,
        "Pygmy_Rainforest": 0.00, "Horn_of_Africa_Semitic": 0.01,
    }),
    # rs1799945: HFE H63D (more cosmopolitan than C282Y).
    "rs1799945": ("C", "G", 1.0, {
        "West_African_Niger_Congo": 0.03, "Central_African_Bantu": 0.04,
        "East_African_Cushitic": 0.06, "East_African_Nilotic": 0.03,
        "Southern_African_Khoisan": 0.03, "North_African_Berber": 0.12,
        "Pygmy_Rainforest": 0.02, "Horn_of_Africa_Semitic": 0.07,
    }),
    # --- MALARIA: ATP2B4 ---
    # rs10900585: ATP2B4 malaria resistance variant.
    "rs10900585": ("T", "G", 2.0, {
        "West_African_Niger_Congo": 0.78, "Central_African_Bantu": 0.75,
        "East_African_Cushitic": 0.55, "East_African_Nilotic": 0.70,
        "Southern_African_Khoisan": 0.40, "North_African_Berber": 0.45,
        "Pygmy_Rainforest": 0.72, "Horn_of_Africa_Semitic": 0.52,
    }),
    # --- MCM6 / LACTASE REGULATORY REGION EAST AFRICAN G-13915 ---
    "rs41380347": ("T", "G", 2.5, {
        "West_African_Niger_Congo": 0.00, "Central_African_Bantu": 0.00,
        "East_African_Cushitic": 0.12, "East_African_Nilotic": 0.08,
        "Southern_African_Khoisan": 0.00, "North_African_Berber": 0.04,
        "Pygmy_Rainforest": 0.00, "Horn_of_Africa_Semitic": 0.10,
    }),
    # --- CCR5 delta32 (ABSENT in Africans, Eurasian marker) ---
    "rs333": ("I", "D", 2.0, {
        "West_African_Niger_Congo": 0.00, "Central_African_Bantu": 0.00,
        "East_African_Cushitic": 0.01, "East_African_Nilotic": 0.00,
        "Southern_African_Khoisan": 0.00, "North_African_Berber": 0.05,
        "Pygmy_Rainforest": 0.00, "Horn_of_Africa_Semitic": 0.01,
    }),
    # --- Y-CHROMOSOME A/B-BASAL MARKER ---
    # rs2032652: tags basal haplogroup A/B (oldest Y lineages, Khoisan/Pygmy).
    "rs2032652": ("G", "A", 1.5, {
        "West_African_Niger_Congo": 0.05, "Central_African_Bantu": 0.08,
        "East_African_Cushitic": 0.15, "East_African_Nilotic": 0.25,
        "Southern_African_Khoisan": 0.60, "North_African_Berber": 0.04,
        "Pygmy_Rainforest": 0.45, "Horn_of_Africa_Semitic": 0.12,
    }),
    # --- APOBEC3B / TRIM5 (HIV restriction) ---
    "rs11465804": ("T", "G", 1.0, {
        "West_African_Niger_Congo": 0.05, "Central_African_Bantu": 0.06,
        "East_African_Cushitic": 0.10, "East_African_Nilotic": 0.08,
        "Southern_African_Khoisan": 0.12, "North_African_Berber": 0.14,
        "Pygmy_Rainforest": 0.05, "Horn_of_Africa_Semitic": 0.10,
    }),
    # --- FY GENE REGION (DUFFY CONTEXT SNP) ---
    "rs12075": ("A", "G", 2.0, {
        "West_African_Niger_Congo": 0.92, "Central_African_Bantu": 0.90,
        "East_African_Cushitic": 0.60, "East_African_Nilotic": 0.80,
        "Southern_African_Khoisan": 0.50, "North_African_Berber": 0.42,
        "Pygmy_Rainforest": 0.88, "Horn_of_Africa_Semitic": 0.55,
    }),
    # --- IL10 INFLAMMATORY ---
    "rs1800872": ("C", "A", 1.0, {
        "West_African_Niger_Congo": 0.22, "Central_African_Bantu": 0.20,
        "East_African_Cushitic": 0.18, "East_African_Nilotic": 0.20,
        "Southern_African_Khoisan": 0.25, "North_African_Berber": 0.22,
        "Pygmy_Rainforest": 0.24, "Horn_of_Africa_Semitic": 0.18,
    }),
    # --- VDR VITAMIN D (relevant in low UV N. Europeans, high UV Africans) ---
    "rs2228570": ("G", "A", 1.0, {
        "West_African_Niger_Congo": 0.65, "Central_African_Bantu": 0.62,
        "East_African_Cushitic": 0.50, "East_African_Nilotic": 0.58,
        "Southern_African_Khoisan": 0.55, "North_African_Berber": 0.38,
        "Pygmy_Rainforest": 0.68, "Horn_of_Africa_Semitic": 0.48,
    }),
    # --- TLR1 Population-specific variant ---
    "rs4833095": ("T", "C", 1.0, {
        "West_African_Niger_Congo": 0.25, "Central_African_Bantu": 0.28,
        "East_African_Cushitic": 0.30, "East_African_Nilotic": 0.28,
        "Southern_African_Khoisan": 0.35, "North_African_Berber": 0.42,
        "Pygmy_Rainforest": 0.28, "Horn_of_Africa_Semitic": 0.32,
    }),
    # --- CYP3A5*1 (high expresser allele frequent in Africans) ---
    "rs776746": ("C", "T", 1.5, {
        "West_African_Niger_Congo": 0.80, "Central_African_Bantu": 0.78,
        "East_African_Cushitic": 0.70, "East_African_Nilotic": 0.75,
        "Southern_African_Khoisan": 0.72, "North_African_Berber": 0.28,
        "Pygmy_Rainforest": 0.80, "Horn_of_Africa_Semitic": 0.65,
    }),
    # --- KRT77 hair keratin variant ---
    "rs11803731": ("T", "A", 1.0, {
        "West_African_Niger_Congo": 0.00, "Central_African_Bantu": 0.00,
        "East_African_Cushitic": 0.05, "East_African_Nilotic": 0.00,
        "Southern_African_Khoisan": 0.02, "North_African_Berber": 0.20,
        "Pygmy_Rainforest": 0.00, "Horn_of_Africa_Semitic": 0.04,
    }),
    # --- BHLHE40 hair morphology ---
    "rs2326171": ("C", "T", 1.0, {
        "West_African_Niger_Congo": 0.10, "Central_African_Bantu": 0.12,
        "East_African_Cushitic": 0.28, "East_African_Nilotic": 0.15,
        "Southern_African_Khoisan": 0.22, "North_African_Berber": 0.42,
        "Pygmy_Rainforest": 0.10, "Horn_of_Africa_Semitic": 0.25,
    }),
    # --- MHC HLA-DQA1 ---
    "rs2647046": ("A", "G", 1.0, {
        "West_African_Niger_Congo": 0.40, "Central_African_Bantu": 0.38,
        "East_African_Cushitic": 0.32, "East_African_Nilotic": 0.35,
        "Southern_African_Khoisan": 0.48, "North_African_Berber": 0.30,
        "Pygmy_Rainforest": 0.42, "Horn_of_Africa_Semitic": 0.33,
    }),
    # --- UGT1A1*28 (pigmentation/bilirubin) ---
    "rs8175347": ("T", "A", 1.0, {
        "West_African_Niger_Congo": 0.42, "Central_African_Bantu": 0.40,
        "East_African_Cushitic": 0.35, "East_African_Nilotic": 0.38,
        "Southern_African_Khoisan": 0.32, "North_African_Berber": 0.30,
        "Pygmy_Rainforest": 0.45, "Horn_of_Africa_Semitic": 0.33,
    }),
    # --- HP haptoglobin (malaria response) ---
    "rs2000999": ("G", "A", 1.0, {
        "West_African_Niger_Congo": 0.08, "Central_African_Bantu": 0.10,
        "East_African_Cushitic": 0.15, "East_African_Nilotic": 0.12,
        "Southern_African_Khoisan": 0.10, "North_African_Berber": 0.25,
        "Pygmy_Rainforest": 0.08, "Horn_of_Africa_Semitic": 0.18,
    }),
    # --- PARK2 selection ---
    "rs10945895": ("C", "T", 0.8, {
        "West_African_Niger_Congo": 0.55, "Central_African_Bantu": 0.52,
        "East_African_Cushitic": 0.48, "East_African_Nilotic": 0.50,
        "Southern_African_Khoisan": 0.45, "North_African_Berber": 0.42,
        "Pygmy_Rainforest": 0.55, "Horn_of_Africa_Semitic": 0.46,
    }),
    # --- AHR metabolism ---
    "rs2066853": ("G", "A", 0.8, {
        "West_African_Niger_Congo": 0.05, "Central_African_Bantu": 0.06,
        "East_African_Cushitic": 0.08, "East_African_Nilotic": 0.06,
        "Southern_African_Khoisan": 0.08, "North_African_Berber": 0.18,
        "Pygmy_Rainforest": 0.04, "Horn_of_Africa_Semitic": 0.10,
    }),
    # --- KITLG pigmentation locus (Berber-enriched) ---
    "rs12821256": ("T", "C", 1.0, {
        "West_African_Niger_Congo": 0.01, "Central_African_Bantu": 0.01,
        "East_African_Cushitic": 0.05, "East_African_Nilotic": 0.02,
        "Southern_African_Khoisan": 0.02, "North_African_Berber": 0.30,
        "Pygmy_Rainforest": 0.02, "Horn_of_Africa_Semitic": 0.05,
    }),
    # --- ALDH1A2 ---
    "rs4646404": ("G", "A", 0.8, {
        "West_African_Niger_Congo": 0.35, "Central_African_Bantu": 0.32,
        "East_African_Cushitic": 0.28, "East_African_Nilotic": 0.30,
        "Southern_African_Khoisan": 0.35, "North_African_Berber": 0.25,
        "Pygmy_Rainforest": 0.38, "Horn_of_Africa_Semitic": 0.28,
    }),
    # --- LARGE Lassa fever resistance ---
    "rs72646792": ("A", "G", 1.5, {
        "West_African_Niger_Congo": 0.18, "Central_African_Bantu": 0.08,
        "East_African_Cushitic": 0.03, "East_African_Nilotic": 0.05,
        "Southern_African_Khoisan": 0.02, "North_African_Berber": 0.02,
        "Pygmy_Rainforest": 0.05, "Horn_of_Africa_Semitic": 0.03,
    }),
    # --- MCM6 CC additional East African LP variant ---
    "rs41525747": ("C", "G", 1.5, {
        "West_African_Niger_Congo": 0.00, "Central_African_Bantu": 0.00,
        "East_African_Cushitic": 0.08, "East_African_Nilotic": 0.05,
        "Southern_African_Khoisan": 0.00, "North_African_Berber": 0.02,
        "Pygmy_Rainforest": 0.00, "Horn_of_Africa_Semitic": 0.06,
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
    Run weighted NNLS decomposition into African sub-population components.

    Returns (proportions_dict, used_snps_list, residual).
    """
    n_pops = len(AFRICAN_POPS)
    rows_A = []
    rows_f = []
    rows_w = []
    used_snps = []

    for rsid, (ref, alt, weight, pop_freqs) in AIMS_AFRICAN.items():
        variant_data = variants.get(rsid.lower())
        if variant_data is None:
            continue
        _chrom, _pos, genotype = variant_data
        dosage = _get_dosage(genotype, alt)
        if dosage is None:
            continue
        obs_freq = dosage / 2.0
        row = [pop_freqs.get(p, 0.0) for p in AFRICAN_POPS]
        rows_A.append(row)
        rows_f.append(obs_freq)
        rows_w.append(weight)
        used_snps.append(rsid)

    if len(used_snps) < 5:
        return {p: 1.0 / n_pops for p in AFRICAN_POPS}, used_snps, 0.0

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

    props = {p: float(x[i]) for i, p in enumerate(AFRICAN_POPS)}
    return props, used_snps, float(resid)


# ---------------------------------------------------------------------------
# AFFINITY SCORES (min-max normalized distances)
# ---------------------------------------------------------------------------

def _compute_affinity_scores(
    variants: Dict[str, Tuple[str, str, str]],
) -> Dict[str, float]:
    """
    Compute relative affinity to each African population using min-max
    normalization. Closest population scores 100, furthest scores 0.
    """
    raw_distances: Dict[str, float] = {}

    for pop in AFRICAN_POPS:
        weighted_sq_diff = 0.0
        w_sum = 0.0
        for rsid, (ref, alt, weight, pop_freqs) in AIMS_AFRICAN.items():
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
        return {p: 0.0 for p in AFRICAN_POPS}

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
# KEY AFRICAN MARKERS SUMMARY
# ---------------------------------------------------------------------------

def _summarize_key_markers(
    variants: Dict[str, Tuple[str, str, str]],
) -> Dict[str, Dict[str, str]]:
    """Report the key markers distinguishing African sub-populations."""
    results: Dict[str, Dict[str, str]] = {}

    # Duffy-null (rs2814778) - THE core West/Central African malaria adaptation
    variant_data = variants.get("rs2814778")
    if variant_data is not None:
        _c, _p, geno = variant_data
        if geno and geno not in ("--", "00", "NN"):
            c_count = geno.upper().count("C")
            if c_count == 2:
                results["duffy_null"] = {
                    "status": "CC",
                    "detail": "CC -- Duffy-null homozygote; strong West/Central African "
                              "signature. Confers resistance to Plasmodium vivax malaria.",
                }
            elif c_count == 1:
                results["duffy_null"] = {
                    "status": "TC",
                    "detail": "TC -- one Duffy-null allele; admixed African ancestry "
                              "or East African/North African.",
                }
            else:
                results["duffy_null"] = {
                    "status": "TT",
                    "detail": "TT -- Duffy positive; Eurasian or Southern African Khoisan pattern.",
                }

    # HbS sickle cell (rs334) - malaria heterozygote advantage
    variant_data = variants.get("rs334")
    if variant_data is not None:
        _c, _p, geno = variant_data
        if geno and geno not in ("--", "00", "NN"):
            t_count = geno.upper().count("T")
            if t_count == 2:
                results["sickle_cell"] = {
                    "status": "TT",
                    "detail": "TT -- sickle cell anemia homozygote (HbSS). Medical condition; "
                              "West African malaria-endemic ancestry.",
                }
            elif t_count == 1:
                results["sickle_cell"] = {
                    "status": "AT",
                    "detail": "AT -- sickle cell TRAIT carrier (HbAS); protective against "
                              "P. falciparum malaria. Strong West/Central African signal.",
                }
            else:
                results["sickle_cell"] = {
                    "status": "AA",
                    "detail": "AA -- no HbS; normal hemoglobin.",
                }

    # APOL1 G1 (rs73885319) - kidney risk + trypanosomiasis resistance
    variant_data = variants.get("rs73885319")
    if variant_data is not None:
        _c, _p, geno = variant_data
        if geno and geno not in ("--", "00", "NN"):
            g_count = geno.upper().count("G")
            if g_count > 0:
                status = "AG" if g_count == 1 else "GG"
                results["apol1_g1"] = {
                    "status": status,
                    "detail": f"{status} -- APOL1 G1 risk allele present. "
                              "West African specific. Confers kidney disease risk but "
                              "protects against African sleeping sickness.",
                }
            else:
                results["apol1_g1"] = {
                    "status": "AA",
                    "detail": "AA -- no APOL1 G1 risk allele.",
                }

    # East African lactase persistence (rs145946881) - independent LP evolution
    variant_data = variants.get("rs145946881")
    if variant_data is not None:
        _c, _p, geno = variant_data
        if geno and geno not in ("--", "00", "NN"):
            g_count = geno.upper().count("G")
            if g_count == 2:
                results["lactase_african"] = {
                    "status": "GG",
                    "detail": "GG -- homozygous for African-specific lactase persistence "
                              "(C-14010G). Strong Cushitic/Nilotic pastoralist signal.",
                }
            elif g_count == 1:
                results["lactase_african"] = {
                    "status": "CG",
                    "detail": "CG -- heterozygous for African LP variant. East African "
                              "pastoralist ancestry.",
                }
            else:
                results["lactase_african"] = {
                    "status": "CC",
                    "detail": "CC -- ancestral; no African LP variant.",
                }

    # SLC24A5 (rs1426654) - N African / Ethio-Semitic lighter skin
    variant_data = variants.get("rs1426654")
    if variant_data is not None:
        _c, _p, geno = variant_data
        if geno and geno not in ("--", "00", "NN"):
            g_count = geno.upper().count("G")
            if g_count == 2:
                results["slc24a5"] = {
                    "status": "GG",
                    "detail": "GG -- homozygous derived (light-skin). Very atypical for "
                              "sub-Saharan Africans; suggests Berber/Ethio-Semitic or "
                              "Eurasian admixture.",
                }
            elif g_count == 1:
                results["slc24a5"] = {
                    "status": "AG",
                    "detail": "AG -- heterozygous; North African Berber, Horn of Africa, "
                              "or admixed sub-Saharan pattern.",
                }
            else:
                results["slc24a5"] = {
                    "status": "AA",
                    "detail": "AA -- ancestral allele; typical sub-Saharan African pattern.",
                }

    # G6PD A- (rs1050828) - African malaria adaptive deficiency
    variant_data = variants.get("rs1050828")
    if variant_data is not None:
        _c, _p, geno = variant_data
        if geno and geno not in ("--", "00", "NN"):
            t_count = geno.upper().count("T")
            if t_count > 0:
                status = "CT" if t_count == 1 else "TT"
                results["g6pd_a_minus"] = {
                    "status": status,
                    "detail": f"{status} -- G6PD A- variant present. African-specific "
                              "deficiency allele conferring P. falciparum protection. "
                              "(X-linked -- males hemizygous.)",
                }

    return results


# ---------------------------------------------------------------------------
# PUBLIC API
# ---------------------------------------------------------------------------

def analyze_african_breakdown(
    variants: Dict[str, Tuple[str, str, str]],
) -> Dict[str, Any]:
    """
    Analyze variants for fine-grained African ancestry breakdown.

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


def generate_african_breakdown_json(result: Dict[str, Any]) -> Dict[str, Any]:
    """
    Convert raw analysis result into the structured JSON for the frontend.

    Args:
        result: Output from analyze_african_breakdown()

    Returns:
        Structured JSON dict matching the report schema.
    """
    props = result["proportions"]
    affinity_scores = result["affinity_scores"]
    used_snps = result["used_snps"]
    resid = result["residual"]
    key_markers = result["key_markers"]

    panel_size = len(AIMS_AFRICAN)
    snps_used = len(used_snps)

    # Build population list sorted by affinity score descending
    sorted_pops = sorted(AFRICAN_POPS, key=lambda p: -affinity_scores.get(p, 0.0))

    populations = []
    for pop in sorted_pops:
        proportion = round(props.get(pop, 0.0), 4)
        populations.append({
            "code": pop,
            "label": AFRICAN_POP_LABELS[pop],
            "affinityScore": affinity_scores.get(pop, 0.0),
            "proportion": proportion,
            "percentage": round(proportion * 100, 1),
            "description": AFRICAN_POP_DESCRIPTIONS[pop],
            "color": AFRICAN_POP_COLORS[pop],
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
            "label": AFRICAN_POP_LABELS[top_pop],
            "description": AFRICAN_POP_DESCRIPTIONS[top_pop],
        },
        "keyMarkers": key_markers,
    }
