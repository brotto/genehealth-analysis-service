"""
Native American Heritage Analyzer
Estimates genetic affinity to Indigenous peoples of the Americas using
Ancestry-Informative Markers (AIMs) with allele frequencies calibrated from
published Amerindian population genomics.

KEY REFERENCE:
  Reich et al. 2012 (Nature): "Reconstructing Native American population history"
    - Identified 3 streams of migration from Asia to the Americas:
        (1) First Americans (majority of ancestry in most Native populations)
        (2) Eskimo-Aleut (Arctic coast)
        (3) Na-Dene speakers (northwest North America)
    - Used 52 Native American populations and ~365,000 SNPs

SUPPORTING REFERENCES:
  - Moreno-Estrada et al. 2013 (PLoS Genetics): Mexican Indigenous population structure
  - Moreno-Estrada et al. 2014 (Science): Caribbean and Mesoamerican structure
  - HGDP: Human Genome Diversity Project (Cann et al. 2002; Bergstrom et al. 2020)
  - Skoglund et al. 2015 (Nature): Amazonian "Population Y" Australasian signal
  - Raghavan et al. 2014 (Nature): Upward Sun River genome, Beringian ancestry
  - Posth et al. 2018 (Cell): Paleogenomics of the Americas, South American structure

SCIENTIFIC LIMITATIONS:
  Native American populations share substantial ancestry derived from a single
  ancestral Beringian source, so individual sub-populations (e.g. Mesoamerican
  vs. Andean) differ only subtly at most loci. A 40-50 SNP AIM panel can
  identify broad Native American ancestry robustly, but sub-regional assignment
  within the Americas carries wider uncertainty. Results reflect relative
  affinity scores, not direct descent or tribal affiliation.

Returns structured JSON for frontend visualization.
"""

from typing import Dict, Tuple, Any, Optional
import numpy as np
from scipy.optimize import nnls

# ---------------------------------------------------------------------------
# POPULATION DEFINITIONS
# ---------------------------------------------------------------------------

NATIVE_POPS = [
    "North_American_Algonquian",
    "Mesoamerican",
    "Andean",
    "Amazonian",
    "Eskimo_Aleut",
    "Patagonian_Southern",
]

NATIVE_POP_LABELS = {
    "North_American_Algonquian": "North American (Algonquian / Iroquoian / Eastern Woodlands)",
    "Mesoamerican":              "Mesoamerican (Nahua / Maya / Zapotec)",
    "Andean":                    "Andean (Quechua / Aymara / Inca Descendants)",
    "Amazonian":                 "Amazonian (Tupi / Arawak / Amazon Basin)",
    "Eskimo_Aleut":              "Eskimo-Aleut (Inuit / Yupik / Aleut)",
    "Patagonian_Southern":       "Southern Cone (Mapuche / Selknam / Tehuelche)",
}

NATIVE_POP_COLORS = {
    "North_American_Algonquian": "#B7410E",  # rust red
    "Mesoamerican":              "#C89B3C",  # ochre
    "Andean":                    "#1F9B9B",  # turquoise
    "Amazonian":                 "#3E2A16",  # deep brown
    "Eskimo_Aleut":              "#6B8E9B",  # cold grey-blue
    "Patagonian_Southern":       "#D4B483",  # sand
}

NATIVE_POP_DESCRIPTIONS = {
    "North_American_Algonquian": (
        "Eastern Woodlands peoples (Algonquian, Iroquoian and related nations) "
        "descend from the First Americans migration (Reich et al. 2012), the "
        "founding stream that populated most of the Americas south of the Arctic. "
        "This cluster includes ancestors of the Ojibwe, Cree, Mi'kmaq, Lenape, "
        "Mohawk, Cherokee and dozens of other nations. Genetically, they carry "
        "high frequencies of Amerindian-specific mtDNA haplogroups A2, B2, C1, "
        "D1, and X2a -- the last of which is uniquely concentrated in northeastern "
        "North America. The EDAR rs3827760 'shovel-incisor' allele is essentially "
        "fixed, and Duffy-null (rs2814778-C) is absent. Lactase non-persistence "
        "is nearly universal, reflecting thousands of years without dairy farming."
    ),
    "Mesoamerican": (
        "Mesoamerican peoples include the Nahua (Aztec descendants), Maya, "
        "Zapotec, Mixtec, Purepecha and many other Indigenous Mexican and Central "
        "American nations. Moreno-Estrada et al. 2014 showed Mesoamerica has one "
        "of the richest Indigenous population structures in the Americas, with "
        "distinct Maya, Nahua, and Zapotec genetic clusters. The EDAR 'shovel "
        "tooth' allele is near fixation, and the iconic hairless-skin pigmentation "
        "profile reflects selection for UV tolerance at tropical latitudes. "
        "Mesoamerican populations contribute the largest Indigenous component "
        "of Mexican mestizo ancestry (~45-60% on average)."
    ),
    "Andean": (
        "Andean peoples -- Quechua, Aymara, and descendants of the Inca, Wari, "
        "Tiwanaku and Chavin civilizations -- inhabit the high-altitude spine of "
        "South America from Ecuador to northern Chile and Argentina. They show "
        "strong high-altitude adaptation signatures, including selection on "
        "EGLN1, PRKAA1 and other hypoxia-response genes distinct from the "
        "Tibetan EPAS1 pathway. Posth et al. 2018 identified a deep North-South "
        "split among early South Americans, with Andean populations retaining "
        "the clearest 'Southern Native American' ancestry. Their lactase "
        "non-persistence rate is essentially 100%, and they show the highest "
        "frequency of certain Native-specific Y-chromosome haplogroup Q lineages."
    ),
    "Amazonian": (
        "Amazonian peoples include Tupi, Arawak, Carib, Ge and many other "
        "language families spanning the Amazon Basin from the Andes to the "
        "Atlantic. Skoglund et al. 2015 famously detected an unexpected "
        "Australasian-related genetic signal (nicknamed 'Population Y') in "
        "some Amazonian groups (Surui, Karitiana, Xavante), suggesting an "
        "additional ancient migration or structured source population. "
        "Amazonian populations carry high frequencies of mtDNA haplogroups "
        "C1 and D1, and show strong signatures of adaptation to tropical "
        "pathogens, including pathogen-response HLA alleles. The ABCC11 "
        "wet-earwax allele is almost absent."
    ),
    "Eskimo_Aleut": (
        "Eskimo-Aleut peoples (Inuit, Yupik, Inupiat, Aleut) represent the "
        "second of three migration streams into the Americas identified by "
        "Reich et al. 2012, arriving ~5,000-6,000 years ago across the "
        "Bering Strait. They carry distinctive mtDNA haplogroups A2a, A2b "
        "and D2a, and show strong cold-adaptation signals in CPT1A (high "
        "frequency of the rare P479L variant unique to Arctic peoples), "
        "FADS genes (fatty-acid metabolism adapted to marine-mammal diet), "
        "and TBX15/WARS2 (body-fat distribution). They are also the only "
        "Native American group with substantial ABCC11 'dry earwax' frequencies, "
        "reflecting their distinct Northeast Asian affinity."
    ),
    "Patagonian_Southern": (
        "Southern Cone peoples -- Mapuche, Tehuelche, Selk'nam (Ona), Yaghan "
        "and Kaweskar -- inhabited the southernmost Americas, from central "
        "Chile and Argentina down to Tierra del Fuego. Ancient DNA from "
        "Patagonian hunter-gatherers (Posth et al. 2018; de la Fuente et al. "
        "2018) shows they descend from the 'Southern Native American' branch "
        "that split from other South American populations ~14,000 years ago. "
        "Their genomes carry signatures of adaptation to cold, wind, and "
        "marine foraging. The Selk'nam were tragically driven to extinction "
        "in the late 19th century, but Mapuche and other groups remain "
        "genetically and culturally vibrant across southern Chile and Argentina."
    ),
}

# ---------------------------------------------------------------------------
# AIM ALLELE FREQUENCIES FOR NATIVE AMERICAN POPULATIONS
# Frequencies based on HGDP Native samples (Karitiana, Surui, Pima, Maya,
# Colombian), 1000 Genomes admixture references, and published Amerindian
# population-genetics literature (Reich 2012; Moreno-Estrada 2014;
# Raghavan 2014; Posth 2018; Bergstrom 2020 HGDP).
#
# For each rsid: (ref, alt, weight, {pop: alt_freq})
# ---------------------------------------------------------------------------

AIMS_NATIVE = {
    # --- High-weight Native-American diagnostic markers --------------------

    # EDAR V370A (rs3827760) -- "shovel teeth", thick hair, sweat glands.
    # Nearly fixed (~0.95-1.00) in Native Americans / East Asians. Near zero
    # in West Eurasians and Africans. Reference: Kamberov et al. 2013.
    "rs3827760": ("A", "G", 5.0, {
        "North_American_Algonquian": 0.96, "Mesoamerican": 0.97, "Andean": 0.98,
        "Amazonian": 0.99, "Eskimo_Aleut": 0.95, "Patagonian_Southern": 0.97,
    }),

    # DARC/FY rs2814778 -- Duffy-null. FIXED ancestral (T) in Native Americans;
    # the derived (C) allele that hides Duffy is ~1.0 in West Africans.
    "rs2814778": ("T", "C", 5.0, {
        "North_American_Algonquian": 0.00, "Mesoamerican": 0.00, "Andean": 0.00,
        "Amazonian": 0.00, "Eskimo_Aleut": 0.00, "Patagonian_Southern": 0.00,
    }),

    # ABCC11 rs17822931 -- C=wet earwax, T=dry earwax. Most Native Americans
    # carry ~10-30% dry-earwax allele (low vs. East Asians) except Eskimo-Aleut
    # who show higher frequencies (~60-75%) reflecting Northeast Asian affinity.
    "rs17822931": ("C", "T", 4.5, {
        "North_American_Algonquian": 0.20, "Mesoamerican": 0.15, "Andean": 0.10,
        "Amazonian": 0.10, "Eskimo_Aleut": 0.70, "Patagonian_Southern": 0.12,
    }),

    # LCT -13910 rs4988235 -- lactase persistence. Essentially ABSENT in
    # pre-contact Native populations (no dairy cattle). Very low frequencies.
    "rs4988235": ("A", "G", 4.0, {
        "North_American_Algonquian": 0.01, "Mesoamerican": 0.01, "Andean": 0.00,
        "Amazonian": 0.00, "Eskimo_Aleut": 0.02, "Patagonian_Southern": 0.00,
    }),

    # SLC24A5 rs1426654 -- A=light skin European allele. Near-zero in
    # un-admixed Native Americans (Norton et al. 2007; Canfield et al. 2013).
    "rs1426654": ("A", "G", 4.0, {
        "North_American_Algonquian": 0.05, "Mesoamerican": 0.04, "Andean": 0.03,
        "Amazonian": 0.01, "Eskimo_Aleut": 0.06, "Patagonian_Southern": 0.04,
    }),

    # SLC45A2 rs16891982 -- G=European light-skin variant. Very low in
    # un-admixed Native populations.
    "rs16891982": ("C", "G", 4.0, {
        "North_American_Algonquian": 0.02, "Mesoamerican": 0.02, "Andean": 0.01,
        "Amazonian": 0.01, "Eskimo_Aleut": 0.02, "Patagonian_Southern": 0.02,
    }),

    # HERC2 rs12913832 -- G=blue eyes. Virtually absent in un-admixed Native
    # Americans (eye-color gradient is ~zero blue-eye allele frequency).
    "rs12913832": ("A", "G", 3.5, {
        "North_American_Algonquian": 0.02, "Mesoamerican": 0.02, "Andean": 0.01,
        "Amazonian": 0.01, "Eskimo_Aleut": 0.03, "Patagonian_Southern": 0.02,
    }),

    # MC1R rs1805007 (R151C) -- European red-hair allele. Absent in Native Americans.
    "rs1805007": ("C", "T", 2.5, {
        "North_American_Algonquian": 0.00, "Mesoamerican": 0.00, "Andean": 0.00,
        "Amazonian": 0.00, "Eskimo_Aleut": 0.00, "Patagonian_Southern": 0.00,
    }),

    # OCA2 rs1800414 (H615R) -- "East Asian/Native American light-skin" allele.
    # HIGH in East Asians/Native Americans (~0.4-0.6), near zero in Europeans.
    "rs1800414": ("T", "C", 3.5, {
        "North_American_Algonquian": 0.45, "Mesoamerican": 0.50, "Andean": 0.55,
        "Amazonian": 0.55, "Eskimo_Aleut": 0.35, "Patagonian_Southern": 0.50,
    }),

    # CPT1A P479L rs80356779 -- Arctic-specific fatty acid oxidation variant.
    # FIXED near 1.0 in Inuit/Yupik, essentially zero in all other populations.
    "rs80356779": ("G", "A", 4.5, {
        "North_American_Algonquian": 0.01, "Mesoamerican": 0.00, "Andean": 0.00,
        "Amazonian": 0.00, "Eskimo_Aleut": 0.85, "Patagonian_Southern": 0.00,
    }),

    # EGLN1 rs186996510 -- Andean high-altitude adaptation (Bigham et al. 2010).
    # Near-fixed derived allele in Andean populations (~0.85-0.95); low elsewhere.
    # Note: some sources use proxy SNPs; rs480902 is a common Andean-enriched tag.
    "rs480902": ("T", "C", 3.0, {
        "North_American_Algonquian": 0.15, "Mesoamerican": 0.20, "Andean": 0.78,
        "Amazonian": 0.20, "Eskimo_Aleut": 0.15, "Patagonian_Southern": 0.30,
    }),

    # FADS1 rs174546 -- fatty-acid desaturase. Higher derived frequency in
    # Inuit (marine-diet adaptation, Fumagalli et al. 2015).
    "rs174546": ("C", "T", 3.0, {
        "North_American_Algonquian": 0.40, "Mesoamerican": 0.40, "Andean": 0.45,
        "Amazonian": 0.50, "Eskimo_Aleut": 0.95, "Patagonian_Southern": 0.45,
    }),

    # TBX15 rs2298080 -- Inuit body-fat distribution (Racimo et al. 2017).
    "rs2298080": ("G", "A", 2.5, {
        "North_American_Algonquian": 0.15, "Mesoamerican": 0.15, "Andean": 0.15,
        "Amazonian": 0.15, "Eskimo_Aleut": 0.60, "Patagonian_Southern": 0.15,
    }),

    # --- Amerindian-elevated metabolic / immune markers --------------------

    # ABCA1 Arg230Cys rs9282541 -- Native-American-specific variant ("Amerindian
    # allele"), strongly associated with low HDL cholesterol. Up to ~25-30%
    # in Mesoamerican populations, very low elsewhere in the world.
    "rs9282541": ("C", "T", 4.0, {
        "North_American_Algonquian": 0.12, "Mesoamerican": 0.28, "Andean": 0.18,
        "Amazonian": 0.15, "Eskimo_Aleut": 0.03, "Patagonian_Southern": 0.15,
    }),

    # HNF1A rs1169288 -- MODY3 / "SIGMA" diabetes risk variant; common
    # at moderate frequencies in Native American populations.
    "rs1169288": ("C", "A", 2.0, {
        "North_American_Algonquian": 0.30, "Mesoamerican": 0.35, "Andean": 0.30,
        "Amazonian": 0.30, "Eskimo_Aleut": 0.25, "Patagonian_Southern": 0.30,
    }),

    # SLC16A11 rs13342232 -- Native American diabetes-risk variant (SIGMA
    # consortium 2014). Introgressed from Neanderthals. High in Mesoamericans.
    "rs13342232": ("C", "T", 3.5, {
        "North_American_Algonquian": 0.30, "Mesoamerican": 0.50, "Andean": 0.35,
        "Amazonian": 0.30, "Eskimo_Aleut": 0.20, "Patagonian_Southern": 0.30,
    }),

    # --- Pigmentation gradient (secondary) ---------------------------------

    # TYR rs1042602 -- pigmentation, high in Europeans.
    "rs1042602": ("C", "A", 1.5, {
        "North_American_Algonquian": 0.02, "Mesoamerican": 0.02, "Andean": 0.01,
        "Amazonian": 0.01, "Eskimo_Aleut": 0.02, "Patagonian_Southern": 0.02,
    }),

    # KITLG rs12821256 -- European blonde-hair variant, near zero.
    "rs12821256": ("T", "C", 1.5, {
        "North_American_Algonquian": 0.00, "Mesoamerican": 0.00, "Andean": 0.00,
        "Amazonian": 0.00, "Eskimo_Aleut": 0.00, "Patagonian_Southern": 0.00,
    }),

    # TYRP1 rs2733832 -- brown pigmentation variant.
    "rs2733832": ("T", "C", 1.5, {
        "North_American_Algonquian": 0.50, "Mesoamerican": 0.55, "Andean": 0.55,
        "Amazonian": 0.60, "Eskimo_Aleut": 0.45, "Patagonian_Southern": 0.55,
    }),

    # IRF4 rs12203592 -- light skin/hair, European-enriched.
    "rs12203592": ("C", "T", 1.5, {
        "North_American_Algonquian": 0.00, "Mesoamerican": 0.00, "Andean": 0.00,
        "Amazonian": 0.00, "Eskimo_Aleut": 0.00, "Patagonian_Southern": 0.00,
    }),

    # --- Bitter taste receptors (East Asian / Amerindian profile) ----------

    # TAS2R38 rs713598 -- PAV/AVI bitter taste. Amerindian populations
    # show elevated PAV (taster) frequencies.
    "rs713598": ("G", "C", 1.5, {
        "North_American_Algonquian": 0.60, "Mesoamerican": 0.55, "Andean": 0.55,
        "Amazonian": 0.58, "Eskimo_Aleut": 0.58, "Patagonian_Southern": 0.55,
    }),

    # TAS2R38 rs1726866 -- linked to rs713598.
    "rs1726866": ("G", "A", 1.0, {
        "North_American_Algonquian": 0.58, "Mesoamerican": 0.56, "Andean": 0.56,
        "Amazonian": 0.57, "Eskimo_Aleut": 0.58, "Patagonian_Southern": 0.56,
    }),

    # --- Alcohol metabolism (East Asian / Amerindian signal) ----------------

    # ADH1B rs1229984 -- "alcohol flush". High in East Asians, low but present
    # in some Native American groups (~5-15%).
    "rs1229984": ("G", "A", 2.0, {
        "North_American_Algonquian": 0.05, "Mesoamerican": 0.08, "Andean": 0.06,
        "Amazonian": 0.10, "Eskimo_Aleut": 0.12, "Patagonian_Southern": 0.06,
    }),

    # ALDH2 rs671 -- East Asian flushing variant. Very low in Native Americans.
    "rs671": ("G", "A", 2.0, {
        "North_American_Algonquian": 0.00, "Mesoamerican": 0.01, "Andean": 0.00,
        "Amazonian": 0.00, "Eskimo_Aleut": 0.03, "Patagonian_Southern": 0.00,
    }),

    # --- Immune / HLA-linked markers ---------------------------------------

    # HLA-B rs2523608 -- broad region tag, elevated in some Native populations.
    "rs2523608": ("G", "A", 1.5, {
        "North_American_Algonquian": 0.30, "Mesoamerican": 0.30, "Andean": 0.32,
        "Amazonian": 0.32, "Eskimo_Aleut": 0.25, "Patagonian_Southern": 0.30,
    }),

    # IL13 rs20541 -- asthma/atopy. Elevated in Native Americans and Asians.
    "rs20541": ("C", "T", 1.0, {
        "North_American_Algonquian": 0.40, "Mesoamerican": 0.42, "Andean": 0.45,
        "Amazonian": 0.45, "Eskimo_Aleut": 0.30, "Patagonian_Southern": 0.42,
    }),

    # CCR5 rs333 (delta32) -- essentially absent in Native Americans and Asians.
    "rs333": ("I", "D", 2.0, {
        "North_American_Algonquian": 0.00, "Mesoamerican": 0.00, "Andean": 0.00,
        "Amazonian": 0.00, "Eskimo_Aleut": 0.00, "Patagonian_Southern": 0.00,
    }),

    # --- Metabolic / CYP450 pharmacogenetic markers ------------------------

    # CYP3A5 rs776746 -- high in non-Europeans; moderate in Native Americans.
    "rs776746": ("A", "G", 1.5, {
        "North_American_Algonquian": 0.60, "Mesoamerican": 0.55, "Andean": 0.60,
        "Amazonian": 0.60, "Eskimo_Aleut": 0.70, "Patagonian_Southern": 0.55,
    }),

    # CYP2C19 rs4244285 -- poor-metabolizer allele, moderate in Amerindians.
    "rs4244285": ("G", "A", 1.0, {
        "North_American_Algonquian": 0.15, "Mesoamerican": 0.18, "Andean": 0.15,
        "Amazonian": 0.15, "Eskimo_Aleut": 0.20, "Patagonian_Southern": 0.15,
    }),

    # CYP2D6 rs1065852 -- metabolizer variant, elevated in East Asian / Amerindian.
    "rs1065852": ("G", "A", 1.0, {
        "North_American_Algonquian": 0.55, "Mesoamerican": 0.50, "Andean": 0.55,
        "Amazonian": 0.55, "Eskimo_Aleut": 0.60, "Patagonian_Southern": 0.52,
    }),

    # --- Continental-distinguishing markers (low in Africa / Europe) --------

    # MCM6 rs182549 -- linked to LCT, European lactase persistence.
    "rs182549": ("C", "T", 1.5, {
        "North_American_Algonquian": 0.01, "Mesoamerican": 0.01, "Andean": 0.00,
        "Amazonian": 0.00, "Eskimo_Aleut": 0.02, "Patagonian_Southern": 0.00,
    }),

    # HFE C282Y rs1800562 -- Northern European hemochromatosis variant.
    "rs1800562": ("G", "A", 1.0, {
        "North_American_Algonquian": 0.00, "Mesoamerican": 0.00, "Andean": 0.00,
        "Amazonian": 0.00, "Eskimo_Aleut": 0.00, "Patagonian_Southern": 0.00,
    }),

    # APOL1 G1 (rs73885319) -- African kidney-disease allele, zero in NA.
    "rs73885319": ("A", "G", 1.5, {
        "North_American_Algonquian": 0.00, "Mesoamerican": 0.00, "Andean": 0.00,
        "Amazonian": 0.00, "Eskimo_Aleut": 0.00, "Patagonian_Southern": 0.00,
    }),

    # G6PD A- (rs1050828) -- African/Mediterranean, absent in NA.
    "rs1050828": ("C", "T", 1.0, {
        "North_American_Algonquian": 0.00, "Mesoamerican": 0.00, "Andean": 0.00,
        "Amazonian": 0.00, "Eskimo_Aleut": 0.00, "Patagonian_Southern": 0.00,
    }),

    # HBB sickle rs334 -- African, absent in NA.
    "rs334": ("A", "T", 1.0, {
        "North_American_Algonquian": 0.00, "Mesoamerican": 0.00, "Andean": 0.00,
        "Amazonian": 0.00, "Eskimo_Aleut": 0.00, "Patagonian_Southern": 0.00,
    }),

    # --- mtDNA-linked / admixture proxy markers ----------------------------

    # NADSYN1/DHCR7 rs12785878 -- Vitamin D-related, gradient across populations.
    "rs12785878": ("T", "G", 1.0, {
        "North_American_Algonquian": 0.50, "Mesoamerican": 0.52, "Andean": 0.55,
        "Amazonian": 0.55, "Eskimo_Aleut": 0.50, "Patagonian_Southern": 0.55,
    }),

    # GC (vitamin D binding) rs2282679 -- Amerindians carry elevated allele.
    "rs2282679": ("T", "G", 1.0, {
        "North_American_Algonquian": 0.55, "Mesoamerican": 0.55, "Andean": 0.58,
        "Amazonian": 0.58, "Eskimo_Aleut": 0.50, "Patagonian_Southern": 0.55,
    }),

    # --- Broad AIM panel markers (standard HGDP set) -----------------------

    "rs1834640": ("A", "G", 1.5, {
        "North_American_Algonquian": 0.95, "Mesoamerican": 0.95, "Andean": 0.97,
        "Amazonian": 0.97, "Eskimo_Aleut": 0.90, "Patagonian_Southern": 0.96,
    }),

    "rs1876482": ("T", "C", 1.5, {
        "North_American_Algonquian": 0.85, "Mesoamerican": 0.90, "Andean": 0.92,
        "Amazonian": 0.93, "Eskimo_Aleut": 0.80, "Patagonian_Southern": 0.90,
    }),

    "rs1800498": ("C", "T", 1.0, {
        "North_American_Algonquian": 0.40, "Mesoamerican": 0.42, "Andean": 0.42,
        "Amazonian": 0.43, "Eskimo_Aleut": 0.40, "Patagonian_Southern": 0.42,
    }),

    "rs3814134": ("C", "T", 1.0, {
        "North_American_Algonquian": 0.65, "Mesoamerican": 0.68, "Andean": 0.72,
        "Amazonian": 0.72, "Eskimo_Aleut": 0.62, "Patagonian_Southern": 0.70,
    }),

    "rs2237717": ("C", "T", 1.0, {
        "North_American_Algonquian": 0.35, "Mesoamerican": 0.38, "Andean": 0.42,
        "Amazonian": 0.42, "Eskimo_Aleut": 0.32, "Patagonian_Southern": 0.40,
    }),

    "rs10496971": ("G", "T", 1.0, {
        "North_American_Algonquian": 0.75, "Mesoamerican": 0.78, "Andean": 0.80,
        "Amazonian": 0.80, "Eskimo_Aleut": 0.70, "Patagonian_Southern": 0.78,
    }),

    "rs7554936": ("C", "T", 1.0, {
        "North_American_Algonquian": 0.80, "Mesoamerican": 0.80, "Andean": 0.82,
        "Amazonian": 0.82, "Eskimo_Aleut": 0.75, "Patagonian_Southern": 0.80,
    }),

    "rs2814720": ("A", "G", 1.0, {
        "North_American_Algonquian": 0.60, "Mesoamerican": 0.62, "Andean": 0.65,
        "Amazonian": 0.65, "Eskimo_Aleut": 0.58, "Patagonian_Southern": 0.62,
    }),

    "rs1800404": ("C", "T", 1.0, {
        "North_American_Algonquian": 0.55, "Mesoamerican": 0.55, "Andean": 0.60,
        "Amazonian": 0.60, "Eskimo_Aleut": 0.50, "Patagonian_Southern": 0.58,
    }),

    "rs2065160": ("C", "T", 1.0, {
        "North_American_Algonquian": 0.30, "Mesoamerican": 0.30, "Andean": 0.32,
        "Amazonian": 0.32, "Eskimo_Aleut": 0.28, "Patagonian_Southern": 0.30,
    }),

    "rs12075": ("A", "G", 1.0, {
        "North_American_Algonquian": 0.08, "Mesoamerican": 0.08, "Andean": 0.05,
        "Amazonian": 0.05, "Eskimo_Aleut": 0.10, "Patagonian_Southern": 0.06,
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
    Run weighted NNLS decomposition into Native American sub-population components.

    Returns (proportions_dict, used_snps_list, residual).
    """
    n_pops = len(NATIVE_POPS)
    rows_A = []
    rows_f = []
    rows_w = []
    used_snps = []

    for rsid, (ref, alt, weight, pop_freqs) in AIMS_NATIVE.items():
        variant_data = variants.get(rsid.lower())
        if variant_data is None:
            continue
        _chrom, _pos, genotype = variant_data
        dosage = _get_dosage(genotype, alt)
        if dosage is None:
            continue
        obs_freq = dosage / 2.0
        row = [pop_freqs.get(p, 0.0) for p in NATIVE_POPS]
        rows_A.append(row)
        rows_f.append(obs_freq)
        rows_w.append(weight)
        used_snps.append(rsid)

    if len(used_snps) < 5:
        return {p: 1.0 / n_pops for p in NATIVE_POPS}, used_snps, 0.0

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

    props = {p: float(x[i]) for i, p in enumerate(NATIVE_POPS)}
    return props, used_snps, float(resid)


# ---------------------------------------------------------------------------
# AFFINITY SCORES (min-max normalized distances)
# ---------------------------------------------------------------------------

def _compute_affinity_scores(
    variants: Dict[str, Tuple[str, str, str]],
) -> Dict[str, float]:
    """
    Compute relative affinity to each Native American sub-population using
    min-max normalization. Since Amerindian sub-populations are closely related,
    we normalize within the observed range so the closest scores 100 and
    furthest 0.
    """
    raw_distances: Dict[str, float] = {}

    for pop in NATIVE_POPS:
        weighted_sq_diff = 0.0
        w_sum = 0.0
        for rsid, (ref, alt, weight, pop_freqs) in AIMS_NATIVE.items():
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
        return {p: 0.0 for p in NATIVE_POPS}

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
# KEY NATIVE AMERICAN MARKERS SUMMARY
# ---------------------------------------------------------------------------

def _summarize_key_markers(
    variants: Dict[str, Tuple[str, str, str]],
) -> Dict[str, Dict[str, str]]:
    """Report the key markers distinguishing Native American populations."""
    results: Dict[str, Dict[str, str]] = {}

    # EDAR shovel-tooth (rs3827760)
    variant_data = variants.get("rs3827760")
    if variant_data is not None:
        _c, _p, geno = variant_data
        if geno and geno not in ("--", "00", "NN"):
            g_count = geno.upper().count("G")
            if g_count == 2:
                results["edar_shovel_teeth"] = {
                    "status": "GG",
                    "detail": "GG -- fully derived EDAR 370A; shovel-shaped incisors, thick hair (Native American / East Asian).",
                }
            elif g_count == 1:
                results["edar_shovel_teeth"] = {
                    "status": "AG",
                    "detail": "AG -- one derived EDAR allele (mixed ancestry consistent with Native American contribution).",
                }
            else:
                results["edar_shovel_teeth"] = {
                    "status": "AA",
                    "detail": "AA -- ancestral EDAR; typical of European / African ancestry.",
                }

    # ABCC11 earwax (rs17822931)
    variant_data = variants.get("rs17822931")
    if variant_data is not None:
        _c, _p, geno = variant_data
        if geno and geno not in ("--", "00", "NN"):
            t_count = geno.upper().count("T")
            if t_count == 2:
                results["earwax"] = {
                    "status": "TT",
                    "detail": "TT -- dry earwax. Characteristic of East Asian / Arctic ancestry (Eskimo-Aleut).",
                }
            elif t_count == 1:
                results["earwax"] = {
                    "status": "CT",
                    "detail": "CT -- intermediate; one dry-earwax allele.",
                }
            else:
                results["earwax"] = {
                    "status": "CC",
                    "detail": "CC -- wet earwax (typical of most Americas Indigenous populations outside the Arctic).",
                }

    # Duffy rs2814778 (FY-null)
    variant_data = variants.get("rs2814778")
    if variant_data is not None:
        _c, _p, geno = variant_data
        if geno and geno not in ("--", "00", "NN"):
            c_count = geno.upper().count("C")
            if c_count == 0:
                results["duffy"] = {
                    "status": "TT",
                    "detail": "TT -- Duffy expressed (ancestral). Consistent with Native American ancestry (no African-origin Duffy-null allele).",
                }
            elif c_count == 1:
                results["duffy"] = {
                    "status": "CT",
                    "detail": "CT -- one Duffy-null allele (African contribution).",
                }
            else:
                results["duffy"] = {
                    "status": "CC",
                    "detail": "CC -- Duffy-null (strong African ancestry signal; not typical of un-admixed Native Americans).",
                }

    # LCT lactase persistence (rs4988235)
    variant_data = variants.get("rs4988235")
    if variant_data is not None:
        _c, _p, geno = variant_data
        if geno and geno not in ("--", "00", "NN"):
            a_count = geno.upper().count("A")
            if a_count == 0:
                results["lactase"] = {
                    "status": "GG",
                    "detail": "GG -- lactase non-persistent (ancestral). Universal pattern among pre-contact Native Americans.",
                }
            elif a_count == 1:
                results["lactase"] = {
                    "status": "AG",
                    "detail": "AG -- one European-origin lactase persistence allele.",
                }
            else:
                results["lactase"] = {
                    "status": "AA",
                    "detail": "AA -- lactase persistent. Indicates European admixture.",
                }

    # ABCA1 "Amerindian allele" (rs9282541)
    variant_data = variants.get("rs9282541")
    if variant_data is not None:
        _c, _p, geno = variant_data
        if geno and geno not in ("--", "00", "NN"):
            t_count = geno.upper().count("T")
            if t_count > 0:
                status = "CT" if t_count == 1 else "TT"
                results["abca1_amerindian"] = {
                    "status": status,
                    "detail": "T allele present -- 'Amerindian' ABCA1 R230C variant; characteristic of Native American ancestry (especially Mesoamerican).",
                }
            else:
                results["abca1_amerindian"] = {
                    "status": "CC",
                    "detail": "CC -- ancestral allele (no Amerindian ABCA1 signal).",
                }

    # CPT1A Arctic variant (rs80356779)
    variant_data = variants.get("rs80356779")
    if variant_data is not None:
        _c, _p, geno = variant_data
        if geno and geno not in ("--", "00", "NN"):
            a_count = geno.upper().count("A")
            if a_count > 0:
                status = "GA" if a_count == 1 else "AA"
                results["cpt1a_arctic"] = {
                    "status": status,
                    "detail": "A allele present -- Arctic-specific CPT1A P479L variant (Inuit / Yupik / Aleut marine-diet adaptation).",
                }
            else:
                results["cpt1a_arctic"] = {
                    "status": "GG",
                    "detail": "GG -- ancestral; no Arctic adaptation signal.",
                }

    return results


# ---------------------------------------------------------------------------
# PUBLIC API
# ---------------------------------------------------------------------------

def analyze_native_american(
    variants: Dict[str, Tuple[str, str, str]],
) -> Dict[str, Any]:
    """
    Analyze variants for Native American population affinity.

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


def generate_native_american_json(result: Dict[str, Any]) -> Dict[str, Any]:
    """
    Convert raw analysis result into the structured JSON for the frontend.

    Args:
        result: Output from analyze_native_american()

    Returns:
        Structured JSON dict matching the report schema.
    """
    props = result["proportions"]
    affinity_scores = result["affinity_scores"]
    used_snps = result["used_snps"]
    resid = result["residual"]
    key_markers = result["key_markers"]

    panel_size = len(AIMS_NATIVE)
    snps_used = len(used_snps)

    sorted_pops = sorted(NATIVE_POPS, key=lambda p: -affinity_scores.get(p, 0.0))

    populations = []
    for pop in sorted_pops:
        proportion = round(props.get(pop, 0.0), 4)
        populations.append({
            "code": pop,
            "label": NATIVE_POP_LABELS[pop],
            "affinityScore": affinity_scores.get(pop, 0.0),
            "proportion": proportion,
            "percentage": round(proportion * 100, 1),
            "description": NATIVE_POP_DESCRIPTIONS[pop],
            "color": NATIVE_POP_COLORS[pop],
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
            "label": NATIVE_POP_LABELS[top_pop],
            "description": NATIVE_POP_DESCRIPTIONS[top_pop],
        },
        "keyMarkers": key_markers,
    }
