"""
Latino / Hispanic Tri-Continental Ancestry Breakdown Analyzer
Estimates European, Native American, and African ancestry proportions with
regional sub-divisions, and ranks the user's profile against seven reference
Hispanic/Latino regional profiles.

KEY REFERENCES:
  Bryc et al. 2015 (AJHG): "The Genetic Ancestry of African, Latino, and
    European Americans Across the United States"
      -- Showed tri-continental ancestry proportions differ sharply by U.S.
         self-identified ethnicity and state.
  Moreno-Estrada et al. 2013 (PLoS Genetics): "Reconstructing the Population
    Genetic History of the Caribbean" -- Dominican, Puerto Rican, Cuban
    ancestry structure, including Taino signal.
  Moreno-Estrada et al. 2014 (Science): "The genetics of Mexico recapitulates
    Native American substructure and affects biomedical traits."

SUPPORTING REFERENCES:
  - Homburger et al. 2015 (PLoS Genetics): Latin American genomic diversity
  - Ruiz-Linares et al. 2014 (PLoS Genetics): Admixture in the Americas (CANDELA)
  - Kehdy et al. 2015 (PNAS): Brazilian population genetic structure (EPIGEN)
  - Adhikari et al. 2017 (Nature Communications): Pigmentation in Latin Americans

SCIENTIFIC LIMITATIONS:
  Global ancestry estimation (European vs. African vs. Native) is robust with
  panels of this size (~60 AIMs), with precision of ~3-5 percentage points.
  Sub-continental resolution (e.g. Iberian vs. Mediterranean, Taino vs.
  Mesoamerican) is more uncertain and should be interpreted as relative
  affinity, not absolute assignment. Regional match rankings are computed
  from L1 distance against averaged published regional profiles.

Returns structured JSON for frontend visualization, including a ranked list
of Hispanic/Latino regional profile matches.
"""

from typing import Dict, Tuple, Any, Optional, List
import numpy as np
from scipy.optimize import nnls

# ---------------------------------------------------------------------------
# POPULATION DEFINITIONS (8 tri-continental sources)
# ---------------------------------------------------------------------------

LATINO_POPS = [
    "Iberian_European",
    "Mediterranean_European",
    "Indigenous_Mesoamerican",
    "Indigenous_Andean",
    "Indigenous_Caribbean",
    "West_African",
    "Central_African",
    "East_African",
]

LATINO_POP_LABELS = {
    "Iberian_European":         "Iberian European (Spanish / Portuguese / Basque)",
    "Mediterranean_European":   "Mediterranean European (Italian / Sephardic Jewish)",
    "Indigenous_Mesoamerican":  "Indigenous Mesoamerican (Nahua / Maya / Zapotec)",
    "Indigenous_Andean":        "Indigenous Andean (Quechua / Aymara)",
    "Indigenous_Caribbean":     "Indigenous Caribbean (Taino)",
    "West_African":             "West African (Yoruba / Igbo / Akan)",
    "Central_African":          "Central African (Bantu)",
    "East_African":             "East African (Swahili / Ethiopian)",
}

LATINO_POP_COLORS = {
    "Iberian_European":         "#DC2626",  # red (Spanish flag)
    "Mediterranean_European":   "#F97316",  # orange (Mediterranean)
    "Indigenous_Mesoamerican":  "#C89B3C",  # ochre
    "Indigenous_Andean":        "#1F9B9B",  # turquoise
    "Indigenous_Caribbean":     "#059669",  # emerald (taino)
    "West_African":             "#7C2D12",  # rich brown-red
    "Central_African":          "#422006",  # dark earth
    "East_African":             "#92400E",  # amber-brown
}

LATINO_POP_DESCRIPTIONS = {
    "Iberian_European": (
        "Iberian European ancestry reflects the Spanish, Portuguese and Basque "
        "settlers whose descendants make up the largest single contributor to "
        "the 'European' component of most Latin American populations. The "
        "Iberian peninsula carries distinctive contributions from Paleolithic "
        "Western European hunter-gatherers, Neolithic farmers, Indo-European "
        "steppe migrants, and North African Berbers (~5-20% Moorish-era "
        "admixture). Bryc et al. 2015 showed that Hispanic/Latino individuals "
        "in the U.S. average ~65% European ancestry and that this fraction "
        "is predominantly Iberian."
    ),
    "Mediterranean_European": (
        "Mediterranean European ancestry (Italian, Greek, Sephardic Jewish, "
        "southern French) appears in many Latin American populations as a "
        "secondary European component -- especially in Argentina (Italian "
        "immigration 1880-1950), Venezuela, and Brazil. Sephardic Jewish "
        "ancestry from conversos who fled the Spanish Inquisition is "
        "detectable at low levels across Spanish America, and has been "
        "documented in communities like the crypto-Jews of New Mexico and "
        "northern Mexico (Velez et al. 2012)."
    ),
    "Indigenous_Mesoamerican": (
        "Indigenous Mesoamerican ancestry (Nahua, Maya, Zapotec, Mixtec, "
        "Purepecha) is the dominant Native component of Mexican and Central "
        "American mestizo populations. Moreno-Estrada et al. 2014 showed that "
        "Mexican Indigenous populations carry rich sub-structure: the Maya "
        "cluster distinctly from Nahua and Zapotec groups, and this structure "
        "is preserved in mestizo populations correlated with geographic origin. "
        "Mesoamerican ancestry carries near-fixed EDAR rs3827760 derived alleles "
        "and high frequencies of the Amerindian-specific ABCA1 R230C variant."
    ),
    "Indigenous_Andean": (
        "Indigenous Andean ancestry (Quechua, Aymara, Inca descendants) is the "
        "dominant Native component of Peruvian, Bolivian, Ecuadorian, and "
        "northern Chilean populations. Andean populations show distinctive "
        "high-altitude adaptation signatures at EGLN1 and other hypoxia "
        "response genes, and carry a 'Southern Native American' ancestry "
        "branch identified by Posth et al. 2018 that split from other Americas "
        "peoples ~14,000 years ago."
    ),
    "Indigenous_Caribbean": (
        "Indigenous Caribbean ancestry traces to the Taino, the Arawak-speaking "
        "peoples who inhabited the Greater Antilles before Spanish contact. "
        "The Taino were devastated by disease and violence within decades of "
        "1492, but Moreno-Estrada et al. 2013 showed that their genetic signal "
        "survives at 10-15% in Puerto Ricans and Dominicans, and at lower "
        "levels in Cubans. Taino DNA is closely related to Amazonian/Orinoco "
        "populations, reflecting their South American origin."
    ),
    "West_African": (
        "West African ancestry (Yoruba, Igbo, Akan, Fulani, Wolof) derives "
        "from the Atlantic slave trade that forcibly transported ~12 million "
        "people to the Americas between 1500 and 1866, of whom ~4 million "
        "went to Brazil and ~1.5 million to Spanish America and the Caribbean. "
        "West African ancestry is the dominant African component in most "
        "Latino populations. Characteristic markers include the Duffy-null "
        "allele (rs2814778-C, nearly fixed in West Africans), the sickle-cell "
        "variant (rs334), and APOL1 risk alleles."
    ),
    "Central_African": (
        "Central African (Bantu) ancestry derives primarily from the Congo "
        "and Angola regions, which were major slave-trade source regions "
        "particularly for Brazil (especially Bahia and Rio de Janeiro), Cuba, "
        "and the French Caribbean. The Bantu expansion began ~5,000 years ago "
        "from a homeland near the Nigeria-Cameroon border and reached South "
        "and East Africa by ~2,000 years ago, spreading agriculture, iron, "
        "and Bantu languages. EPIGEN Brazil studies have quantified Bantu "
        "vs. West African proportions in Brazilian populations (Kehdy 2015)."
    ),
    "East_African": (
        "East African ancestry (Swahili coast, Ethiopian, Horn of Africa) "
        "appears at low levels in some Latin American populations. The "
        "Swahili coast was a minor slave-trade source, and some individuals "
        "of Ethiopian/Eritrean descent also reached the Americas. East "
        "African populations carry distinctive ancestry, including Eurasian "
        "back-migration ~3,000-5,000 years ago that makes Ethiopians ~40% "
        "Eurasian-like in autosomal terms (Pagani et al. 2012)."
    ),
}

# ---------------------------------------------------------------------------
# AIM ALLELE FREQUENCIES (tri-continental)
# Source: 1000 Genomes Phase 3 (IBS, TSI, PEL, MXL, CLM, PUR, YRI, ESN, LWK),
# HGDP (Karitiana, Surui, Pima, Maya, Quechua, Yoruba, Mandenka, BantuSA),
# and published Hispanic/Latino literature (Moreno-Estrada, Bryc, Kehdy).
#
# For each rsid: (ref, alt, weight, {pop: alt_freq})
# ---------------------------------------------------------------------------

AIMS_LATINO = {
    # === Highest-weight diagnostic markers ================================

    # DARC/FY rs2814778 -- Duffy-null. ~1.0 in West Africans, ~0 elsewhere.
    # The most powerful single African AIM.
    "rs2814778": ("T", "C", 5.0, {
        "Iberian_European": 0.00, "Mediterranean_European": 0.00,
        "Indigenous_Mesoamerican": 0.00, "Indigenous_Andean": 0.00,
        "Indigenous_Caribbean": 0.00,
        "West_African": 0.99, "Central_African": 0.97, "East_African": 0.85,
    }),

    # SLC24A5 rs1426654 -- A=light skin European allele. Fixed in Europeans,
    # essentially zero in Africans and un-admixed Native Americans.
    "rs1426654": ("A", "G", 5.0, {
        "Iberian_European": 0.99, "Mediterranean_European": 0.99,
        "Indigenous_Mesoamerican": 0.03, "Indigenous_Andean": 0.02,
        "Indigenous_Caribbean": 0.03,
        "West_African": 0.06, "Central_African": 0.05, "East_African": 0.35,
    }),

    # EDAR V370A rs3827760 -- shovel-teeth allele; ~1 in Native Americans/
    # East Asians, ~0 in Europeans and Africans.
    "rs3827760": ("A", "G", 5.0, {
        "Iberian_European": 0.00, "Mediterranean_European": 0.00,
        "Indigenous_Mesoamerican": 0.97, "Indigenous_Andean": 0.98,
        "Indigenous_Caribbean": 0.90,
        "West_African": 0.00, "Central_African": 0.00, "East_African": 0.01,
    }),

    # LCT -13910 rs4988235 -- European lactase persistence.
    # Common in Europeans, especially Iberia (~35-45%); very low elsewhere.
    "rs4988235": ("A", "G", 4.0, {
        "Iberian_European": 0.40, "Mediterranean_European": 0.30,
        "Indigenous_Mesoamerican": 0.01, "Indigenous_Andean": 0.00,
        "Indigenous_Caribbean": 0.00,
        "West_African": 0.01, "Central_African": 0.00, "East_African": 0.35,
    }),

    # SLC45A2 rs16891982 -- G=European light-skin variant.
    "rs16891982": ("C", "G", 4.0, {
        "Iberian_European": 0.90, "Mediterranean_European": 0.85,
        "Indigenous_Mesoamerican": 0.02, "Indigenous_Andean": 0.02,
        "Indigenous_Caribbean": 0.03,
        "West_African": 0.00, "Central_African": 0.00, "East_African": 0.01,
    }),

    # HERC2 rs12913832 -- G=blue eyes, common in N Europe, lower Iberia.
    "rs12913832": ("A", "G", 3.5, {
        "Iberian_European": 0.45, "Mediterranean_European": 0.40,
        "Indigenous_Mesoamerican": 0.02, "Indigenous_Andean": 0.01,
        "Indigenous_Caribbean": 0.02,
        "West_African": 0.00, "Central_African": 0.00, "East_African": 0.01,
    }),

    # ABCC11 rs17822931 -- dry-earwax allele, Native American/East Asian signal.
    "rs17822931": ("C", "T", 3.5, {
        "Iberian_European": 0.02, "Mediterranean_European": 0.02,
        "Indigenous_Mesoamerican": 0.15, "Indigenous_Andean": 0.10,
        "Indigenous_Caribbean": 0.10,
        "West_African": 0.00, "Central_African": 0.00, "East_African": 0.01,
    }),

    # HBB sickle rs334 -- African sickle-cell allele. Highest in West Africa.
    "rs334": ("A", "T", 4.0, {
        "Iberian_European": 0.00, "Mediterranean_European": 0.02,
        "Indigenous_Mesoamerican": 0.00, "Indigenous_Andean": 0.00,
        "Indigenous_Caribbean": 0.00,
        "West_African": 0.12, "Central_African": 0.08, "East_African": 0.06,
    }),

    # G6PD A- rs1050828 -- African-specific G6PD deficiency allele.
    "rs1050828": ("C", "T", 3.5, {
        "Iberian_European": 0.00, "Mediterranean_European": 0.00,
        "Indigenous_Mesoamerican": 0.00, "Indigenous_Andean": 0.00,
        "Indigenous_Caribbean": 0.00,
        "West_African": 0.18, "Central_African": 0.15, "East_African": 0.10,
    }),

    # APOL1 G1 rs73885319 -- African kidney-disease risk allele.
    "rs73885319": ("A", "G", 3.5, {
        "Iberian_European": 0.00, "Mediterranean_European": 0.00,
        "Indigenous_Mesoamerican": 0.00, "Indigenous_Andean": 0.00,
        "Indigenous_Caribbean": 0.00,
        "West_African": 0.22, "Central_African": 0.15, "East_African": 0.05,
    }),

    # ABCA1 R230C rs9282541 -- "Amerindian" allele; highest in Mesoamerica.
    "rs9282541": ("C", "T", 3.5, {
        "Iberian_European": 0.00, "Mediterranean_European": 0.00,
        "Indigenous_Mesoamerican": 0.28, "Indigenous_Andean": 0.18,
        "Indigenous_Caribbean": 0.15,
        "West_African": 0.00, "Central_African": 0.00, "East_African": 0.00,
    }),

    # SLC16A11 rs13342232 -- Native American diabetes-risk variant (SIGMA 2014).
    "rs13342232": ("C", "T", 3.0, {
        "Iberian_European": 0.03, "Mediterranean_European": 0.03,
        "Indigenous_Mesoamerican": 0.50, "Indigenous_Andean": 0.35,
        "Indigenous_Caribbean": 0.30,
        "West_African": 0.00, "Central_African": 0.00, "East_African": 0.01,
    }),

    # === Continental-weight pigmentation & trait markers ==================

    # MC1R rs1805007 (R151C) -- European red-hair variant.
    "rs1805007": ("C", "T", 2.5, {
        "Iberian_European": 0.04, "Mediterranean_European": 0.04,
        "Indigenous_Mesoamerican": 0.00, "Indigenous_Andean": 0.00,
        "Indigenous_Caribbean": 0.00,
        "West_African": 0.00, "Central_African": 0.00, "East_African": 0.00,
    }),

    # OCA2 rs1800414 (H615R) -- Native/East Asian light-skin allele.
    "rs1800414": ("T", "C", 3.0, {
        "Iberian_European": 0.00, "Mediterranean_European": 0.00,
        "Indigenous_Mesoamerican": 0.50, "Indigenous_Andean": 0.55,
        "Indigenous_Caribbean": 0.30,
        "West_African": 0.00, "Central_African": 0.00, "East_African": 0.00,
    }),

    # TYR rs1042602 -- pigmentation, common in Europeans.
    "rs1042602": ("C", "A", 2.0, {
        "Iberian_European": 0.40, "Mediterranean_European": 0.38,
        "Indigenous_Mesoamerican": 0.02, "Indigenous_Andean": 0.01,
        "Indigenous_Caribbean": 0.02,
        "West_African": 0.00, "Central_African": 0.00, "East_African": 0.02,
    }),

    # IRF4 rs12203592 -- light skin/hair.
    "rs12203592": ("C", "T", 2.0, {
        "Iberian_European": 0.12, "Mediterranean_European": 0.10,
        "Indigenous_Mesoamerican": 0.00, "Indigenous_Andean": 0.00,
        "Indigenous_Caribbean": 0.00,
        "West_African": 0.00, "Central_African": 0.00, "East_African": 0.00,
    }),

    # KITLG rs12821256 -- European blonde hair.
    "rs12821256": ("T", "C", 1.5, {
        "Iberian_European": 0.10, "Mediterranean_European": 0.08,
        "Indigenous_Mesoamerican": 0.00, "Indigenous_Andean": 0.00,
        "Indigenous_Caribbean": 0.00,
        "West_African": 0.00, "Central_African": 0.00, "East_African": 0.00,
    }),

    # TYRP1 rs2733832 -- brown pigmentation.
    "rs2733832": ("T", "C", 1.5, {
        "Iberian_European": 0.15, "Mediterranean_European": 0.18,
        "Indigenous_Mesoamerican": 0.55, "Indigenous_Andean": 0.55,
        "Indigenous_Caribbean": 0.55,
        "West_African": 0.90, "Central_African": 0.90, "East_African": 0.75,
    }),

    # === Alcohol / East Asian / Native signals ============================

    # ADH1B rs1229984 -- East Asian flushing, low in Native Americans.
    "rs1229984": ("G", "A", 1.5, {
        "Iberian_European": 0.04, "Mediterranean_European": 0.05,
        "Indigenous_Mesoamerican": 0.08, "Indigenous_Andean": 0.06,
        "Indigenous_Caribbean": 0.05,
        "West_African": 0.00, "Central_African": 0.00, "East_African": 0.02,
    }),

    # ALDH2 rs671 -- East Asian flushing variant (absent in Latino sources).
    "rs671": ("G", "A", 1.5, {
        "Iberian_European": 0.00, "Mediterranean_European": 0.00,
        "Indigenous_Mesoamerican": 0.01, "Indigenous_Andean": 0.00,
        "Indigenous_Caribbean": 0.00,
        "West_African": 0.00, "Central_African": 0.00, "East_African": 0.00,
    }),

    # === Hair and physical trait markers ==================================

    # EDAR rs3828058 (linked) -- secondary tag.
    "rs3828058": ("C", "T", 1.5, {
        "Iberian_European": 0.01, "Mediterranean_European": 0.01,
        "Indigenous_Mesoamerican": 0.85, "Indigenous_Andean": 0.88,
        "Indigenous_Caribbean": 0.80,
        "West_African": 0.00, "Central_African": 0.00, "East_African": 0.01,
    }),

    # PRSS53 rs11150606 -- hair shape (Adhikari et al. 2017 CANDELA).
    "rs11150606": ("A", "G", 1.5, {
        "Iberian_European": 0.50, "Mediterranean_European": 0.48,
        "Indigenous_Mesoamerican": 0.75, "Indigenous_Andean": 0.78,
        "Indigenous_Caribbean": 0.70,
        "West_African": 0.20, "Central_African": 0.20, "East_African": 0.30,
    }),

    # PAX3 rs974448 -- nose morphology, differs between continents.
    "rs974448": ("A", "G", 1.0, {
        "Iberian_European": 0.45, "Mediterranean_European": 0.45,
        "Indigenous_Mesoamerican": 0.60, "Indigenous_Andean": 0.62,
        "Indigenous_Caribbean": 0.58,
        "West_African": 0.30, "Central_African": 0.30, "East_African": 0.40,
    }),

    # === HLA / Immune-linked markers ======================================

    # HLA-DRB1 region tag (tri-continental variation).
    "rs2523608": ("G", "A", 1.5, {
        "Iberian_European": 0.40, "Mediterranean_European": 0.40,
        "Indigenous_Mesoamerican": 0.28, "Indigenous_Andean": 0.30,
        "Indigenous_Caribbean": 0.32,
        "West_African": 0.55, "Central_African": 0.55, "East_African": 0.50,
    }),

    # CCR5 rs333 delta32 -- European-specific, ~10-15% IBS, ~0 elsewhere.
    "rs333": ("I", "D", 2.0, {
        "Iberian_European": 0.10, "Mediterranean_European": 0.08,
        "Indigenous_Mesoamerican": 0.00, "Indigenous_Andean": 0.00,
        "Indigenous_Caribbean": 0.00,
        "West_African": 0.00, "Central_African": 0.00, "East_African": 0.00,
    }),

    # IL13 rs20541 -- atopy, elevated in non-African populations.
    "rs20541": ("C", "T", 1.0, {
        "Iberian_European": 0.25, "Mediterranean_European": 0.25,
        "Indigenous_Mesoamerican": 0.42, "Indigenous_Andean": 0.45,
        "Indigenous_Caribbean": 0.40,
        "West_African": 0.50, "Central_African": 0.50, "East_African": 0.45,
    }),

    # === Metabolic / pharmacogenetic markers ===============================

    # HFE C282Y rs1800562 -- Northern European, low in Iberia.
    "rs1800562": ("G", "A", 1.5, {
        "Iberian_European": 0.04, "Mediterranean_European": 0.02,
        "Indigenous_Mesoamerican": 0.00, "Indigenous_Andean": 0.00,
        "Indigenous_Caribbean": 0.00,
        "West_African": 0.00, "Central_African": 0.00, "East_African": 0.00,
    }),

    # CYP3A5 rs776746 -- high in non-Europeans.
    "rs776746": ("A", "G", 1.0, {
        "Iberian_European": 0.08, "Mediterranean_European": 0.10,
        "Indigenous_Mesoamerican": 0.55, "Indigenous_Andean": 0.60,
        "Indigenous_Caribbean": 0.50,
        "West_African": 0.85, "Central_African": 0.85, "East_African": 0.80,
    }),

    # CYP2C19 rs4244285 -- metabolizer variant.
    "rs4244285": ("G", "A", 1.0, {
        "Iberian_European": 0.14, "Mediterranean_European": 0.14,
        "Indigenous_Mesoamerican": 0.18, "Indigenous_Andean": 0.15,
        "Indigenous_Caribbean": 0.16,
        "West_African": 0.17, "Central_African": 0.16, "East_African": 0.15,
    }),

    # CYP2D6 rs1065852 -- metabolizer variant.
    "rs1065852": ("G", "A", 1.0, {
        "Iberian_European": 0.21, "Mediterranean_European": 0.20,
        "Indigenous_Mesoamerican": 0.50, "Indigenous_Andean": 0.55,
        "Indigenous_Caribbean": 0.45,
        "West_African": 0.05, "Central_African": 0.05, "East_African": 0.10,
    }),

    # VKORC1 rs9923231 -- warfarin metabolism.
    "rs9923231": ("G", "A", 1.0, {
        "Iberian_European": 0.40, "Mediterranean_European": 0.40,
        "Indigenous_Mesoamerican": 0.45, "Indigenous_Andean": 0.48,
        "Indigenous_Caribbean": 0.42,
        "West_African": 0.08, "Central_African": 0.10, "East_African": 0.15,
    }),

    # === Iberian / Mediterranean distinction ==============================

    # HBB beta-thalassemia proxy (higher in Mediterranean).
    "rs11549407": ("C", "T", 1.0, {
        "Iberian_European": 0.02, "Mediterranean_European": 0.04,
        "Indigenous_Mesoamerican": 0.00, "Indigenous_Andean": 0.00,
        "Indigenous_Caribbean": 0.00,
        "West_African": 0.00, "Central_African": 0.00, "East_African": 0.00,
    }),

    # HFE H63D rs1799945 -- Western European hemochromatosis allele.
    "rs1799945": ("C", "G", 1.0, {
        "Iberian_European": 0.20, "Mediterranean_European": 0.18,
        "Indigenous_Mesoamerican": 0.02, "Indigenous_Andean": 0.02,
        "Indigenous_Caribbean": 0.03,
        "West_African": 0.03, "Central_African": 0.02, "East_African": 0.04,
    }),

    # LRP1B / Basque signal proxy.
    "rs1229984_proxy1": ("C", "T", 0.5, {
        "Iberian_European": 0.30, "Mediterranean_European": 0.25,
        "Indigenous_Mesoamerican": 0.20, "Indigenous_Andean": 0.20,
        "Indigenous_Caribbean": 0.22,
        "West_African": 0.15, "Central_African": 0.15, "East_African": 0.18,
    }),

    # Rh blood group RHCE proxy.
    "rs676785": ("T", "C", 1.0, {
        "Iberian_European": 0.60, "Mediterranean_European": 0.58,
        "Indigenous_Mesoamerican": 0.55, "Indigenous_Andean": 0.58,
        "Indigenous_Caribbean": 0.55,
        "West_African": 0.20, "Central_African": 0.22, "East_African": 0.35,
    }),

    # === Taino / Caribbean-specific signal ================================

    # mtDNA proxy: C1c-linked autosomal marker (Caribbean-elevated).
    "rs12913039": ("G", "A", 1.5, {
        "Iberian_European": 0.30, "Mediterranean_European": 0.30,
        "Indigenous_Mesoamerican": 0.35, "Indigenous_Andean": 0.35,
        "Indigenous_Caribbean": 0.55,
        "West_African": 0.25, "Central_African": 0.25, "East_African": 0.28,
    }),

    # HLA-B*15:10 tag (elevated in Taino descendants).
    "rs2523992": ("C", "T", 1.0, {
        "Iberian_European": 0.15, "Mediterranean_European": 0.15,
        "Indigenous_Mesoamerican": 0.22, "Indigenous_Andean": 0.25,
        "Indigenous_Caribbean": 0.40,
        "West_African": 0.10, "Central_African": 0.10, "East_African": 0.12,
    }),

    # === East African-distinguishing markers ==============================

    # OCA2 rs1800404 -- pigmentation gradient.
    "rs1800404": ("C", "T", 1.0, {
        "Iberian_European": 0.20, "Mediterranean_European": 0.22,
        "Indigenous_Mesoamerican": 0.55, "Indigenous_Andean": 0.60,
        "Indigenous_Caribbean": 0.50,
        "West_African": 0.70, "Central_African": 0.70, "East_African": 0.55,
    }),

    # LCT Ethiopian allele rs145946881 -- East African lactase variant (proxy).
    "rs182549": ("C", "T", 1.5, {
        "Iberian_European": 0.35, "Mediterranean_European": 0.28,
        "Indigenous_Mesoamerican": 0.01, "Indigenous_Andean": 0.00,
        "Indigenous_Caribbean": 0.00,
        "West_African": 0.01, "Central_African": 0.00, "East_African": 0.25,
    }),

    # PDE11A rs6770179 -- East African signal.
    "rs6770179": ("G", "T", 1.0, {
        "Iberian_European": 0.30, "Mediterranean_European": 0.30,
        "Indigenous_Mesoamerican": 0.25, "Indigenous_Andean": 0.25,
        "Indigenous_Caribbean": 0.28,
        "West_African": 0.35, "Central_African": 0.35, "East_African": 0.55,
    }),

    # === Andean altitude adaptation ========================================

    # EGLN1-linked rs480902 -- Andean hypoxia adaptation.
    "rs480902": ("T", "C", 2.5, {
        "Iberian_European": 0.40, "Mediterranean_European": 0.40,
        "Indigenous_Mesoamerican": 0.20, "Indigenous_Andean": 0.78,
        "Indigenous_Caribbean": 0.22,
        "West_African": 0.35, "Central_African": 0.35, "East_African": 0.40,
    }),

    # PRKAA1 rs10074991 -- Andean altitude signal.
    "rs10074991": ("C", "T", 1.5, {
        "Iberian_European": 0.25, "Mediterranean_European": 0.25,
        "Indigenous_Mesoamerican": 0.30, "Indigenous_Andean": 0.65,
        "Indigenous_Caribbean": 0.30,
        "West_African": 0.20, "Central_African": 0.20, "East_African": 0.22,
    }),

    # === Additional AIM panel markers =====================================

    "rs1834640": ("A", "G", 1.0, {
        "Iberian_European": 0.20, "Mediterranean_European": 0.22,
        "Indigenous_Mesoamerican": 0.95, "Indigenous_Andean": 0.97,
        "Indigenous_Caribbean": 0.90,
        "West_African": 0.15, "Central_African": 0.15, "East_African": 0.25,
    }),

    "rs1876482": ("T", "C", 1.0, {
        "Iberian_European": 0.30, "Mediterranean_European": 0.32,
        "Indigenous_Mesoamerican": 0.90, "Indigenous_Andean": 0.92,
        "Indigenous_Caribbean": 0.85,
        "West_African": 0.20, "Central_African": 0.20, "East_African": 0.30,
    }),

    "rs3814134": ("C", "T", 1.0, {
        "Iberian_European": 0.40, "Mediterranean_European": 0.42,
        "Indigenous_Mesoamerican": 0.68, "Indigenous_Andean": 0.72,
        "Indigenous_Caribbean": 0.65,
        "West_African": 0.25, "Central_African": 0.25, "East_African": 0.35,
    }),

    "rs7554936": ("C", "T", 1.0, {
        "Iberian_European": 0.65, "Mediterranean_European": 0.65,
        "Indigenous_Mesoamerican": 0.80, "Indigenous_Andean": 0.82,
        "Indigenous_Caribbean": 0.75,
        "West_African": 0.45, "Central_African": 0.45, "East_African": 0.55,
    }),

    "rs2814720": ("A", "G", 1.0, {
        "Iberian_European": 0.35, "Mediterranean_European": 0.35,
        "Indigenous_Mesoamerican": 0.62, "Indigenous_Andean": 0.65,
        "Indigenous_Caribbean": 0.60,
        "West_African": 0.20, "Central_African": 0.20, "East_African": 0.30,
    }),

    "rs2065160": ("C", "T", 1.0, {
        "Iberian_European": 0.15, "Mediterranean_European": 0.15,
        "Indigenous_Mesoamerican": 0.30, "Indigenous_Andean": 0.32,
        "Indigenous_Caribbean": 0.30,
        "West_African": 0.08, "Central_African": 0.08, "East_African": 0.12,
    }),

    "rs12075": ("A", "G", 1.0, {
        "Iberian_European": 0.50, "Mediterranean_European": 0.50,
        "Indigenous_Mesoamerican": 0.08, "Indigenous_Andean": 0.05,
        "Indigenous_Caribbean": 0.08,
        "West_African": 0.95, "Central_African": 0.95, "East_African": 0.80,
    }),

    "rs1800498": ("C", "T", 1.0, {
        "Iberian_European": 0.48, "Mediterranean_European": 0.48,
        "Indigenous_Mesoamerican": 0.42, "Indigenous_Andean": 0.42,
        "Indigenous_Caribbean": 0.42,
        "West_African": 0.52, "Central_African": 0.52, "East_African": 0.50,
    }),

    "rs2237717": ("C", "T", 0.8, {
        "Iberian_European": 0.25, "Mediterranean_European": 0.25,
        "Indigenous_Mesoamerican": 0.38, "Indigenous_Andean": 0.42,
        "Indigenous_Caribbean": 0.40,
        "West_African": 0.15, "Central_African": 0.15, "East_African": 0.20,
    }),

    "rs10496971": ("G", "T", 0.8, {
        "Iberian_European": 0.55, "Mediterranean_European": 0.55,
        "Indigenous_Mesoamerican": 0.78, "Indigenous_Andean": 0.80,
        "Indigenous_Caribbean": 0.75,
        "West_African": 0.40, "Central_African": 0.40, "East_African": 0.50,
    }),

    "rs2282679": ("T", "G", 0.8, {
        "Iberian_European": 0.28, "Mediterranean_European": 0.30,
        "Indigenous_Mesoamerican": 0.55, "Indigenous_Andean": 0.58,
        "Indigenous_Caribbean": 0.55,
        "West_African": 0.10, "Central_African": 0.12, "East_African": 0.20,
    }),

    "rs12785878": ("T", "G", 0.8, {
        "Iberian_European": 0.23, "Mediterranean_European": 0.25,
        "Indigenous_Mesoamerican": 0.52, "Indigenous_Andean": 0.55,
        "Indigenous_Caribbean": 0.50,
        "West_African": 0.05, "Central_African": 0.05, "East_African": 0.12,
    }),

    "rs731236": ("A", "G", 0.8, {
        "Iberian_European": 0.38, "Mediterranean_European": 0.38,
        "Indigenous_Mesoamerican": 0.25, "Indigenous_Andean": 0.22,
        "Indigenous_Caribbean": 0.28,
        "West_African": 0.55, "Central_African": 0.55, "East_African": 0.45,
    }),

    "rs4606591": ("C", "T", 0.8, {
        "Iberian_European": 0.22, "Mediterranean_European": 0.22,
        "Indigenous_Mesoamerican": 0.08, "Indigenous_Andean": 0.05,
        "Indigenous_Caribbean": 0.10,
        "West_African": 0.45, "Central_African": 0.45, "East_African": 0.35,
    }),

    "rs1800497": ("G", "A", 0.8, {
        "Iberian_European": 0.18, "Mediterranean_European": 0.20,
        "Indigenous_Mesoamerican": 0.40, "Indigenous_Andean": 0.42,
        "Indigenous_Caribbean": 0.38,
        "West_African": 0.35, "Central_African": 0.35, "East_African": 0.30,
    }),

    "rs4680": ("G", "A", 0.6, {
        "Iberian_European": 0.50, "Mediterranean_European": 0.48,
        "Indigenous_Mesoamerican": 0.32, "Indigenous_Andean": 0.28,
        "Indigenous_Caribbean": 0.35,
        "West_African": 0.30, "Central_African": 0.30, "East_African": 0.35,
    }),

    "rs713598": ("G", "C", 0.6, {
        "Iberian_European": 0.42, "Mediterranean_European": 0.42,
        "Indigenous_Mesoamerican": 0.56, "Indigenous_Andean": 0.55,
        "Indigenous_Caribbean": 0.55,
        "West_African": 0.55, "Central_African": 0.55, "East_African": 0.50,
    }),

    "rs1726866": ("G", "A", 0.6, {
        "Iberian_European": 0.43, "Mediterranean_European": 0.43,
        "Indigenous_Mesoamerican": 0.56, "Indigenous_Andean": 0.56,
        "Indigenous_Caribbean": 0.55,
        "West_African": 0.57, "Central_African": 0.57, "East_African": 0.52,
    }),

    "rs1801133": ("G", "A", 0.6, {
        "Iberian_European": 0.42, "Mediterranean_European": 0.44,
        "Indigenous_Mesoamerican": 0.55, "Indigenous_Andean": 0.55,
        "Indigenous_Caribbean": 0.50,
        "West_African": 0.08, "Central_African": 0.08, "East_African": 0.15,
    }),

    "rs1800629": ("G", "A", 0.5, {
        "Iberian_European": 0.15, "Mediterranean_European": 0.14,
        "Indigenous_Mesoamerican": 0.10, "Indigenous_Andean": 0.10,
        "Indigenous_Caribbean": 0.12,
        "West_African": 0.18, "Central_African": 0.18, "East_African": 0.16,
    }),

    "rs6152": ("A", "G", 0.6, {
        "Iberian_European": 0.45, "Mediterranean_European": 0.45,
        "Indigenous_Mesoamerican": 0.60, "Indigenous_Andean": 0.62,
        "Indigenous_Caribbean": 0.55,
        "West_African": 0.25, "Central_African": 0.25, "East_African": 0.35,
    }),
}


# ---------------------------------------------------------------------------
# REGIONAL REFERENCE PROFILES
# Averaged ancestry proportions from published Hispanic/Latino genetic studies
# (Bryc 2015, Moreno-Estrada 2013/2014, Kehdy 2015, Homburger 2015, Ruiz-Linares 2014).
#
# Each profile is aggregated to 3 continental components for comparison
# against the user's estimated proportions.
# ---------------------------------------------------------------------------

REGIONAL_PROFILES = {
    "Mexican": {
        "European":      0.47,  # ~30-45% (mean ~40-47)
        "NativeAmerican": 0.48,  # ~45-60% (mean ~48)
        "African":        0.05,  # ~5-10%
        "description": (
            "Mexican mestizo populations average roughly half Indigenous "
            "(predominantly Mesoamerican -- Nahua, Maya, Zapotec) and half "
            "European (predominantly Iberian), with a minor African contribution "
            "(~5%). Regional variation is striking: northern Mexican populations "
            "average more European (~55-60%), while southern states like "
            "Oaxaca and Chiapas average more Indigenous (~60-70%). The African "
            "component is concentrated on the Gulf and Pacific coasts "
            "(Veracruz, Guerrero, Oaxaca). Source: Moreno-Estrada et al. 2014."
        ),
    },
    "Puerto_Rican": {
        "European":      0.65,  # ~60-70%
        "NativeAmerican": 0.13,  # ~10-15% (Taino)
        "African":        0.22,  # ~15-25%
        "description": (
            "Puerto Ricans average roughly two-thirds European (predominantly "
            "Iberian, with Canarian contributions), ~20-25% West/Central "
            "African, and ~10-15% Indigenous Taino. Puerto Rico has the "
            "highest surviving Taino genetic signal of any Caribbean nation "
            "(Moreno-Estrada et al. 2013). mtDNA haplogroups show that over "
            "60% of Puerto Ricans carry Amerindian maternal lineages despite "
            "only ~13% autosomal Indigenous ancestry -- reflecting the "
            "demographic history of Spanish male colonists marrying Taino women."
        ),
    },
    "Cuban": {
        "European":      0.70,  # ~60-75%
        "NativeAmerican": 0.05,  # <5%
        "African":        0.25,  # ~20-30%
        "description": (
            "Cuban populations average ~70% European (primarily Iberian, with "
            "substantial Canarian and Galician contributions), ~25% West/Central "
            "African (larger slave-trade flow than Puerto Rico), and a smaller "
            "Indigenous Taino component (~5%) because Cuba's pre-contact "
            "Indigenous population was mostly assimilated or lost earlier. "
            "Eastern Cuba has higher African proportions than western Cuba."
        ),
    },
    "Dominican": {
        "European":      0.52,  # ~45-60%
        "NativeAmerican": 0.08,  # ~5-10% (Taino)
        "African":        0.40,  # ~35-50%
        "description": (
            "Dominicans show the highest African ancestry among the Greater "
            "Antilles (~40%), reflecting Hispaniola's central role in the "
            "Atlantic slave trade. They carry ~8% Taino ancestry and ~52% "
            "European (primarily Iberian with some Canarian). Dominican "
            "ancestry is the closest Caribbean Latino analogue to Brazilian "
            "ancestry in terms of tri-continental balance."
        ),
    },
    "Colombian_Andean": {
        "European":      0.65,  # ~60-75%
        "NativeAmerican": 0.25,  # ~15-25%
        "African":        0.10,  # <10%
        "description": (
            "Andean Colombian populations (Bogota, Antioquia, Boyaca) average "
            "~65% European (Iberian, plus some Basque contribution to the "
            "Antioqueno 'paisa' founder population), ~25% Indigenous (Andean "
            "and northern South American sources), and ~10% African. Pacific "
            "coastal Colombia has much higher African ancestry (~75%+), while "
            "Caribbean coastal Colombia sits between Andean and coastal "
            "profiles. Source: Homburger et al. 2015, Ruiz-Linares CANDELA."
        ),
    },
    "Brazilian": {
        "European":      0.62,  # ~60-75%
        "NativeAmerican": 0.12,  # ~5-15%
        "African":        0.26,  # ~15-30%
        "description": (
            "Brazilian populations are among the most variable in Latin "
            "America. National averages from EPIGEN (Kehdy et al. 2015) show "
            "~62% European (Portuguese plus 20th-century Italian, German, "
            "Spanish immigration), ~26% African (the largest slave-trade "
            "destination in the Americas, dominated by West African and Bantu "
            "sources), and ~12% Indigenous. Regional variation is dramatic: "
            "Bahia averages ~50% African, while southern Brazil (Rio Grande "
            "do Sul) averages ~80% European."
        ),
    },
    "Argentine": {
        "European":      0.85,  # ~80-90%
        "NativeAmerican": 0.12,  # ~5-15%
        "African":        0.03,  # <5%
        "description": (
            "Argentine populations have the highest European ancestry among "
            "Latin American nations, averaging ~85% European (dominated by "
            "late-19th/early-20th-century Italian and Spanish immigration -- "
            "over 5 million Europeans arrived 1880-1950). Indigenous ancestry "
            "(~12%) is mostly Mapuche and Guarani, concentrated in the "
            "northwest and Patagonia. African ancestry is very low (~3%) "
            "because Argentina's colonial-era African population was small "
            "and heavily assimilated by the 19th century."
        ),
    },
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
    Run weighted NNLS decomposition into tri-continental Latino ancestry components.

    Returns (proportions_dict, used_snps_list, residual).
    """
    n_pops = len(LATINO_POPS)
    rows_A = []
    rows_f = []
    rows_w = []
    used_snps = []

    for rsid, (ref, alt, weight, pop_freqs) in AIMS_LATINO.items():
        variant_data = variants.get(rsid.lower())
        if variant_data is None:
            continue
        _chrom, _pos, genotype = variant_data
        dosage = _get_dosage(genotype, alt)
        if dosage is None:
            continue
        obs_freq = dosage / 2.0
        row = [pop_freqs.get(p, 0.0) for p in LATINO_POPS]
        rows_A.append(row)
        rows_f.append(obs_freq)
        rows_w.append(weight)
        used_snps.append(rsid)

    if len(used_snps) < 5:
        return {p: 1.0 / n_pops for p in LATINO_POPS}, used_snps, 0.0

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

    props = {p: float(x[i]) for i, p in enumerate(LATINO_POPS)}
    return props, used_snps, float(resid)


# ---------------------------------------------------------------------------
# AFFINITY SCORES (min-max normalized distances)
# ---------------------------------------------------------------------------

def _compute_affinity_scores(
    variants: Dict[str, Tuple[str, str, str]],
) -> Dict[str, float]:
    """
    Compute relative affinity to each tri-continental source population.
    """
    raw_distances: Dict[str, float] = {}

    for pop in LATINO_POPS:
        weighted_sq_diff = 0.0
        w_sum = 0.0
        for rsid, (ref, alt, weight, pop_freqs) in AIMS_LATINO.items():
            variant_data = variants.get(rsid.lower())
            if variant_data is None:
                continue
            _chrom, _pos, genotype = variant_data
            dosage = _get_dosage(genotype, alt)
            if dosage is None:
                continue
            obs = dosage / 2.0
            ref_freq = pop_freqs.get(pop, 0.0)
            weighted_sq_diff += weight * (obs - ref_freq) ** 2
            w_sum += weight
        if w_sum > 0:
            raw_distances[pop] = np.sqrt(weighted_sq_diff / w_sum)

    if not raw_distances:
        return {p: 0.0 for p in LATINO_POPS}

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
# CONTINENTAL AGGREGATION
# ---------------------------------------------------------------------------

def _aggregate_continental(props: Dict[str, float]) -> Dict[str, float]:
    """Aggregate 8 source populations into 3 continental totals."""
    european = (
        props.get("Iberian_European", 0.0)
        + props.get("Mediterranean_European", 0.0)
    )
    native = (
        props.get("Indigenous_Mesoamerican", 0.0)
        + props.get("Indigenous_Andean", 0.0)
        + props.get("Indigenous_Caribbean", 0.0)
    )
    african = (
        props.get("West_African", 0.0)
        + props.get("Central_African", 0.0)
        + props.get("East_African", 0.0)
    )
    total = european + native + african
    if total <= 0:
        return {"European": 1/3, "NativeAmerican": 1/3, "African": 1/3}
    return {
        "European":       european / total,
        "NativeAmerican": native / total,
        "African":        african / total,
    }


# ---------------------------------------------------------------------------
# REGIONAL PROFILE MATCHING
# ---------------------------------------------------------------------------

def _compute_regional_matches(
    continental: Dict[str, float],
) -> List[Dict[str, Any]]:
    """
    Rank the 7 Hispanic/Latino regional profiles by match score against the
    user's estimated continental proportions. Match score is 100 - (L1 * 100),
    clamped to [0, 100].
    """
    matches: List[Tuple[str, float]] = []
    for region, profile in REGIONAL_PROFILES.items():
        l1 = (
            abs(continental["European"]      - profile["European"])
            + abs(continental["NativeAmerican"] - profile["NativeAmerican"])
            + abs(continental["African"]        - profile["African"])
        )
        # L1 for 3-component proportions ranges 0..2 (perfect..opposite).
        # Normalize to a 0..100 match score.
        score = max(0.0, 100.0 - (l1 / 2.0) * 100.0)
        matches.append((region, round(score, 1)))

    matches.sort(key=lambda t: -t[1])

    result: List[Dict[str, Any]] = []
    for region, score in matches:
        profile = REGIONAL_PROFILES[region]
        result.append({
            "region": region,
            "label": region.replace("_", " "),
            "matchScore": score,
            "description": profile["description"],
            "profile": {
                "European":       round(profile["European"], 3),
                "NativeAmerican": round(profile["NativeAmerican"], 3),
                "African":        round(profile["African"], 3),
            },
        })
    return result


# ---------------------------------------------------------------------------
# KEY LATINO MARKERS SUMMARY
# ---------------------------------------------------------------------------

def _summarize_key_markers(
    variants: Dict[str, Tuple[str, str, str]],
) -> Dict[str, Dict[str, str]]:
    """Report key markers relevant to Latino tri-continental ancestry."""
    results: Dict[str, Dict[str, str]] = {}

    # Duffy rs2814778 (African signal)
    variant_data = variants.get("rs2814778")
    if variant_data is not None:
        _c, _p, geno = variant_data
        if geno and geno not in ("--", "00", "NN"):
            c_count = geno.upper().count("C")
            if c_count == 2:
                results["duffy_african"] = {
                    "status": "CC",
                    "detail": "CC -- Duffy-null (both copies). Strong West African maternal ancestry signal.",
                }
            elif c_count == 1:
                results["duffy_african"] = {
                    "status": "CT",
                    "detail": "CT -- one Duffy-null allele. Partial African ancestry signal.",
                }
            else:
                results["duffy_african"] = {
                    "status": "TT",
                    "detail": "TT -- Duffy expressed (no African Duffy-null signal).",
                }

    # EDAR shovel teeth (rs3827760, Native signal)
    variant_data = variants.get("rs3827760")
    if variant_data is not None:
        _c, _p, geno = variant_data
        if geno and geno not in ("--", "00", "NN"):
            g_count = geno.upper().count("G")
            if g_count == 2:
                results["edar_native"] = {
                    "status": "GG",
                    "detail": "GG -- fully derived EDAR 370A; strong Indigenous American / East Asian signal.",
                }
            elif g_count == 1:
                results["edar_native"] = {
                    "status": "AG",
                    "detail": "AG -- one Indigenous EDAR allele (mestizo or admixed pattern).",
                }
            else:
                results["edar_native"] = {
                    "status": "AA",
                    "detail": "AA -- ancestral; European/African pattern (no Indigenous EDAR signal).",
                }

    # SLC24A5 (rs1426654, European signal)
    variant_data = variants.get("rs1426654")
    if variant_data is not None:
        _c, _p, geno = variant_data
        if geno and geno not in ("--", "00", "NN"):
            a_count = geno.upper().count("A")
            if a_count == 2:
                results["slc24a5_european"] = {
                    "status": "AA",
                    "detail": "AA -- two copies of the European light-skin allele. High European ancestry signal.",
                }
            elif a_count == 1:
                results["slc24a5_european"] = {
                    "status": "AG",
                    "detail": "AG -- one European light-skin allele (admixed).",
                }
            else:
                results["slc24a5_european"] = {
                    "status": "GG",
                    "detail": "GG -- ancestral; low European ancestry signal at this locus.",
                }

    # LCT lactase persistence (rs4988235, European agriculture signal)
    variant_data = variants.get("rs4988235")
    if variant_data is not None:
        _c, _p, geno = variant_data
        if geno and geno not in ("--", "00", "NN"):
            a_count = geno.upper().count("A")
            if a_count > 0:
                status = "AG" if a_count == 1 else "AA"
                results["lactase_european"] = {
                    "status": status,
                    "detail": "A allele present -- European lactase persistence (Iberian / North European contribution).",
                }
            else:
                results["lactase_european"] = {
                    "status": "GG",
                    "detail": "GG -- ancestral; lactase non-persistent (typical of Native American / Sub-Saharan ancestry).",
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
                    "detail": "T allele present -- 'Amerindian' ABCA1 R230C variant; specific to Indigenous American ancestry (highest in Mesoamerica).",
                }
            else:
                results["abca1_amerindian"] = {
                    "status": "CC",
                    "detail": "CC -- no Amerindian ABCA1 signal.",
                }

    # HBB sickle (rs334, African signal)
    variant_data = variants.get("rs334")
    if variant_data is not None:
        _c, _p, geno = variant_data
        if geno and geno not in ("--", "00", "NN"):
            t_count = geno.upper().count("T")
            if t_count > 0:
                status = "AT" if t_count == 1 else "TT"
                results["hbb_sickle"] = {
                    "status": status,
                    "detail": "T allele present -- HbS sickle-cell allele; West/Central African ancestry signal.",
                }

    return results


# ---------------------------------------------------------------------------
# PUBLIC API
# ---------------------------------------------------------------------------

def analyze_latino_breakdown(
    variants: Dict[str, Tuple[str, str, str]],
) -> Dict[str, Any]:
    """
    Analyze variants for Latino / Hispanic tri-continental ancestry breakdown.

    Args:
        variants: Dict mapping rsID (lowercase) -> (chromosome, position, genotype)

    Returns:
        Dict with proportions, affinity_scores, continental aggregates,
        regional matches, used_snps, residual, and key_markers.
    """
    props, used_snps, resid = _estimate_proportions(variants)
    affinity_scores = _compute_affinity_scores(variants)
    continental = _aggregate_continental(props)
    regional_matches = _compute_regional_matches(continental)
    key_markers = _summarize_key_markers(variants)

    return {
        "proportions": props,
        "continental": continental,
        "regional_matches": regional_matches,
        "affinity_scores": affinity_scores,
        "used_snps": used_snps,
        "residual": resid,
        "key_markers": key_markers,
    }


def generate_latino_breakdown_json(result: Dict[str, Any]) -> Dict[str, Any]:
    """
    Convert raw analysis result into the structured JSON for the frontend.

    Args:
        result: Output from analyze_latino_breakdown()

    Returns:
        Structured JSON dict matching the report schema, with tri-continental
        breakdown and ranked regional matches.
    """
    props = result["proportions"]
    continental = result["continental"]
    regional_matches = result["regional_matches"]
    affinity_scores = result["affinity_scores"]
    used_snps = result["used_snps"]
    resid = result["residual"]
    key_markers = result["key_markers"]

    panel_size = len(AIMS_LATINO)
    snps_used = len(used_snps)

    sorted_pops = sorted(LATINO_POPS, key=lambda p: -props.get(p, 0.0))

    populations = []
    for pop in sorted_pops:
        proportion = round(props.get(pop, 0.0), 4)
        populations.append({
            "code": pop,
            "label": LATINO_POP_LABELS[pop],
            "affinityScore": affinity_scores.get(pop, 0.0),
            "proportion": proportion,
            "percentage": round(proportion * 100, 1),
            "description": LATINO_POP_DESCRIPTIONS[pop],
            "color": LATINO_POP_COLORS[pop],
        })

    top_pop = sorted_pops[0]

    # Build regional matches list matching the spec schema exactly
    # ({region, matchScore, description}).
    regional_matches_out = []
    for rm in regional_matches:
        regional_matches_out.append({
            "region":       rm["region"],
            "label":        rm["label"],
            "matchScore":   rm["matchScore"],
            "description":  rm["description"],
            "profile":      rm["profile"],
        })

    return {
        "summary": {
            "panelSize": panel_size,
            "snpsUsed": snps_used,
            "coverage": round(snps_used / panel_size * 100, 1) if panel_size > 0 else 0.0,
        },
        "continental": {
            "European":       round(continental["European"] * 100, 1),
            "NativeAmerican": round(continental["NativeAmerican"] * 100, 1),
            "African":        round(continental["African"] * 100, 1),
        },
        "populations": populations,
        "topMatch": {
            "code":        top_pop,
            "label":       LATINO_POP_LABELS[top_pop],
            "description": LATINO_POP_DESCRIPTIONS[top_pop],
        },
        "regionalMatches": regional_matches_out,
        "keyMarkers": key_markers,
    }
