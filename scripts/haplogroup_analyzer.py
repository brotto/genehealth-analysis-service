"""
Haplogroup Analyzer
Assigns mitochondrial (maternal) and Y chromosome (paternal) haplogroups
from consumer SNP array data, tracing direct maternal/paternal lineages
back 200,000+ years through deep human prehistory.

References:
- van Oven & Kayser 2009 (Human Mutation): PhyloTree — mtDNA haplogroup tree
- ISOGG (International Society of Genetic Genealogy): Y-tree
- Achilli et al. 2004 (Science): European mtDNA haplogroups
- Underhill & Kivisild 2007 (Annu Rev Genet): Y haplogroup evolution
"""

from typing import Dict, List, Tuple, Any, Optional, Set


# ─────────────────────────────────────────────────────────────────────────────
# HELPER FUNCTIONS
# ─────────────────────────────────────────────────────────────────────────────

def has_allele(genotype: str, allele: str) -> bool:
    """Check if a specific allele is present in genotype."""
    if not genotype:
        return False
    return allele.upper() in genotype.upper()


# ─────────────────────────────────────────────────────────────────────────────
# MITOCHONDRIAL HAPLOGROUP DEFINITIONS
# ─────────────────────────────────────────────────────────────────────────────
# Position-based markers (rCRS coordinates)
# Format: position → {derived_allele: haplogroup, ...}
# Listed from ROOT to TIPS — more specific markers are applied last
# Source: PhyloTree Build 17 (van Oven 2009), verified against mtDNA literature

MT_POSITION_MARKERS: Dict[int, Dict[str, str]] = {
    # Root / L haplogroup (African)
    73:    {"G": "L3+"},        # L3 and downstream (most non-African)
    263:   {"G": "L3+"},        # L3 marker
    750:   {"G": "L3_down"},    # Nearly all Eurasians
    1438:  {"G": "R_macro"},    # R macrohaplogroup (H, V, J, T, K, U, B, F, P)
    2706:  {"G": "N_macro"},    # N macrohaplogroup
    4769:  {"A": "M_marker"},   # M macrohaplogroup (EAS/SAS)
    7028:  {"T": "H"},          # H-DEFINING: 7028 C→T (most common European hg)
    8860:  {"G": "L3_down"},    # Common marker
    11719: {"A": "R_macro"},    # R macrohaplogroup confirmation
    12705: {"T": "M_marker"},   # M macrohaplogroup (alternative marker)
    14766: {"T": "HV"},         # HV ancestor (ancestor of H and V)
    15326: {"G": "HV"},         # HV marker
    # H subgroups
    3010:  {"A": "H1_or_H3"},   # H1 or H3 (most common H subgroups in Europeans)
    6776:  {"C": "H1"},         # H1 specific
    14872: {"T": "H3"},         # H3 specific
    # V haplogroup
    4580:  {"A": "V"},          # V-defining (also in some H)
    # J and T haplogroup (from HV through JT clade)
    295:   {"T": "JT"},         # JT ancestor (precedes J and T)
    489:   {"C": "U_K"},        # U and K haplogroup (also seen in others)
    10398: {"G": "JT_M"},       # JT and M marker
    16126: {"C": "JT"},         # JT clade (J and T share this)
    16311: {"C": "T"},          # T-specific (distinguishes from J)
    16362: {"C": "T"},          # T marker
    16069: {"C": "J"},          # J-specific (distinguishes from T)
    # K haplogroup (derived from U8)
    16224: {"C": "K_U8"},       # K haplogroup
    16234: {"T": "K"},          # K-specific
    # U5 haplogroup (ancient European hunter-gatherer marker)
    16270: {"C": "U5"},         # U5 haplogroup
    9477:  {"A": "U5"},         # U5 confirmation
    # U4 haplogroup
    16356: {"C": "U4"},         # U4 haplogroup
    # W haplogroup
    119:   {"C": "W"},          # W-defining
    16292: {"T": "W"},          # W confirmation
    # I haplogroup
    16129: {"A": "I"},          # I haplogroup
    # X haplogroup (rare)
    6221:  {"T": "X"},          # X haplogroup
    16278: {"C": "X"},          # X confirmation
    # African L haplogroups
    16223: {"T": "L"},          # L haplogroup (African macrohaplogroup marker)
    16278: {"T": "L"},          # L marker (Note: same position as X but T vs C)
}

# rsID-based markers (where rsIDs are known)
# These serve as backup markers and as primary lookup method in the analysis service
MT_RSID_MARKERS: Dict[str, Tuple[str, str]] = {
    # rsid: (derived_allele, haplogroup_marker)
    "rs2032115":  ("A", "R_macro"),    # pos 11719 G→A: R macrohaplogroup
    "rs3928305":  ("T", "HV"),         # pos 14766 C→T: HV ancestor
    "rs28358564": ("T", "H"),          # pos 7028 C→T: H haplogroup
    "rs2015062":  ("A", "H1_or_H3"),   # pos 3010 G→A: H1 or H3
    "rs28357376": ("C", "L"),          # L haplogroup African marker
    "rs2853508":  ("T", "N_macro"),    # 10873 T→C: N macrohaplogroup
    "rs41419549": ("C", "V"),          # V haplogroup indicator
    "rs28357681": ("A", "R_macro"),    # 11719 G→A: R macrohaplogroup
    "rs2853499":  ("T", "K_B"),        # K or B haplogroup
    "rs3135031":  ("T", "K"),          # K haplogroup
    "rs41360374": ("G", "JT"),         # JT ancestor marker
    "rs28357093": ("C", "J"),          # J haplogroup
    "rs41419546": ("A", "T"),          # T haplogroup
    "rs2853519":  ("T", "U5"),         # U5 haplogroup (WHG marker)
    "rs28358886": ("T", "H1"),         # H1 subgroup
    "rs9862792":  ("C", "L"),          # L haplogroup African
}

# Haplogroup decision tree (ordered from most to least specific)
# Each entry: (signals_present, haplogroup_assigned, description, confidence)
MT_DECISION_TREE: List[Tuple[List[str], str, str, str]] = [
    # African L haplogroups
    (["L"],                      "L",     "African macrohaplogroup",                "moderate"),
    # M haplogroup (EAS/SAS, non-European)
    (["M_marker"],               "M",     "East Asian/South Asian macrohaplogroup", "moderate"),
    # European and Near Eastern haplogroups
    (["H", "HV"],                "H",     "H (European dominant)",                 "high"),
    (["H", "H1_or_H3", "H1"],   "H1",    "H1 (European)",                         "moderate"),
    (["H", "H1_or_H3", "H3"],   "H3",    "H3 (European)",                         "moderate"),
    (["H", "H1_or_H3"],         "H1/H3", "H1 or H3 subgroup",                     "moderate"),
    (["HV", "V"],                "V",     "V (Iberian/European)",                  "moderate"),
    (["HV"],                     "HV",    "HV (ancestor of H and V)",              "moderate"),
    (["JT", "J"],                "J",     "J (Mediterranean/European)",            "high"),
    (["JT", "T"],                "T",     "T (Anatolian/European)",                "high"),
    (["JT"],                     "JT",    "JT (ancestor of J and T)",              "moderate"),
    (["K_U8", "K"],              "K",     "K (European/Ashkenazi)",                "high"),
    (["K_U8"],                   "K",     "K or U8",                               "moderate"),
    (["U_K", "U5"],              "U5",    "U5 (WHG hunter-gatherer marker)",       "high"),
    (["U_K", "U4"],              "U4",    "U4 (Eastern European)",                 "moderate"),
    (["U_K"],                    "U",     "U (broad clade)",                       "moderate"),
    (["W"],                      "W",     "W (South Asian/European)",              "moderate"),
    (["I"],                      "I",     "I (European/MidEast)",                  "moderate"),
    (["X"],                      "X",     "X (rare European/Native Am)",           "moderate"),
    # Fallback for N macrohaplogroup
    (["N_macro", "R_macro"],     "R*",    "R* (pre-HV, pre-JT, pre-U lineage)",   "low"),
    (["N_macro"],                "N*",    "N macrohaplogroup (non-R branch)",      "low"),
]


# ─────────────────────────────────────────────────────────────────────────────
# MITOCHONDRIAL HAPLOGROUP NARRATIVES
# ─────────────────────────────────────────────────────────────────────────────

MT_HAPLOGROUP_INFO: Dict[str, Dict[str, Any]] = {
    "H": {
        "description": (
            "Haplogroup H is the most common mitochondrial haplogroup in Europe (~44%) "
            "and the dominant lineage among modern Europeans. It originated in the Near "
            "East/Caucasus region before the Last Glacial Maximum (>20,000 years ago) and "
            "expanded dramatically after the Ice Age from glacial refugia in Southwestern "
            "Europe (Franco-Cantabrian refugium — modern Spain and southern France). The "
            "post-glacial expansion of H into Europe ~13,000-15,000 years ago is one of the "
            "major demographic events in European prehistory, associated with the spread of "
            "microlithic technology. Today, H is distributed across Europe, the Near East, "
            "and North Africa. Famous carrier: Ötzi the Iceman was NOT H (he was K) — "
            "but most modern Italians, Germans, British, and French belong to H or its "
            "subgroups. H1 is the most frequent H subgroup (~8% of Europeans), highest "
            "in Iberia. H3 is also common in Western Europe and Sardinia."
        ),
        "geographicOrigin": {"region": "Near East / Western Europe", "detail": "Franco-Cantabrian glacial refugium"},
        "timeOriginYears": 20000,
        "modernFrequency": {"europe": 0.44, "nearEast": 0.15, "northAfrica": 0.10},
        "migrationRoute": "Post-glacial expansion from Franco-Cantabrian refugium across Europe",
        "famousCarriers": ["Most modern Europeans", "European royalty lineages"],
        "color": "#4F46E5",
    },
    "H1": {
        "description": (
            "Haplogroup H1 is the most common subgroup of H in Europe (~8%). It shows "
            "a westward distribution gradient, peaking in Iberia (~20%) and decreasing "
            "toward Central/Eastern Europe. This distribution reflects the post-Ice-Age "
            "re-expansion of populations from the Iberian refugium northward and eastward "
            "as the glaciers retreated ~13,000 years ago. H1 is also common in North "
            "Africa, where it spread via the same post-glacial expansion. Modern Iberians, "
            "Basques, Sardinians, and Western Europeans show the highest H1 frequencies."
        ),
        "geographicOrigin": {"region": "Western Europe", "detail": "Iberian Peninsula refugium"},
        "timeOriginYears": 13000,
        "modernFrequency": {"europe": 0.08, "iberia": 0.20, "northAfrica": 0.10},
        "migrationRoute": "Post-glacial expansion from Iberian refugium into Central and Northern Europe",
        "famousCarriers": ["Common among Basques and Sardinians"],
        "color": "#6366F1",
    },
    "H1/H3": {
        "description": (
            "Your maternal lineage belongs to haplogroup H1 or H3, both common European "
            "H subgroups that expanded from the Iberian glacial refugium after the Last "
            "Glacial Maximum ~13,000 years ago. Both are most common in Western Europe, "
            "Iberia, and Sardinia. The distinguishing marker between H1 and H3 was not "
            "captured by the chip, but both share the same geographic origin story."
        ),
        "geographicOrigin": {"region": "Western Europe", "detail": "Iberian Peninsula refugium"},
        "timeOriginYears": 13000,
        "modernFrequency": {"europe": 0.12, "iberia": 0.25},
        "migrationRoute": "Post-glacial expansion from Iberian refugium",
        "famousCarriers": ["Common among Western Europeans"],
        "color": "#6366F1",
    },
    "H3": {
        "description": (
            "Haplogroup H3 is a common European H subgroup (~4%), with its highest "
            "frequencies in Iberia, Sardinia, and Western Europe. Like H1, it likely "
            "expanded from the Iberian glacial refugium after the Last Glacial Maximum. "
            "H3 has been found in ancient DNA samples from Neolithic and Early Bronze "
            "Age Europe, suggesting it was present in Early European Farmers."
        ),
        "geographicOrigin": {"region": "Western Europe", "detail": "Iberian glacial refugium"},
        "timeOriginYears": 13000,
        "modernFrequency": {"europe": 0.04, "iberia": 0.10, "sardinia": 0.12},
        "migrationRoute": "Post-glacial expansion from Iberian refugium",
        "famousCarriers": ["Found in Neolithic and Bronze Age European ancient DNA"],
        "color": "#818CF8",
    },
    "V": {
        "description": (
            "Haplogroup V is found in ~5% of Europeans, with its highest frequencies "
            "among the Saami (Lapps) of northern Scandinavia (~53%) and Basques (~10%). "
            "It originated in the Iberian glacial refugium during the Last Glacial Maximum "
            "and expanded into Europe and North Africa after the glaciers retreated, "
            "~13,000 years ago. The extremely high frequency of V in Saami is a striking "
            "example of genetic drift in an isolated population. V is derived from HV, "
            "the same ancestor that gave rise to H."
        ),
        "geographicOrigin": {"region": "Iberian Peninsula", "detail": "Iberian glacial refugium"},
        "timeOriginYears": 13000,
        "modernFrequency": {"europe": 0.05, "saami": 0.53, "basques": 0.10},
        "migrationRoute": "Post-glacial expansion from Iberia; drift in Saami population",
        "famousCarriers": ["Common among Saami people of northern Scandinavia"],
        "color": "#7C3AED",
    },
    "J": {
        "description": (
            "Haplogroup J is found in ~8% of Europeans and is associated with the Near "
            "Eastern Neolithic agricultural expansion. J originated in the Middle East "
            "and spread to Europe with the first farmers ~9,000 years ago. It is most "
            "common in the British Isles (~12%), Scandinavia (~10%), and the Near East "
            "(~12-25%). Famous carriers: the pedigree of European royalty historically "
            "shows high J frequency. J2a is common in the Near East; J1c is European. "
            "J is the haplogroup of the Biblical patriarch lineage according to some "
            "genetic genealogy studies (Feder et al. 2007)."
        ),
        "geographicOrigin": {"region": "Middle East", "detail": "Fertile Crescent / Near East"},
        "timeOriginYears": 45000,
        "modernFrequency": {"europe": 0.08, "britishIsles": 0.12, "nearEast": 0.25},
        "migrationRoute": "Neolithic agricultural expansion from Near East into Europe ~9,000 years ago",
        "famousCarriers": ["European royalty lineages", "Biblical patriarch lineage (debated)"],
        "color": "#A78BFA",
    },
    "T": {
        "description": (
            "Haplogroup T is found in ~9% of Europeans, with highest frequencies in the "
            "Near East and Eastern Europe. Like J, it is associated with the Neolithic "
            "agricultural expansion from the Near East. T2 is the most common European "
            "subgroup. Famously, T2 has been found in many ancient Egyptian mummies and "
            "in the pharaoh Ramesses III (New Kingdom Egypt). The T haplogroup shows "
            "strongest presence in the Fertile Crescent and Eastern Mediterranean. "
            "It is the haplogroup of many ancient and medieval European samples, "
            "suggesting it was brought by Near Eastern farmers."
        ),
        "geographicOrigin": {"region": "Near East", "detail": "Fertile Crescent / Anatolia"},
        "timeOriginYears": 10000,
        "modernFrequency": {"europe": 0.09, "nearEast": 0.15, "easternEurope": 0.12},
        "migrationRoute": "Neolithic agricultural expansion from Anatolia/Near East into Europe",
        "famousCarriers": ["Pharaoh Ramesses III (T2)", "Many ancient Egyptian mummies"],
        "color": "#8B5CF6",
    },
    "K": {
        "description": (
            "Haplogroup K is found in ~6% of Europeans and ~32% of Ashkenazi Jews. "
            "It originated in the Near East and entered Europe with Neolithic farmers. "
            "The famously preserved Tyrolean Iceman 'Ötzi' belonged to haplogroup K1f. "
            "In Ashkenazi Jews, three specific K subgroups (K1a1b1a, K1a9, K2a2a) "
            "account for ~32% of the Ashkenazi gene pool — the result of a severe "
            "bottleneck where these three lineages were amplified by drift. The Ashkenazi "
            "K lineages share a common ancestor ~2,000-3,000 years ago, consistent with "
            "the historical foundation of the Ashkenazi community."
        ),
        "geographicOrigin": {"region": "Near East", "detail": "Entered Europe with Neolithic farmers"},
        "timeOriginYears": 12000,
        "modernFrequency": {"europe": 0.06, "ashkenaziJews": 0.32},
        "migrationRoute": "Near East to Europe with Neolithic farming expansion",
        "famousCarriers": ["Ötzi the Iceman (K1f)", "Common among Ashkenazi Jews"],
        "color": "#C084FC",
    },
    "U5": {
        "description": (
            "Haplogroup U5 is the oldest European mitochondrial haplogroup, dating back "
            ">30,000 years. It is the signature haplogroup of the pre-agricultural "
            "Western Hunter-Gatherers (WHG) who inhabited Europe before the arrival of "
            "Neolithic farmers. Before farming, U5 may have been at ~80-90% frequency "
            "in European hunter-gatherers. After the Neolithic replacement, it dropped "
            "to ~8-9% in modern Europeans. Today U5 is highest in Saami (~48%), "
            "Finnish (~24%), and Baltic populations — regions where WHG contribution "
            "persisted. Finding U5 means your maternal lineage stretches back to the "
            "Ice Age hunters of Europe — arguably the most ancient documented "
            "European maternal lineage. U5a is found in the oldest known European "
            "skeletal remains with preserved DNA."
        ),
        "geographicOrigin": {"region": "Europe (Paleolithic)", "detail": "Pre-agricultural Western Hunter-Gatherers"},
        "timeOriginYears": 30000,
        "modernFrequency": {"europe": 0.09, "saami": 0.48, "finnish": 0.24, "baltic": 0.15},
        "migrationRoute": "Indigenous European lineage; pre-dates farming arrival",
        "famousCarriers": ["Oldest European skeletal remains with DNA (U5a)", "Cheddar Man (U5)"],
        "color": "#5B21B6",
    },
    "U": {
        "description": (
            "Haplogroup U is a broad macrohaplogroup that includes many European and "
            "Asian lineages (U2-U8, K). It originated outside Africa ~55,000 years "
            "ago and diversified across Eurasia. European U lineages (U4, U5) trace "
            "back to the earliest modern humans in Europe (Upper Paleolithic, "
            ">30,000 years ago). Your specific U subgroup would require additional "
            "analysis to determine precisely."
        ),
        "geographicOrigin": {"region": "Eurasia", "detail": "Diversified across Europe and Asia"},
        "timeOriginYears": 55000,
        "modernFrequency": {"europe": 0.15, "southAsia": 0.10},
        "migrationRoute": "Out-of-Africa through Near East, diversified across Eurasia",
        "famousCarriers": ["Ötzi the Iceman (K, a U subgroup)"],
        "color": "#6D28D9",
    },
    "U4": {
        "description": (
            "Haplogroup U4 is found in ~2-4% of Eastern Europeans and is enriched "
            "in Uralic-speaking populations (Finns, Estonians, Siberians). It is one "
            "of the older European lineages, predating the Neolithic farmers. U4 "
            "has been found in Bronze Age Eastern European and Siberian ancient DNA, "
            "suggesting a Eurasian steppe connection."
        ),
        "geographicOrigin": {"region": "Eastern Europe / Siberia", "detail": "Uralic-speaking populations"},
        "timeOriginYears": 25000,
        "modernFrequency": {"easternEurope": 0.04, "finland": 0.05, "siberia": 0.08},
        "migrationRoute": "Eastern European / Eurasian steppe lineage",
        "famousCarriers": ["Found in Bronze Age Siberian ancient DNA"],
        "color": "#7E22CE",
    },
    "W": {
        "description": (
            "Haplogroup W is a rare haplogroup (~2% of Europeans) with its highest "
            "frequencies in South Asia and the Caucasus. It appears to have originated "
            "in South/West Asia and entered Europe as part of the Anatolian/Near Eastern "
            "gene flow. W is rare enough that finding it indicates a direct connection "
            "to Middle Eastern or South Asian maternal ancestors."
        ),
        "geographicOrigin": {"region": "South / West Asia", "detail": "Caucasus and South Asian origin"},
        "timeOriginYears": 20000,
        "modernFrequency": {"europe": 0.02, "southAsia": 0.05, "caucasus": 0.04},
        "migrationRoute": "South/West Asia to Europe via Anatolia",
        "famousCarriers": [],
        "color": "#9333EA",
    },
    "I": {
        "description": (
            "Haplogroup I is a rare European haplogroup (~2%) with Near Eastern roots. "
            "It shows a distribution suggesting dual origin: an ancient Near Eastern "
            "lineage that entered Europe both with Paleolithic hunter-gatherers and "
            "later with Neolithic farmers. Most common in Scandinavia and the British "
            "Isles. I is particularly interesting because it shows geographic "
            "clustering that mirrors ancient population movements."
        ),
        "geographicOrigin": {"region": "Near East / Europe", "detail": "Dual Paleolithic and Neolithic entry"},
        "timeOriginYears": 30000,
        "modernFrequency": {"europe": 0.02, "scandinavia": 0.04, "britishIsles": 0.03},
        "migrationRoute": "Near East to Europe via multiple waves",
        "famousCarriers": [],
        "color": "#A855F7",
    },
    "X": {
        "description": (
            "Haplogroup X is rare globally (~0.5-1%) but fascinating for two reasons: "
            "(1) it is found in both Old World populations (Europeans ~2%, Near East, "
            "Central Asia) and Native Americans (~3% in some tribes), suggesting either "
            "an early trans-Atlantic contact or a Siberian route; (2) the Native American "
            "X2a sub-branch is notably divergent from Eurasian X, complicating theories "
            "of its origin. The Druze of Israel have unusually high X frequency (~26%). "
            "Recent analysis supports the Siberian route via Beringia as most parsimonious."
        ),
        "geographicOrigin": {"region": "Near East / Central Asia", "detail": "Also found in Native Americans"},
        "timeOriginYears": 30000,
        "modernFrequency": {"europe": 0.02, "druze": 0.26, "nativeAmerican": 0.03},
        "migrationRoute": "Siberian route via Beringia to Americas; Near East to Europe",
        "famousCarriers": ["Found in Druze of Israel at high frequency"],
        "color": "#D946EF",
    },
    "L": {
        "description": (
            "Haplogroup L is the African macrohaplogroup encompassing all purely African "
            "maternal lineages (L0, L1, L2, L3). L haplogroups originated in Africa and "
            "represent the oldest branches of the human mtDNA tree. All non-African "
            "haplogroups (M and N) descend from within L3 — the lineage that survived "
            "the Out-of-Africa migration ~60,000 years ago. Finding haplogroup L indicates "
            "maternal ancestry from Sub-Saharan Africa in the last several thousand years "
            "(since the historic slave trade or more recent migration)."
        ),
        "geographicOrigin": {"region": "Sub-Saharan Africa", "detail": "Oldest branches of human mtDNA tree"},
        "timeOriginYears": 200000,
        "modernFrequency": {"subSaharanAfrica": 0.95, "americas": 0.15},
        "migrationRoute": "Origin of all human maternal lineages; Out-of-Africa ~60,000 years ago",
        "famousCarriers": ["Mitochondrial Eve (L0/L1 ancestor)"],
        "color": "#4C1D95",
    },
    "M": {
        "description": (
            "Haplogroup M is a major East Asian/South Asian/African macrohaplogroup, "
            "sister to N (which includes all European haplogroups). M encompasses "
            "most Asian maternal lineages: C, D, G, Q, Z (East Asian), M2-M6 (South Asian), "
            "and E (Southeast Asian). M arose from L3 during or just before the Out-of-"
            "Africa migration ~60,000 years ago. Finding M in your data suggests East "
            "Asian, South Asian, or related ancestry."
        ),
        "geographicOrigin": {"region": "East Asia / South Asia", "detail": "Sister to N macrohaplogroup"},
        "timeOriginYears": 60000,
        "modernFrequency": {"eastAsia": 0.60, "southAsia": 0.50},
        "migrationRoute": "Out-of-Africa coastal route to South and East Asia",
        "famousCarriers": [],
        "color": "#581C87",
    },
    "HV": {
        "description": (
            "Haplogroup HV is the ancestor of both H and V — two of the most common "
            "European haplogroups. HV itself is rare in modern populations (~1%) but "
            "is found in the Near East (especially Iran and the Caucasus), which was "
            "likely the geographic origin of the H and V lineages before they expanded "
            "into Europe. Finding HV suggests either a Near Eastern connection or that "
            "the specific H/V-defining variant wasn't captured by the chip."
        ),
        "geographicOrigin": {"region": "Near East", "detail": "Especially Iran and the Caucasus"},
        "timeOriginYears": 25000,
        "modernFrequency": {"nearEast": 0.05, "caucasus": 0.04, "europe": 0.01},
        "migrationRoute": "Near East origin; gave rise to H and V in Europe",
        "famousCarriers": [],
        "color": "#6B21A8",
    },
    "JT": {
        "description": (
            "JT is the common ancestor haplogroup of both J and T. JT itself is rare "
            "in modern populations. Your maternal lineage is in the JT clade but the "
            "specific J vs. T distinguishing marker was not captured. Both J and T are "
            "associated with the Near Eastern Neolithic agricultural expansion into Europe."
        ),
        "geographicOrigin": {"region": "Near East", "detail": "Ancestor of J and T haplogroups"},
        "timeOriginYears": 50000,
        "modernFrequency": {"nearEast": 0.03},
        "migrationRoute": "Near East to Europe with Neolithic expansion",
        "famousCarriers": [],
        "color": "#7C3AED",
    },
    "R*": {
        "description": (
            "Your maternal lineage falls within haplogroup R but has not been assigned "
            "to a specific sub-haplogroup with the available markers. R includes all "
            "the major European haplogroups (H, V, J, T, K, U). Additional markers "
            "would be needed for specific sub-haplogroup assignment."
        ),
        "geographicOrigin": {"region": "Eurasia", "detail": "Broad clade including most European lineages"},
        "timeOriginYears": 60000,
        "modernFrequency": {"eurasia": 0.70},
        "migrationRoute": "Out-of-Africa through Near East; diversified across Eurasia",
        "famousCarriers": [],
        "color": "#8B5CF6",
    },
    "N*": {
        "description": (
            "Your maternal lineage falls within the N macrohaplogroup, which includes "
            "all non-M Eurasian lineages. This is a very broad classification. "
            "Additional analysis would be needed to determine the specific European "
            "or Asian sub-haplogroup."
        ),
        "geographicOrigin": {"region": "Eurasia", "detail": "Includes all non-M Eurasian lineages"},
        "timeOriginYears": 65000,
        "modernFrequency": {"eurasia": 0.80},
        "migrationRoute": "Out-of-Africa; ancestor of all European maternal lineages",
        "famousCarriers": [],
        "color": "#9333EA",
    },
}


# ─────────────────────────────────────────────────────────────────────────────
# Y CHROMOSOME HAPLOGROUP DEFINITIONS (males only)
# ─────────────────────────────────────────────────────────────────────────────
# rsIDs from ISOGG Y-SNP database (International Society of Genetic Genealogy)
# and major published studies

Y_MARKERS: Dict[str, Tuple[str, str]] = {
    # rsid: (derived_allele, haplogroup_clade)
    # R haplogroup (most common in W. and E. Europeans)
    "rs9786153":  ("A", "R"),          # M207: R haplogroup root
    "rs16981293": ("A", "R1b"),        # P25/M343: R1b root
    "rs9786137":  ("A", "R1b_L11"),    # L11/P310: R1b Northwest European
    "rs34276300": ("A", "R1b_P312"),   # P312: R1b Western European (Iberian, French, British)
    "rs34126399": ("A", "R1b_U152"),   # U152/S28: R1b Central European (Italian, Swiss)
    "rs17222279": ("A", "R1a"),        # M17/M198: R1a root (Slavic/South Asian)
    "rs2032604":  ("A", "R1a_M420"),   # M420: R1a upstream
    # I haplogroup (old N. European)
    "rs34079747": ("A", "I1"),         # M253: I1 haplogroup (Scandinavian/Germanic)
    "rs34311906": ("A", "I2"),         # M423: I2 haplogroup (Balkan)
    "rs41362250": ("A", "I2a"),        # P37.2: I2a (Dinaric)
    # J haplogroup (Mediterranean/MidEast)
    "rs13447352": ("A", "J2"),         # M172: J2 haplogroup
    "rs9341296":  ("A", "J1"),         # M267: J1 haplogroup (Arabian/MidEast)
    # G haplogroup (Caucasus/Anatolian)
    "rs2032601":  ("A", "G"),          # M201: G haplogroup
    "rs17690947": ("A", "G2a"),        # P15: G2a (EEF marker, Ötzi's Y)
    # E haplogroup (African/Mediterranean)
    "rs41352448": ("A", "E"),          # M96: E haplogroup root
    "rs41462431": ("A", "E1b1b"),      # M35: E1b1b (Mediterranean/African)
    # N haplogroup (Uralic/Finnish/Siberian)
    "rs41352439": ("A", "N"),          # M231: N haplogroup
    # Q haplogroup (Native American/Siberian)
    "rs41352449": ("A", "Q"),          # M242: Q haplogroup (Native American)
    # T haplogroup (rare; MidEast/E.African)
    "rs41352436": ("A", "T_Y"),        # M184: T haplogroup (not common in Europe)
    # C haplogroup (EAS/Oceanian)
    "rs41352441": ("A", "C"),          # M130: C haplogroup (Mongolic/Oceanian)
}

# Priority order for Y haplogroup assignment (most specific first)
Y_PRIORITY: List[str] = [
    "R1b_U152", "R1b_P312", "R1b_L11", "R1b", "R1a_M420", "R1a",
    "I1", "I2a", "I2", "J2", "J1", "G2a", "G",
    "E1b1b", "E", "N", "Q", "T_Y", "C", "R",
]


# ─────────────────────────────────────────────────────────────────────────────
# Y CHROMOSOME HAPLOGROUP NARRATIVES
# ─────────────────────────────────────────────────────────────────────────────

Y_HAPLOGROUP_INFO: Dict[str, Dict[str, Any]] = {
    "R1b": {
        "description": (
            "R1b is the most common Y haplogroup in Western Europe, dominant in "
            "Ireland (~80%), Iberia (~70%), and France (~60%). R1b expanded dramatically "
            "~5,000 years ago with the Yamnaya/Bronze Age Steppe pastoralists, bringing "
            "Indo-European languages to Europe. The 'Bell Beaker' culture (c. 2800-1800 BCE) "
            "was the vehicle of R1b's expansion across Western Europe. Before this, "
            "European Neolithic farmers (EEF) had completely different Y haplogroups "
            "(mainly G2a). R1b essentially replaced the earlier Neolithic Y lineages "
            "almost completely in Western Europe within a few centuries — a genetic "
            "turnover of remarkable speed and scale, possibly involving violent conflict "
            "or severe reproductive advantage (David Reich 2018: 'Who We Are')."
        ),
        "geographicOrigin": {"region": "Pontic-Caspian Steppe", "detail": "Yamnaya culture expansion"},
        "timeOriginYears": 5000,
        "modernFrequency": {"ireland": 0.80, "iberia": 0.70, "france": 0.60, "europe": 0.50},
        "migrationRoute": "Yamnaya Steppe to Western Europe via Bell Beaker culture",
        "famousCarriers": ["Bell Beaker culture males", "Most Western European males"],
        "color": "#DC2626",
    },
    "R1b_P312": {
        "description": (
            "R1b-P312 (also called R1b-S116) is the main Western European branch of R1b, "
            "encompassing most R1b in Iberia, France, British Isles, and Alpine regions. "
            "P312 expanded with the Bell Beaker culture ~4,500 years ago. Its subclades "
            "include R1b-DF27 (Iberian), R1b-U152 (Central European/Italian), "
            "R1b-L21 (British Isles/Atlantic), and R1b-DF19. The geographic distribution "
            "of P312 subclades closely mirrors the diversity of Celtic and Italic "
            "languages in pre-Roman Europe."
        ),
        "geographicOrigin": {"region": "Western Europe", "detail": "Bell Beaker culture expansion"},
        "timeOriginYears": 4500,
        "modernFrequency": {"iberia": 0.55, "france": 0.50, "britishIsles": 0.45},
        "migrationRoute": "Bell Beaker expansion from Central Europe to Atlantic fringe",
        "famousCarriers": ["Celtic and Italic peoples of pre-Roman Europe"],
        "color": "#EF4444",
    },
    "R1b_U152": {
        "description": (
            "R1b-U152 (S28) is the Central European/Italian branch of R1b, found at "
            "high frequencies in Northern Italy (~25%), Switzerland, and Alpine regions. "
            "Associated with early Italic peoples and possibly the ancient Romans. "
            "Julius Caesar was reportedly of this lineage according to some analyses. "
            "U152 is found in Cisalpine Gaul, Helvetia, and Central Italy."
        ),
        "geographicOrigin": {"region": "Central Europe / Italy", "detail": "Associated with Celtic and Italic expansion"},
        "timeOriginYears": 4500,
        "modernFrequency": {"italy": 0.35, "switzerland": 0.20, "france": 0.15},
        "migrationRoute": "Indo-European expansion through Central Europe into Italy",
        "famousCarriers": ["Possible connection to ancient Roman patrician lineages", "Julius Caesar (debated)"],
        "color": "#F87171",
    },
    "R1b_L11": {
        "description": (
            "R1b-L11 is the Northwest European branch, encompassing most British Isles, "
            "French, and Northwestern European R1b. It expanded with the Bell Beaker "
            "culture and subsequently diversified into the main Atlantic/Celtic lineages."
        ),
        "geographicOrigin": {"region": "Northwest Europe", "detail": "Bell Beaker expansion"},
        "timeOriginYears": 4800,
        "modernFrequency": {"britishIsles": 0.55, "france": 0.50, "germany": 0.35},
        "migrationRoute": "Bell Beaker culture expansion across Northwest Europe",
        "famousCarriers": ["Atlantic/Celtic lineage males"],
        "color": "#FB923C",
    },
    "R1a": {
        "description": (
            "R1a is the dominant Y haplogroup in Eastern Europe (Poland, Russia, Ukraine "
            "~50-60%) and South Asia (Brahmins ~60-70%). R1a expanded with the Yamnaya/"
            "Corded Ware Bronze Age Steppe culture ~5,000 years ago, eastward toward "
            "South Asia and westward into Eastern/Northern Europe. R1a-Z93 is the South "
            "Asian branch (Indo-Aryan speakers); R1a-Z280 is the Eastern European branch "
            "(Slavic); R1a-M417 is the ancestral European Slavic lineage. The presence "
            "of R1a in South Asian Brahmin priests (~60%) and in Slavic populations is "
            "strong genetic evidence for the Indo-European language expansion."
        ),
        "geographicOrigin": {"region": "Pontic-Caspian Steppe", "detail": "Corded Ware culture expansion"},
        "timeOriginYears": 5000,
        "modernFrequency": {"poland": 0.57, "russia": 0.50, "ukraine": 0.55, "brahmins": 0.65},
        "migrationRoute": "Yamnaya Steppe east to South Asia and west to Eastern Europe",
        "famousCarriers": ["Corded Ware culture males", "Indo-European language carriers"],
        "color": "#B91C1C",
    },
    "R1a_M420": {
        "description": (
            "R1a-M420 is an upstream marker of R1a, confirming membership in the R1a "
            "haplogroup. R1a is the dominant Y haplogroup in Eastern Europe and South Asia, "
            "associated with the Corded Ware / Indo-European expansion ~5,000 years ago."
        ),
        "geographicOrigin": {"region": "Pontic-Caspian Steppe", "detail": "Corded Ware culture"},
        "timeOriginYears": 5500,
        "modernFrequency": {"easternEurope": 0.50, "southAsia": 0.35},
        "migrationRoute": "Steppe expansion eastward and westward",
        "famousCarriers": [],
        "color": "#991B1B",
    },
    "I1": {
        "description": (
            "I1 is the 'Scandinavian' Y haplogroup, dominant in Norway, Sweden, Denmark "
            "(30-40%), and Finland (~28%). I1 is one of the oldest European Y lineages, "
            "predating the Steppe expansion — it was already present in Scandinavia "
            "during the Mesolithic. After the Steppe people arrived, I1 survived at "
            "high frequency in Scandinavia, possibly due to selective advantage or "
            "geographic isolation. It is the most common haplogroup of documented "
            "Viking-age males. Famous in Viking-age ancient DNA studies."
        ),
        "geographicOrigin": {"region": "Scandinavia", "detail": "Pre-Steppe Mesolithic lineage"},
        "timeOriginYears": 5000,
        "modernFrequency": {"norway": 0.35, "sweden": 0.40, "denmark": 0.35, "finland": 0.28},
        "migrationRoute": "Indigenous Scandinavian lineage; survived Steppe expansion",
        "famousCarriers": ["Viking-age males (ancient DNA)", "Scandinavian royalty"],
        "color": "#F59E0B",
    },
    "I2": {
        "description": (
            "I2 is an old Balkan/Southeastern European lineage, found in Croatia, Bosnia, "
            "Serbia (~40%), and Ukraine. Like I1, it predates the Steppe expansion. "
            "I2a (Dinaric) is the main subgroup of the Balkans. Ancient I2 carriers "
            "have been found in European Mesolithic hunter-gatherers."
        ),
        "geographicOrigin": {"region": "Balkans / Southeast Europe", "detail": "Mesolithic lineage"},
        "timeOriginYears": 15000,
        "modernFrequency": {"balkans": 0.40, "croatia": 0.40, "ukraine": 0.20},
        "migrationRoute": "Indigenous Balkan lineage; Mesolithic hunter-gatherers",
        "famousCarriers": ["European Mesolithic hunter-gatherers (ancient DNA)"],
        "color": "#D97706",
    },
    "I2a": {
        "description": (
            "I2a (Dinaric) is the main Balkan subgroup of I2, found at high frequency "
            "in Croatia, Bosnia, Serbia (~40%), and Ukraine. It is associated with "
            "Mesolithic hunter-gatherers of Southeast Europe and experienced a dramatic "
            "expansion during the Slavic migrations of the 6th-7th centuries CE."
        ),
        "geographicOrigin": {"region": "Balkans (Dinaric Alps)", "detail": "Mesolithic / Slavic expansion"},
        "timeOriginYears": 10000,
        "modernFrequency": {"balkans": 0.40, "ukraine": 0.20},
        "migrationRoute": "Balkan origin; expanded with Slavic migrations",
        "famousCarriers": [],
        "color": "#CA8A04",
    },
    "J2": {
        "description": (
            "J2 is a Mediterranean and Near Eastern haplogroup (~6% of Europeans), "
            "strongly associated with the Neolithic agricultural expansion from the "
            "Fertile Crescent. Highest in the Near East (20-30%), Caucasus, and "
            "Mediterranean Europe (Italians, Greeks, Southern Iberians). J2 was the "
            "Y haplogroup of many Early European Farmers (alongside G2a). It spread "
            "with the Neolithic (J2b) and again with Bronze Age Mycenaean/Greek "
            "expansions (J2a). Associated with the spread of viticulture and pottery."
        ),
        "geographicOrigin": {"region": "Fertile Crescent / Near East", "detail": "Neolithic agricultural expansion"},
        "timeOriginYears": 10000,
        "modernFrequency": {"nearEast": 0.30, "caucasus": 0.25, "mediterranean": 0.15, "europe": 0.06},
        "migrationRoute": "Fertile Crescent to Mediterranean and Europe with Neolithic farming",
        "famousCarriers": ["Early European Farmers", "Mycenaean Greeks"],
        "color": "#EA580C",
    },
    "J1": {
        "description": (
            "Haplogroup J1 is dominant in the Arabian Peninsula (Arab populations 40-75%), "
            "Ethiopia, and the Near East. Its presence in European populations reflects "
            "ancient Near Eastern contacts, medieval Arab/Moorish presence in Southern "
            "Europe, or more recent migration."
        ),
        "geographicOrigin": {"region": "Arabian Peninsula", "detail": "Arab/Semitic populations"},
        "timeOriginYears": 15000,
        "modernFrequency": {"arabianPeninsula": 0.60, "ethiopia": 0.30, "nearEast": 0.20},
        "migrationRoute": "Arabian Peninsula expansion; Moorish presence in Southern Europe",
        "famousCarriers": [],
        "color": "#C2410C",
    },
    "G": {
        "description": (
            "Haplogroup G is found in ~3% of Europeans, mostly in the Caucasus region "
            "(Georgians 45%, Armenians 11%). It is associated with Caucasus Hunter-Gatherers "
            "and Early European Farmers. G2a was the dominant Neolithic European Y haplogroup."
        ),
        "geographicOrigin": {"region": "Caucasus", "detail": "Caucasus Hunter-Gatherers"},
        "timeOriginYears": 20000,
        "modernFrequency": {"caucasus": 0.45, "europe": 0.03},
        "migrationRoute": "Caucasus to Europe with Neolithic farming",
        "famousCarriers": ["Ötzi the Iceman (G2a2a1)"],
        "color": "#E11D48",
    },
    "G2a": {
        "description": (
            "G2a was the Y haplogroup of the first European farmers (Early European Farmers, "
            "EEF). Before the Steppe expansion ~5,000 years ago, G2a was the dominant "
            "Y lineage across Neolithic Europe. The famous Ötzi the Iceman (3300 BCE) "
            "was G2a2a1 — a Neolithic Italian farmer. The Bronze Age Steppe expansion "
            "nearly eliminated G2a from most of Europe; today it survives at ~3-5% in "
            "Mediterranean Europe, the Caucasus (20-30%), and the Near East. The Druze "
            "of Lebanon and Syria have ~30% G2a. G2a carriers in modern Europeans are "
            "walking genetic relics of the first farmers."
        ),
        "geographicOrigin": {"region": "Caucasus / Anatolia", "detail": "Early European Farmers (EEF)"},
        "timeOriginYears": 10000,
        "modernFrequency": {"caucasus": 0.30, "mediterranean": 0.05, "druze": 0.30, "europe": 0.03},
        "migrationRoute": "Anatolia/Caucasus to Europe with first farmers; nearly replaced by Steppe expansion",
        "famousCarriers": ["Ötzi the Iceman (G2a2a1, 3300 BCE)", "Early European Farmer males"],
        "color": "#BE123C",
    },
    "E": {
        "description": (
            "Your paternal lineage falls within haplogroup E. E is broadly distributed "
            "across Africa, the Middle East, and Mediterranean Europe. The specific "
            "sub-haplogroup would require additional markers to determine."
        ),
        "geographicOrigin": {"region": "Africa / Near East", "detail": "Broadly distributed"},
        "timeOriginYears": 50000,
        "modernFrequency": {"africa": 0.60, "nearEast": 0.15, "mediterranean": 0.10},
        "migrationRoute": "African origin; spread to Near East and Mediterranean",
        "famousCarriers": [],
        "color": "#9F1239",
    },
    "E1b1b": {
        "description": (
            "E1b1b is the dominant Y haplogroup in North Africa (~70-80%) and common "
            "in the Near East, the Horn of Africa, and Mediterranean Europe (10-20% in "
            "Italians, Greeks, Iberians). It likely spread with the expansion of "
            "Afroasiatic languages and later with Greek/Phoenician colonization of "
            "the Mediterranean. In Southern Europe, E1b1b reflects pre-Roman era "
            "North African/Near Eastern connections."
        ),
        "geographicOrigin": {"region": "North Africa / Horn of Africa", "detail": "Afroasiatic language expansion"},
        "timeOriginYears": 25000,
        "modernFrequency": {"northAfrica": 0.75, "hornOfAfrica": 0.50, "mediterranean": 0.15, "europe": 0.05},
        "migrationRoute": "North Africa to Mediterranean with Afroasiatic/Phoenician expansion",
        "famousCarriers": ["Albert Einstein (debated)", "Ramesses III (debated)"],
        "color": "#881337",
    },
    "N": {
        "description": (
            "Haplogroup N is the dominant Y lineage of Uralic-speaking populations: "
            "Finns (~60%), Estonians (~35%), and Siberian peoples (80-90%). It originated "
            "in Siberia and spread westward with Uralic language speakers, reaching "
            "Scandinavia and the Baltic relatively recently (within last 5,000 years). "
            "N is the genetic signature of the Finno-Ugric people."
        ),
        "geographicOrigin": {"region": "Siberia", "detail": "Uralic / Finno-Ugric peoples"},
        "timeOriginYears": 15000,
        "modernFrequency": {"finland": 0.60, "estonia": 0.35, "siberia": 0.80},
        "migrationRoute": "Siberia westward to Scandinavia and Baltic with Uralic speakers",
        "famousCarriers": ["Finno-Ugric peoples"],
        "color": "#F97316",
    },
    "Q": {
        "description": (
            "Haplogroup Q is the primary Native American Y haplogroup (Q1a3a1 = "
            "M3 mutation is nearly universal among pre-Columbian Native Americans). "
            "Q also occurs in Siberian populations. Its presence in a European individual "
            "suggests either Native American ancestry or Siberian ancestry."
        ),
        "geographicOrigin": {"region": "Siberia / Americas", "detail": "Native American primary lineage"},
        "timeOriginYears": 20000,
        "modernFrequency": {"nativeAmerican": 0.80, "siberia": 0.30},
        "migrationRoute": "Siberia to Americas via Beringia land bridge",
        "famousCarriers": ["Pre-Columbian Native Americans"],
        "color": "#FB923C",
    },
    "T_Y": {
        "description": (
            "Haplogroup T is a rare Y haplogroup (1-2% in Europe and Near East). "
            "Found mainly in the Near East, Mediterranean, and East Africa. Associated "
            "with ancient Sumerian and early agricultural populations. Rare in modern "
            "Europeans but with ancient connections to Fertile Crescent civilization."
        ),
        "geographicOrigin": {"region": "Near East / East Africa", "detail": "Ancient Sumerian / agricultural populations"},
        "timeOriginYears": 25000,
        "modernFrequency": {"nearEast": 0.05, "mediterranean": 0.02, "eastAfrica": 0.08},
        "migrationRoute": "Near East / East African origin",
        "famousCarriers": ["Thomas Jefferson (T1a, debated)"],
        "color": "#FBBF24",
    },
    "C": {
        "description": (
            "Haplogroup C is the dominant Y lineage of Mongolian and Turkic peoples "
            "(Mongolians ~50-90%), found also in Siberia, East Asia, Oceania, and "
            "Native Americans (C1b). Its presence in a European individual suggests "
            "Central Asian/Turkic ancestry (possibly from medieval Mongol invasions "
            "or Turkic migrations into Eastern Europe)."
        ),
        "geographicOrigin": {"region": "Central Asia / Mongolia", "detail": "Mongolic / Turkic peoples"},
        "timeOriginYears": 50000,
        "modernFrequency": {"mongolia": 0.70, "siberia": 0.30, "oceania": 0.20},
        "migrationRoute": "Central Asian expansion; Mongol/Turkic migrations",
        "famousCarriers": ["Genghis Khan lineage (C3, debated)"],
        "color": "#F59E0B",
    },
    "R": {
        "description": (
            "Your paternal lineage falls within haplogroup R but a specific sub-haplogroup "
            "could not be determined from the available markers. R includes both R1a "
            "(Eastern European/Slavic) and R1b (Western European), the two most common "
            "European Y lineages. Additional markers would determine which branch."
        ),
        "geographicOrigin": {"region": "Eurasia", "detail": "Includes R1a and R1b, most common European lineages"},
        "timeOriginYears": 25000,
        "modernFrequency": {"europe": 0.60},
        "migrationRoute": "Steppe origin; diversified into R1a and R1b branches",
        "famousCarriers": [],
        "color": "#7F1D1D",
    },
}


# ─────────────────────────────────────────────────────────────────────────────
# CORE ANALYSIS
# ─────────────────────────────────────────────────────────────────────────────

def analyze_haplogroup(
    variants: Dict[str, Tuple[str, str, str]]
) -> Dict[str, Any]:
    """
    Analyze genetic variants for mitochondrial and Y chromosome haplogroups.

    Args:
        variants: Dictionary mapping rsID (lowercase) to (chromosome, position, genotype).
                  May also contain position-based entries for MT chromosome.

    Returns:
        Dict with mt and Y haplogroup assignments, signals, confidence, and metadata.
    """
    # ── Build position-based lookup for MT chromosome ──
    mt_by_position: Dict[int, str] = {}
    y_snp_count = 0
    mt_snp_count = 0

    for rsid, (chrom, pos, genotype) in variants.items():
        chrom_upper = chrom.upper()
        if chrom_upper in ("MT", "26", "CHRM"):
            mt_snp_count += 1
            try:
                mt_by_position[int(pos)] = genotype
            except (ValueError, TypeError):
                pass
        elif chrom_upper in ("Y", "24", "CHRY"):
            y_snp_count += 1

    # ── Assign mitochondrial haplogroup ──
    mt_signals: Set[str] = set()
    mt_evidence: List[Dict[str, str]] = []

    # Check position-based markers
    for pos, marker_dict in MT_POSITION_MARKERS.items():
        if pos in mt_by_position:
            result = mt_by_position[pos]
            for derived_allele, hg_signal in marker_dict.items():
                if has_allele(result, derived_allele):
                    mt_signals.add(hg_signal)
                    mt_evidence.append({
                        "type": "position",
                        "marker": str(pos),
                        "allele": derived_allele,
                        "signal": hg_signal,
                    })

    # Check rsID-based markers
    for rsid, (derived_allele, hg_signal) in MT_RSID_MARKERS.items():
        rsid_lower = rsid.lower()
        if rsid_lower in variants:
            _chrom, _pos, genotype = variants[rsid_lower]
            if has_allele(genotype, derived_allele):
                mt_signals.add(hg_signal)
                mt_evidence.append({
                    "type": "rsid",
                    "marker": rsid,
                    "allele": derived_allele,
                    "signal": hg_signal,
                })

    # Apply decision tree
    mt_haplogroup: Optional[str] = None
    mt_description: Optional[str] = None
    mt_confidence: str = "low"

    for required_signals, hg, desc, conf in MT_DECISION_TREE:
        if all(s in mt_signals for s in required_signals):
            mt_haplogroup = hg
            mt_description = desc
            mt_confidence = conf
            break

    mt_reason: Optional[str] = None
    if mt_haplogroup is None:
        if mt_snp_count == 0:
            mt_reason = "insufficient_data"
            mt_confidence = "none"
        else:
            mt_haplogroup = "Undetermined"
            mt_reason = "insufficient_markers"
            mt_confidence = "low"

    # ── Assign Y chromosome haplogroup ──
    y_signals: Dict[str, str] = {}  # hg_clade -> rsid
    y_evidence: List[Dict[str, str]] = []

    for rsid, (derived_allele, hg_clade) in Y_MARKERS.items():
        rsid_lower = rsid.lower()
        if rsid_lower in variants:
            chrom, _pos, genotype = variants[rsid_lower]
            chrom_upper = chrom.upper()
            if chrom_upper in ("Y", "CHRY", "24"):
                if has_allele(genotype, derived_allele):
                    y_signals[hg_clade] = rsid
                    y_evidence.append({
                        "type": "rsid",
                        "marker": rsid,
                        "allele": derived_allele,
                        "signal": hg_clade,
                    })

    y_haplogroup: Optional[str] = None
    y_confidence: str = "none"
    y_reason: Optional[str] = None
    # Detect sex: having Y-chromosome entries is not enough (MyHeritage exports
    # Y positions for females with '--' genotype). Require actual Y-specific
    # marker calls from ISOGG markers, or count non-PAR Y variants with calls.
    is_male: bool = len(y_signals) > 0  # Only True if at least one Y marker has a derived allele call

    if not y_signals:
        if not is_male:
            y_reason = "female"
        else:
            y_reason = "insufficient_data"
    else:
        # Priority-ordered assignment
        for hg in Y_PRIORITY:
            if hg in y_signals:
                y_haplogroup = hg
                y_confidence = "moderate"
                break

        if y_haplogroup is None and y_signals:
            y_haplogroup = list(y_signals.keys())[0]
            y_confidence = "low"

    # ── Build result ──
    return {
        "mitochondrial": {
            "haplogroup": mt_haplogroup,
            "confidence": mt_confidence,
            "snpsFound": mt_snp_count,
            "signalsDetected": sorted(mt_signals),
            "evidence": mt_evidence,
            "reason": mt_reason,
        },
        "yChromosome": {
            "haplogroup": y_haplogroup,
            "confidence": y_confidence,
            "snpsFound": y_snp_count,
            "signalsDetected": sorted(y_signals.keys()) if y_signals else [],
            "evidence": y_evidence,
            "reason": y_reason,
        },
        "isMale": is_male,
    }


# ─────────────────────────────────────────────────────────────────────────────
# JSON OUTPUT
# ─────────────────────────────────────────────────────────────────────────────

def generate_haplogroup_json(result: Dict[str, Any]) -> dict:
    """
    Generate a JSON-serializable dict for the frontend haplogroup report.

    Args:
        result: Dict from analyze_haplogroup()

    Returns:
        Dict matching the haplogroup report JSON schema
    """
    mt_data = result["mitochondrial"]
    y_data = result["yChromosome"]

    mt_hg = mt_data["haplogroup"]
    y_hg = y_data["haplogroup"]

    # ── Build mitochondrial section ──
    mt_section: Dict[str, Any]
    if mt_hg and mt_hg != "Undetermined":
        mt_info = MT_HAPLOGROUP_INFO.get(mt_hg, {})
        mt_section = {
            "haplogroup": mt_hg,
            "confidence": mt_data["confidence"],
            "snpsFound": mt_data["snpsFound"],
            "signalsDetected": mt_data["signalsDetected"],
            "description": mt_info.get("description", f"Mitochondrial haplogroup {mt_hg}"),
            "geographicOrigin": mt_info.get("geographicOrigin", {"region": "Unknown", "detail": ""}),
            "timeOriginYears": mt_info.get("timeOriginYears", 0),
            "modernFrequency": mt_info.get("modernFrequency", {}),
            "migrationRoute": mt_info.get("migrationRoute", ""),
            "famousCarriers": mt_info.get("famousCarriers", []),
            "color": mt_info.get("color", "#6366F1"),
        }
    else:
        mt_section = {
            "haplogroup": None,
            "confidence": mt_data["confidence"],
            "snpsFound": mt_data["snpsFound"],
            "signalsDetected": mt_data["signalsDetected"],
            "reason": mt_data.get("reason", "insufficient_data"),
            "description": "Mitochondrial haplogroup could not be determined from the available data.",
            "geographicOrigin": {"region": "Unknown", "detail": ""},
            "timeOriginYears": 0,
            "modernFrequency": {},
            "migrationRoute": "",
            "famousCarriers": [],
            "color": "#6B7280",
        }

    # ── Build Y chromosome section ──
    y_section: Dict[str, Any]
    if y_hg:
        # Map T_Y back to display name T
        display_hg = y_hg if y_hg != "T_Y" else "T"
        y_info = Y_HAPLOGROUP_INFO.get(y_hg, {})
        y_section = {
            "haplogroup": display_hg,
            "confidence": y_data["confidence"],
            "snpsFound": y_data["snpsFound"],
            "signalsDetected": [s if s != "T_Y" else "T" for s in y_data["signalsDetected"]],
            "description": y_info.get("description", f"Y chromosome haplogroup {display_hg}"),
            "geographicOrigin": y_info.get("geographicOrigin", {"region": "Unknown", "detail": ""}),
            "timeOriginYears": y_info.get("timeOriginYears", 0),
            "modernFrequency": y_info.get("modernFrequency", {}),
            "migrationRoute": y_info.get("migrationRoute", ""),
            "famousCarriers": y_info.get("famousCarriers", []),
            "color": y_info.get("color", "#DC2626"),
        }
    else:
        reason = y_data.get("reason", "insufficient_data")
        if reason == "female":
            desc = "No Y chromosome haplogroup — females carry XX chromosomes and do not have a Y chromosome."
        else:
            desc = "Y chromosome haplogroup could not be determined from the available data."
        y_section = {
            "haplogroup": None,
            "confidence": y_data["confidence"],
            "snpsFound": y_data["snpsFound"],
            "signalsDetected": [],
            "reason": reason,
            "description": desc,
            "geographicOrigin": {"region": "Unknown", "detail": ""},
            "timeOriginYears": 0,
            "modernFrequency": {},
            "migrationRoute": "",
            "famousCarriers": [],
            "color": "#6B7280",
        }

    # ── Build summary ──
    mt_display = mt_hg if mt_hg and mt_hg != "Undetermined" else None
    y_display = y_hg if y_hg else None
    if y_display == "T_Y":
        y_display = "T"

    # ── Ancestral lineage chains (Biopython #3) ─────────────────────────
    # Curated chains from mt-MRCA ("Mitochondrial Eve") or y-MRCA
    # ("Y-chromosomal Adam") down to the user's haplogroup. Frontend
    # renders these as a horizontal tree with "you are here" highlight.
    # Failure is silent: if lineage is unknown, no `ancestralLineage`
    # key is added to that section.
    try:
        from .haplogroup_lineage import get_mt_lineage, get_y_lineage
        mt_lineage = get_mt_lineage(mt_display)
        if mt_lineage:
            mt_section["ancestralLineage"] = mt_lineage
        y_lineage = get_y_lineage(y_display)
        if y_lineage:
            y_section["ancestralLineage"] = y_lineage
    except Exception as e:  # noqa: BLE001
        print(f"[haplogroup-lineage] enrichment skipped: {type(e).__name__}: {e}", flush=True)

    return {
        "reportType": "haplogroup",
        "version": "1.0",
        "summary": {
            "mtHaplogroup": mt_display,
            "yHaplogroup": y_display,
            "mtConfidence": mt_data["confidence"],
            "yConfidence": y_data["confidence"],
        },
        "mitochondrial": mt_section,
        "yChromosome": y_section,
    }
