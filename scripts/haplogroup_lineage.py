"""
Haplogroup ancestral lineage — curated chains from the root of the tree
(Mitochondrial Eve for mtDNA, Y-chromosomal Adam for Y-DNA) down to the
user's haplogroup.

Produces a structured payload the frontend renders as a horizontal tree
with "you are here" highlight. Approximate dates are in years Before
Present (BP), sourced from peer-reviewed population genetics (Behar 2012,
Karmin 2015, Fu 2018, Jobling & Tyler-Smith 2017, ISOGG 2024, Phylotree
Build 17). Ages are approximations — exact dating still debated.

The key insight for the user: their haplogroup is not a dead-end
classification — it's a waypoint in a ~200,000-year chain going back to a
single woman (mtDNA) or man (Y-DNA).

Public API:
    from scripts.haplogroup_lineage import get_mt_lineage, get_y_lineage

    payload = get_mt_lineage("H1")
    # -> {
    #   "kind": "mitochondrial",
    #   "rootLabel": "Mitochondrial Eve",
    #   "totalSpanYears": 200000,
    #   "nodes": [
    #       {"code": "mt-MRCA", "label": "Mitochondrial Eve", "ageBP": 200000, "region": "Africa", "summary": "..."},
    #       {"code": "L",      "label": "L", "ageBP": 170000, "region": "Africa",   "summary": "..."},
    #       {"code": "L3",     "label": "L3", "ageBP": 70000,  "region": "East Africa", "summary": "..."},
    #       {"code": "N",      "label": "N",  "ageBP": 65000,  "region": "Near East",   "summary": "..."},
    #       {"code": "R",      "label": "R",  "ageBP": 60000,  "region": "West Asia",   "summary": "..."},
    #       {"code": "HV",     "label": "HV", "ageBP": 28000,  "region": "Caucasus",    "summary": "..."},
    #       {"code": "H",      "label": "H",  "ageBP": 25000,  "region": "West Asia / Europe", "summary": "..."},
    #       {"code": "H1",     "label": "H1", "ageBP": 13000,  "region": "Iberia refugium",    "summary": "...", "isUser": true},
    #   ]
    # }

Failure mode: unknown haplogroup returns None.
"""

from __future__ import annotations

from typing import Any, Dict, List, Optional

# ─── Node metadata (shared across lineages) ──────────────────────────────

_NODE_META: Dict[str, Dict[str, Any]] = {
    # ── mtDNA root + African ──
    "mt-MRCA": {
        "label": "Mitochondrial Eve",
        "ageBP": 200000,
        "region": "East Africa",
        "summary": "The most recent common matrilineal ancestor of all living humans. Not the only woman alive then — just the one whose mtDNA lineage survived to today through an unbroken chain of daughters.",
    },
    "L": {"label": "L", "ageBP": 170000, "region": "Africa", "summary": "The deepest split. All modern humans descend from L. Subclades L0–L6 still dominate Africa."},
    "L0": {"label": "L0", "ageBP": 160000, "region": "Southern Africa", "summary": "The oldest surviving mtDNA branch. Common in Khoisan peoples — living descendants of the most divergent living human lineage."},
    "L1": {"label": "L1", "ageBP": 140000, "region": "Central Africa", "summary": "Second-oldest branch. Common among Central African foragers including the Mbuti and Biaka."},
    "L2": {"label": "L2", "ageBP": 90000, "region": "West/Central Africa", "summary": "Major lineage across West Africa. Brought to the Americas via the Atlantic slave trade — the most common African mtDNA among African-diaspora communities."},
    "L3": {"label": "L3", "ageBP": 70000, "region": "East Africa", "summary": "The mother lineage of everyone outside Africa. A single L3 woman's descendants left Africa ~60,000 years ago and populated the world."},

    # ── Eurasian trunk ──
    "M": {"label": "M", "ageBP": 60000, "region": "Arabia / South Asia", "summary": "One of the two great out-of-Africa branches. M spread east along the southern Asian coast and dominates South and East Asia today."},
    "N": {"label": "N", "ageBP": 65000, "region": "Near East", "summary": "The other out-of-Africa branch. Ancestor of almost all West Eurasian mtDNA, plus East Asian A, N9a, Y."},
    "R": {"label": "R", "ageBP": 60000, "region": "West Asia", "summary": "Major sub-branch of N. Spread rapidly ~50,000 BP — ancestor of most European and Near-Eastern mtDNA lineages."},
    "U": {"label": "U", "ageBP": 55000, "region": "West Asia", "summary": "One of the oldest European lineages. Peak frequency in European hunter-gatherers before the Neolithic."},
    "U4": {"label": "U4", "ageBP": 25000, "region": "Northeast Europe / Siberia", "summary": "Classic hunter-gatherer lineage of post-Ice-Age Europe. Peaks in Volga-Uralic populations and among hunter-fisher remains."},
    "U5": {"label": "U5", "ageBP": 30000, "region": "Pre-glacial Europe", "summary": "The signature mtDNA of European Paleolithic hunter-gatherers. Found at ~11% in modern Europeans and up to 50% in Saami populations."},

    "HV": {"label": "HV", "ageBP": 28000, "region": "Caucasus / Near East", "summary": "Cradle of the H and V branches. Expanded out of Southwest Asia during and after the Last Glacial Maximum."},
    "H": {"label": "H", "ageBP": 25000, "region": "West Asia / Europe", "summary": "The single most common European mtDNA — carried by ~40% of Europeans. Expanded explosively after the last Ice Age from Franco-Cantabrian refugia."},
    "H1": {"label": "H1", "ageBP": 13000, "region": "Iberian refugium", "summary": "Post-glacial expansion from the Iberian peninsula as ice sheets retreated. Peak frequency today in the Basques, Berbers, and Iberians."},
    "H3": {"label": "H3", "ageBP": 11000, "region": "Iberia / Atlantic Europe", "summary": "Sister clade of H1. Strong in Western Europe, especially Portugal, Sardinia, and Britain. Associated with Mesolithic re-population of Europe."},
    "V": {"label": "V", "ageBP": 15000, "region": "Iberia / Caucasus", "summary": "A post-glacial European lineage centered on Iberia. Reaches up to 30%+ among the Saami of northern Scandinavia."},

    "JT": {"label": "JT", "ageBP": 50000, "region": "Near East", "summary": "The common mother of J and T — two of the major Neolithic lineages that spread farming across Europe."},
    "J": {"label": "J", "ageBP": 45000, "region": "Near East", "summary": "Classic Neolithic farmer lineage. Spread with agriculture from the Fertile Crescent into Europe ~8,000 BP."},
    "T": {"label": "T", "ageBP": 30000, "region": "Near East / Mesopotamia", "summary": "Another Neolithic expansion marker. High in modern Near Eastern populations, spread with farming into Europe and northern Africa."},

    "I": {"label": "I", "ageBP": 30000, "region": "West Asia / Europe", "summary": "A rare European lineage carried by ~2-3%. Present in several ancient Steppe samples."},
    "W": {"label": "W", "ageBP": 23000, "region": "South Caucasus", "summary": "Low-frequency European and South-Asian lineage. Found in Indian Tharu populations and Eastern Europeans."},
    "X": {"label": "X", "ageBP": 30000, "region": "Near East", "summary": "Unusual distribution: Europe, North America, and pockets of Near East. One of the five founding mtDNA lineages of Native Americans."},

    "K": {"label": "K", "ageBP": 16000, "region": "Near East", "summary": "Sub-branch of U8. Common among Ashkenazi Jews and Druze. The famous Ötzi the Iceman carried K1a."},

    # ── Y-DNA root ──
    "y-MRCA": {
        "label": "Y-chromosomal Adam",
        "ageBP": 275000,
        "region": "Africa",
        "summary": "The most recent common patrilineal ancestor of all living men. Not the only man alive then — the one whose Y-chromosome lineage made it through ~10,000 generations of unbroken father-to-son transmission.",
    },
    "A": {"label": "A", "ageBP": 270000, "region": "Africa", "summary": "The deepest Y-DNA lineage. Carried mainly by Khoisan foragers of southern Africa."},
    "BT": {"label": "BT", "ageBP": 130000, "region": "Africa", "summary": "The paternal trunk of every man outside of deep African A lineages."},
    "CT": {"label": "CT", "ageBP": 88000, "region": "Africa", "summary": "Carried by nearly every non-African man and many African men. The moment the Y-tree expanded outward."},
    "CF": {"label": "CF", "ageBP": 68000, "region": "Near East / Arabia", "summary": "The bottleneck lineage that produced all Eurasian, American, and Oceanian males outside the DE haplogroup."},
    "DE": {"label": "DE", "ageBP": 76000, "region": "Africa / Asia", "summary": "A deep paternal split — D became dominant in Tibet and Japan, E in Africa."},
    "F": {"label": "F", "ageBP": 55000, "region": "South Asia / Near East", "summary": "F is the ancestor of all non-African Y-DNA except D. Over 90% of men outside Africa trace to F."},
    "IJK": {"label": "IJK", "ageBP": 48000, "region": "West Asia", "summary": "Trunk of the great Eurasian patrilines I, J, K — which split into the lineages that later covered Europe, the Near East, and the Steppe."},
    "K": {"label": "K", "ageBP": 47000, "region": "South Asia", "summary": "The Eurasian super-haplogroup. Ancestor of N, O, P — covering Siberia, East Asia, Americas."},
    "P": {"label": "P", "ageBP": 45000, "region": "Central Asia / Siberia", "summary": "The Steppe-and-beyond ancestor of Q and R — carried by early Paleolithic colonizers of Siberia and the New World."},
    "R": {"label": "R", "ageBP": 27000, "region": "Central Asia", "summary": "One of the dominant Y lineages in modern Europe and South Asia. Split into R1 (~22,000 BP) and R2."},
    "R1": {"label": "R1", "ageBP": 22000, "region": "Western Eurasia", "summary": "Ancestor of R1a (Slavic / Indo-Iranian) and R1b (Western European) — the two most frequent Y-haplogroups in modern Europe."},
    "R1b": {"label": "R1b", "ageBP": 18000, "region": "West Asia → Europe", "summary": "Carried by ~55% of Western European men. Expanded dramatically with the Bronze Age Yamnaya and Bell Beaker cultures ~4,500 BP."},
    "R1b1a": {"label": "R1b1a", "ageBP": 13000, "region": "Pontic Steppe", "summary": "Pre-Yamnaya ancestor of most Western European R1b. Expanded out of the Pontic-Caspian steppe carrying early Indo-European languages."},
    "R1a": {"label": "R1a", "ageBP": 22000, "region": "Pontic / Caspian", "summary": "Dominant in Eastern Europe, Scandinavia, and among Indo-Aryan-speaking populations in India. Reached 50%+ frequency in modern Slavs."},

    "I": {"label": "I", "ageBP": 27000, "region": "Europe", "summary": "One of the earliest European Y-DNA lineages. Carried by paleolithic European hunter-gatherers before farmers and Steppe pastoralists arrived."},
    "I1": {"label": "I1", "ageBP": 4600, "region": "Scandinavia", "summary": "Young but dominant lineage of Scandinavia and Northern Europe. Expanded from a single man ~4,600 years ago — possibly through the Nordic Bronze Age."},
    "I2": {"label": "I2", "ageBP": 21000, "region": "Europe", "summary": "Ancient European paternal lineage. High frequencies today in the Balkans and Sardinia — thought to represent surviving pockets of Mesolithic hunter-gatherer ancestry."},

    "J": {"label": "J", "ageBP": 32000, "region": "Mesopotamia / Caucasus", "summary": "The Near Eastern and Semitic lineage. Split into J1 (Arabia) and J2 (Anatolia / Levant / Mediterranean)."},
    "J1": {"label": "J1", "ageBP": 18000, "region": "Arabian Peninsula", "summary": "Dominant lineage of Arabia and the Semitic world. High frequencies in Bedouin, Cohanim, and Saudi populations."},
    "J2": {"label": "J2", "ageBP": 27000, "region": "Anatolia / Levant", "summary": "Classic Neolithic farmer lineage. Spread with agriculture around the Mediterranean and into South Asia."},

    "G": {"label": "G", "ageBP": 48000, "region": "Near East / Caucasus", "summary": "Carried by Ötzi the Iceman. Marker of the Neolithic farmer expansion into Europe ~9,000 BP."},

    "E": {"label": "E", "ageBP": 70000, "region": "East Africa / Near East", "summary": "Rose from DE in Africa. Dominant in modern Africa and widespread in the Mediterranean through E1b1b."},

    "N": {"label": "N", "ageBP": 20000, "region": "Siberia", "summary": "Siberian lineage carried by Uralic-speaking peoples (Finns, Saami, Yakuts). Associated with the spread of reindeer herding across the taiga."},

    "Q": {"label": "Q", "ageBP": 28000, "region": "Central Asia / Siberia", "summary": "Classic founder lineage of the Americas. Q1a-M3 is carried by ~90% of Native American men — tracing a single ancestor who crossed Beringia."},

    "T_Y": {"label": "T", "ageBP": 28000, "region": "Near East", "summary": "Uncommon lineage centered on the Horn of Africa and parts of Anatolia and the Iberian peninsula."},

    "C": {"label": "C", "ageBP": 50000, "region": "Central Asia / Oceania", "summary": "One of the first lineages outside Africa. Reached Oceania 45,000 BP, Siberia via C2, and the Americas in small numbers. C2 is associated with Genghis Khan's male-line descendants."},
}


# ─── Chains (leaf → root, as a shorthand list of node codes) ──────────────
#
# Each chain is the sequence from oldest ancestor to the user's haplogroup.
# Keeping them flat (as lists) lets us share intermediate nodes across
# chains without maintaining a proper tree structure.

_MT_CHAINS: Dict[str, List[str]] = {
    "L": ["mt-MRCA", "L"],
    "L0": ["mt-MRCA", "L", "L0"],
    "L1": ["mt-MRCA", "L", "L1"],
    "L2": ["mt-MRCA", "L", "L2"],
    "L3": ["mt-MRCA", "L", "L3"],
    "M": ["mt-MRCA", "L", "L3", "M"],
    "N": ["mt-MRCA", "L", "L3", "N"],
    "R": ["mt-MRCA", "L", "L3", "N", "R"],
    "U": ["mt-MRCA", "L", "L3", "N", "R", "U"],
    "U4": ["mt-MRCA", "L", "L3", "N", "R", "U", "U4"],
    "U5": ["mt-MRCA", "L", "L3", "N", "R", "U", "U5"],
    "HV": ["mt-MRCA", "L", "L3", "N", "R", "HV"],
    "H": ["mt-MRCA", "L", "L3", "N", "R", "HV", "H"],
    "H1": ["mt-MRCA", "L", "L3", "N", "R", "HV", "H", "H1"],
    "H3": ["mt-MRCA", "L", "L3", "N", "R", "HV", "H", "H3"],
    "V": ["mt-MRCA", "L", "L3", "N", "R", "HV", "V"],
    "JT": ["mt-MRCA", "L", "L3", "N", "R", "JT"],
    "J": ["mt-MRCA", "L", "L3", "N", "R", "JT", "J"],
    "T": ["mt-MRCA", "L", "L3", "N", "R", "JT", "T"],
    "I": ["mt-MRCA", "L", "L3", "N", "I"],
    "W": ["mt-MRCA", "L", "L3", "N", "W"],
    "X": ["mt-MRCA", "L", "L3", "N", "X"],
    "K": ["mt-MRCA", "L", "L3", "N", "R", "U", "K"],
}

_Y_CHAINS: Dict[str, List[str]] = {
    "A": ["y-MRCA", "A"],
    "BT": ["y-MRCA", "A", "BT"],
    "CT": ["y-MRCA", "A", "BT", "CT"],
    "DE": ["y-MRCA", "A", "BT", "CT", "DE"],
    "E": ["y-MRCA", "A", "BT", "CT", "DE", "E"],
    "CF": ["y-MRCA", "A", "BT", "CT", "CF"],
    "F": ["y-MRCA", "A", "BT", "CT", "CF", "F"],
    "G": ["y-MRCA", "A", "BT", "CT", "CF", "F", "G"],
    "IJK": ["y-MRCA", "A", "BT", "CT", "CF", "F", "IJK"],
    "I": ["y-MRCA", "A", "BT", "CT", "CF", "F", "IJK", "I"],
    "I1": ["y-MRCA", "A", "BT", "CT", "CF", "F", "IJK", "I", "I1"],
    "I2": ["y-MRCA", "A", "BT", "CT", "CF", "F", "IJK", "I", "I2"],
    "J": ["y-MRCA", "A", "BT", "CT", "CF", "F", "IJK", "J"],
    "J1": ["y-MRCA", "A", "BT", "CT", "CF", "F", "IJK", "J", "J1"],
    "J2": ["y-MRCA", "A", "BT", "CT", "CF", "F", "IJK", "J", "J2"],
    "K": ["y-MRCA", "A", "BT", "CT", "CF", "F", "IJK", "K"],
    "N": ["y-MRCA", "A", "BT", "CT", "CF", "F", "IJK", "K", "N"],
    "P": ["y-MRCA", "A", "BT", "CT", "CF", "F", "IJK", "K", "P"],
    "Q": ["y-MRCA", "A", "BT", "CT", "CF", "F", "IJK", "K", "P", "Q"],
    "R": ["y-MRCA", "A", "BT", "CT", "CF", "F", "IJK", "K", "P", "R"],
    "R1": ["y-MRCA", "A", "BT", "CT", "CF", "F", "IJK", "K", "P", "R", "R1"],
    "R1a": ["y-MRCA", "A", "BT", "CT", "CF", "F", "IJK", "K", "P", "R", "R1", "R1a"],
    "R1b": ["y-MRCA", "A", "BT", "CT", "CF", "F", "IJK", "K", "P", "R", "R1", "R1b"],
    "R1b1a": ["y-MRCA", "A", "BT", "CT", "CF", "F", "IJK", "K", "P", "R", "R1", "R1b", "R1b1a"],
    "C": ["y-MRCA", "A", "BT", "CT", "CF", "C"],
    "T_Y": ["y-MRCA", "A", "BT", "CT", "CF", "F", "IJK", "K", "T_Y"],
}


def _chain_to_payload(
    chain: List[str],
    user_code: str,
    kind: str,
    root_code: str,
) -> Dict[str, Any]:
    nodes: List[Dict[str, Any]] = []
    for code in chain:
        meta = _NODE_META.get(code, {})
        node = {
            "code": code,
            "label": meta.get("label", code),
            "ageBP": meta.get("ageBP"),
            "region": meta.get("region", ""),
            "summary": meta.get("summary", ""),
        }
        if code == user_code:
            node["isUser"] = True
        nodes.append(node)

    root_age = _NODE_META.get(root_code, {}).get("ageBP") or 0
    return {
        "kind": kind,  # "mitochondrial" | "yChromosome"
        "rootCode": root_code,
        "rootLabel": _NODE_META.get(root_code, {}).get("label", root_code),
        "userCode": user_code,
        "totalSpanYears": root_age,
        "nodeCount": len(nodes),
        "nodes": nodes,
    }


def get_mt_lineage(haplogroup: Optional[str]) -> Optional[Dict[str, Any]]:
    """Return the mitochondrial ancestral lineage for `haplogroup`, or None."""
    if not haplogroup:
        return None
    # Strip variant suffixes — try full match first, then progressively shorter.
    for probe in (haplogroup, haplogroup.split("+")[0], haplogroup.rstrip("abcdefghijklmnop0123456789")):
        chain = _MT_CHAINS.get(probe)
        if chain:
            return _chain_to_payload(chain, user_code=probe, kind="mitochondrial", root_code="mt-MRCA")
    return None


def get_y_lineage(haplogroup: Optional[str]) -> Optional[Dict[str, Any]]:
    """Return the Y-chromosome ancestral lineage for `haplogroup`, or None."""
    if not haplogroup:
        return None
    # Map display name "T" back to our internal "T_Y".
    canonical = "T_Y" if haplogroup == "T" else haplogroup
    for probe in (canonical, canonical.split("+")[0], canonical.rstrip("abcdefghijklmnop0123456789")):
        chain = _Y_CHAINS.get(probe)
        if chain:
            return _chain_to_payload(chain, user_code=probe, kind="yChromosome", root_code="y-MRCA")
    return None
