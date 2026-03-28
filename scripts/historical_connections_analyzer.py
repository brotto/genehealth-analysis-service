"""
Historical Connections Analyzer
Traces genetic connections between user haplogroups and historically confirmed
ancient DNA of notable historical figures.

METHODOLOGY:
  1. Y-DNA haplogroup matching (patrilineal line, males only)
  2. mtDNA haplogroup matching (matrilineal line, all individuals)
  3. Otzi phenotypic comparison (published genome, Keller et al. 2012)

THREE MATCH TIERS:
  5 stars  Confirmed   - direct ancient DNA sequencing from verified remains
  4 stars  High        - multiple confirmations or from authenticated relatives
  3 stars  Inferred    - population genetics study, not from individual remains
  2 stars  Contested   - reported but methodologically questioned
  1 star   Speculative - hypothesis only, no published data

IMPORTANT SCIENTIFIC CAVEAT:
  Shared haplogroup does NOT mean direct descent from that individual.
  It means sharing a common patrilineal or matrilineal ancestor somewhere
  in deep history. The shared ancestor could be 500 or 50,000 years back.

References:
  - King et al. 2014, Nature Communications - Richard III
  - Gill et al. 1994, Nature Genetics - Romanovs
  - Keller et al. 2012, Nature Communications - Otzi (full genome)
  - Begg et al. 2023, Current Biology - Beethoven
  - Hagelberg et al. 2011 - Napoleon (inferred from relatives)
  - Hawass et al. 2012, JAMA - Egyptian pharaohs (contested)
  - Brace et al. 2019, Nature Ecology & Evolution - Cheddar Man
  - Lazaridis et al. 2014, Nature - Loschbour Man
  - Zerjal et al. 2003, AJHG - Mongol imperial lineage (inferred)
"""

from typing import Dict, List, Tuple, Any, Optional
from collections import defaultdict


# ---------------------------------------------------------------------------
# HISTORICAL FIGURES DATABASE
# ---------------------------------------------------------------------------

HISTORICAL_FIGURES: List[Dict[str, Any]] = [
    {
        "name": "Richard III",
        "dates": "1452-1485 AD",
        "title": "Last Plantagenet King of England",
        "category": "royalty",
        "bio": (
            "Last Plantagenet king of England, killed at the Battle of Bosworth "
            "Field in 1485. His skeleton was discovered under a Leicester parking "
            "lot in 2012."
        ),
        "origin": "English (Plantagenet dynasty)",
        "y_hg": "G2a",
        "mt_hg": "J1c2c3",
        "y_matches": [
            ("G2a2b", 4, "Close G2a branch - lineage of early European farmers"),
            ("G2a", 5, "Exact match - G2a (G-P287), confirmed from skeleton"),
            ("G2", 4, "G2 clade"),
            ("G", 3, "Macrohaplogroup G - Caucasus/European Neolithic"),
        ],
        "mt_matches": [
            ("J1c2c3", 5, "Exact maternal haplogroup match"),
            ("J1c2c", 4, "Close branch - recent common ancestor"),
            ("J1c2", 3, "J1c2 clade"),
            ("J1c", 3, "J1c clade"),
            ("J1", 2, "J1 branch"),
            ("J", 1, "Macrohaplogroup J - Near Eastern/Neolithic origin"),
        ],
        "evidence_y": 5,
        "evidence_mt": 5,
        "source": "King et al. 2014, Nature Communications",
        "verified": True,
        "story": (
            "Last English king to die in battle (Battle of Bosworth Field, 1485). "
            "His skeleton was identified by DNA in 2014 - one of the greatest cases "
            "of historical forensic genetics. NOTE: The skeleton's Y-DNA is G2a "
            "(G-P287), NOT R1b. Living descendants of the Somerset family (expected "
            "paternal lineage) carry R1b-L21 - the divergence confirms a non-paternity "
            "event in the Plantagenet lineage between Edward III and Richard III. "
            "The mtDNA J1c2c3 was confirmed by matching with Joy Ibsen, a modern "
            "descendant through the maternal line of Anne of York."
        ),
    },
    {
        "name": "Nicholas II - Romanov Family",
        "dates": "1868-1918 AD",
        "title": "Last Tsar of Russia",
        "category": "royalty",
        "bio": (
            "Last Tsar of Russia, executed with his family in 1918. DNA identification "
            "resolved decades of mystery about the Romanov fate."
        ),
        "origin": "Russian-German (Romanov dynasty)",
        "y_hg": "R1b",
        "mt_hg": "T2",
        "y_matches": [
            ("R1b-U152", 3, "Distinct R1b branches - common ancestor ~5,000 years"),
            ("R1b-P312", 3, "Distinct R1b branches"),
            ("R1b-L11", 3, "Distinct R1b branches"),
            ("R1b", 2, "Shared macrohaplogroup R1b - Steppe Bronze Age lineage"),
            ("R1", 1, "R1 superlineage"),
        ],
        "mt_matches": [
            ("T2", 5, "Exact match - Romanov maternal haplogroup T2"),
            ("T", 4, "T clade - Anatolian origin"),
            ("JT", 2, "JT ancestor"),
        ],
        "evidence_y": 5,
        "evidence_mt": 5,
        "source": "Gill et al. 1994, Nature Genetics; Rogaev et al. 2009",
        "verified": True,
        "story": (
            "Nicholas II and his family were executed in 1918 in Yekaterinburg. "
            "DNA identification of the remains in 1994 was a milestone in forensic "
            "genetics. The mtDNA T2 was confirmed by matching with Prince Philip of "
            "the United Kingdom (great-grandson of Tsarina Alexandra's sister through "
            "the maternal line). The 2009 results (Rogaev et al.) confirmed the "
            "children's remains, including Anastasia - definitively debunking Anna "
            "Anderson's fraud."
        ),
    },
    {
        "name": "Otzi the Iceman",
        "dates": "c. 3300 BC",
        "title": "Natural mummy from the Italian Alps - 5,300 years old",
        "category": "ancient",
        "bio": (
            "The 5,300-year-old natural mummy found in the Italian Alps in 1991 - "
            "the oldest human with a fully published genome."
        ),
        "origin": "Neolithic farmer (Alps between Italy and Austria)",
        "y_hg": "G2a2b",
        "mt_hg": "K1f",
        "y_matches": [
            ("G2a2b", 5, "Match in G2a2b clade - lineage of early farmers"),
            ("G2a2", 4, "G2a2 branch - Neolithic European farmers"),
            ("G2a", 5, "G2a clade - haplogroup of Europe's first farmers"),
            ("G2", 4, "G2 clade"),
            ("G", 3, "Macrohaplogroup G - Caucasus/Anatolian origin"),
        ],
        "mt_matches": [
            ("K1f", 5, "Exact match - K1f (Otzi's precise haplogroup)"),
            ("K1", 4, "K1 clade"),
            ("K", 4, "Macrohaplogroup K - associated with Neolithic farmers and Ashkenazi Jews"),
            ("UK", 2, "U/K superlineage"),
        ],
        "evidence_y": 5,
        "evidence_mt": 5,
        "source": "Keller et al. 2012, Nature Communications (full genome published)",
        "verified": True,
        "has_otzi_compare": True,
        "story": (
            "Otzi is the best-preserved natural mummy in history, found in the Otztal "
            "Alps in 1991. His complete genome was the first prehistoric human genome "
            "to be published (2012). He was a Neolithic farmer of Anatolian descent - "
            "dark skin, brown eyes, lactose intolerant, carrier of atherosclerosis and "
            "Borrelia (Lyme disease). In 2016, living descendants were identified in "
            "Austrian Tyrol by sharing the same G2a haplotype (Niederstatter et al.) - "
            "making Otzi the oldest historical figure with verifiable autosomal kinship."
        ),
    },
    {
        "name": "Ludwig van Beethoven",
        "dates": "1770-1827 AD",
        "title": "Composer - one of the greatest figures in Western music",
        "category": "artist",
        "bio": (
            "One of the greatest composers in Western music history. His hair was "
            "sequenced in 2023, revealing surprising genetic secrets."
        ),
        "origin": "German-Flemish, Bonn",
        "y_hg": "I1",
        "mt_hg": "H17",
        "y_matches": [
            ("I1", 5, "Exact match - I1 haplogroup (Germanic/Scandinavian)"),
            ("I", 3, "Macrohaplogroup I - European Paleolithic hunters"),
            ("IJ", 2, "IJ superlineage"),
        ],
        "mt_matches": [
            ("H17", 5, "Exact match - H17"),
            ("H1", 3, "H1 clade - close branch"),
            ("H", 4, "Macrohaplogroup H - the most common in Europe (~44%)"),
            ("HV", 3, "HV ancestor"),
        ],
        "evidence_y": 5,
        "evidence_mt": 4,
        "source": "Begg et al. 2023, Current Biology",
        "verified": True,
        "story": (
            "In 2023, five authenticated locks of Beethoven's hair were sequenced. "
            "Results: (1) Y haplogroup = I1 - typically Germanic/Scandinavian; (2) a "
            "paternity discontinuity was detected: Beethoven was likely not biologically "
            "descended from the official 'Beethoven' lineage; (3) genetic predisposition "
            "to liver disease - consistent with historical suspicions of cirrhosis as "
            "cause of death. Beethoven's deafness still has no identified genetic "
            "explanation."
        ),
    },
    {
        "name": "Napoleon Bonaparte",
        "dates": "1769-1821 AD",
        "title": "Emperor of France",
        "category": "conqueror",
        "bio": (
            "Emperor of France who conquered most of Europe. His Y-DNA was inferred "
            "from living paternal relatives."
        ),
        "origin": "Corsican (family of Italian-Genoese origin)",
        "y_hg": "E1b1b1c",
        "mt_hg": None,
        "y_matches": [
            ("E1b1b1c", 5, "Same haplogroup inferred for the Bonaparte lineage"),
            ("E1b1b1", 4, "E1b1b1 clade"),
            ("E1b1b", 3, "Macrohaplogroup E1b1b - Mediterranean/North African"),
            ("E1b1", 2, "E1b1 clade"),
            ("E", 1, "Macrohaplogroup E"),
        ],
        "mt_matches": [],
        "evidence_y": 4,
        "evidence_mt": 0,
        "source": "Hagelberg et al. 2011; Lucotte & Thomasset 2011 (paternal relatives)",
        "verified": False,
        "story": (
            "Napoleon's haplogroup was inferred from tests on genealogically confirmed "
            "paternal-line male relatives - not from Napoleon himself. E1b1b is common "
            "in Corsica (~20%) and southern Italy, consistent with the Corsican origin "
            "of the Bonapartes. This haplogroup reflects pre-Roman Mediterranean "
            "connections - it arrived in Europe mainly with Greek, Phoenician, and "
            "Berber migrations in the first millennium BC."
        ),
    },
    {
        "name": "Ramesses III",
        "dates": "c. 1186-1155 BC",
        "title": "Pharaoh of Egypt (20th Dynasty, New Kingdom)",
        "category": "ancient",
        "bio": (
            "Pharaoh of Egypt's 20th Dynasty who defended against the Sea Peoples "
            "invasion (~1178 BC)."
        ),
        "origin": "Ancient Egyptian",
        "y_hg": "E1b1a",
        "mt_hg": None,
        "y_matches": [
            ("E1b1a", 5, "Same reported haplogroup - E1b1a"),
            ("E1b1", 3, "E1b1 clade"),
            ("E", 2, "Macrohaplogroup E"),
        ],
        "mt_matches": [],
        "evidence_y": 2,
        "evidence_mt": 0,
        "source": "Hawass et al. 2012, JAMA (methodology contested by scientific community)",
        "verified": False,
        "story": (
            "WARNING: This result is strongly contested. The Hawass et al. (2012) study "
            "never underwent full peer review, raw data were not published, and "
            "independent researchers question both methodology and conclusions. E1b1a "
            "is predominantly sub-Saharan - if confirmed, it would suggest significant "
            "African ancestry in the 20th Dynasty pharaohs. Treat as hypothesis, not "
            "established fact."
        ),
    },
    {
        "name": "Cheddar Man",
        "dates": "c. 7150 BC",
        "title": "Mesolithic British hunter-gatherer (9,100 years old)",
        "category": "ancient",
        "bio": (
            "Britain's oldest complete skeleton (9,100 years old). Had dark skin and "
            "blue eyes - rewrote our understanding of European appearance."
        ),
        "origin": "WHG (Western Hunter-Gatherer), Somerset, England",
        "y_hg": "I2a",
        "mt_hg": "U5b1",
        "y_matches": [
            ("I2a", 5, "Exact match - I2a (WHG signature)"),
            ("I2", 4, "I2 clade - European Paleolithic hunters"),
            ("I", 3, "Macrohaplogroup I"),
            ("IJ", 2, "IJ superlineage"),
        ],
        "mt_matches": [
            ("U5b1", 5, "Exact match - U5b1 (WHG signature)"),
            ("U5b", 4, "U5b clade"),
            ("U5", 5, "U5 - Europe's oldest haplogroup (>30,000 years)"),
            ("U", 3, "Macrohaplogroup U"),
        ],
        "evidence_y": 5,
        "evidence_mt": 5,
        "source": "Brace et al. 2019, Nature Ecology & Evolution",
        "verified": True,
        "story": (
            "Cheddar Man is the best-preserved human skeleton in Britain. His results "
            "generated worldwide headlines: despite being a 'native' Briton living "
            "9,100 years ago in Somerset, he had dark skin, blue eyes, and dark curly "
            "hair. This demonstrates that light European skin is a recent evolutionary "
            "adaptation (last 6,000-8,000 years) - not an 'original' trait. Cheddar Man "
            "belonged to the Western Hunter-Gatherers (WHG), who inhabited pre-agricultural "
            "Europe and were largely replaced by Neolithic farmers from Anatolia around "
            "6,000 BC."
        ),
    },
    {
        "name": "Loschbour Man",
        "dates": "c. 6000 BC",
        "title": "Mesolithic hunter-gatherer, Luxembourg",
        "category": "ancient",
        "bio": (
            "A Mesolithic hunter-gatherer from Luxembourg (~8,000 years old). His "
            "genome defines the Western Hunter-Gatherer reference profile."
        ),
        "origin": "WHG (Western Hunter-Gatherer), Western Europe",
        "y_hg": "I2a",
        "mt_hg": "U5b1a",
        "y_matches": [
            ("I2a", 5, "Exact match - I2a WHG"),
            ("I2", 4, "I2 clade"),
            ("I", 3, "Macrohaplogroup I"),
        ],
        "mt_matches": [
            ("U5b1a", 5, "Exact match - U5b1a"),
            ("U5b1", 4, "U5b1 clade"),
            ("U5b", 4, "U5b clade"),
            ("U5", 4, "U5 - Europe's oldest haplogroup"),
            ("U", 3, "Macrohaplogroup U"),
        ],
        "evidence_y": 5,
        "evidence_mt": 5,
        "source": "Lazaridis et al. 2014, Nature",
        "verified": True,
        "story": (
            "Loschbour Man is one of the fundamental references for Western "
            "Hunter-Gatherers (WHG). His genome was sequenced by Lazaridis et al. "
            "(2014) - the study that established the three ancestral components of "
            "modern Europe (WHG, EEF, Steppe). Like Cheddar Man, he had dark skin "
            "and light eyes. His DNA defines the 'WHG profile' used as reference in "
            "all ancient European ancestry studies."
        ),
    },
    {
        "name": "Mongol Imperial Clan (Genghis Khan lineage)",
        "dates": "c. 1162-1227 AD",
        "title": "Founder of the largest contiguous empire in history",
        "category": "conqueror",
        "bio": (
            "Founder of the Mongol Empire, the largest contiguous empire in history. "
            "His tomb has never been found."
        ),
        "origin": "Mongolia (Borjigin clan)",
        "y_hg": "C2b",
        "mt_hg": None,
        "y_matches": [
            ("C2b", 5, "Same C2b branch (star cluster - Mongol imperial lineage)"),
            ("C2", 4, "C2 clade"),
            ("C", 3, "Macrohaplogroup C - Mongolic/Siberian/Oceanic"),
        ],
        "mt_matches": [],
        "evidence_y": 3,
        "evidence_mt": 0,
        "source": "Zerjal et al. 2003, AJHG (population study - no direct Genghis Khan DNA)",
        "verified": False,
        "story": (
            "CRITICAL NOTE: No confirmed DNA exists from Genghis Khan - his tomb has "
            "never been found. The 'C2b star cluster' was inferred by Zerjal et al. "
            "(2003): ~8% of men in the region from Mongolia to the Pacific share a Y "
            "branch that expanded explosively ~1,000 years ago - coinciding in time and "
            "place with the Mongol Empire. Estimated: ~16 million living men would be "
            "descendants of the Borjigin clan. What we can say: you share a Y lineage "
            "with the estimated Mongol imperial family - nothing more."
        ),
    },
    {
        "name": "Johann Sebastian Bach",
        "dates": "1685-1750 AD",
        "title": "Composer - father of Western classical music",
        "category": "artist",
        "bio": (
            "Father of Western classical music. His Y-DNA was inferred from living "
            "male Bach family members."
        ),
        "origin": "German, Eisenach, Thuringia",
        "y_hg": "I1",
        "mt_hg": None,
        "y_matches": [
            ("I1", 5, "Same inferred haplogroup - I1 (Germanic/Scandinavian)"),
            ("I", 3, "Macrohaplogroup I"),
        ],
        "mt_matches": [],
        "evidence_y": 3,
        "evidence_mt": 0,
        "source": "Bach Family Y-DNA genealogical project (no direct DNA from JSB)",
        "verified": False,
        "story": (
            "The I1 haplogroup for Bach was inferred by Y-DNA testing of living male "
            "Bach family members with documented paternal lineage - not from the "
            "composer himself. I1 is the quintessential Scandinavian/Germanic "
            "haplogroup, consistent with the Bach family origin in central Thuringia "
            "(a region with high I1 frequency). Interestingly, Beethoven is also I1 - "
            "but this does not imply close kinship, as I1 encompasses ~16% of German men."
        ),
    },
    {
        "name": "Tutankhamun and the 18th Egyptian Dynasty",
        "dates": "c. 1341-1323 BC",
        "title": "Pharaoh of Egypt - New Kingdom",
        "category": "ancient",
        "bio": (
            "The boy king of Egypt's 18th Dynasty, famous for his intact tomb "
            "discovered by Howard Carter in 1922."
        ),
        "origin": "Ancient Egyptian (Upper Egypt)",
        "y_hg": "R1b",
        "mt_hg": "K",
        "y_matches": [
            ("R1b-U152", 2, "Distinct R1b branches - common ancestor >10,000 years"),
            ("R1b-P312", 2, "Distinct R1b branches"),
            ("R1b-L11", 2, "Distinct R1b branches"),
            ("R1b", 4, "Same macrohaplogroup R1b - broad ancestral lineage"),
            ("R1", 2, "R1 superlineage"),
        ],
        "mt_matches": [
            ("K1f", 4, "Close K branch (K1f is Otzi's)"),
            ("K", 5, "Exact macrohaplogroup match - mtDNA K of the 18th Dynasty"),
        ],
        "evidence_y": 3,
        "evidence_mt": 3,
        "source": "Gad, Hawass et al. 2021, Human Molecular Genetics 30(R1):R24-R28",
        "verified": True,
        "story": (
            "Gad & Hawass (2021) published haplogroups for Tutankhamun, Amenhotep III, "
            "and Akhenaten: all carry Y-DNA R1b (R-M343, the ancestral root) and mtDNA K. "
            "IMPORTANT: R-M343 is the ancestral trunk of R1b - present in North Africa "
            "and the Middle East for at least 20,000 years, long before being associated "
            "with Europe. It does not imply European ancestry. The specific subclade was "
            "not determined. Caveats: raw data were never deposited in a public "
            "repository; no independent replication."
        ),
    },
    {
        "name": "Nicolaus Copernicus",
        "dates": "1473-1543 AD",
        "title": "Astronomer - founder of the heliocentric model",
        "category": "scientist",
        "bio": (
            "Polish astronomer who proposed the heliocentric model. His skeleton was "
            "identified by DNA from hair in his own book."
        ),
        "origin": "Polish (Frombork, Royal Prussia)",
        "y_hg": None,
        "mt_hg": "H27",
        "y_matches": [],
        "mt_matches": [
            ("H27", 5, "Exact match - H27 (confirmed identification)"),
            ("H", 4, "Macrohaplogroup H - most common in Europe (~44%)"),
            ("HV", 3, "HV ancestor - root of H and V"),
        ],
        "evidence_y": 0,
        "evidence_mt": 4,
        "source": "Gierczyński et al. 2009, PNAS (skeleton + hair from authenticated book)",
        "verified": True,
        "story": (
            "The skeleton exhumed in 2005 from Frombork Cathedral (where Copernicus "
            "served as canon for decades) was identified by DNA in 2009 with elegant "
            "methodology: the same mtDNA H27 profile - with two private mutations "
            "defining the subclade - was found in both the skeleton's molars/femur and "
            "two hairs recovered from a calendar belonging to Copernicus, held at the "
            "Museum Gustavianum in Uppsala, Sweden."
        ),
    },
    {
        "name": "Marie Antoinette",
        "dates": "1755-1793 AD",
        "title": "Queen of France (Habsburg-Lorraine)",
        "category": "royalty",
        "bio": (
            "Queen of France, executed during the French Revolution. DNA extracted "
            "from a hair locket preserved by her mother."
        ),
        "origin": "Austrian, daughter of Empress Maria Theresa",
        "y_hg": None,
        "mt_hg": "H",
        "y_matches": [],
        "mt_matches": [
            ("H27", 3, "Close H branch (Copernicus is also H)"),
            ("H17", 3, "Close H branch (Beethoven is H17)"),
            ("H1", 3, "Close H branch"),
            ("H", 4, "Match in macrohaplogroup H - Habsburg-Lorraine lineage"),
            ("HV", 3, "HV ancestor of H and V"),
        ],
        "evidence_y": 0,
        "evidence_mt": 3,
        "source": "Jehaes et al. 1998, Eur J Hum Genet; refined in 2020, Int J Sciences",
        "verified": False,
        "story": (
            "DNA was extracted from a hair locket preserved by Empress Maria Theresa "
            "(Marie Antoinette's mother). The heart of Louis XVII (her son) - which "
            "survived decades in family custody before being analyzed - matched the "
            "maternal profile, confirming the lineage. Haplogroup H with specific "
            "private mutations (152C, 194T, 263G, 315.1C in HVR2) were identified. "
            "Caveats: the studies predate modern aDNA authentication standards."
        ),
    },
    {
        "name": "Rurikid Dynasty - Dmitry Alexandrovich",
        "dates": "9th-13th century AD",
        "title": "Founding dynasty of Russia (son of Alexander Nevsky)",
        "category": "royalty",
        "bio": (
            "The founding dynasty of Russia, from Rurik (830 AD) through Alexander "
            "Nevsky to Ivan the Terrible."
        ),
        "origin": "Scandinavian-Slavic (Novgorod/Kiev)",
        "y_hg": "N1a",
        "mt_hg": None,
        "y_matches": [
            ("N1a", 5, "Exact match - N1a (confirmed Rurikid haplogroup)"),
            ("N", 3, "Macrohaplogroup N - Uralic/Finno-Ugric/Siberian"),
        ],
        "mt_matches": [],
        "evidence_y": 4,
        "evidence_mt": 0,
        "source": "Paleogenomic study 2023, PMC10615192 (Dmitry Alexandrovich, son of Alexander Nevsky)",
        "verified": True,
        "story": (
            "Ancient DNA from Prince Dmitry Alexandrovich (d. 1294), son of Alexander "
            "Nevsky, confirmed Y-DNA N1a - consistent with living male Rurikid "
            "descendants today (N1c1-L550/Y4343). The Rurikids were Russia's founding "
            "dynasty, descended from Rurik (c. 830 AD), a Varangian (Scandinavian "
            "Viking) leader who ruled Novgorod. N1a is predominantly Uralic/Finno-Ugric "
            "- its presence in the Varangian-Slavic elite is unexpected and suggests "
            "the Rurikid founders carried paternal lineage from Euro-Siberian populations "
            "of northern Scandinavia."
        ),
    },
    {
        "name": "Corvinus - Hunyadi Family",
        "dates": "15th century AD",
        "title": "John and Christopher Corvinus (sons of King Matthias Corvinus of Hungary)",
        "category": "royalty",
        "bio": (
            "Hungary's most powerful medieval royal family. Sons of King Matthias "
            "Corvinus had their full genomes sequenced."
        ),
        "origin": "Hungarian (Wallachian of Transylvanian origin)",
        "y_hg": "E1b1b",
        "mt_hg": "T2",
        "y_matches": [
            ("E1b1b1c", 4, "Close E1b1b branch (Napoleon is also E1b1b)"),
            ("E1b1b", 5, "Match in E1b1b clade - Corvinus/Hunyadi lineage"),
            ("E1b1", 3, "E1b1 clade"),
            ("E", 2, "Macrohaplogroup E"),
        ],
        "mt_matches": [
            ("T2", 5, "Exact match - T2 (same as the Romanovs!)"),
            ("T", 4, "T clade - Anatolian origin"),
        ],
        "evidence_y": 4,
        "evidence_mt": 4,
        "source": "Neparaczki et al. 2022, Heliyon (complete genomes from two confirmed sons of Matthias Corvinus)",
        "verified": True,
        "story": (
            "Matthias Corvinus (1443-1490) was one of Hungary's greatest medieval kings. "
            "Neparaczki et al. (2022) sequenced the complete genomes of John Corvinus "
            "and Christopher Corvinus - confirmed illegitimate sons of Matthias. Result: "
            "Y-DNA E1b1b-V13, mtDNA T2b (John) and T2c1+146 (Christopher). E1b1b-V13 "
            "is common in the Balkans and Romania (~15-20%), consistent with the "
            "Transylvanian origin of the Hunyadis."
        ),
    },
    {
        "name": "Birger Jarl",
        "dates": "c. 1210-1266 AD",
        "title": "Founder of Stockholm and regent of Sweden",
        "category": "royalty",
        "bio": (
            "Founder of Stockholm (~1252 AD) and architect of Sweden as a feudal "
            "kingdom."
        ),
        "origin": "Swedish (Bjalbo clan)",
        "y_hg": "I1",
        "mt_hg": None,
        "y_matches": [
            ("I1", 5, "Exact match - I1 (the quintessential Viking haplogroup)"),
            ("I", 3, "Macrohaplogroup I - European Paleolithic hunters"),
        ],
        "mt_matches": [],
        "evidence_y": 4,
        "evidence_mt": 0,
        "source": "National Board of Forensic Medicine, Sweden (identified skeletal remains)",
        "verified": True,
        "story": (
            "Birger Jarl was the statesman who founded Stockholm (~1252 AD) and "
            "transformed Sweden into an organized feudal kingdom. His Y-DNA I1 was "
            "confirmed from skeletal remains by Sweden's National Board of Forensic "
            "Medicine. I1 (I-M253) is the dominant haplogroup in Scandinavia (30-40% "
            "of Swedes, Danes, and Norwegians) and the most common in Viking Age "
            "individuals."
        ),
    },
    {
        "name": "Hunnic Elite - Lineage Associated with Attila",
        "dates": "4th-5th century AD",
        "title": "Huns - conquerors from Central Asia to Western Europe",
        "category": "conqueror",
        "bio": (
            "The Hunnic warrior elite who swept from Central Asia to Europe in the "
            "5th century AD."
        ),
        "origin": "Eurasian (Xiongnu origin -> Central Asia -> Carpathian Basin)",
        "y_hg": "R1a",
        "mt_hg": None,
        "y_matches": [
            ("R1a", 4, "Match in R1a clade (R1a-Z93 - central Hunnic haplotype)"),
            ("R1", 2, "R1 superlineage"),
        ],
        "mt_matches": [],
        "evidence_y": 2,
        "evidence_mt": 0,
        "source": "Gnecchi-Ruscone et al. 2025, PNAS; Maroti et al. 2022, Current Biology",
        "verified": False,
        "story": (
            "NOTE: Attila was never found - his tomb remains unknown. This profile "
            "represents the 5th-century Hunnic elite, based on ancient DNA from Hunnic "
            "cemeteries in the Carpathian Basin (Hungary/Romania). Dominant haplogroups "
            "in the Hunnic elite are R1a-Z93 (~43%) and Q-M242 (~39%). R1a-Z93 is also "
            "the dominant haplogroup among Scythians and Indo-Iranians of the steppes - "
            "it is not exclusive to the Huns."
        ),
    },
    {
        "name": "Charlemagne (Carolingians)",
        "dates": "747-814 AD",
        "title": "Emperor of the Holy Roman Empire - unified Western Europe",
        "category": "royalty",
        "bio": (
            "Emperor who united Western Europe under the Carolingian dynasty (~800 AD)."
        ),
        "origin": "Frankish (Rhine region, modern Germany/Belgium)",
        "y_hg": "R1b",
        "mt_hg": None,
        "y_matches": [
            ("R1b-U152", 3, "Speculated haplogroup for the Carolingian/Habsburg lineage"),
            ("R1b-P312", 3, "R1b branch - Frankish/Carolingian ancestor"),
            ("R1b-L11", 3, "R1b-L11 clade - Western Europe"),
            ("R1b", 2, "Macrohaplogroup R1b"),
            ("R1", 1, "R1 superlineage"),
        ],
        "mt_matches": [],
        "evidence_y": 2,
        "evidence_mt": 0,
        "source": "Inferred from living descendants of the Habsburg-Lorraine lineage (no direct aDNA from Aachen relics)",
        "verified": False,
        "story": (
            "WARNING: No ancient DNA has been extracted from Charlemagne's relics at "
            "Aachen Cathedral. The R1b-U152/L2 haplogroup is speculative, inferred from "
            "~20 men with the Habsburg surname who share a close haplotype. Problems: "
            "(1) the unbroken male line from Charlemagne likely went extinct within a "
            "few generations; (2) mathematical analyses (Ralph & Coop 2013) show that "
            "virtually all Europeans descend from Charlemagne through SOME line - but "
            "Y-DNA requires an unbroken paternal line, which likely broke centuries ago."
        ),
    },
    {
        "name": "Galilean Jewish Men - 1st century AD (historical context)",
        "dates": "c. 1-100 AD",
        "title": "Population reference - Galilee and Second Temple Judea",
        "category": "historical",
        "bio": (
            "Population reference for 1st century AD Galilean Jewish men, based on "
            "ancient DNA from contemporary Levantine sites."
        ),
        "origin": "Ancient Levant (northern Israel/Palestine)",
        "y_hg": "J1",
        "mt_hg": None,
        "y_matches": [
            ("J1", 5, "Most common haplogroup (~30-40%) in 1st century Galilean Jewish men"),
            ("J2", 4, "Second most common haplogroup (~15-20%) in Galilean Jews"),
            ("J", 4, "J clade - predominant in the ancient Levant"),
        ],
        "mt_matches": [],
        "evidence_y": 3,
        "evidence_mt": 0,
        "source": "Haber et al. 2020, Current Biology; Harney et al. 2018, Nature Communications; Feldman et al. 2019, Nature Communications",
        "verified": False,
        "story": (
            "POPULATION REFERENCE - NOT AN IDENTIFIED INDIVIDUAL. This entry represents "
            "the most likely Y haplogroups for 1st century AD Galilean Jewish men, based "
            "on ancient DNA from contemporary archaeological sites in the Levant. "
            "Documented haplogroups: J1-P58 (~30-40%), J2a-M410 (~15-20%), E1b1b "
            "(~10-15%), T1a (~5-8%), G2a (~3-6%). No confirmed DNA exists from any "
            "specific historical figure of this period and region."
        ),
    },
    # ── INFAMOUS ───────────────────────────────────────────────────────────────
    {
        "name": "Adolf Hitler",
        "dates": "1889-1945 AD",
        "title": "Dictator of Nazi Germany",
        "category": "infamous",
        "bio": (
            "Leader of Nazi Germany responsible for the Holocaust and World War II. "
            "His Y-DNA was determined from 39 patrilineal relatives in 2010."
        ),
        "origin": "Austrian (Braunau am Inn)",
        "y_hg": "E1b1b",
        "mt_hg": None,
        "y_matches": [
            ("E1b1b1c", 4, "Close E1b1b branch (also shared by Napoleon)"),
            ("E1b1b", 5, "Exact match - E1b1b (E-M35)"),
            ("E1b1", 3, "E1b1 clade"),
            ("E", 2, "Macrohaplogroup E"),
        ],
        "mt_matches": [],
        "evidence_y": 4,
        "evidence_mt": 0,
        "source": "Mulders & Vermeeren 2010, Knack magazine (39 patrilineal relatives tested)",
        "verified": False,
        "story": (
            "In 2010, journalist Jean-Paul Mulders and historian Marc Vermeeren collected "
            "saliva samples from 39 patrilineal relatives of Hitler across Austria and the "
            "United States. All carried Y-DNA haplogroup E1b1b (E-M35). Ironically, E1b1b "
            "is most common in North Africa, the Horn of Africa, and among some Jewish "
            "populations - contradicting Nazi racial ideology. However, Professor Michael "
            "Hammer cautioned that no conclusions about Jewish ancestry can be drawn from "
            "haplogroup alone, as E1b1b is widespread across many populations."
        ),
    },
    {
        "name": "Joseph Stalin",
        "dates": "1878-1953 AD",
        "title": "Dictator of the Soviet Union",
        "category": "infamous",
        "bio": (
            "General Secretary of the Soviet Communist Party who ruled through terror. "
            "Responsible for the deaths of millions through purges, gulags, and famine."
        ),
        "origin": "Georgian (Gori, Georgia)",
        "y_hg": "G2a",
        "mt_hg": None,
        "y_matches": [
            ("G2a2b", 4, "Close G2a branch (Otzi is also G2a2b)"),
            ("G2a", 5, "Exact match - G2a (G-FGC595/Z6553)"),
            ("G2", 4, "G2 clade"),
            ("G", 3, "Macrohaplogroup G - Caucasus/European Neolithic"),
        ],
        "mt_matches": [],
        "evidence_y": 4,
        "evidence_mt": 0,
        "source": "DNA test on grandson Alexander Burdonsky (son of Vasily Stalin)",
        "verified": False,
        "story": (
            "Stalin's Y-DNA haplogroup G2a was determined through testing of his grandson "
            "Alexander Burdonsky, son of Vasily Stalin. G2a1a (G-FGC595/Z6553) is found at "
            "high frequencies in the Caucasus Mountains, consistent with Stalin's Georgian "
            "heritage. The Jughashvili family originated from the village of Didi Lilo near "
            "Tbilisi. G2a is also the haplogroup of Otzi the Iceman and Richard III - "
            "connecting a 20th-century dictator to a 5,300-year-old mummy and a medieval king "
            "through deep patrilineal ancestry."
        ),
    },
    # ── ROYALTY (EXPANDED) ─────────────────────────────────────────────────────
    {
        "name": "House of Habsburg",
        "dates": "13th-20th century AD",
        "title": "Europe's most powerful dynasty for 600 years",
        "category": "royalty",
        "bio": (
            "The Habsburg dynasty ruled the Holy Roman Empire, Spain, Austria, and much of "
            "Europe for over six centuries. Known for strategic marriages and genetic effects "
            "of consanguinity."
        ),
        "origin": "Swiss-Austrian (Habsburg Castle, Aargau, Switzerland)",
        "y_hg": "R1b",
        "mt_hg": None,
        "y_matches": [
            ("R1b-U152", 5, "Exact match - R1b-U152 (L2+ branch, Habsburg signature)"),
            ("R1b-P312", 4, "R1b-P312 clade"),
            ("R1b-L11", 3, "R1b-L11 branch"),
            ("R1b", 3, "Macrohaplogroup R1b"),
            ("R1", 2, "R1 superlineage"),
        ],
        "mt_matches": [],
        "evidence_y": 4,
        "evidence_mt": 0,
        "source": "FamilyTreeDNA Habsburg Family Project (20 documented members tested)",
        "verified": False,
        "story": (
            "The Habsburg Y-DNA was determined through the FamilyTreeDNA Habsburg Family "
            "Project, where 20 living members with documented Habsburg surname all share "
            "R1b-U152 (L2+ branch). This subclade is consistent with the family's origin in "
            "northern Switzerland, a known hotspot for U152. The Habsburgs are infamous for "
            "extensive inbreeding - Charles II of Spain (1661-1700) had an inbreeding "
            "coefficient higher than a child of siblings, contributing to the dynasty's "
            "eventual biological decline."
        ),
    },
    {
        "name": "House of Bourbon - Louis XVI",
        "dates": "16th-21st century AD",
        "title": "Royal house of France, Spain, and Luxembourg",
        "category": "royalty",
        "bio": (
            "The Bourbon dynasty ruled France from Henry IV to Louis XVI, and continues "
            "in Spain today. Their Y-DNA was confirmed from three living descendants."
        ),
        "origin": "French (House of Capet branch)",
        "y_hg": "R1b",
        "mt_hg": None,
        "y_matches": [
            ("R1b-U152", 3, "Distinct R1b branches"),
            ("R1b-P312", 3, "Distinct R1b branches"),
            ("R1b-L11", 3, "R1b-L11 branch"),
            ("R1b", 4, "Same macrohaplogroup R1b (R-Z381* subclade)"),
            ("R1", 2, "R1 superlineage"),
        ],
        "mt_matches": [],
        "evidence_y": 5,
        "evidence_mt": 0,
        "source": "Larmuseau et al. 2014, European Journal of Human Genetics",
        "verified": True,
        "story": (
            "A 2014 peer-reviewed study by Larmuseau et al. tested three living male-line "
            "Bourbon descendants and confirmed Y-DNA R1b (R-Z381*). Importantly, this "
            "contradicted a previous claim that a blood-stained handkerchief from Louis XVI's "
            "execution carried haplogroup G2a - proving the relic was likely inauthentic. The "
            "Bourbons descend from the House of Capet, making their Y-DNA the confirmed "
            "patrilineal signature of the French royal house from Hugh Capet (987 AD) onward."
        ),
    },
    {
        "name": "Queen Victoria",
        "dates": "1819-1901 AD",
        "title": "Queen of the United Kingdom and Empress of India",
        "category": "royalty",
        "bio": (
            "Queen Victoria reigned for 63 years and became the 'grandmother of Europe' "
            "through her descendants' marriages into royal houses across the continent."
        ),
        "origin": "British-German (House of Hanover)",
        "y_hg": None,
        "mt_hg": "H",
        "y_matches": [],
        "mt_matches": [
            ("H27", 3, "Close H branch (Copernicus is also H)"),
            ("H17", 3, "Close H branch (Beethoven is H17)"),
            ("H1", 3, "Close H branch"),
            ("H", 5, "Exact match - Haplogroup H (confirmed via Prince Philip of Edinburgh)"),
            ("HV", 3, "HV ancestor"),
        ],
        "evidence_y": 0,
        "evidence_mt": 5,
        "source": "Confirmed via Prince Philip matching Romanov remains (Gill et al. 1994)",
        "verified": True,
        "story": (
            "Queen Victoria's mtDNA haplogroup H was confirmed through a remarkable chain of "
            "evidence: Prince Philip, Duke of Edinburgh, shared Victoria's maternal lineage "
            "(both descend from Princess Alice of the UK through the female line). When Prince "
            "Philip's mtDNA was used to identify the Romanov remains in 1994, it also "
            "confirmed the mtDNA lineage of Victoria herself. Victoria's maternal line "
            "connects to virtually every surviving European royal house."
        ),
    },
    # ── SCIENTISTS ─────────────────────────────────────────────────────────────
    {
        "name": "Charles Darwin",
        "dates": "1809-1882 AD",
        "title": "Naturalist - author of On the Origin of Species",
        "category": "scientist",
        "bio": (
            "Father of the theory of evolution by natural selection. His great-great-grandson "
            "was DNA tested to determine the Darwin patrilineal haplogroup."
        ),
        "origin": "English (Shrewsbury, Shropshire)",
        "y_hg": "R1b",
        "mt_hg": None,
        "y_matches": [
            ("R1b-U152", 3, "Distinct R1b branches"),
            ("R1b-P312", 3, "Distinct R1b branches"),
            ("R1b-L11", 3, "R1b-L11 branch"),
            ("R1b", 4, "Same macrohaplogroup R1b"),
            ("R1", 2, "R1 superlineage"),
        ],
        "mt_matches": [],
        "evidence_y": 4,
        "evidence_mt": 0,
        "source": "Great-great-grandson DNA tested via FamilyTreeDNA",
        "verified": False,
        "story": (
            "Darwin's Y-DNA R1b was determined by testing a great-great-grandson through the "
            "paternal line. R1b is the most common Y-DNA haplogroup in Western Europe (~60-70% "
            "of men in the British Isles), consistent with the Darwin family's deep English "
            "roots in Shropshire. Ironically, the man who revolutionized our understanding of "
            "heredity and descent carries the most common Western European patrilineal lineage "
            "- a haplogroup that expanded massively during the Bronze Age with Indo-European "
            "migrations ~4,500 years ago."
        ),
    },
    {
        "name": "Nikola Tesla",
        "dates": "1856-1943 AD",
        "title": "Inventor and electrical engineer",
        "category": "scientist",
        "bio": (
            "Serbian-American inventor who pioneered alternating current (AC) electrical "
            "systems, the Tesla coil, and wireless energy transmission."
        ),
        "origin": "Serbian (Smiljan, Austrian Empire - modern Croatia)",
        "y_hg": "R1a",
        "mt_hg": None,
        "y_matches": [
            ("R1a", 5, "Exact match - R1a (Slavic/Indo-European signature)"),
            ("R1", 2, "R1 superlineage"),
        ],
        "mt_matches": [],
        "evidence_y": 3,
        "evidence_mt": 0,
        "source": "Serbian DNA Project (relative/descendant testing)",
        "verified": False,
        "story": (
            "Tesla's Y-DNA R1a-M458 (L1029 subclade) was determined through the Serbian DNA "
            "Project via relative testing. R1a is the most common haplogroup among Slavic "
            "peoples (~40-50% in Serbia, Poland, Russia) and traces back to Indo-European "
            "migrations from the Pontic Steppe ~5,000 years ago. Tesla's R1a-M458 subclade "
            "is specifically associated with West Slavic populations and is one of the "
            "largest Slavic-specific branches of R1a."
        ),
    },
    # ── POLITICAL LEADERS ──────────────────────────────────────────────────────
    {
        "name": "Abraham Lincoln",
        "dates": "1809-1865 AD",
        "title": "16th President of the United States",
        "category": "political",
        "bio": (
            "President who preserved the Union during the American Civil War and abolished "
            "slavery. His haplogroup was determined through multiple paternal descendants."
        ),
        "origin": "American (English ancestry, Kentucky)",
        "y_hg": "R1b",
        "mt_hg": None,
        "y_matches": [
            ("R1b-U152", 5, "Exact match - R1b-U152 (Lincoln lineage confirmed)"),
            ("R1b-P312", 4, "R1b-P312 clade"),
            ("R1b-L11", 3, "R1b-L11 branch"),
            ("R1b", 3, "Macrohaplogroup R1b"),
            ("R1", 2, "R1 superlineage"),
        ],
        "mt_matches": [],
        "evidence_y": 4,
        "evidence_mt": 0,
        "source": "Lincoln DNA Project at FamilyTreeDNA (multiple descendants of Samuel 'the weaver' Lincoln)",
        "verified": False,
        "story": (
            "Lincoln's Y-DNA was determined through the Lincoln DNA Project at FamilyTreeDNA, "
            "where multiple paternal descendants of Samuel 'the weaver' Lincoln (Abraham's "
            "earliest documented ancestor, who emigrated from Norfolk, England to Massachusetts "
            "in 1637) were tested. Result: R1b-U152 (R1b-S20376, within U152>L2>Z142>Z150). "
            "Remarkably, this places Lincoln in the same R1b-U152 branch as the Habsburg "
            "dynasty and possibly George Washington - a ~4,500-year-old lineage connecting "
            "an American president to European royal houses."
        ),
    },
    {
        "name": "Thomas Jefferson",
        "dates": "1743-1826 AD",
        "title": "3rd President of the United States",
        "category": "political",
        "bio": (
            "Principal author of the Declaration of Independence. DNA testing of his "
            "descendants famously confirmed the Sally Hemings paternity controversy."
        ),
        "origin": "American (Welsh/English ancestry, Virginia)",
        "y_hg": "T",
        "mt_hg": None,
        "y_matches": [],
        "mt_matches": [],
        "evidence_y": 5,
        "evidence_mt": 0,
        "source": "Foster et al. 1998, Nature",
        "verified": True,
        "story": (
            "Jefferson's Y-DNA haplogroup T (T-M184, formerly K2) was determined through "
            "direct male-line descendant testing, published in Nature in 1998. This landmark "
            "study also confirmed that Jefferson fathered children with Sally Hemings, an "
            "enslaved woman. Haplogroup T is extremely rare in Europe (<1% of men), making "
            "it a powerful forensic identifier. The same rare haplogroup was found in both "
            "the Jefferson and Hemings descendants, with odds of a coincidental match "
            "effectively zero. NOTE: Your DNA chip may not detect haplogroup T markers, "
            "so this connection may not appear even if you carry haplogroup T."
        ),
    },
    # ── WARRIORS & CULTURAL FIGURES ────────────────────────────────────────────
    {
        "name": "Somerled - King of the Isles",
        "dates": "c. 1100-1164 AD",
        "title": "Norse-Gaelic King of the Hebrides and Argyll",
        "category": "conqueror",
        "bio": (
            "Norse-Gaelic leader who drove the Vikings from western Scotland and founded "
            "Clan Donald (MacDonald), the most powerful Highland clan."
        ),
        "origin": "Norse-Gaelic (Scotland/Hebrides)",
        "y_hg": "R1a",
        "mt_hg": None,
        "y_matches": [
            ("R1a", 5, "Exact match - R1a (Norse-Gaelic warrior lineage)"),
            ("R1", 2, "R1 superlineage"),
        ],
        "mt_matches": [],
        "evidence_y": 3,
        "evidence_mt": 0,
        "source": "Oxford University 2003; Bryan Sykes 2005",
        "verified": False,
        "story": (
            "Somerled is the progenitor of Clan Donald (MacDonald, MacDougall, MacAlister) - "
            "the most powerful clan in Scottish history. His Y-DNA R1a1 was determined by "
            "testing clan descendants, estimated at ~500,000 living men today. R1a is unusual "
            "for a 'Celtic' leader, suggesting Norse paternal ancestry - consistent with his "
            "Norse-Gaelic (Gall-Ghaidheal) heritage. Somerled expelled the Vikings from the "
            "Hebrides around 1156 and established a maritime kingdom that lasted centuries."
        ),
    },
    {
        "name": "Niall of the Nine Hostages",
        "dates": "c. 400 AD",
        "title": "Legendary High King of Ireland",
        "category": "conqueror",
        "bio": (
            "Semi-legendary Irish High King whose descendants became the most powerful "
            "dynasty in Ireland for a millennium."
        ),
        "origin": "Irish (Connacht/Tara)",
        "y_hg": "R1b",
        "mt_hg": None,
        "y_matches": [
            ("R1b-U152", 3, "Distinct R1b branches"),
            ("R1b-P312", 4, "R1b-P312 clade (Niall is within P312)"),
            ("R1b-L11", 3, "R1b-L11 branch"),
            ("R1b", 3, "Macrohaplogroup R1b"),
            ("R1", 2, "R1 superlineage"),
        ],
        "mt_matches": [],
        "evidence_y": 3,
        "evidence_mt": 0,
        "source": "Moore et al. 2006, Trinity College Dublin (population study)",
        "verified": False,
        "story": (
            "A 2006 study from Trinity College Dublin found that R1b-M222 occurs at very "
            "high frequency (~21%) in northwest Ireland, with a common ancestor estimated "
            "at ~1,700 years ago - consistent with Niall's era. The association with Niall "
            "specifically is genealogical inference (Ui Neill dynasty), not confirmed ancient "
            "DNA. Today, an estimated 2-3 million men worldwide carry the M222 marker. "
            "NOTE: Consumer DNA chips typically cannot resolve R1b-M222 specifically, so "
            "this match is at the broader R1b level."
        ),
    },
    {
        "name": "Martin Luther",
        "dates": "1483-1546 AD",
        "title": "Theologian who sparked the Protestant Reformation",
        "category": "religious",
        "bio": (
            "German theologian whose Ninety-Five Theses (1517) triggered the Protestant "
            "Reformation, fundamentally reshaping Christianity and European history."
        ),
        "origin": "German (Eisleben, Saxony)",
        "y_hg": "I2a",
        "mt_hg": None,
        "y_matches": [
            ("I2a", 5, "Exact match - I2a (European Paleolithic hunter lineage)"),
            ("I2", 4, "I2 clade"),
            ("I", 3, "Macrohaplogroup I - European Paleolithic"),
            ("IJ", 2, "IJ superlineage"),
        ],
        "mt_matches": [],
        "evidence_y": 3,
        "evidence_mt": 0,
        "source": "Tested relatives via genealogical DNA projects",
        "verified": False,
        "story": (
            "Luther's Y-DNA I2a was determined through tested relatives with documented "
            "paternal lineage to the Luther family. I2a is one of Europe's oldest "
            "haplogroups, carried by Paleolithic hunter-gatherers including Cheddar Man "
            "(9,100 years old). Its presence in Luther connects the Protestant reformer "
            "to the oldest surviving male lineage in Europe - a lineage that was largely "
            "replaced by Neolithic farmers and Bronze Age steppe migrants but persists "
            "today at ~15-20% in central and southeastern Europe."
        ),
    },
    {
        "name": "Jesse James",
        "dates": "1847-1882 AD",
        "title": "American outlaw and bank robber",
        "category": "infamous",
        "bio": (
            "America's most famous outlaw. His identity was confirmed by DNA extracted "
            "from his exhumed remains in 1995."
        ),
        "origin": "American (Missouri, Scots-Irish ancestry)",
        "y_hg": None,
        "mt_hg": "T",
        "y_matches": [],
        "mt_matches": [
            ("T2", 4, "Close T branch (T2 is shared by the Romanovs)"),
            ("T", 5, "Exact match - mtDNA T (confirmed from direct remains)"),
            ("JT", 3, "JT ancestor"),
        ],
        "evidence_y": 0,
        "evidence_mt": 5,
        "source": "Stone et al. 2001, Journal of Forensic Sciences (teeth from 1995 exhumation)",
        "verified": True,
        "story": (
            "Jesse James' mtDNA was extracted from teeth recovered during the 1995 exhumation "
            "of his grave. The profile (T: 16126C, 16274A, 16294T, 16296T, 16304C) matched "
            "living maternal descendants of Jesse's sister Susan, confirming the remains were "
            "indeed James - debunking persistent conspiracy theories that he had faked his "
            "death. This was one of the first uses of ancient DNA to resolve a historical "
            "identity dispute in the United States."
        ),
    },
    {
        "name": "Fath Ali Shah Qajar",
        "dates": "1772-1834 AD",
        "title": "Shah of Persia (Qajar Dynasty)",
        "category": "royalty",
        "bio": (
            "Second ruler of the Qajar dynasty of Iran, known for his exceptionally long "
            "beard and hundreds of children."
        ),
        "origin": "Persian-Turkic (Tehran, Iran)",
        "y_hg": "J1",
        "mt_hg": None,
        "y_matches": [
            ("J1", 5, "Exact match - J1 (J-M267, Qajar dynasty lineage)"),
            ("J", 3, "Macrohaplogroup J - Semitic/Middle Eastern"),
        ],
        "mt_matches": [],
        "evidence_y": 4,
        "evidence_mt": 0,
        "source": "Descendant testing of multiple lineages from Fath Ali Shah's sons",
        "verified": False,
        "story": (
            "Fath Ali Shah's Y-DNA J1 (J-M267) was determined through testing descendants "
            "from multiple sons' lineages. Fath Ali Shah reportedly fathered over 100 sons "
            "and is believed to have approximately 5,000 descendants alive today. J1 is the "
            "most common haplogroup in the Arabian Peninsula and Levant (~40-70% of men in "
            "some regions). Its presence in the Turkic Qajar dynasty reflects the deep "
            "genetic connections between Persian, Turkic, and Semitic populations in the "
            "Middle East."
        ),
    },
]


# ---------------------------------------------------------------------------
# OTZI PHENOTYPIC COMPARISON
# ---------------------------------------------------------------------------

OTZI_GENOTYPES: Dict[str, Dict[str, str]] = {
    "rs12913832": {
        "gt": "GG",
        "trait": "Eye color (HERC2/OCA2)",
        "otzi_desc": "Brown/dark eyes - no blue eye allele",
        "evidence": "Strong (confirmed)",
    },
    "rs4988235": {
        "gt": "GG",
        "trait": "Lactase persistence (LCT)",
        "otzi_desc": "Lactose intolerant - ancestral pre-pastoral genotype",
        "evidence": "Strong (confirmed)",
    },
    "rs1426654": {
        "gt": "AA",
        "trait": "Skin pigmentation (SLC24A5)",
        "otzi_desc": "European light skin allele (EEF derived variant)",
        "evidence": "Strong (confirmed)",
    },
    "rs1805007": {
        "gt": "CC",
        "trait": "Red hair (MC1R R151C)",
        "otzi_desc": "No red hair - ancestral allele",
        "evidence": "Strong (confirmed)",
    },
    "rs1815739": {
        "gt": "TT",
        "trait": "Muscle fiber type (ACTN3 R577X)",
        "otzi_desc": "Endurance type (XX) - consistent with high-altitude hiking",
        "evidence": "Strong (confirmed)",
    },
    "rs1800562": {
        "gt": "GG",
        "trait": "Hemochromatosis HFE C282Y",
        "otzi_desc": "No C282Y variant - ancestral allele",
        "evidence": "Moderate",
    },
}

OTZI_TRAIT_DISCLAIMERS: Dict[str, str] = {
    "rs12913832": (
        "Eye color prediction has ~93% accuracy for brown and blue, but only ~65% "
        "for green and hazel - colors that emerge from rare variants not captured "
        "in consumer genotyping chips. If you have green or hazel eyes, this result "
        "may be incorrect."
    ),
    "rs4988235": (
        "Lactase persistence testing is highly reliable (>99% accuracy) as it "
        "depends on a single well-characterized variant."
    ),
    "rs1426654": (
        "This SNP explains most of the variation between European and non-European "
        "skin pigmentation but doesn't capture the full spectrum of skin tones "
        "within populations."
    ),
    "rs1805007": (
        "MC1R testing captures the most common red hair variant but misses rarer "
        "MC1R variants (R160W, D294H) that also contribute."
    ),
    "rs1815739": (
        "ACTN3 R577X is well-validated but muscle performance depends on hundreds "
        "of genes, training, and environment."
    ),
    "rs1800562": (
        "C282Y is the most clinically significant HFE variant. Homozygosity is "
        "required for disease risk."
    ),
}


# ---------------------------------------------------------------------------
# Y-DNA AND mtDNA HAPLOGROUP DETECTION
# ---------------------------------------------------------------------------

# Y-DNA markers: (rsid, derived_allele, clade_name, specificity)
# Higher specificity = more downstream (more specific subclade)
Y_MARKERS_ORDERED = [
    ("rs34126399", "A", "R1b-U152", 9),
    ("rs34276300", "A", "R1b-P312", 8),
    ("rs9786137",  "A", "R1b-L11",  7),
    ("rs16981293", "A", "R1b",      6),
    ("rs9786153",  "A", "R",        5),
    ("rs17222279", "A", "R1a",      7),
    ("rs34079747", "A", "I1",       7),
    ("rs34311906", "A", "I2",       7),
    ("rs41362250", "A", "I2a",      8),
    ("rs2032601",  "A", "G",        5),
    ("rs17690947", "A", "G2a",      7),
    ("rs41462431", "A", "E1b1b",    7),
    ("rs41352448", "A", "E",        5),
    ("rs13447352", "A", "J2",       7),
    ("rs9341296",  "A", "J1",       7),
    ("rs41352439", "A", "N",        5),
    ("rs41352441", "A", "C",        5),
    ("rs41352449", "A", "Q",        5),
]

# mtDNA position-based markers (rCRS coordinates)
MT_POS_MARKERS: Dict[int, Dict[str, str]] = {
    14766: {"T": "HV"},
    15326: {"G": "HV"},
    7028:  {"T": "H"},
    4580:  {"A": "V"},
    295:   {"T": "JT"},
    16069: {"C": "J"},
    16311: {"C": "T"},
    16224: {"C": "K"},
    16270: {"C": "U5"},
    16356: {"C": "U4"},
    119:   {"C": "W"},
    16129: {"A": "I"},
    16223: {"T": "L"},
    16126: {"C": "JT"},
}

MT_DECISION_TREE = [
    ({"H", "HV"},  "H",  "H"),
    ({"HV", "V"},  "V",  "V"),
    ({"HV"},       "HV", "HV"),
    ({"JT", "J"},  "J",  "J"),
    ({"JT", "T"},  "T",  "T"),
    ({"JT"},       "JT", "JT"),
    ({"K"},        "K",  "K"),
    ({"U5"},       "U5", "U5"),
    ({"U4"},       "U4", "U4"),
    ({"W"},        "W",  "W"),
    ({"I"},        "I",  "I"),
    ({"L"},        "L",  "L"),
]

# rsID-based backup markers for mtDNA
RSID_MT_MARKERS: Dict[str, Tuple[str, str]] = {
    "rs28358564": ("T", "H"),
    "rs41419549": ("C", "V"),
    "rs28357093": ("C", "J"),
    "rs41419546": ("A", "T"),
    "rs3135031":  ("T", "K"),
    "rs2853519":  ("T", "U5"),
    "rs41360374": ("G", "JT"),
    "rs28357376": ("C", "L"),
}


# ---------------------------------------------------------------------------
# SEX DETECTION
# ---------------------------------------------------------------------------

def _detect_is_male(variants: Dict[str, Tuple[str, str, str]]) -> bool:
    """Detect if the user is male by checking for Y chromosome markers."""
    for rsid, derived, clade, spec in Y_MARKERS_ORDERED:
        rsid_lower = rsid.lower()
        if rsid_lower in variants:
            chrom, pos, gt = variants[rsid_lower]
            # If the chromosome is Y or the genotype is not empty/null
            if chrom.upper() in ("Y", "24"):
                return True
            # If we find a derived allele on a Y marker, user is male
            if derived in gt.upper():
                return True
    # Also check if any variant is on chromosome Y
    for rsid, (chrom, pos, gt) in variants.items():
        if chrom.upper() in ("Y", "24"):
            return True
    return False


# ---------------------------------------------------------------------------
# HAPLOGROUP DETECTION (adapted for pipeline variants dict)
# ---------------------------------------------------------------------------

def _detect_y_haplogroup(
    variants: Dict[str, Tuple[str, str, str]],
    is_male: bool,
) -> Optional[str]:
    """Returns best Y haplogroup string, or None if female / not detected."""
    if not is_male:
        return None

    best: Optional[str] = None
    best_spec = -1
    for rsid, derived, clade, spec in Y_MARKERS_ORDERED:
        rsid_lower = rsid.lower()
        if rsid_lower not in variants:
            continue
        _chrom, _pos, gt = variants[rsid_lower]
        if derived in gt.upper():
            if spec > best_spec:
                best_spec = spec
                best = clade
    return best


def _detect_mt_haplogroup(
    variants: Dict[str, Tuple[str, str, str]],
) -> Optional[str]:
    """Returns best mtDNA haplogroup string."""
    signals: set = set()

    # Build position lookup for mitochondrial chromosome
    mt_positions: Dict[int, str] = {}
    for rsid, (chrom, pos, gt) in variants.items():
        if chrom.upper() in ("MT", "M", "26"):
            try:
                mt_positions[int(pos)] = gt
            except (ValueError, TypeError):
                continue

    # Position-based markers on chrM
    for pos, allele_map in MT_POS_MARKERS.items():
        gt = mt_positions.get(pos, "")
        for allele, marker in allele_map.items():
            if allele in gt.upper():
                signals.add(marker)

    # rsID-based backup markers
    for rsid, (allele, marker) in RSID_MT_MARKERS.items():
        rsid_lower = rsid.lower()
        if rsid_lower in variants:
            _chrom, _pos, gt = variants[rsid_lower]
            if allele in gt.upper():
                signals.add(marker)

    # Decision tree
    for required, hg, label in MT_DECISION_TREE:
        if required.issubset(signals):
            return hg

    return None if not signals else "R*"


# ---------------------------------------------------------------------------
# HAPLOGROUP MATCHING ENGINE
# ---------------------------------------------------------------------------

def _normalize_hg(hg: Optional[str]) -> str:
    """Standardize haplogroup name for prefix matching."""
    if hg is None:
        return ""
    return hg.upper().replace("_", "-")


def _match_figure(
    user_hg: Optional[str],
    match_list: List[Tuple[str, int, str]],
) -> Tuple[int, Optional[str]]:
    """
    Find best match for user_hg in match_list.
    match_list: [(prefix, score, label), ...]
    Returns: (score, label) or (0, None)
    """
    if not user_hg:
        return 0, None
    user_hg_norm = _normalize_hg(user_hg)
    best_score = 0
    best_label: Optional[str] = None
    for prefix, score, label in match_list:
        p = _normalize_hg(prefix)
        if user_hg_norm == p or user_hg_norm.startswith(p):
            if score > best_score:
                best_score = score
                best_label = label
    return best_score, best_label


def _score_to_label(score: int) -> str:
    labels = {
        5: "Exact Match",
        4: "Close Branch",
        3: "Secondary Clade",
        2: "Shared Macrohaplogroup",
        1: "Distant Lineage",
        0: "No Match",
    }
    return labels.get(score, "Unknown")


# ---------------------------------------------------------------------------
# OTZI COMPARISON
# ---------------------------------------------------------------------------

def _otzi_comparison(
    variants: Dict[str, Tuple[str, str, str]],
) -> Tuple[List[Dict[str, Any]], int, int, int]:
    """Compare user SNPs against Otzi's published genotypes."""
    results: List[Dict[str, Any]] = []
    matches = 0
    partial = 0
    total = 0

    for rsid, info in OTZI_GENOTYPES.items():
        otzi_gt = info["gt"]
        trait = info["trait"]
        odesc = info["otzi_desc"]
        ev = info["evidence"]
        disclaimer = OTZI_TRAIT_DISCLAIMERS.get(rsid, "")

        rsid_lower = rsid.lower()
        if rsid_lower in variants:
            _chrom, _pos, user_gt = variants[rsid_lower]
            user_gt = user_gt.upper().strip()
            total += 1

            if user_gt == otzi_gt or set(user_gt) == set(otzi_gt):
                status = "exact"
                status_label = "Identical to Otzi"
                matches += 1
            elif len(set(user_gt) & set(otzi_gt)) > 0:
                status = "partial"
                status_label = "Partial match (1 shared allele)"
                partial += 1
            else:
                status = "different"
                status_label = "Different from Otzi"

            results.append({
                "rsid": rsid,
                "trait": trait,
                "otziGenotype": otzi_gt,
                "userGenotype": user_gt,
                "otziDescription": odesc,
                "status": status,
                "statusLabel": status_label,
                "evidence": ev,
                "disclaimer": disclaimer,
            })
        else:
            results.append({
                "rsid": rsid,
                "trait": trait,
                "otziGenotype": otzi_gt,
                "userGenotype": None,
                "otziDescription": odesc,
                "status": "not_typed",
                "statusLabel": "Not typed on chip",
                "evidence": ev,
                "disclaimer": disclaimer,
            })

    return results, matches, partial, total


# ---------------------------------------------------------------------------
# CORE ANALYSIS
# ---------------------------------------------------------------------------

def analyze_historical_connections(
    variants: Dict[str, Tuple[str, str, str]],
) -> Dict[str, Any]:
    """
    Analyze genetic variants for historical connections via haplogroup matching
    and Otzi phenotypic comparison.

    Args:
        variants: Dictionary mapping rsID (lowercase) to (chromosome, position, genotype)

    Returns:
        Dictionary with analysis results
    """
    # Detect sex
    is_male = _detect_is_male(variants)

    # Detect haplogroups
    user_y_hg = _detect_y_haplogroup(variants, is_male)
    user_mt_hg = _detect_mt_haplogroup(variants)

    # Y-DNA connections
    y_connections: List[Dict[str, Any]] = []
    if is_male and user_y_hg:
        for fig in HISTORICAL_FIGURES:
            if not fig["y_matches"]:
                continue
            score, label = _match_figure(user_y_hg, fig["y_matches"])
            if score > 0:
                y_connections.append({
                    "figure": fig,
                    "matchScore": score,
                    "matchLabel": label,
                    "evidenceScore": fig["evidence_y"],
                })
        y_connections.sort(key=lambda x: -x["matchScore"])

    # mtDNA connections
    mt_connections: List[Dict[str, Any]] = []
    if user_mt_hg:
        for fig in HISTORICAL_FIGURES:
            if not fig["mt_matches"]:
                continue
            score, label = _match_figure(user_mt_hg, fig["mt_matches"])
            if score > 0:
                mt_connections.append({
                    "figure": fig,
                    "matchScore": score,
                    "matchLabel": label,
                    "evidenceScore": fig["evidence_mt"],
                })
        mt_connections.sort(key=lambda x: -x["matchScore"])

    # Otzi comparison
    otzi_results, exact_matches, partial_matches, total_compared = _otzi_comparison(variants)

    # Build top connections list
    top_connections: List[str] = []
    seen_names: set = set()
    for conn in sorted(y_connections + mt_connections, key=lambda x: -x["matchScore"]):
        name = conn["figure"]["name"]
        if name not in seen_names:
            seen_names.add(name)
            top_connections.append(name)
        if len(top_connections) >= 5:
            break

    return {
        "is_male": is_male,
        "user_y_hg": user_y_hg,
        "user_mt_hg": user_mt_hg,
        "y_connections": y_connections,
        "mt_connections": mt_connections,
        "otzi_results": otzi_results,
        "otzi_exact_matches": exact_matches,
        "otzi_partial_matches": partial_matches,
        "otzi_total_compared": total_compared,
        "top_connections": top_connections,
    }


# ---------------------------------------------------------------------------
# JSON OUTPUT
# ---------------------------------------------------------------------------

def _evidence_stars(score: int) -> str:
    """Return star string for evidence score."""
    stars = {5: "★★★★★", 4: "★★★★☆", 3: "★★★☆☆", 2: "★★☆☆☆", 1: "★☆☆☆☆", 0: "—"}
    return stars.get(score, "—")


def _build_connection_json(conn: Dict[str, Any], line_type: str, user_hg: Optional[str]) -> Dict[str, Any]:
    """Build a single historical connection JSON entry."""
    fig = conn["figure"]
    score = conn["matchScore"]
    label = conn["matchLabel"]
    ev_score = conn["evidenceScore"]

    return {
        "name": fig["name"],
        "dates": fig["dates"],
        "title": fig["title"],
        "category": fig.get("category", "historical"),
        "bio": fig["bio"],
        "origin": fig["origin"],
        "haplogroup": fig["y_hg"] if line_type == "Y" else fig["mt_hg"],
        "userHaplogroup": user_hg,
        "matchScore": score,
        "matchLabel": _score_to_label(score),
        "evidenceScore": ev_score,
        "evidenceStars": _evidence_stars(ev_score),
        "verified": fig["verified"],
        "source": fig["source"],
        "story": fig["story"],
    }


def generate_historical_connections_json(result: Dict[str, Any]) -> Dict[str, Any]:
    """
    Generate a JSON-serializable dict for the frontend historical connections report.

    Args:
        result: Dict from analyze_historical_connections()

    Returns:
        Dict matching the HistoricalConnectionsReport JSON schema
    """
    user_y_hg = result["user_y_hg"]
    user_mt_hg = result["user_mt_hg"]
    is_male = result["is_male"]

    # Build Y-DNA connections JSON
    y_dna_json: List[Dict[str, Any]] = []
    for conn in result["y_connections"]:
        y_dna_json.append(_build_connection_json(conn, "Y", user_y_hg))

    # Build mtDNA connections JSON
    mt_dna_json: List[Dict[str, Any]] = []
    for conn in result["mt_connections"]:
        mt_dna_json.append(_build_connection_json(conn, "mt", user_mt_hg))

    # Otzi comparison
    otzi_results = result["otzi_results"]
    exact = result["otzi_exact_matches"]
    partial = result["otzi_partial_matches"]
    total = result["otzi_total_compared"]

    # Build Otzi interpretation
    if total > 0:
        pct_exact = 100 * exact / total
        if pct_exact >= 80:
            interpretation = (
                "High degree of phenotypic similarity with Otzi - typical profile "
                "of well-preserved Neolithic Mediterranean (EEF) ancestry."
            )
        elif pct_exact >= 50:
            interpretation = (
                "Moderate similarity with Otzi - mix of ancestral Neolithic traits "
                "with more recent derived variants (post-Bronze Age)."
            )
        else:
            interpretation = (
                "Phenotypic profile distinct from Otzi - higher proportion of derived "
                "variants than in the Neolithic, which is expected in populations with "
                "significant Nordic/Steppe ancestry."
            )
    else:
        interpretation = "No Otzi traits could be compared - SNPs not available in genotyping data."

    otzi_comparison_json = {
        "traits": otzi_results,
        "exactMatches": exact,
        "partialMatches": partial,
        "totalCompared": total,
        "interpretation": interpretation,
        "disclaimer": (
            "Trait similarity does NOT equal direct kinship. Otzi lived 5,300 years ago. "
            "Real autosomal kinship has only been confirmed for G2a carriers with an "
            "identical haplotype (Niederstatter et al. 2016)."
        ),
    }

    # Disclaimers
    disclaimers = [
        (
            "Shared haplogroup does NOT mean direct descent from a historical individual. "
            "It means sharing a common patrilineal or matrilineal ancestor at some point "
            "in history - which could be 500 years ago or 50,000 years ago."
        ),
        (
            "Match levels indicate proximity on the human family tree: "
            "5 stars = same haplogroup branch; "
            "4 stars = common ancestor ~1,000-5,000 years; "
            "3 stars = common ancestor ~5,000-10,000 years; "
            "2 stars = common ancestor >10,000 years; "
            "1 star = very ancestral branch, no genealogical information."
        ),
        (
            "This analysis is educational. It does not constitute confirmed genealogy "
            "and does not replace high-resolution genealogical DNA testing."
        ),
        (
            "Evidence quality reflects the source of the historical figure's DNA: "
            "5 stars = confirmed from direct ancient DNA sequencing of verified remains; "
            "4 stars = multiple confirmations or from authenticated relatives; "
            "3 stars = population genetics study, not from individual remains; "
            "2 stars = reported but methodologically questioned; "
            "1 star = hypothesis only, no published data."
        ),
    ]

    return {
        "reportType": "historical_connections",
        "version": "2.0",
        "summary": {
            "userYHaplogroup": user_y_hg,
            "userMtHaplogroup": user_mt_hg,
            "isMale": is_male,
            "totalFigures": len(HISTORICAL_FIGURES),
            "yMatches": len(y_dna_json),
            "mtMatches": len(mt_dna_json),
            "topConnections": result["top_connections"],
        },
        "yDnaConnections": y_dna_json,
        "mtDnaConnections": mt_dna_json,
        "otziComparison": otzi_comparison_json,
        "disclaimers": disclaimers,
    }
