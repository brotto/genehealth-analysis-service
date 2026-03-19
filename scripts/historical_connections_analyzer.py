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
        "version": "1.0",
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
