"""
Fun Traits Analyzer
Analyzes lighthearted, shareable genetic traits across five categories:
Taste & Smell, Physical Quirks, Sensory Perception, Sleep & Daily Rhythm,
and Miscellaneous Fun.

Unlike scored analyzers, fun traits are descriptive — each variant maps
to a specific outcome (e.g., "cilantro tastes like soap") rather than a
score on a dimension.

References:
- Kim et al. (2003) Science - TAS2R38 and PTC bitter taste
- Eriksson et al. (2012) Flavour - OR6A2 and cilantro
- Pelchat et al. (2011) Chem Senses - Asparagus urine smell
- Lindpaintner (2008) - ABCC11 earwax and body odor
- Kimura et al. (2009) Am J Hum Genet - EDAR and East Asian hair/teeth
- Sulem et al. (2007) Nat Genet - HERC2/OCA2 eye color
- Valverde et al. (1995) Nat Genet - MC1R red hair/freckles
- Han et al. (2020) TRPV1 spicy tolerance
- Higuchi et al. (2004) - ALDH2 alcohol flush
- Bhatti et al. (2013) - ADH1B and alcohol metabolism
- Cornelis et al. (2011) - CYP1A2 caffeine metabolism
- Archer et al. (2003) Sleep - PER3 VNTR and chronotype
- Katzenberg et al. (1998) Sleep - CLOCK and morningness
- He et al. (2009) Science - DEC2 and short sleeper phenotype
- Bachner-Melman et al. (2005) - AVPR1A and dance creativity
- Ebstein et al. (1996) - DRD4 and novelty seeking
- Jonsson et al. (2010) Nature - IRF4 and hair graying
- Eriksson et al. (2010) PLoS Genet - Photic sneeze reflex
"""

from typing import Dict, List, Tuple, Any
from dataclasses import dataclass, field


# ─────────────────────────────────────────────────────────────────────────────
# COMPLEMENT MAP (for strand-flip handling)
# ─────────────────────────────────────────────────────────────────────────────

COMPLEMENT = {"A": "T", "T": "A", "C": "G", "G": "C"}


def complement_allele(allele: str) -> str:
    """Return the complement of a single nucleotide."""
    return COMPLEMENT.get(allele.upper(), allele.upper())


def count_allele(genotype: str, allele: str) -> int:
    """
    Count occurrences of an allele in a genotype, handling complement strands.
    """
    genotype = genotype.upper().replace("-", "")
    allele = allele.upper()
    comp = complement_allele(allele)

    count = genotype.count(allele)

    if count == 0 and comp != allele:
        count = genotype.count(comp)

    return min(count, 2)


# ─────────────────────────────────────────────────────────────────────────────
# CATEGORY DEFINITIONS
# ─────────────────────────────────────────────────────────────────────────────

CATEGORIES = [
    "Taste & Smell",
    "Physical Quirks",
    "Sensory Perception",
    "Sleep & Daily Rhythm",
    "Miscellaneous Fun",
]

CATEGORY_EMOJI = {
    "Taste & Smell": "\U0001f445",          # 👅
    "Physical Quirks": "\u2728",            # ✨
    "Sensory Perception": "\U0001f441\ufe0f", # 👁️
    "Sleep & Daily Rhythm": "\U0001f634",   # 😴
    "Miscellaneous Fun": "\U0001f3b2",      # 🎲
}


# ─────────────────────────────────────────────────────────────────────────────
# FUN TRAIT VARIANT DATABASE
# ─────────────────────────────────────────────────────────────────────────────

@dataclass
class FunTraitSNP:
    rsid: str
    gene: str
    category: str
    trait: str                     # short trait name
    # map of allele-count (0, 1, 2 of target_allele) -> outcome text
    target_allele: str
    result_by_count: Dict[int, str]
    description: str               # what this SNP does in general
    fun_fact: str                  # interesting science tidbit


FUN_VARIANTS: Dict[str, FunTraitSNP] = {
    # ── TASTE & SMELL ──
    "rs713598": FunTraitSNP(
        rsid="rs713598", gene="TAS2R38", category="Taste & Smell",
        trait="Bitter taste perception (PTC)",
        target_allele="G",
        result_by_count={
            2: "You are likely a SUPER-TASTER — bitter foods like broccoli, "
               "grapefruit and dark chocolate taste very intense to you.",
            1: "You are likely a medium taster of bitter compounds — you "
               "detect bitterness but it is not overwhelming.",
            0: "You are likely a NON-TASTER — bitter compounds barely "
               "register for you, so brussels sprouts may taste almost sweet.",
        },
        description="TAS2R38 encodes a bitter taste receptor. This variant "
                    "(P49A) determines whether PTC and related compounds taste bitter.",
        fun_fact="About 25% of people are non-tasters, 50% medium tasters, "
                 "and 25% super-tasters. Supertasters also tend to dislike "
                 "coffee, beer and spinach.",
    ),
    "rs1726866": FunTraitSNP(
        rsid="rs1726866", gene="TAS2R38", category="Taste & Smell",
        trait="Bitter taste (A262V variant)",
        target_allele="A",
        result_by_count={
            2: "Reinforces a sensitive bitter-taste phenotype — strong "
               "reaction to brassica vegetables.",
            1: "Intermediate bitter sensitivity.",
            0: "Weaker bitter sensitivity; less aversion to greens.",
        },
        description="Second major TAS2R38 SNP contributing to PTC haplotype.",
        fun_fact="Combined with rs713598, the PAV/AVI haplotype explains "
                 "most of the variance in PTC bitter tasting.",
    ),
    "rs72921001": FunTraitSNP(
        rsid="rs72921001", gene="OR6A2", category="Taste & Smell",
        trait="Cilantro soap taste",
        target_allele="C",
        result_by_count={
            2: "You probably think cilantro tastes like SOAP or metal — "
               "your olfactory receptor is tuned to its aldehydes.",
            1: "You may occasionally notice a soapy note in cilantro, but "
               "likely tolerate it fine.",
            0: "You probably enjoy cilantro without any soapy taste.",
        },
        description="OR6A2 is an olfactory receptor that detects the "
                    "aldehydes in fresh cilantro (coriander).",
        fun_fact="The 'cilantro tastes like soap' phenomenon affects roughly "
                 "4 to 14% of the population depending on ancestry.",
    ),
    "rs4481887": FunTraitSNP(
        rsid="rs4481887", gene="OR2M7", category="Taste & Smell",
        trait="Asparagus urine smell detection",
        target_allele="A",
        result_by_count={
            2: "You can smell 'asparagus pee' — your nose detects the sulfur "
               "compounds produced after eating asparagus.",
            1: "You can likely smell asparagus urine, though possibly less "
               "strongly.",
            0: "You probably cannot smell asparagus pee at all — you are an "
               "'asparagus anosmic'.",
        },
        description="OR2M7 is an olfactory receptor gene involved in "
                    "detecting asparagus urinary metabolites.",
        fun_fact="Everyone who eats asparagus produces the smelly "
                 "compounds, but only ~40% of people can actually smell them.",
    ),
    "rs7501331": FunTraitSNP(
        rsid="rs7501331", gene="BCMO1", category="Taste & Smell",
        trait="Beta-carotene conversion efficiency",
        target_allele="T",
        result_by_count={
            2: "You convert beta-carotene from vegetables into vitamin A "
               "less efficiently — you may need more dietary vitamin A.",
            1: "You have reduced but still functional beta-carotene conversion.",
            0: "You convert beta-carotene efficiently into vitamin A.",
        },
        description="BCMO1 encodes the enzyme that converts plant "
                    "beta-carotene into active vitamin A (retinol).",
        fun_fact="Low converters may benefit more from animal-source vitamin "
                 "A (liver, eggs, dairy) than from carrots.",
    ),

    # ── PHYSICAL QUIRKS ──
    "rs3827760": FunTraitSNP(
        rsid="rs3827760", gene="EDAR", category="Physical Quirks",
        trait="Thick hair & shovel-shaped incisors",
        target_allele="C",
        result_by_count={
            2: "You likely have THICKER hair shafts and possibly "
               "shovel-shaped front teeth — the classic EDAR 'V370A' variant.",
            1: "You may have moderately thick hair and mild shoveling of incisors.",
            0: "You have typical-diameter hair and flat incisors.",
        },
        description="EDAR V370A is an East Asian-associated variant "
                    "affecting hair, sweat glands and tooth shape.",
        fun_fact="This variant arose ~30,000 years ago in East Asia and is "
                 "found in over 90% of modern East Asians.",
    ),
    "rs1805007": FunTraitSNP(
        rsid="rs1805007", gene="MC1R", category="Physical Quirks",
        trait="Red hair & freckles (R151C)",
        target_allele="T",
        result_by_count={
            2: "Very high chance of red hair, fair skin, and freckles — "
               "you carry two copies of this classic MC1R red variant.",
            1: "Carrier — you might have red undertones, freckles, or "
               "redheaded children.",
            0: "No copies of this red-hair variant, though others exist.",
        },
        description="MC1R R151C is one of the three main red-hair variants "
                    "(R, R151C; R160W; D294H).",
        fun_fact="Redheads need ~20% more anesthetic than average because "
                 "MC1R also affects pain perception.",
    ),
    "rs1805008": FunTraitSNP(
        rsid="rs1805008", gene="MC1R", category="Physical Quirks",
        trait="Red hair & UV sensitivity (R160W)",
        target_allele="T",
        result_by_count={
            2: "Strong red hair and very fair, UV-sensitive skin likely.",
            1: "Carrier of the R160W variant — possible freckles and "
               "sun sensitivity.",
            0: "No R160W copies present.",
        },
        description="MC1R R160W is another major red-hair-associated variant.",
        fun_fact="Multiple 'red' MC1R variants exist; the more red alleles "
                 "you carry, the more red-haired the phenotype.",
    ),
    "rs12913832": FunTraitSNP(
        rsid="rs12913832", gene="HERC2", category="Physical Quirks",
        trait="Blue vs brown eyes",
        target_allele="G",
        result_by_count={
            2: "Very likely BLUE or green eyes — you carry two copies of "
               "the main eye-color variant.",
            1: "Likely mixed eye color (hazel, green, or intermediate).",
            0: "Most likely BROWN eyes.",
        },
        description="HERC2/OCA2 rs12913832 is the single biggest genetic "
                    "predictor of blue vs brown eyes in Europeans.",
        fun_fact="All blue-eyed people likely share a common ancestor from "
                 "6,000 to 10,000 years ago who first carried this mutation.",
    ),
    "rs8065080": FunTraitSNP(
        rsid="rs8065080", gene="TRPV1", category="Physical Quirks",
        trait="Spicy food tolerance",
        target_allele="T",
        result_by_count={
            2: "You probably tolerate spicy food WELL — your capsaicin "
               "receptor is less sensitive.",
            1: "Moderate spicy tolerance.",
            0: "You probably find spicy food more intense — your TRPV1 "
               "receptor is fully sensitive.",
        },
        description="TRPV1 is the capsaicin (chili pepper) receptor. "
                    "Variants alter how strongly spicy food registers as pain.",
        fun_fact="Capsaicin is not a taste — it binds the TRPV1 heat/pain "
                 "receptor. Your brain interprets it as burning.",
    ),
    "rs17822931": FunTraitSNP(
        rsid="rs17822931", gene="ABCC11", category="Physical Quirks",
        trait="Earwax type & body odor",
        target_allele="T",
        result_by_count={
            2: "You likely have DRY, flaky earwax and reduced underarm "
               "body odor. Very common in East Asian populations.",
            1: "Mixed — possibly intermediate earwax and mild body odor.",
            0: "You likely have WET (sticky) earwax and typical body odor.",
        },
        description="ABCC11 G180R determines whether earwax is wet or dry "
                    "and influences apocrine sweat (body odor).",
        fun_fact="This is one of the few human genes where a single SNP "
                 "perfectly determines a visible trait — wet vs dry earwax.",
    ),

    # ── SENSORY PERCEPTION ──
    "rs671": FunTraitSNP(
        rsid="rs671", gene="ALDH2", category="Sensory Perception",
        trait="Alcohol flush reaction",
        target_allele="A",
        result_by_count={
            2: "Severe alcohol flush likely — you lack functional ALDH2 "
               "and accumulate toxic acetaldehyde. Drinking carries higher "
               "cancer risk for you.",
            1: "You probably experience the 'Asian flush' — face reddening, "
               "fast heartbeat, nausea after alcohol.",
            0: "Typical alcohol metabolism without flushing from this variant.",
        },
        description="ALDH2*2 disables aldehyde dehydrogenase 2, preventing "
                    "breakdown of acetaldehyde.",
        fun_fact="~540 million people, mostly East Asian, carry this variant. "
                 "The WHO classifies acetaldehyde as a Group 1 carcinogen.",
    ),
    "rs1229984": FunTraitSNP(
        rsid="rs1229984", gene="ADH1B", category="Sensory Perception",
        trait="Alcohol metabolism speed",
        target_allele="A",
        result_by_count={
            2: "You metabolize alcohol VERY FAST into acetaldehyde — often "
               "causing unpleasant symptoms after drinking.",
            1: "Faster-than-average alcohol metabolism.",
            0: "Typical alcohol metabolism rate.",
        },
        description="ADH1B*2 (Arg48His) is a super-active alcohol "
                    "dehydrogenase variant.",
        fun_fact="ADH1B*2 carriers have ~60% lower alcoholism risk because "
                 "drinking quickly becomes unpleasant.",
    ),
    "rs762551": FunTraitSNP(
        rsid="rs762551", gene="CYP1A2", category="Sensory Perception",
        trait="Caffeine metabolism",
        target_allele="A",
        result_by_count={
            2: "You are a FAST caffeine metabolizer — coffee has minimal "
               "lingering effects and you can drink it late without trouble.",
            1: "Intermediate caffeine metabolism.",
            0: "You are a SLOW caffeine metabolizer — caffeine stays in "
               "your system longer and late coffee may wreck your sleep.",
        },
        description="CYP1A2 -163C>A (1F variant) alters the activity of "
                    "the enzyme that clears ~95% of caffeine.",
        fun_fact="Slow metabolizers have ~2x higher risk of heart attack "
                 "from heavy coffee; fast metabolizers are protected.",
    ),
    "rs4481887_sensory_placeholder": FunTraitSNP(
        rsid="rs11642015", gene="FTO", category="Sensory Perception",
        trait="Perception of food reward",
        target_allele="T",
        result_by_count={
            2: "You may experience stronger food reward signals — higher "
               "hedonic response to palatable foods.",
            1: "Intermediate food reward perception.",
            0: "Typical food reward response.",
        },
        description="FTO variants modulate reward perception related to "
                    "food cues in the brain.",
        fun_fact="FTO was the first 'obesity gene' discovered via GWAS; "
                 "it influences how rewarding fatty foods feel.",
    ),
    "rs11562975": FunTraitSNP(
        rsid="rs11562975", gene="TRPM8", category="Sensory Perception",
        trait="Cold & menthol sensitivity",
        target_allele="T",
        result_by_count={
            2: "You may be especially sensitive to cold temperatures and "
               "menthol sensations.",
            1: "Intermediate cold sensitivity.",
            0: "Typical cold perception.",
        },
        description="TRPM8 is the cold/menthol receptor; variants alter "
                    "thermal perception.",
        fun_fact="TRPM8 is why menthol feels cool — it activates the same "
                 "receptor as actual cold temperatures.",
    ),

    # ── SLEEP & DAILY RHYTHM ──
    "rs228697": FunTraitSNP(
        rsid="rs228697", gene="PER3", category="Sleep & Daily Rhythm",
        trait="Morningness vs eveningness",
        target_allele="C",
        result_by_count={
            2: "You probably lean toward being a NIGHT OWL — natural late "
               "bedtime and late wake time.",
            1: "Intermediate chronotype — flexible sleep preference.",
            0: "You likely lean toward being a MORNING LARK — early bed, "
               "early rise.",
        },
        description="PER3 is a core circadian clock gene; variants shift "
                    "preferred sleep timing.",
        fun_fact="About 40% of people are clearly morning or evening types; "
                 "the rest are somewhere in between.",
    ),
    "rs10462021": FunTraitSNP(
        rsid="rs10462021", gene="PER3", category="Sleep & Daily Rhythm",
        trait="PER3 sleep homeostasis",
        target_allele="A",
        result_by_count={
            2: "You may be more vulnerable to sleep deprivation — needing "
               "full sleep to perform well.",
            1: "Intermediate sleep resilience.",
            0: "You may tolerate sleep loss relatively well.",
        },
        description="PER3 length polymorphism affects sleep homeostasis "
                    "and cognitive response to sleep loss.",
        fun_fact="People with long PER3 variants show bigger cognitive "
                 "deficits when sleep-deprived than short variants.",
    ),
    "rs1801260": FunTraitSNP(
        rsid="rs1801260", gene="CLOCK", category="Sleep & Daily Rhythm",
        trait="CLOCK 3111T>C (chronotype)",
        target_allele="C",
        result_by_count={
            2: "More likely evening preference, later bedtime, and delayed "
               "sleep phase tendency.",
            1: "Slight evening preference.",
            0: "More likely morning preference.",
        },
        description="CLOCK is the master circadian transcription factor; "
                    "3111T>C shifts diurnal preference.",
        fun_fact="The CLOCK gene was discovered in mice in 1994 and is "
                 "conserved across almost every living organism.",
    ),
    "rs11932595": FunTraitSNP(
        rsid="rs11932595", gene="ARNTL", category="Sleep & Daily Rhythm",
        trait="BMAL1 circadian amplitude",
        target_allele="G",
        result_by_count={
            2: "You may have a strong, stable circadian rhythm.",
            1: "Intermediate circadian amplitude.",
            0: "You may have a flatter circadian rhythm, possibly leading "
               "to irregular sleep patterns.",
        },
        description="ARNTL (BMAL1) partners with CLOCK to drive the "
                    "master circadian oscillator.",
        fun_fact="BMAL1 knockout mice lose all circadian rhythms and age "
                 "prematurely — it is essential for biological timing.",
    ),
    "rs11046205": FunTraitSNP(
        rsid="rs11046205", gene="ABCC9", category="Sleep & Daily Rhythm",
        trait="Habitual sleep duration",
        target_allele="A",
        result_by_count={
            2: "You may naturally sleep LESS than average (shorter "
               "habitual sleep).",
            1: "Slightly shorter sleep duration than average.",
            0: "You likely need average-to-longer sleep duration.",
        },
        description="ABCC9 was identified in a GWAS for habitual sleep "
                    "duration in the general population.",
        fun_fact="This variant was linked to ~30 minutes of difference in "
                 "nightly sleep duration across populations.",
    ),

    # ── MISCELLANEOUS FUN ──
    "rs3021529": FunTraitSNP(
        rsid="rs3021529", gene="AVPR1A", category="Miscellaneous Fun",
        trait="Dance & creative movement",
        target_allele="T",
        result_by_count={
            2: "You may have stronger natural inclination toward dance "
               "and expressive movement — AVPR1A is linked with dancing ability.",
            1: "Intermediate predisposition for dance/music response.",
            0: "Typical response — though anyone can learn to dance!",
        },
        description="AVPR1A (vasopressin receptor 1A) has been associated "
                    "with creative dance and musical behavior.",
        fun_fact="A 2005 study found professional dancers carry specific "
                 "AVPR1A variants more often than non-dancers.",
    ),
    "rs1800955": FunTraitSNP(
        rsid="rs1800955", gene="DRD4", category="Miscellaneous Fun",
        trait="Novelty seeking",
        target_allele="T",
        result_by_count={
            2: "You may score higher on NOVELTY SEEKING — craving new "
               "experiences and excitement.",
            1: "Moderate novelty seeking.",
            0: "You may prefer routine and familiarity.",
        },
        description="DRD4 -521C>T is a dopamine D4 receptor promoter "
                    "variant linked with novelty seeking behavior.",
        fun_fact="DRD4 variants are sometimes called the 'adventure gene' "
                 "and are more common in nomadic populations.",
    ),
    "rs12203592": FunTraitSNP(
        rsid="rs12203592", gene="IRF4", category="Miscellaneous Fun",
        trait="Hair graying",
        target_allele="T",
        result_by_count={
            2: "You may GRAY EARLIER than average — IRF4 is the first "
               "gene clearly linked with premature graying.",
            1: "Slightly earlier graying than average.",
            0: "Typical graying timeline.",
        },
        description="IRF4 rs12203592 was identified in the first GWAS "
                    "for hair graying, conducted in Latin Americans.",
        fun_fact="This was the first-ever gene identified for going gray. "
                 "It also affects freckles, eye color and hair color.",
    ),
    "rs10427255": FunTraitSNP(
        rsid="rs10427255", gene="ZEB2", category="Miscellaneous Fun",
        trait="Photic sneeze reflex (ACHOO)",
        target_allele="C",
        result_by_count={
            2: "You probably SNEEZE in bright sunlight — the photic "
               "sneeze reflex (ACHOO syndrome).",
            1: "You may sneeze in bright light occasionally.",
            0: "Typical — no sunlight sneezing.",
        },
        description="The ACHOO reflex (Autosomal Dominant Compelling "
                    "Helio-Ophthalmic Outburst) is linked to variants near ZEB2.",
        fun_fact="Around 18-35% of people sneeze when stepping into bright "
                 "sunlight. Aristotle was the first to describe it.",
    ),
    "rs1815739_fun": FunTraitSNP(
        rsid="rs2229616", gene="MC4R", category="Miscellaneous Fun",
        trait="Appetite & snacking tendency",
        target_allele="T",
        result_by_count={
            2: "You may have slightly reduced appetite signals — MC4R "
               "variants protect against obesity.",
            1: "Intermediate appetite regulation.",
            0: "Typical appetite regulation.",
        },
        description="MC4R V103I is a well-studied variant of the "
                    "melanocortin-4 receptor affecting hunger signaling.",
        fun_fact="MC4R is the most common monogenic cause of obesity, but "
                 "the V103I variant is actually protective.",
    ),
}


# ─────────────────────────────────────────────────────────────────────────────
# ANALYSIS RESULT DATACLASSES
# ─────────────────────────────────────────────────────────────────────────────

@dataclass
class FunVariantResult:
    rsid: str
    gene: str
    category: str
    trait: str
    genotype: str
    result: str          # the outcome text for this person
    description: str
    fun_fact: str


@dataclass
class FunTraitsAnalysisResult:
    total_checked: int
    found: int
    not_found: int
    results_by_category: Dict[str, List[FunVariantResult]]
    missing_rsids: List[str]


# ─────────────────────────────────────────────────────────────────────────────
# CORE ANALYSIS
# ─────────────────────────────────────────────────────────────────────────────

def analyze_fun_traits(
    variants: Dict[str, Tuple[str, str, str]]
) -> FunTraitsAnalysisResult:
    """
    Analyze genetic variants for fun, shareable traits.

    Args:
        variants: Dictionary mapping rsID (lowercase) to (chromosome, position, genotype)

    Returns:
        FunTraitsAnalysisResult with per-category descriptive results
    """
    results_by_category: Dict[str, List[FunVariantResult]] = {
        cat: [] for cat in CATEGORIES
    }
    missing_rsids: List[str] = []

    found = 0
    not_found = 0

    for _key, snp_info in FUN_VARIANTS.items():
        rsid = snp_info.rsid
        rsid_lower = rsid.lower()

        if rsid_lower in variants:
            found += 1
            _chrom, _pos, genotype = variants[rsid_lower]
            genotype = genotype.upper().replace("-", "")

            target_count = count_allele(genotype, snp_info.target_allele)
            result_text = snp_info.result_by_count.get(
                target_count,
                "Indeterminate — genotype not recognized.",
            )

            results_by_category[snp_info.category].append(
                FunVariantResult(
                    rsid=rsid,
                    gene=snp_info.gene,
                    category=snp_info.category,
                    trait=snp_info.trait,
                    genotype=genotype,
                    result=result_text,
                    description=snp_info.description,
                    fun_fact=snp_info.fun_fact,
                )
            )
        else:
            not_found += 1
            missing_rsids.append(f"{rsid} ({snp_info.gene})")

    return FunTraitsAnalysisResult(
        total_checked=len(FUN_VARIANTS),
        found=found,
        not_found=not_found,
        results_by_category=results_by_category,
        missing_rsids=missing_rsids,
    )


# ─────────────────────────────────────────────────────────────────────────────
# JSON OUTPUT
# ─────────────────────────────────────────────────────────────────────────────

def generate_fun_traits_json(result: FunTraitsAnalysisResult) -> dict:
    """
    Generate a JSON-serializable dict for the frontend fun traits report.

    Args:
        result: FunTraitsAnalysisResult from analyze_fun_traits()

    Returns:
        Dict matching the fun traits report JSON schema
    """
    categories_list = []

    for cat in CATEGORIES:
        findings = result.results_by_category.get(cat, [])
        variants_json = []

        for f in findings:
            variants_json.append({
                "rsid": f.rsid,
                "gene": f.gene,
                "trait": f.trait,
                "genotype": f.genotype,
                "result": f.result,
                "description": f.description,
                "funFact": f.fun_fact,
            })

        categories_list.append({
            "name": cat,
            "emoji": CATEGORY_EMOJI.get(cat, ""),
            "variants": variants_json,
        })

    return {
        "summary": {
            "totalChecked": result.total_checked,
            "found": result.found,
            "notFound": result.not_found,
        },
        "categories": categories_list,
    }


# ─────────────────────────────────────────────────────────────────────────────
# MARKDOWN REPORT (optional, for compatibility)
# ─────────────────────────────────────────────────────────────────────────────

def generate_fun_traits_report(
    result: FunTraitsAnalysisResult, subject_name: str = "Subject"
) -> str:
    """Generate a Markdown fun traits report."""
    from datetime import datetime

    report = f"""# Fun Genetic Traits Report

**Subject:** {subject_name}
**Generated:** {datetime.now().strftime("%B %d, %Y")}

---

## Overview

These are lighthearted, shareable genetic traits — the quirky stuff people
love to discover about themselves. They are not health information.

| Metric | Value |
|--------|-------|
| Total Fun SNPs Checked | {result.total_checked} |
| SNPs Found in Your Data | {result.found} |
| SNPs Not Available | {result.not_found} |

---

"""

    for cat in CATEGORIES:
        findings = result.results_by_category.get(cat, [])
        if not findings:
            continue

        emoji = CATEGORY_EMOJI.get(cat, "")
        report += f"## {emoji} {cat}\n\n"

        for f in findings:
            report += f"### {f.gene} ({f.rsid}) — {f.trait}\n\n"
            report += f"- **Your genotype:** {f.genotype}\n"
            report += f"- **Result:** {f.result}\n"
            report += f"- **About this gene:** {f.description}\n"
            report += f"- **Fun fact:** {f.fun_fact}\n\n"

        report += "---\n\n"

    report += """## Important Notes

1. **For fun only** — None of these traits are medical predictions. They are
   lighthearted associations from peer-reviewed research.

2. **Probabilistic, not deterministic** — Even the strongest association
   genes only increase the likelihood of a trait. Environment and chance
   still play major roles.

3. **Share away** — These are the kinds of results people love to compare
   with friends and family. Have fun!

---
*Generated by GeneHealth Analysis Platform*
"""

    return report
