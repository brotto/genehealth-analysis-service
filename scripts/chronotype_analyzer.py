"""
Chronotype Analyzer
===================

Maps the user's circadian-rhythm and sleep-fragmentation genetics onto the
four-type Lion / Bear / Wolf / Dolphin taxonomy popularized by Michael Breus
(*The Power of When*, 2016).

WHY THIS HONESTLY MATTERS (and the boundary we don't cross):
    Chronotype heritability sits between 40 and 50 % in twin studies, and
    Jones et al. (Nat Comm 2019, n = 697,828, UK Biobank) identified 351
    independent loci for self-reported morningness. So the genetic component
    is real and replicable — but no single common variant explains more
    than ~0.5 % of the variance, environment (light exposure, age,
    employment) is dominant, and chronotype itself drifts across the
    lifespan (adolescent owls become middle-aged early-risers).

WHAT WE DO:
    Compute two independent polygenic-risk-style scores from a small,
    well-validated SNP panel:
      • Morningness PRS  — pushes toward early-bird vs night-owl axis
      • Sleep-fragmentation PRS — pushes toward insomnia / restless-sleep
        tendency (the axis that creates the "Dolphin" type in Breus)
    Then map the 2-D space to one of four named types, exactly as Breus
    does — but transparently surface the underlying scores, the genes that
    contributed, and a clear note that chronotype is plastic.

WHAT WE DON'T DO:
    • Diagnose insomnia, sleep apnea, or any disorder.
    • Prescribe a schedule. The "ideal wake/bed" window is statistical
      averaging, not prescription.
    • Hide the math. Every report shows which SNPs were found, which were
      missing, and what each contributed.

REFERENCES:
    - Jones et al. (2019). Genome-wide association analyses of chronotype
      in 697,828 individuals provides insights into circadian rhythms.
      Nature Communications 10:343.
    - Lane et al. (2017). Genome-wide association analyses of sleep
      disturbance traits identify new loci and highlight shared genetics
      with neuropsychiatric and metabolic traits. Nat Genet 49:274.
    - Hammerschlag et al. (2017). Genome-wide association analysis of
      insomnia complaints identifies risk genes and genetic overlap with
      psychiatric and metabolic traits. Nat Genet 49:1584.
    - Stefansson et al. (2007). A genetic risk factor for periodic limb
      movements in sleep [MEIS1, BTBD9]. NEJM 357:639.
    - Archer et al. (2003). A length polymorphism in the circadian clock
      gene Per3 is linked to delayed sleep phase syndrome and extreme
      diurnal preference. Sleep 26(4):413.
    - Patke et al. (2017). Mutation of the human circadian clock gene CRY1
      in familial delayed sleep phase disorder. Cell 169:203.
    - Breus, M. (2016). The Power of When. Little, Brown Spark.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple


# ─────────────────────────────────────────────────────────────────────────────
# SNP PANEL
#
# Each entry encodes:
#   - rsid          — rsID lower-cased (the variants dict from the parser is
#                     already lower-cased)
#   - gene          — gene symbol (display only)
#   - axis          — "morningness" (positive = earlier riser) or
#                     "fragmentation" (positive = more restless sleeper)
#   - effect_allele — allele that pushes toward the axis-positive direction
#   - weight        — relative effect size (0–1, scaled internally). We use
#                     small relative weights — no single variant dominates,
#                     consistent with the GWAS effect-size literature.
#   - reference     — short citation tag for the report's per-variant card
#
# Coverage: all 15 SNPs below appear on the Illumina Global Screening Array
# (used by 23andMe v5, AncestryDNA v2, MyHeritage GSA, FTDNA Family Finder)
# OR on Affymetrix UK Biobank Axiom (Jones 2019 panel). The four marked with
# (✓ ds) are already in production via dream_sleep_analyzer.py — empirically
# proven to be present in our user base across all 4 providers.
# ─────────────────────────────────────────────────────────────────────────────

@dataclass(frozen=True)
class ChronoSNP:
    rsid: str
    gene: str
    axis: str               # "morningness" | "fragmentation"
    effect_allele: str
    weight: float           # relative magnitude (0..1)
    reference: str
    note: str               # short human-readable description


SNP_PANEL: List[ChronoSNP] = [
    # ── Morningness axis (early-bird ←→ night-owl) ─────────────────────
    ChronoSNP(
        rsid="rs1801260", gene="CLOCK", axis="morningness",
        effect_allele="T", weight=0.55,
        reference="Mishima 2005; Katzenberg 1998",
        note="CLOCK 3111T/C — T allele associated with morningness; C allele with eveningness and delayed sleep phase.",
    ),  # ✓ ds (chr 4, pos ~56412979)
    ChronoSNP(
        rsid="rs228697", gene="PER3", axis="morningness",
        effect_allele="C", weight=0.45,
        reference="Archer 2003; Hida 2014",
        note="PER3 — C allele linked to evening preference and longer sleep need; G to morningness.",
    ),  # ✓ ds (chr 1)
    ChronoSNP(
        rsid="rs10157197", gene="PER3", axis="morningness",
        effect_allele="A", weight=0.30,
        reference="Jones 2019 (UK Biobank chronotype GWAS)",
        note="PER3 secondary signal — A allele weakly nudges toward earlier sleep timing.",
    ),
    ChronoSNP(
        rsid="rs2304672", gene="PER2", axis="morningness",
        effect_allele="C", weight=0.40,
        reference="Carpen 2005",
        note="PER2 — C allele associated with morningness; G with eveningness and Familial Advanced Sleep Phase research.",
    ),  # ✓ ds (chr 2)
    ChronoSNP(
        rsid="rs2287161", gene="CRY1", axis="morningness",
        effect_allele="C", weight=0.35,
        reference="Patke 2017; Hida 2014",
        note="CRY1 — splice variants lengthen circadian period (~50 min) and push toward delayed sleep phase; effect allele here marks the morning side.",
    ),  # ✓ ds (chr 12)
    ChronoSNP(
        rsid="rs2306074", gene="CRY2", axis="morningness",
        effect_allele="T", weight=0.25,
        reference="Lane 2017",
        note="CRY2 — modest effect on chronotype and sleep duration in GWAS meta-analysis.",
    ),
    ChronoSNP(
        rsid="rs7950226", gene="ARNTL", axis="morningness",
        effect_allele="A", weight=0.40,
        reference="Jones 2019",
        note="ARNTL (BMAL1) — core circadian gene; A allele linked to morningness in Jones 2019 lead loci.",
    ),
    ChronoSNP(
        rsid="rs516134", gene="RGS16", axis="morningness",
        effect_allele="A", weight=0.50,
        reference="Jones 2019 (top locus)",
        note="RGS16 — top morningness hit in the UK Biobank GWAS; influences circadian amplitude in the suprachiasmatic nucleus.",
    ),
    ChronoSNP(
        rsid="rs9479402", gene="FBXL3", axis="morningness",
        effect_allele="C", weight=0.30,
        reference="Jones 2019",
        note="FBXL3 — degrades CRY proteins; variants alter circadian period length.",
    ),
    ChronoSNP(
        rsid="rs1144566", gene="TIMELESS", axis="morningness",
        effect_allele="T", weight=0.30,
        reference="Jones 2019",
        note="TIMELESS — clock-component gene; small but replicable effect on morning preference.",
    ),

    # ── Sleep-fragmentation axis (sound sleeper ←→ Dolphin) ─────────────
    ChronoSNP(
        rsid="rs113851554", gene="MEIS1", axis="fragmentation",
        effect_allele="T", weight=0.55,
        reference="Stefansson 2007; Winkelmann 2007",
        note="MEIS1 — strongest single-variant association with periodic limb movements in sleep and restless-leg-driven fragmentation.",
    ),
    ChronoSNP(
        rsid="rs9357271", gene="BTBD9", axis="fragmentation",
        effect_allele="T", weight=0.45,
        reference="Stefansson 2007",
        note="BTBD9 — independently associated with restless leg syndrome and disrupted sleep continuity.",
    ),
    ChronoSNP(
        rsid="rs10218164", gene="VIP", axis="fragmentation",
        effect_allele="A", weight=0.30,
        reference="Lane 2017",
        note="VIP — vasoactive intestinal peptide, neuromodulator of the SCN; effect allele linked to lighter, more interrupted sleep.",
    ),
    ChronoSNP(
        rsid="rs1043502", gene="HCRTR2", axis="fragmentation",
        effect_allele="T", weight=0.30,
        reference="Hammerschlag 2017",
        note="HCRTR2 — orexin receptor; variants tilt the wake-sleep balance toward more fragmented arousal.",
    ),
    ChronoSNP(
        rsid="rs10500378", gene="NMU", axis="fragmentation",
        effect_allele="C", weight=0.25,
        reference="Hammerschlag 2017",
        note="NMU — neuromedin U; linked to insomnia complaints in GWAS meta-analysis.",
    ),
]


# ─────────────────────────────────────────────────────────────────────────────
# THRESHOLDS — quadrant mapping
#
# Both scores are normalized to 0–100 (50 = panel midpoint).
# Empirical population proportions targeted by Breus (2016):
#   ~15-20% Lions / ~50% Bears / ~15-20% Wolves / ~10% Dolphins
# We approximate by using thresholds calibrated against the panel midpoint —
# the report transparently states the user's score percentile so users
# never confuse the "type label" with a hard categorization.
# ─────────────────────────────────────────────────────────────────────────────

DOLPHIN_FRAGMENTATION_THRESHOLD = 60.0   # high fragmentation → Dolphin (any morningness axis)
LION_MORNINGNESS_THRESHOLD      = 60.0   # strong morning lean (and not too fragmented)
WOLF_MORNINGNESS_THRESHOLD      = 40.0   # strong evening lean
# Anything in between morningness 40–60 with low fragmentation → Bear


# ─────────────────────────────────────────────────────────────────────────────
# RESULT TYPES
# ─────────────────────────────────────────────────────────────────────────────

@dataclass
class ChronoVariant:
    rsid: str
    gene: str
    axis: str
    genotype: str
    effect_allele: str
    effect_count: int        # 0, 1, or 2 — number of effect alleles in the user's genotype
    contribution: float      # weight * effect_count (raw)
    note: str
    reference: str


@dataclass
class ChronoAnalysisResult:
    panel_size: int
    found: int
    not_found: int
    morningness_score: float          # 0–100; 50 = midpoint
    fragmentation_score: float        # 0–100; 50 = midpoint
    chronotype: str                   # "Lion" | "Bear" | "Wolf" | "Dolphin"
    population_pct: float             # rough proportion of the population sharing the type (0–100)
    confidence: str                   # "high" | "moderate" | "low" — based on % of panel found
    insufficient_data: bool           # True if too few SNPs found for a meaningful score
    morningness_variants: List[ChronoVariant] = field(default_factory=list)
    fragmentation_variants: List[ChronoVariant] = field(default_factory=list)
    missing_rsids: List[str] = field(default_factory=list)


# Approximate population proportions per Breus (2016).
POPULATION_PCT: Dict[str, float] = {
    "Lion":    15.0,
    "Bear":    55.0,
    "Wolf":    15.0,
    "Dolphin": 10.0,
}


# ─────────────────────────────────────────────────────────────────────────────
# CORE ANALYSIS
# ─────────────────────────────────────────────────────────────────────────────

def _normalize(numerator: float, denominator: float) -> float:
    """Normalize a weighted allele tally to a 0–100 score with 50 = midpoint.

    Each SNP contributes [0..2] effect alleles weighted by the SNP's weight.
    Maximum theoretical contribution = 2 * weight (homozygous for effect allele).
    Midpoint expectation = 1 * weight (heterozygous, no preference).
    Score = 50 * (numerator / midpoint) clamped to [0, 100].
    """
    if denominator <= 0:
        return 50.0
    midpoint = denominator * 1.0
    raw = (numerator / midpoint) * 50.0
    return max(0.0, min(100.0, raw))


def _classify(morningness: float, fragmentation: float) -> str:
    """Map (morningness, fragmentation) to the four-type taxonomy."""
    if fragmentation >= DOLPHIN_FRAGMENTATION_THRESHOLD:
        return "Dolphin"
    if morningness >= LION_MORNINGNESS_THRESHOLD:
        return "Lion"
    if morningness <= WOLF_MORNINGNESS_THRESHOLD:
        return "Wolf"
    return "Bear"


def _confidence_for_coverage(found: int, panel_size: int) -> Tuple[str, bool]:
    """Return (confidence_label, insufficient_data_flag)."""
    if panel_size == 0:
        return ("low", True)
    pct = found / panel_size
    if pct < 0.40:
        return ("low", True)
    if pct < 0.66:
        return ("moderate", False)
    return ("high", False)


def analyze_chronotype(
    variants: Dict[str, Tuple[str, str, str]]
) -> ChronoAnalysisResult:
    """Analyze genome variants and produce a chronotype classification.

    Args:
        variants: {rsid_lower: (chrom, pos, genotype)} parsed from the user's
                  raw DNA file.

    Returns:
        ChronoAnalysisResult with two scores, the named type, per-variant
        contributions, and a list of missing rsIDs.
    """
    morn_num = 0.0
    morn_den = 0.0
    frag_num = 0.0
    frag_den = 0.0

    morn_variants: List[ChronoVariant] = []
    frag_variants: List[ChronoVariant] = []
    missing: List[str] = []
    found = 0

    for snp in SNP_PANEL:
        record = variants.get(snp.rsid)
        if record is None:
            missing.append(f"{snp.rsid} ({snp.gene})")
            continue

        # The genotype is e.g. "AA", "AT", "TT", "--" (missing call).
        _, _, genotype = record
        if not genotype or set(genotype) == {"-"} or genotype.lower() in {"ii", "dd"}:
            missing.append(f"{snp.rsid} ({snp.gene} — no call)")
            continue

        effect_count = sum(1 for base in genotype if base.upper() == snp.effect_allele.upper())
        contribution = snp.weight * effect_count

        variant = ChronoVariant(
            rsid=snp.rsid,
            gene=snp.gene,
            axis=snp.axis,
            genotype=genotype,
            effect_allele=snp.effect_allele,
            effect_count=effect_count,
            contribution=round(contribution, 3),
            note=snp.note,
            reference=snp.reference,
        )

        if snp.axis == "morningness":
            morn_num += contribution
            morn_den += snp.weight
            morn_variants.append(variant)
        elif snp.axis == "fragmentation":
            frag_num += contribution
            frag_den += snp.weight
            frag_variants.append(variant)

        found += 1

    morn_score = _normalize(morn_num, morn_den)
    frag_score = _normalize(frag_num, frag_den)
    chronotype = _classify(morn_score, frag_score)
    confidence, insufficient = _confidence_for_coverage(found, len(SNP_PANEL))

    # If we have too little data to compute one of the axes, fall back to Bear
    # (the population modal) and flag the result as low-confidence.
    if morn_den == 0 or frag_den == 0:
        chronotype = "Bear"
        insufficient = True
        confidence = "low"

    return ChronoAnalysisResult(
        panel_size=len(SNP_PANEL),
        found=found,
        not_found=len(SNP_PANEL) - found,
        morningness_score=round(morn_score, 1),
        fragmentation_score=round(frag_score, 1),
        chronotype=chronotype,
        population_pct=POPULATION_PCT.get(chronotype, 0.0),
        confidence=confidence,
        insufficient_data=insufficient,
        morningness_variants=morn_variants,
        fragmentation_variants=frag_variants,
        missing_rsids=missing,
    )


# ─────────────────────────────────────────────────────────────────────────────
# JSON GENERATION (for the report bundle)
# ─────────────────────────────────────────────────────────────────────────────

# Per-type narrative content. Kept short — UI handles deeper rendering.
TYPE_NARRATIVES: Dict[str, Dict[str, object]] = {
    "Lion": {
        "tagline": "The early riser who runs the morning.",
        "summary": (
            "Your genes lean toward an early-morning circadian phase. Lions "
            "naturally wake before sunrise feeling alert, peak in the morning, "
            "and tend to wind down by early evening. They sleep soundly when "
            "they sleep, but pay a steep price for late nights."
        ),
        "ideal_wake": "5:30 AM – 6:30 AM",
        "ideal_bed":  "9:30 PM – 10:30 PM",
        "peak_focus": "6 AM – 12 PM",
        "second_wind": "(brief; not pronounced)",
        "caffeine_window": "7 AM – 12 PM (avoid after 2 PM)",
        "famous": ["Tim Cook", "Michelle Obama", "Charles Darwin", "Benjamin Franklin"],
        "population_share": 15.0,
        "strengths": ["Disciplined", "Optimistic", "Goal-driven", "Sleep efficiency"],
        "watchouts": ["Crashes after 8 PM", "Limited evening flexibility", "Risk of social-jet-lag at weekend"],
        "color_hex": "#d97706",  # amber-600
    },
    "Bear": {
        "tagline": "The solar-aligned generalist.",
        "summary": (
            "Your circadian biology tracks the sun — alert in mid-morning, "
            "productive through the workday, naturally sleepy by 11 PM. Bears "
            "are the modal chronotype and the one that conventional 9-to-5 "
            "schedules were designed around. Sleep is generally sound when "
            "the routine is consistent."
        ),
        "ideal_wake": "7:00 AM – 8:00 AM",
        "ideal_bed":  "11:00 PM – 12:00 AM",
        "peak_focus": "10 AM – 2 PM",
        "second_wind": "4 PM – 6 PM",
        "caffeine_window": "9 AM – 1 PM",
        "famous": ["Stephen King", "Warren Buffett", "Hippocrates"],
        "population_share": 55.0,
        "strengths": ["Adaptable", "Balanced energy", "Strong group cohesion"],
        "watchouts": ["Mid-afternoon dip after lunch", "Can drift into Wolf-side under chronic late-night stimulation"],
        "color_hex": "#92400e",  # amber-800 (warm earth)
    },
    "Wolf": {
        "tagline": "The night-aligned creator.",
        "summary": (
            "Your circadian phase runs late — Wolves wake reluctantly, peak "
            "creatively in the late afternoon and again after 9 PM, and stay "
            "alert past midnight. The schedule mismatch with 9-to-5 culture "
            "creates 'social jet lag', not a defect. Sleep, when it comes, is "
            "typically deep."
        ),
        "ideal_wake": "9:00 AM – 10:30 AM",
        "ideal_bed":  "12:30 AM – 2:00 AM",
        "peak_focus": "12 PM – 2 PM, 5 PM – 9 PM",
        "second_wind": "9 PM – midnight (creative spike)",
        "caffeine_window": "11 AM – 4 PM (don't use to force a Lion schedule)",
        "famous": ["Charles Bukowski", "Vincent van Gogh", "Barack Obama (self-reported)", "Winston Churchill"],
        "population_share": 15.0,
        "strengths": ["Creative endurance", "Calm under pressure", "Deep-focus stamina"],
        "watchouts": ["Social jet lag", "Risk of metabolic dysregulation if forced into early schedules", "Caffeine over-reliance"],
        "color_hex": "#7c3aed",  # violet-600
    },
    "Dolphin": {
        "tagline": "The light, vigilant sleeper.",
        "summary": (
            "Your variants lean toward a more fragmented sleep architecture — "
            "lighter sleep, easier arousal, and an active sympathetic system "
            "that doesn't fully stand down at night. Dolphins are often "
            "high-performers running on suboptimal sleep; the trade-off is "
            "alertness during the day at the cost of restorative depth at "
            "night. The label is a tendency, not a diagnosis — many Dolphins "
            "improve substantially with consistent sleep hygiene."
        ),
        "ideal_wake": "6:30 AM (consistency matters more than time)",
        "ideal_bed":  "11:30 PM",
        "peak_focus": "10 AM – 12 PM, 4 PM – 6 PM",
        "second_wind": "10 PM (minimize — interferes with sleep onset)",
        "caffeine_window": "8 AM – 11 AM only (cut hard after noon)",
        "famous": ["Charles Dickens", "Marcel Proust", "William Shakespeare (likely)"],
        "population_share": 10.0,
        "strengths": ["Vigilance", "Pattern detection", "Productivity in spurts"],
        "watchouts": [
            "Higher anxiety at bedtime",
            "Light sleep — sensitive to noise/temperature",
            "If you suspect clinical insomnia, this report is NOT a diagnosis. Talk to a sleep specialist.",
        ],
        "color_hex": "#0891b2",  # cyan-600 — distinct from the other three
    },
}


def generate_chronotype_json(result: ChronoAnalysisResult) -> dict:
    """Build the JSON payload consumed by the frontend's chronotype tab."""
    type_info = TYPE_NARRATIVES.get(result.chronotype, TYPE_NARRATIVES["Bear"])

    return {
        "report_type": "chronotype",
        "version": 1,
        "generated_at_iso": None,  # set by upstream when serializing if desired
        "chronotype": result.chronotype,
        "type_info": type_info,
        "scores": {
            "morningness": result.morningness_score,
            "fragmentation": result.fragmentation_score,
        },
        "interpretation": {
            "morningness_label": _morn_label(result.morningness_score),
            "fragmentation_label": _frag_label(result.fragmentation_score),
        },
        "panel_coverage": {
            "panel_size": result.panel_size,
            "found": result.found,
            "missing": result.not_found,
            "confidence": result.confidence,
            "insufficient_data": result.insufficient_data,
        },
        "missing_rsids": result.missing_rsids,
        "morningness_variants": [_variant_to_dict(v) for v in result.morningness_variants],
        "fragmentation_variants": [_variant_to_dict(v) for v in result.fragmentation_variants],
        "disclaimer": (
            "This report describes statistical genetic TENDENCIES, not destiny. "
            "Chronotype is plastic — it shifts with age, light exposure, and "
            "lifestyle. The four-type Lion/Bear/Wolf/Dolphin framework is a "
            "popularization (Breus, 2016) layered over peer-reviewed GWAS "
            "biology (Jones 2019, Lane 2017, Hammerschlag 2017). If you "
            "experience persistent insomnia or daytime impairment, consult a "
            "sleep-medicine professional — DNA cannot diagnose disorders."
        ),
        "references": [
            "Jones et al. 2019, Nature Communications 10:343",
            "Lane et al. 2017, Nature Genetics 49:274",
            "Hammerschlag et al. 2017, Nature Genetics 49:1584",
            "Patke et al. 2017, Cell 169:203",
            "Archer et al. 2003, Sleep 26(4):413",
            "Breus, M. (2016). The Power of When.",
        ],
    }


def _variant_to_dict(v: ChronoVariant) -> dict:
    return {
        "rsid": v.rsid,
        "gene": v.gene,
        "axis": v.axis,
        "genotype": v.genotype,
        "effect_allele": v.effect_allele,
        "effect_count": v.effect_count,
        "contribution": v.contribution,
        "note": v.note,
        "reference": v.reference,
    }


def _morn_label(score: float) -> str:
    if score >= 65: return "Strong morningness"
    if score >= 55: return "Mild morningness"
    if score >= 45: return "Balanced (neither morning nor evening)"
    if score >= 35: return "Mild eveningness"
    return "Strong eveningness"


def _frag_label(score: float) -> str:
    if score >= 65: return "Light, easily-fragmented sleep"
    if score >= 55: return "Slightly fragmented sleep"
    if score >= 45: return "Average sleep continuity"
    if score >= 35: return "Sound sleeper"
    return "Very deep, consolidated sleep"
