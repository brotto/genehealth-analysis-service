"""
Career Genes Analyzer
=====================

Maps personality and cognitive genetics to Holland's RIASEC career-style
taxonomy. Outputs ranked work-style tendencies with example career fits.

CRITICAL DISCLAIMER (also embedded in JSON output):
Genes do NOT determine your career. This report describes statistical
TENDENCIES toward broad work-style categories observed in research
populations. Education, opportunity, training, interest, family context,
luck, and personal choice are decisive — not the variants below.

We use a small, well-studied SNP panel (the same loci used in the
Big Five Personality and Intelligence reports), then map underlying
trait scores to RIASEC dimensions using Holland's framework.

Why RIASEC and not "your gene says you should be a doctor"?
- RIASEC categories (Realistic / Investigative / Artistic / Social /
  Enterprising / Conventional) are validated work-style dimensions, not
  job titles. Two careers can sit in the same dimension.
- Vocational interest twin studies (Schermer & Vernon 2008; Lykken &
  Bouchard) report heritability in the 0.30–0.50 range, BUT the variance
  explained by any single common variant is well under 0.5%. We frame
  outputs as gentle tendencies, never as prescriptions.
- Educational-attainment GWAS (Lee et al. 2018, n = 1.1M) supplies
  cognitive-component loci; risk-taking GWAS (Sanchez-Roige et al. 2018)
  supplies the enterprising / artistic component.

References (further reading inside the report):
- Holland, J. L. (1997). Making Vocational Choices, 3rd ed.
- Plomin & Deary (2015). Genetics and intelligence differences. Mol Psych.
- Lee et al. (2018). Gene discovery for educational attainment. Nat Genet.
- Sanchez-Roige et al. (2018). Genome-wide association of risk-taking. Nat Hum Behav.
- Schermer & Vernon (2008). The heritability of vocational interests. Behav Genet.
- Lykken et al. (1993). Heritability of interests: a twin study. J Appl Psych.
- Rounds & Su (2014). Vocational interests. J Counseling Psych.
- Munafo et al. (2008). DRD4 and novelty seeking meta-analysis. Biol Psychiatry.
- Kosfeld et al. (2005). OXTR rs53576 and trust. Nature.
"""

from typing import Dict, List, Tuple, Any
from dataclasses import dataclass, field


CAREER_DISCLAIMER = (
    "Genes do NOT determine your career. This report describes statistical "
    "TENDENCIES toward broad work-style categories observed in research "
    "populations. Education, opportunity, training, interest, family context, "
    "luck, and personal choice are decisive — not the variants below. Use this "
    "as a curiosity-driven self-reflection prompt, not a vocational test."
)


# ─────────────────────────────────────────────────────────────────────────────
# COMPLEMENT MAP (strand-flip handling, mirrors other analyzers)
# ─────────────────────────────────────────────────────────────────────────────

COMPLEMENT = {"A": "T", "T": "A", "C": "G", "G": "C"}


def complement_allele(allele: str) -> str:
    return COMPLEMENT.get(allele.upper(), allele.upper())


def count_allele(genotype: str, allele: str) -> int:
    """Count occurrences of an allele in a genotype, handling complement strands."""
    genotype = genotype.upper().replace("-", "")
    allele = allele.upper()
    comp = complement_allele(allele)
    count = genotype.count(allele)
    if count == 0 and comp != allele:
        count = genotype.count(comp)
    return min(count, 2)


# ─────────────────────────────────────────────────────────────────────────────
# UNDERLYING TRAIT DIMENSIONS
# ─────────────────────────────────────────────────────────────────────────────
# Each SNP pushes one of these intermediate trait scores HIGHER (favoured allele
# present) or LOWER (allele absent). Trait scores are then mapped to RIASEC.

TRAIT_DIMS = [
    "openness",
    "conscientiousness",
    "extraversion",
    "agreeableness",
    "emotional_stability",   # inverse of neuroticism
    "cognitive_ability",
    "working_memory",
    "processing_speed",
    "verbal_ability",
    "risk_tolerance",
    "social_bonding",
    "creativity",
]


@dataclass
class CareerSNP:
    """A SNP whose genotype contributes to one trait dimension."""
    rsid: str
    gene: str
    trait: str               # which TRAIT_DIM this contributes to
    higher_allele: str       # allele that pushes the trait HIGHER
    lower_allele: str
    weight: float            # how strongly this SNP contributes (typically 0.5–1.0)
    note: str                # short description for the report
    evidence: str            # "Strong" | "Moderate" | "Preliminary"


# Curated SNPs — overlapping with Big Five + Intelligence panels for consistency,
# plus a few specifically chosen for risk-taking and social-bonding signals
# relevant to RIASEC.
CAREER_SNPS: List[CareerSNP] = [
    # ── OPENNESS (curiosity, exploration) ──────────────────────────────────
    CareerSNP("rs1800955", "DRD4",  "openness",  "T", "C", 1.0,
              "DRD4 promoter variant linked to novelty seeking",
              "Moderate"),
    CareerSNP("rs17518584", "CADM2", "openness", "T", "C", 1.0,
              "CADM2 GWAS variant for risk-taking and exploration",
              "Strong"),
    CareerSNP("rs4818",    "COMT",  "openness",  "C", "G", 0.5,
              "COMT haplotype tied to cognitive flexibility",
              "Moderate"),

    # ── CONSCIENTIOUSNESS (self-discipline, planning) ──────────────────────
    CareerSNP("rs1076560", "DRD2",     "conscientiousness", "G", "T", 1.0,
              "DRD2 splicing variant tied to self-control",
              "Moderate"),
    CareerSNP("rs1800497", "DRD2/ANKK1", "conscientiousness", "G", "A", 1.0,
              "Taq1A — D2 receptor density and impulse control",
              "Moderate"),

    # ── EXTRAVERSION (sociability, assertiveness) ──────────────────────────
    CareerSNP("rs1611115", "DBH",   "extraversion", "C", "T", 1.0,
              "Dopamine beta-hydroxylase activity & arousal",
              "Moderate"),
    CareerSNP("rs6277",    "DRD2",  "extraversion", "T", "C", 0.5,
              "DRD2 C957T — striatal reward signaling",
              "Moderate"),

    # ── AGREEABLENESS (empathy, prosociality) ──────────────────────────────
    CareerSNP("rs53576",   "OXTR",   "agreeableness", "G", "A", 1.0,
              "Oxytocin receptor — empathy and trust",
              "Strong"),
    CareerSNP("rs2254298", "OXTR",   "agreeableness", "G", "A", 0.5,
              "OXTR intronic variant — attachment behaviour",
              "Moderate"),

    # ── EMOTIONAL STABILITY (low neuroticism — note inverted polarity) ────
    # Higher-stability allele = lower neuroticism
    CareerSNP("rs6265",    "BDNF",   "emotional_stability", "G", "A", 1.0,
              "BDNF Val66Met — Val/Val supports stress resilience",
              "Moderate"),
    CareerSNP("rs4680",    "COMT",   "emotional_stability", "G", "A", 0.5,
              "COMT Val158Met — Val variant linked to stress resilience",
              "Moderate"),

    # ── COGNITIVE ABILITY (general g) ──────────────────────────────────────
    CareerSNP("rs2251499", "GWAS",   "cognitive_ability", "C", "T", 1.0,
              "Intelligence GWAS locus (Savage et al. 2018, n≈248k)",
              "Moderate"),
    CareerSNP("rs4851266", "LINC01104","cognitive_ability","T","C",1.0,
              "Educational-attainment GWAS variant",
              "Preliminary"),
    CareerSNP("rs429358",  "APOE",   "cognitive_ability", "T", "C", 0.5,
              "APOE non-E4 — preserved cognition with age",
              "Moderate"),

    # ── WORKING MEMORY ──────────────────────────────────────────────────────
    CareerSNP("rs4680",    "COMT",   "working_memory", "A", "G", 1.0,
              "COMT Val158Met — Met variant boosts prefrontal dopamine",
              "Moderate"),
    CareerSNP("rs17070145","KIBRA",  "working_memory", "T", "C", 1.0,
              "KIBRA/WWC1 — episodic and working memory",
              "Moderate"),

    # ── PROCESSING SPEED ───────────────────────────────────────────────────
    CareerSNP("rs1800497", "DRD2/ANKK1","processing_speed", "G", "A", 1.0,
              "D2 receptor density and reaction time",
              "Moderate"),

    # ── VERBAL ABILITY ─────────────────────────────────────────────────────
    CareerSNP("rs2253478", "FOXP2",  "verbal_ability", "G", "A", 1.0,
              "FOXP2 — speech and language regulator",
              "Moderate"),
    CareerSNP("rs7794745", "CNTNAP2","verbal_ability", "A", "T", 1.0,
              "CNTNAP2 — language acquisition trajectory",
              "Moderate"),

    # ── RISK TOLERANCE (sensation seeking, exploration) ────────────────────
    CareerSNP("rs17518584","CADM2",  "risk_tolerance", "T", "C", 1.0,
              "CADM2 — risk-taking GWAS (Sanchez-Roige 2018)",
              "Strong"),
    CareerSNP("rs1800955", "DRD4",   "risk_tolerance", "T", "C", 0.5,
              "DRD4 — novelty seeking and behavioral risk",
              "Moderate"),

    # ── SOCIAL BONDING (attachment, cooperation) ────────────────────────────
    CareerSNP("rs53576",   "OXTR",   "social_bonding", "G", "A", 1.0,
              "OXTR — empathic accuracy and pair bonding",
              "Strong"),
    CareerSNP("rs11174811","AVPR1A", "social_bonding", "A", "C", 1.0,
              "AVPR1A — pair-bonding & altruistic behaviour",
              "Moderate"),

    # ── CREATIVITY (divergent thinking proxy) ──────────────────────────────
    CareerSNP("rs1800955", "DRD4",   "creativity", "T", "C", 1.0,
              "DRD4 — exploration and novelty motivation",
              "Moderate"),
    CareerSNP("rs6265",    "BDNF",   "creativity", "A", "G", 0.5,
              "BDNF Met — divergent thinking studies (mixed)",
              "Preliminary"),
]


# ─────────────────────────────────────────────────────────────────────────────
# RIASEC (Holland) DIMENSIONS
# ─────────────────────────────────────────────────────────────────────────────

RIASEC = ["Realistic", "Investigative", "Artistic", "Social",
          "Enterprising", "Conventional"]

RIASEC_EMOJI = {
    "Realistic":     "\U0001f527",   # 🔧
    "Investigative": "\U0001f52c",   # 🔬
    "Artistic":      "\U0001f3a8",   # 🎨
    "Social":        "\U0001f91d",   # 🤝
    "Enterprising":  "\U0001f680",   # 🚀
    "Conventional":  "\U0001f4ca",   # 📊
}

RIASEC_TAGLINE = {
    "Realistic":     "The Doer — hands-on, practical, results-oriented",
    "Investigative": "The Thinker — analytical, curious, evidence-driven",
    "Artistic":      "The Creator — expressive, original, idea-rich",
    "Social":        "The Helper — empathic, supportive, people-first",
    "Enterprising":  "The Persuader — bold, leadership-driven, opportunity-seeking",
    "Conventional":  "The Organizer — methodical, detail-oriented, structured",
}

# Weights mapping each underlying trait to each RIASEC score.
# Positive = trait pushes RIASEC dimension UP. Negative = pushes DOWN.
# Weights are normalized internally so absolute magnitudes don't matter much —
# what matters is the relative profile across RIASEC.
RIASEC_WEIGHTS: Dict[str, Dict[str, float]] = {
    "Realistic": {
        "conscientiousness":   0.30,
        "openness":           -0.15,
        "extraversion":        0.05,
        "emotional_stability": 0.20,
        "processing_speed":    0.20,
        "risk_tolerance":      0.20,
        "verbal_ability":     -0.05,
    },
    "Investigative": {
        "cognitive_ability":   0.40,
        "openness":            0.25,
        "conscientiousness":   0.15,
        "working_memory":      0.25,
        "agreeableness":      -0.05,
        "emotional_stability": 0.05,
    },
    "Artistic": {
        "openness":            0.40,
        "creativity":          0.30,
        "conscientiousness":  -0.20,
        "verbal_ability":      0.10,
        "emotional_stability":-0.10,
        "risk_tolerance":      0.10,
    },
    "Social": {
        "agreeableness":       0.35,
        "extraversion":        0.25,
        "social_bonding":      0.30,
        "verbal_ability":      0.10,
        "openness":            0.05,
    },
    "Enterprising": {
        "extraversion":        0.35,
        "risk_tolerance":      0.30,
        "conscientiousness":   0.10,
        "emotional_stability": 0.20,
        "agreeableness":      -0.10,
        "verbal_ability":      0.10,
    },
    "Conventional": {
        "conscientiousness":   0.40,
        "openness":           -0.20,
        "emotional_stability": 0.20,
        "working_memory":      0.20,
        "agreeableness":       0.05,
    },
}

# Example careers per dimension (kept generic; not deterministic suggestions)
CAREER_EXAMPLES: Dict[str, List[str]] = {
    "Realistic": [
        "Engineer", "Surgeon", "Mechanic", "Architect", "Pilot",
        "Athlete / Coach", "Veterinarian", "Field biologist",
        "Skilled tradesperson", "Firefighter",
    ],
    "Investigative": [
        "Scientist", "Physician (diagnostic specialties)", "Data analyst",
        "Software engineer (R&D)", "Research economist", "Mathematician",
        "Forensic analyst", "Epidemiologist", "Quant trader",
        "University researcher",
    ],
    "Artistic": [
        "Writer", "Designer", "Musician", "Filmmaker", "Photographer",
        "Creative director", "UX/UI designer", "Architect (creative)",
        "Game designer", "Art therapist",
    ],
    "Social": [
        "Teacher", "Counselor / Therapist", "Nurse", "Coach",
        "Social worker", "HR partner", "Pediatrician", "Community manager",
        "Speech-language pathologist", "Non-profit leader",
    ],
    "Enterprising": [
        "Entrepreneur / Founder", "Sales executive", "Lawyer",
        "Investment banker", "Marketing director", "Politician",
        "Realtor", "Talent agent", "Business development",
        "Strategy consultant",
    ],
    "Conventional": [
        "Accountant", "Auditor", "Financial analyst", "Operations manager",
        "Project manager", "Database administrator", "Compliance officer",
        "Logistics coordinator", "Tax specialist", "Quality assurance",
    ],
}


# ─────────────────────────────────────────────────────────────────────────────
# RESULT TYPES
# ─────────────────────────────────────────────────────────────────────────────

@dataclass
class CareerVariantResult:
    rsid: str
    gene: str
    trait: str
    genotype: str
    higher_allele_count: int
    note: str
    evidence: str


@dataclass
class CareerAnalysisResult:
    total_checked: int
    found: int
    not_found: int
    trait_scores: Dict[str, float]           # 0–100 per underlying trait
    riasec_scores: Dict[str, float]          # 0–100 per RIASEC dimension
    riasec_ranked: List[Tuple[str, float]]   # sorted descending
    holland_code: str                        # top 3 letters, e.g. "IRE"
    variants: List[CareerVariantResult]
    missing_rsids: List[str]


# ─────────────────────────────────────────────────────────────────────────────
# CORE ANALYSIS
# ─────────────────────────────────────────────────────────────────────────────

def analyze_career(
    variants: Dict[str, Tuple[str, str, str]]
) -> CareerAnalysisResult:
    """Analyze genome data and produce RIASEC career-style scores.

    Args:
        variants: {rsid_lower: (chrom, pos, genotype)} as parsed from genome file.

    Returns:
        CareerAnalysisResult with trait + RIASEC scores and per-variant detail.
    """
    # 1. Tally allele counts per trait dimension, weighted by SNP weight.
    trait_numerator: Dict[str, float] = {t: 0.0 for t in TRAIT_DIMS}
    trait_denominator: Dict[str, float] = {t: 0.0 for t in TRAIT_DIMS}

    variant_results: List[CareerVariantResult] = []
    missing: List[str] = []
    found = 0

    for snp in CAREER_SNPS:
        rsid_lc = snp.rsid.lower()
        # Each SNP can contribute up to 2 favourable allele copies.
        trait_denominator[snp.trait] += 2.0 * snp.weight

        if rsid_lc in variants:
            found += 1
            _chrom, _pos, gt = variants[rsid_lc]
            gt = gt.upper().replace("-", "")
            high_count = count_allele(gt, snp.higher_allele)
            trait_numerator[snp.trait] += high_count * snp.weight

            variant_results.append(CareerVariantResult(
                rsid=snp.rsid,
                gene=snp.gene,
                trait=snp.trait,
                genotype=gt,
                higher_allele_count=high_count,
                note=snp.note,
                evidence=snp.evidence,
            ))
        else:
            missing.append(f"{snp.rsid} ({snp.gene})")

    # 2. Compute trait scores 0–100. If no data found for a trait → 50 (neutral).
    trait_scores: Dict[str, float] = {}
    for t in TRAIT_DIMS:
        denom = trait_denominator[t]
        if denom > 0:
            num = trait_numerator[t]
            # Only report a real score if at least ONE SNP for this trait was
            # found in the genome — otherwise it's all "missing", default to
            # neutral 50.
            any_found = any(v.trait == t for v in variant_results)
            trait_scores[t] = round(num / denom * 100, 1) if any_found else 50.0
        else:
            trait_scores[t] = 50.0

    # 3. Compute RIASEC scores. Each dimension is the weighted sum of trait
    #    scores (centered around 50), normalized back to a 0–100 scale.
    riasec_scores: Dict[str, float] = {}
    for dim in RIASEC:
        weights = RIASEC_WEIGHTS[dim]
        numerator = 0.0
        weight_sum = 0.0
        for trait, w in weights.items():
            # Center trait score around 50 (i.e. positive deviation = above avg)
            # then weight signed contribution.
            delta = trait_scores.get(trait, 50.0) - 50.0
            numerator += delta * w
            weight_sum += abs(w)
        # Normalize so an "all maxed" profile maps to ~100, "all min" to ~0.
        # delta range is [-50, +50] per trait; weighted sum range is
        # [-50*weight_sum, +50*weight_sum]. Map back to [0, 100].
        if weight_sum > 0:
            normalized = 50.0 + (numerator / weight_sum)
        else:
            normalized = 50.0
        # Clamp to [0,100]
        riasec_scores[dim] = round(max(0.0, min(100.0, normalized)), 1)

    # 4. Sort RIASEC descending and form Holland 3-letter code.
    riasec_ranked = sorted(riasec_scores.items(), key=lambda x: -x[1])
    holland_code = "".join(d[0] for d, _ in riasec_ranked[:3])

    return CareerAnalysisResult(
        total_checked=len(CAREER_SNPS),
        found=found,
        not_found=len(CAREER_SNPS) - found,
        trait_scores=trait_scores,
        riasec_scores=riasec_scores,
        riasec_ranked=riasec_ranked,
        holland_code=holland_code,
        variants=variant_results,
        missing_rsids=missing,
    )


# ─────────────────────────────────────────────────────────────────────────────
# JSON OUTPUT (consumed by frontend React component)
# ─────────────────────────────────────────────────────────────────────────────

def _riasec_interpretation(dim: str, score: float) -> str:
    if score >= 70:
        level = "strong"
    elif score >= 60:
        level = "moderately strong"
    elif score >= 45:
        level = "balanced"
    elif score >= 30:
        level = "moderately weak"
    else:
        level = "weak"

    descriptions = {
        "Realistic": {
            "strong": "Your profile leans clearly toward hands-on, results-driven work where outcomes are tangible.",
            "moderately strong": "Your profile gently favours practical, action-oriented work over abstract or socially-driven roles.",
            "balanced": "Hands-on practical work is neither strongly indicated nor contraindicated by your profile.",
            "moderately weak": "Your profile leans away from hands-on, outcome-driven work in favor of more abstract or interpersonal styles.",
            "weak": "Your profile suggests hands-on, mechanical work is a poorer match than other styles.",
        },
        "Investigative": {
            "strong": "Your profile strongly favours analytical, evidence-driven work — research, diagnosis, deep technical problem-solving.",
            "moderately strong": "You lean toward analytical work that rewards curiosity and rigorous thinking.",
            "balanced": "Analytical work is one option among several that fit your profile equally.",
            "moderately weak": "Pure analytical / research work is a slightly weaker fit than other styles.",
            "weak": "Your profile suggests deep-dive analytical work may feel less energising than other domains.",
        },
        "Artistic": {
            "strong": "Your profile strongly favours creative, expressive, original work where there is no single right answer.",
            "moderately strong": "You lean toward creative work that rewards originality and personal voice.",
            "balanced": "Creative work is a viable but not standout direction.",
            "moderately weak": "Pure creative work is a slightly weaker fit than other styles.",
            "weak": "Your profile suggests open-ended creative work may feel less motivating than structured or social work.",
        },
        "Social": {
            "strong": "Your profile strongly favours people-centered work — teaching, healing, advising, building communities.",
            "moderately strong": "You lean toward people-centered roles where empathy and relationships are central.",
            "balanced": "Social work is a viable option among several that fit equally.",
            "moderately weak": "Highly people-centric roles are a slightly weaker fit than other styles.",
            "weak": "Your profile suggests heavily people-centric work may feel more draining than energising.",
        },
        "Enterprising": {
            "strong": "Your profile strongly favours leadership, persuasion, and opportunity-seeking work — sales, founding, strategy, advocacy.",
            "moderately strong": "You lean toward bold, persuasive, leadership-driven roles.",
            "balanced": "Leadership / persuasion work is a viable option among several.",
            "moderately weak": "Highly persuasive or leadership-heavy work is a slightly weaker fit.",
            "weak": "Your profile suggests opportunity-driven, persuasion-heavy work may feel less natural than other styles.",
        },
        "Conventional": {
            "strong": "Your profile strongly favours structured, detail-oriented, methodical work — finance, ops, compliance, logistics.",
            "moderately strong": "You lean toward methodical, organized, detail-oriented work.",
            "balanced": "Structured / methodical work is a viable option among several.",
            "moderately weak": "Highly structured rule-based work is a slightly weaker fit than other styles.",
            "weak": "Your profile suggests highly procedural work may feel constraining compared to other styles.",
        },
    }
    return descriptions.get(dim, {}).get(level, "Balanced profile.")


def generate_career_genes_json(result: CareerAnalysisResult) -> dict:
    """Build the JSON payload consumed by the frontend Career Genes report."""
    riasec_list = []
    for dim in RIASEC:
        score = result.riasec_scores.get(dim, 50.0)
        riasec_list.append({
            "name": dim,
            "emoji": RIASEC_EMOJI[dim],
            "tagline": RIASEC_TAGLINE[dim],
            "score": score,
            "interpretation": _riasec_interpretation(dim, score),
            "exampleCareers": CAREER_EXAMPLES[dim],
        })

    # Top fits — first 2 from ranked
    top_fits = []
    for dim, score in result.riasec_ranked[:2]:
        top_fits.append({
            "name": dim,
            "emoji": RIASEC_EMOJI[dim],
            "tagline": RIASEC_TAGLINE[dim],
            "score": score,
            "exampleCareers": CAREER_EXAMPLES[dim][:6],
        })

    # Variant breakdown grouped by trait
    variants_json = []
    for v in result.variants:
        variants_json.append({
            "rsid": v.rsid,
            "gene": v.gene,
            "trait": v.trait,
            "genotype": v.genotype,
            "higherAlleleCount": v.higher_allele_count,
            "note": v.note,
            "evidence": v.evidence,
        })

    return {
        "summary": {
            "hollandCode": result.holland_code,
            "topFits": top_fits,
            "totalChecked": result.total_checked,
            "found": result.found,
            "notFound": result.not_found,
        },
        "riasec": riasec_list,
        "traitScores": result.trait_scores,
        "variants": variants_json,
        "missingRsids": result.missing_rsids,
        "disclaimer": CAREER_DISCLAIMER,
        "references": [
            "Holland, J. L. (1997). Making Vocational Choices, 3rd ed.",
            "Plomin & Deary (2015). Genetics and intelligence differences. Mol Psychiatry.",
            "Lee et al. (2018). Educational attainment GWAS (n≈1.1M). Nat Genet.",
            "Sanchez-Roige et al. (2018). Risk-taking GWAS. Nat Hum Behav.",
            "Schermer & Vernon (2008). Heritability of vocational interests. Behav Genet.",
            "Munafo et al. (2008). DRD4 and novelty seeking meta-analysis. Biol Psychiatry.",
            "Kosfeld et al. (2005). OXTR and trust. Nature.",
        ],
    }


# ─────────────────────────────────────────────────────────────────────────────
# MARKDOWN OUTPUT (optional — used by some downstream pipelines)
# ─────────────────────────────────────────────────────────────────────────────

def generate_career_genes_report(
    result: CareerAnalysisResult, subject_name: str = "Subject"
) -> str:
    from datetime import datetime
    out = []
    out.append(f"# Career Genes Report\n")
    out.append(f"**Subject:** {subject_name}  \n")
    out.append(f"**Generated:** {datetime.now().strftime('%B %d, %Y')}\n\n")
    out.append("---\n\n")
    out.append(f"> {CAREER_DISCLAIMER}\n\n")
    out.append(f"## Your Holland Code: **{result.holland_code}**\n\n")

    out.append("### RIASEC Profile\n\n")
    out.append("| Style | Score |\n|---|---|\n")
    for dim, score in result.riasec_ranked:
        out.append(f"| {RIASEC_EMOJI[dim]} {dim} — {RIASEC_TAGLINE[dim]} | {score}/100 |\n")

    out.append("\n---\n\n")
    out.append("### Top Fits\n\n")
    for dim, score in result.riasec_ranked[:2]:
        out.append(f"**{RIASEC_EMOJI[dim]} {dim} ({score}/100)** — {RIASEC_TAGLINE[dim]}\n\n")
        out.append(f"{_riasec_interpretation(dim, score)}\n\n")
        examples = ", ".join(CAREER_EXAMPLES[dim][:6])
        out.append(f"_Example careers in this dimension:_ {examples}.\n\n")

    out.append("\n### Important Notes\n")
    out.append("1. RIASEC describes work styles, not job titles.\n")
    out.append("2. Heritability of vocational interests is ~0.30–0.50 — environment dominates the rest.\n")
    out.append("3. Per-SNP effects are tiny; this is a coarse, curiosity-driven sketch.\n")
    out.append("4. Validated career counseling instruments (Strong Interest Inventory, Self-Directed Search) remain the gold standard for vocational guidance.\n")

    return "".join(out)
