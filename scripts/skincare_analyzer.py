"""
Skincare & Dermatogenomics Analyzer
Analyzes genetic variants for skin pigmentation, UV sensitivity,
collagen integrity, skin aging, and antioxidant defense.

Returns structured JSON for frontend visualization.

References:
- Raimondi et al. (2008) Eur J Cancer 44:2275-84
- Lamason et al. (2005) Science 310:1782-6
- Sturm (2009) Hum Genet: Molecular genetics of human pigmentation
- Liu et al. (2015) Forensic Sci Int Genet: HIrisPlex-S skin color prediction
- Cole et al. (2019) Aging: Genomic predictors of skin aging

DISCLAIMER:
  This is NOT medical dermatological advice. Results reflect genetic
  tendencies based on population studies. Consult a dermatologist for
  medical skin concerns.
"""

from typing import Dict, Tuple, Optional, List, Any


# ─────────────────────────────────────────────────────────────────────────────
# COMPLEMENT STRAND MAP
# ─────────────────────────────────────────────────────────────────────────────

COMPLEMENT = {"A": "T", "T": "A", "C": "G", "G": "C"}


# ─────────────────────────────────────────────────────────────────────────────
# SKINCARE GENETIC VARIANT PANEL
# ─────────────────────────────────────────────────────────────────────────────

SKINCARE_VARIANTS = {

    # ── PIGMENTATION & UV SENSITIVITY ────────────────────────────────────────

    "rs1805007": {
        "gene": "MC1R", "name": "R151C — Major red hair/fair skin variant",
        "category": "Pigmentation & UV Sensitivity",
        "ref": "C", "alt": "T",
        "genotypes": {
            "CC": {"label": "Wild-type", "uv_risk": 0, "aging_risk": 0,
                   "meaning": "No R151C variant. Normal melanocortin-1 receptor function."},
            "CT": {"label": "Carrier (R/r)", "uv_risk": 2, "aging_risk": 0,
                   "meaning": "One copy of R151C. Partial MC1R loss-of-function. "
                              "Increased fair skin, freckling, and sun sensitivity. "
                              "~2-4x increased melanoma risk per allele."},
            "TT": {"label": "Homozygous variant (r/r)", "uv_risk": 3, "aging_risk": 0,
                   "meaning": "Both copies of R151C. Strong shift toward pheomelanin "
                              "over eumelanin. Very fair skin, red/auburn hair, "
                              "high UV sensitivity, ~6-10x melanoma risk."},
        },
        "evidence": "Strong",
    },

    "rs1805008": {
        "gene": "MC1R", "name": "R160W — Red hair/fair skin variant",
        "category": "Pigmentation & UV Sensitivity",
        "ref": "C", "alt": "T",
        "genotypes": {
            "CC": {"label": "Wild-type", "uv_risk": 0, "aging_risk": 0,
                   "meaning": "Normal MC1R function at this position."},
            "CT": {"label": "Carrier (R160W heterozygous)", "uv_risk": 1, "aging_risk": 0,
                   "meaning": "One copy of R160W. Mild increase in sun sensitivity and freckling."},
            "TT": {"label": "Homozygous R160W", "uv_risk": 3, "aging_risk": 0,
                   "meaning": "Both copies. Strong fair skin and red hair association. "
                              "High UV sensitivity."},
        },
        "evidence": "Strong",
    },

    "rs1805009": {
        "gene": "MC1R", "name": "D294H — Red hair variant",
        "category": "Pigmentation & UV Sensitivity",
        "ref": "G", "alt": "C",
        "genotypes": {
            "GG": {"label": "Wild-type", "uv_risk": 0, "aging_risk": 0,
                   "meaning": "Normal MC1R function at this position."},
            "GC": {"label": "Carrier D294H", "uv_risk": 2, "aging_risk": 0,
                   "meaning": "One copy. Associated with red hair, fair skin, UV sensitivity."},
            "CC": {"label": "Homozygous D294H", "uv_risk": 3, "aging_risk": 0,
                   "meaning": "Both copies. Strong red hair/fair skin association. "
                              "Very high UV sensitivity."},
        },
        "evidence": "Strong",
    },

    "rs1126809": {
        "gene": "TYR", "name": "Tyrosinase R402Q",
        "category": "Pigmentation & UV Sensitivity",
        "ref": "G", "alt": "A",
        "genotypes": {
            "GG": {"label": "Wild-type", "uv_risk": 0, "aging_risk": 0,
                   "meaning": "Normal tyrosinase activity. Standard melanin production."},
            "GA": {"label": "Heterozygous R402Q", "uv_risk": 1, "aging_risk": 0,
                   "meaning": "Slightly reduced tyrosinase activity. Modestly lighter "
                              "skin pigmentation."},
            "AA": {"label": "Homozygous R402Q", "uv_risk": 2, "aging_risk": 0,
                   "meaning": "Reduced tyrosinase activity. Lighter pigmentation, "
                              "increased sun sensitivity."},
        },
        "evidence": "Strong",
    },

    "rs16891982": {
        "gene": "SLC45A2", "name": "Membrane-Associated Transporter Protein L374F",
        "category": "Pigmentation & UV Sensitivity",
        "ref": "C", "alt": "G",
        "genotypes": {
            "GG": {"label": "European light pigmentation (L374)", "uv_risk": 1, "aging_risk": 0,
                   "meaning": "Derived allele near-fixed in Europeans (~97%). Associated with "
                              "lighter skin, hair, and eye color. Moderate UV sensitivity."},
            "CG": {"label": "Heterozygous", "uv_risk": 0, "aging_risk": 0,
                   "meaning": "Mixed pigmentation effect. Intermediate skin tone."},
            "CC": {"label": "Ancestral (F374)", "uv_risk": -1, "aging_risk": 0,
                   "meaning": "Original darker pigmentation allele. Better natural UV protection."},
        },
        "evidence": "Strong",
    },

    "rs1426654": {
        "gene": "SLC24A5", "name": "Skin pigmentation (A111T)",
        "category": "Pigmentation & UV Sensitivity",
        "ref": "G", "alt": "A",
        "genotypes": {
            "AA": {"label": "Derived — Light skin (European)", "uv_risk": 1, "aging_risk": 0,
                   "meaning": "Near-fixed in Europeans (~99%). Major skin-lightening allele. "
                              "Explains ~25-38% of skin color difference between Europeans "
                              "and Africans."},
            "GA": {"label": "Heterozygous", "uv_risk": 0, "aging_risk": 0,
                   "meaning": "Intermediate pigmentation effect."},
            "GG": {"label": "Ancestral — Darker pigmentation", "uv_risk": -1, "aging_risk": 0,
                   "meaning": "Original allele. Better natural UV protection. Near-fixed in "
                              "African populations."},
        },
        "evidence": "Strong",
    },

    "rs12203592": {
        "gene": "IRF4", "name": "Interferon Regulatory Factor 4",
        "category": "Pigmentation & UV Sensitivity",
        "ref": "C", "alt": "T",
        "genotypes": {
            "CC": {"label": "No freckling allele", "uv_risk": 0, "aging_risk": 0,
                   "meaning": "Standard pigmentation response."},
            "CT": {"label": "Freckling carrier", "uv_risk": 1, "aging_risk": 0,
                   "meaning": "T allele associated with freckling, lighter hair, sun "
                              "sensitivity. IRF4 regulates melanocyte development."},
            "TT": {"label": "Homozygous freckling", "uv_risk": 2, "aging_risk": 0,
                   "meaning": "Strong freckling tendency, lighter hair, increased sun "
                              "sensitivity."},
        },
        "evidence": "Strong",
    },

    # ── SKIN AGING & COLLAGEN ────────────────────────────────────────────────

    "rs2228145": {
        "gene": "IL6R", "name": "Interleukin-6 Receptor (D358A)",
        "category": "Skin Aging & Collagen",
        "ref": "A", "alt": "C",
        "genotypes": {
            "AA": {"label": "Asp/Asp — Standard IL-6 signaling", "uv_risk": 0, "aging_risk": 0,
                   "meaning": "Normal inflammatory signaling through IL-6 receptor."},
            "AC": {"label": "Asp/Ala — Reduced membrane IL-6R", "uv_risk": 0, "aging_risk": -1,
                   "meaning": "One copy of D358A. Higher soluble IL-6R levels. May reduce "
                              "local skin inflammation and slow inflammatory aging."},
            "CC": {"label": "Ala/Ala — Reduced membrane IL-6R", "uv_risk": 0, "aging_risk": -1,
                   "meaning": "Both copies. Significantly higher soluble IL-6R. Reduced "
                              "classic IL-6 signaling. May protect against inflammation-driven "
                              "skin aging."},
        },
        "evidence": "Moderate",
    },

    "rs1800629": {
        "gene": "TNF", "name": "TNF-alpha Promoter (-308 G/A)",
        "category": "Skin Aging & Collagen",
        "ref": "G", "alt": "A",
        "genotypes": {
            "GG": {"label": "G/G — Standard TNF-alpha", "uv_risk": 0, "aging_risk": 0,
                   "meaning": "Normal TNF-alpha production. Standard inflammatory response."},
            "GA": {"label": "G/A — Higher TNF-alpha", "uv_risk": 0, "aging_risk": 1,
                   "meaning": "One A allele. ~2-fold higher TNF-alpha production. May "
                              "experience more skin inflammation, redness, and reactive "
                              "skin conditions."},
            "AA": {"label": "A/A — Highest TNF-alpha", "uv_risk": 0, "aging_risk": 2,
                   "meaning": "Both A alleles. Highest TNF-alpha production. Increased risk "
                              "of inflammatory skin conditions, faster inflamm-aging."},
        },
        "evidence": "Moderate",
    },

    # ── SKIN HEALTH & VITAMIN D ──────────────────────────────────────────────

    "rs1544410": {
        "gene": "VDR", "name": "Vitamin D Receptor (BsmI)",
        "category": "Skin Health & Vitamin D",
        "ref": "C", "alt": "T",
        "genotypes": {
            "CC": {"label": "b/b — Higher VDR expression", "uv_risk": 0, "aging_risk": -1,
                   "meaning": "Higher vitamin D receptor expression in skin. Better vitamin D "
                              "signaling may support skin barrier function and immune defense."},
            "CT": {"label": "B/b — Intermediate", "uv_risk": 0, "aging_risk": 0,
                   "meaning": "Intermediate VDR expression."},
            "TT": {"label": "B/B — Lower VDR expression", "uv_risk": 0, "aging_risk": 1,
                   "meaning": "Lower VDR activity. May need more sun/supplementation for "
                              "adequate vitamin D. Potentially reduced skin barrier support."},
        },
        "evidence": "Moderate",
    },

    "rs2282679": {
        "gene": "GC", "name": "Vitamin D Binding Protein",
        "category": "Skin Health & Vitamin D",
        "ref": "A", "alt": "C",
        "genotypes": {
            "AA": {"label": "Higher vitamin D levels", "uv_risk": 0, "aging_risk": -1,
                   "meaning": "Associated with higher circulating 25(OH)D levels. Better "
                              "vitamin D status supports skin health."},
            "AC": {"label": "Intermediate vitamin D", "uv_risk": 0, "aging_risk": 0,
                   "meaning": "Intermediate vitamin D levels."},
            "CC": {"label": "Lower vitamin D levels", "uv_risk": 0, "aging_risk": 1,
                   "meaning": "Associated with lower circulating vitamin D. May need "
                              "supplementation, especially in northern climates."},
        },
        "evidence": "Strong",
    },

    # ── ANTIOXIDANT DEFENSE ──────────────────────────────────────────────────

    "rs4880": {
        "gene": "SOD2", "name": "Superoxide Dismutase 2 (Ala16Val)",
        "category": "Antioxidant Defense",
        "ref": "C", "alt": "T",
        "genotypes": {
            "CC": {"label": "Ala/Ala — Efficient mitochondrial SOD", "uv_risk": 0, "aging_risk": -1,
                   "meaning": "Optimal mitochondrial targeting of SOD2. Better defense against "
                              "reactive oxygen species (ROS). May slow oxidative skin aging."},
            "CT": {"label": "Ala/Val — Intermediate", "uv_risk": 0, "aging_risk": 0,
                   "meaning": "Intermediate SOD2 efficiency."},
            "TT": {"label": "Val/Val — Reduced SOD transport", "uv_risk": 0, "aging_risk": 1,
                   "meaning": "Val16 variant reduces mitochondrial import efficiency of SOD2. "
                              "Less efficient ROS neutralization. May benefit more from "
                              "antioxidant-rich skincare and diet."},
        },
        "evidence": "Moderate",
    },

    "rs1001179": {
        "gene": "CAT", "name": "Catalase Promoter (-262 C/T)",
        "category": "Antioxidant Defense",
        "ref": "C", "alt": "T",
        "genotypes": {
            "CC": {"label": "C/C — Standard catalase", "uv_risk": 0, "aging_risk": 0,
                   "meaning": "Normal catalase expression. Standard hydrogen peroxide "
                              "detoxification."},
            "CT": {"label": "C/T — Higher catalase activity", "uv_risk": 0, "aging_risk": -1,
                   "meaning": "T allele associated with higher catalase promoter activity. "
                              "Better H2O2 neutralization. May protect against oxidative "
                              "skin damage."},
            "TT": {"label": "T/T — Highest catalase", "uv_risk": 0, "aging_risk": -1,
                   "meaning": "Highest catalase expression. Best antioxidant defense against "
                              "hydrogen peroxide. May offer better protection against "
                              "photo-aging."},
        },
        "evidence": "Emerging",
    },
}

# Ordered category list for consistent output
SKINCARE_CATEGORIES = [
    "Pigmentation & UV Sensitivity",
    "Skin Aging & Collagen",
    "Skin Health & Vitamin D",
    "Antioxidant Defense",
]

CATEGORY_EMOJIS = {
    "Pigmentation & UV Sensitivity": "☀️",
    "Skin Aging & Collagen": "🧬",
    "Skin Health & Vitamin D": "💊",
    "Antioxidant Defense": "🛡️",
}

# ─────────────────────────────────────────────────────────────────────────────
# SKINCARE ADVICE PROFILES
# ─────────────────────────────────────────────────────────────────────────────

ADVICE = {
    "sun_protection": {
        "title": "Sun Protection",
        "low": "Your pigmentation genetics suggest moderate natural UV protection. "
               "Standard sun protection (SPF 30+, reapply every 2h) is recommended.",
        "moderate": "You carry one or more variants associated with increased UV sensitivity. "
                    "Use SPF 50+ broad-spectrum daily, especially on face and hands. "
                    "Seek shade during peak UV hours (10am-4pm). Regular skin checks recommended.",
        "high": "Multiple variants indicate high UV sensitivity and fair skin tendency. "
                "DAILY SPF 50+ is essential, even on cloudy days. UPF clothing recommended "
                "for extended outdoor exposure. Annual dermatologist skin check strongly advised. "
                "Vitamin D supplementation may be needed if avoiding sun exposure.",
    },
    "antiaging": {
        "title": "Anti-Aging Strategy",
        "low": "Your collagen/inflammation genetics appear favorable. Focus on prevention: "
               "sunscreen, retinol (start with 0.025%), and adequate hydration.",
        "moderate": "Some genetic factors may accelerate skin aging. Consider adding vitamin C "
                    "serum (antioxidant + collagen stimulation), retinol (0.025-0.05%), "
                    "and peptide-based products to your routine.",
        "high": "Genetic markers suggest accelerated collagen breakdown and/or inflammation. "
                "Prioritize: daily retinol, vitamin C serum (L-ascorbic acid 10-20%), "
                "niacinamide (vitamin B3) for barrier repair, and peptide treatments. "
                "Sun protection is even more critical for you.",
    },
    "anti_inflammation": {
        "title": "Anti-Inflammation Care",
        "low": "No significant inflammatory genetic markers detected. Standard gentle "
               "skincare routine is appropriate.",
        "moderate": "You carry variants associated with elevated inflammatory response. "
                    "Use gentle, fragrance-free products. Consider niacinamide (3-5%) "
                    "to calm skin. Avoid harsh exfoliants.",
        "high": "Multiple inflammatory markers detected. Your skin may be more reactive. "
                "Use barrier-supporting ceramide products, centella asiatica extracts, "
                "and anti-inflammatory ingredients (licorice root, green tea extract). "
                "Patch-test new products. Avoid fragrance and essential oils.",
    },
    "vitamin_d": {
        "title": "Vitamin D & Skin Health",
        "low": "Your vitamin D genetics appear favorable. Moderate sun exposure + diet "
               "should maintain adequate levels.",
        "moderate": "Some variants suggest lower vitamin D status. Consider supplementation "
                    "(1000-2000 IU/day) especially in winter months.",
        "high": "Genetic variants suggest lower vitamin D production/transport. "
                "Supplementation (2000-4000 IU/day with vitamin K2) recommended, "
                "especially in winter or if using high-SPF sunscreen daily. "
                "Have 25(OH)D levels checked annually.",
    },
    "antioxidant": {
        "title": "Antioxidant Support",
        "low": "Antioxidant defense genetics appear adequate. Include antioxidant-rich "
               "foods in diet (berries, green vegetables, nuts).",
        "moderate": "Some variants suggest reduced antioxidant capacity. Add vitamin C "
                    "serum to skincare. Increase dietary antioxidants.",
        "high": "Genetic variants suggest reduced antioxidant defense. Prioritize: "
                "topical vitamin C + E serum (synergistic), CoQ10, resveratrol. "
                "Diet rich in colorful vegetables, omega-3 fatty acids. "
                "Consider oral antioxidant supplements (consult healthcare provider).",
    },
}

# Map each variant to its advice key based on category
VARIANT_ADVICE_KEYS = {
    "rs1805007": "sun_protection",
    "rs1805008": "sun_protection",
    "rs1805009": "sun_protection",
    "rs1126809": "sun_protection",
    "rs16891982": "sun_protection",
    "rs1426654": "sun_protection",
    "rs12203592": "sun_protection",
    "rs2228145": "antiaging",
    "rs1800629": "anti_inflammation",
    "rs1544410": "vitamin_d",
    "rs2282679": "vitamin_d",
    "rs4880": "antioxidant",
    "rs1001179": "antioxidant",
}


# ─────────────────────────────────────────────────────────────────────────────
# GENOTYPE NORMALIZATION
# ─────────────────────────────────────────────────────────────────────────────

def _normalize_genotype(genotype: str, ref: str, alt: str) -> Optional[str]:
    """
    Normalize a genotype to ref/alt alleles, handling complement strand.

    Tries the genotype as-is first, then tries the complement of each allele.
    Returns a sorted genotype string (ref allele first) or None if unrecognizable.
    """
    if not genotype or len(genotype) < 2:
        return None

    genotype = genotype.upper().replace("-", "")
    if len(genotype) != 2:
        return None

    a1, a2 = genotype[0], genotype[1]
    alleles = {ref, alt}

    # Direct match
    if a1 in alleles and a2 in alleles:
        return "".join(sorted([a1, a2], key=lambda x: (x != ref, x)))

    # Complement strand match
    a1c = COMPLEMENT.get(a1, a1)
    a2c = COMPLEMENT.get(a2, a2)
    if a1c in alleles and a2c in alleles:
        return "".join(sorted([a1c, a2c], key=lambda x: (x != ref, x)))

    return None


# ─────────────────────────────────────────────────────────────────────────────
# ANALYSIS
# ─────────────────────────────────────────────────────────────────────────────

def analyze_skincare(
    variants: Dict[str, Tuple[str, str, str]]
) -> Dict[str, Any]:
    """
    Analyze genetic variants for skincare-related traits.

    Args:
        variants: Dictionary mapping rsID -> (chromosome, position, genotype).
                  Keys should be lowercase rsIDs.

    Returns:
        Dictionary with analysis results (variant interpretations, scores,
        advice scores) suitable for passing to generate_skincare_json().
    """
    results: List[Dict[str, Any]] = []
    advice_scores: Dict[str, int] = {k: 0 for k in ADVICE}
    found = 0
    not_found = 0

    for rsid, var_info in SKINCARE_VARIANTS.items():
        rsid_lower = rsid.lower()
        entry: Dict[str, Any] = {
            "rsid": rsid,
            "var": var_info,
            "found": False,
            "genotype": None,
            "interp": None,
        }

        if rsid_lower not in variants:
            not_found += 1
            results.append(entry)
            continue

        _chrom, _pos, raw_genotype = variants[rsid_lower]
        norm = _normalize_genotype(raw_genotype, var_info["ref"], var_info["alt"])

        if norm is None:
            not_found += 1
            entry["genotype"] = raw_genotype
            results.append(entry)
            continue

        found += 1
        entry["found"] = True
        entry["genotype"] = norm

        # Look up interpretation
        interp = var_info["genotypes"].get(norm)
        if interp is None:
            # Try reversed order
            rev = norm[1] + norm[0]
            interp = var_info["genotypes"].get(rev)

        entry["interp"] = interp

        # Accumulate advice scores
        if interp:
            advice_key = VARIANT_ADVICE_KEYS.get(rsid)
            if advice_key:
                uv = interp.get("uv_risk", 0)
                aging = interp.get("aging_risk", 0)
                advice_scores[advice_key] += uv + aging

        results.append(entry)

    return {
        "results": results,
        "advice_scores": advice_scores,
        "found": found,
        "not_found": not_found,
        "total": len(SKINCARE_VARIANTS),
    }


def _risk_level(score: int) -> str:
    """Convert numeric score to risk level string."""
    if score <= 0:
        return "low"
    elif score <= 3:
        return "moderate"
    else:
        return "high"


def _compute_uv_risk_score(results: List[Dict[str, Any]]) -> int:
    """
    Compute a 0-100 UV risk score based on pigmentation variants.

    The raw score is derived from uv_risk values across all found variants.
    Maximum theoretical raw score is ~17 (all worst-case pigmentation alleles).
    """
    raw = 0
    max_possible = 0

    for entry in results:
        if not entry["found"] or entry["interp"] is None:
            continue
        uv = entry["interp"].get("uv_risk", 0)
        # Only count pigmentation category for UV score
        if entry["var"]["category"] == "Pigmentation & UV Sensitivity":
            raw += uv
            # Max possible per variant: find the highest uv_risk in its genotypes
            max_uv = max(
                g.get("uv_risk", 0) for g in entry["var"]["genotypes"].values()
            )
            max_possible += max_uv

    if max_possible <= 0:
        return 0

    # Normalize to 0-100, clamped
    score = int(round((raw / max_possible) * 100))
    return max(0, min(100, score))


def _compute_aging_risk_score(results: List[Dict[str, Any]]) -> int:
    """
    Compute a 0-100 aging risk score based on aging-related variants.

    aging_risk values can be negative (protective) or positive (risk).
    We map the range so that 0 raw = 50 (neutral), negative = lower, positive = higher.
    """
    raw = 0
    min_possible = 0
    max_possible = 0

    for entry in results:
        if not entry["found"] or entry["interp"] is None:
            continue
        aging = entry["interp"].get("aging_risk", 0)
        raw += aging

        genotype_values = [
            g.get("aging_risk", 0) for g in entry["var"]["genotypes"].values()
        ]
        min_possible += min(genotype_values)
        max_possible += max(genotype_values)

    # Map raw from [min_possible, max_possible] to [0, 100]
    total_range = max_possible - min_possible
    if total_range <= 0:
        return 50

    score = int(round(((raw - min_possible) / total_range) * 100))
    return max(0, min(100, score))


def _determine_skin_type(results: List[Dict[str, Any]]) -> str:
    """Determine Fitzpatrick skin type tendency from MC1R and pigmentation variants."""
    mc1r_risk = 0
    for entry in results:
        if not entry["found"] or entry["interp"] is None:
            continue
        if entry["var"]["gene"].startswith("MC1R"):
            uv = entry["interp"].get("uv_risk", 0)
            if uv > 0:
                mc1r_risk += uv

    if mc1r_risk >= 4:
        return "Very Fair I-II"
    elif mc1r_risk >= 2:
        return "Fair II-III"
    elif mc1r_risk >= 1:
        return "Light-Medium III"
    else:
        return "Medium III+"


def _determine_uv_sensitivity(uv_risk_score: int) -> str:
    """Map UV risk score to a sensitivity label."""
    if uv_risk_score >= 75:
        return "Very High"
    elif uv_risk_score >= 50:
        return "High"
    elif uv_risk_score >= 25:
        return "Moderate"
    else:
        return "Standard"


# ─────────────────────────────────────────────────────────────────────────────
# JSON OUTPUT
# ─────────────────────────────────────────────────────────────────────────────

def generate_skincare_json(result: Dict[str, Any]) -> Dict[str, Any]:
    """
    Generate a JSON-serializable dictionary for the skincare report.

    Args:
        result: Output from analyze_skincare().

    Returns:
        Structured dict with summary, categories, and recommendations.
    """
    analysis_results = result["results"]

    uv_risk_score = _compute_uv_risk_score(analysis_results)
    aging_risk_score = _compute_aging_risk_score(analysis_results)
    skin_type = _determine_skin_type(analysis_results)
    uv_sensitivity = _determine_uv_sensitivity(uv_risk_score)

    # Build categories
    categories: List[Dict[str, Any]] = []
    for cat_name in SKINCARE_CATEGORIES:
        cat_variants = []
        for entry in analysis_results:
            if entry["var"]["category"] != cat_name:
                continue
            if not entry["found"]:
                continue

            interp = entry["interp"]
            cat_variants.append({
                "rsid": entry["rsid"],
                "gene": entry["var"]["gene"],
                "name": entry["var"]["name"],
                "genotype": entry["genotype"] or "",
                "label": interp["label"] if interp else "Unknown",
                "meaning": interp["meaning"] if interp else "Interpretation not available.",
                "evidence": entry["var"]["evidence"],
            })

        if cat_variants:
            categories.append({
                "name": cat_name,
                "emoji": CATEGORY_EMOJIS.get(cat_name, "📊"),
                "variants": cat_variants,
            })

    # Build recommendations
    recommendations: List[Dict[str, str]] = []
    for advice_key, advice_info in ADVICE.items():
        score = result["advice_scores"].get(advice_key, 0)
        level = _risk_level(score)
        recommendations.append({
            "title": advice_info["title"],
            "level": level,
            "text": advice_info[level],
        })

    return {
        "summary": {
            "totalChecked": result["total"],
            "found": result["found"],
            "notFound": result["not_found"],
            "skinType": skin_type,
            "uvSensitivity": uv_sensitivity,
            "uvRiskScore": uv_risk_score,
            "agingRiskScore": aging_risk_score,
        },
        "categories": categories,
        "recommendations": recommendations,
    }
