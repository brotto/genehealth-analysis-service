"""
Steroid Pharmacogenomics Analyzer
Analyzes genetic markers that determine how the body metabolizes,
absorbs, and responds to anabolic steroids and testosterone.

Based on peer-reviewed pharmacogenomic literature:
- PMC4760435, PMC2877023, PMC4382380 (UGT2B17)
- PMC4572434, PubMed 15377647 (AR CAG)
- PMC6722847, PubMed 38652149 (SHBG)
- NCBI NBK539904, PubMed 10898110 (SRD5A2)
- PMC12383265, PMC11060041 (ACTN3)
- PMC11625447 (CYP3A4/CYP3A5)
- PMC3098758, PMC5635765 (ABCB1)
"""

from typing import Dict, Tuple, Any, List, Optional
from datetime import datetime, timezone


# ─── SNP Database ────────────────────────────────────────────────────────────────

STEROID_SNP_DATABASE = {
    # ── ACTN3 — Muscle Fiber Type + Androgen Receptor Density ──
    "rs1815739": {
        "gene": "ACTN3",
        "full_name": "Alpha-Actinin-3 (R577X)",
        "category": "muscle_fiber",
        "chromosome": "11",
        "position": "66328095",
        "ref_allele": "C",
        "alt_allele": "T",
        "description": "Alpha-actinin-3 is a structural protein exclusive to fast-twitch (type II) muscle fibers. The R577X polymorphism determines whether you produce functional alpha-actinin-3.",
        "what_it_does": "Controls the presence of alpha-actinin-3 protein in fast-twitch muscle fibers and directly influences androgen receptor (AR) density in skeletal muscle.",
        "genotypes": {
            "CC": {
                "label": "RR — Full Alpha-Actinin-3",
                "status": "optimal",
                "fiber_type": "Predominantly fast-twitch (Type II)",
                "ar_density": "High",
                "summary": "You produce full alpha-actinin-3 protein. Your muscles have 65-72% MORE androgen receptors compared to XX genotype. This means your muscle tissue has significantly more \"docking stations\" for testosterone and other androgens — each molecule of hormone produces a stronger anabolic signal in your muscles.",
                "practical": [
                    "Your muscles are highly responsive to androgens — you get more anabolic effect per unit of hormone",
                    "You have a natural advantage for explosive strength, power, and muscle hypertrophy",
                    "Fast-twitch dominant fiber composition favors rapid force production",
                    "Testosterone and derivatives will produce strong anabolic responses in your muscle tissue"
                ],
                "response_rating": 5,
            },
            "CT": {
                "label": "RX — Partial Alpha-Actinin-3",
                "status": "moderate",
                "fiber_type": "Mixed fiber type (Type I + Type II)",
                "ar_density": "Moderate",
                "summary": "You have one working copy of ACTN3. Your androgen receptor density is intermediate — about 30-35% more than XX genotype. You have a balanced muscle fiber profile with good responsiveness to androgens.",
                "practical": [
                    "Good hormonal responsiveness — your muscles respond well to androgens, though not at the maximum level",
                    "Balanced fiber composition gives versatility for both strength and endurance",
                    "You can expect moderate-to-good anabolic response from testosterone",
                    "Training stimulus combined with androgens produces solid hypertrophy"
                ],
                "response_rating": 3,
            },
            "TT": {
                "label": "XX — No Alpha-Actinin-3",
                "status": "reduced",
                "fiber_type": "Predominantly slow-twitch (Type I)",
                "ar_density": "Low",
                "summary": "You do not produce alpha-actinin-3. Your muscles have 65-72% FEWER androgen receptors compared to RR genotype. This means your muscle tissue has significantly fewer \"docking stations\" for testosterone — you need more hormone to achieve the same anabolic signal.",
                "practical": [
                    "Your muscles are less responsive to androgens — you need higher concentrations to achieve the same effect",
                    "Endurance-oriented fiber profile with lower explosive power capacity",
                    "Testosterone and derivatives will produce a weaker anabolic response per unit of hormone",
                    "Higher volume training may partially compensate for lower receptor density"
                ],
                "response_rating": 1,
            },
        },
        "references": ["PMC12383265", "PMC11060041"],
        "population_note": "RR is found in ~30% of Europeans, ~52% of East Africans. XX is ~18% in Europeans.",
    },

    # ── CYP3A4*22 — Hepatic Metabolism (Major) ──
    "rs35599367": {
        "gene": "CYP3A4",
        "full_name": "CYP3A4*22 — Hepatic Metabolism",
        "category": "hepatic_metabolism",
        "chromosome": "7",
        "position": "99366316",
        "ref_allele": "G",  # *1 (normal)
        "alt_allele": "A",  # *22 (reduced)
        "description": "CYP3A4 is the primary liver enzyme responsible for metabolizing oral anabolic steroids including Stanozolol (Winstrol), Oxandrolone (Anavar), Methandienone (Dianabol), and Turinabol.",
        "what_it_does": "Determines how quickly your liver breaks down oral steroids. This directly affects how much of each oral dose actually reaches your bloodstream (bioavailability).",
        "genotypes": {
            "GG": {
                "label": "*1/*1 — Normal Metabolizer",
                "status": "normal",
                "metabolism_speed": "Normal",
                "summary": "Your CYP3A4 enzyme works at standard speed. Oral steroids are metabolized at a normal rate through your liver. You get the expected bioavailability from standard oral doses.",
                "practical": [
                    "Standard oral steroid metabolism — no dose adjustment needed for liver processing speed",
                    "Normal hepatic clearance rate for Stanozolol, Oxandrolone, Dianabol, Turinabol",
                    "Standard liver stress from oral 17-alpha-alkylated compounds",
                    "Blood levels follow typical pharmacokinetic curves for oral compounds"
                ],
                "response_rating": 3,
            },
            "GA": {
                "label": "*1/*22 — Slow Metabolizer",
                "status": "favorable",
                "metabolism_speed": "Slow (~30-50% reduced activity)",
                "summary": "You carry one copy of the CYP3A4*22 allele, which reduces enzyme activity by approximately 30-50%. Oral steroids stay in your system longer and reach higher blood concentrations. You get MORE effect per milligram of oral compound.",
                "practical": [
                    "You have HIGHER bioavailability of oral steroids — each mg produces more effect than average",
                    "Oral compounds like Oxandrolone and Stanozolol will circulate longer in your bloodstream",
                    "Lower doses may achieve the same blood levels as standard doses in normal metabolizers",
                    "Increased liver exposure time — be mindful of hepatotoxicity with 17-alpha-alkylated orals"
                ],
                "response_rating": 4,
            },
            "AA": {
                "label": "*22/*22 — Very Slow Metabolizer",
                "status": "favorable_with_caution",
                "metabolism_speed": "Very Slow (~50-70% reduced activity)",
                "summary": "You carry two copies of CYP3A4*22. Your liver metabolizes oral steroids very slowly — compounds remain active in your bloodstream significantly longer. You get substantially MORE effect per milligram, but also increased hepatic stress.",
                "practical": [
                    "You have SIGNIFICANTLY higher bioavailability of oral steroids — strong effect per dose",
                    "Oral compounds persist much longer in your system than in most people",
                    "Much lower doses needed to achieve therapeutic/effective blood levels",
                    "CAUTION: Higher hepatotoxicity risk — liver is exposed to active compounds for longer periods"
                ],
                "response_rating": 5,
            },
        },
        "references": ["PMC11625447"],
        "population_note": "CYP3A4*22 carriers are ~5-7% of European populations. Homozygous *22/*22 is very rare (<1%).",
    },

    # ── CYP3A4*1B — Hepatic Metabolism (Secondary) ──
    "rs2740574": {
        "gene": "CYP3A4",
        "full_name": "CYP3A4*1B — Hepatic Expression",
        "category": "hepatic_metabolism",
        "chromosome": "7",
        "position": "99382096",
        "ref_allele": "T",  # *1A (normal)
        "alt_allele": "C",  # *1B (modest reduction)
        "description": "CYP3A4*1B is a promoter variant that modestly reduces CYP3A4 expression. It contributes to overall hepatic metabolism capacity alongside CYP3A4*22.",
        "what_it_does": "Modestly reduces the amount of CYP3A4 enzyme your liver produces, complementing the effect of other CYP3A4 variants.",
        "genotypes": {
            "TT": {
                "label": "*1A/*1A — Normal Expression",
                "status": "normal",
                "metabolism_speed": "Normal",
                "summary": "Standard CYP3A4 enzyme expression. This allele does not alter your oral steroid metabolism.",
                "practical": [
                    "Normal CYP3A4 expression — no additional effect on oral steroid metabolism"
                ],
                "response_rating": 3,
            },
            "TC": {
                "label": "*1A/*1B — Slightly Reduced Expression",
                "status": "slightly_favorable",
                "metabolism_speed": "Slightly reduced",
                "summary": "One copy of *1B slightly reduces CYP3A4 production. Minor increase in oral steroid bioavailability.",
                "practical": [
                    "Marginally slower metabolism of oral steroids",
                    "Small increase in bioavailability — effect is modest but additive with other CYP3A4 variants"
                ],
                "response_rating": 3,
            },
            "CC": {
                "label": "*1B/*1B — Reduced Expression",
                "status": "favorable",
                "metabolism_speed": "Moderately reduced",
                "summary": "Two copies of *1B reduce CYP3A4 expression. Noticeable increase in oral steroid bioavailability when combined with other metabolic variants.",
                "practical": [
                    "Moderately slower oral steroid metabolism",
                    "Better bioavailability from oral compounds",
                    "Effect is additive with CYP3A4*22 if both are present"
                ],
                "response_rating": 4,
            },
        },
        "references": ["PMC11625447"],
        "population_note": "*1B is found in ~35-67% of African populations and ~2-9% of Europeans.",
    },

    # ── CYP3A5*3 — Hepatic Metabolism (Auxiliary) ──
    "rs776746": {
        "gene": "CYP3A5",
        "full_name": "CYP3A5*3 — Auxiliary Hepatic Metabolism",
        "category": "hepatic_metabolism",
        "chromosome": "7",
        "position": "99270539",
        "ref_allele": "C",  # *3 (loss of function — most common in Europeans)
        "alt_allele": "T",  # *1 (functional)
        "description": "CYP3A5 is a secondary enzyme that assists CYP3A4 in metabolizing oral steroids. Most Europeans lack functional CYP3A5 due to the *3 allele.",
        "what_it_does": "Determines whether you have an additional metabolic pathway for breaking down oral steroids alongside CYP3A4.",
        "genotypes": {
            "CC": {
                "label": "*3/*3 — Non-Expressor (No CYP3A5)",
                "status": "favorable",
                "metabolism_speed": "One less metabolic pathway",
                "summary": "You do NOT produce functional CYP3A5 enzyme. You rely solely on CYP3A4 for oral steroid metabolism. This means oral compounds are cleared more slowly since you have one fewer metabolic pathway.",
                "practical": [
                    "Only CYP3A4 handles oral steroid metabolism — slower overall clearance",
                    "Higher effective bioavailability of oral steroids",
                    "This is the most common genotype in Europeans (~85-95%)"
                ],
                "response_rating": 4,
            },
            "CT": {
                "label": "*1/*3 — Partial Expressor",
                "status": "moderate",
                "metabolism_speed": "Additional pathway present",
                "summary": "You produce some CYP3A5 enzyme, providing an additional metabolic pathway for oral steroid clearance. Oral compounds are cleared somewhat faster than in non-expressors.",
                "practical": [
                    "You have an extra metabolic pathway that clears oral steroids faster",
                    "Slightly reduced bioavailability of oral compounds compared to non-expressors",
                    "May need slightly higher oral doses to achieve the same blood levels"
                ],
                "response_rating": 3,
            },
            "TT": {
                "label": "*1/*1 — Full Expressor",
                "status": "reduced",
                "metabolism_speed": "Full additional pathway",
                "summary": "You produce full CYP3A5 enzyme activity. Both CYP3A4 and CYP3A5 actively metabolize oral steroids, leading to faster clearance and lower bioavailability per dose.",
                "practical": [
                    "Oral steroids are metabolized faster through two active pathways (CYP3A4 + CYP3A5)",
                    "Lower bioavailability per oral dose — compounds are cleared more quickly",
                    "May need higher or more frequent oral doses to maintain blood levels",
                    "Benefit: faster liver recovery since compounds are cleared quicker"
                ],
                "response_rating": 2,
            },
        },
        "references": ["PMC11625447"],
        "population_note": "*3/*3 (non-expressor) is ~85-95% of Europeans but only ~15-30% of Africans.",
    },

    # ── ABCB1/MDR1 — Intestinal Absorption (Oral Compounds) ──
    "rs1045642": {
        "gene": "ABCB1",
        "full_name": "ABCB1/MDR1 C3435T — P-glycoprotein (Intestinal Absorption)",
        "category": "absorption",
        "chromosome": "7",
        "position": "87138645",
        "ref_allele": "A",  # Note: strand-dependent, A=T allele in some references
        "alt_allele": "G",  # C allele
        "description": "ABCB1 encodes P-glycoprotein (P-gp), an efflux pump in the intestinal wall that actively pushes ingested compounds back into the gut lumen, reducing oral absorption.",
        "what_it_does": "Controls how much of an oral steroid actually gets absorbed through your intestines into your bloodstream. More P-gp = less absorption. Only affects oral compounds — injectable steroids bypass this entirely.",
        "genotypes": {
            "AA": {
                "label": "TT — Low P-gp (High Absorption)",
                "status": "favorable",
                "absorption_level": "High",
                "summary": "You produce LESS P-glycoprotein in your intestinal wall. Oral steroids pass through your gut lining more efficiently. You absorb a higher percentage of each oral dose into your bloodstream.",
                "practical": [
                    "You have GOOD oral absorption — a higher percentage of each oral dose reaches your bloodstream",
                    "Oral compounds like Oxandrolone, Stanozolol, Dianabol, Turinabol are well absorbed",
                    "You get more effect per milligram of oral compound compared to CC genotype",
                    "This only affects ORAL compounds — injectables (enanthate, cypionate, propionate) bypass the gut entirely"
                ],
                "response_rating": 5,
            },
            "AG": {
                "label": "CT — Moderate P-gp (Moderate Absorption)",
                "status": "moderate",
                "absorption_level": "Moderate",
                "summary": "You have intermediate P-glycoprotein levels. Oral steroid absorption is average — you get a standard percentage of each oral dose into your bloodstream.",
                "practical": [
                    "Average oral absorption — standard bioavailability from oral compounds",
                    "No specific advantage or disadvantage for oral steroid absorption",
                    "Standard doses should produce expected blood levels"
                ],
                "response_rating": 3,
            },
            "GG": {
                "label": "CC — High P-gp (Reduced Absorption)",
                "status": "reduced",
                "absorption_level": "Low",
                "summary": "You produce MORE P-glycoprotein, which actively pumps oral steroids back out of your intestinal cells. A significant portion of each oral dose never reaches your bloodstream.",
                "practical": [
                    "You have REDUCED oral absorption — the P-gp pump pushes oral compounds back into your gut",
                    "You may need higher oral doses to achieve the same blood levels as TT genotype",
                    "Consider that injectable compounds completely bypass this limitation",
                    "Oral-only cycles may be less effective for you compared to injectable alternatives"
                ],
                "response_rating": 1,
            },
        },
        "references": ["PMC3098758", "PMC5635765"],
        "population_note": "TT is ~25% of Europeans, CT ~50%, CC ~25%.",
    },

    # ── ABCB1/MDR1 G2677T — Intestinal Absorption (Secondary) ──
    "rs2032582": {
        "gene": "ABCB1",
        "full_name": "ABCB1/MDR1 G2677T — P-glycoprotein (Secondary Marker)",
        "category": "absorption",
        "chromosome": "7",
        "position": "87160618",
        "ref_allele": "A",
        "alt_allele": "C",
        "description": "A secondary ABCB1 marker that further modulates P-glycoprotein activity. Works together with rs1045642 to determine overall intestinal absorption capacity.",
        "what_it_does": "Fine-tunes P-glycoprotein function alongside C3435T. The combined haplotype of both SNPs determines your overall oral absorption efficiency.",
        "genotypes": {
            "AA": {
                "label": "TT — Lower P-gp Activity",
                "status": "favorable",
                "absorption_level": "Higher absorption",
                "summary": "This variant reduces P-glycoprotein activity, improving oral compound absorption. Combined with rs1045642 TT, you have excellent oral absorption.",
                "practical": [
                    "Additional reduction in P-gp activity — supports better oral absorption",
                    "Combined with C3435T results, gives a complete picture of your oral absorption profile"
                ],
                "response_rating": 4,
            },
            "AC": {
                "label": "GT — Intermediate",
                "status": "moderate",
                "absorption_level": "Average",
                "summary": "Intermediate P-glycoprotein activity from this marker.",
                "practical": [
                    "Standard contribution to oral absorption — no major effect alone"
                ],
                "response_rating": 3,
            },
            "CC": {
                "label": "GG — Higher P-gp Activity",
                "status": "reduced",
                "absorption_level": "Reduced absorption",
                "summary": "This variant supports higher P-glycoprotein activity, contributing to reduced oral absorption.",
                "practical": [
                    "Contributes to higher P-gp activity and lower oral absorption",
                    "Combined with rs1045642 CC, suggests consistently reduced oral bioavailability"
                ],
                "response_rating": 2,
            },
        },
        "references": ["PMC3098758", "PMC5635765"],
        "population_note": "Allele frequencies vary widely by population.",
    },

    # ── SHBG — Sex Hormone-Binding Globulin ──
    "rs1799941": {
        "gene": "SHBG",
        "full_name": "SHBG Promoter Variant — Free Testosterone Availability",
        "category": "hormone_binding",
        "chromosome": "17",
        "position": "7533423",
        "ref_allele": "G",
        "alt_allele": "A",
        "description": "SHBG (Sex Hormone-Binding Globulin) is a protein that binds to testosterone in your blood, making it inactive. Only UNBOUND (free) testosterone can enter muscle cells and activate androgen receptors.",
        "what_it_does": "Determines how much of your total testosterone is actually FREE and available to build muscle. Higher SHBG = more testosterone locked up and unavailable. Lower SHBG = more testosterone actively working in your muscles.",
        "genotypes": {
            "GG": {
                "label": "GG — Low SHBG (More Free Testosterone)",
                "status": "favorable",
                "shbg_level": "Lower",
                "summary": "You have the genotype associated with LOWER SHBG production. This means a larger fraction of your total testosterone circulates in its FREE, bioactive form. More testosterone is available to enter muscle cells and activate androgen receptors.",
                "practical": [
                    "You have MORE free testosterone available — more hormone reaching your muscles",
                    "Better anabolic response per unit of total testosterone (endogenous or exogenous)",
                    "Your testosterone-to-SHBG ratio is genetically favorable for muscle building",
                    "Both natural and exogenous testosterone will have higher effective bioavailability"
                ],
                "response_rating": 5,
            },
            "GA": {
                "label": "GA — Moderate SHBG",
                "status": "moderate",
                "shbg_level": "Average",
                "summary": "You have one copy of the high-SHBG allele. Your SHBG levels are average, with a standard proportion of free vs. bound testosterone.",
                "practical": [
                    "Average free testosterone availability — standard hormonal bioavailability",
                    "No specific advantage or disadvantage for testosterone utilization"
                ],
                "response_rating": 3,
            },
            "AA": {
                "label": "AA — High SHBG (Less Free Testosterone)",
                "status": "reduced",
                "shbg_level": "Higher (+26.9% SHBG)",
                "summary": "You carry two copies of the high-SHBG allele, associated with ~26.9% higher SHBG levels. More of your testosterone is BOUND and inactive. Less free testosterone reaches your muscle cells.",
                "practical": [
                    "You have LESS free testosterone — more hormone is locked up by SHBG protein",
                    "You may need higher total testosterone levels to achieve the same free testosterone as GG genotype",
                    "Endogenous and exogenous testosterone are both affected — more gets bound before reaching muscles",
                    "Consider that some compounds (e.g., mesterolone) can lower SHBG and increase free T"
                ],
                "response_rating": 1,
            },
        },
        "references": ["PMC6722847", "PubMed 38652149", "PubMed 24327369"],
        "population_note": "Allele A frequency is ~35-40% in European populations.",
    },

    # ── SRD5A2 V89L — 5-Alpha Reductase Type 2 ──
    "rs523349": {
        "gene": "SRD5A2",
        "full_name": "SRD5A2 V89L — 5-Alpha Reductase Type 2 (Testosterone → DHT Conversion)",
        "category": "dht_conversion",
        "chromosome": "2",
        "position": "31805706",
        "ref_allele": "C",  # V89 (Valine — wild type)
        "alt_allele": "G",  # L89 (Leucine — reduced activity)
        "description": "SRD5A2 converts testosterone into DHT (Dihydrotestosterone), which is 3-5x more potent at the androgen receptor. DHT drives androgenic effects: body/facial hair, prostate stimulation, sebaceous gland activity, and some muscle effects.",
        "what_it_does": "Controls the rate of conversion of testosterone into its more potent form (DHT). Higher conversion = more androgenic effects (both desirable and undesirable).",
        "genotypes": {
            "CC": {
                "label": "VV — Full 5α-Reductase Activity",
                "status": "high_conversion",
                "conversion_rate": "Normal/High",
                "summary": "You have full SRD5A2 enzyme activity. Testosterone is efficiently converted to DHT at the standard rate. You experience the full spectrum of androgenic effects from both endogenous and exogenous testosterone.",
                "practical": [
                    "Full testosterone → DHT conversion — maximum androgenic potency in target tissues",
                    "Higher DHT means stronger effects on: facial/body hair, prostate, scalp (hair loss), skin oil",
                    "Compounds like testosterone will produce significant DHT-mediated side effects",
                    "Nandrolone (Deca) converts to DHN instead (weaker), so this enzyme matters less for nandrolone"
                ],
                "response_rating": 4,
            },
            "CG": {
                "label": "VL — Moderately Reduced Activity",
                "status": "moderate_conversion",
                "conversion_rate": "Moderately Reduced (~15%)",
                "summary": "One copy of the L89 variant reduces your 5α-reductase activity by approximately 15%. You convert somewhat less testosterone to DHT.",
                "practical": [
                    "Slightly less DHT production — moderately reduced androgenic side effects",
                    "Some protection against DHT-driven hair loss and prostate stimulation",
                    "Still significant DHT conversion occurs — this is a moderate reduction, not elimination"
                ],
                "response_rating": 3,
            },
            "GG": {
                "label": "LL — Significantly Reduced Activity (~30%)",
                "status": "low_conversion",
                "conversion_rate": "Reduced by ~30%",
                "summary": "Two copies of L89 reduce your 5α-reductase activity by approximately 30%. You produce significantly less DHT from testosterone. This affects both androgenic side effects and DHT-specific tissue responses.",
                "practical": [
                    "Notably less DHT production — reduced risk of DHT-driven hair loss and acne",
                    "Less prostate stimulation from testosterone-based compounds",
                    "The purely anabolic (muscle-building) effects of testosterone are less affected",
                    "You may experience fewer androgenic side effects but also less facial/body hair growth"
                ],
                "response_rating": 2,
            },
        },
        "references": ["NCBI NBK539904", "PubMed 10898110"],
        "population_note": "L89 allele frequency varies: ~37% in East Asian, ~17% in European, ~10% in African populations.",
    },

    # ── SRD5A2 A49T — 5-Alpha Reductase (Secondary) ──
    "rs9282858": {
        "gene": "SRD5A2",
        "full_name": "SRD5A2 A49T — Enhanced 5α-Reductase Activity",
        "category": "dht_conversion",
        "chromosome": "2",
        "position": "31805826",
        "ref_allele": "C",  # A49 (normal)
        "alt_allele": "T",  # T49 (increased activity)
        "description": "The A49T variant of SRD5A2 INCREASES enzyme activity, leading to higher DHT conversion. This is the opposite direction from V89L.",
        "what_it_does": "When present, this variant increases the conversion of testosterone to DHT beyond normal levels.",
        "genotypes": {
            "CC": {
                "label": "AA — Normal Activity",
                "status": "normal",
                "conversion_rate": "Normal",
                "summary": "Standard 5α-reductase activity at this position. No enhanced DHT conversion from this variant.",
                "practical": [
                    "No additional DHT conversion effect from this variant"
                ],
                "response_rating": 3,
            },
            "CT": {
                "label": "AT — Increased Activity",
                "status": "enhanced",
                "conversion_rate": "Increased",
                "summary": "One copy of T49 increases 5α-reductase activity. More testosterone is converted to DHT than average.",
                "practical": [
                    "Higher DHT conversion — increased androgenic potency from testosterone-based compounds",
                    "Greater risk of DHT-related side effects: hair loss, acne, prostate stimulation",
                    "Combined with rs523349, gives a complete picture of your DHT conversion profile"
                ],
                "response_rating": 4,
            },
            "TT": {
                "label": "TT — Significantly Increased Activity",
                "status": "highly_enhanced",
                "conversion_rate": "Significantly Increased",
                "summary": "Two copies of T49 significantly increase 5α-reductase activity. You convert substantially more testosterone to DHT than most people.",
                "practical": [
                    "Significantly higher DHT production from testosterone — strong androgenic effects",
                    "Elevated risk of hair loss, acne, and prostate stimulation on testosterone",
                    "Consider DHT-sparing compounds (Nandrolone, Boldenone) if androgenic sides are a concern"
                ],
                "response_rating": 5,
            },
        },
        "references": ["PubMed 10898110"],
        "population_note": "T49 allele is rare in most populations (<5%).",
    },

    # ── UGT2B17 Proxy — Elimination Speed ──
    "rs7436962": {
        "gene": "UGT2B17",
        "full_name": "UGT2B17 Proxy — Testosterone Elimination Speed",
        "category": "elimination",
        "chromosome": "4",
        "position": "69400000",
        "ref_allele": "G",
        "alt_allele": "A",
        "description": "UGT2B17 is the primary enzyme responsible for glucuronidation of testosterone — the process that makes testosterone water-soluble so it can be excreted by the kidneys. This SNP is a proxy for UGT2B17 gene deletion.",
        "what_it_does": "Controls how quickly your body eliminates testosterone through the kidneys. Slower elimination = testosterone stays active longer. Faster elimination = testosterone is cleared quickly.",
        "genotypes": {
            "GG": {
                "label": "Full UGT2B17 (Fast Elimination)",
                "status": "fast_elimination",
                "elimination_speed": "Fast",
                "summary": "You likely have two functional copies of UGT2B17. Testosterone is glucuronidated and excreted efficiently. Blood levels drop faster after administration.",
                "practical": [
                    "Testosterone is cleared from your system relatively quickly",
                    "You may need more frequent dosing or higher doses to maintain stable blood levels",
                    "Benefit: lower risk of metabolite accumulation in kidneys"
                ],
                "response_rating": 2,
            },
            "GA": {
                "label": "Partial UGT2B17 (Moderate Elimination)",
                "status": "moderate_elimination",
                "elimination_speed": "Moderate",
                "summary": "You likely have one functional copy of UGT2B17. Testosterone elimination is at an intermediate rate.",
                "practical": [
                    "Intermediate testosterone clearance — balanced elimination speed",
                    "Standard dosing intervals should work well"
                ],
                "response_rating": 3,
            },
            "AA": {
                "label": "UGT2B17 Deletion (Slow Elimination)",
                "status": "slow_elimination",
                "elimination_speed": "Slow",
                "summary": "You likely have reduced or absent UGT2B17 activity. Testosterone circulates in your body much longer before being cleared. Each dose has a longer effective duration.",
                "practical": [
                    "Testosterone stays active in your bloodstream significantly longer",
                    "You may need less frequent dosing or lower doses to maintain target blood levels",
                    "CAUTION: increased risk of renal stress with prolonged use — unconjugated metabolites accumulate",
                    "This is common in East Asian populations (~80% deletion)"
                ],
                "response_rating": 4,
            },
        },
        "references": ["PMC4760435", "PMC2877023", "PMC4382380"],
        "population_note": "Full deletion (del/del) is ~12% in Europeans and ~80% in East Asians.",
    },
}


# ─── Analysis Functions ──────────────────────────────────────────────────────────

def analyze_steroid_pharmacogenomics(variants: Dict[str, Tuple[str, str, str]]) -> Dict[str, Any]:
    """
    Analyze genetic markers related to steroid metabolism and androgen response.

    Args:
        variants: Dictionary mapping rsID to (chromosome, position, genotype)

    Returns:
        Complete analysis result dictionary
    """
    markers_found = []
    markers_not_found = []

    for rsid, snp_info in STEROID_SNP_DATABASE.items():
        rsid_lower = rsid.lower()

        if rsid_lower in variants:
            chrom, pos, genotype = variants[rsid_lower]
            genotype_upper = genotype.upper().replace('-', '').replace(' ', '')

            # Sort genotype for consistent lookup (e.g., AG -> AG, GA -> AG)
            if len(genotype_upper) == 2:
                sorted_genotype = ''.join(sorted(genotype_upper))
            else:
                sorted_genotype = genotype_upper

            # Try to find matching genotype info
            geno_info = snp_info["genotypes"].get(genotype_upper)
            if not geno_info:
                geno_info = snp_info["genotypes"].get(sorted_genotype)
            if not geno_info:
                # Try reverse
                if len(genotype_upper) == 2:
                    reverse = genotype_upper[1] + genotype_upper[0]
                    geno_info = snp_info["genotypes"].get(reverse)

            if geno_info:
                markers_found.append({
                    "rsid": rsid,
                    "gene": snp_info["gene"],
                    "full_name": snp_info["full_name"],
                    "category": snp_info["category"],
                    "genotype": genotype_upper,
                    "description": snp_info["description"],
                    "what_it_does": snp_info["what_it_does"],
                    "label": geno_info["label"],
                    "status": geno_info["status"],
                    "summary": geno_info["summary"],
                    "practical": geno_info["practical"],
                    "response_rating": geno_info["response_rating"],
                    "references": snp_info.get("references", []),
                    "population_note": snp_info.get("population_note", ""),
                    "extra": {k: v for k, v in geno_info.items() if k not in ("label", "status", "summary", "practical", "response_rating")},
                })
            else:
                markers_found.append({
                    "rsid": rsid,
                    "gene": snp_info["gene"],
                    "full_name": snp_info["full_name"],
                    "category": snp_info["category"],
                    "genotype": genotype_upper,
                    "description": snp_info["description"],
                    "what_it_does": snp_info["what_it_does"],
                    "label": f"Genotype {genotype_upper}",
                    "status": "indeterminate",
                    "summary": f"Your genotype ({genotype_upper}) was detected but does not match the studied variants for this marker. This may be due to a rare allele or strand orientation difference.",
                    "practical": ["Genotype detected but classification is unavailable for this specific combination"],
                    "response_rating": 3,
                    "references": snp_info.get("references", []),
                    "population_note": snp_info.get("population_note", ""),
                    "extra": {},
                })
        else:
            markers_not_found.append({
                "rsid": rsid,
                "gene": snp_info["gene"],
                "full_name": snp_info["full_name"],
                "category": snp_info["category"],
            })

    # Compute composite scores
    composite = _compute_composite_scores(markers_found)

    # Generate category summaries
    categories = _generate_category_analysis(markers_found, markers_not_found)

    # Generate practical profile
    profile = _generate_practical_profile(markers_found, composite)

    return {
        "markers_found": markers_found,
        "markers_not_found": markers_not_found,
        "composite_scores": composite,
        "categories": categories,
        "practical_profile": profile,
    }


def _compute_composite_scores(markers: List[Dict]) -> Dict[str, Any]:
    """Compute composite sensitivity and response scores."""
    if not markers:
        return {
            "overall_androgen_sensitivity": 50,
            "oral_bioavailability": 50,
            "dht_conversion_tendency": 50,
            "muscle_receptor_density": 50,
            "free_testosterone_availability": 50,
            "elimination_speed": 50,
            "overall_label": "Insufficient Data",
        }

    category_scores = {}
    for m in markers:
        cat = m["category"]
        rating = m["response_rating"]
        if cat not in category_scores:
            category_scores[cat] = []
        category_scores[cat].append(rating)

    def avg_score(cat: str, default: float = 3.0) -> float:
        vals = category_scores.get(cat, [default])
        return sum(vals) / len(vals)

    def to_100(score: float) -> int:
        return max(0, min(100, int((score / 5.0) * 100)))

    muscle_score = avg_score("muscle_fiber")
    hepatic_score = avg_score("hepatic_metabolism")
    absorption_score = avg_score("absorption")
    hormone_score = avg_score("hormone_binding")
    dht_score = avg_score("dht_conversion")
    elimination_score = avg_score("elimination")

    # Overall = weighted average
    weights = {
        "muscle_fiber": 0.25,
        "hepatic_metabolism": 0.20,
        "absorption": 0.15,
        "hormone_binding": 0.20,
        "dht_conversion": 0.10,
        "elimination": 0.10,
    }

    total_weight = 0
    weighted_sum = 0
    for cat, w in weights.items():
        if cat in category_scores:
            weighted_sum += avg_score(cat) * w
            total_weight += w

    overall = weighted_sum / total_weight if total_weight > 0 else 3.0
    overall_100 = to_100(overall)

    if overall_100 >= 75:
        overall_label = "High Androgen Sensitivity"
    elif overall_100 >= 55:
        overall_label = "Moderate-High Sensitivity"
    elif overall_100 >= 40:
        overall_label = "Moderate Sensitivity"
    elif overall_100 >= 25:
        overall_label = "Moderate-Low Sensitivity"
    else:
        overall_label = "Low Androgen Sensitivity"

    return {
        "overall_androgen_sensitivity": overall_100,
        "oral_bioavailability": to_100((hepatic_score + absorption_score) / 2 if (category_scores.get("hepatic_metabolism") or category_scores.get("absorption")) else 3.0),
        "dht_conversion_tendency": to_100(dht_score),
        "muscle_receptor_density": to_100(muscle_score),
        "free_testosterone_availability": to_100(hormone_score),
        "elimination_speed": to_100(elimination_score),
        "overall_label": overall_label,
    }


def _generate_category_analysis(found: List[Dict], not_found: List[Dict]) -> List[Dict]:
    """Group markers into functional categories with summaries."""
    CATEGORY_META = {
        "muscle_fiber": {
            "name": "Muscle Fiber Composition & AR Density",
            "icon": "💪",
            "description": "How your muscle fibers are composed and how many androgen receptors they contain",
        },
        "hepatic_metabolism": {
            "name": "Liver Metabolism (Oral Compounds)",
            "icon": "🧬",
            "description": "How quickly your liver breaks down oral steroids — affects bioavailability of oral compounds",
        },
        "absorption": {
            "name": "Intestinal Absorption (Oral Compounds)",
            "icon": "💊",
            "description": "How efficiently your gut absorbs oral steroids — only relevant for oral administration",
        },
        "hormone_binding": {
            "name": "Hormone Binding & Free Testosterone",
            "icon": "🔗",
            "description": "How much of your testosterone is free (active) vs. bound to SHBG (inactive)",
        },
        "dht_conversion": {
            "name": "DHT Conversion (Testosterone → DHT)",
            "icon": "⚡",
            "description": "How much testosterone gets converted to DHT — affects androgenic side effects and potency",
        },
        "elimination": {
            "name": "Elimination & Clearance",
            "icon": "🔄",
            "description": "How quickly your body eliminates testosterone through the kidneys",
        },
    }

    categories = []
    for cat_key, meta in CATEGORY_META.items():
        cat_found = [m for m in found if m["category"] == cat_key]
        cat_missing = [m for m in not_found if m["category"] == cat_key]

        if not cat_found and not cat_missing:
            continue

        avg_rating = sum(m["response_rating"] for m in cat_found) / len(cat_found) if cat_found else None

        categories.append({
            "key": cat_key,
            "name": meta["name"],
            "icon": meta["icon"],
            "description": meta["description"],
            "markers_found": cat_found,
            "markers_not_found": cat_missing,
            "average_rating": avg_rating,
            "total_markers": len(cat_found) + len(cat_missing),
            "found_count": len(cat_found),
        })

    return categories


def _generate_practical_profile(markers: List[Dict], composite: Dict) -> Dict[str, Any]:
    """Generate a user-friendly practical profile based on all findings."""
    strengths = []
    considerations = []
    warnings = []

    for m in markers:
        rating = m["response_rating"]
        gene = m["gene"]
        label = m["label"]

        if rating >= 4:
            if m["category"] == "muscle_fiber":
                strengths.append(f"Your {gene} genotype ({label.split('—')[0].strip()}) gives you high androgen receptor density — your muscles are highly responsive to hormones")
            elif m["category"] == "hepatic_metabolism":
                strengths.append(f"Your {gene} genotype means slower liver metabolism — oral steroids have higher bioavailability in your body")
            elif m["category"] == "absorption":
                strengths.append(f"Your {gene} genotype ({label.split('—')[0].strip()}) means good intestinal absorption of oral compounds")
            elif m["category"] == "hormone_binding":
                strengths.append(f"Your {gene} genotype favors lower SHBG — more of your testosterone circulates in its free, active form")
            elif m["category"] == "elimination":
                strengths.append(f"Your {gene} genotype means slower testosterone elimination — each dose has a longer effective duration")

        elif rating <= 2:
            if m["category"] == "muscle_fiber":
                considerations.append(f"Your {gene} genotype ({label.split('—')[0].strip()}) means fewer androgen receptors in muscle — you may need higher hormone levels for the same anabolic effect")
            elif m["category"] == "hepatic_metabolism":
                considerations.append(f"Your {gene} enzyme is highly active — oral steroids are cleared quickly, reducing their effective duration")
            elif m["category"] == "absorption":
                considerations.append(f"Your {gene} genotype means higher P-glycoprotein — oral compounds are less efficiently absorbed")
            elif m["category"] == "hormone_binding":
                considerations.append(f"Your {gene} genotype favors higher SHBG — more testosterone is bound and inactive")
            elif m["category"] == "elimination":
                considerations.append(f"Your {gene} genotype means faster testosterone clearance — you may need more frequent dosing")

        # Check for caution combinations
        if m.get("status") == "favorable_with_caution":
            warnings.append(f"{gene}: {m['summary']}")

    # Check for slow elimination + oral caution
    has_slow_metabolism = any(m["category"] == "hepatic_metabolism" and m["response_rating"] >= 4 for m in markers)
    has_slow_elimination = any(m["category"] == "elimination" and m["response_rating"] >= 4 for m in markers)
    if has_slow_metabolism and has_slow_elimination:
        warnings.append("You have both slow hepatic metabolism AND slow elimination — oral compounds will have significantly extended duration and higher accumulation. Monitor liver and kidney markers closely.")

    overall = composite["overall_androgen_sensitivity"]
    if overall >= 70:
        overall_summary = "Your genetic profile suggests HIGH overall androgen sensitivity. Your body is well-equipped to respond to both endogenous and exogenous androgens. You likely experience strong anabolic responses and may need lower doses to achieve target effects."
    elif overall >= 50:
        overall_summary = "Your genetic profile suggests MODERATE androgen sensitivity. You have a balanced response profile with some favorable and some neutral genetic factors. Standard approaches should work well for you."
    else:
        overall_summary = "Your genetic profile suggests LOWER androgen sensitivity. You may need higher hormone levels or longer exposure times to achieve the same effects as someone with higher genetic sensitivity. Consider optimizing training, nutrition, and recovery to maximize your response."

    return {
        "overall_summary": overall_summary,
        "strengths": strengths,
        "considerations": considerations,
        "warnings": warnings,
    }


# ─── JSON Report Generator ──────────────────────────────────────────────────────

def generate_steroid_pharmacogenomics_json(result: Dict[str, Any]) -> Dict[str, Any]:
    """Generate the JSON report for frontend consumption."""
    markers_found = result["markers_found"]
    markers_not_found = result["markers_not_found"]
    composite = result["composite_scores"]
    categories = result["categories"]
    profile = result["practical_profile"]

    return {
        "report_type": "steroid_pharmacogenomics",
        "version": "1.0",
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "disclaimer": "IMPORTANT DISCLAIMER: This report presents pharmacogenomic information based on published scientific studies. It does NOT constitute medical advice, prescription, or encouragement to use controlled substances. The use of anabolic steroids without a medical prescription is illegal in many countries and may cause serious health risks. Always consult an endocrinologist or qualified physician before making any decisions related to hormonal therapy. This report is for educational and informational purposes only.",
        "summary": {
            "total_markers_analyzed": len(STEROID_SNP_DATABASE),
            "markers_found": len(markers_found),
            "markers_not_found": len(markers_not_found),
            "overall_score": composite["overall_androgen_sensitivity"],
            "overall_label": composite["overall_label"],
            "oral_bioavailability_score": composite["oral_bioavailability"],
            "dht_conversion_score": composite["dht_conversion_tendency"],
            "muscle_receptor_score": composite["muscle_receptor_density"],
            "free_testosterone_score": composite["free_testosterone_availability"],
            "elimination_score": composite["elimination_speed"],
        },
        "categories": categories,
        "practical_profile": profile,
        "composite_scores": composite,
    }
