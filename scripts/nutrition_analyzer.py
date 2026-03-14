"""
Nutrition Analyzer (Nutrigenomics)
Analyzes genetic variants related to nutrition, metabolism, and dietary needs.
Based on peer-reviewed research for personalized nutrition insights.
"""

from typing import Dict, List, Tuple, Optional
from dataclasses import dataclass, field


# =============================================================================
# Strand complement mapping
# =============================================================================

COMPLEMENT = {"A": "T", "T": "A", "C": "G", "G": "C"}


def complement_genotype(genotype: str) -> str:
    """Return the complement of a genotype string (e.g., 'AG' -> 'TC')."""
    return "".join(COMPLEMENT.get(b, b) for b in genotype.upper())


def normalize_genotype(genotype: str) -> str:
    """Sort alleles alphabetically so AG == GA for comparison."""
    g = genotype.upper().replace("-", "")
    if len(g) == 2:
        return "".join(sorted(g))
    return g


# =============================================================================
# Data model
# =============================================================================

@dataclass
class GenotypeInterpretation:
    label: str
    meaning: str
    recommendation: str


@dataclass
class NutritionSNP:
    rsid: str
    gene: str
    category: str
    nutrient: str
    ref_allele: str
    alt_allele: str
    genotype_interpretations: Dict[str, GenotypeInterpretation]
    references: List[str] = field(default_factory=list)


@dataclass
class NutritionVariantResult:
    rsid: str
    gene: str
    category: str
    nutrient: str
    genotype: str
    label: str
    meaning: str
    recommendation: str


@dataclass
class NutritionAnalysisResult:
    total_checked: int
    found: int
    not_found: int
    results_by_category: Dict[str, List[NutritionVariantResult]]


# =============================================================================
# Categories
# =============================================================================

NUTRITION_CATEGORIES = [
    "Macronutrients",
    "Micronutrients",
    "Substance Metabolism",
    "Taste & Sensitivity",
]

CATEGORY_EMOJIS = {
    "Macronutrients": "\U0001F95B",       # 🥛
    "Micronutrients": "\U0001F34E",       # 🍎
    "Substance Metabolism": "\u2615",     # ☕
    "Taste & Sensitivity": "\U0001F445", # 👅
}


# =============================================================================
# SNP Database
# =============================================================================

NUTRITION_DATABASE: Dict[str, NutritionSNP] = {

    # =========================================================================
    # MACRONUTRIENTS
    # =========================================================================

    "rs4988235": NutritionSNP(
        rsid="rs4988235",
        gene="LCT/MCM6",
        category="Macronutrients",
        nutrient="Lactose Tolerance",
        ref_allele="G",
        alt_allele="A",
        genotype_interpretations={
            "GG": GenotypeInterpretation(
                label="Likely Lactose Intolerant",
                meaning="You likely do not produce lactase in adulthood. Most adults worldwide carry this genotype. Dairy consumption may cause bloating, gas, or digestive discomfort.",
                recommendation="Consider lactose-free dairy products, plant-based milks, or a lactase enzyme supplement before consuming dairy. Ensure adequate calcium from non-dairy sources such as leafy greens, fortified foods, and sardines.",
            ),
            "AG": GenotypeInterpretation(
                label="Likely Lactose Tolerant",
                meaning="You carry one copy of the lactase persistence allele. You most likely continue to produce lactase into adulthood and can digest dairy without issue.",
                recommendation="Dairy is generally well-tolerated. Include a variety of dairy and non-dairy calcium sources for optimal nutrition.",
            ),
            "AA": GenotypeInterpretation(
                label="Lactose Tolerant",
                meaning="You carry two copies of the lactase persistence allele. You continue to produce lactase and can digest lactose efficiently.",
                recommendation="Dairy is well-tolerated. Enjoy a balanced intake of dairy products as desired.",
            ),
        },
        references=["PMID:11788828", "PMID:15114531"],
    ),

    "rs9939609": NutritionSNP(
        rsid="rs9939609",
        gene="FTO",
        category="Macronutrients",
        nutrient="Obesity & Appetite Regulation",
        ref_allele="T",
        alt_allele="A",
        genotype_interpretations={
            "TT": GenotypeInterpretation(
                label="Typical Risk",
                meaning="You do not carry the FTO risk allele. Your genetic predisposition for increased appetite and weight gain from this variant is not elevated.",
                recommendation="Maintain a balanced diet and regular physical activity. Standard dietary guidelines apply.",
            ),
            "AT": GenotypeInterpretation(
                label="Slightly Increased Risk",
                meaning="You carry one copy of the FTO risk allele. This is associated with roughly 1.2x increased risk of obesity and slightly higher appetite/caloric intake.",
                recommendation="Be mindful of portion sizes and satiety cues. Regular physical activity can offset much of the genetic predisposition. A higher-protein diet may help with satiety.",
            ),
            "AA": GenotypeInterpretation(
                label="Increased Risk",
                meaning="You carry two copies of the FTO risk allele. This is associated with roughly 1.7x increased obesity risk and higher appetite drive. On average, carriers eat ~100 kcal/day more.",
                recommendation="Prioritize high-protein, high-fiber meals to promote satiety. Regular exercise (especially resistance training) significantly mitigates this risk. Consider working with a dietitian for personalized guidance.",
            ),
        },
        references=["PMID:17434869", "PMID:18454148"],
    ),

    "rs7903146": NutritionSNP(
        rsid="rs7903146",
        gene="TCF7L2",
        category="Macronutrients",
        nutrient="Carbohydrate Metabolism & T2D Risk",
        ref_allele="C",
        alt_allele="T",
        genotype_interpretations={
            "CC": GenotypeInterpretation(
                label="Typical Risk",
                meaning="You do not carry the TCF7L2 risk allele. Your carbohydrate metabolism and type 2 diabetes risk from this locus are not elevated.",
                recommendation="Follow standard dietary guidelines. A balanced diet with whole grains, vegetables, and moderate refined carbohydrate intake is recommended.",
            ),
            "CT": GenotypeInterpretation(
                label="Moderately Increased T2D Risk",
                meaning="You carry one copy of the T risk allele, associated with roughly 1.4x increased risk of type 2 diabetes. This variant affects insulin secretion in response to carbohydrates.",
                recommendation="Favour complex carbohydrates with a low glycaemic index. Limit refined sugars and processed carbs. Regular blood glucose monitoring is advisable, especially if other risk factors are present.",
            ),
            "TT": GenotypeInterpretation(
                label="Higher T2D Risk",
                meaning="You carry two copies of the T risk allele, associated with roughly 2x increased risk of type 2 diabetes. Insulin secretion in response to carbohydrate intake may be impaired.",
                recommendation="Strongly consider a low-glycaemic, fibre-rich diet. Limit refined carbohydrates and sugary beverages. Regular glucose screening is recommended. Physical activity significantly reduces risk.",
            ),
        },
        references=["PMID:16415884", "PMID:17463246"],
    ),

    "rs1801282": NutritionSNP(
        rsid="rs1801282",
        gene="PPARG",
        category="Macronutrients",
        nutrient="Fat Metabolism",
        ref_allele="C",
        alt_allele="G",
        genotype_interpretations={
            "CC": GenotypeInterpretation(
                label="Typical Fat Metabolism",
                meaning="You carry the common Pro12Pro genotype. Your fat metabolism and insulin sensitivity from this locus are typical.",
                recommendation="Follow standard dietary guidelines for fat intake. Include healthy fats (olive oil, nuts, fatty fish) and limit saturated fats.",
            ),
            "CG": GenotypeInterpretation(
                label="Improved Insulin Sensitivity",
                meaning="You carry one copy of the Ala allele (Pro12Ala). This is associated with improved insulin sensitivity and a modest protective effect against type 2 diabetes.",
                recommendation="A Mediterranean-style diet rich in unsaturated fats may complement your genetic advantage. Maintain an active lifestyle to maximise insulin sensitivity.",
            ),
            "GG": GenotypeInterpretation(
                label="Enhanced Insulin Sensitivity",
                meaning="You carry two copies of the Ala allele (Ala12Ala). This rare genotype is linked to higher insulin sensitivity, but may also be associated with slightly higher BMI in some populations.",
                recommendation="Focus on a balanced fat intake with emphasis on unsaturated sources. Monitor weight and metabolic markers. A Mediterranean dietary pattern is advisable.",
            ),
        },
        references=["PMID:9614169", "PMID:12829785"],
    ),

    # =========================================================================
    # MICRONUTRIENTS
    # =========================================================================

    "rs1801133": NutritionSNP(
        rsid="rs1801133",
        gene="MTHFR",
        category="Micronutrients",
        nutrient="Folate Metabolism (C677T)",
        ref_allele="G",
        alt_allele="A",
        genotype_interpretations={
            "GG": GenotypeInterpretation(
                label="Normal Folate Metabolism",
                meaning="You have the common 677CC genotype (reported on plus strand as GG). MTHFR enzyme activity is normal, and folate metabolism is efficient.",
                recommendation="Standard folate intake through leafy greens, legumes, and fortified foods is sufficient.",
            ),
            "AG": GenotypeInterpretation(
                label="Mildly Reduced MTHFR Activity",
                meaning="You are heterozygous for the C677T variant (roughly 65% normal enzyme activity). Folate metabolism is mildly impaired, and homocysteine may be slightly elevated.",
                recommendation="Ensure adequate dietary folate (400-600 mcg DFE daily). Consider foods rich in methylfolate (leafy greens, lentils). Supplementation with methylfolate (L-5-MTHF) rather than folic acid may be more efficient.",
            ),
            "AA": GenotypeInterpretation(
                label="Significantly Reduced MTHFR Activity",
                meaning="You are homozygous for the C677T variant (roughly 30% normal enzyme activity). Folate metabolism is substantially impaired and homocysteine levels may be elevated, which is a cardiovascular risk factor.",
                recommendation="Prioritise methylfolate-rich foods and consider supplementation with L-5-MTHF (methylfolate) rather than folic acid. Ensure adequate B6 and B12 intake. Have homocysteine levels monitored by your doctor.",
            ),
        },
        references=["PMID:7647779", "PMID:9545397"],
    ),

    "rs1801131": NutritionSNP(
        rsid="rs1801131",
        gene="MTHFR",
        category="Micronutrients",
        nutrient="Folate Metabolism (A1298C)",
        ref_allele="T",
        alt_allele="G",
        genotype_interpretations={
            "TT": GenotypeInterpretation(
                label="Normal MTHFR (1298 position)",
                meaning="You have the common 1298AA genotype (reported on plus strand as TT). Enzyme activity at this position is not impaired.",
                recommendation="Standard folate intake is generally sufficient from this variant alone.",
            ),
            "GT": GenotypeInterpretation(
                label="Mildly Reduced MTHFR Activity (A1298C)",
                meaning="You are heterozygous for the A1298C variant. This causes a mild reduction in MTHFR activity. Impact on folate metabolism is modest in isolation but may be significant if combined with the C677T variant.",
                recommendation="Ensure adequate dietary folate. If you also carry the C677T variant, consider methylfolate supplementation and homocysteine testing.",
            ),
            "GG": GenotypeInterpretation(
                label="Reduced MTHFR Activity (A1298C)",
                meaning="You are homozygous for the A1298C variant. MTHFR enzyme activity is moderately reduced. The clinical significance is generally less than C677T homozygosity, but the combination can be impactful.",
                recommendation="Increase folate-rich food intake. Consider methylfolate supplementation, especially if combined with the C677T variant. Monitor homocysteine levels.",
            ),
        },
        references=["PMID:10444342", "PMID:11814284"],
    ),

    "rs7041": NutritionSNP(
        rsid="rs7041",
        gene="GC",
        category="Micronutrients",
        nutrient="Vitamin D Levels",
        ref_allele="G",
        alt_allele="T",
        genotype_interpretations={
            "GG": GenotypeInterpretation(
                label="Higher Vitamin D Binding",
                meaning="You carry the GC*1F variant associated with higher vitamin D binding protein levels and generally adequate circulating vitamin D.",
                recommendation="Standard sun exposure and dietary vitamin D intake are usually sufficient. Monitor levels annually.",
            ),
            "GT": GenotypeInterpretation(
                label="Moderate Vitamin D Levels",
                meaning="You are heterozygous at this locus. Vitamin D binding protein levels and circulating 25(OH)D may be intermediate.",
                recommendation="Ensure adequate sun exposure and dietary vitamin D (fatty fish, fortified foods). Consider supplementation during winter months, especially at higher latitudes.",
            ),
            "TT": GenotypeInterpretation(
                label="Lower Vitamin D Levels",
                meaning="You carry two copies of the T allele (GC*2 variant), associated with lower vitamin D binding protein and reduced circulating 25(OH)D levels.",
                recommendation="Regular vitamin D supplementation (1000-2000 IU daily) is advisable, especially in winter or with limited sun exposure. Have 25(OH)D levels tested and supplement accordingly.",
            ),
        },
        references=["PMID:20541252", "PMID:22573406"],
    ),

    "rs12934922": NutritionSNP(
        rsid="rs12934922",
        gene="BCMO1",
        category="Micronutrients",
        nutrient="Beta-Carotene to Vitamin A Conversion",
        ref_allele="A",
        alt_allele="T",
        genotype_interpretations={
            "AA": GenotypeInterpretation(
                label="Normal Conversion",
                meaning="You efficiently convert beta-carotene (from plant foods) into active vitamin A (retinol). Your BCMO1 enzyme functions normally.",
                recommendation="A diet rich in colourful vegetables and fruits (carrots, sweet potatoes, spinach) provides adequate vitamin A through beta-carotene conversion.",
            ),
            "AT": GenotypeInterpretation(
                label="Reduced Conversion",
                meaning="You carry one copy of the variant allele, which reduces beta-carotene to vitamin A conversion by approximately 32%. You may have slightly lower vitamin A status if relying solely on plant sources.",
                recommendation="Include both plant sources of beta-carotene and preformed vitamin A (retinol) from eggs, dairy, or liver. Consider vitamin A supplementation if following a strict vegan diet.",
            ),
            "TT": GenotypeInterpretation(
                label="Significantly Reduced Conversion",
                meaning="You carry two copies of the variant allele, reducing beta-carotene conversion by approximately 69%. You are a 'low converter' and may not obtain sufficient vitamin A from plant sources alone.",
                recommendation="Include preformed vitamin A (retinol) from animal sources such as eggs, dairy, and liver. If vegan, consider retinol supplementation or monitor vitamin A status. Do not rely solely on beta-carotene.",
            ),
        },
        references=["PMID:19103647", "PMID:22113863"],
    ),

    "rs602662": NutritionSNP(
        rsid="rs602662",
        gene="FUT2",
        category="Micronutrients",
        nutrient="Vitamin B12 Absorption",
        ref_allele="G",
        alt_allele="A",
        genotype_interpretations={
            "GG": GenotypeInterpretation(
                label="Normal B12 Absorption",
                meaning="You are a 'secretor' with normal FUT2 function. Your gut microbiome composition supports typical vitamin B12 absorption.",
                recommendation="Standard dietary B12 intake from animal products, fortified foods, or supplementation (for vegetarians/vegans) is adequate.",
            ),
            "AG": GenotypeInterpretation(
                label="Slightly Reduced B12 Levels",
                meaning="You carry one copy of the non-secretor allele. B12 levels may be slightly lower than average due to altered gut microbiome composition.",
                recommendation="Ensure regular B12 intake. If vegetarian or vegan, supplementation is particularly important. Monitor B12 levels periodically.",
            ),
            "AA": GenotypeInterpretation(
                label="Lower B12 Levels Likely",
                meaning="You are a 'non-secretor' with reduced FUT2 function. This is associated with lower circulating vitamin B12 levels and altered gut microbiome.",
                recommendation="Consider B12 supplementation (methylcobalamin or hydroxocobalamin form). Regular B12 level monitoring is recommended. Include B12-rich foods or fortified products consistently.",
            ),
        },
        references=["PMID:19744961", "PMID:23754956"],
    ),

    "rs1800562": NutritionSNP(
        rsid="rs1800562",
        gene="HFE",
        category="Micronutrients",
        nutrient="Iron Absorption (C282Y)",
        ref_allele="G",
        alt_allele="A",
        genotype_interpretations={
            "GG": GenotypeInterpretation(
                label="Normal Iron Regulation",
                meaning="You do not carry the C282Y variant. Iron absorption and regulation are typical.",
                recommendation="Follow standard dietary guidelines for iron intake. Include iron-rich foods with vitamin C to enhance absorption.",
            ),
            "AG": GenotypeInterpretation(
                label="Carrier - Mildly Increased Iron Absorption",
                meaning="You carry one copy of the C282Y variant. You are a carrier for hereditary hemochromatosis. Iron absorption may be slightly elevated, but clinical iron overload is uncommon in heterozygotes.",
                recommendation="Periodic serum ferritin monitoring is prudent. Avoid unnecessary iron supplementation. Moderate red meat and alcohol intake. If combined with H63D, discuss with your doctor.",
            ),
            "AA": GenotypeInterpretation(
                label="High Risk - Hereditary Hemochromatosis",
                meaning="You are homozygous for C282Y, the primary variant for hereditary hemochromatosis. You are at significant risk of iron overload, which can damage the liver, heart, and pancreas if untreated.",
                recommendation="Seek medical evaluation and regular serum ferritin and transferrin saturation testing. Avoid iron supplements and vitamin C supplements with meals. Limit red meat. Therapeutic phlebotomy may be recommended.",
            ),
        },
        references=["PMID:8696333", "PMID:11196646"],
    ),

    "rs1799945": NutritionSNP(
        rsid="rs1799945",
        gene="HFE",
        category="Micronutrients",
        nutrient="Iron Absorption (H63D)",
        ref_allele="C",
        alt_allele="G",
        genotype_interpretations={
            "CC": GenotypeInterpretation(
                label="Normal Iron Regulation",
                meaning="You do not carry the H63D variant. Iron absorption at this locus is typical.",
                recommendation="Standard dietary guidelines for iron apply.",
            ),
            "CG": GenotypeInterpretation(
                label="Carrier - Mild Effect on Iron",
                meaning="You carry one copy of the H63D variant. This has a mild effect on iron absorption. Clinical significance is limited in isolation but may be relevant if combined with C282Y.",
                recommendation="No specific action needed unless also a C282Y carrier. If compound heterozygous (C282Y/H63D), consult your doctor for iron studies.",
            ),
            "GG": GenotypeInterpretation(
                label="Mildly Increased Iron Absorption",
                meaning="You are homozygous for H63D. This may modestly increase iron absorption. Risk of clinical hemochromatosis from H63D alone is low but not negligible.",
                recommendation="Consider periodic ferritin monitoring. Avoid unnecessary iron supplementation. Discuss with your doctor if ferritin levels are elevated.",
            ),
        },
        references=["PMID:8696333", "PMID:12472866"],
    ),

    # =========================================================================
    # SUBSTANCE METABOLISM
    # =========================================================================

    "rs762551": NutritionSNP(
        rsid="rs762551",
        gene="CYP1A2",
        category="Substance Metabolism",
        nutrient="Caffeine Metabolism",
        ref_allele="A",
        alt_allele="C",
        genotype_interpretations={
            "AA": GenotypeInterpretation(
                label="Fast Caffeine Metaboliser",
                meaning="You are a rapid caffeine metaboliser (CYP1A2*1A/*1A). You clear caffeine quickly and may tolerate higher amounts without adverse effects. Moderate coffee intake (3-4 cups) has been associated with reduced cardiovascular risk in fast metabolisers.",
                recommendation="Moderate coffee consumption (up to 3-4 cups/day) is generally safe and may be protective. Be aware you may need more caffeine for the desired effect due to rapid clearance.",
            ),
            "AC": GenotypeInterpretation(
                label="Intermediate Caffeine Metaboliser",
                meaning="You are an intermediate caffeine metaboliser (heterozygous CYP1A2*1A/*1F). Caffeine clearance is moderate.",
                recommendation="Moderate coffee intake (2-3 cups/day) is generally fine. Avoid caffeine late in the day if it affects your sleep. Monitor how you feel after coffee to calibrate intake.",
            ),
            "CC": GenotypeInterpretation(
                label="Slow Caffeine Metaboliser",
                meaning="You are a slow caffeine metaboliser (CYP1A2*1F/*1F). Caffeine stays in your system longer. High coffee intake (>2-3 cups/day) has been associated with increased risk of hypertension and heart attack in slow metabolisers.",
                recommendation="Limit coffee to 1-2 cups/day. Avoid caffeine in the afternoon and evening. Consider switching to decaf or tea. Be aware of caffeine in chocolate, energy drinks, and medications.",
            ),
        },
        references=["PMID:16522833", "PMID:18226225"],
    ),

    "rs4244285": NutritionSNP(
        rsid="rs4244285",
        gene="CYP2C19",
        category="Substance Metabolism",
        nutrient="Drug & Substance Metabolism",
        ref_allele="G",
        alt_allele="A",
        genotype_interpretations={
            "GG": GenotypeInterpretation(
                label="Normal Metaboliser (CYP2C19)",
                meaning="You are an extensive (normal) metaboliser for CYP2C19 substrates. This enzyme metabolises several medications (e.g., clopidogrel, omeprazole, some antidepressants) as well as certain dietary compounds.",
                recommendation="Standard dosing of CYP2C19-metabolised medications is appropriate. No special dietary considerations from this variant.",
            ),
            "AG": GenotypeInterpretation(
                label="Intermediate Metaboliser (CYP2C19)",
                meaning="You carry one loss-of-function allele (*2). You are an intermediate metaboliser. Some medications metabolised by CYP2C19 may have altered efficacy or side effect profiles.",
                recommendation="Inform your healthcare provider about this variant when prescribing CYP2C19 substrates (clopidogrel, PPIs, certain SSRIs). Pharmacogenomic-guided dosing may be beneficial.",
            ),
            "AA": GenotypeInterpretation(
                label="Poor Metaboliser (CYP2C19)",
                meaning="You carry two loss-of-function alleles (*2/*2). You are a poor metaboliser. CYP2C19-dependent medications may accumulate or fail to activate properly (e.g., clopidogrel may be less effective).",
                recommendation="Inform your doctor of this result. Alternative medications or adjusted dosing may be needed for CYP2C19 substrates. Carry this information on a pharmacogenomics card if available.",
            ),
        },
        references=["PMID:16625025", "PMID:21270786"],
    ),

    "rs671": NutritionSNP(
        rsid="rs671",
        gene="ALDH2",
        category="Substance Metabolism",
        nutrient="Alcohol Metabolism (Acetaldehyde Clearance)",
        ref_allele="G",
        alt_allele="A",
        genotype_interpretations={
            "GG": GenotypeInterpretation(
                label="Normal Aldehyde Metabolism",
                meaning="You have functional ALDH2 enzyme. Acetaldehyde (the toxic intermediate of alcohol metabolism) is cleared efficiently.",
                recommendation="Standard alcohol guidelines apply (moderation). No increased sensitivity from this variant.",
            ),
            "AG": GenotypeInterpretation(
                label="Reduced Aldehyde Metabolism (Alcohol Flush)",
                meaning="You carry one copy of the ALDH2*2 variant. Your ability to clear acetaldehyde is reduced. You likely experience facial flushing, nausea, and rapid heart rate after alcohol (Asian flush reaction). Acetaldehyde accumulation increases esophageal cancer risk with alcohol use.",
                recommendation="Limit alcohol consumption significantly. If you drink, do so minimally and slowly. Be aware of the elevated cancer risk from acetaldehyde accumulation. Avoid combining alcohol with medications.",
            ),
            "AA": GenotypeInterpretation(
                label="Severely Impaired Aldehyde Metabolism",
                meaning="You are homozygous for ALDH2*2. Your ALDH2 enzyme is essentially non-functional. Even small amounts of alcohol cause severe flushing, nausea, and discomfort. Acetaldehyde accumulation poses a significant cancer risk.",
                recommendation="Avoid alcohol entirely. Even small amounts can cause severe reactions and long-term health risks (esophageal and head/neck cancers). Be cautious with fermented foods and alcohol in cooking.",
            ),
        },
        references=["PMID:19276327", "PMID:20304485"],
    ),

    "rs1229984": NutritionSNP(
        rsid="rs1229984",
        gene="ADH1B",
        category="Substance Metabolism",
        nutrient="Alcohol Metabolism (Ethanol Oxidation)",
        ref_allele="T",
        alt_allele="C",
        genotype_interpretations={
            "TT": GenotypeInterpretation(
                label="Typical Alcohol Metabolism",
                meaning="You carry the common ADH1B*1 variant. Alcohol (ethanol) is converted to acetaldehyde at a normal rate.",
                recommendation="Standard alcohol moderation guidelines apply. Be mindful of overall intake and its health effects.",
            ),
            "CT": GenotypeInterpretation(
                label="Faster Alcohol Metabolism",
                meaning="You carry one copy of the ADH1B*2 (His47Arg) variant. You convert ethanol to acetaldehyde more rapidly, which can cause faster accumulation of acetaldehyde after drinking. This variant is protective against alcoholism.",
                recommendation="You may experience unpleasant effects from alcohol more quickly. This can be protective against excessive drinking. Moderate your alcohol intake accordingly.",
            ),
            "CC": GenotypeInterpretation(
                label="Rapid Alcohol Metabolism",
                meaning="You are homozygous for the ADH1B*2 variant. Ethanol is converted to acetaldehyde very rapidly. You likely experience stronger reactions to alcohol. This variant is strongly protective against alcohol dependence.",
                recommendation="Alcohol tolerance is likely low. Respect your body's signals and limit consumption. This variant is protective against alcohol dependence but does not protect against alcohol-related organ damage if you do drink.",
            ),
        },
        references=["PMID:15057824", "PMID:19261174"],
    ),

    # =========================================================================
    # TASTE & SENSITIVITY
    # =========================================================================

    "rs713598": NutritionSNP(
        rsid="rs713598",
        gene="TAS2R38",
        category="Taste & Sensitivity",
        nutrient="Bitter Taste Perception (Position 49)",
        ref_allele="G",
        alt_allele="C",
        genotype_interpretations={
            "GG": GenotypeInterpretation(
                label="Likely Non-Taster",
                meaning="You carry the AVI haplotype alleles at this position. Combined with other TAS2R38 variants, you likely have reduced sensitivity to bitter compounds such as PTC, PROP, and glucosinolates found in cruciferous vegetables.",
                recommendation="You may find bitter vegetables more palatable and can include them easily in your diet. Enjoy broccoli, kale, Brussels sprouts, and other bitter greens freely.",
            ),
            "CG": GenotypeInterpretation(
                label="Moderate Bitter Taster",
                meaning="You are heterozygous at this position. You likely have intermediate sensitivity to bitter compounds. Some bitter foods may be mildly unpleasant.",
                recommendation="Try preparing bitter vegetables with seasoning, olive oil, or roasting to reduce bitterness while retaining nutritional benefits.",
            ),
            "CC": GenotypeInterpretation(
                label="Likely Strong Bitter Taster",
                meaning="You carry the PAV haplotype alleles at this position. Combined with other TAS2R38 variants, you likely have heightened sensitivity to bitter compounds. Cruciferous vegetables and some foods may taste very bitter.",
                recommendation="Prepare bitter vegetables with flavour-masking techniques (roasting, cheese sauce, spices). Despite taste aversion, aim to include these cancer-protective vegetables in your diet through creative preparation.",
            ),
        },
        references=["PMID:12595690", "PMID:16051168"],
    ),

    "rs1726866": NutritionSNP(
        rsid="rs1726866",
        gene="TAS2R38",
        category="Taste & Sensitivity",
        nutrient="Bitter Taste Perception (Position 262)",
        ref_allele="G",
        alt_allele="A",
        genotype_interpretations={
            "GG": GenotypeInterpretation(
                label="Likely Non-Taster (Position 262)",
                meaning="You carry the AVI haplotype allele at position 262. This contributes to reduced bitter taste perception when combined with concordant variants at other TAS2R38 positions.",
                recommendation="Bitter vegetables are likely palatable. Include a variety of cruciferous vegetables for their cancer-protective properties.",
            ),
            "AG": GenotypeInterpretation(
                label="Intermediate Bitter Taster (Position 262)",
                meaning="You are heterozygous at this position, contributing to intermediate bitter taste sensitivity.",
                recommendation="Experiment with different preparation methods for bitter vegetables. Roasting, sauteing with garlic, or adding citrus can reduce perceived bitterness.",
            ),
            "AA": GenotypeInterpretation(
                label="Likely Strong Bitter Taster (Position 262)",
                meaning="You carry the PAV haplotype allele at position 262. Combined with concordant variants at other TAS2R38 positions, this contributes to strong bitter taste perception.",
                recommendation="Use flavour-masking strategies for cruciferous vegetables. Smoothies, stir-fries with strong sauces, or raw preparations with dips can help. These vegetables are important for cancer prevention.",
            ),
        },
        references=["PMID:12595690", "PMID:16051168"],
    ),

    "rs17822931": NutritionSNP(
        rsid="rs17822931",
        gene="ABCC11",
        category="Taste & Sensitivity",
        nutrient="Earwax Type & Body Odour",
        ref_allele="C",
        alt_allele="T",
        genotype_interpretations={
            "CC": GenotypeInterpretation(
                label="Wet Earwax / Typical Body Odour",
                meaning="You have the common variant producing wet, sticky earwax. You also have typical apocrine gland secretion, meaning standard body odour production.",
                recommendation="Standard hygiene practices apply. Antiperspirant/deodorant use is recommended as needed. This has no direct dietary implication but is an interesting genetic trait.",
            ),
            "CT": GenotypeInterpretation(
                label="Intermediate Earwax / Reduced Body Odour",
                meaning="You carry one copy of the dry earwax allele. Earwax consistency may be intermediate. Body odour production is somewhat reduced due to decreased apocrine gland secretion.",
                recommendation="You may need less antiperspirant than average. No specific dietary considerations.",
            ),
            "TT": GenotypeInterpretation(
                label="Dry Earwax / Minimal Body Odour",
                meaning="You are homozygous for the dry earwax variant. You produce dry, flaky earwax and have significantly reduced apocrine gland secretion, resulting in minimal body odour. This variant is common in East Asian populations.",
                recommendation="You likely need little to no antiperspirant. No dietary implications, but this is a useful pharmacogenomic marker.",
            ),
        },
        references=["PMID:16444273", "PMID:19710689"],
    ),
}


# =============================================================================
# Analysis functions
# =============================================================================

def _match_genotype(observed: str, snp: NutritionSNP) -> Optional[GenotypeInterpretation]:
    """
    Match an observed genotype against the known interpretations.
    Handles strand complement matching: if the observed genotype (or its
    normalized form) is not found directly, try the complement.
    """
    norm = normalize_genotype(observed)

    # Build a lookup with normalized keys
    lookup: Dict[str, GenotypeInterpretation] = {}
    for gt, interp in snp.genotype_interpretations.items():
        lookup[normalize_genotype(gt)] = interp

    # Direct match
    if norm in lookup:
        return lookup[norm]

    # Complement match
    comp = normalize_genotype(complement_genotype(observed))
    if comp in lookup:
        return lookup[comp]

    return None


def analyze_nutrition(
    variants: Dict[str, Tuple[str, str, str]],
) -> NutritionAnalysisResult:
    """
    Analyze genetic variants for nutrigenomics insights.

    Args:
        variants: Dictionary mapping rsID (lowercase) to
                  (chromosome, position, genotype).

    Returns:
        NutritionAnalysisResult with categorized findings.
    """
    results_by_category: Dict[str, List[NutritionVariantResult]] = {
        cat: [] for cat in NUTRITION_CATEGORIES
    }
    found = 0
    not_found = 0

    for rsid, snp_info in NUTRITION_DATABASE.items():
        rsid_lower = rsid.lower()

        if rsid_lower not in variants:
            not_found += 1
            continue

        found += 1
        _chrom, _pos, genotype = variants[rsid_lower]
        genotype = genotype.upper().replace("-", "")

        interp = _match_genotype(genotype, snp_info)

        if interp is not None:
            results_by_category[snp_info.category].append(
                NutritionVariantResult(
                    rsid=rsid,
                    gene=snp_info.gene,
                    category=snp_info.category,
                    nutrient=snp_info.nutrient,
                    genotype=genotype,
                    label=interp.label,
                    meaning=interp.meaning,
                    recommendation=interp.recommendation,
                )
            )
        else:
            # Genotype present but not in our interpretations —
            # still count as found but report with a generic message.
            results_by_category[snp_info.category].append(
                NutritionVariantResult(
                    rsid=rsid,
                    gene=snp_info.gene,
                    category=snp_info.category,
                    nutrient=snp_info.nutrient,
                    genotype=genotype,
                    label="Uncommon Genotype",
                    meaning=f"Your genotype ({genotype}) at this locus is not among the most commonly studied combinations. Further research may be needed.",
                    recommendation="No specific recommendation for this genotype. Consult a healthcare professional for personalised advice.",
                )
            )

    return NutritionAnalysisResult(
        total_checked=len(NUTRITION_DATABASE),
        found=found,
        not_found=not_found,
        results_by_category=results_by_category,
    )


def generate_nutrition_json(result: NutritionAnalysisResult) -> dict:
    """
    Generate a JSON-serializable dictionary for the nutrition report.

    Returns:
        A dict with 'summary' and 'categories' keys suitable for
        serialisation and frontend rendering.
    """
    categories = []

    for cat_name in NUTRITION_CATEGORIES:
        findings = result.results_by_category.get(cat_name, [])
        if not findings:
            continue

        variants_list = []
        for v in findings:
            variants_list.append({
                "rsid": v.rsid,
                "gene": v.gene,
                "nutrient": v.nutrient,
                "genotype": v.genotype,
                "label": v.label,
                "meaning": v.meaning,
                "recommendation": v.recommendation,
            })

        categories.append({
            "name": cat_name,
            "emoji": CATEGORY_EMOJIS.get(cat_name, "\U0001F4CA"),
            "variants": variants_list,
        })

    categories_with_findings = len(categories)

    return {
        "summary": {
            "totalChecked": result.total_checked,
            "found": result.found,
            "notFound": result.not_found,
            "categoriesWithFindings": categories_with_findings,
        },
        "categories": categories,
    }
