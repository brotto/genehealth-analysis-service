"""
Comprehensive Traits & Wellness SNP Database
SNPs for cognitive traits, physical characteristics, metabolism, and more
"""

from typing import Dict, List
from dataclasses import dataclass


@dataclass
class TraitSNP:
    rsid: str
    gene: str
    category: str
    trait: str
    risk_allele: str  # Allele associated with the trait
    effect: str  # What having this allele means
    description: str
    references: List[str]  # PubMed or study references


# Comprehensive database of trait-related SNPs
TRAITS_DATABASE: Dict[str, TraitSNP] = {

    # ============================================
    # COGNITIVE & NEUROLOGICAL
    # ============================================

    # Autism Spectrum
    "rs2710102": TraitSNP(
        rsid="rs2710102",
        gene="CNTNAP2",
        category="Cognitive",
        trait="Autism Spectrum Risk",
        risk_allele="C",
        effect="Slightly increased risk of autism spectrum traits",
        description="CNTNAP2 is involved in neural development and language. Variants associated with autism and language delays.",
        references=["PMID:18179900"]
    ),
    "rs7794745": TraitSNP(
        rsid="rs7794745",
        gene="CNTNAP2",
        category="Cognitive",
        trait="Autism Spectrum Risk",
        risk_allele="T",
        effect="Associated with autism spectrum traits",
        description="Another CNTNAP2 variant linked to autism and social communication differences.",
        references=["PMID:18179900"]
    ),
    "rs1858830": TraitSNP(
        rsid="rs1858830",
        gene="MET",
        category="Cognitive",
        trait="Autism Spectrum Risk",
        risk_allele="C",
        effect="Increased autism risk (2x with CC genotype)",
        description="MET gene regulates brain development. C allele associated with autism in multiple studies.",
        references=["PMID:17053076"]
    ),
    "rs4307059": TraitSNP(
        rsid="rs4307059",
        gene="CDH9/CDH10",
        category="Cognitive",
        trait="Autism Spectrum Risk",
        risk_allele="T",
        effect="Associated with autism spectrum",
        description="Cadherin genes involved in neuronal connectivity.",
        references=["PMID:19404256"]
    ),

    # ADHD
    "rs1800955": TraitSNP(
        rsid="rs1800955",
        gene="DRD4",
        category="Cognitive",
        trait="ADHD Risk",
        risk_allele="C",
        effect="Increased ADHD risk, novelty seeking",
        description="Dopamine receptor D4 variant. Associated with attention and novelty-seeking behavior.",
        references=["PMID:11943377"]
    ),
    "rs27072": TraitSNP(
        rsid="rs27072",
        gene="DAT1/SLC6A3",
        category="Cognitive",
        trait="ADHD Risk",
        risk_allele="A",
        effect="Associated with ADHD and dopamine transport",
        description="Dopamine transporter gene. Affects dopamine reuptake in the brain.",
        references=["PMID:12893112"]
    ),
    "rs5569": TraitSNP(
        rsid="rs5569",
        gene="NET/SLC6A2",
        category="Cognitive",
        trait="ADHD Risk",
        risk_allele="T",
        effect="Associated with ADHD symptoms",
        description="Norepinephrine transporter. Affects attention and arousal systems.",
        references=["PMID:16642436"]
    ),
    "rs1801260": TraitSNP(
        rsid="rs1801260",
        gene="CLOCK",
        category="Cognitive",
        trait="ADHD Risk / Circadian Rhythm",
        risk_allele="C",
        effect="Evening preference, possible ADHD link",
        description="Circadian rhythm gene affecting sleep patterns and potentially ADHD.",
        references=["PMID:17346975"]
    ),

    # Intelligence / Cognitive Performance
    "rs363050": TraitSNP(
        rsid="rs363050",
        gene="SNAP25",
        category="Cognitive",
        trait="Cognitive Performance",
        risk_allele="G",
        effect="Associated with higher IQ scores",
        description="SNAP25 is crucial for neurotransmitter release. G allele linked to better cognitive performance.",
        references=["PMID:19197363"]
    ),
    "rs17070145": TraitSNP(
        rsid="rs17070145",
        gene="KIBRA",
        category="Cognitive",
        trait="Memory Performance",
        risk_allele="T",
        effect="Better episodic memory",
        description="T allele carriers show enhanced memory recall in multiple studies.",
        references=["PMID:16741276"]
    ),
    "rs6265": TraitSNP(
        rsid="rs6265",
        gene="BDNF",
        category="Cognitive",
        trait="Memory & Learning",
        risk_allele="T",  # Met allele
        effect="Reduced memory performance, anxiety risk",
        description="BDNF Val66Met. Met carriers may have reduced hippocampal volume and memory. Also affects mood.",
        references=["PMID:12805289", "PMID:15728823"]
    ),
    "rs4680": TraitSNP(
        rsid="rs4680",
        gene="COMT",
        category="Cognitive",
        trait="Cognitive Style / Stress Response",
        risk_allele="A",  # Met/Met
        effect="Warrior vs Worrier: A=better cognition under low stress, G=better under high stress",
        description="COMT Val158Met affects dopamine in prefrontal cortex. Met/Met (AA): better working memory but more stress sensitive. Val/Val (GG): more stress resilient but lower baseline cognition.",
        references=["PMID:11182883"]
    ),
    "rs1800497": TraitSNP(
        rsid="rs1800497",
        gene="ANKK1/DRD2",
        category="Cognitive",
        trait="Learning & Reward Processing",
        risk_allele="T",  # A1 allele
        effect="Reduced dopamine receptors, different learning style",
        description="Taq1A polymorphism. T allele (A1) linked to fewer D2 receptors, affecting reward-based learning.",
        references=["PMID:9449002"]
    ),

    # ============================================
    # METABOLISM & NUTRITION
    # ============================================

    # Caffeine
    "rs762551": TraitSNP(
        rsid="rs762551",
        gene="CYP1A2",
        category="Metabolism",
        trait="Caffeine Metabolism",
        risk_allele="C",
        effect="Slow caffeine metabolizer",
        description="AA genotype = fast metabolizer (can drink more coffee). AC/CC = slow metabolizer (caffeine stays longer, higher cardiovascular risk with high intake).",
        references=["PMID:16522833"]
    ),
    "rs2472297": TraitSNP(
        rsid="rs2472297",
        gene="CYP1A2",
        category="Metabolism",
        trait="Caffeine Consumption",
        risk_allele="T",
        effect="Higher caffeine consumption tendency",
        description="Affects how much caffeine you naturally tend to consume.",
        references=["PMID:21490707"]
    ),
    "rs4410790": TraitSNP(
        rsid="rs4410790",
        gene="AHR",
        category="Metabolism",
        trait="Caffeine Sensitivity",
        risk_allele="T",
        effect="Higher caffeine sensitivity",
        description="Affects caffeine's effects on the body, including anxiety and sleep disruption.",
        references=["PMID:25288136"]
    ),

    # Alcohol
    "rs671": TraitSNP(
        rsid="rs671",
        gene="ALDH2",
        category="Metabolism",
        trait="Alcohol Metabolism",
        risk_allele="A",
        effect="Alcohol flush reaction, poor alcohol tolerance",
        description="Common in East Asian populations. A allele causes acetaldehyde buildup, flushing, and nausea.",
        references=["PMID:20626054"]
    ),
    "rs1229984": TraitSNP(
        rsid="rs1229984",
        gene="ADH1B",
        category="Metabolism",
        trait="Alcohol Metabolism",
        risk_allele="T",
        effect="Fast alcohol metabolism, protective against alcoholism",
        description="Processes alcohol faster, making intoxication less pleasant. Protective against alcohol dependence.",
        references=["PMID:18385738"]
    ),

    # Lactose
    "rs4988235": TraitSNP(
        rsid="rs4988235",
        gene="MCM6/LCT",
        category="Metabolism",
        trait="Lactose Tolerance",
        risk_allele="G",
        effect="Lactose intolerance (GG genotype)",
        description="G/G = likely lactose intolerant. A/G or A/A = lactose tolerant (lactase persistence).",
        references=["PMID:11788828"]
    ),

    # Vitamin & Nutrient Metabolism
    "rs12934922": TraitSNP(
        rsid="rs12934922",
        gene="BCMO1",
        category="Metabolism",
        trait="Beta-Carotene Conversion",
        risk_allele="T",
        effect="Poor converter of beta-carotene to Vitamin A",
        description="May need preformed Vitamin A (retinol) instead of relying on plant sources.",
        references=["PMID:19103647"]
    ),
    "rs7041": TraitSNP(
        rsid="rs7041",
        gene="GC",
        category="Metabolism",
        trait="Vitamin D Levels",
        risk_allele="T",
        effect="Lower vitamin D levels",
        description="Affects vitamin D binding protein. May need more sun or supplementation.",
        references=["PMID:20418485"]
    ),
    "rs602662": TraitSNP(
        rsid="rs602662",
        gene="FUT2",
        category="Metabolism",
        trait="Vitamin B12 Absorption",
        risk_allele="A",
        effect="Lower B12 levels",
        description="Affects B12 absorption. May benefit from monitoring B12 status.",
        references=["PMID:18779456"]
    ),

    # ============================================
    # SLEEP & ENERGY
    # ============================================

    "rs57875989": TraitSNP(
        rsid="rs57875989",
        gene="DEC2/BHLHE41",
        category="Sleep",
        trait="Sleep Duration",
        risk_allele="G",
        effect="Short sleeper (needs less sleep)",
        description="Rare variant allowing some people to function well on 4-6 hours of sleep.",
        references=["PMID:19679812"]
    ),
    "rs12649507": TraitSNP(
        rsid="rs12649507",
        gene="ADA",
        category="Sleep",
        trait="Sleep Depth",
        risk_allele="A",
        effect="Deeper sleep, more slow-wave sleep",
        description="Affects adenosine metabolism which regulates sleep pressure.",
        references=["PMID:15931224"]
    ),
    "rs1801260_clock": TraitSNP(  # Note: same as ADHD one, different effect
        rsid="rs1801260",
        gene="CLOCK",
        category="Sleep",
        trait="Chronotype (Morning/Evening)",
        risk_allele="C",
        effect="Evening preference (night owl)",
        description="C allele associated with being a night owl. T allele with morning preference.",
        references=["PMID:12957359"]
    ),
    "rs10830963": TraitSNP(
        rsid="rs10830963",
        gene="MTNR1B",
        category="Sleep",
        trait="Melatonin Sensitivity",
        risk_allele="G",
        effect="Higher fasting glucose, sleep disruption effects on metabolism",
        description="Affects melatonin receptor. Night eating/late meals may particularly affect blood sugar.",
        references=["PMID:19060906"]
    ),

    # Fatigue & Energy
    "rs9939609": TraitSNP(
        rsid="rs9939609",
        gene="FTO",
        category="Metabolism",
        trait="Energy & Weight Regulation",
        risk_allele="A",
        effect="Increased appetite, obesity risk",
        description="FTO affects hunger hormones. A allele carriers may feel less satiated after eating.",
        references=["PMID:17434869"]
    ),

    # ============================================
    # PHYSICAL TRAITS
    # ============================================

    # Hair
    "rs1805007": TraitSNP(
        rsid="rs1805007",
        gene="MC1R",
        category="Physical",
        trait="Red Hair",
        risk_allele="T",
        effect="Red hair, fair skin, increased sun sensitivity",
        description="One of several MC1R variants causing red hair phenotype.",
        references=["PMID:8651275"]
    ),
    "rs1805008": TraitSNP(
        rsid="rs1805008",
        gene="MC1R",
        category="Physical",
        trait="Red Hair",
        risk_allele="T",
        effect="Red hair, fair skin",
        description="Another MC1R variant for red hair.",
        references=["PMID:8651275"]
    ),
    "rs12821256": TraitSNP(
        rsid="rs12821256",
        gene="KITLG",
        category="Physical",
        trait="Blonde Hair",
        risk_allele="C",
        effect="Blonde hair",
        description="Affects hair color through melanocyte development.",
        references=["PMID:25283252"]
    ),

    # Baldness
    "rs2180439": TraitSNP(
        rsid="rs2180439",
        gene="Chr20p11",
        category="Physical",
        trait="Male Pattern Baldness",
        risk_allele="T",
        effect="Increased baldness risk",
        description="One of several variants associated with androgenetic alopecia.",
        references=["PMID:18849991"]
    ),
    "rs6152": TraitSNP(
        rsid="rs6152",
        gene="AR",
        category="Physical",
        trait="Male Pattern Baldness",
        risk_allele="A",
        effect="Increased baldness risk",
        description="Androgen receptor gene variant. Key predictor of male pattern baldness.",
        references=["PMID:15902657"]
    ),
    "rs1385699": TraitSNP(
        rsid="rs1385699",
        gene="Chr7p21.1",
        category="Physical",
        trait="Male Pattern Baldness",
        risk_allele="T",
        effect="Increased baldness risk",
        description="Genome-wide association with baldness.",
        references=["PMID:18849991"]
    ),

    # Eyes
    "rs12913832": TraitSNP(
        rsid="rs12913832",
        gene="HERC2/OCA2",
        category="Physical",
        trait="Eye Color",
        risk_allele="G",
        effect="Blue eyes (GG), Brown eyes (AA)",
        description="Major determinant of blue vs brown eyes. GG = ~80% blue. AA = brown.",
        references=["PMID:18172690"]
    ),
    "rs1800407": TraitSNP(
        rsid="rs1800407",
        gene="OCA2",
        category="Physical",
        trait="Eye Color",
        risk_allele="T",
        effect="Green/hazel eyes",
        description="Modifies eye color, associated with green and hazel.",
        references=["PMID:18172690"]
    ),
    "rs12896399": TraitSNP(
        rsid="rs12896399",
        gene="SLC24A4",
        category="Physical",
        trait="Eye Color",
        risk_allele="T",
        effect="Lighter eye color",
        description="Contributes to lighter eye pigmentation.",
        references=["PMID:18488028"]
    ),

    # Skin
    "rs1426654": TraitSNP(
        rsid="rs1426654",
        gene="SLC24A5",
        category="Physical",
        trait="Skin Pigmentation",
        risk_allele="A",
        effect="Lighter skin",
        description="Major gene for European skin lightening. A allele = lighter skin.",
        references=["PMID:16357253"]
    ),

    # ============================================
    # ATHLETIC PERFORMANCE
    # ============================================

    "rs1815739": TraitSNP(
        rsid="rs1815739",
        gene="ACTN3",
        category="Athletic",
        trait="Muscle Fiber Type",
        risk_allele="T",
        effect="Endurance athlete type (TT = no fast-twitch protein)",
        description="ACTN3 R577X. CC = sprint/power. TT = endurance. CT = mixed. Elite sprinters rarely have TT.",
        references=["PMID:12879365"]
    ),
    "rs4994": TraitSNP(
        rsid="rs4994",
        gene="ADRB3",
        category="Athletic",
        trait="Exercise Recovery",
        risk_allele="T",
        effect="May have slower recovery from exercise",
        description="Affects metabolic rate and fat oxidation during exercise.",
        references=["PMID:9771751"]
    ),
    "rs8192678": TraitSNP(
        rsid="rs8192678",
        gene="PPARGC1A",
        category="Athletic",
        trait="Aerobic Capacity",
        risk_allele="A",
        effect="Better response to endurance training",
        description="PGC-1alpha affects mitochondrial biogenesis. Important for endurance.",
        references=["PMID:14557860"]
    ),
    "rs1042713": TraitSNP(
        rsid="rs1042713",
        gene="ADRB2",
        category="Athletic",
        trait="Exercise Performance",
        risk_allele="G",
        effect="Better endurance performance",
        description="Beta-2 adrenergic receptor. Affects bronchodilation and metabolic response to exercise.",
        references=["PMID:17456243"]
    ),

    # ============================================
    # MENTAL HEALTH & PERSONALITY
    # ============================================

    "rs6295": TraitSNP(
        rsid="rs6295",
        gene="HTR1A",
        category="Mental Health",
        trait="Depression/Anxiety Risk",
        risk_allele="G",
        effect="Increased anxiety and depression risk",
        description="Serotonin receptor variant. G allele associated with mood disorders.",
        references=["PMID:16631126"]
    ),
    "rs25531": TraitSNP(
        rsid="rs25531",
        gene="SLC6A4",
        category="Mental Health",
        trait="Stress Sensitivity",
        risk_allele="G",
        effect="Lower serotonin transporter expression, more stress sensitive",
        description="Part of 5-HTTLPR. Affects how you respond to life stress.",
        references=["PMID:12869766"]
    ),
    "rs4570625": TraitSNP(
        rsid="rs4570625",
        gene="TPH2",
        category="Mental Health",
        trait="Emotional Processing",
        risk_allele="T",
        effect="Enhanced emotional reactivity",
        description="Affects serotonin synthesis in the brain. Influences emotional responses.",
        references=["PMID:16402131"]
    ),
    "rs1800532": TraitSNP(
        rsid="rs1800532",
        gene="TPH1",
        category="Mental Health",
        trait="Mood Regulation",
        risk_allele="A",
        effect="Associated with mood disorders",
        description="Tryptophan hydroxylase variant affecting serotonin production.",
        references=["PMID:10369410"]
    ),
    "rs6313": TraitSNP(
        rsid="rs6313",
        gene="HTR2A",
        category="Mental Health",
        trait="SSRI Response",
        risk_allele="T",
        effect="May respond differently to SSRI antidepressants",
        description="Serotonin receptor 2A. Affects response to serotonergic medications.",
        references=["PMID:16642437"]
    ),

    # ============================================
    # LONGEVITY & AGING
    # ============================================

    "rs429358": TraitSNP(  # Note: also in disease database
        rsid="rs429358",
        gene="APOE",
        category="Longevity",
        trait="Longevity / Alzheimer's",
        risk_allele="C",
        effect="APOE4 - reduced longevity, increased Alzheimer's risk",
        description="APOE4 carriers have higher cardiovascular and Alzheimer's risk. APOE2 is protective.",
        references=["PMID:8346443"]
    ),
    "rs2802292": TraitSNP(
        rsid="rs2802292",
        gene="FOXO3",
        category="Longevity",
        trait="Longevity",
        risk_allele="G",
        effect="Associated with longevity",
        description="FOXO3 is a key longevity gene. G allele found more often in centenarians.",
        references=["PMID:18765803"]
    ),
    "rs1042522": TraitSNP(
        rsid="rs1042522",
        gene="TP53",
        category="Longevity",
        trait="Cancer Risk / Longevity",
        risk_allele="C",
        effect="Pro72 variant - different cancer profile",
        description="p53 codon 72. Affects cancer risk profile and cellular aging.",
        references=["PMID:12474142"]
    ),

    # ============================================
    # PAIN & SENSITIVITY
    # ============================================

    "rs1799971": TraitSNP(
        rsid="rs1799971",
        gene="OPRM1",
        category="Sensitivity",
        trait="Pain Sensitivity / Opioid Response",
        risk_allele="G",
        effect="Higher pain sensitivity, may need more pain medication",
        description="Mu-opioid receptor. G allele associated with higher pain sensitivity and opioid requirements.",
        references=["PMID:15583379"]
    ),
    "rs4680_pain": TraitSNP(
        rsid="rs4680",
        gene="COMT",
        category="Sensitivity",
        trait="Pain Sensitivity",
        risk_allele="A",
        effect="Met/Met - higher pain sensitivity",
        description="COMT affects pain processing. Met/Met individuals may be more pain sensitive.",
        references=["PMID:14985765"]
    ),
}


def get_all_trait_rsids() -> set:
    """Return all rsIDs in the traits database"""
    return set(TRAITS_DATABASE.keys())


def get_traits_by_category(category: str) -> Dict[str, TraitSNP]:
    """Get all traits in a specific category"""
    return {
        rsid: trait for rsid, trait in TRAITS_DATABASE.items()
        if trait.category == category
    }


TRAIT_CATEGORIES = [
    "Cognitive",
    "Metabolism",
    "Sleep",
    "Physical",
    "Athletic",
    "Mental Health",
    "Longevity",
    "Sensitivity"
]
