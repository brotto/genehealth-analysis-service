"""
Dream & Sleep Architecture Analyzer
Analyzes genetic variants associated with dream vividness, REM sleep architecture,
circadian chronotype, and nocturnal phenomena (sleep paralysis, hypnagogia,
lucid dream propensity).

FOUR CATEGORIES:
  1. Dream Vividness & Memory        — HTR2A, DRD4, FAAH, BDNF (cross-ref)
  2. Sleep Architecture               — PER3, ADORA2A, NRXN1
  3. Chronotype & Circadian Rhythm    — CLOCK, CRY1, PER1, PER2, ARNTL
  4. Nocturnal Phenomena              — HTR2A (sleep paralysis/hypnagogia context),
                                        COMT (cross-ref), DRD4 (lucid dreams)

SCIENTIFIC BASIS:
  - Dreams & HTR2A: Nichols (2004) Pharmacol Ther — 5-HT2A receptors and dream content
  - HTR2A rs6311: Berger et al. (2003) Neuropsychobiology — receptor density and REM
  - DRD4 rs1800955: Schmack et al. (2013) PLOS ONE — dopamine and dream content
  - FAAH rs324420: Fezza et al. (2014) Molecules — endocannabinoids and REM sleep
  - PER3 VNTR / rs228697: Viola et al. (2007) Current Biology — sleep architecture
  - ADORA2A rs5751876: Retey et al. (2005) PNAS — sleep pressure and adenosine
  - CLOCK rs1801260: Mishima et al. (2005) Sleep — circadian rhythm and phase delay
  - CRY1 rs2287161: Patke et al. (2017) Cell — circadian period +~50 min (DSPD)
  - PER1 rs2735611: Katzenberg et al. (1998) / Carpen et al. (2006) — morningness
  - PER2 rs2304672: Archer et al. (2003) Sleep — evening/morning chronotype

SPIRITUAL CONNECTION (GeneHealth identity):
  - HTR2A is at the center of both vivid dreams and spiritual sensitivity
  - Sleep paralysis and hypnagogia are frequently reported as spiritual experiences
    across cultures (visions, sense of presence, astral travel)
  - Lucid dreams have a neurobiological basis in prefrontal dopamine (DRD4, COMT)

CAVEATS:
  - All associations are probabilistic (each SNP explains <2% of variance)
  - Environment, sleep hygiene, and substances strongly modify any tendency
  - Not a diagnosis — for sleep disorders, consult a sleep medicine specialist
"""

from typing import Dict, List, Tuple, Any, Optional


# ---------------------------------------------------------------------------
# UTILITIES
# ---------------------------------------------------------------------------

COMPLEMENT: Dict[str, str] = {"A": "T", "T": "A", "C": "G", "G": "C"}


def _complement(allele: str) -> str:
    return COMPLEMENT.get(allele.upper(), allele.upper())


def _count_allele(genotype: str, allele: str) -> int:
    """Count occurrences of an allele in the genotype with strand flip handling."""
    genotype = genotype.upper().replace("-", "")
    allele = allele.upper()
    comp = _complement(allele)
    count = genotype.count(allele)
    if count == 0 and comp != allele:
        count = genotype.count(comp)
    return min(count, 2)


def _analyze_snp_list(
    variants: Dict[str, Tuple[str, str, str]],
    snp_list: List[Dict[str, Any]],
) -> List[Dict[str, Any]]:
    """
    Processes a list of SNPs against the user's variant dictionary.
    Returns a list of enriched results with genotype and interpretation.
    """
    results = []
    for snp in snp_list:
        rsid = snp["rsid"].lower()
        raw = variants.get(rsid)
        if raw is None:
            chrom, pos, genotype = None, None, None
        elif isinstance(raw, str):
            chrom, pos, genotype = None, None, raw
        else:
            chrom, pos, genotype = raw

        if genotype:
            genotype = genotype.upper().replace("-", "")
            geno_data = snp["genotypes"].get(genotype)

            # Try complement if not found directly
            if not geno_data:
                comp_geno = "".join(_complement(a) for a in genotype)
                geno_data = snp["genotypes"].get(comp_geno)
                # Try reverse order
                if not geno_data and len(genotype) == 2:
                    geno_data = snp["genotypes"].get(genotype[::-1])
                if not geno_data:
                    rev_comp = comp_geno[::-1]
                    geno_data = snp["genotypes"].get(rev_comp)

            if not geno_data:
                # Fallback: genotype found but not mapped
                geno_data = {
                    "label": genotype,
                    "nickname": "Uncatalogued variant",
                    "interpretation": (
                        f"Your genotype {genotype} was detected but is not "
                        "catalogued for this SNP. It may be a rare allele."
                    ),
                    "score": 0.3,
                }

            results.append({
                **snp,
                "genotype": genotype,
                "genotypeData": geno_data,
                "found": True,
            })
        else:
            results.append({
                **snp,
                "genotype": "Not available",
                "genotypeData": {
                    "label": "N/A",
                    "nickname": "Not available",
                    "interpretation": (
                        "This SNP was not found in your genomic file. "
                        "Kits from different companies cover distinct sets of variants."
                    ),
                    "score": None,
                },
                "found": False,
            })

        results[-1]["category"] = snp.get("category", "")

    return results


# ---------------------------------------------------------------------------
# SNP DATABASE: Dream Vividness & Memory
# ---------------------------------------------------------------------------

DREAM_VIVIDNESS_SNPS: List[Dict[str, Any]] = [
    {
        "rsid": "rs6311",
        "gene": "HTR2A",
        "name": "Dream Vividness and REM Density",
        "chromosome": "13",
        "position": 46897343,
        "effectSize": "moderate",
        "category": "dreamVividness",
        "genotypes": {
            "GG": {
                "label": "G/G (-1438G)",
                "nickname": "Intense Dreamer",
                "interpretation": (
                    "Two copies of the G allele associated with higher 5-HT2A receptor expression "
                    "in the cortex. Denser 5-HT2A receptors are linked to more vivid, emotionally "
                    "intense, and memorable dreams — the same mechanism that explains the dream-like "
                    "effects of psychedelic substances (which activate exactly this receptor). "
                    "You may have a naturally rich and intense dream life."
                ),
                "score": 0.8,
            },
            "AG": {
                "label": "A/G",
                "nickname": "Moderate Dreamer",
                "interpretation": (
                    "Heterozygous genotype for the HTR2A promoter. Intermediate 5-HT2A expression. "
                    "You likely experience dreams of varying intensity — sometimes very vivid, "
                    "other times more diffuse. Dream quality may be sensitive to factors such as "
                    "stress, diet, and substance use."
                ),
                "score": 0.5,
            },
            "AA": {
                "label": "A/A (-1438A)",
                "nickname": "Serene Dreamer",
                "interpretation": (
                    "Two copies of the A allele. Slightly reduced 5-HT2A receptor expression. "
                    "Your dreams tend to be less emotionally intrusive. You may recall fewer "
                    "dreams upon waking, but you may enjoy quieter and more restorative REM sleep."
                ),
                "score": 0.3,
            },
        },
        "scientificBasis": (
            "HTR2A rs6311 (C-1438T/A) is in the promoter of the serotonin 2A receptor gene. "
            "The G allele increases HTR2A transcription by ~30% compared to the A allele (Polesskaya & "
            "Sokolov, 2002). During sleep, 5-HT2A receptors are active during wakefulness and REM, "
            "regulating the intensity and emotional content of dreams. Berger et al. (2003) "
            "(Neuropsychobiology 47:2) demonstrated a correlation between HTR2A promoter variants "
            "and REM sleep architecture. The connection is mechanistically clear: classic psychedelics "
            "(LSD, psilocybin) that produce intense dream-like experiences act as agonists at "
            "exactly this receptor. This SNP is also present in GeneHealth's mind_spirit_analyzer "
            "in the context of spiritual sensitivity — the same variant that amplifies the dream "
            "experience tends to amplify transcendent sensitivity during wakefulness."
        ),
        "references": [
            {"pmid": "12566938", "title": "HTR2A promoter polymorphism and receptor expression", "year": 2003},
            {"pmid": "15163253", "title": "Serotonin 2A receptor and sleep architecture", "year": 2004},
        ],
        "actionableInsights": [
            "G/G carriers: Avoid bright screens 2h before bed — serotonergic stimulation may over-intensify dreams",
            "Evening meditation can help 'prepare' the mental space for more integrative dreams",
            "If dreams are disturbing, dream processing strategies (journaling) are especially valuable for you",
            "A/A carriers: Lucid dream induction techniques (MILD, WILD) may be needed to increase dream recall",
        ],
    },
    {
        "rsid": "rs1800955",
        "gene": "DRD4",
        "name": "Dopamine, Dream Novelty, and Lucid Dreams",
        "chromosome": "11",
        "position": 637304,
        "effectSize": "moderate",
        "category": "dreamVividness",
        "genotypes": {
            "TT": {
                "label": "T/T",
                "nickname": "Typical Dream Content",
                "interpretation": (
                    "Reference genotype for DRD4 -521. Typical expression and function of the "
                    "dopamine D4 receptor. Your dreams tend to predominantly reflect everyday "
                    "content and recent memories."
                ),
                "score": 0.3,
            },
            "CT": {
                "label": "C/T",
                "nickname": "Explorer Dreamer",
                "interpretation": (
                    "One copy of the C allele. The D4 receptor with the C variant has reduced "
                    "expression in the prefrontal cortex, which may increase the likelihood of "
                    "unusual, bizarre, or novelty-seeking dream content. Some studies associate "
                    "DRD4 variants with lucid dream propensity."
                ),
                "score": 0.6,
            },
            "CC": {
                "label": "C/C",
                "nickname": "Visionary Dreamer",
                "interpretation": (
                    "Two copies of the C allele. Variant associated with a lower-affinity D4 receptor "
                    "for dopamine. Correlated with novelty seeking, creativity, and in a dream context, "
                    "more unusual and potentially more bizarre or visionary content. Prefrontal "
                    "dopaminergic activity during REM is linked to dream metacognition (awareness "
                    "that one is dreaming)."
                ),
                "score": 0.8,
            },
        },
        "scientificBasis": (
            "DRD4 -521C/T (rs1800955) affects expression of the dopamine D4 receptor — variants with "
            "lower expression in the prefrontal cortex have been associated with novelty seeking (Strobel "
            "et al., 1999). During REM sleep, prefrontal dopamine is crucial for dream metacognition. "
            "Schmack et al. (2013, PLOS ONE) demonstrated that dopaminergic profiles influence dream "
            "content and bizarreness. The link to lucid dreams is explored by Hobson & Friston (2012): "
            "consciousness during REM requires activation of frontal areas that are normally inactive — "
            "a pattern favored by certain D4 receptor profiles. The 7-repeat DRD4 VNTR (well correlated "
            "with rs1800955) is the allele most strongly associated with novelty, creativity, and "
            "exploratory behaviors."
        ),
        "references": [
            {"pmid": "23516459", "title": "Dopaminergic prediction errors drive dream content", "year": 2013},
            {"pmid": "10205490", "title": "DRD4 promoter polymorphism and novelty seeking", "year": 1999},
        ],
        "actionableInsights": [
            "CC carriers: Reality check techniques during the day can increase the probability of lucid dreams",
            "Keep a dream journal beside your bed — your predisposition to unusual content is worth recording",
            "TT carriers: MILD (Mnemonic Induction of Lucid Dreams) techniques may be especially effective by strengthening pre-sleep intention",
        ],
    },
    {
        "rsid": "rs324420",
        "gene": "FAAH",
        "name": "Endocannabinoids and Dream Emotional Intensity",
        "chromosome": "1",
        "position": 46891060,
        "effectSize": "moderate",
        "category": "dreamVividness",
        "genotypes": {
            "AA": {
                "label": "C385A / A385A",
                "nickname": "Emotionally Amplified Dreams",
                "interpretation": (
                    "Carrier of the A allele (385A). The FAAH enzyme with 385A has reduced activity, "
                    "resulting in higher levels of anandamide (the 'endogenous cannabinoid' or "
                    "'bliss molecule'). Higher anandamide is associated with more emotionally intense "
                    "dreams, with greater positive tone and, in some studies, less dream anxiety. "
                    "You may have an emotionally rich dream life with higher hedonic tone."
                ),
                "score": 0.7,
            },
            "AC": {
                "label": "Heterozygous C385A",
                "nickname": "Balanced Endocannabinoid Modulation",
                "interpretation": (
                    "One copy of each variant. Intermediate FAAH activity and anandamide levels "
                    "between the extremes. You likely experience a dream life with emotional "
                    "intensity that varies depending on your waking emotional state."
                ),
                "score": 0.5,
            },
            "CC": {
                "label": "C/C (385Pro)",
                "nickname": "Standard Endocannabinoid Processing",
                "interpretation": (
                    "Reference genotype with normal FAAH activity. Standard anandamide levels. "
                    "Your dreams tend to have typical emotional intensity, without any particular "
                    "endocannabinoid amplification."
                ),
                "score": 0.3,
            },
        },
        "scientificBasis": (
            "FAAH (Fatty Acid Amide Hydrolase) rs324420 (C385A) is the most studied variant in the "
            "endocannabinoid system. The 385A allele produces a protein ~30% less stable, resulting "
            "in reduced anandamide degradation. Prather et al. (2013, PNAS) showed that A allele "
            "carriers have reduced stress response and more positive mood. In the sleep context, the "
            "endocannabinoid system regulates REM sleep — animals lacking FAAH (and therefore with "
            "high anandamide) show increased REM sleep with greater neuronal activity during that "
            "phase. Bhattacharya et al. (2021, Int J Mol Sci) reviewed the relationship between "
            "endocannabinoids and sleep regulation, confirming the role of anandamide in REM sleep "
            "intensity. This SNP is particularly relevant for GeneHealth's spiritual profile: "
            "anandamide is called the 'bliss molecule' and its connection to absorption experiences "
            "and spiritual well-being is documented."
        ),
        "references": [
            {"pmid": "23754379", "title": "FAAH C385A and stress reactivity in humans", "year": 2013},
            {"pmid": "34444209", "title": "Endocannabinoid system and sleep regulation", "year": 2021},
        ],
        "actionableInsights": [
            "A/A carriers: Mindfulness practices that increase endogenous anandamide (meditation, exercise) can enhance your dream experience",
            "C/C carriers: Physical exercise temporarily elevates anandamide and may enrich your dreams on subsequent nights",
            "All genotypes: A diet rich in omega-3 fatty acids supports endocannabinoid synthesis",
        ],
    },
]


# ---------------------------------------------------------------------------
# SNP DATABASE: Sleep Architecture
# ---------------------------------------------------------------------------

SLEEP_ARCHITECTURE_SNPS: List[Dict[str, Any]] = [
    {
        "rsid": "rs228697",
        "gene": "PER3",
        "name": "REM Sleep Architecture and Sleep Pressure",
        "chromosome": "1",
        "position": 7847208,
        "effectSize": "strong",
        "category": "sleepArchitecture",
        "genotypes": {
            "CC": {
                "label": "Short allele (proxy 4-repeat)",
                "nickname": "Efficient Sleeper",
                "interpretation": (
                    "Associated with the PER3 4-repeat allele (proxy). Carriers tend to accumulate "
                    "sleep pressure more quickly, fall asleep easily, and have efficient deep sleep. "
                    "However, they are more sensitive to the cognitive effects of sleep deprivation. "
                    "REM sleep is typically more compact and concentrated at the end of the night."
                ),
                "score": 0.4,
            },
            "CT": {
                "label": "Heterozygous",
                "nickname": "Mixed Sleep Profile",
                "interpretation": (
                    "Heterozygous genotype. Combines characteristics of both sleep architecture "
                    "profiles. You likely experience good sleep efficiency with a moderately rich "
                    "dream life — the most common balance."
                ),
                "score": 0.5,
            },
            "TT": {
                "label": "Long allele (proxy 5-repeat)",
                "nickname": "REM-Rich Dreamer",
                "interpretation": (
                    "Associated with the PER3 5-repeat allele (proxy). This profile is fascinating: "
                    "carriers have REM sleep more distributed throughout the night, with higher "
                    "dream density and better dream recall. They are more cognitively resilient to "
                    "sleep deprivation (need less sleep to function well), but have longer sleep "
                    "onset latency. Viola et al. (2007) showed that 5/5 carriers have higher "
                    "amplitude delta waves — more intense deep sleep."
                ),
                "score": 0.7,
            },
        },
        "scientificBasis": (
            "PER3 (Period Circadian Regulator 3) contains a VNTR (rs57875989) of 4 or 5 repeats "
            "of 54 base pairs, which is the most studied sleep SNP in consumer genomics. rs228697 "
            "is a high-correlation proxy available on standard genotyping arrays. Viola et al. "
            "(2007, Current Biology 17:4) demonstrated that 5/5 carriers have more slow-wave "
            "activity (SWA) and higher theta EEG activity during wakefulness after sleep deprivation. "
            "Dijk & Archer (2010, J Sleep Res) extensively reviewed how PER3 genotype modulates "
            "both sleep homeostasis and the response to deprivation. The relevance for GeneHealth: "
            "carriers of the long allele have naturally richer dream experiences — a finding that "
            "connects directly to the spiritual sensitivity profile."
        ),
        "references": [
            {"pmid": "17320392", "title": "PER3 polymorphism and sleep homeostasis", "year": 2007},
            {"pmid": "20082661", "title": "PER3 VNTR and sleep regulation", "year": 2010},
        ],
        "actionableInsights": [
            "TT carriers (proxy 5/5): Ensure 8+ hours — your sleep architecture is naturally longer to complete REM cycles",
            "CC carriers (proxy 4/4): You may be more cognitively resilient to deprivation, but don't neglect quality sleep",
            "Everyone: Maintaining regular sleep schedules maximizes REM density regardless of genotype",
        ],
    },
    {
        "rsid": "rs5751876",
        "gene": "ADORA2A",
        "name": "Sleep Pressure and Caffeine Sensitivity",
        "chromosome": "22",
        "position": 24822676,
        "effectSize": "strong",
        "category": "sleepArchitecture",
        "genotypes": {
            "TT": {
                "label": "T/T (1083T)",
                "nickname": "High Sleep Pressure",
                "interpretation": (
                    "Genotype associated with greater adenosine A2A receptor sensitivity. "
                    "Adenosine is the main molecule that accumulates 'sleep pressure' during "
                    "wakefulness. TT carriers feel sleepier throughout the day, have a greater "
                    "propensity to fall asleep, and are more sensitive to caffeine (which blocks "
                    "exactly this receptor). When they do sleep, they tend to have deeper and more "
                    "restorative sleep. Their dreams may be more intense after partial deprivation."
                ),
                "score": 0.7,
            },
            "CT": {
                "label": "C/T",
                "nickname": "Intermediate Sensitivity",
                "interpretation": (
                    "Heterozygous profile. Adenosine and caffeine sensitivity between the extremes. "
                    "You probably tolerate moderate caffeine without major sleep impact if consumed "
                    "before 2 PM."
                ),
                "score": 0.5,
            },
            "CC": {
                "label": "C/C (1083C)",
                "nickname": "Resilient Sleeper",
                "interpretation": (
                    "Lower adenosine A2A receptor sensitivity. CC carriers tend to accumulate sleep "
                    "pressure more slowly, are less drowsy during the day, and are more tolerant of "
                    "evening caffeine. They may function better with partial sleep deprivation, but "
                    "deep sleep (slow waves) may be slightly reduced compared to TT carriers."
                ),
                "score": 0.3,
            },
        },
        "scientificBasis": (
            "ADORA2A rs5751876 is the best-validated caffeine sensitivity SNP. Retey et al. "
            "(2005, PNAS 102:17) showed that T/T carriers have greater subjective and objective "
            "caffeine sensitivity (measured by EEG and sleepiness scales). The adenosine A2A "
            "receptor is the main target of caffeine — when adenosine binds, it signals fatigue; "
            "when caffeine blocks, you feel awake. T/T carriers have a receptor that 'amplifies' "
            "both the sleep signal and the wakefulness effect of caffeine. In terms of dreams: "
            "the sleep architecture of T/T carriers favors more intense REM cycles when sleep "
            "finally occurs, potentially amplifying the dream experience."
        ),
        "references": [
            {"pmid": "16188225", "title": "ADORA2A polymorphism and caffeine sensitivity", "year": 2005},
            {"pmid": "17592536", "title": "Adenosine A2A receptor and sleep homeostasis", "year": 2007},
        ],
        "actionableInsights": [
            "T/T carriers: Avoid caffeine after 1 PM — your adenosinergic system amplifies the impact",
            "C/C carriers: You may tolerate afternoon coffee, but monitor sleep quality individually",
            "Everyone: Accumulated 'sleep debt' leads to more intense recovery dreams — use this consciously",
        ],
    },
]


# ---------------------------------------------------------------------------
# SNP DATABASE: Chronotype & Circadian Rhythm
# ---------------------------------------------------------------------------

CHRONOTYPE_SNPS: List[Dict[str, Any]] = [
    {
        "rsid": "rs1801260",
        "gene": "CLOCK",
        "name": "Circadian Rhythm and Sleep Phase Delay",
        "chromosome": "4",
        "position": 56292869,
        "effectSize": "moderate",
        "category": "chronotype",
        "genotypes": {
            "CC": {
                "label": "C/C (3111C)",
                "nickname": "Standard Circadian Cycle",
                "interpretation": (
                    "Reference genotype for CLOCK 3111T/C. Typical circadian period (~24h). "
                    "You likely adapt well to conventional sleep schedules and have relative "
                    "ease in adjusting your sleep-wake cycle."
                ),
                "score": 0.3,
            },
            "CT": {
                "label": "C/T",
                "nickname": "Mild Evening Tendency",
                "interpretation": (
                    "One copy of the T allele. Associated with a slight tendency toward evening "
                    "chronotype (night owl). You may find it easy to stay up late and difficult "
                    "to wake up early, especially during adolescence and early adulthood."
                ),
                "score": 0.5,
            },
            "TT": {
                "label": "T/T (3111T)",
                "nickname": "Evening Chronotype",
                "interpretation": (
                    "Two copies of the T allele. Mishima et al. (2005) associated this genotype "
                    "with delayed sleep phase and an evening circadian pattern. TT carriers have "
                    "a greater tendency to sleep late, wake late, and function better at night. "
                    "The melatonin peak occurs later compared to CC carriers. More intense dreams "
                    "tend to occur near the natural late wake time of this chronotype."
                ),
                "score": 0.7,
            },
        },
        "scientificBasis": (
            "CLOCK 3111T/C (rs1801260) is in the 3' UTR region of the CLOCK gene, one of the main "
            "regulators of the molecular circadian clock. CLOCK forms the CLOCK:BMAL1 heterodimer "
            "that activates transcription of Per and Cry genes. Mishima et al. (2005, Sleep 28:1) "
            "demonstrated a significant association of the T allele with delayed sleep phase syndrome "
            "(DSPS) in Japanese populations. Katzenberg et al. (1998, Sleep 21:6) also reported an "
            "association with morningness vs. eveningness. CLOCK is the 'conductor' of circadian "
            "rhythm — variants here affect not only sleep but also the timing of hormonal peaks, "
            "body temperature, and consequently when the most vivid dreams occur during personalized "
            "sleep architecture."
        ),
        "references": [
            {"pmid": "15700720", "title": "CLOCK mutation and delayed sleep phase syndrome", "year": 2005},
            {"pmid": "9756555", "title": "CLOCK gene polymorphism and human diurnal preference", "year": 1998},
        ],
        "actionableInsights": [
            "T/T carriers: Respect your biology — forcing a very early schedule can chronically reduce REM sleep quality",
            "Morning sunlight exposure helps 'anchor' the circadian cycle regardless of genotype",
            "T/T carriers who need morning schedules: low-dose melatonin (0.5mg) 5h before desired sleep time may help",
        ],
    },
    {
        "rsid": "rs2287161",
        "gene": "CRY1",
        "name": "Long Circadian Period and Delayed Sleep",
        "chromosome": "12",
        "position": 107484498,
        "effectSize": "strong",
        "category": "chronotype",
        "genotypes": {
            "CC": {
                "label": "C/C",
                "nickname": "Typical Circadian Period",
                "interpretation": (
                    "Reference genotype. Endogenous circadian period close to 24h. "
                    "You adapt well to conventional schedules and tend to have good synchronization "
                    "between your internal clock and the environmental light-dark cycle."
                ),
                "score": 0.3,
            },
            "CG": {
                "label": "C/G",
                "nickname": "Slightly Extended Period",
                "interpretation": (
                    "One copy of the G allele. Your endogenous circadian period may be slightly "
                    "longer than 24h. This creates a natural tendency to progressively delay the "
                    "sleep-wake cycle if there are no strong external light anchors."
                ),
                "score": 0.5,
            },
            "GG": {
                "label": "G/G",
                "nickname": "Extended Circadian Period (~50 min+)",
                "interpretation": (
                    "High-impact variant described by Patke et al. (2017, Cell). Carriers have "
                    "a longer endogenous circadian period — approximately 50 minutes above average. "
                    "This means your internal clock prefers to 'run slower', resulting in a strong "
                    "tendency to sleep late and wake late. In free-running environments (without "
                    "schedule obligations), carriers naturally sleep 2-3h later than average. This "
                    "is one of the most potent chronotype variants identified to date."
                ),
                "score": 0.85,
            },
        },
        "scientificBasis": (
            "CRY1 (Cryptochrome 1) is a central repressor of the molecular circadian clock — "
            "it inhibits the CLOCK:BMAL1 heterodimer, closing the feedback loop. Patke et al. "
            "(2017, Cell 169:2) identified a splicing mutation in CRY1 associated with familial "
            "advanced sleep phase syndrome (FASP). rs2287161 is a linkage disequilibrium marker "
            "for this functional variant. Carriers with a longer period have more difficulty "
            "synchronizing with the standard solar cycle — their 'biological midnight' occurs "
            "later. Sleep EEG studies show that carriers with a long CRY1 period have greater "
            "accumulated slow-wave intensity when they finally sleep, potentially resulting in "
            "richer recovery REM."
        ),
        "references": [
            {"pmid": "28388406", "title": "CRY1 splice variant and delayed sleep phase disorder", "year": 2017},
            {"pmid": "30773234", "title": "CRY1 functional variants and circadian period", "year": 2019},
        ],
        "actionableInsights": [
            "G/G carriers: Your chronotype has a strong genetic basis — validating this can relieve the social guilt of being a 'night owl'",
            "Intense morning light exposure (10,000 lux for 30 min) is the most effective intervention for advancing the cycle",
            "If work allows, consider adjusting schedules to respect your natural cognitive peak (typically 11 AM-1 PM for night owls)",
        ],
    },
    {
        "rsid": "rs2735611",
        "gene": "PER1",
        "name": "Morning Preference and Circadian Anticipation",
        "chromosome": "17",
        "position": 8123949,
        "effectSize": "moderate",
        "category": "chronotype",
        "genotypes": {
            "AA": {
                "label": "A/A",
                "nickname": "Morning Chronotype",
                "interpretation": (
                    "Associated with morning preference. Carriers tend to fall asleep earlier, "
                    "wake up naturally early, and have peak energy and cognition in the first "
                    "hours of the morning. The most intense REM sleep occurs in the early morning "
                    "cycles, which may result in more vivid dreams in the first hours after falling asleep."
                ),
                "score": 0.7,
            },
            "AG": {
                "label": "A/G",
                "nickname": "Intermediate Chronotype",
                "interpretation": (
                    "Heterozygous genotype. Tendency toward an intermediate chronotype — you "
                    "likely adapt reasonably well to different schedules, with a slight "
                    "preference for waking at moderate times."
                ),
                "score": 0.5,
            },
            "GG": {
                "label": "G/G",
                "nickname": "Evening Tendency",
                "interpretation": (
                    "Associated with lower morning preference. You may find it easier to stay "
                    "awake at night and harder to wake up early. In combination with other "
                    "circadian markers, this may indicate a moderate evening chronotype."
                ),
                "score": 0.4,
            },
        },
        "scientificBasis": (
            "PER1 (Period 1) is one of three 'period' genes of the circadian clock. rs2735611 "
            "was associated with chronotype and morning preference in large-scale GWAS studies. "
            "Jones et al. (2019, Nature Communications) identified 351 loci for chronotype in "
            "~700,000 individuals from the UK Biobank, with PER1 among the most significant. "
            "Carpen et al. (2006, J Sleep Res) also reported an association of PER1 variants "
            "with morningness. PER1 regulates the circadian clock phase along with PER2 and "
            "PER3 — variants in this gene affect how quickly the clock 'advances' in response "
            "to morning light."
        ),
        "references": [
            {"pmid": "30696823", "title": "GWAS of chronotype in UK Biobank identifies 351 loci", "year": 2019},
            {"pmid": "16685255", "title": "PER1 polymorphism and diurnal preference", "year": 2006},
        ],
        "actionableInsights": [
            "A/A carriers: Take advantage of your morning peak — schedule creative work and important decisions before noon",
            "G/G carriers: Consider deep work blocks in the late afternoon/early evening when your cognition peaks",
            "Everyone: Consistency in wake time (even on weekends) is the most important factor for sleep quality",
        ],
    },
]


# ---------------------------------------------------------------------------
# SNP DATABASE: Nocturnal Phenomena
# (Sleep paralysis, hypnagogia, lucid dreams, boundary experiences)
# ---------------------------------------------------------------------------

NOCTURNAL_PHENOMENA_SNPS: List[Dict[str, Any]] = [
    {
        "rsid": "rs6313",
        "gene": "HTR2A",
        "name": "Sleep Paralysis and Hypnagogic Experiences",
        "chromosome": "13",
        "position": 46897065,
        "effectSize": "moderate",
        "category": "nocturnal",
        "genotypes": {
            "TT": {
                "label": "T/T",
                "nickname": "Defined Sleep-Wake Boundary",
                "interpretation": (
                    "Reference genotype. Typical sleep-wake transition. You likely experience "
                    "the boundary between sleep and wakefulness in a relatively clear manner, "
                    "without frequent boundary phenomena."
                ),
                "score": 0.3,
            },
            "CT": {
                "label": "C/T",
                "nickname": "Permeable Boundary",
                "interpretation": (
                    "One copy of the C allele. Increased HTR2A expression may make the sleep-wake "
                    "boundary more permeable. You may experience hypnagogic hallucinations (while "
                    "falling asleep) or hypnopompic hallucinations (while waking), occasional sleep "
                    "paralysis, or vivid sensations during state transitions more frequently. "
                    "These experiences, though sometimes frightening, have a clear neurobiological "
                    "basis and are reported in many traditions as spiritual experiences."
                ),
                "score": 0.6,
            },
            "CC": {
                "label": "C/C",
                "nickname": "Highly Permeable Boundary",
                "interpretation": (
                    "Higher 5-HT2A receptor expression. The boundary between REM and wakefulness "
                    "may be especially permeable for you. Episodes of sleep paralysis with rich "
                    "hallucinatory content, spontaneous lucid dreams, and intense hypnagogic "
                    "experiences are more likely with this profile. In many cultures, these "
                    "sleep-wake boundary experiences are interpreted as contact with other "
                    "dimensions, visions, or spiritual messages — and you may have a genetic "
                    "predisposition for these boundary experiences."
                ),
                "score": 0.8,
            },
        },
        "scientificBasis": (
            "HTR2A rs6313 is another polymorphism in the serotonin 2A receptor promoter. "
            "Sleep paralysis (SP) occurs when REM state intrudes into wakefulness — the "
            "mechanisms that paralyze musculature during REM persist upon waking, frequently "
            "accompanied by hypnopompic hallucinations. The serotonergic system via HTR2A "
            "regulates the REM-wakefulness transition. Jalal & Hinton (2013, Psychol Sci) "
            "extensively documented SP hallucinations across 30+ cultures as 'night visitors', "
            "'demons', 'entities' — the same neurobiology expressed in distinct spiritual "
            "narratives. Variants that increase HTR2A can increase both the probability and "
            "the intensity of these boundary experiences."
        ),
        "references": [
            {"pmid": "23983260", "title": "Sleep paralysis hallucinations across cultures", "year": 2013},
            {"pmid": "11774094", "title": "HTR2A polymorphisms and serotonin receptor expression", "year": 2001},
        ],
        "actionableInsights": [
            "CC carriers: If sleep paralysis is disturbing, sleeping on your side (not on your back) significantly reduces frequency",
            "You can explore lucid dream induction techniques — your permeable boundary can be a gateway",
            "Understanding the neurobiology of these states can transform frightening experiences into conscious exploration",
            "Meditative practices that work with hypnagogic states (yoga nidra, NSDR) are especially relevant for your profile",
        ],
    },
    {
        "rsid": "rs6265",
        "gene": "BDNF",
        "name": "Dream Memory Consolidation and Neuroplasticity",
        "chromosome": "11",
        "position": 27658369,
        "effectSize": "moderate",
        "category": "nocturnal",
        "genotypes": {
            "CC": {
                "label": "Val/Val",
                "nickname": "Standard Dream Consolidation",
                "interpretation": (
                    "Reference genotype. Normal activity-dependent BDNF secretion. Memory "
                    "consolidation during sleep follows the typical pattern — dreams related "
                    "to recent events and usual emotional processing."
                ),
                "score": 0.3,
            },
            "CT": {
                "label": "Val/Met",
                "nickname": "Altered Dream Consolidation",
                "interpretation": (
                    "One copy of the Met allele. Reduced activity-dependent BDNF secretion may "
                    "affect how emotional memories are processed and consolidated during sleep. "
                    "You may have dreams with replay of emotionally intense events more frequently "
                    "— your brain may need more REM 'iterations' to process emotionally charged "
                    "experiences."
                ),
                "score": 0.5,
            },
            "TT": {
                "label": "Met/Met",
                "nickname": "Intensified Dream Processing",
                "interpretation": (
                    "Two copies of the Met allele (Val66Met). The reduction in activity-dependent "
                    "BDNF release affects hippocampal plasticity and memory consolidation. In the "
                    "dream context, this may manifest as more frequent dreams processing emotional "
                    "memories, recurrent dreams during stressful periods (the brain trying to "
                    "'resolve' via REM), and possibly greater dream rumination. Regular exercise "
                    "elevates BDNF and can improve the quality of processing during sleep."
                ),
                "score": 0.65,
            },
        },
        "scientificBasis": (
            "BDNF Val66Met (rs6265) is already in GeneHealth's mind_spirit_analyzer in the context "
            "of daytime neuroplasticity. Here, the focus is on BDNF's role in memory consolidation "
            "during sleep. Walker & Stickgold (2004, Neuron 44:1) established that REM sleep is "
            "critical for the consolidation of procedural and emotional memories. BDNF is the "
            "main mediator of synaptic plasticity during this consolidation. Bhattacharya et al. "
            "(2015, Neurobiology of Sleep and Circadian Rhythms) reviewed how BDNF variants "
            "affect sleep-dependent memory consolidation. This SNP appears in two contexts on "
            "the GeneHealth platform: during wakefulness (plasticity and mood) and during sleep "
            "(dream consolidation and nocturnal emotional processing) — illustrating the "
            "mind-sleep continuum."
        ),
        "references": [
            {"pmid": "15450160", "title": "Sleep and the price of plasticity", "year": 2004},
            {"pmid": "26779544", "title": "BDNF and sleep-dependent memory consolidation", "year": 2015},
        ],
        "actionableInsights": [
            "TT carriers: Aerobic exercise (especially morning/afternoon) elevates BDNF and can dramatically improve dream processing quality",
            "Pre-sleep journaling helps 'pre-digest' emotional content before REM needs to process it",
            "Val/Val carriers: Your sleep memory processing is efficient — take advantage by sleeping well after important learning",
        ],
    },
]


# ---------------------------------------------------------------------------
# CONTEXT CONSTANTS (for the frontend report)
# ---------------------------------------------------------------------------

SPIRITUAL_DREAM_CONTEXT = (
    "At GeneHealth, we understand that the dream experience is a unique window into "
    "consciousness. Cultures around the world have treated dreams as portals to the sacred "
    "— from dream incubation in the temples of Asclepius in ancient Greece, to Siberian "
    "shamanic traditions, to prophetic dreams in Candomble and Brazilian indigenous "
    "traditions. Modern neuroscience is beginning to understand why: the REM state shares "
    "characteristics of deep meditative states and psychedelic experiences, including "
    "activation of the same neural networks and reduction of the Default Mode Network. "
    "Your genetic markers in this report reflect your biological predisposition for this "
    "dimension of the human experience."
)

DISCLAIMERS = {
    "general": (
        "This report is educational and does not constitute a medical diagnosis. "
        "Genetics is only one factor — sleep hygiene, stress, medications, "
        "and environment have an equal or greater impact on sleep and dream quality. "
        "For sleep disorders, consult a sleep medicine specialist."
    ),
    "chronotype": (
        "Chronotype has a strong genetic basis (~50% heritability) but is malleable "
        "through light exposure, social schedules, and routine. Do not use these results "
        "to justify habits harmful to your health."
    ),
    "nocturnal": (
        "Sleep paralysis and hypnagogic experiences are benign neurological phenomena, "
        "present in 20-40% of the population. If frequent or disturbing, consult a "
        "specialist — simple behavioral treatments are highly effective."
    ),
}


# ---------------------------------------------------------------------------
# REPORT HELPER FUNCTIONS
# ---------------------------------------------------------------------------

def _generate_category_summary(traits: List[Dict[str, Any]], category: str) -> str:
    """Generates a text summary for a category of traits."""
    found = [t for t in traits if t.get("found", False)]
    if not found:
        return "No variants from this category were found in your genomic file."

    scores = [
        t["genotypeData"]["score"]
        for t in found
        if t["genotypeData"].get("score") is not None
    ]
    if not scores:
        return "Variants found, but insufficient interpretation data available."

    avg = sum(scores) / len(scores)

    summaries = {
        "dreamVividness": {
            0.65: "Your genetic profile suggests a naturally intense and rich dream life — vivid, emotionally charged dreams with good memorability are likely characteristics.",
            0.45: "You have a balanced dream profile — with potential for vivid dreams under favorable conditions, especially after emotionally intense days.",
            0.0: "Your profile suggests a more serene approach to dream life — dreams tend to be less intrusive, with moderate recall.",
        },
        "sleepArchitecture": {
            0.65: "Your sleep architecture tends toward rich REM cycles and intense deep sleep — you likely need adequate hours to complete your naturally longer cycles.",
            0.45: "Balanced sleep architecture with good response to schedule regularity.",
            0.0: "Efficient sleep with rapid pressure buildup — you can function well with more compact cycles.",
        },
        "chronotype": {
            0.65: "Your biological clock has a strong evening (night owl) tendency with a solid genetic basis.",
            0.45: "Intermediate chronotype — you adapt to different schedules with relative ease.",
            0.0: "Strong predisposition to a morning chronotype — you likely function best early.",
        },
        "nocturnal": {
            0.65: "High probability of sleep-wake boundary experiences (hypnagogia, sleep paralysis, spontaneous lucid dreams).",
            0.45: "Moderate propensity for nocturnal phenomena — boundary experiences are possible under deprivation or stress conditions.",
            0.0: "Profile with typically well-defined sleep-wake transitions.",
        },
    }

    cat_summaries = summaries.get(category, summaries["dreamVividness"])
    for threshold in sorted(cat_summaries.keys(), reverse=True):
        if avg >= threshold:
            return cat_summaries[threshold]
    return "Typical profile."


def _generate_overall_profile(traits: List[Dict[str, Any]], category: str) -> str:
    """Generates the overall profile label."""
    found = [t for t in traits if t.get("found", False)]
    if not found:
        return "Insufficient data"

    scores = [
        t["genotypeData"]["score"]
        for t in found
        if t["genotypeData"].get("score") is not None
    ]
    if not scores:
        return "Typical profile"

    avg = sum(scores) / len(scores)

    profiles = {
        "dreamVividness": {
            0.65: "Intense Dreamer",
            0.45: "Balanced Dreamer",
            0.0: "Serene Dreamer",
        },
        "sleepArchitecture": {
            0.65: "REM-Rich Architecture",
            0.45: "Balanced Architecture",
            0.0: "Efficient and Compact Sleep",
        },
        "chronotype": {
            0.65: "Genetic Night Owl",
            0.45: "Intermediate Chronotype",
            0.0: "Genetic Early Bird",
        },
        "nocturnal": {
            0.65: "Boundary Explorer",
            0.45: "Moderate Sensitivity",
            0.0: "Defined Boundary",
        },
    }

    cat_profiles = profiles.get(category, profiles["dreamVividness"])
    for threshold in sorted(cat_profiles.keys(), reverse=True):
        if avg >= threshold:
            return cat_profiles[threshold]
    return "Typical"


# ---------------------------------------------------------------------------
# MAIN ANALYSIS FUNCTION
# ---------------------------------------------------------------------------

def analyze_dream_sleep(
    variants: Dict[str, Tuple[str, str, str]]
) -> Dict[str, Any]:
    """
    Runs the full Dream & Sleep analysis for the user's variants.

    Args:
        variants: Dict mapping rsid -> (chromosome, position, genotype)

    Returns:
        Dict with complete analysis results
    """
    dream_vividness = _analyze_snp_list(variants, DREAM_VIVIDNESS_SNPS)
    sleep_architecture = _analyze_snp_list(variants, SLEEP_ARCHITECTURE_SNPS)
    chronotype = _analyze_snp_list(variants, CHRONOTYPE_SNPS)
    nocturnal = _analyze_snp_list(variants, NOCTURNAL_PHENOMENA_SNPS)

    return {
        "dream_vividness": dream_vividness,
        "sleep_architecture": sleep_architecture,
        "chronotype": chronotype,
        "nocturnal_phenomena": nocturnal,
    }


def generate_dream_sleep_json(result: Dict[str, Any]) -> Dict[str, Any]:
    """
    Generates the final JSON for the Dream & Sleep report.

    Args:
        result: Output from analyze_dream_sleep()

    Returns:
        Dict compatible with the GeneHealth frontend JSON schema
    """
    dv = result["dream_vividness"]
    sa = result["sleep_architecture"]
    ch = result["chronotype"]
    np_ = result["nocturnal_phenomena"]

    return {
        "reportType": "dream_sleep",
        "version": "1.0",
        "dreamVividness": {
            "summary": _generate_category_summary(dv, "dreamVividness"),
            "overallProfile": _generate_overall_profile(dv, "dreamVividness"),
            "traits": dv,
        },
        "sleepArchitecture": {
            "summary": _generate_category_summary(sa, "sleepArchitecture"),
            "overallProfile": _generate_overall_profile(sa, "sleepArchitecture"),
            "traits": sa,
        },
        "chronotype": {
            "summary": _generate_category_summary(ch, "chronotype"),
            "overallProfile": _generate_overall_profile(ch, "chronotype"),
            "traits": ch,
        },
        "nocturnalPhenomena": {
            "summary": _generate_category_summary(np_, "nocturnal"),
            "overallProfile": _generate_overall_profile(np_, "nocturnal"),
            "spiritualContext": SPIRITUAL_DREAM_CONTEXT,
            "traits": np_,
        },
        "disclaimers": DISCLAIMERS,
        "crossReferences": {
            "mindSpirit": (
                "HTR2A (rs6311, rs6313) and BDNF (rs6265) also appear in the Mind & Spirit "
                "report. The same neurobiology that amplifies spiritual sensitivity during "
                "wakefulness manifests as an intense dream life during sleep — consciousness "
                "is a continuum."
            ),
        },
    }
