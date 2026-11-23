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
  - HTR2A & hallucinogens: Nichols (2004) Pharmacol Ther — 5-HT2A receptors and hallucinogens
  - HTR2A rs6311: Parsons et al. (2004) Biol Psychiatry — promoter activity (dream link exploratory)
  - DRD4 rs1800955: Okuyama et al. (2000) Mol Psychiatry — novelty seeking (dream link exploratory)
  - FAAH rs324420: Sipe et al. (2002) PNAS; Dincheva et al. (2015) Nat Commun (dream link exploratory)
  - PER3 rs228697: Hida et al. (2014) Sci Rep — minor allele and evening preference; PER3 VNTR: Viola et al. (2007)
  - ADORA2A rs5751876: Rétey et al. (2007) Clin Pharmacol Ther; Alsene et al. (2003) Neuropsychopharmacology
  - CLOCK rs1801260: Katzenberg et al. (1998) Sleep; Mishima et al. (2005) Am J Med Genet B — 3111C allele and evening preference
  - CRY1 rs2287161: Patke et al. (2017) Cell — rare CRY1 variant in familial DSPD; rs2287161 itself has no established chronotype effect
  - PER1 rs2735611: Jones et al. (2019) Nat Commun; Carpen et al. (2006) J Hum Genet — direction for this SNP not established

SPIRITUAL CONNECTION (GeneHealth identity):
  - HTR2A is at the center of both vivid dreams and spiritual sensitivity
  - Sleep paralysis and hypnagogia are frequently reported as spiritual experiences
    across cultures (visions, sense of presence, astral travel)
  - Links between these variants and dream content are hypotheses; none has been tested directly

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
        "name": "Serotonin 2A Receptor Expression and Dreams (exploratory)",
        "chromosome": "13",
        "position": 46897343,
        "effectSize": "moderate",
        "category": "dreamVividness",
        "genotypes": {
            "AA": {
                "label": "A/A (-1438A)",
                "nickname": "Higher expected expression",
                "interpretation": (
                    "Two copies of the -1438A allele, which showed higher promoter activity "
                    "than the G allele in laboratory cell experiments and travels with the 102T "
                    "allele of rs6313, which is expressed more in human cortex. This points to "
                    "somewhat higher 5-HT2A receptor expression. If 5-HT2A receptors do shape "
                    "dream experience, this is the genotype where an effect would be most "
                    "expected, but no study has yet measured dream vividness by HTR2A genotype, "
                    "so this link is a hypothesis, not an established finding."
                ),
                "score": 0.8,
            },
            "AG": {
                "label": "A/G",
                "nickname": "Intermediate expected expression",
                "interpretation": (
                    "One copy of each allele. Intermediate promoter activity is expected. "
                    "Whether this has any bearing on your dreams has not been tested; dream "
                    "quality is strongly shaped by stress, diet, sleep timing and substance use."
                ),
                "score": 0.5,
            },
            "GG": {
                "label": "G/G (-1438G)",
                "nickname": "Lower expected expression",
                "interpretation": (
                    "Two copies of the -1438G allele, which showed lower promoter activity than "
                    "the A allele in laboratory cell experiments, pointing to somewhat lower "
                    "5-HT2A receptor expression. Whether this relates to dream experience has "
                    "not been tested; it remains a hypothesis."
                ),
                "score": 0.3,
            },
        },
        "scientificBasis": (
            "HTR2A rs6311 (-1438A/G) lies in the promoter of the serotonin 2A receptor gene. "
            "Parsons et al. (2004, Biol Psychiatry) found that the -1438A allele drove significantly "
            "greater promoter activity than G in cell lines that express HTR2A. rs6311 is in strong "
            "linkage disequilibrium with rs6313 (102T/C), and -1438A travels with 102T. In human "
            "temporal cortex, Polesskaya & Sokolov (2002, J Neurosci Res) found the 102C allele was "
            "expressed about 20% less than 102T, and C/C individuals had lower 5-HT2A mRNA and "
            "protein than T/T. Taken together, A (and T) corresponds to higher expression. The link "
            "to dreams is a hypothesis: classic psychedelics (LSD, psilocybin) act as agonists at "
            "this receptor and produce dream-like experiences, which motivates the idea, but no "
            "study has measured dream vividness by HTR2A genotype. It is also a hypothesis, not a "
            "finding, that the same variant could relate to transcendent experiences during "
            "wakefulness. DNA files report this SNP on the forward strand, where -1438A appears as "
            "T and -1438G as C. The score reflects expected receptor expression, not dream intensity."
        ),
        "references": [
            {"pmid": "15364038", "title": "The -1438A/G polymorphism in the 5-hydroxytryptamine type 2A receptor gene affects promoter activity.", "year": 2004},
            {"pmid": "11891796", "title": "Differential expression of the \"C\" and \"T\" alleles of the 5-HT2A receptor gene in the temporal cortex of normal individuals and schizophrenics.", "year": 2002},
        ],
        "actionableInsights": [
            "A/A carriers: Limiting bright screens in the 2 hours before bed is good sleep hygiene for everyone",
            "Evening relaxation or meditation can help prepare for sleep",
            "If dreams are disturbing, journaling can help you process them",
            "Keeping a dream journal beside your bed can help with dream recall, whatever your genotype",
        ],
    },
    {
        "rsid": "rs1800955",
        "gene": "DRD4",
        "name": "Dopamine Signalling and Novelty Seeking (dream link exploratory)",
        "chromosome": "11",
        "position": 637304,
        "effectSize": "moderate",
        "category": "dreamVividness",
        "genotypes": {
            "TT": {
                "label": "T/T",
                "nickname": "Lower Novelty-Seeking Tendency",
                "interpretation": (
                    "Two copies of the T allele, which reduces transcriptional efficiency of DRD4. "
                    "In the original study, T/T individuals had the lowest novelty-seeking scores. "
                    "No study has measured dream content or lucid dreaming by DRD4 genotype, so any "
                    "link to dreams is a hypothesis."
                ),
                "score": 0.3,
            },
            "CT": {
                "label": "C/T",
                "nickname": "Intermediate Novelty-Seeking",
                "interpretation": (
                    "One copy of the C allele. Novelty-seeking scores in the original study were "
                    "intermediate between the two homozygous groups. Whether this has any bearing "
                    "on dream content has not been tested."
                ),
                "score": 0.6,
            },
            "CC": {
                "label": "C/C",
                "nickname": "Higher Novelty-Seeking Tendency",
                "interpretation": (
                    "Two copies of the C allele, associated with higher DRD4 transcriptional "
                    "efficiency. In the original study, C/C individuals had the highest "
                    "novelty-seeking scores. Whether this relates to dream content, such as more "
                    "unusual dreams or lucid dreaming, has not been studied; it is a hypothesis "
                    "built on the novelty-seeking finding."
                ),
                "score": 0.8,
            },
        },
        "scientificBasis": (
            "DRD4 -521C/T (rs1800955) lies in the promoter of the dopamine D4 receptor gene. "
            "Okuyama et al. (2000, Mol Psychiatry) found it significantly associated with Novelty "
            "Seeking (P = 0.0001): C/C individuals scored highest and T/T lowest, and the T variant "
            "reduces transcriptional efficiency. This SNP is distinct from the DRD4 exon III "
            "repeat polymorphism. No study has measured dream content or lucid dreaming by DRD4 "
            "genotype, so any connection to dreams is exploratory."
        ),
        "references": [
            {"pmid": "10673770", "title": "Identification of a polymorphism in the promoter region of DRD4 associated with the human novelty seeking personality trait.", "year": 2000},
        ],
        "actionableInsights": [
            "Keep a dream journal beside your bed to notice your own patterns",
            "Regular sleep and wake times support dream recall and sleep quality for everyone",
            "If you are curious about lucid dreaming, techniques such as reality checks and MILD (Mnemonic Induction of Lucid Dreams) are widely discussed; individual results vary",
        ],
    },
    {
        "rsid": "rs324420",
        "gene": "FAAH",
        "name": "Endocannabinoid Enzyme FAAH (dream link exploratory)",
        "chromosome": "1",
        "position": 46891060,
        "effectSize": "moderate",
        "category": "dreamVividness",
        "genotypes": {
            "AA": {
                "label": "C385A / A385A",
                "nickname": "Reduced FAAH Expression",
                "interpretation": (
                    "Two copies of the A allele (385A, Pro129Thr). This variant gives an FAAH "
                    "enzyme with normal catalytic properties but enhanced sensitivity to "
                    "degradation, and in a mouse and human study it was linked to reduced FAAH "
                    "expression, stronger fronto-amygdala connectivity, better fear-extinction "
                    "learning and lower anxiety-like behaviour. No study has measured dream "
                    "intensity or emotional tone by FAAH genotype, so any link to your dreams is "
                    "exploratory."
                ),
                "score": 0.7,
            },
            "AC": {
                "label": "Heterozygous C385A",
                "nickname": "Balanced Endocannabinoid Modulation",
                "interpretation": (
                    "One copy of each variant. FAAH function is expected to be intermediate. "
                    "Whether this has any bearing on dreams has not been tested."
                ),
                "score": 0.5,
            },
            "CC": {
                "label": "C/C (385Pro)",
                "nickname": "Standard Endocannabinoid Processing",
                "interpretation": (
                    "Reference genotype with the common form of the FAAH enzyme. No particular "
                    "effect on dreams is expected or has been shown."
                ),
                "score": 0.3,
            },
        },
        "scientificBasis": (
            "FAAH (Fatty Acid Amide Hydrolase) rs324420 (C385A) is a well-studied variant in the "
            "endocannabinoid system. Sipe et al. (2002, PNAS) showed that the 385A (Pro129Thr) "
            "variant gives an FAAH with normal catalytic properties but enhanced sensitivity to "
            "proteolytic degradation. Dincheva et al. (2015, Nat Commun) found, in a knock-in mouse "
            "and in humans, that the variant allele reduced FAAH expression, enhanced fronto-amygdala "
            "connectivity and fear-extinction learning, and decreased anxiety-like behaviours. "
            "Huitron-Resendiz et al. (2004, Sleep) found that mice lacking FAAH had more and more "
            "intense slow-wave sleep at baseline. No study has measured dream intensity by FAAH "
            "genotype, so any link to dreams is exploratory."
        ),
        "references": [
            {"pmid": "12060782", "title": "A missense mutation in human fatty acid amide hydrolase associated with problem drug use.", "year": 2002},
            {"pmid": "25731744", "title": "FAAH genetic variation enhances fronto-amygdala function in mouse and human.", "year": 2015},
            {"pmid": "15453543", "title": "Characterization of the sleep-wake patterns in mice lacking fatty acid amide hydrolase.", "year": 2004},
        ],
        "actionableInsights": [
            "Regular physical exercise supports sleep quality for everyone",
            "Consistent sleep and wake times support restorative sleep whatever your genotype",
            "Keeping a dream journal can help you notice your own dream patterns",
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
        "name": "PER3 and Sleep Timing Preference",
        "chromosome": "1",
        "position": 7847208,
        "effectSize": "strong",
        "category": "sleepArchitecture",
        "genotypes": {
            "CC": {
                "label": "C/C",
                "nickname": "Typical Timing",
                "interpretation": (
                    "Two copies of the common C allele. The minor G allele of this marker was "
                    "more prevalent in evening types, so C/C is the genotype least associated "
                    "with eveningness in that study. This single marker explains only a small "
                    "part of sleep-timing preference."
                ),
                "score": 0.4,
            },
            "CG": {
                "label": "C/G",
                "nickname": "Mixed Profile",
                "interpretation": (
                    "One copy of the minor G allele. Intermediate between the two homozygous "
                    "groups for evening preference. Your actual sleep timing depends on many "
                    "genes, light exposure and daily routine."
                ),
                "score": 0.5,
            },
            "GG": {
                "label": "G/G",
                "nickname": "Evening-Leaning",
                "interpretation": (
                    "Two copies of the minor G allele, which was more common in evening types "
                    "than in morning types (sex-adjusted odds ratio 2.48) and in people with "
                    "free-running circadian rhythm disorder (odds ratio 2.02) in Hida et al. "
                    "(2014). This points to a lean toward later sleep timing, not a certainty."
                ),
                "score": 0.7,
            },
        },
        "scientificBasis": (
            "PER3 (Period Circadian Regulator 3) is a circadian clock gene. Hida et al. (2014, "
            "Sci Rep) screened clock gene polymorphisms and found that the minor allele of rs228697 "
            "(G) was more prevalent in evening types than morning types (sex-adjusted OR 2.48) and "
            "in people with free-running circadian rhythm disorder (OR 2.02). A separate, much-studied "
            "PER3 length polymorphism (VNTR) is linked to sleep structure and sleep-loss response "
            "(Viola et al., 2007, Curr Biol: 5/5 carriers had more slow-wave activity and larger "
            "cognitive deficits after sleep loss) and to delayed sleep phase syndrome and extreme "
            "diurnal preference (Archer et al., 2003, Sleep). That VNTR is not on standard "
            "genotyping arrays, and rs228697 is not an established proxy for it, so those findings "
            "are listed for context only."
        ),
        "references": [
            {"pmid": "25201053", "title": "Screening of clock gene polymorphisms demonstrates association of a PER3 polymorphism with morningness-eveningness preference and circadian rhythm sleep disorder.", "year": 2014},
            {"pmid": "17346965", "title": "PER3 polymorphism predicts sleep structure and waking performance.", "year": 2007},
            {"pmid": "12841365", "title": "A length polymorphism in the circadian clock gene Per3 is linked to delayed sleep phase syndrome and extreme diurnal preference.", "year": 2003},
        ],
        "actionableInsights": [
            "G/G carriers: If you lean toward evenings, bright light in the morning can help anchor your schedule",
            "Everyone: Maintaining regular sleep schedules supports sleep quality regardless of genotype",
        ],
    },
    {
        "rsid": "rs5751876",
        "gene": "ADORA2A",
        "name": "Caffeine Sensitivity and Sleep",
        "chromosome": "22",
        "position": 24822676,
        "effectSize": "moderate",
        "category": "sleepArchitecture",
        "genotypes": {
            "TT": {
                "label": "T/T (1083T)",
                "nickname": "Mixed Evidence",
                "interpretation": (
                    "Two copies of the T allele (c.1083T>C, also called 1976T>C). Studies point in "
                    "different directions: T/T carriers reported more anxiety after 150 mg of caffeine "
                    "in one study (Alsene et al. 2003), while in 1,023 workers the T allele was linked to "
                    "fewer sleep complaints and less insomnia (Erblang et al. 2019). Limiting caffeine "
                    "late in the day is sensible advice for everyone."
                ),
                "score": 0.5,
            },
            "CT": {
                "label": "C/T",
                "nickname": "Mixed Evidence",
                "interpretation": (
                    "One copy of each allele. Findings for this variant point in different directions "
                    "for anxiety and for sleep, so no clear tendency can be given. Limiting caffeine "
                    "late in the day is good practice for sleep quality."
                ),
                "score": 0.5,
            },
            "CC": {
                "label": "C/C (1083C)",
                "nickname": "Mixed Evidence",
                "interpretation": (
                    "Two copies of the C allele. Studies point in different directions: the caffeine "
                    "anxiety increase was reported for T/T rather than C/C carriers (Alsene et al. 2003), "
                    "while T-allele carriers reported fewer sleep complaints than C/C carriers in a "
                    "workplace study (Erblang et al. 2019). Limiting late-day caffeine remains good advice."
                ),
                "score": 0.5,
            },
        },
        "scientificBasis": (
            "ADORA2A encodes the adenosine A2A receptor, one of the targets of caffeine. Retey et al. "
            "(2007, Clin Pharmacol Ther) found that c.1083T>C (rs5751876) genotype distributions differed "
            "between self-rated caffeine-sensitive and -insensitive adults, and that genotype determined how "
            "caffeine changed the sleep EEG. Byrne et al. (2012, Sleep), in a genome-wide study of "
            "caffeine-induced insomnia, replicated an association at ADORA2A. For individual alleles the "
            "findings conflict: T/T carriers reported more anxiety after caffeine (Alsene et al. 2003), "
            "while the T allele was linked to fewer sleep complaints and less insomnia, but also shorter "
            "sleep in low caffeine consumers (Erblang et al. 2019). Because the directions conflict, this "
            "marker does not shift your sleep score."
        ),
        "references": [
            {"pmid": "17329997", "title": "A genetic variation in the adenosine A2A receptor gene (ADORA2A) contributes to individual sensitivity to caffeine effects on sleep.", "year": 2007},
            {"pmid": "12825092", "title": "Association between A2a receptor gene polymorphisms and caffeine-induced anxiety.", "year": 2003},
            {"pmid": "22754043", "title": "A genome-wide association study of caffeine-related sleep disturbance: confirmation of a role for a common variant in the adenosine receptor.", "year": 2012},
            {"pmid": "31817803", "title": "The Impact of Genetic Variations in ADORA2A in the Association between Caffeine Consumption and Sleep.", "year": 2019},
        ],
        "actionableInsights": [
            "Everyone: Limiting caffeine late in the day supports sleep",
            "If caffeine makes you jittery or disturbs your sleep, try cutting back after early afternoon",
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
        "name": "Circadian Rhythm and Evening Preference",
        "chromosome": "4",
        "position": 56292869,
        "effectSize": "moderate",
        "category": "chronotype",
        "genotypes": {
            "CC": {
                "label": "C/C (3111C)",
                "nickname": "Evening Chronotype",
                "interpretation": (
                    "Two copies of the 3111C allele (reported as G/G on the forward strand "
                    "that most DNA files use). This genotype was associated with evening "
                    "preference and later sleep timing in Katzenberg et al. (1998) and Mishima "
                    "et al. (2005). You may find it easier to stay up late and harder to wake "
                    "early."
                ),
                "score": 0.7,
            },
            "CT": {
                "label": "C/T",
                "nickname": "Mild Evening Tendency",
                "interpretation": (
                    "One copy of the 3111C allele. Carriers in the original study leaned "
                    "slightly toward eveningness."
                ),
                "score": 0.5,
            },
            "TT": {
                "label": "T/T (3111T)",
                "nickname": "Standard Circadian Cycle",
                "interpretation": (
                    "Two copies of the 3111T allele (A/A on the forward strand). This is the "
                    "most common genotype in the original studies, with relatively earlier "
                    "sleep timing than C carriers."
                ),
                "score": 0.3,
            },
        },
        "scientificBasis": (
            "CLOCK 3111T/C (rs1801260) is in the 3' UTR of the CLOCK gene, which encodes a core "
            "clock transcription factor that forms the CLOCK:BMAL1 heterodimer activating "
            "transcription of Per and Cry genes. In two studies the 3111C allele was the evening allele. "
            "Katzenberg et al. (1998, Sleep) found in 410 adults that carriers of 3111C had lower "
            "Horne-Ostberg scores (a shift toward eveningness) and a 10-44 minute delay in preferred "
            "timing of activity and sleep. Mishima et al. (2005, Am J Med Genet B) found in 421 "
            "Japanese subjects that 3111C/C homozygotes had the strongest eveningness, delayed sleep "
            "onset, shorter sleep and more daytime sleepiness compared with T carriers. The evidence "
            "is mixed: Iwase et al. (2002, Psychiatry Res) found the 3111C allele was less frequent, "
            "not more, in patients with delayed sleep phase syndrome."
        ),
        "references": [
            {"pmid": "9779516", "title": "A CLOCK polymorphism associated with human diurnal preference.", "year": 1998},
            {"pmid": "15578592", "title": "The 3111T/C polymorphism of hClock is associated with evening preference and delayed sleep timing in a Japanese population sample.", "year": 2005},
            {"pmid": "11927136", "title": "Mutation screening of the human Clock gene in circadian rhythm sleep disorders.", "year": 2002},
        ],
        "actionableInsights": [
            "C/C carriers: Respect your natural timing where you can, and avoid chronic sleep loss from forcing very early schedules",
            "Morning sunlight exposure helps anchor the circadian cycle regardless of genotype",
            "C/C carriers who need morning schedules: gradually shifting bedtime earlier and getting bright light soon after waking may help",
        ],
    },
    {
        "rsid": "rs2287161",
        "gene": "CRY1",
        "name": "CRY1 Region Variant (context only)",
        "chromosome": "12",
        "position": 107484498,
        "effectSize": "strong",
        "category": "chronotype",
        "genotypes": {
            "CC": {
                "label": "C/C",
                "nickname": "Common Genotype",
                "interpretation": (
                    "No established effect on circadian period or chronotype is known for this "
                    "genotype. It is shown for context only."
                ),
                "score": 0.5,
            },
            "CG": {
                "label": "C/G",
                "nickname": "Heterozygous",
                "interpretation": (
                    "No established effect on circadian period or chronotype is known for this "
                    "genotype. It is shown for context only."
                ),
                "score": 0.5,
            },
            "GG": {
                "label": "G/G",
                "nickname": "Variant Genotype",
                "interpretation": (
                    "No established effect on circadian period or chronotype is known for this "
                    "genotype. It is shown for context only."
                ),
                "score": 0.5,
            },
        },
        "scientificBasis": (
            "CRY1 (Cryptochrome 1) is a central repressor of the molecular circadian clock: it "
            "inhibits the CLOCK:BMAL1 heterodimer, closing the feedback loop. Patke et al. (2017, "
            "Cell) found that a rare dominant CRY1 variant causes familial delayed sleep phase "
            "disorder (DSPD). rs2287161 is a different, common variant near CRY1 and has not been "
            "shown to tag Patke's variant. It has been studied mainly for an interaction with "
            "carbohydrate intake on insulin resistance (Dashti et al., 2014, Chronobiol Int), and "
            "has no established effect on chronotype or circadian period. This marker is shown for "
            "context only and does not shift the chronotype score."
        ),
        "references": [
            {"pmid": "28388406", "title": "Mutation of the Human Circadian Clock Gene CRY1 in Familial Delayed Sleep Phase Disorder.", "year": 2017},
            {"pmid": "24548145", "title": "CRY1 circadian gene variant interacts with carbohydrate intake for insulin resistance in two independent populations: Mediterranean and North American.", "year": 2014},
        ],
        "actionableInsights": [
            "Morning light exposure and regular schedules support a healthy circadian rhythm for everyone",
            "If work allows, consider scheduling demanding tasks around your own natural cognitive peak",
        ],
    },
    {
        "rsid": "rs2735611",
        "gene": "PER1",
        "name": "PER1 Region Variant (context only)",
        "chromosome": "17",
        "position": 8123949,
        "effectSize": "moderate",
        "category": "chronotype",
        "genotypes": {
            "AA": {
                "label": "A/A",
                "nickname": "Common Genotype",
                "interpretation": (
                    "No established direction of effect on chronotype is known for rs2735611. "
                    "It is shown for context only."
                ),
                "score": 0.5,
            },
            "AG": {
                "label": "A/G",
                "nickname": "Heterozygous",
                "interpretation": (
                    "No established direction of effect on chronotype is known for rs2735611. "
                    "It is shown for context only."
                ),
                "score": 0.5,
            },
            "GG": {
                "label": "G/G",
                "nickname": "Variant Genotype",
                "interpretation": (
                    "No established direction of effect on chronotype is known for rs2735611. "
                    "It is shown for context only."
                ),
                "score": 0.5,
            },
        },
        "scientificBasis": (
            "PER1 (Period 1) is one of three 'period' genes of the circadian clock. Chronotype is "
            "polygenic: Jones et al. (2019, Nat Commun) increased the number of known chronotype "
            "loci from 24 to 351 using 697,828 UK Biobank and 23andMe participants. Carpen et al. "
            "(2006, J Hum Genet) reported a synonymous PER1 variant (T2434C) linked to extreme "
            "morning preference, but whether rs2735611 is that variant is not established. This "
            "marker is shown for context and does not shift the chronotype score."
        ),
        "references": [
            {"pmid": "30696823", "title": "Genome-wide association analyses of chronotype in 697,828 individuals provides insights into circadian rhythms.", "year": 2019},
            {"pmid": "17051316", "title": "A silent polymorphism in the PER1 gene associates with extreme diurnal preference in humans.", "year": 2006},
        ],
        "actionableInsights": [
            "Consider scheduling demanding work around the time of day when you personally feel sharpest",
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
        "name": "Sleep Paralysis and Hypnagogic Experiences (exploratory)",
        "chromosome": "13",
        "position": 46897065,
        "effectSize": "moderate",
        "category": "nocturnal",
        "genotypes": {
            "TT": {
                "label": "T/T",
                "nickname": "Higher expected expression",
                "interpretation": (
                    "Two copies of the 102T allele, which was expressed more than the 102C allele "
                    "in human cortex, suggesting somewhat higher 5-HT2A receptor expression. "
                    "Whether genotype relates to sleep paralysis or hypnagogic experiences has "
                    "not been tested; this link is exploratory."
                ),
                "score": 0.8,
            },
            "CT": {
                "label": "C/T",
                "nickname": "Intermediate expected expression",
                "interpretation": (
                    "One copy of each allele; expression is expected to be intermediate. Whether "
                    "this relates to sleep paralysis or hypnagogic experiences has not been "
                    "tested; this link is exploratory."
                ),
                "score": 0.5,
            },
            "CC": {
                "label": "C/C",
                "nickname": "Lower expected expression",
                "interpretation": (
                    "Two copies of the 102C allele, which was expressed about 20% less than 102T "
                    "in human cortex, suggesting somewhat lower 5-HT2A receptor expression. "
                    "Whether genotype relates to sleep paralysis or hypnagogic experiences has "
                    "not been tested; this link is exploratory."
                ),
                "score": 0.3,
            },
        },
        "scientificBasis": (
            "HTR2A rs6313 (102T/C) is a synonymous variant in exon 1 of the serotonin 2A receptor "
            "gene (not a promoter polymorphism) and is in strong linkage disequilibrium with rs6311. "
            "Polesskaya & Sokolov (2002, J Neurosci Res) found the 102T allele expressed more than "
            "102C in human temporal cortex. Sleep paralysis (SP) occurs when REM-state muscle "
            "paralysis persists into wakefulness. A systematic review (Denis et al., 2018, Sleep Med "
            "Rev) found it is associated with substance use, stress and trauma, genetic influences, "
            "physical illness, personality, anomalous beliefs, sleep problems or disruption, and "
            "anxiety symptoms. No study has linked HTR2A genotype to sleep paralysis, so this "
            "marker is shown as an exploratory hypothesis. DNA files report this SNP on the "
            "forward strand, where 102T appears as A and 102C as G. The score reflects expected "
            "receptor expression, not the likelihood of any nocturnal experience."
        ),
        "references": [
            {"pmid": "11891796", "title": "Differential expression of the \"C\" and \"T\" alleles of the 5-HT2A receptor gene in the temporal cortex of normal individuals and schizophrenics.", "year": 2002},
            {"pmid": "28735779", "title": "A systematic review of variables associated with sleep paralysis.", "year": 2018},
        ],
        "actionableInsights": [
            "If sleep paralysis is disturbing, keeping a regular sleep schedule and reducing stress and sleep disruption are sensible first steps",
            "Sleeping on your side rather than your back is a commonly suggested tip for some people",
            "Meditative practices such as yoga nidra or NSDR can be helpful for relaxation around sleep",
            "If episodes are frequent or distressing, consult a sleep medicine specialist",
        ],
    },
    {
        "rsid": "rs6265",
        "gene": "BDNF",
        "name": "Sleep Intensity and Neuroplasticity",
        "chromosome": "11",
        "position": 27658369,
        "effectSize": "moderate",
        "category": "nocturnal",
        "genotypes": {
            "CC": {
                "label": "Val/Val",
                "nickname": "Standard Sleep Profile",
                "interpretation": (
                    "Reference genotype with normal activity-dependent BDNF secretion. In a small "
                    "study (Bachmann et al., 2012), Val/Val participants had higher deep (stage 4) "
                    "sleep and NREM slow-wave activity than Val/Met participants."
                ),
                "score": 0.3,
            },
            "CT": {
                "label": "Val/Met",
                "nickname": "Intermediate Profile",
                "interpretation": (
                    "One copy of the Met allele, which impairs activity-dependent secretion of "
                    "mature BDNF. In Bachmann et al. (2012), Val/Met participants showed lower "
                    "deep-sleep intensity than Val/Val. No study has measured dream recall or "
                    "content by BDNF genotype, so any dream link is exploratory."
                ),
                "score": 0.5,
            },
            "TT": {
                "label": "Met/Met",
                "nickname": "Altered Sleep Intensity",
                "interpretation": (
                    "Two copies of the Met allele (Val66Met), which impairs activity-dependent "
                    "secretion of mature BDNF. In Bachmann et al. (2012), Met carriers showed "
                    "lower deep-sleep intensity than Val/Val. No study has measured dream recall "
                    "or content by BDNF genotype, so any dream link is exploratory. Regular "
                    "exercise supports overall brain health and sleep."
                ),
                "score": 0.65,
            },
        },
        "scientificBasis": (
            "BDNF Val66Met (rs6265) is also in GeneHealth's mind_spirit_analyzer in the context of "
            "daytime neuroplasticity. The Met allele impairs activity-dependent secretion of mature "
            "BDNF. Tononi & Cirelli (2014, Neuron) review how sleep supports synaptic homeostasis and "
            "memory consolidation. Bachmann et al. (2012, Sleep) compared 11 Val/Met with 11 Val/Val "
            "participants and found deep stage-4 sleep and NREM slow-wave activity were higher in "
            "Val/Val than Val/Met in both baseline and recovery nights, and Val/Val had better 2-back "
            "accuracy. No study has measured dream recall or content by BDNF genotype, so any "
            "dream-related interpretation is exploratory."
        ),
        "references": [
            {"pmid": "24411729", "title": "Sleep and the price of plasticity: from synaptic and cellular homeostasis to memory consolidation and integration.", "year": 2014},
            {"pmid": "22379239", "title": "The BDNF Val66Met polymorphism modulates sleep intensity: EEG frequency- and state-specificity.", "year": 2012},
        ],
        "actionableInsights": [
            "Regular aerobic exercise supports sleep and brain health for everyone",
            "Pre-sleep journaling can help quiet a busy mind before bed",
            "Sleeping well after important learning is a good habit for every genotype",
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
            0.65: "Your markers lean toward higher serotonin 2A receptor expression and a higher novelty-seeking tendency. Whether this relates to dream life is an untested hypothesis.",
            0.45: "Your markers fall in an intermediate range. Any link to dream life is an untested hypothesis.",
            0.0: "Your markers lean toward lower expression or lower novelty-seeking scores. Any link to dream life is an untested hypothesis.",
        },
        "sleepArchitecture": {
            0.65: "Your markers here lean toward evening timing preference. Regular schedules and limiting late caffeine are sensible for everyone.",
            0.45: "Intermediate profile for sleep timing; regular schedules help everyone.",
            0.0: "Your markers here lean toward typical sleep timing.",
        },
        "chronotype": {
            0.65: "Your biological clock has a strong evening (night owl) tendency with a solid genetic basis.",
            0.45: "Intermediate chronotype — you adapt to different schedules with relative ease.",
            0.0: "Strong predisposition to a morning chronotype — you likely function best early.",
        },
        "nocturnal": {
            0.65: "Your markers lean toward higher HTR2A expression and altered sleep intensity. Links to sleep paralysis or hypnagogia are untested hypotheses.",
            0.45: "Intermediate markers; links to nocturnal phenomena are untested hypotheses. Stress and sleep disruption are known contributors to sleep paralysis.",
            0.0: "Markers lean toward lower expression; links to nocturnal phenomena are untested hypotheses.",
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
        # All dream-vividness markers are exploratory (no study has measured dreams by
        # genotype): describe the marker pattern, not a kind of dreamer.
        "dreamVividness": {
            0.65: "Higher-Score Marker Pattern (exploratory)",
            0.45: "Intermediate Marker Pattern (exploratory)",
            0.0: "Lower-Score Marker Pattern (exploratory)",
        },
        "sleepArchitecture": {
            0.65: "Evening-Leaning Timing",
            0.45: "Intermediate Profile",
            0.0: "Typical Timing",
        },
        "chronotype": {
            0.65: "Genetic Night Owl",
            0.45: "Intermediate Chronotype",
            0.0: "Genetic Early Bird",
        },
        "nocturnal": {
            0.65: "Higher HTR2A Expression",
            0.45: "Intermediate Profile",
            0.0: "Lower HTR2A Expression",
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
                "report. Whether the same variants relate to experiences during wakefulness "
                "and during sleep is a hypothesis that has not been tested."
            ),
        },
    }
