"""
Mind & Spirit Analyzer
Analyzes genetic variants associated with personality traits, mental health
predispositions, and spiritual sensitivity (including mediumship genetics).

THREE CATEGORIES:
  1. Personality Traits - COMT, CADM2, MSRA, OPRM1, OXTR
  2. Mental Health & Stress Response - FKBP5, BDNF, MTHFR, CACNA1C, RGS2, SNAP25
  3. Spiritual Sensitivity - HTR2A, TPH2 (exploratory), OXTR (absorption context)

SCIENTIFIC BASIS:
  - Personality: GWAS-validated SNPs (Nature Human Behaviour 2024, Translational Psychiatry 2023)
  - Mental Health: PGC consortium studies (Cell 2024, Nature Genetics 2025)
  - Spiritual: Gattaz et al. 2025, Braz J Psychiatry (PMID 39874024), HTR2A mystical experience research

IMPORTANT CAVEATS:
  - All associations are probabilistic, not deterministic
  - Effect sizes are generally small (each SNP explains <1% of trait variance)
  - Environment, culture, and personal choice strongly modify any genetic tendency
  - Mental health reports are educational only, not diagnostic
  - Spiritual sensitivity findings are preliminary and exploratory

References:
  - PGC MDD Working Group 2024, Cell (PMC11092713) - 697 depression loci
  - Gattaz et al. 2025, Braz J Psychiatry (PMID 39874024) - mediumship candidate genes
  - Pasternak & Zangari 2025, Questao de Ciencia - critical review
  - Sanchez-Roige et al. 2023, Translational Psychiatry - CADM2 impulsivity
"""

from typing import Dict, List, Tuple, Any, Optional


# ---------------------------------------------------------------------------
# SNP DATABASE: Personality Traits
# ---------------------------------------------------------------------------

PERSONALITY_SNPS: List[Dict[str, Any]] = [
    {
        "rsid": "rs4680",
        "gene": "COMT",
        "name": "Dopamine Regulation and Cognitive Style",
        "chromosome": "22",
        "position": 19963748,
        "effectSize": "moderate",
        "genotypes": {
            "GG": {
                "label": "Val/Val",
                "nickname": "Warrior",
                "interpretation": (
                    "Higher COMT enzyme activity leads to faster dopamine clearance in the "
                    "prefrontal cortex. Associated with better stress tolerance, higher pain threshold, "
                    "and advantage under pressure. May show slightly lower baseline cognitive performance, "
                    "but better performance under stress."
                ),
                "score": 0.7,
            },
            "AG": {
                "label": "Val/Met",
                "nickname": "Balanced",
                "interpretation": (
                    "Intermediate COMT activity provides a balanced dopamine profile. "
                    "You likely have a versatile cognitive style — adequate stress tolerance "
                    "with good baseline cognitive flexibility. Most common genotype."
                ),
                "score": 0.5,
            },
            "AA": {
                "label": "Met/Met",
                "nickname": "Worrier",
                "interpretation": (
                    "Lower COMT enzyme activity means slower dopamine clearance, resulting "
                    "in higher prefrontal dopamine levels. Associated with superior baseline cognitive "
                    "performance (better working memory, attention), but potentially greater "
                    "vulnerability to stress and anxiety under pressure."
                ),
                "score": 0.3,
            },
        },
        "scientificBasis": (
            "COMT Val158Met (rs4680) is one of the most studied functional SNPs in behavioral genetics. "
            "The Val allele (G) is the predominant determinant of higher COMT enzyme activity in "
            "human prefrontal cortex (Chen et al. 2004), compared with the Met allele (A). A "
            "meta-analysis of cognitive studies found little if any consistent effect on cognition "
            "(Barnett et al. 2008), so this marker should be read as a tendency, not a prediction. "
            "The 'Warrior vs Worrier' model (Stein et al. 2006) captures the proposed trade-off: Val/Val excels under stress but has lower "
            "baseline prefrontal dopamine, while Met/Met excels in baseline cognition but is more "
            "sensitive to stress."
        ),
        "references": [
            {"pmid": "17008817", "title": "Warriors versus worriers: the role of COMT gene variants.", "year": 2006},
            {"pmid": "18339359", "title": "Meta-analysis of the cognitive effects of the catechol-O-methyltransferase gene Val158/108Met polymorphism.", "year": 2008},
            {"pmid": "15457404", "title": "Functional analysis of genetic variation in catechol-O-methyltransferase (COMT): effects on mRNA, protein, and enzyme activity in postmortem human brain.", "year": 2004},
        ],
        "actionableInsights": [
            "Val/Val: You may thrive in high-pressure environments — consider leveraging this in career choices",
            "Met/Met: Mindfulness and stress-reduction practices may be especially beneficial for you",
            "Regular physical exercise helps optimize dopamine levels regardless of genotype",
        ],
    },
    {
        "rsid": "rs17518584",
        "gene": "CADM2",
        "name": "Risk Tolerance and Impulsivity",
        "chromosome": "3",
        "position": 85890189,
        "effectSize": "moderate",
        "genotypes": {
            "CC": {
                "label": "Typical",
                "nickname": "Cautious",
                "interpretation": (
                    "Reference genotype associated with typical levels of risk tolerance. "
                    "You likely have a balanced approach to risk decisions."
                ),
                "score": 0.3,
            },
            "CT": {
                "label": "Intermediate",
                "nickname": "Moderate",
                "interpretation": (
                    "One copy of the risk-associated allele. Slightly elevated tendency "
                    "toward risk-taking behavior and novelty seeking compared to CC carriers."
                ),
                "score": 0.5,
            },
            "TT": {
                "label": "Elevated",
                "nickname": "Risk-taker",
                "interpretation": (
                    "Two copies of the risk-associated allele. Associated with higher risk "
                    "tolerance, novelty seeking, and impulsivity in GWAS studies. "
                    "Also linked to substance use tendencies in population studies."
                ),
                "score": 0.8,
            },
        },
        "scientificBasis": (
            "CADM2 (Cell Adhesion Molecule 2) rs17518584 is one of the best-validated GWAS findings "
            "for risk-taking behavior. Replicated in UK Biobank (n>400,000) and 23andMe research "
            "cohorts. CADM2 is expressed in the brain and plays a role in synaptic development. "
            "Effect sizes are small but consistent (Sanchez-Roige et al. 2023, Translational "
            "Psychiatry). A GWAS of over 1 million individuals identified 99 loci for general "
            "risk tolerance (Karlsson Linnér et al. 2019, Nat Genet)."
        ),
        "references": [
            {"pmid": "37173343", "title": "CADM2 is implicated in impulsive personality and numerous other traits by genome- and phenome-wide association studies in humans and mice.", "year": 2023},
            {"pmid": "30643258", "title": "Genome-wide association analyses of risk tolerance and risky behaviors in over 1 million individuals identify hundreds of loci and shared genetic influences.", "year": 2019},
        ],
        "actionableInsights": [
            "TT carriers: Being aware of elevated impulsivity tendencies can help with financial and health decisions",
            "Consider using decision-making frameworks for important life choices",
            "Channel risk-taking tendencies constructively through entrepreneurship or sports",
        ],
    },
    {
        "rsid": "rs4925638",
        "gene": "MSRA",
        "name": "Irritability and Behavioral Inhibition",
        "chromosome": "8",
        "position": 10040702,
        "effectSize": "moderate",
        "genotypes": {
            "AA": {
                "label": "Typical",
                "nickname": "Steady temperament",
                "interpretation": (
                    "Reference genotype associated with typical levels of emotional reactivity "
                    "and behavioral inhibition."
                ),
                "score": 0.3,
            },
            "AG": {
                "label": "Intermediate",
                "nickname": "Moderate reactivity",
                "interpretation": (
                    "One copy of the variant allele. Slightly elevated emotional reactivity "
                    "compared to AA carriers."
                ),
                "score": 0.5,
            },
            "GG": {
                "label": "Elevated",
                "nickname": "Elevated reactivity",
                "interpretation": (
                    "Two copies of the variant allele. Associated with higher irritability "
                    "and reduced behavioral inhibition in GWAS studies. MSRA plays a role "
                    "in protection against oxidative stress in the brain."
                ),
                "score": 0.7,
            },
        },
        "scientificBasis": (
            "MSRA (Methionine Sulfoxide Reductase A) rs4925638 has been associated with irritability, "
            "risk tolerance, and impulsivity in UK Biobank GWAS analyses. MSRA protects proteins "
            "against oxidative damage — variants may affect neuronal stress resilience and emotional "
            "regulation pathways. Replicated by Boutwell et al. 2017 (Heliyon)."
        ),
        "references": [
            {"pmid": "28795158", "title": "Replication and characterization of CADM2 and MSRA genes on human behavior.", "year": 2017},
        ],
        "actionableInsights": [
            "GG carriers: An antioxidant-rich diet (berries, leafy greens) may support MSRA function",
            "Regular physical exercise is one of the strongest moderators of irritability",
            "Cognitive-behavioral techniques can help manage elevated emotional reactivity",
        ],
    },
    {
        "rsid": "rs1799971",
        "gene": "OPRM1",
        "name": "Emotional Sensitivity and Reward Processing",
        "chromosome": "6",
        "position": 154360797,
        "effectSize": "moderate",
        "genotypes": {
            "AA": {
                "label": "Asp/Asp",
                "nickname": "Typical sensitivity",
                "interpretation": (
                    "Reference genotype for the mu-opioid receptor (A118). Most common "
                    "genotype worldwide."
                ),
                "score": 0.4,
            },
            "AG": {
                "label": "Asp/Asn",
                "nickname": "Increased sensitivity",
                "interpretation": (
                    "One copy of the G allele (A118G). In one study, G allele carriers reported "
                    "higher sensitivity to social rejection and showed stronger brain responses "
                    "to it (Way et al. 2009)."
                ),
                "score": 0.6,
            },
            "GG": {
                "label": "Asn/Asn",
                "nickname": "High sensitivity",
                "interpretation": (
                    "Two copies of the G allele. G allele carriers reported higher sensitivity to "
                    "social rejection in one study (Way et al. 2009). Response to opioid pain "
                    "medication may also differ; discuss pain management with your clinician."
                ),
                "score": 0.8,
            },
        },
        "scientificBasis": (
            "OPRM1 A118G (rs1799971) affects the mu-opioid receptor, which mediates pain, reward, and "
            "social bonding. The G allele yielded lower OPRM1 mRNA and protein in human brain tissue and "
            "cell models (Zhang et al. 2005), and binds beta-endorphin about three times more tightly "
            "in vitro (Bond et al. 1998). G allele carriers reported higher dispositional sensitivity "
            "to rejection and showed greater dACC and anterior insula responses to social "
            "rejection (Way et al. 2009, PNAS)."
        ),
        "references": [
            {"pmid": "19706472", "title": "Variation in the mu-opioid receptor gene (OPRM1) is associated with dispositional and neural sensitivity to social rejection.", "year": 2009},
            {"pmid": "16046395", "title": "Allelic expression imbalance of human mu opioid receptor (OPRM1) caused by variant A118G.", "year": 2005},
            {"pmid": "9689128", "title": "Single-nucleotide polymorphism in the human mu opioid receptor gene alters beta-endorphin binding and activity: possible implications for opiate addiction.", "year": 1998},
        ],
        "actionableInsights": [
            "G carriers: Your elevated emotional sensitivity can be a strength in empathic relationships",
            "Be aware that pain perception and medication response may differ from population averages",
            "Social connection and support networks may be especially important for G carriers",
        ],
    },
    {
        "rsid": "rs53576",
        "gene": "OXTR",
        "name": "Social Sensitivity and Empathy",
        "chromosome": "3",
        "position": 8804371,
        "effectSize": "weak",
        "genotypes": {
            "GG": {
                "label": "GG",
                "nickname": "Higher empathic tendency",
                "interpretation": (
                    "Associated with higher empathy scores, greater social sensitivity, "
                    "and stronger emotional reactivity to social cues in early studies. "
                    "Note: large-scale replications have shown mixed results — effect sizes "
                    "are very small."
                ),
                "score": 0.7,
            },
            "AG": {
                "label": "AG",
                "nickname": "Intermediate",
                "interpretation": (
                    "Intermediate oxytocin receptor sensitivity. Typical social bonding "
                    "tendencies. Most common genotype."
                ),
                "score": 0.5,
            },
            "AA": {
                "label": "AA",
                "nickname": "Lower social sensitivity",
                "interpretation": (
                    "Associated with slightly lower social sensitivity and empathy scores "
                    "in some studies. However, a broader meta-analysis across social "
                    "outcomes found effects indistinguishable from zero."
                ),
                "score": 0.3,
            },
        },
        "scientificBasis": (
            "OXTR rs53576 is the most studied oxytocin receptor variant for social behavior. "
            "Early studies linked the GG genotype to higher empathy and prosocial behavior. "
            "A population study plus meta-analysis found a significant association with empathy "
            "in Europeans and Asians (Gong et al. 2017, PMID 28343138). However, a broader "
            "meta-analysis across social outcomes found combined effects indistinguishable from "
            "zero (Bakermans-Kranenburg & van IJzendoorn 2014). This SNP "
            "is included for completeness, but results should be interpreted with caution."
        ),
        "references": [
            {"pmid": "28343138", "title": "Revisiting the impact of OXTR rs53576 on empathy: A population-based study and a meta-analysis.", "year": 2017},
            {"pmid": "23921259", "title": "A sociability gene? Meta-analysis of oxytocin receptor genotype effects in humans.", "year": 2014},
        ],
        "actionableInsights": [
            "Social behavior is predominantly shaped by environment and personal experiences",
            "Empathy is a skill that can be developed regardless of genotype",
            "Strong social connections benefit everyone — invest in relationships",
        ],
    },
    {
        "rsid": "rs2254298",
        "gene": "OXTR",
        "name": "Social Bonding and Attachment",
        "chromosome": "3",
        "position": 8798609,
        "effectSize": "weak",
        "genotypes": {
            "GG": {
                "label": "GG",
                "nickname": "Reference",
                "interpretation": (
                    "Reference genotype for this OXTR variant. Typical attachment "
                    "and social bonding patterns."
                ),
                "score": 0.4,
            },
            "GA": {
                "label": "GA",
                "nickname": "Intermediate",
                "interpretation": (
                    "One copy of the A allele. Some studies link it to altered social "
                    "bonding patterns, but findings are inconsistent across populations."
                ),
                "score": 0.5,
            },
            "AA": {
                "label": "AA",
                "nickname": "Variant",
                "interpretation": (
                    "Two copies of the A allele. Linked to attachment or social bonding in "
                    "some individual studies, but a meta-analysis found no significant "
                    "association with social outcomes."
                ),
                "score": 0.6,
            },
        },
        "scientificBasis": (
            "OXTR rs2254298 has been studied for associations with social bonding, attachment security, "
            "and prosocial behavior. Minor-allele frequency differs by ancestry, but a "
            "meta-analysis of 34 effect sizes (N = 13,547) found no significant association "
            "with social outcomes (Bakermans-Kranenburg & van IJzendoorn 2014). Limited "
            "predictive value for individual social behavior."
        ),
        "references": [
            {"pmid": "23921259", "title": "A sociability gene? Meta-analysis of oxytocin receptor genotype effects in humans.", "year": 2014},
        ],
        "actionableInsights": [
            "Attachment styles are primarily shaped by childhood caregiving experiences",
            "Secure attachment can be developed at any age through therapy and relationships",
        ],
    },
]


# ---------------------------------------------------------------------------
# SNP DATABASE: Mental Health and Stress Response
# ---------------------------------------------------------------------------

MENTAL_HEALTH_SNPS: List[Dict[str, Any]] = [
    {
        "rsid": "rs1360780",
        "gene": "FKBP5",
        "name": "Stress Response and HPA Axis Regulation",
        "chromosome": "6",
        "position": 35639794,
        "effectSize": "strong",
        "genotypes": {
            "CC": {
                "label": "CC",
                "nickname": "Standard stress response",
                "interpretation": (
                    "Reference genotype. Normal FKBP5 expression and glucocorticoid receptor "
                    "feedback. Standard cortisol stress response and recovery."
                ),
                "score": 0.3,
            },
            "CT": {
                "label": "CT",
                "nickname": "Moderately increased stress sensitivity",
                "interpretation": (
                    "One copy of the risk T allele. Moderately increased FKBP5 expression "
                    "after stress, which may slow cortisol feedback recovery. "
                    "Mildly elevated stress sensitivity, particularly in the context of "
                    "adverse experiences."
                ),
                "score": 0.5,
            },
            "TT": {
                "label": "TT",
                "nickname": "Increased stress sensitivity",
                "interpretation": (
                    "Two copies of the T allele (high-induction genotype). Significant increase "
                    "in FKBP5 mRNA expression after stress exposure, impairing "
                    "glucocorticoid receptor feedback and prolonging cortisol response. "
                    "Strongest association with PTSD and stress-related depression "
                    "in the context of childhood adversity."
                ),
                "score": 0.8,
            },
        },
        "scientificBasis": (
            "FKBP5 rs1360780 is one of the most studied single-SNP findings in stress-related "
            "psychiatry. In people exposed to childhood trauma, the risk allele enables "
            "allele-specific DNA demethylation in FKBP5, increasing FKBP5 responsiveness "
            "(Klengel et al. 2013); FKBP5 protein binds and inhibits the glucocorticoid receptor, "
            "disrupting cortisol feedback. FKBP5 polymorphisms interacted with childhood abuse to "
            "predict adult PTSD symptoms (Binder et al. 2008, JAMA), an interaction that has been "
            "reported in several cohorts. Effect sizes are modest."
        ),
        "references": [
            {"pmid": "18349090", "title": "Association of FKBP5 polymorphisms and childhood abuse with risk of posttraumatic stress disorder symptoms in adults.", "year": 2008},
            {"pmid": "23201972", "title": "Allele-specific FKBP5 DNA demethylation mediates gene-childhood trauma interactions.", "year": 2013},
        ],
        "actionableInsights": [
            "TT carriers: Stress management practices (meditation, exercise) may be especially important",
            "Early intervention for stress-related symptoms is valuable for T allele carriers",
            "Therapeutic approaches targeting stress reactivity (CBT, EMDR) are evidence-based options",
            "This variant does NOT determine your mental health — environment and coping strategies matter enormously",
        ],
    },
    {
        "rsid": "rs6265",
        "gene": "BDNF",
        "name": "Neuroplasticity and Mood Regulation",
        "chromosome": "11",
        "position": 27658369,
        "effectSize": "moderate",
        "genotypes": {
            "CC": {
                "label": "Val/Val",
                "nickname": "Standard neuroplasticity",
                "interpretation": (
                    "Reference genotype (Val66). Normal BDNF secretion and activity-dependent "
                    "release. Standard neuroplasticity and hippocampal function."
                ),
                "score": 0.3,
            },
            "CT": {
                "label": "Val/Met",
                "nickname": "Moderately altered",
                "interpretation": (
                    "One copy of the Met allele. Reduced activity-dependent BDNF secretion. "
                    "Subtle effects on memory and hippocampal volume in some studies. The association "
                    "with depression per se is weak at the population level."
                ),
                "score": 0.5,
            },
            "TT": {
                "label": "Met/Met",
                "nickname": "Reduced BDNF secretion",
                "interpretation": (
                    "Two copies of the Met allele. Significant reduction in activity-dependent "
                    "BDNF release. Associated with lower hippocampal volume, altered "
                    "memory consolidation, and in some studies, modestly elevated depression risk — "
                    "in a meta-analysis, the effect was seen in men."
                ),
                "score": 0.7,
            },
        },
        "scientificBasis": (
            "BDNF Val66Met (rs6265) affects activity-dependent secretion of brain-derived neurotrophic "
            "factor, crucial for synaptic plasticity. A meta-analysis of 14 studies (2,812 cases, "
            "10,843 controls) found no overall association with major depression, but a "
            "significant effect in men (Met allele OR 1.27) (Verhagen et al. 2010). "
            "BDNF is also studied in exercise and cognition research."
        ),
        "references": [
            {"pmid": "12553913", "title": "The BDNF val66met polymorphism affects activity-dependent secretion of BDNF and human memory and hippocampal function.", "year": 2003},
            {"pmid": "18852698", "title": "Meta-analysis of the BDNF Val66Met polymorphism in major depressive disorder: effects of gender and ethnicity.", "year": 2010},
        ],
        "actionableInsights": [
            "Exercise is the best-known way to increase BDNF levels — especially important for Met carriers",
            "Aerobic exercise (running, swimming, cycling) has the most evidence for BDNF elevation",
            "Met carriers may particularly benefit from cognitive training and novel learning experiences",
        ],
    },
    {
        "rsid": "rs1801133",
        "gene": "MTHFR",
        "name": "Folate Metabolism and Mood",
        "chromosome": "1",
        "position": 11856378,
        "effectSize": "moderate",
        "genotypes": {
            "GG": {
                "label": "CC (normal)",
                "nickname": "Full enzyme activity",
                "interpretation": (
                    "Normal MTHFR enzyme activity (~100%). Efficient folate metabolism "
                    "and methylation cycle. No folate-related mood concerns."
                ),
                "score": 0.2,
            },
            "AG": {
                "label": "CT (heterozygous)",
                "nickname": "Reduced activity (~65%)",
                "interpretation": (
                    "One copy of the T allele (C677T). MTHFR enzyme activity reduced to "
                    "approximately 65% of normal. Generally clinically insignificant, but "
                    "adequate folate intake is advisable."
                ),
                "score": 0.5,
            },
            "AA": {
                "label": "TT (homozygous)",
                "nickname": "Significantly reduced (~30%)",
                "interpretation": (
                    "Two copies of the T allele. MTHFR enzyme activity reduced to approximately "
                    "30% of normal. In a HuGE meta-analysis, TT homozygotes had higher odds of "
                    "unipolar depression than CC (OR 1.36, 95% CI 1.11-1.67) (Gilbody et al. 2007). "
                    "TT prevalence in Brazil is approximately 10%. Adjunctive L-methylfolate has "
                    "been tested in SSRI-resistant depression (Papakostas et al. 2012), although "
                    "that trial did not select patients by MTHFR genotype."
                ),
                "score": 0.8,
            },
        },
        "scientificBasis": (
            "MTHFR C677T (rs1801133) reduces methylenetetrahydrofolate reductase activity, "
            "affecting the one-carbon metabolism cycle that produces SAMe (the primary methyl "
            "donor for neurotransmitter synthesis). TT homozygotes have ~30% enzyme activity "
            "and may develop mild hyperhomocysteinemia. In a HuGE meta-analysis, TT homozygotes "
            "had higher odds of unipolar depression than CC (OR 1.36, 95% CI 1.11-1.67) "
            "(Gilbody et al. 2007). Adjunctive L-methylfolate has been tested in SSRI-resistant "
            "depression (Papakostas et al. 2012), although that trial did not select patients "
            "by MTHFR genotype."
        ),
        "references": [
            {"pmid": "7647779", "title": "A candidate genetic risk factor for vascular disease: a common mutation in methylenetetrahydrofolate reductase.", "year": 1995},
            {"pmid": "17074966", "title": "Methylenetetrahydrofolate reductase (MTHFR) genetic polymorphisms and psychiatric disorders: a HuGE review.", "year": 2007},
            {"pmid": "23212058", "title": "L-methylfolate as adjunctive therapy for SSRI-resistant major depression: results of two randomized, double-blind, parallel-sequential trials.", "year": 2012},
        ],
        "actionableInsights": [
            "TT carriers: Discuss folate status and any supplementation with a healthcare professional",
            "Ensure adequate dietary folate intake: leafy greens, legumes, fortified foods",
            "TT carriers should monitor homocysteine levels with their doctor",
            "B12 and B6 are important cofactors — ensure adequate intake",
        ],
    },
    {
        "rsid": "rs1006737",
        "gene": "CACNA1C",
        "name": "Calcium Channel and Mood Stability",
        "chromosome": "12",
        "position": 2345295,
        "effectSize": "strong",
        "genotypes": {
            "GG": {
                "label": "GG",
                "nickname": "Reference",
                "interpretation": (
                    "Reference genotype. Standard L-type calcium channel function and "
                    "neuronal excitability."
                ),
                "score": 0.2,
            },
            "AG": {
                "label": "AG",
                "nickname": "One risk allele",
                "interpretation": (
                    "One copy of the risk A allele. Slightly altered calcium channel "
                    "expression. Modestly elevated risk for bipolar disorder."
                ),
                "score": 0.5,
            },
            "AA": {
                "label": "AA",
                "nickname": "Two risk alleles",
                "interpretation": (
                    "Two copies of the A allele. One of the earliest and most studied genetic "
                    "risk variants for bipolar disorder. CACNA1C has also been implicated in other "
                    "psychiatric conditions. Most carriers never develop any condition."
                ),
                "score": 0.7,
            },
        },
        "scientificBasis": (
            "CACNA1C rs1006737 is one of the earliest and most studied bipolar disorder GWAS "
            "signals (combined P = 7.0 × 10^-8; Ferreira et al. 2008). It encodes the alpha-1C subunit "
            "of the voltage-dependent L-type calcium channel, critical for neuronal excitability "
            "and synaptic plasticity. CACNA1C has also been implicated in other psychiatric "
            "conditions. The 2021 Psychiatric Genomics Consortium study of more than 40,000 "
            "bipolar disorder cases identified 64 loci (Mullins et al. 2021)."
        ),
        "references": [
            {"pmid": "18711365", "title": "Collaborative genome-wide association analysis supports a role for ANK3 and CACNA1C in bipolar disorder.", "year": 2008},
            {"pmid": "34002096", "title": "Genome-wide association study of more than 40,000 bipolar disorder cases provides new insights into the underlying biology.", "year": 2021},
        ],
        "actionableInsights": [
            "AA carriers: Mood tracking can help identify patterns early",
            "Regular sleep schedules are particularly important for mood stability",
            "This variant increases risk slightly — most carriers never develop mood disorders",
            "If you have a family history of bipolar disorder, discuss with your doctor",
        ],
    },
    {
        "rsid": "rs4606",
        "gene": "RGS2",
        "name": "Anxiety Sensitivity and Stress Reactivity",
        "chromosome": "1",
        "position": 192839582,
        "effectSize": "moderate",
        "genotypes": {
            "CC": {
                "label": "CC",
                "nickname": "Elevated anxiety sensitivity",
                "interpretation": (
                    "Two copies of the C allele. In studies of RGS2, the C allele has been "
                    "associated with higher odds of generalized anxiety disorder and with "
                    "behavioral inhibition and greater amygdala/insula activation to emotional "
                    "faces. RGS2 dampens G-protein signaling, so this is a tendency, "
                    "not a prediction."
                ),
                "score": 0.7,
            },
            "CG": {
                "label": "CG",
                "nickname": "Intermediate",
                "interpretation": (
                    "One copy of the C allele. An intermediate profile; any effect on "
                    "anxiety sensitivity is expected to be subtle."
                ),
                "score": 0.5,
            },
            "GG": {
                "label": "GG",
                "nickname": "Reference",
                "interpretation": (
                    "Reference genotype. Typical anxiety response levels; no anxiety-associated "
                    "C allele copies."
                ),
                "score": 0.3,
            },
        },
        "scientificBasis": (
            "RGS2 (Regulator of G-protein Signaling 2) rs4606 modulates signaling downstream "
            "of corticotropin-releasing hormone (CRH) receptors. RGS2 markers including rs4606 "
            "were associated with childhood behavioral inhibition, introversion and increased "
            "amygdala/insula activation to emotional faces; rs4606 genotype explained 10-15% of the "
            "variance in that brain activation (Smoller et al. 2008, Arch Gen Psychiatry). Each C "
            "allele was associated with higher odds of generalized anxiety disorder (Koenen et al. "
            "2009, Depress Anxiety), consistent with an earlier report of a higher C-allele "
            "prevalence in panic disorder patients. Effect sizes are small."
        ),
        "references": [
            {"pmid": "18316676", "title": "Influence of RGS2 on anxiety-related temperament, personality, and brain function.", "year": 2008},
            {"pmid": "18833580", "title": "RGS2 and generalized anxiety disorder in an epidemiologic sample of hurricane-exposed adults.", "year": 2009},
        ],
        "actionableInsights": [
            "CC carriers: Structured relaxation techniques (progressive muscle relaxation, deep breathing) may be especially helpful",
            "Regular physical exercise is one of the strongest anxiolytic interventions",
            "Consider whether anxiety management strategies would benefit your daily routine",
        ],
    },
    {
        "rsid": "rs3746544",
        "gene": "SNAP25",
        "name": "Synaptic Signaling and Attention",
        "chromosome": "20",
        "position": 10284862,
        "effectSize": "weak",
        "genotypes": {
            "TT": {
                "label": "TT",
                "nickname": "Reference",
                "interpretation": (
                    "Reference genotype. Standard SNAP25 protein function and synaptic "
                    "vesicle release."
                ),
                "score": 0.3,
            },
            "TG": {
                "label": "TG",
                "nickname": "Intermediate",
                "interpretation": (
                    "One copy of the G allele. Subtle effects on synaptic transmission efficiency. "
                    "Some ADHD association studies show modest risk elevation."
                ),
                "score": 0.5,
            },
            "GG": {
                "label": "GG",
                "nickname": "Variant",
                "interpretation": (
                    "Two copies of the G allele. Associated with ADHD tendencies in some "
                    "studies, although effect sizes are small. SNAP25 is involved in "
                    "synaptic vesicle exocytosis and neurotransmitter release."
                ),
                "score": 0.7,
            },
        },
        "scientificBasis": (
            "SNAP25 (Synaptosomal-Associated Protein 25kDa) rs3746544 is a candidate gene "
            "for ADHD. SNAP25 is essential for synaptic vesicle fusion and neurotransmitter "
            "release. A meta-analysis of ADHD candidate-gene studies found a significant but "
            "heterogeneous association for SNAP25 (Gizer et al. 2009). Effect sizes are small. The main classic ADHD candidate genes "
            "(DRD4, DAT1 VNTRs) are not available on consumer DNA chips. The PGC ADHD GWAS "
            "(2023, PMC10914347) identified 27 genome-wide significant loci through "
            "GWAS approaches."
        ),
        "references": [
            {"pmid": "19506906", "title": "Candidate gene studies of ADHD: a meta-analytic review.", "year": 2009},
        ],
        "actionableInsights": [
            "Attention and focus are highly trainable skills regardless of genotype",
            "Structured routines and environmental organization can support attention",
            "If you suspect attention difficulties, professional evaluation is recommended",
        ],
    },
]


# ---------------------------------------------------------------------------
# SNP DATABASE: Spiritual Sensitivity and Mediumship
# ---------------------------------------------------------------------------

SPIRITUAL_SNPS: List[Dict[str, Any]] = [
    {
        "rsid": "rs6313",
        "gene": "HTR2A",
        "name": "Serotonin 2A Receptor and Mystical Experience Intensity",
        "chromosome": "13",
        "position": 47471478,
        "effectSize": "moderate",
        "genotypes": {
            "CC": {
                "label": "C/C (102C)",
                "nickname": "Lower receptor expression",
                "interpretation": (
                    "In human cortex, C/C carriers had lower 5-HT2A mRNA and protein than T/T "
                    "carriers. In a PET study, lower 5-HT2A binding predicted higher "
                    "mystical-experience scores after psilocybin. Linking this genotype to "
                    "mystical experience is an untested inference: that study measured receptors, "
                    "not genotype."
                ),
                "score": 0.8,
            },
            "CT": {
                "label": "C/T (heterozygous)",
                "nickname": "Intermediate receptor expression",
                "interpretation": (
                    "One copy of the C allele (T102C). Expected to have intermediate 5-HT2A "
                    "expression. Any link to the intensity of mystical or altered states is an "
                    "untested inference."
                ),
                "score": 0.5,
            },
            "TT": {
                "label": "T/T (102T)",
                "nickname": "Higher receptor expression",
                "interpretation": (
                    "Two copies of the T allele. In human cortex, T/T carriers had higher 5-HT2A "
                    "mRNA and protein than C/C carriers. Because lower receptor binding (not "
                    "genotype) predicted more intense psilocybin experiences in PET research, "
                    "the genotype-to-experience link is an untested inference."
                ),
                "score": 0.3,
            },
        },
        "scientificBasis": (
            "HTR2A (serotonin 5-HT2A receptor) is the primary target of classic psychedelics "
            "and plays a central role in consciousness, perception, and mystical experiences. "
            "A 2021 PET study of 16 healthy volunteers (Stenbæk et al., PMID 33501857) found that "
            "lower pre-drug neocortical 5-HT2A binding predicted higher mystical-experience "
            "scores after psilocybin. That study measured receptors, not genotype. In human "
            "cortex, the 102C allele was expressed about 20% less than 102T, and C/C carriers had "
            "lower 5-HT2A mRNA and protein than T/T carriers (Polesskaya & Sokolov 2002), so any "
            "link between this genotype and mystical experience is an untested inference. "
            "rs6313 (T102C) is a synonymous variant in strong linkage disequilibrium "
            "with the promoter variant rs6311 (A-1438G). The 5-HT2A receptor has been "
            "hypothesised to relate to absorption."
        ),
        "references": [
            {"pmid": "33501857", "title": "Brain serotonin 2A receptor binding predicts subjective temporal and mystical effects of psilocybin in healthy humans.", "year": 2021},
            {"pmid": "11891796", "title": "Differential expression of the \"C\" and \"T\" alleles of the 5-HT2A receptor gene in the temporal cortex of normal individuals and schizophrenics.", "year": 2002},
        ],
        "actionableInsights": [
            "Contemplative practices (meditation, prayer) are open to everyone regardless of genotype",
            "Whether this variant relates to the vividness of contemplative experiences has not been tested",
            "This variant relates to consciousness and perception — not to any specific spiritual ability",
        ],
    },
    {
        "rsid": "rs6311",
        "gene": "HTR2A",
        "name": "Serotonin 2A Promoter and Transcendent Experiences",
        "chromosome": "13",
        "position": 47471867,
        "effectSize": "moderate",
        "genotypes": {
            "GG": {
                "label": "G/G (-1438G)",
                "nickname": "Lower promoter activity",
                "interpretation": (
                    "In cells expressing HTR2A, the G allele drove lower promoter activity than "
                    "the A allele. This is the same direction as rs6313 C/C (lower expression), "
                    "but any link to mystical or transcendent experience is untested."
                ),
                "score": 0.7,
            },
            "GA": {
                "label": "G/A (heterozygous)",
                "nickname": "Intermediate",
                "interpretation": (
                    "One copy of the A allele. Expected to have intermediate HTR2A promoter "
                    "activity. Any link to transcendent experience is untested."
                ),
                "score": 0.5,
            },
            "AA": {
                "label": "A/A (-1438A)",
                "nickname": "Higher promoter activity",
                "interpretation": (
                    "Two copies of the A allele. In cells expressing HTR2A, the A allele drove "
                    "greater promoter activity than G. This is the same direction as rs6313 T/T "
                    "(higher expression), but any link to mystical or transcendent experience "
                    "is untested."
                ),
                "score": 0.3,
            },
        },
        "scientificBasis": (
            "HTR2A rs6311 (A-1438G) is a promoter polymorphism in strong linkage disequilibrium "
            "with rs6313 (T102C). In cells expressing HTR2A, the A allele drove greater promoter activity than G "
            "(Parsons et al. 2004, Biol Psychiatry). The "
            "serotonin 2A receptor is the primary mediator of psychedelic and mystical experiences, and "
            "variations in this receptor system have been hypothesised to relate to individual differences in "
            "absorption, openness to experience, and self-transcendence. In older Australian twins, "
            "self-transcendence heritability was 0.37 (men) and 0.41 (women) (Kirk et al. 1999, Twin Res)."
        ),
        "references": [
            {"pmid": "15364038", "title": "The -1438A/G polymorphism in the 5-hydroxytryptamine type 2A receptor gene affects promoter activity.", "year": 2004},
            {"pmid": "10480742", "title": "Self-transcendence as a measure of spirituality in a sample of older Australian twins.", "year": 1999},
        ],
        "actionableInsights": [
            "Combined with rs6313, these variants shape your serotonergic system's sensitivity profile",
            "Contemplative traditions across cultures describe phenomena consistent with serotonergic variation",
        ],
    },
    {
        "rsid": "rs4570625",
        "gene": "TPH2",
        "name": "Serotonin Synthesis (TPH2) and Self-Transcendence (Exploratory)",
        "chromosome": "12",
        "position": 71938143,
        "effectSize": "preliminary",
        "genotypes": {
            "GG": {
                "label": "GG",
                "nickname": "Reference",
                "interpretation": (
                    "Reference genotype for this TPH2 promoter variant. No known effect on "
                    "amygdala reactivity or self-transcendence."
                ),
                "score": 0.5,
            },
            "GT": {
                "label": "GT",
                "nickname": "Intermediate",
                "interpretation": (
                    "One copy of the T allele. This TPH2 regulatory variant has been associated "
                    "with greater amygdala reactivity (T allele). No study has linked rs4570625 to "
                    "self-transcendence. Interpret with extreme caution."
                ),
                "score": 0.5,
            },
            "TT": {
                "label": "TT",
                "nickname": "Variant",
                "interpretation": (
                    "Two copies of the T allele. This TPH2 regulatory variant has been associated "
                    "with greater amygdala reactivity (Brown et al. 2005). No study has linked "
                    "rs4570625 to self-transcendence, so any link to spirituality remains "
                    "scientifically unvalidated."
                ),
                "score": 0.5,
            },
        },
        "scientificBasis": (
            "rs4570625 is a promoter variant of TPH2 (tryptophan hydroxylase 2, the G-703T variant), "
            "not a variant in VMAT2, the gene Dean Hamer popularised as the 'God Gene'. A regulatory "
            "TPH2 variant biases amygdala reactivity, with the T allele associated with greater "
            "amygdala reactivity (Brown et al. 2005). No study has linked rs4570625 to "
            "self-transcendence. Hamer's VMAT2 hypothesis has NOT been replicated in large-scale "
            "GWAS. This marker is included as exploratory context only, clearly marked as "
            "preliminary evidence."
        ),
        "references": [
            {"pmid": "16044172", "title": "A regulatory variant of the human tryptophan hydroxylase-2 gene biases amygdala reactivity.", "year": 2005},
        ],
        "actionableInsights": [
            "Self-transcendence is a validated personality dimension regardless of its genetic basis",
            "Contemplative practices (meditation, prayer, nature immersion) cultivate transcendence in everyone",
            "The 'God Gene' is a cultural concept, not an established scientific finding; this marker is not in that gene",
        ],
    },
]


# ---------------------------------------------------------------------------
# USP MEDIUMSHIP STUDY CONTEXT
# ---------------------------------------------------------------------------

USP_STUDY_CONTEXT = (
    "In January 2025, researchers affiliated with USP (University of Sao Paulo) and "
    "IPq (Institute of Psychiatry) published a landmark study in the Brazilian Journal "
    "of Psychiatry (PMID 39874024): 'Candidate genes related to spiritual mediumship: "
    "a whole-exome sequencing analysis of highly gifted mediums.' The study compared "
    "whole-exome sequencing of 54 highly experienced mediums (10+ years of practice) "
    "against 53 non-medium first-degree relatives. Key findings: 15,669 genetic variants "
    "were found exclusively in mediums, with 33 genes altered in >=1/3 "
    "of mediums but in none of their relatives. The most affected biological pathway was the "
    "inflammatory and immune system (43.9%), with ZAP-70 translocation to the immunological "
    "synapse as the main finding.\n\n"
    "IMPORTANT SCIENTIFIC CONTEXT: This study was criticized by Natalia Pasternak "
    "(Questao de Ciencia Institute) and Wellington Zangari (USP Institute of Psychology) for "
    "a fundamental design flaw — using first-degree relatives as controls inflates "
    "false positives, since they share 50% of DNA. The phenotype definition (self-reported "
    "mediumship) lacks objective biological criteria. No independent replication exists. "
    "The authors themselves state that replication is 'necessary and indispensable.' This research "
    "is preliminary and should be interpreted as an intriguing initial exploration, not as "
    "validated science."
)

NEUROIMAGING_CONTEXT = (
    "Brazilian research groups have produced significant neuroimaging studies of Spiritist "
    "mediums during trance states:\n\n"
    "* SPECT studies (PLoS ONE 2012) of experienced Brazilian Spiritist mediums showed "
    "frontal lobe deactivation during psychography (automatic writing) — consistent with "
    "intentional inhibition of executive control, not pathological dissociation.\n\n"
    "* EEG studies (SciELO 2016) found greater theta and beta power during anomalous "
    "sensory experiences vs. controls. Gamma and beta waves distinguished mediumistic "
    "communication states from ordinary mental tasks.\n\n"
    "* A consistent finding across studies: mediumistic trance is a measurable altered state of "
    "consciousness, distinct from ordinary wakefulness — not faked. Mental health studies "
    "of Brazilian Spiritist mediums consistently find below-average rates of "
    "psychiatric disorders and above-average socio-educational levels, suggesting that "
    "mediumistic practice in structured religious contexts is psychologically adaptive.\n\n"
    "* The serotonin 2A receptor (HTR2A) — analyzed in this report — is the primary "
    "mediator of altered states of consciousness in both pharmacological (psilocybin) "
    "and non-pharmacological (meditation, trance) contexts."
)

CULTURAL_CONTEXT = (
    "Brazil has a uniquely favorable cultural context for exploring the intersection "
    "between genetics and spiritual experiences:\n\n"
    "* IBGE Census 2022: 1.84% of Brazilians identify as Spiritists (~3.2 million), but "
    "a much larger population across Catholic, Umbanda, Candomble, and Evangelical traditions "
    "engages with mediumship culturally.\n\n"
    "* Brazil is widely described as the largest Spiritist country in the world, with over 12,000 "
    "Spiritist institutions.\n\n"
    "* Spiritist demographic profile: 60.6% female, 48% with university degree, 96.6% "
    "with internet access — a premium and digitally literate consumer segment.\n\n"
    "* The NUPES research center (Prof. Alexander Moreira-Almeida, UFJF) has worked for "
    "decades to legitimize the scientific study of spiritual experiences in Brazil.\n\n"
    "* The 2025 USP study generated major media coverage (CNN Brasil, Estado de Minas), "
    "demonstrating significant public appetite for content bridging science and spirituality."
)


# ---------------------------------------------------------------------------
# DISCLAIMERS
# ---------------------------------------------------------------------------

DISCLAIMERS = {
    "personality": (
        "This report provides educational information about genetic variants that scientific research "
        "has associated with personality tendencies. Personality is shaped by a complex "
        "interplay of genetic, environmental, cultural, and experiential factors. "
        "These results describe probabilistic tendencies based on population-level research "
        "— they do not determine who you are or how you will behave. Effect sizes "
        "are generally small, and genetic variants explain only a fraction of "
        "variation in personality traits. This is not a medical or psychological diagnosis."
    ),
    "mentalHealth": (
        "This report describes genetic variants that scientific research has associated with "
        "risk factors for certain mood and anxiety-related conditions. Having a "
        "variant associated with elevated risk does NOT mean you will develop any "
        "condition. Most people with these variants never develop the associated "
        "condition. This report is NOT a diagnostic tool and does not replace "
        "evaluation by a qualified mental health professional. If you have concerns "
        "about your mental health, consult a healthcare professional. This report is provided "
        "for educational purposes only and does not constitute medical advice. "
        "If you are in crisis, please contact a crisis helpline in your country."
    ),
    "spiritualSensitivity": (
        "This report explores emerging scientific research on the genetics of spiritual "
        "experiences, absorption, and transcendent states. The 2025 USP study published in the Brazilian Journal of Psychiatry "
        "(PMID 39874024) identified candidate genes potentially associated with mediumship "
        "in Spiritist practitioners — this research is preliminary and has not been independently "
        "replicated. No validated genetic test for mediumship or spiritual ability "
        "exists. This report reflects current scientific hypotheses and should be interpreted "
        "as a cultural and scientific exploration, not as a predictive or "
        "diagnostic tool. GeneHealth respects all spiritual traditions and presents this "
        "information in the spirit of scientific curiosity and cultural appreciation."
    ),
}


# ---------------------------------------------------------------------------
# ANALYSIS FUNCTIONS
# ---------------------------------------------------------------------------

def _lookup_snp(variants: Dict[str, Tuple[str, str, str]], rsid: str) -> Optional[str]:
    """Look up a SNP genotype from the variants dict.

    Returns the genotype string (e.g., 'AG') or None if not found.
    The variants dict maps rsid -> (chromosome, position, genotype).
    """
    if rsid in variants:
        val = variants[rsid]
        # Support both formats: plain genotype string or (chr, pos, genotype) tuple
        if isinstance(val, str):
            return val
        return val[2]  # genotype is the third element
    return None


def _analyze_snp_category(
    variants: Dict[str, Tuple[str, str, str]],
    snp_database: List[Dict[str, Any]]
) -> List[Dict[str, Any]]:
    """Analyze a category of SNPs against user variants."""
    results = []

    for snp_info in snp_database:
        rsid = snp_info["rsid"]
        genotype = _lookup_snp(variants, rsid)

        if genotype is None:
            # SNP not found in user data
            results.append({
                "name": snp_info["name"],
                "gene": snp_info["gene"],
                "rsid": rsid,
                "genotype": "Not available",
                "interpretation": (
                    f"This SNP ({rsid}) was not found in your DNA file. "
                    f"This may occur because your provider's genotyping chip does not "
                    f"include this variant."
                ),
                "effectSize": snp_info["effectSize"],
                "scientificBasis": snp_info["scientificBasis"],
                "references": snp_info["references"],
                "actionableInsights": [],
                "score": None,
            })
            continue

        # Normalize genotype for lookup (try both orientations)
        genotype_upper = genotype.upper().strip()
        reversed_gt = genotype_upper[::-1] if len(genotype_upper) == 2 else genotype_upper

        genotype_info = snp_info["genotypes"].get(genotype_upper)
        if genotype_info is None:
            genotype_info = snp_info["genotypes"].get(reversed_gt)

        if genotype_info is None:
            # Unknown genotype combination
            results.append({
                "name": snp_info["name"],
                "gene": snp_info["gene"],
                "rsid": rsid,
                "genotype": genotype_upper,
                "interpretation": (
                    f"Your genotype ({genotype_upper}) for {snp_info['gene']} {rsid} "
                    f"is not in our standard interpretation database. This may represent "
                    f"a rare variant or a strand orientation difference."
                ),
                "effectSize": snp_info["effectSize"],
                "scientificBasis": snp_info["scientificBasis"],
                "references": snp_info["references"],
                "actionableInsights": [],
                "score": None,
            })
            continue

        results.append({
            "name": snp_info["name"],
            "gene": snp_info["gene"],
            "rsid": rsid,
            "genotype": genotype_upper,
            "genotypeLabel": genotype_info.get("label", genotype_upper),
            "genotypeNickname": genotype_info.get("nickname", ""),
            "interpretation": genotype_info["interpretation"],
            "effectSize": snp_info["effectSize"],
            "scientificBasis": snp_info["scientificBasis"],
            "references": snp_info["references"],
            "actionableInsights": snp_info["actionableInsights"],
            "score": genotype_info.get("score"),
        })

    return results


def _generate_category_summary(traits: List[Dict[str, Any]], category: str) -> str:
    """Generate a summary paragraph for a category based on analyzed traits."""
    analyzed = [t for t in traits if t["genotype"] != "Not available" and t.get("score") is not None]

    if not analyzed:
        return f"No {category} SNPs were found in your DNA data. Your genotyping chip may not include these variants."

    avg_score = sum(t["score"] for t in analyzed) / len(analyzed)
    count = len(analyzed)
    total = len(traits)

    if category == "personality":
        if avg_score >= 0.6:
            return (
                f"Based on {count} of {total} personality variants analyzed, your genetic profile "
                f"suggests elevated tendencies for emotional sensitivity and risk-taking. "
                f"Remember that personality is shaped primarily by life experience, "
                f"culture, and personal choice — genetics provides only a small part of the picture."
            )
        elif avg_score >= 0.4:
            return (
                f"Based on {count} of {total} personality variants analyzed, your genetic profile "
                f"shows a balanced combination of traits. You carry a typical combination "
                f"of personality-related variants. Personality is predominantly "
                f"shaped by environment and experience."
            )
        else:
            return (
                f"Based on {count} of {total} personality variants analyzed, your genetic profile "
                f"suggests typical emotional regulation and moderate risk tolerance. "
                f"These are population-level tendencies — your actual personality is shaped "
                f"primarily by your unique life experiences."
            )

    elif category == "mentalHealth":
        if avg_score >= 0.6:
            return (
                f"Based on {count} of {total} stress response variants analyzed, your "
                f"genetic profile shows some variants associated with increased stress "
                f"sensitivity. This is informational only — most people with these variants "
                f"never develop mental health conditions. Proactive stress management, regular "
                f"exercise, and strong social support are beneficial for everyone. Consult a healthcare "
                f"professional for any mental health concerns."
            )
        elif avg_score >= 0.4:
            return (
                f"Based on {count} of {total} variants analyzed, your genetic stress "
                f"response profile is within typical ranges. You carry a common combination "
                f"of protective and sensitivity variants. Mental health is shaped by "
                f"many factors beyond genetics."
            )
        else:
            return (
                f"Based on {count} of {total} variants analyzed, your genetic profile "
                f"shows predominantly protective variants for stress response. This is "
                f"encouraging, but does not guarantee mental health outcomes — environment, lifestyle, "
                f"and life events play important roles."
            )

    else:  # spiritualSensitivity
        htr2a_traits = [t for t in analyzed if t["gene"] == "HTR2A"]
        htr2a_score = sum(t["score"] for t in htr2a_traits) / len(htr2a_traits) if htr2a_traits else 0.5

        if htr2a_score >= 0.6:
            return (
                f"Based on {count} of {total} variants analyzed, your serotonergic profile "
                f"suggests potentially increased sensitivity to transcendent and "
                f"contemplative experiences. Your HTR2A (serotonin 2A receptor) variants are associated "
                f"with altered receptor expression; a PET study found that lower receptor binding "
                f"predicted more intense psilocybin experiences, but the link from genotype to "
                f"experience is untested. This is a genetic tendency, not a spiritual "
                f"diagnosis. The 2025 USP study explored similar themes through a different "
                f"lens (immune system genetics in mediums), and both lines of research "
                f"point to the complex biology underlying spiritual experiences."
            )
        elif htr2a_score >= 0.4:
            return (
                f"Based on {count} of {total} variants analyzed, your serotonergic profile "
                f"shows an intermediate sensitivity pattern. Your HTR2A variants suggest "
                f"intermediate receptor expression. Spiritual "
                f"experiences are shaped by practice, culture, belief, and biology together."
            )
        else:
            return (
                f"Based on {count} of {total} variants analyzed, your serotonergic profile "
                f"shows a standard sensitivity pattern. Spiritual and contemplative experiences "
                f"are accessible to everyone regardless of genetic profile — practice and "
                f"intention are the primary drivers."
            )


def _generate_overall_profile(traits: List[Dict[str, Any]], category: str) -> str:
    """Generate a one-line overall profile label."""
    analyzed = [t for t in traits if t.get("score") is not None]
    if not analyzed:
        return "Insufficient data"

    avg = sum(t["score"] for t in analyzed) / len(analyzed)

    profiles = {
        "personality": {
            0.6: "Sensitive and Exploratory",
            0.4: "Balanced and Adaptable",
            0.0: "Steady and Measured",
        },
        "mentalHealth": {
            0.6: "Increased Stress Sensitivity",
            0.4: "Typical Stress Response",
            0.0: "Resilient Profile",
        },
        "spiritualSensitivity": {
            0.6: "Elevated Transcendent Sensitivity",
            0.4: "Intermediate Sensitivity",
            0.0: "Standard Sensitivity",
        },
    }

    cat_profiles = profiles.get(category, profiles["personality"])
    for threshold in sorted(cat_profiles.keys(), reverse=True):
        if avg >= threshold:
            return cat_profiles[threshold]
    return "Typical"


# ---------------------------------------------------------------------------
# MAIN ANALYSIS FUNCTION
# ---------------------------------------------------------------------------

def analyze_mind_spirit(
    variants: Dict[str, Tuple[str, str, str]]
) -> Dict[str, Any]:
    """
    Run full Mind & Spirit analysis on user variants.

    Args:
        variants: Dict mapping rsid -> (chromosome, position, genotype)

    Returns:
        Dict with complete analysis results
    """
    # Analyze each category
    personality_traits = _analyze_snp_category(variants, PERSONALITY_SNPS)
    mental_health_traits = _analyze_snp_category(variants, MENTAL_HEALTH_SNPS)
    spiritual_traits = _analyze_snp_category(variants, SPIRITUAL_SNPS)

    # Also include OXTR rs53576 in spiritual context (absorption/empathy)
    # It's already analyzed in personality — reference it
    oxtr_in_spiritual = None
    for t in personality_traits:
        if t["rsid"] == "rs53576":
            oxtr_in_spiritual = {
                **t,
                "name": "Absorption and Empathic Sensitivity (Spiritual Context)",
                "scientificBasis": (
                    t["scientificBasis"] + " In the context of spiritual experiences, "
                    "OXTR variants have been associated with the psychological construct of "
                    "'absorption' — the tendency to become deeply immersed in experiences. "
                    "Absorption has ~50% heritability in twin studies and is the personality "
                    "trait most consistently associated with reported mediumistic and "
                    "mystical experiences across cultures."
                ),
            }
            break

    if oxtr_in_spiritual:
        spiritual_traits.append(oxtr_in_spiritual)

    return {
        "personality_traits": personality_traits,
        "mental_health_traits": mental_health_traits,
        "spiritual_traits": spiritual_traits,
    }


def generate_mind_spirit_json(result: Dict[str, Any]) -> Dict[str, Any]:
    """
    Generate the final JSON output for the Mind & Spirit report.

    Args:
        result: Output from analyze_mind_spirit()

    Returns:
        Dict matching the MindSpiritReport frontend JSON schema
    """
    personality_traits = result["personality_traits"]
    mental_health_traits = result["mental_health_traits"]
    spiritual_traits = result["spiritual_traits"]

    return {
        "reportType": "mind_spirit",
        "version": "1.0",
        "personality": {
            "summary": _generate_category_summary(personality_traits, "personality"),
            "overallProfile": _generate_overall_profile(personality_traits, "personality"),
            "traits": personality_traits,
        },
        "mentalHealth": {
            "summary": _generate_category_summary(mental_health_traits, "mentalHealth"),
            "overallProfile": _generate_overall_profile(mental_health_traits, "mentalHealth"),
            "disclaimer": DISCLAIMERS["mentalHealth"],
            "traits": mental_health_traits,
        },
        "spiritualSensitivity": {
            "summary": _generate_category_summary(spiritual_traits, "spiritualSensitivity"),
            "overallProfile": _generate_overall_profile(spiritual_traits, "spiritualSensitivity"),
            "uspStudyContext": USP_STUDY_CONTEXT,
            "neuroimagingContext": NEUROIMAGING_CONTEXT,
            "culturalContext": CULTURAL_CONTEXT,
            "traits": spiritual_traits,
        },
        "disclaimers": DISCLAIMERS,
    }
