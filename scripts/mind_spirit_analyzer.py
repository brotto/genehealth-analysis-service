"""
Mind & Spirit Analyzer
Analyzes genetic variants associated with personality traits, mental health
predispositions, and spiritual sensitivity (including mediumship genetics).

THREE CATEGORIES:
  1. Personality Traits - COMT, CADM2, MSRA, OPRM1, OXTR
  2. Mental Health & Stress Response - FKBP5, BDNF, MTHFR, CACNA1C, RGS2, SNAP25
  3. Spiritual Sensitivity - HTR2A, VMAT2/SLC18A2, OXTR (absorption context)

SCIENTIFIC BASIS:
  - Personality: GWAS-validated SNPs (Nature Human Behaviour 2024, Translational Psychiatry 2023)
  - Mental Health: PGC consortium studies (Cell 2024, Nature Genetics 2025)
  - Spiritual: USP/BJPsych 2025 mediumship study (PMID 39874024), HTR2A mystical experience research

IMPORTANT CAVEATS:
  - All associations are probabilistic, not deterministic
  - Effect sizes are generally small (each SNP explains <1% of trait variance)
  - Environment, culture, and personal choice strongly modify any genetic tendency
  - Mental health reports are educational only, not diagnostic
  - Spiritual sensitivity findings are preliminary and exploratory

References:
  - PGC MDD Working Group 2024, Cell (PMC11092713) - 697 depression loci
  - Anxiety GWAS 2025, Nature Genetics - 58 anxiety loci
  - Wagner et al. 2025, BJPsych (PMID 39874024) - mediumship candidate genes
  - Pasternak & Zangari 2025, Questao de Ciencia - critical review
  - Singh et al. 2023, Translational Psychiatry - CADM2 impulsivity
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
            "The Val allele (G) produces a high-activity enzyme that clears dopamine ~4x faster "
            "than the Met allele (A). Meta-analyses confirm small but real effects on stress "
            "response (d ~ 0.15-0.25) and cognitive flexibility. The 'Warrior vs Worrier' model "
            "(Goldman et al. 2005) captures the trade-off: Val/Val excels under stress but has lower "
            "baseline prefrontal dopamine, while Met/Met excels in baseline cognition but is more "
            "sensitive to stress."
        ),
        "references": [
            {"pmid": "16151010", "title": "The 'warrior' and 'worrier' model for COMT", "year": 2005},
            {"pmid": "29520078", "title": "COMT Val158Met and cognition meta-analysis", "year": 2018},
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
            "Effect sizes are small but consistent (Singh et al. 2023, Translational "
            "Psychiatry). A 2023 Nature Human Behaviour study identified 27+ independent loci "
            "for risk tolerance, with CADM2 among the strongest."
        ),
        "references": [
            {"pmid": "36658148", "title": "CADM2 implicated in impulsive personality by GWAS", "year": 2023},
            {"pmid": "30643258", "title": "Risk tolerance GWAS in UK Biobank", "year": 2019},
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
            "regulation pathways. Replicated in Davies et al. 2017 (PMC5537199)."
        ),
        "references": [
            {"pmid": "28696412", "title": "Replication of CADM2 and MSRA on human behavior", "year": 2017},
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
                    "Reference genotype for the mu-opioid receptor. Standard pain sensitivity "
                    "and typical placebo response. Most common genotype worldwide."
                ),
                "score": 0.4,
            },
            "AG": {
                "label": "Asp/Asn",
                "nickname": "Increased sensitivity",
                "interpretation": (
                    "One copy of the G allele (A118G). Associated with altered pain sensitivity, "
                    "reduced placebo response, and potentially greater frustration sensitivity. "
                    "May experience social rejection more intensely."
                ),
                "score": 0.6,
            },
            "GG": {
                "label": "Asn/Asn",
                "nickname": "High sensitivity",
                "interpretation": (
                    "Two copies of the variant allele. Greater emotional sensitivity to social "
                    "interactions and potentially altered reward processing. Less responsive "
                    "to placebo effects and may require different approaches for pain management."
                ),
                "score": 0.8,
            },
        },
        "scientificBasis": (
            "OPRM1 A118G (rs1799971) affects the mu-opioid receptor, which mediates pain, reward, and "
            "social bonding. The G allele reduces receptor expression and has been associated with greater "
            "frustration sensitivity, reduced placebo response, and altered social "
            "pain processing. Multiple neuroimaging studies show that the G allele is associated with "
            "greater neural responses to social rejection (Way et al. 2009, PNAS)."
        ),
        "references": [
            {"pmid": "19622738", "title": "OPRM1 A118G and social rejection sensitivity", "year": 2009},
            {"pmid": "15583718", "title": "Functional effects of the OPRM1 A118G polymorphism", "year": 2005},
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
                    "in some studies. However, large-scale replications in UK Biobank "
                    "(n>100,000) found weak and inconsistent effects."
                ),
                "score": 0.3,
            },
        },
        "scientificBasis": (
            "OXTR rs53576 is the most studied oxytocin receptor variant for social behavior. "
            "Early studies linked the GG genotype to higher empathy and prosocial behavior. "
            "However, a 2017 meta-analysis (n>10,000) found the association weak and "
            "inconsistent (PMID 28343138). UK Biobank analyses show minimal effect. This SNP "
            "is included for completeness, but results should be interpreted with caution."
        ),
        "references": [
            {"pmid": "28343138", "title": "Meta-analysis of OXTR rs53576 and empathy", "year": 2017},
            {"pmid": "21775986", "title": "OXTR and social behavior review", "year": 2011},
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
                    "Two copies of the A allele. Associated with altered attachment styles "
                    "in some studies, particularly in East Asian populations. Limited "
                    "replication in European-ancestry cohorts."
                ),
                "score": 0.6,
            },
        },
        "scientificBasis": (
            "OXTR rs2254298 has been studied for associations with social bonding, attachment security, "
            "and prosocial behavior. Results vary significantly by population — the "
            "A allele shows stronger effects in East Asian population studies than in "
            "European-ancestry samples. Limited predictive value for individual social behavior."
        ),
        "references": [
            {"pmid": "21775986", "title": "OXTR polymorphisms and social behavior", "year": 2011},
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
            "FKBP5 rs1360780 is one of the most robust single-SNP findings in stress-related "
            "psychiatry. The T allele creates an allele-specific DNA methylation response to stress "
            "that increases FKBP5 protein expression, which binds and inhibits the glucocorticoid receptor, "
            "disrupting cortisol feedback. A 2025 MDPI Genes review (Zannas et al.) "
            "characterized behavioral phenotypes. The gene-by-environment interaction with childhood "
            "maltreatment is well-replicated across multiple independent cohorts. Effect sizes are "
            "modest but mechanistically compelling."
        ),
        "references": [
            {"pmid": "39856776", "title": "FKBP5 rs1360780: genetic variation and behavioral phenotypes", "year": 2025},
            {"pmid": "24029109", "title": "FKBP5 epigenetics and stress vulnerability", "year": 2013},
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
                    "particularly in males and elderly populations."
                ),
                "score": 0.7,
            },
        },
        "scientificBasis": (
            "BDNF Val66Met (rs6265) affects activity-dependent secretion of brain-derived neurotrophic "
            "factor, crucial for synaptic plasticity. A 2023 Frontiers in Psychiatry meta-analysis "
            "found near-null effect for MDD (OR ~0.96). However, the "
            "Met allele shows stronger associations in males (OR 1.27) and geriatric depression "
            "(OR ~1.48). The greatest value lies in neuroplasticity and exercise biology — BDNF "
            "mediates many exercise-related cognitive benefits."
        ),
        "references": [
            {"pmid": "37398592", "title": "BDNF Val66Met and depression meta-analysis", "year": 2023},
            {"pmid": "15457404", "title": "BDNF Val66Met and hippocampal function", "year": 2004},
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
                "label": "CT (heterozigoto)",
                "nickname": "Reduced activity (~65%)",
                "interpretation": (
                    "One copy of the T allele (C677T). MTHFR enzyme activity reduced to "
                    "approximately 65% of normal. Generally clinically insignificant, but "
                    "adequate folate intake is advisable."
                ),
                "score": 0.5,
            },
            "AA": {
                "label": "TT (homozigoto)",
                "nickname": "Significantly reduced (~30%)",
                "interpretation": (
                    "Two copies of the T allele. MTHFR enzyme activity reduced to approximately "
                    "30% of normal. Meta-analysis of 30 studies finds OR = 1.20 for depression "
                    "(TT vs CC). Effect concentrated in Asian populations. TT prevalence in "
                    "Brazil is approximately 10%. Actionable value is high: methylfolate "
                    "supplementation for TT carriers is widely accepted clinically."
                ),
                "score": 0.8,
            },
        },
        "scientificBasis": (
            "MTHFR C677T (rs1801133) reduces methylenetetrahydrofolate reductase activity, "
            "affecting the one-carbon metabolism cycle that produces SAMe (the primary methyl "
            "donor for neurotransmitter synthesis). TT homozygotes have ~30% enzyme activity "
            "and may develop mild hyperhomocysteinemia. Meta-analysis shows modest association "
            "with depression (OR = 1.20, p=0.0004) concentrated in Asian populations. The greatest "
            "clinical value is actionability: folate/methylfolate supplementation."
        ),
        "references": [
            {"pmid": "23212058", "title": "MTHFR C677T and depression meta-analysis", "year": 2013},
            {"pmid": "15166018", "title": "MTHFR and psychiatric disorders review", "year": 2004},
        ],
        "actionableInsights": [
            "TT carriers: Consider L-methylfolate supplementation (consult a healthcare professional)",
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
                    "expression. Modestly elevated transdiagnostic risk for mood conditions "
                    "(bipolar, depression, schizophrenia) — OR approximately 1.07."
                ),
                "score": 0.5,
            },
            "AA": {
                "label": "AA",
                "nickname": "Two risk alleles",
                "interpretation": (
                    "Two copies of the A allele. The most consistently replicated genetic "
                    "risk variant for bipolar disorder (OR ~1.14). This is a transdiagnostic variant "
                    "— it confers small risk increases across multiple mood and "
                    "psychotic conditions. Most carriers never develop any condition."
                ),
                "score": 0.7,
            },
        },
        "scientificBasis": (
            "CACNA1C rs1006737 is the most consistently replicated bipolar disorder risk SNP "
            "across all PGC meta-analyses (OR ~1.14, p = 9.78x10^-10). It encodes the alpha-1C subunit "
            "of the voltage-dependent L-type calcium channel, critical for neuronal excitability "
            "and synaptic plasticity. Important: this is a transdiagnostic variant — it also increases "
            "risk for schizophrenia and recurrent depression. The 2025 PGC4 trans-ancestral study "
            "found 93 bipolar loci (23 new), with CACNA1C remaining among the strongest."
        ),
        "references": [
            {"pmid": "34002096", "title": "PGC Bipolar GWAS identifies 64 loci", "year": 2021},
            {"pmid": "21926972", "title": "CACNA1C and psychiatric phenotypes", "year": 2011},
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
                "nickname": "Reference",
                "interpretation": (
                    "Reference genotype. Standard RGS2 expression and G-protein signaling "
                    "regulation. Typical anxiety response levels."
                ),
                "score": 0.3,
            },
            "CG": {
                "label": "CG",
                "nickname": "Intermediate",
                "interpretation": (
                    "One copy of the G allele. Potentially reduced RGS2 expression, which "
                    "may subtly increase stress hormone signaling downstream "
                    "of CRH receptors."
                ),
                "score": 0.5,
            },
            "GG": {
                "label": "GG",
                "nickname": "Elevated anxiety sensitivity",
                "interpretation": (
                    "Two copies of the G allele. Associated with higher anxiety sensitivity "
                    "and stress reactivity in multiple studies. RGS2 dampens G-protein "
                    "signaling — reduced RGS2 may lead to increased stress hormone "
                    "response."
                ),
                "score": 0.7,
            },
        },
        "scientificBasis": (
            "RGS2 (Regulator of G-protein Signaling 2) rs4606 modulates signaling downstream "
            "of corticotropin-releasing hormone (CRH) receptors. The G allele has been associated "
            "with anxiety and stress reactivity in multiple independent studies (2008-2024). "
            "Effect sizes are small but the finding has been replicated. A landmark 2025 anxiety "
            "GWAS from Nature Genetics (n=122,341) identified 58 loci, further supporting "
            "the involvement of GABAergic and stress signaling pathways."
        ),
        "references": [
            {"pmid": "18266781", "title": "RGS2 and anxiety disorders", "year": 2008},
        ],
        "actionableInsights": [
            "GG carriers: Structured relaxation techniques (progressive muscle relaxation, deep breathing) may be especially helpful",
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
            "release. Effect sizes are small. The main classic ADHD candidate genes "
            "(DRD4, DAT1 VNTRs) are not available on consumer DNA chips. The PGC ADHD GWAS "
            "(2023, PMC10914347) identified 27 genome-wide significant loci through "
            "GWAS approaches."
        ),
        "references": [
            {"pmid": "15742474", "title": "SNAP25 and ADHD association", "year": 2005},
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
                "nickname": "Higher receptor density",
                "interpretation": (
                    "Reference genotype associated with higher 5-HT2A receptor density in the "
                    "cortex. Standard intensity of transcendent or mystical experiences. "
                    "This is the most common genotype."
                ),
                "score": 0.3,
            },
            "CT": {
                "label": "C/T (heterozigoto)",
                "nickname": "Intermediate receptor density",
                "interpretation": (
                    "One copy of the T allele (T102C). Intermediate 5-HT2A receptor density. "
                    "Potentially moderate sensitivity to transcendent or altered states of "
                    "consciousness. Psilocybin research shows intermediate intensity of "
                    "mystical experiences with this genotype."
                ),
                "score": 0.5,
            },
            "TT": {
                "label": "T/T (102T)",
                "nickname": "Lower receptor density",
                "interpretation": (
                    "Two copies of the T allele. Associated with lower 5-HT2A receptor density "
                    "in the cortex. Paradoxically, studies show that individuals with fewer "
                    "5-HT2A receptors report MORE intense mystical experiences during psilocybin "
                    "administration. This suggests a compensatory sensitivity mechanism — fewer "
                    "receptors, but greater individual sensitivity per receptor."
                ),
                "score": 0.8,
            },
        },
        "scientificBasis": (
            "HTR2A (serotonin 5-HT2A receptor) is the primary target of classic psychedelics "
            "and plays a central role in consciousness, perception, and mystical experiences. "
            "A 2021 study (PMID 33501857) found that brain 5-HT2A receptor binding "
            "predicts the intensity of mystical effects during psilocybin administration — "
            "people with fewer receptors (TT genotype) paradoxically report MORE intense "
            "experiences. rs6313 (T102C) is a synonymous variant in strong linkage disequilibrium "
            "with the promoter variant rs6311 (A-1438G). The 5-HT2A receptor is also implicated in "
            "absorption — the psychological trait most consistently associated with reported spiritual "
            "and mediumistic experiences."
        ),
        "references": [
            {"pmid": "33501857", "title": "Brain 5-HT2A receptor binding predicts mystical experiences", "year": 2021},
            {"pmid": "17606772", "title": "HTR2A polymorphisms and hallucinogen response", "year": 2007},
        ],
        "actionableInsights": [
            "TT carriers: You may have increased sensitivity to contemplative practices (meditation, prayer)",
            "Mindfulness and meditation practices may produce particularly vivid experiences for TT carriers",
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
                "label": "G/G (A-1438)",
                "nickname": "Reference promoter",
                "interpretation": (
                    "Reference genotype for the HTR2A promoter. Standard levels of "
                    "5-HT2A receptor transcription and expression."
                ),
                "score": 0.3,
            },
            "GA": {
                "label": "G/A (heterozigoto)",
                "nickname": "Intermediate expression",
                "interpretation": (
                    "One copy of the A allele. Altered promoter activity affecting 5-HT2A "
                    "expression. Intermediate phenotype for serotonergic sensitivity."
                ),
                "score": 0.5,
            },
            "AA": {
                "label": "A/A (-1438G>A)",
                "nickname": "Altered expression",
                "interpretation": (
                    "Two copies of the A allele. This promoter variant alters HTR2A "
                    "transcription. In linkage disequilibrium with rs6313 — together, they modulate "
                    "the serotonergic system's role in consciousness and perception. Associated with "
                    "altered response to serotonergic stimuli."
                ),
                "score": 0.7,
            },
        },
        "scientificBasis": (
            "HTR2A rs6311 (A-1438G) is a promoter polymorphism in strong linkage disequilibrium "
            "with rs6313 (T102C). Together, they modulate 5-HT2A receptor density and function. The "
            "serotonin 2A receptor is the primary mediator of psychedelic and mystical experiences, and "
            "variations in this receptor system have been associated with individual differences in absorption, "
            "openness to experience, and self-transcendence. Twin studies show 40-50% "
            "heritability for self-transcendence (Bouchard et al.)."
        ),
        "references": [
            {"pmid": "11779785", "title": "HTR2A promoter polymorphism and receptor density", "year": 2002},
        ],
        "actionableInsights": [
            "Combined with rs6313, these variants shape your serotonergic system's sensitivity profile",
            "Contemplative traditions across cultures describe phenomena consistent with serotonergic variation",
        ],
    },
    {
        "rsid": "rs4570625",
        "gene": "SLC18A2",
        "name": "Vesicular Monoamine Transport and Self-Transcendence",
        "chromosome": "10",
        "position": 119003566,
        "effectSize": "preliminary",
        "genotypes": {
            "GG": {
                "label": "GG",
                "nickname": "Reference",
                "interpretation": (
                    "Reference genotype for the VMAT2 region. Standard vesicular monoamine "
                    "transporter function and neurotransmitter packaging."
                ),
                "score": 0.3,
            },
            "GT": {
                "label": "GT",
                "nickname": "Intermediate",
                "interpretation": (
                    "One copy of the T allele. In Dean Hamer's 'God Gene' hypothesis (2004), "
                    "the T allele was associated with higher self-transcendence scores on Cloninger's "
                    "TCI personality inventory. This finding has NOT been replicated in large-scale "
                    "GWAS. Interpret with extreme caution."
                ),
                "score": 0.5,
            },
            "TT": {
                "label": "TT",
                "nickname": "Variant",
                "interpretation": (
                    "Two copies of the T allele. Hamer's original study (n~1,000) suggested "
                    "greater self-transcendence. However, the 'God Gene' hypothesis has NOT "
                    "been replicated in large-scale studies. VMAT2 is a real monoamine transporter "
                    "involved in packaging dopamine, serotonin, and norepinephrine "
                    "into vesicles — but linking it specifically to spirituality remains "
                    "scientifically unvalidated."
                ),
                "score": 0.6,
            },
        },
        "scientificBasis": (
            "VMAT2 (SLC18A2) was proposed as the 'God Gene' by Dean Hamer in 2004. The variant "
            "rs4570625 was associated with self-transcendence scores in a small sample "
            "(n~1,000). VMAT2 packages monoamines (dopamine, serotonin, norepinephrine) into synaptic "
            "vesicles — it has a genuine role in neurotransmission. However, the 'God Gene' "
            "hypothesis has NOT been replicated in large-scale GWAS. Scientific consensus considers this "
            "speculative. We include it for completeness and cultural interest, clearly marked as "
            "preliminary evidence."
        ),
        "references": [
            {"pmid": "15520385", "title": "Hamer: The God Gene (book review)", "year": 2004},
        ],
        "actionableInsights": [
            "Self-transcendence is a validated personality dimension regardless of its genetic basis",
            "Contemplative practices (meditation, prayer, nature immersion) cultivate transcendence in everyone",
            "The 'God Gene' is a cultural concept, not an established scientific finding",
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
        "experiences, absorption, and transcendent states. The 2025 USP/BJPsych study "
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
                f"with altered receptor density — research shows this correlates with the "
                f"intensity of mystical experiences. This is a genetic tendency, not a spiritual "
                f"diagnosis. The 2025 USP study explored similar themes through a different "
                f"lens (immune system genetics in mediums), and both lines of research "
                f"point to the complex biology underlying spiritual experiences."
            )
        elif htr2a_score >= 0.4:
            return (
                f"Based on {count} of {total} variants analyzed, your serotonergic profile "
                f"shows an intermediate sensitivity pattern. Your HTR2A variants suggest "
                f"typical receptor density and sensitivity to altered states. Spiritual "
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
