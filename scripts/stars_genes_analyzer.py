"""
Stars & Genes Analyzer
Compares astrological natal chart predictions against genetic personality markers.

INPUTS:
  - variants: Dict[str, Tuple[str, str, str]]  (rsid -> (chrom, pos, genotype))
  - birth_data: Dict with date, time, lat, lng, tz

OUTPUTS:
  - Structured JSON with natal chart + DNA comparison per trait category

ASTROLOGY ENGINE: pyswisseph (Swiss Ephemeris)
GENETIC DATA: Cross-references mind_spirit, traits, dream_sleep SNPs already analyzed

References:
  - Swiss Ephemeris: astro.com/swisseph (NASA JPL DE431)
  - COMT rs4680: Egan et al. 2001, PNAS (Warrior/Worrier)
  - OXTR rs53576: Rodrigues et al. 2009, PNAS (Empathy)
  - CADM2 rs17518584: Day et al. 2016, Nature Genetics (Risk-taking)
  - DRD4 rs1800955: Ebstein et al. 1996, Nature Genetics (Novelty seeking)
"""

from typing import Dict, List, Tuple, Any, Optional
from datetime import datetime, timezone
import math

try:
    import swisseph as swe
    HAS_SWE = True
except ImportError:
    HAS_SWE = False


# ─────────────────────────────────────────────────────────────────────────────
# ZODIAC CONSTANTS
# ─────────────────────────────────────────────────────────────────────────────

SIGNS = [
    "Aries", "Taurus", "Gemini", "Cancer", "Leo", "Virgo",
    "Libra", "Scorpio", "Sagittarius", "Capricorn", "Aquarius", "Pisces",
]

SIGN_ELEMENTS = {
    "Aries": "fire", "Leo": "fire", "Sagittarius": "fire",
    "Taurus": "earth", "Virgo": "earth", "Capricorn": "earth",
    "Gemini": "air", "Libra": "air", "Aquarius": "air",
    "Cancer": "water", "Scorpio": "water", "Pisces": "water",
}

SIGN_MODALITIES = {
    "Aries": "cardinal", "Cancer": "cardinal", "Libra": "cardinal", "Capricorn": "cardinal",
    "Taurus": "fixed", "Leo": "fixed", "Scorpio": "fixed", "Aquarius": "fixed",
    "Gemini": "mutable", "Virgo": "mutable", "Sagittarius": "mutable", "Pisces": "mutable",
}

SIGN_SYMBOLS = {
    "Aries": "\u2648", "Taurus": "\u2649", "Gemini": "\u264A", "Cancer": "\u264B",
    "Leo": "\u264C", "Virgo": "\u264D", "Libra": "\u264E", "Scorpio": "\u264F",
    "Sagittarius": "\u2650", "Capricorn": "\u2651", "Aquarius": "\u2652", "Pisces": "\u2653",
}

PLANET_IDS = {
    "sun": 0,       # swe.SUN
    "moon": 1,      # swe.MOON
    "mercury": 2,   # swe.MERCURY
    "venus": 3,     # swe.VENUS
    "mars": 4,      # swe.MARS
    "jupiter": 5,   # swe.JUPITER
    "saturn": 6,    # swe.SATURN
}


# ─────────────────────────────────────────────────────────────────────────────
# GENETIC MARKERS FOR PERSONALITY DIMENSIONS
# ─────────────────────────────────────────────────────────────────────────────
# Each marker maps to a personality dimension with interpretations per genotype

GENETIC_MARKERS: Dict[str, Dict[str, Any]] = {
    # COMT — Warrior/Worrier (stress resilience vs baseline cognition)
    "rs4680": {
        "gene": "COMT",
        "dimension": "stress_resilience",
        "label": "Stress Response (Warrior/Worrier)",
        "interpretations": {
            "GG": {"label": "Warrior (Val/Val)", "score": 1.0, "desc": "Higher stress resilience, fast dopamine clearance — thrives under pressure"},
            "AG": {"label": "Balanced (Val/Met)", "score": 0.5, "desc": "Mixed stress response — adaptable to both calm and intense situations"},
            "AA": {"label": "Worrier (Met/Met)", "score": 0.0, "desc": "Higher baseline cognition, slower dopamine clearance — excels in calm, predictable settings"},
        },
    },
    # CADM2 — Risk tolerance / impulsivity
    "rs17518584": {
        "gene": "CADM2",
        "dimension": "risk_tolerance",
        "label": "Risk Tolerance",
        "interpretations": {
            "TT": {"label": "Risk-taker", "score": 1.0, "desc": "Higher genetic risk tolerance and novelty-seeking tendency"},
            "AT": {"label": "Moderate risk tolerance", "score": 0.5, "desc": "Balanced approach to risk — situationally adventurous"},
            "AA": {"label": "Risk-averse", "score": 0.0, "desc": "Genetically inclined toward caution and careful decision-making"},
        },
    },
    # OPRM1 — Emotional sensitivity / social bonding
    "rs1799971": {
        "gene": "OPRM1",
        "dimension": "emotional_sensitivity",
        "label": "Emotional Sensitivity",
        "interpretations": {
            "AA": {"label": "Standard sensitivity", "score": 0.3, "desc": "Typical emotional and social pain sensitivity"},
            "AG": {"label": "Enhanced sensitivity", "score": 0.7, "desc": "Heightened sensitivity to social rejection and emotional experiences"},
            "GG": {"label": "High sensitivity", "score": 1.0, "desc": "Markedly elevated emotional sensitivity — deeply affected by social bonds"},
        },
    },
    # OXTR — Empathy and social bonding
    "rs53576": {
        "gene": "OXTR",
        "dimension": "empathy",
        "label": "Empathic Tendency",
        "interpretations": {
            "GG": {"label": "High empathy", "score": 1.0, "desc": "Strongly empathic — sensitive to others' emotional states"},
            "AG": {"label": "Moderate empathy", "score": 0.5, "desc": "Balanced empathic response"},
            "AA": {"label": "Lower empathy reactivity", "score": 0.0, "desc": "More analytical approach to social situations — less emotionally reactive"},
        },
    },
    # FKBP5 — Stress reactivity / cortisol regulation
    "rs1360780": {
        "gene": "FKBP5",
        "dimension": "stress_reactivity",
        "label": "Stress Reactivity",
        "interpretations": {
            "CC": {"label": "Resilient", "score": 0.0, "desc": "Efficient cortisol regulation — lower stress reactivity"},
            "CT": {"label": "Moderate reactivity", "score": 0.5, "desc": "Moderate stress response amplification"},
            "TT": {"label": "High reactivity", "score": 1.0, "desc": "Amplified stress response — stronger cortisol reactions to adversity"},
        },
    },
    # BDNF — Neuroplasticity / mood
    "rs6265": {
        "gene": "BDNF",
        "dimension": "neuroplasticity",
        "label": "Brain Plasticity",
        "interpretations": {
            "CC": {"label": "High plasticity (Val/Val)", "score": 1.0, "desc": "Optimal BDNF secretion — robust neuroplasticity and mood regulation"},
            "CT": {"label": "Moderate plasticity (Val/Met)", "score": 0.5, "desc": "Slightly reduced neuroplasticity — still within normal range"},
            "TT": {"label": "Reduced plasticity (Met/Met)", "score": 0.0, "desc": "Lower BDNF secretion — may benefit from extra exercise and environmental enrichment"},
        },
    },
    # HTR2A — Spiritual sensitivity / dream vividness
    "rs6313": {
        "gene": "HTR2A",
        "dimension": "spiritual_sensitivity",
        "label": "Spiritual / Mystical Sensitivity",
        "interpretations": {
            "TT": {"label": "High sensitivity", "score": 1.0, "desc": "Enhanced serotonin 2A receptor activity — higher mystical/spiritual experience tendency"},
            "CT": {"label": "Moderate sensitivity", "score": 0.5, "desc": "Moderate serotonin receptor sensitivity"},
            "CC": {"label": "Standard sensitivity", "score": 0.0, "desc": "Typical serotonin 2A function — more grounded, less prone to mystical states"},
        },
    },
    # DRD4 — Novelty seeking
    "rs1800955": {
        "gene": "DRD4",
        "dimension": "novelty_seeking",
        "label": "Novelty Seeking",
        "interpretations": {
            "CC": {"label": "High novelty seeker", "score": 1.0, "desc": "Enhanced dopamine D4 receptor activity — drawn to new experiences and exploration"},
            "CT": {"label": "Moderate novelty seeking", "score": 0.5, "desc": "Balanced curiosity — enjoys new experiences but also values routine"},
            "TT": {"label": "Routine-preferring", "score": 0.0, "desc": "Genetically inclined toward familiar environments and established routines"},
        },
    },
    # CLOCK — Chronotype (morning/evening)
    "rs1801260": {
        "gene": "CLOCK",
        "dimension": "chronotype",
        "label": "Chronotype (Morning/Evening)",
        "interpretations": {
            "TT": {"label": "Night owl tendency", "score": 1.0, "desc": "CLOCK variant associated with evening preference and later sleep onset"},
            "TC": {"label": "Intermediate", "score": 0.5, "desc": "No strong genetic preference — chronotype shaped mainly by habits"},
            "CC": {"label": "Morning lark tendency", "score": 0.0, "desc": "CLOCK variant associated with earlier wake times and morning alertness"},
        },
    },
    # RGS2 — Anxiety susceptibility
    "rs4606": {
        "gene": "RGS2",
        "dimension": "anxiety",
        "label": "Anxiety Sensitivity",
        "interpretations": {
            "CC": {"label": "Higher anxiety tendency", "score": 1.0, "desc": "RGS2 variant linked to increased anxiety susceptibility and panic sensitivity"},
            "CG": {"label": "Moderate anxiety tendency", "score": 0.5, "desc": "Intermediate anxiety susceptibility"},
            "GG": {"label": "Lower anxiety tendency", "score": 0.0, "desc": "RGS2 variant associated with lower anxiety susceptibility"},
        },
    },
    # SNAP25 — Cognitive performance
    "rs363050": {
        "gene": "SNAP25",
        "dimension": "cognitive_performance",
        "label": "Cognitive Performance",
        "interpretations": {
            "GG": {"label": "Enhanced cognition", "score": 1.0, "desc": "SNAP25 variant associated with higher IQ and cognitive performance"},
            "AG": {"label": "Typical cognition", "score": 0.5, "desc": "Standard cognitive performance"},
            "AA": {"label": "Standard cognition", "score": 0.3, "desc": "Standard cognitive performance — no deficit, just population average"},
        },
    },
    # KIBRA — Memory
    "rs17070145": {
        "gene": "KIBRA",
        "dimension": "memory",
        "label": "Episodic Memory",
        "interpretations": {
            "TT": {"label": "Enhanced memory", "score": 1.0, "desc": "KIBRA variant associated with superior episodic memory formation"},
            "CT": {"label": "Good memory", "score": 0.6, "desc": "Moderately enhanced memory capacity"},
            "CC": {"label": "Standard memory", "score": 0.3, "desc": "Population-average memory capacity"},
        },
    },
    # SLC18A2 (VMAT2) — Self-transcendence / "God Gene"
    "rs1390938": {
        "gene": "SLC18A2",
        "dimension": "transcendence",
        "label": "Self-Transcendence",
        "interpretations": {
            "CC": {"label": "Higher transcendence", "score": 1.0, "desc": "VMAT2 variant associated with higher self-transcendence and spiritual inclination"},
            "CT": {"label": "Moderate transcendence", "score": 0.5, "desc": "Intermediate self-transcendence tendency"},
            "TT": {"label": "Lower transcendence", "score": 0.0, "desc": "More materialist/pragmatic orientation — lower self-transcendence scores"},
        },
    },
    # MSRA — Irritability / behavioral inhibition
    "rs4925638": {
        "gene": "MSRA",
        "dimension": "impulsivity",
        "label": "Behavioral Control",
        "interpretations": {
            "TT": {"label": "Higher impulsivity", "score": 1.0, "desc": "MSRA variant linked to lower behavioral inhibition and higher irritability"},
            "CT": {"label": "Moderate control", "score": 0.5, "desc": "Balanced behavioral control"},
            "CC": {"label": "Strong control", "score": 0.0, "desc": "Strong behavioral inhibition — naturally measured and controlled"},
        },
    },
}


# ─────────────────────────────────────────────────────────────────────────────
# ZODIAC SIGN PERSONALITY PROFILES
# ─────────────────────────────────────────────────────────────────────────────
# Maps each sign to trait expectations with corresponding genetic dimensions

SUN_SIGN_PROFILES: Dict[str, Dict[str, Any]] = {
    "Aries": {
        "element": "fire",
        "profile": "Bold, competitive, impulsive, and fearlessly direct. Aries charges headfirst into challenges with unstoppable energy and infectious confidence.",
        "traits": [
            {"trait": "Risk Tolerance", "prediction": "Aries is the warrior of the zodiac — fearless, competitive, and drawn to bold action.", "dimension": "risk_tolerance", "expect_high": True},
            {"trait": "Stress Resilience", "prediction": "Aries thrives under pressure, using stress as fuel rather than succumbing to it.", "dimension": "stress_resilience", "expect_high": True},
            {"trait": "Behavioral Control", "prediction": "Known for impulsivity — Aries acts first and thinks later.", "dimension": "impulsivity", "expect_high": True},
            {"trait": "Novelty Seeking", "prediction": "Constantly seeking new adventures and challenges — boredom is the enemy.", "dimension": "novelty_seeking", "expect_high": True},
        ],
    },
    "Taurus": {
        "element": "earth",
        "profile": "Steady, sensual, and deeply grounded. Taurus values stability, comfort, and the finer things in life — but watch out for legendary stubbornness.",
        "traits": [
            {"trait": "Behavioral Control", "prediction": "Taurus is patient and measured — rarely acts on impulse.", "dimension": "impulsivity", "expect_high": False},
            {"trait": "Stress Reactivity", "prediction": "Earth signs absorb stress slowly — Taurus is emotionally stable and hard to rattle.", "dimension": "stress_reactivity", "expect_high": False},
            {"trait": "Novelty Seeking", "prediction": "Taurus prefers familiar comforts over risky new experiences.", "dimension": "novelty_seeking", "expect_high": False},
            {"trait": "Brain Plasticity", "prediction": "Steady and persistent — Taurus builds knowledge gradually through repetition.", "dimension": "neuroplasticity", "expect_high": True},
        ],
    },
    "Gemini": {
        "element": "air",
        "profile": "Curious, adaptable, and endlessly communicative. Gemini's dual nature makes them the social chameleon of the zodiac — quick-witted and always learning.",
        "traits": [
            {"trait": "Novelty Seeking", "prediction": "Gemini is the most curious sign — constantly seeking new information and experiences.", "dimension": "novelty_seeking", "expect_high": True},
            {"trait": "Cognitive Performance", "prediction": "Mercury-ruled Gemini excels in mental agility and information processing.", "dimension": "cognitive_performance", "expect_high": True},
            {"trait": "Behavioral Control", "prediction": "Easily bored and restless — Gemini can be scattered and impulsive.", "dimension": "impulsivity", "expect_high": True},
            {"trait": "Empathy", "prediction": "Gemini connects intellectually more than emotionally — empathy is cerebral.", "dimension": "empathy", "expect_high": False},
        ],
    },
    "Cancer": {
        "element": "water",
        "profile": "Nurturing, emotional, and fiercely protective of loved ones. Cancer feels everything deeply and creates emotional safe havens for those they love.",
        "traits": [
            {"trait": "Empathy", "prediction": "Cancer is the most nurturing sign — deeply empathic and attuned to others' emotions.", "dimension": "empathy", "expect_high": True},
            {"trait": "Emotional Sensitivity", "prediction": "Water signs feel everything — Cancer absorbs emotional energy from their environment.", "dimension": "emotional_sensitivity", "expect_high": True},
            {"trait": "Stress Reactivity", "prediction": "Cancer's deep emotional nature means stress hits hard and recovery takes time.", "dimension": "stress_reactivity", "expect_high": True},
            {"trait": "Anxiety Sensitivity", "prediction": "The protective instinct creates a tendency toward worry and anxiety.", "dimension": "anxiety", "expect_high": True},
        ],
    },
    "Leo": {
        "element": "fire",
        "profile": "Confident, dramatic, and magnetically generous. Leo is the performer of the zodiac — radiating warmth, creativity, and an unshakeable belief in themselves.",
        "traits": [
            {"trait": "Stress Resilience", "prediction": "Leo's fire burns bright under pressure — natural confidence creates stress resilience.", "dimension": "stress_resilience", "expect_high": True},
            {"trait": "Risk Tolerance", "prediction": "Leo takes bold risks — the spotlight rewards those who dare.", "dimension": "risk_tolerance", "expect_high": True},
            {"trait": "Brain Plasticity", "prediction": "Creative and adaptable — Leo's brain thrives on new creative challenges.", "dimension": "neuroplasticity", "expect_high": True},
            {"trait": "Empathy", "prediction": "Leo is generous but self-focused — empathy can be overshadowed by ego.", "dimension": "empathy", "expect_high": False},
        ],
    },
    "Virgo": {
        "element": "earth",
        "profile": "Analytical, perfectionist, and quietly anxious. Virgo sees every detail others miss — a natural problem-solver who holds themselves to impossibly high standards.",
        "traits": [
            {"trait": "Cognitive Performance", "prediction": "Mercury-ruled like Gemini, but Virgo channels it into meticulous analysis.", "dimension": "cognitive_performance", "expect_high": True},
            {"trait": "Anxiety Sensitivity", "prediction": "Virgo's perfectionism feeds a constant undercurrent of anxiety and self-criticism.", "dimension": "anxiety", "expect_high": True},
            {"trait": "Behavioral Control", "prediction": "Disciplined and methodical — Virgo controls impulses through analysis.", "dimension": "impulsivity", "expect_high": False},
            {"trait": "Stress Reactivity", "prediction": "Internalizes stress through overthinking — body holds the tension the mind creates.", "dimension": "stress_reactivity", "expect_high": True},
        ],
    },
    "Libra": {
        "element": "air",
        "profile": "Harmonious, social, and endlessly diplomatic. Libra seeks balance in all things — a natural mediator who craves beauty, partnership, and justice.",
        "traits": [
            {"trait": "Empathy", "prediction": "Libra is the diplomat — highly attuned to others' needs and emotions.", "dimension": "empathy", "expect_high": True},
            {"trait": "Anxiety Sensitivity", "prediction": "The desire to please everyone creates decision paralysis and social anxiety.", "dimension": "anxiety", "expect_high": True},
            {"trait": "Emotional Sensitivity", "prediction": "Libra absorbs disharmony — conflict is physically uncomfortable.", "dimension": "emotional_sensitivity", "expect_high": True},
            {"trait": "Risk Tolerance", "prediction": "Prefers safe, balanced choices — avoids rocking the boat.", "dimension": "risk_tolerance", "expect_high": False},
        ],
    },
    "Scorpio": {
        "element": "water",
        "profile": "Intense, perceptive, and powerfully magnetic. Scorpio sees through facades and experiences life at maximum emotional depth — nothing is surface-level.",
        "traits": [
            {"trait": "Emotional Sensitivity", "prediction": "Scorpio feels everything at maximum intensity — the deepest emotional sign.", "dimension": "emotional_sensitivity", "expect_high": True},
            {"trait": "Spiritual Sensitivity", "prediction": "Drawn to the hidden and mystical — Scorpio has a natural affinity for the unseen.", "dimension": "spiritual_sensitivity", "expect_high": True},
            {"trait": "Stress Resilience", "prediction": "Scorpio transforms through crisis — uses stress as a catalyst for rebirth.", "dimension": "stress_resilience", "expect_high": True},
            {"trait": "Behavioral Control", "prediction": "Scorpio is strategic and controlled — emotions run deep but are rarely displayed.", "dimension": "impulsivity", "expect_high": False},
        ],
    },
    "Sagittarius": {
        "element": "fire",
        "profile": "Adventurous, philosophical, and refreshingly honest. Sagittarius is the explorer of the zodiac — driven by an insatiable hunger for meaning, freedom, and truth.",
        "traits": [
            {"trait": "Risk Tolerance", "prediction": "The eternal adventurer — Sagittarius leaps before looking.", "dimension": "risk_tolerance", "expect_high": True},
            {"trait": "Novelty Seeking", "prediction": "Driven by wanderlust and intellectual curiosity — always seeking the next horizon.", "dimension": "novelty_seeking", "expect_high": True},
            {"trait": "Self-Transcendence", "prediction": "The philosopher of the zodiac — drawn to questions of meaning and the divine.", "dimension": "transcendence", "expect_high": True},
            {"trait": "Behavioral Control", "prediction": "Impulsive and brutally honest — tact is not Sagittarius' strong suit.", "dimension": "impulsivity", "expect_high": True},
        ],
    },
    "Capricorn": {
        "element": "earth",
        "profile": "Ambitious, disciplined, and strategically patient. Capricorn is the master builder — playing the long game with unmatched determination and practical wisdom.",
        "traits": [
            {"trait": "Behavioral Control", "prediction": "Capricorn is the most disciplined sign — iron willpower and emotional control.", "dimension": "impulsivity", "expect_high": False},
            {"trait": "Stress Resilience", "prediction": "Saturn-ruled resilience — Capricorn treats obstacles as stepping stones.", "dimension": "stress_resilience", "expect_high": True},
            {"trait": "Risk Tolerance", "prediction": "Calculated risks only — Capricorn never gambles without a backup plan.", "dimension": "risk_tolerance", "expect_high": False},
            {"trait": "Anxiety Sensitivity", "prediction": "The weight of responsibility creates an undercurrent of worry.", "dimension": "anxiety", "expect_high": True},
        ],
    },
    "Aquarius": {
        "element": "air",
        "profile": "Independent, innovative, and delightfully eccentric. Aquarius marches to their own drumbeat — a visionary who challenges convention and champions the collective.",
        "traits": [
            {"trait": "Novelty Seeking", "prediction": "The innovator — Aquarius actively seeks unconventional ideas and experiences.", "dimension": "novelty_seeking", "expect_high": True},
            {"trait": "Cognitive Performance", "prediction": "Brilliant and analytical — Aquarius excels at abstract, systems-level thinking.", "dimension": "cognitive_performance", "expect_high": True},
            {"trait": "Empathy", "prediction": "Aquarius cares about humanity but can be emotionally detached from individuals.", "dimension": "empathy", "expect_high": False},
            {"trait": "Behavioral Control", "prediction": "Rebellious and non-conformist — conventional rules don't apply.", "dimension": "impulsivity", "expect_high": True},
        ],
    },
    "Pisces": {
        "element": "water",
        "profile": "Intuitive, empathic, and boundlessly creative. Pisces swims between reality and dreams — a natural mystic who feels the collective unconscious as vividly as the waking world.",
        "traits": [
            {"trait": "Empathy", "prediction": "The empath of the zodiac — Pisces absorbs others' emotions like a sponge.", "dimension": "empathy", "expect_high": True},
            {"trait": "Spiritual Sensitivity", "prediction": "Most spiritually attuned sign — Pisces has a natural connection to the transcendent.", "dimension": "spiritual_sensitivity", "expect_high": True},
            {"trait": "Self-Transcendence", "prediction": "Drawn to dissolving ego boundaries — Pisces seeks union with something greater.", "dimension": "transcendence", "expect_high": True},
            {"trait": "Emotional Sensitivity", "prediction": "Deeply sensitive — Pisces feels everything and struggles with emotional boundaries.", "dimension": "emotional_sensitivity", "expect_high": True},
        ],
    },
}

# Moon sign profiles map to emotional patterns by element
MOON_ELEMENT_PROFILES: Dict[str, Dict[str, Any]] = {
    "fire": {
        "profile": "Quick emotional reactions, optimistic recovery, passionate intensity. Fire Moons process emotions through action — they need to DO something about how they feel.",
        "traits": [
            {"trait": "Stress Resilience", "prediction": "Fire Moons bounce back quickly — emotions burn hot but don't linger.", "dimension": "stress_resilience", "expect_high": True},
            {"trait": "Behavioral Control", "prediction": "Emotions are expressed immediately and explosively — no poker face here.", "dimension": "impulsivity", "expect_high": True},
        ],
    },
    "earth": {
        "profile": "Steady emotional processing, practical coping, slow to anger. Earth Moons are emotionally grounded — they feel deeply but express through practical action.",
        "traits": [
            {"trait": "Stress Reactivity", "prediction": "Earth Moons absorb stress slowly — emotionally stable and hard to rattle.", "dimension": "stress_reactivity", "expect_high": False},
            {"trait": "Behavioral Control", "prediction": "Measured emotional responses — Earth Moons think before they react.", "dimension": "impulsivity", "expect_high": False},
        ],
    },
    "air": {
        "profile": "Intellectual emotional processing, analytical detachment, social regulation. Air Moons think about their feelings rather than just feeling them.",
        "traits": [
            {"trait": "Cognitive Performance", "prediction": "Air Moons process emotions cerebrally — high emotional intelligence through analysis.", "dimension": "cognitive_performance", "expect_high": True},
            {"trait": "Emotional Sensitivity", "prediction": "Emotionally present but intellectually buffered — Air Moons maintain perspective.", "dimension": "emotional_sensitivity", "expect_high": False},
        ],
    },
    "water": {
        "profile": "Deep emotional intensity, empathic absorption, psychic sensitivity. Water Moons feel the emotional undercurrents that others miss entirely.",
        "traits": [
            {"trait": "Empathy", "prediction": "Water Moons are emotional sponges — they absorb the feelings of everyone around them.", "dimension": "empathy", "expect_high": True},
            {"trait": "Emotional Sensitivity", "prediction": "Extremely sensitive inner emotional world — Water Moons feel everything deeply.", "dimension": "emotional_sensitivity", "expect_high": True},
        ],
    },
}

# Rising sign profiles by element
RISING_ELEMENT_PROFILES: Dict[str, Dict[str, Any]] = {
    "fire": {
        "profile": "Energetic, direct, and confidently assertive first impression. Fire Rising people light up a room with their natural charisma and action-oriented presence.",
        "traits": [
            {"trait": "Risk Tolerance", "prediction": "Fire Rising projects boldness — outwardly fearless and action-oriented.", "dimension": "risk_tolerance", "expect_high": True},
            {"trait": "Stress Resilience", "prediction": "Projects composure under pressure — Fire Rising's external armor.", "dimension": "stress_resilience", "expect_high": True},
        ],
    },
    "earth": {
        "profile": "Calm, reliable, and grounded presence. Earth Rising people project stability and competence — others instinctively trust their steadiness.",
        "traits": [
            {"trait": "Behavioral Control", "prediction": "Earth Rising projects measured composure — outwardly disciplined.", "dimension": "impulsivity", "expect_high": False},
            {"trait": "Brain Plasticity", "prediction": "Projects a practical, learning-oriented approach to life.", "dimension": "neuroplasticity", "expect_high": True},
        ],
    },
    "air": {
        "profile": "Friendly, curious, and intellectually engaging presence. Air Rising people draw others in with their wit, charm, and natural communication skills.",
        "traits": [
            {"trait": "Cognitive Performance", "prediction": "Air Rising projects sharp intellect — first impression is bright and articulate.", "dimension": "cognitive_performance", "expect_high": True},
            {"trait": "Memory", "prediction": "Projects excellent recall and conversational depth.", "dimension": "memory", "expect_high": True},
        ],
    },
    "water": {
        "profile": "Sensitive, mysterious, and deeply intuitive presence. Water Rising people seem to see right through you — projecting emotional depth and psychic sensitivity.",
        "traits": [
            {"trait": "Emotional Sensitivity", "prediction": "Water Rising radiates emotional awareness — others sense your depth.", "dimension": "emotional_sensitivity", "expect_high": True},
            {"trait": "Spiritual Sensitivity", "prediction": "Projects mystical, otherworldly quality — drawn to the hidden.", "dimension": "spiritual_sensitivity", "expect_high": True},
        ],
    },
}

# Planetary insights (non-scored, narrative only)
PLANET_NARRATIVES: Dict[str, Dict[str, str]] = {
    "mercury": {
        "Aries": "Your communication style is direct, fast, and unfiltered. Mercury in Aries speaks first and edits later — a mind built for quick decisions, not careful diplomacy.",
        "Taurus": "Your thinking is methodical and sensual. Mercury in Taurus processes information slowly but deeply — ideas are tested against practical reality before being accepted.",
        "Gemini": "Mercury is at home in Gemini — your mind is a high-speed information processor. Quick-witted, verbally gifted, and endlessly curious. The risk: intellectual scatter.",
        "Cancer": "Your thinking is colored by emotion and memory. Mercury in Cancer processes information through feelings — intuition often outperforms logic.",
        "Leo": "Your communication style is dramatic, warm, and persuasive. Mercury in Leo speaks to inspire — every conversation has a storytelling quality.",
        "Virgo": "Mercury is exalted in Virgo — your mind is a precision instrument. Analytical, detail-oriented, and methodical. The risk: paralysis by analysis.",
        "Libra": "Your thinking seeks balance and fairness. Mercury in Libra weighs every perspective — excellent at diplomacy but sometimes paralyzed by seeing all sides.",
        "Scorpio": "Your mind penetrates surfaces to find hidden truths. Mercury in Scorpio is the detective of the zodiac — nothing escapes your probing intellect.",
        "Sagittarius": "Your thinking is philosophical and expansive. Mercury in Sagittarius connects ideas across cultures and disciplines — the big picture always matters more than details.",
        "Capricorn": "Your communication is strategic and purposeful. Mercury in Capricorn thinks in terms of goals and outcomes — every word serves a plan.",
        "Aquarius": "Your mind operates in unconventional patterns. Mercury in Aquarius thinks in systems and innovations — you see solutions others can't imagine.",
        "Pisces": "Your thinking is intuitive and imaginative. Mercury in Pisces processes information through imagery and feeling — logic takes a back seat to creative intuition.",
    },
    "venus": {
        "Aries": "You love with passionate intensity and directness. Venus in Aries falls fast, pursues boldly, and values excitement in relationships over comfort.",
        "Taurus": "Venus is at home in Taurus — you love deeply, sensually, and loyally. Beauty, comfort, and physical affection are essential to your happiness.",
        "Gemini": "You love through communication and intellectual connection. Venus in Gemini needs mental stimulation in relationships — flirtation is an art form.",
        "Cancer": "You love through nurturing and emotional security. Venus in Cancer creates deep domestic bonds — home and family are the center of your heart.",
        "Leo": "You love with dramatic generosity and warmth. Venus in Leo wants to be adored — and gives lavishly in return. Romance should be a grand production.",
        "Virgo": "You love through acts of service and attention to detail. Venus in Virgo shows care through practical gestures — devotion expressed through doing.",
        "Libra": "Venus is at home in Libra — you love through harmony, beauty, and partnership. Relationships are the central organizing principle of your life.",
        "Scorpio": "You love with transformative intensity. Venus in Scorpio bonds at the deepest level — casual connections hold no interest. All or nothing.",
        "Sagittarius": "You love through shared adventure and philosophical connection. Venus in Sagittarius needs freedom within partnership — love should expand your world.",
        "Capricorn": "You love through commitment and building together. Venus in Capricorn values long-term partnership — love is an investment that compounds over time.",
        "Aquarius": "You love through intellectual connection and shared ideals. Venus in Aquarius values friendship within romance — conventional relationship rules don't apply.",
        "Pisces": "Venus is exalted in Pisces — you love with boundless compassion and romantic idealism. Soulmate connections and transcendent love are your deepest longing.",
    },
    "mars": {
        "Aries": "Mars is at home in Aries — your drive is explosive, competitive, and direct. You attack challenges head-on with relentless energy.",
        "Taurus": "Your drive is steady and unstoppable. Mars in Taurus moves slowly but with immense force — once you start, nothing can stop you.",
        "Gemini": "Your energy is mental and versatile. Mars in Gemini fights with words and wit — you win through intelligence, not brute force.",
        "Cancer": "Your drive is protective and emotionally motivated. Mars in Cancer fights for family and emotional security — the mama/papa bear energy.",
        "Leo": "Your drive is creative and attention-seeking. Mars in Leo channels energy through performance and creative expression.",
        "Virgo": "Your drive is precise and methodical. Mars in Virgo channels energy through efficient systems and detailed planning.",
        "Libra": "Your drive seeks justice and balance. Mars in Libra fights for fairness — but can be passive-aggressive when direct conflict feels too disharmonious.",
        "Scorpio": "Mars is powerful in Scorpio — your drive is strategic, intense, and transformative. You play the long game with laser focus.",
        "Sagittarius": "Your drive is expansive and philosophical. Mars in Sagittarius channels energy through exploration, adventure, and the pursuit of truth.",
        "Capricorn": "Mars is exalted in Capricorn — your drive is disciplined, strategic, and goal-oriented. You climb mountains methodically and relentlessly.",
        "Aquarius": "Your drive is innovative and collective. Mars in Aquarius channels energy toward humanitarian goals and systemic change.",
        "Pisces": "Your drive is intuitive and compassionate. Mars in Pisces channels energy through creative and spiritual pursuits — action guided by inspiration.",
    },
    "jupiter": {
        "Aries": "Jupiter in Aries expands through bold action and pioneering initiatives. Growth comes through courage and taking the first step.",
        "Taurus": "Jupiter in Taurus expands through material abundance and sensory pleasure. Growth comes through building tangible value.",
        "Gemini": "Jupiter in Gemini expands through learning and communication. Growth comes through intellectual exploration and diverse connections.",
        "Cancer": "Jupiter in Cancer is exalted — emotional growth and family connections bring your greatest blessings.",
        "Leo": "Jupiter in Leo expands through creative expression and generous leadership. Growth comes through sharing your gifts.",
        "Virgo": "Jupiter in Virgo expands through service and mastery of craft. Growth comes through meticulous skill development.",
        "Libra": "Jupiter in Libra expands through partnership and social harmony. Growth comes through collaborative relationships.",
        "Scorpio": "Jupiter in Scorpio expands through emotional depth and transformation. Growth comes through confronting hidden truths.",
        "Sagittarius": "Jupiter is at home in Sagittarius — expansion through travel, philosophy, and higher learning is your natural path.",
        "Capricorn": "Jupiter in Capricorn expands through disciplined achievement. Growth comes through patient, strategic long-term building.",
        "Aquarius": "Jupiter in Aquarius expands through innovation and collective vision. Growth comes through embracing the unconventional.",
        "Pisces": "Jupiter is at home in Pisces — spiritual expansion and compassionate connection bring your greatest blessings.",
    },
    "saturn": {
        "Aries": "Saturn in Aries teaches patience and strategic action. Your life lesson: channel impulsive energy into sustained, purposeful direction.",
        "Taurus": "Saturn in Taurus teaches the difference between security and hoarding. Your life lesson: build real value without clinging to material comfort.",
        "Gemini": "Saturn in Gemini teaches depth of thought. Your life lesson: go deep instead of wide — master one thing instead of sampling everything.",
        "Cancer": "Saturn in Cancer teaches emotional boundaries. Your life lesson: nurture others without losing yourself in their needs.",
        "Leo": "Saturn in Leo teaches authentic confidence. Your life lesson: lead through substance, not performance — real authority earned, not demanded.",
        "Virgo": "Saturn in Virgo teaches self-compassion. Your life lesson: excellence without perfectionism — good enough is sometimes the highest wisdom.",
        "Libra": "Saturn is exalted in Libra — teaching the structure of fair relationships. Your life lesson: justice requires hard choices, not just diplomacy.",
        "Scorpio": "Saturn in Scorpio teaches the discipline of transformation. Your life lesson: control without manipulation — power wielded with integrity.",
        "Sagittarius": "Saturn in Sagittarius teaches grounded wisdom. Your life lesson: freedom within structure — philosophy tested against reality.",
        "Capricorn": "Saturn is at home in Capricorn — teaching mastery through patience. Your life lesson: the mountain is climbed one step at a time.",
        "Aquarius": "Saturn is at home in Aquarius — teaching innovation within structure. Your life lesson: revolution requires organization.",
        "Pisces": "Saturn in Pisces teaches spiritual discipline. Your life lesson: transcendence through practice, not escape — grounded mysticism.",
    },
}


# ─────────────────────────────────────────────────────────────────────────────
# NATAL CHART CALCULATION
# ─────────────────────────────────────────────────────────────────────────────

def _longitude_to_sign(longitude: float) -> Tuple[str, float]:
    """Convert ecliptic longitude to zodiac sign and degree within sign."""
    sign_index = int(longitude / 30) % 12
    degree = longitude % 30
    return SIGNS[sign_index], round(degree, 2)


def _datetime_to_jd(date_str: str, time_str: Optional[str], tz_offset_hours: float = 0.0) -> float:
    """Convert date/time to Julian Day for Swiss Ephemeris.

    Args:
        date_str: ISO date 'YYYY-MM-DD'
        time_str: Time 'HH:MM' or None (defaults to 12:00)
        tz_offset_hours: Timezone offset from UTC in hours
    """
    parts = date_str.split("-")
    year = int(parts[0])
    month = int(parts[1])
    day = int(parts[2])

    if time_str:
        t_parts = time_str.split(":")
        hour = int(t_parts[0])
        minute = int(t_parts[1]) if len(t_parts) > 1 else 0
    else:
        hour = 12
        minute = 0

    # Convert to UT (Universal Time)
    ut_hour = hour + minute / 60.0 - tz_offset_hours

    if HAS_SWE:
        jd = swe.julday(year, month, day, ut_hour)
    else:
        # Fallback: approximate Julian Day calculation
        # Meeus algorithm
        if month <= 2:
            year -= 1
            month += 12
        a = year // 100
        b = 2 - a + a // 4
        jd = int(365.25 * (year + 4716)) + int(30.6001 * (month + 1)) + day + ut_hour / 24.0 + b - 1524.5
    return jd


def _get_tz_offset(tz_name: str, date_str: str) -> float:
    """Get UTC offset in hours for a timezone name.

    Simple lookup for common timezones. Falls back to 0 (UTC) for unknown.
    """
    # Common timezone offsets (standard time — doesn't account for DST)
    # For astrology, +-1 hour DST error shifts the chart by ~15 degrees for fast planets
    # but Sun/Moon are stable enough that this rarely changes the sign
    TZ_OFFSETS = {
        "America/Sao_Paulo": -3, "America/New_York": -5, "America/Chicago": -6,
        "America/Denver": -7, "America/Los_Angeles": -8, "America/Anchorage": -9,
        "America/Honolulu": -10, "America/Buenos_Aires": -3, "America/Santiago": -4,
        "America/Bogota": -5, "America/Mexico_City": -6, "America/Lima": -5,
        "America/Caracas": -4, "America/Toronto": -5, "America/Vancouver": -8,
        "America/Manaus": -4, "America/Recife": -3, "America/Fortaleza": -3,
        "America/Belem": -3, "America/Cuiaba": -4, "America/Porto_Velho": -4,
        "America/Campo_Grande": -4, "America/Rio_Branco": -5, "America/Maceio": -3,
        "Europe/London": 0, "Europe/Paris": 1, "Europe/Berlin": 1,
        "Europe/Rome": 1, "Europe/Madrid": 1, "Europe/Lisbon": 0,
        "Europe/Amsterdam": 1, "Europe/Brussels": 1, "Europe/Vienna": 1,
        "Europe/Zurich": 1, "Europe/Warsaw": 1, "Europe/Prague": 1,
        "Europe/Budapest": 1, "Europe/Bucharest": 2, "Europe/Helsinki": 2,
        "Europe/Athens": 2, "Europe/Moscow": 3, "Europe/Istanbul": 3,
        "Europe/Kiev": 2, "Europe/Stockholm": 1, "Europe/Oslo": 1,
        "Europe/Copenhagen": 1, "Europe/Dublin": 0,
        "Asia/Tokyo": 9, "Asia/Shanghai": 8, "Asia/Hong_Kong": 8,
        "Asia/Seoul": 9, "Asia/Kolkata": 5.5, "Asia/Mumbai": 5.5,
        "Asia/Dubai": 4, "Asia/Riyadh": 3, "Asia/Singapore": 8,
        "Asia/Bangkok": 7, "Asia/Jakarta": 7, "Asia/Karachi": 5,
        "Asia/Tehran": 3.5, "Asia/Baghdad": 3, "Asia/Beirut": 2,
        "Asia/Jerusalem": 2, "Asia/Taipei": 8, "Asia/Manila": 8,
        "Africa/Cairo": 2, "Africa/Lagos": 1, "Africa/Johannesburg": 2,
        "Africa/Nairobi": 3, "Africa/Casablanca": 1,
        "Australia/Sydney": 10, "Australia/Melbourne": 10,
        "Australia/Perth": 8, "Australia/Brisbane": 10,
        "Pacific/Auckland": 12, "Pacific/Fiji": 12,
        "UTC": 0, "GMT": 0,
    }
    return TZ_OFFSETS.get(tz_name, 0)


def calculate_natal_chart(
    date_str: str,
    time_str: Optional[str],
    lat: float,
    lng: float,
    tz_name: str,
) -> Dict[str, Any]:
    """Calculate natal chart placements using Swiss Ephemeris."""

    tz_offset = _get_tz_offset(tz_name, date_str)
    jd = _datetime_to_jd(date_str, time_str, tz_offset)

    chart: Dict[str, Any] = {}

    if HAS_SWE:
        swe.set_ephe_path(None)  # Use built-in ephemeris

        for planet_name, planet_id in PLANET_IDS.items():
            result = swe.calc_ut(jd, planet_id)
            longitude = result[0][0]
            sign, degree = _longitude_to_sign(longitude)
            chart[planet_name] = {"sign": sign, "degree": degree}

        # Calculate houses (Placidus system)
        if time_str:
            houses_data = swe.houses(jd, lat, lng, b'P')
            cusps = houses_data[0]
            asc_longitude = houses_data[1][0]
            asc_sign, asc_degree = _longitude_to_sign(asc_longitude)
            chart["rising"] = {"sign": asc_sign, "degree": asc_degree}

            # Assign houses to planets
            for planet_name in PLANET_IDS:
                p_long = chart[planet_name]["degree"] + SIGNS.index(chart[planet_name]["sign"]) * 30
                for h in range(12):
                    next_h = (h + 1) % 12
                    start = cusps[h]
                    end = cusps[next_h]
                    if start > end:
                        if p_long >= start or p_long < end:
                            chart[planet_name]["house"] = h + 1
                            break
                    else:
                        if start <= p_long < end:
                            chart[planet_name]["house"] = h + 1
                            break
        else:
            # No time — Rising is approximate (use noon Ascendant)
            houses_data = swe.houses(jd, lat, lng, b'P')
            asc_longitude = houses_data[1][0]
            asc_sign, asc_degree = _longitude_to_sign(asc_longitude)
            chart["rising"] = {"sign": asc_sign, "degree": asc_degree}
    else:
        # Fallback without pyswisseph — use approximate formula
        # This gives reasonable Sun sign accuracy but Moon/Rising will be off
        chart = _fallback_chart(date_str, time_str, lat, lng, tz_offset)

    return chart


def _fallback_chart(
    date_str: str,
    time_str: Optional[str],
    lat: float,
    lng: float,
    tz_offset: float,
) -> Dict[str, Any]:
    """Approximate natal chart without pyswisseph.
    Uses basic astronomical formulas. Sun sign is accurate;
    Moon and other planets are approximate.
    """
    parts = date_str.split("-")
    year, month, day = int(parts[0]), int(parts[1]), int(parts[2])
    hour = 12
    if time_str:
        t = time_str.split(":")
        hour = int(t[0]) + (int(t[1]) if len(t) > 1 else 0) / 60.0

    # Approximate Sun longitude (accurate to ~1 degree)
    # Days since J2000.0 (Jan 1.5, 2000)
    jd = _datetime_to_jd(date_str, time_str, tz_offset)
    n = jd - 2451545.0
    L = (280.460 + 0.9856474 * n) % 360
    g = math.radians((357.528 + 0.9856003 * n) % 360)
    sun_long = (L + 1.915 * math.sin(g) + 0.020 * math.sin(2 * g)) % 360

    # Approximate Moon longitude (accurate to ~5 degrees)
    moon_long = (218.32 + 13.1764 * n) % 360

    # Other planets — very rough approximations
    merc_long = (sun_long + 20 * math.sin(math.radians(n * 4.09))) % 360
    venus_long = (sun_long + 46 * math.sin(math.radians(n * 1.62))) % 360
    mars_long = (355.45 + 0.5240 * n) % 360
    jup_long = (34.40 + 0.0831 * n) % 360
    sat_long = (50.08 + 0.0335 * n) % 360

    # Rising sign — very rough (based on local sidereal time)
    lst = (100.46 + 0.985647 * n + lng + (hour - tz_offset) * 15) % 360
    asc_long = lst  # Very approximate

    chart = {}
    for name, lon in [("sun", sun_long), ("moon", moon_long), ("mercury", merc_long),
                       ("venus", venus_long), ("mars", mars_long), ("jupiter", jup_long),
                       ("saturn", sat_long), ("rising", asc_long)]:
        sign, degree = _longitude_to_sign(lon)
        chart[name] = {"sign": sign, "degree": degree}

    return chart


# ─────────────────────────────────────────────────────────────────────────────
# DNA ANALYSIS
# ─────────────────────────────────────────────────────────────────────────────

_COMPLEMENT = str.maketrans("ACGT", "TGCA")


def _read_genetic_markers(
    variants: Dict[str, Tuple[str, str, str]],
) -> Dict[str, Dict[str, Any]]:
    """Read all personality-relevant genetic markers from the user's genome."""
    results: Dict[str, Dict[str, Any]] = {}

    for rsid, marker_info in GENETIC_MARKERS.items():
        rsid_lower = rsid.lower()
        if rsid_lower in variants:
            chrom, pos, gt = variants[rsid_lower]
            gt_clean = gt.strip().upper().replace("-", "")
            if not gt_clean:
                continue

            # Sort alleles for consistent matching
            gt_sorted = "".join(sorted(gt_clean))
            # Also try complement strand (A↔T, C↔G) for providers that report opposite strand
            gt_comp = gt_clean.translate(_COMPLEMENT)
            gt_comp_sorted = "".join(sorted(gt_comp))

            # Find matching interpretation (try direct first, then complement)
            for genotype, interp in marker_info["interpretations"].items():
                gs = "".join(sorted(genotype))
                if gs == gt_sorted or gs == gt_comp_sorted:
                    results[rsid] = {
                        "gene": marker_info["gene"],
                        "dimension": marker_info["dimension"],
                        "label": marker_info["label"],
                        "genotype": gt_clean,
                        "interpretation": interp["label"],
                        "score": interp["score"],
                        "description": interp["desc"],
                    }
                    break

    return results


# ─────────────────────────────────────────────────────────────────────────────
# COMPARISON ENGINE
# ─────────────────────────────────────────────────────────────────────────────

def _compare_trait(
    trait_def: Dict[str, Any],
    genetic_results: Dict[str, Dict[str, Any]],
) -> Optional[Dict[str, Any]]:
    """Compare one astrological trait prediction against genetic data."""
    dimension = trait_def["dimension"]
    expect_high = trait_def["expect_high"]

    # Find the genetic marker for this dimension
    for rsid, result in genetic_results.items():
        if result["dimension"] == dimension:
            score = result["score"]

            # Determine concordance
            if expect_high:
                if score >= 0.7:
                    verdict = "aligned"
                elif score >= 0.4:
                    verdict = "partially_aligned"
                elif score <= 0.1:
                    verdict = "contradicted"
                else:
                    verdict = "inconclusive"
            else:
                if score <= 0.3:
                    verdict = "aligned"
                elif score <= 0.6:
                    verdict = "partially_aligned"
                elif score >= 0.9:
                    verdict = "contradicted"
                else:
                    verdict = "inconclusive"

            marker_info = GENETIC_MARKERS.get(rsid, {})

            return {
                "trait": trait_def["trait"],
                "astrologyPrediction": trait_def["prediction"],
                "geneticFinding": result["description"],
                "gene": result["gene"],
                "rsid": rsid,
                "genotype": result["genotype"],
                "verdict": verdict,
                "confidenceNote": f"{result['gene']} {rsid} — {result['interpretation']}",
            }

    return None


def _build_sign_comparison(
    sign: str,
    profile_data: Dict[str, Any],
    genetic_results: Dict[str, Dict[str, Any]],
) -> Dict[str, Any]:
    """Build comparison for a zodiac sign (sun, moon element, or rising element)."""
    trait_comparisons = []
    aligned_count = 0
    total_count = 0

    for trait_def in profile_data.get("traits", []):
        comparison = _compare_trait(trait_def, genetic_results)
        if comparison:
            trait_comparisons.append(comparison)
            total_count += 1
            if comparison["verdict"] == "aligned":
                aligned_count += 1
            elif comparison["verdict"] == "partially_aligned":
                aligned_count += 0.5

    alignment_score = round(aligned_count / total_count * 100, 1) if total_count > 0 else 50.0

    if alignment_score >= 70:
        verdict = "aligned"
        narrative = f"Your DNA strongly confirms the {sign} archetype. The genetic markers for {sign}'s key traits show remarkable alignment with astrological predictions."
    elif alignment_score >= 45:
        verdict = "partially_aligned"
        narrative = f"Your DNA partially aligns with the {sign} archetype. Some traits match the astrological predictions while others tell a different genetic story."
    elif alignment_score < 30:
        verdict = "contradicted"
        narrative = f"Your DNA tells a different story than {sign} would predict. Your genetic personality markers diverge significantly from the {sign} archetype."
    else:
        verdict = "inconclusive"
        narrative = f"The genetic evidence for {sign}'s traits is mixed — some markers align, others are neutral. Your personality is more nuanced than any single archetype."

    return {
        "sign": sign,
        "element": SIGN_ELEMENTS.get(sign, "unknown"),
        "astrologyProfile": profile_data.get("profile", ""),
        "traitComparisons": trait_comparisons,
        "overallVerdict": verdict,
        "alignmentScore": alignment_score,
        "narrative": narrative,
    }


# ─────────────────────────────────────────────────────────────────────────────
# MAIN ANALYSIS
# ─────────────────────────────────────────────────────────────────────────────

def analyze_stars_genes(
    variants: Dict[str, Tuple[str, str, str]],
    birth_data: Dict[str, Any],
) -> Dict[str, Any]:
    """
    Main analysis function: calculate natal chart + compare against DNA.

    Args:
        variants: rsID -> (chrom, pos, genotype)
        birth_data: {"birthDate": "1990-03-15", "birthTime": "14:30" or null,
                     "birthTimeKnown": bool, "birthLat": float, "birthLng": float,
                     "birthTz": "America/Sao_Paulo", "birthCity": "São Paulo, Brazil"}

    Returns:
        Dict with full report data
    """
    date_str = birth_data["birthDate"]
    time_str = birth_data.get("birthTime")
    birth_time_known = birth_data.get("birthTimeKnown", bool(time_str))
    lat = birth_data["birthLat"]
    lng = birth_data["birthLng"]
    tz_name = birth_data["birthTz"]
    city = birth_data.get("birthCity", "Unknown")

    # 1. Calculate natal chart
    chart = calculate_natal_chart(date_str, time_str, lat, lng, tz_name)

    # 2. Read genetic markers
    genetic_results = _read_genetic_markers(variants)

    # 3. Build Sun Sign comparison
    sun_sign = chart["sun"]["sign"]
    sun_profile = SUN_SIGN_PROFILES.get(sun_sign, {})
    sun_comparison = _build_sign_comparison(sun_sign, sun_profile, genetic_results)

    # 4. Build Moon Sign comparison (by element)
    moon_sign = chart["moon"]["sign"]
    moon_element = SIGN_ELEMENTS.get(moon_sign, "water")
    moon_profile = MOON_ELEMENT_PROFILES.get(moon_element, {})
    moon_comparison = _build_sign_comparison(moon_sign, moon_profile, genetic_results)

    # 5. Build Rising Sign comparison (by element)
    rising_sign = chart["rising"]["sign"]
    rising_element = SIGN_ELEMENTS.get(rising_sign, "fire")
    rising_profile = RISING_ELEMENT_PROFILES.get(rising_element, {})
    rising_comparison = _build_sign_comparison(rising_sign, rising_profile, genetic_results)

    # 6. Build planetary insights (narrative only, no scoring)
    planetary_insights = []
    for planet in ["mercury", "venus", "mars", "jupiter", "saturn"]:
        if planet in chart:
            p_sign = chart[planet]["sign"]
            narratives = PLANET_NARRATIVES.get(planet, {})
            narrative = narratives.get(p_sign, f"Your {planet.title()} is in {p_sign}.")

            # Find genetic parallel if one exists
            genetic_parallel = None
            if planet == "mercury":
                for rsid, r in genetic_results.items():
                    if r["dimension"] == "cognitive_performance":
                        genetic_parallel = f"Your {r['gene']} {rsid} shows {r['interpretation'].lower()} — {'supporting' if r['score'] >= 0.5 else 'contrasting with'} Mercury in {p_sign}'s intellectual style."
                        break
            elif planet == "venus":
                for rsid, r in genetic_results.items():
                    if r["dimension"] == "empathy":
                        genetic_parallel = f"Your {r['gene']} {rsid} shows {r['interpretation'].lower()} — {'resonating with' if r['score'] >= 0.5 else 'adding nuance to'} Venus in {p_sign}'s love language."
                        break
            elif planet == "mars":
                for rsid, r in genetic_results.items():
                    if r["dimension"] == "risk_tolerance":
                        genetic_parallel = f"Your {r['gene']} {rsid} shows {r['interpretation'].lower()} — {'fueling' if r['score'] >= 0.5 else 'tempering'} Mars in {p_sign}'s drive."
                        break

            planetary_insights.append({
                "planet": planet.title(),
                "sign": p_sign,
                "meaning": narrative,
                "geneticParallel": genetic_parallel,
            })

    # 7. Calculate Cosmic Alignment Score
    sun_score = sun_comparison["alignmentScore"]
    moon_score = moon_comparison["alignmentScore"]
    rising_score = rising_comparison["alignmentScore"]
    overall_score = round((sun_score * 0.5 + moon_score * 0.3 + rising_score * 0.2), 1)

    if overall_score >= 70:
        score_interp = (
            f"Your DNA and your stars are remarkably aligned! With a {overall_score}% Cosmic Alignment Score, "
            f"your genetic personality profile closely matches what astrology predicts for a {sun_sign} Sun, "
            f"{moon_sign} Moon, {rising_sign} Rising individual. Whether the stars knew something all along "
            f"or it's a beautiful coincidence — your cosmic and genetic identities resonate deeply."
        )
    elif overall_score >= 45:
        score_interp = (
            f"Your DNA and stars show partial alignment at {overall_score}%. Some astrological predictions "
            f"are confirmed by your genetics while others diverge — suggesting your personality is shaped by "
            f"both nature and experiences that transcend any single framework. The most interesting insights "
            f"come from where astrology and DNA disagree."
        )
    else:
        score_interp = (
            f"At {overall_score}%, your DNA tells a different story than your stars. This isn't a bad thing — "
            f"it means your genetic personality is uniquely your own, not easily captured by astrological archetypes. "
            f"The divergences are fascinating: where {sun_sign} predicts one tendency, your DNA reveals another. "
            f"You are, genetically speaking, a cosmic rebel."
        )

    # 8. Element and modality distribution
    all_signs = [chart[p]["sign"] for p in ["sun", "moon", "rising", "mercury", "venus", "mars", "jupiter", "saturn"]]
    elements = {"fire": 0, "earth": 0, "air": 0, "water": 0}
    modalities = {"cardinal": 0, "fixed": 0, "mutable": 0}
    for s in all_signs:
        elements[SIGN_ELEMENTS.get(s, "fire")] += 1
        modalities[SIGN_MODALITIES.get(s, "cardinal")] += 1

    # 9. Build birth chart section
    birth_chart = {
        "birthDate": date_str,
        "birthTime": time_str,
        "birthTimeKnown": birth_time_known,
        "birthCity": city,
        "sun": chart["sun"],
        "moon": chart["moon"],
        "rising": chart["rising"],
        "mercury": chart.get("mercury", {"sign": "Unknown", "degree": 0}),
        "venus": chart.get("venus", {"sign": "Unknown", "degree": 0}),
        "mars": chart.get("mars", {"sign": "Unknown", "degree": 0}),
        "jupiter": chart.get("jupiter", {"sign": "Unknown", "degree": 0}),
        "saturn": chart.get("saturn", {"sign": "Unknown", "degree": 0}),
        "elements": elements,
        "modalities": modalities,
    }

    disclaimers = [
        (
            "Astrology is not a science. Zodiac sign personality descriptions are based on "
            "cultural traditions, not empirical evidence. Large-scale studies have consistently "
            "found no correlation between birth date and personality traits (Shawn Carlson, Nature, 1985)."
        ),
        (
            "Genetic personality markers explain only a small fraction of personality variation. "
            "Environment, culture, upbringing, and personal choices shape who you are far more than "
            "any individual SNP. The genetic findings presented here are population-level associations, "
            "not individual predictions."
        ),
        (
            "The 'Cosmic Alignment Score' is an entertainment metric — it measures how well astrological "
            "archetypes happen to match your genetic markers. A high score does not validate astrology, "
            "and a low score does not mean your personality is wrong. It's a fun lens, not a diagnostic tool."
        ),
        (
            "This report is for entertainment and educational purposes only. It does not constitute "
            "psychological assessment, medical advice, or scientific personality profiling."
        ),
    ]
    if not birth_time_known:
        disclaimers.append(
            "Your birth time was not provided. The Rising Sign (Ascendant) and planetary house "
            "placements are estimated using noon as default. Sun sign and Moon sign are still "
            "accurate for the day."
        )

    return {
        "reportType": "stars_genes",
        "version": "1.0",
        "generatedAt": datetime.now(timezone.utc).isoformat(),
        "birthChart": birth_chart,
        "comparisons": {
            "sunSign": sun_comparison,
            "moonSign": moon_comparison,
            "risingSign": rising_comparison,
        },
        "planetaryInsights": planetary_insights,
        "cosmicAlignmentScore": {
            "overall": overall_score,
            "sunAlignment": sun_score,
            "moonAlignment": moon_score,
            "risingAlignment": rising_score,
            "interpretation": score_interp,
        },
        "disclaimers": disclaimers,
    }
