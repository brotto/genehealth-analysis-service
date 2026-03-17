"""
Evolution Analyzer
Analyzes genetic variants related to recent human evolution, including adaptive sweeps
(natural selection events in the last ~50,000 years) and archaic introgression
(DNA inherited from Neanderthals and Denisovans).

References:
- Bersaglieri et al. (2004) AJHG - LCT lactase persistence sweep
- Norton et al. (2007) Mol Biol Evol - SLC24A5 skin pigmentation
- Kamberov et al. (2013) Cell - EDAR hair/teeth/sweat glands
- MacArthur et al. (2007) Nat Genet - ACTN3 R577X muscle fiber
- Stephens et al. (1998) AJHG - CCR5-delta32 HIV resistance
- Valverde et al. (1995) Nat Genet - MC1R red hair/fair skin
- Mendez et al. (2012) Mol Biol Evol - OAS1 Neanderthal introgression
- Zeberg & Paabo (2020) Nature - Chr3 COVID risk Neanderthal haplotype
- Huerta-Sanchez et al. (2014) Nature - EPAS1 Denisovan altitude adaptation
- Mendez et al. (2013) PLoS Genet - STAT2 Neanderthal introgression
- Abi-Rached et al. (2011) Science - HLA Neanderthal introgression
- Vernot & Akey (2014) Science - BNC2 Neanderthal skin pigmentation
"""

from typing import Dict, List, Tuple, Any
from dataclasses import dataclass, field


# ---------------------------------------------------------------------------
# COMPLEMENT MAP (for strand-flip handling)
# ---------------------------------------------------------------------------

COMPLEMENT = {"A": "T", "T": "A", "C": "G", "G": "C"}


def complement_allele(allele: str) -> str:
    """Return the complement of a single nucleotide."""
    return COMPLEMENT.get(allele.upper(), allele.upper())


def count_allele(genotype: str, allele: str) -> int:
    """
    Count occurrences of an allele in a genotype, handling complement strands.
    For ambiguous SNPs (A/T or C/G), checks both the allele and its complement.
    """
    genotype = genotype.upper().replace("-", "")
    allele = allele.upper()
    comp = complement_allele(allele)

    count = genotype.count(allele)

    # If the allele is A/T or C/G (ambiguous), complement is the same pair
    # so we should NOT double-count. Only count complement if it differs
    # from the original allele and the genotype doesn't contain the original.
    if count == 0 and comp != allele:
        count = genotype.count(comp)

    return min(count, 2)


# ---------------------------------------------------------------------------
# LOCUS DEFINITIONS
# ---------------------------------------------------------------------------

@dataclass
class EvolutionSNP:
    rsid: str
    derived_allele: str


@dataclass
class EvolutionLocus:
    id: str
    name: str
    category: str               # "adaptive_sweep" or "archaic_introgression"
    snps: List[EvolutionSNP]
    description: str
    phenotype_icon: str
    geographic_origin: str
    evidence_level: str         # "Strong", "Moderate", "Emerging"
    # Sweep-specific
    age_years: int = 0
    # Introgression-specific
    source: str = ""            # "Neanderthal", "Denisovan", or ""
    note: str = ""              # optional note (e.g., array limitations)


# ---------------------------------------------------------------------------
# ADAPTIVE SWEEPS DATABASE
# ---------------------------------------------------------------------------

ADAPTIVE_SWEEPS: List[EvolutionLocus] = [
    EvolutionLocus(
        id="LCT",
        name="Lactase Persistence",
        category="adaptive_sweep",
        snps=[
            EvolutionSNP(rsid="rs4988235", derived_allele="T"),
            EvolutionSNP(rsid="rs182549", derived_allele="T"),
        ],
        age_years=8000,
        geographic_origin="Northern Europe",
        evidence_level="Strong",
        phenotype_icon="milk",
        description=(
            "Lactase persistence allows adults to digest lactose in milk, a trait "
            "that arose with the spread of dairy farming in Northern Europe around "
            "8,000 years ago. This is one of the strongest and most recent signals "
            "of positive selection in the human genome."
        ),
    ),
    EvolutionLocus(
        id="SLC24A5",
        name="Skin Lightening",
        category="adaptive_sweep",
        snps=[
            EvolutionSNP(rsid="rs1426654", derived_allele="A"),
        ],
        age_years=10000,
        geographic_origin="Europe",
        evidence_level="Strong",
        phenotype_icon="sun",
        description=(
            "The SLC24A5 gene variant is responsible for a major portion of skin "
            "color difference between European and African populations. The derived "
            "allele swept to near-fixation in Europeans between 8,000 and 10,000 "
            "years ago, likely driven by the need for vitamin D synthesis at "
            "northern latitudes."
        ),
    ),
    EvolutionLocus(
        id="EDAR",
        name="Hair, Teeth & Sweat Glands",
        category="adaptive_sweep",
        snps=[
            EvolutionSNP(rsid="rs3827760", derived_allele="C"),
        ],
        age_years=30000,
        geographic_origin="East Asia",
        evidence_level="Strong",
        phenotype_icon="sparkles",
        description=(
            "The EDAR V370A variant affects hair thickness, tooth shape (shovel-shaped "
            "incisors), and sweat gland density. It rose to high frequency in East Asian "
            "populations around 30,000 years ago. This single variant has pleiotropic "
            "effects on multiple ectodermal structures."
        ),
    ),
    EvolutionLocus(
        id="ACTN3",
        name="Muscle Fiber Type (R577X)",
        category="adaptive_sweep",
        snps=[
            EvolutionSNP(rsid="rs1815739", derived_allele="T"),
        ],
        age_years=50000,
        geographic_origin="Global",
        evidence_level="Moderate",
        phenotype_icon="muscle",
        description=(
            "The ACTN3 R577X variant causes loss of alpha-actinin-3 protein in "
            "fast-twitch muscle fibers. The X allele (T) has increased in frequency "
            "globally, especially in populations that migrated to colder climates, "
            "possibly because it improves metabolic efficiency and cold tolerance."
        ),
    ),
    EvolutionLocus(
        id="CCR5",
        name="HIV Resistance (CCR5-delta32)",
        category="adaptive_sweep",
        snps=[
            EvolutionSNP(rsid="rs333", derived_allele="del32"),
        ],
        age_years=5000,
        geographic_origin="Northern Europe",
        evidence_level="Moderate",
        phenotype_icon="shield",
        description=(
            "The CCR5-delta32 deletion confers strong resistance to HIV infection "
            "in homozygous carriers. It rose to ~10% frequency in Northern Europeans, "
            "possibly selected by historical plague or smallpox epidemics. Note: most "
            "genotyping arrays cannot detect this deletion directly."
        ),
        note="Most genotyping arrays cannot reliably detect this deletion.",
    ),
    EvolutionLocus(
        id="MC1R",
        name="Red Hair & Fair Skin",
        category="adaptive_sweep",
        snps=[
            EvolutionSNP(rsid="rs1805007", derived_allele="T"),
            EvolutionSNP(rsid="rs1805008", derived_allele="C"),
        ],
        age_years=50000,
        geographic_origin="Europe",
        evidence_level="Strong",
        phenotype_icon="palette",
        description=(
            "MC1R variants reduce eumelanin production, leading to red hair, fair "
            "skin, and freckling. Multiple independent MC1R variants have risen in "
            "frequency in European populations over the last 50,000 years, likely "
            "due to relaxed selection on dark pigmentation at higher latitudes."
        ),
    ),
]

# ---------------------------------------------------------------------------
# ARCHAIC INTROGRESSION DATABASE
# ---------------------------------------------------------------------------

ARCHAIC_INTROGRESSION: List[EvolutionLocus] = [
    EvolutionLocus(
        id="OAS1",
        name="OAS1 Antiviral Immunity",
        category="archaic_introgression",
        snps=[
            EvolutionSNP(rsid="rs10774671", derived_allele="G"),
        ],
        source="Neanderthal",
        geographic_origin="Eurasia",
        evidence_level="Strong",
        phenotype_icon="shield",
        description=(
            "The Neanderthal-derived OAS1 haplotype enhances antiviral defense by "
            "producing a more active form of the 2'-5'-oligoadenylate synthetase "
            "enzyme. This introgressed variant has been positively selected in "
            "modern humans, likely due to its protective effects against RNA viruses."
        ),
    ),
    EvolutionLocus(
        id="CHR3_COVID",
        name="Chr3 COVID-19 Risk Locus",
        category="archaic_introgression",
        snps=[
            EvolutionSNP(rsid="rs35044562", derived_allele="A"),
            EvolutionSNP(rsid="rs10490770", derived_allele="T"),
        ],
        source="Neanderthal",
        geographic_origin="South Asia / Europe",
        evidence_level="Strong",
        phenotype_icon="virus",
        description=(
            "A Neanderthal-derived haplotype on chromosome 3 is the strongest "
            "genetic risk factor for severe COVID-19. Carriers have approximately "
            "doubled risk of requiring hospitalization. This haplotype is most "
            "common in South Asian populations (~50%) and present at ~16% in Europeans."
        ),
    ),
    EvolutionLocus(
        id="EPAS1",
        name="EPAS1 Altitude Adaptation",
        category="archaic_introgression",
        snps=[
            EvolutionSNP(rsid="rs1868092", derived_allele="A"),
        ],
        source="Denisovan",
        geographic_origin="Tibet",
        evidence_level="Strong",
        phenotype_icon="mountain",
        description=(
            "The EPAS1 variant inherited from Denisovans helps Tibetans thrive at "
            "high altitude by regulating hemoglobin production. This is the most "
            "celebrated example of adaptive introgression, where archaic DNA provided "
            "a crucial survival advantage in extreme environments."
        ),
    ),
    EvolutionLocus(
        id="STAT2",
        name="STAT2 Immune Regulation",
        category="archaic_introgression",
        snps=[
            EvolutionSNP(rsid="rs2066807", derived_allele="T"),
        ],
        source="Neanderthal",
        geographic_origin="Eurasia",
        evidence_level="Moderate",
        phenotype_icon="dna",
        description=(
            "A Neanderthal-derived STAT2 haplotype affects innate immune signaling "
            "through the JAK-STAT pathway. This variant modulates the interferon "
            "response and may influence susceptibility to viral infections, "
            "representing a trade-off between immune activation and autoimmunity."
        ),
    ),
    EvolutionLocus(
        id="HLA_MHC",
        name="MHC/HLA Immune Region",
        category="archaic_introgression",
        snps=[
            EvolutionSNP(rsid="rs2596988", derived_allele="T"),
        ],
        source="Neanderthal",
        geographic_origin="Eurasia",
        evidence_level="Moderate",
        phenotype_icon="antibody",
        description=(
            "Neanderthal-derived HLA alleles in the major histocompatibility complex "
            "region have been retained at high frequency in modern Eurasians. These "
            "variants help the immune system recognize pathogens and may have provided "
            "early modern humans with pre-adapted immunity to Eurasian infections."
        ),
    ),
    EvolutionLocus(
        id="BNC2",
        name="BNC2 Skin Pigmentation",
        category="archaic_introgression",
        snps=[
            EvolutionSNP(rsid="rs10756819", derived_allele="G"),
        ],
        source="Neanderthal",
        geographic_origin="Europe",
        evidence_level="Moderate",
        phenotype_icon="sun",
        description=(
            "A Neanderthal-derived variant near BNC2 influences skin pigmentation "
            "and freckling in modern Europeans. This suggests that introgression "
            "from Neanderthals contributed to the lighter skin pigmentation that "
            "evolved in European populations after they left Africa."
        ),
    ),
]


# ---------------------------------------------------------------------------
# ANALYSIS RESULT DATACLASSES
# ---------------------------------------------------------------------------

@dataclass
class LocusSNPResult:
    rsid: str
    genotype: str
    derived_copies: int     # 0, 1, or 2
    found: bool


@dataclass
class LocusResult:
    locus: EvolutionLocus
    snp_results: List[LocusSNPResult]
    status: str             # "homozygous", "heterozygous", "absent", "not_typed"
    interpretation: str


@dataclass
class EvolutionAnalysisResult:
    sweeps: List[LocusResult]
    introgression: List[LocusResult]
    total_snps_checked: int
    total_snps_found: int


# ---------------------------------------------------------------------------
# INTERPRETATION BUILDERS
# ---------------------------------------------------------------------------

def _build_sweep_interpretation(locus: EvolutionLocus, status: str) -> str:
    """Build a human-readable interpretation for an adaptive sweep locus."""
    name = locus.name
    age_label = f"~{locus.age_years:,}" if locus.age_years else "unknown"

    if status == "homozygous":
        return (
            f"You carry two copies of the derived allele for {name}. "
            f"This adaptive variant, which arose approximately {age_label} years ago "
            f"in {locus.geographic_origin}, is fully present in your genome."
        )
    elif status == "heterozygous":
        return (
            f"You carry one copy of the derived allele for {name}. "
            f"This adaptive variant, selected approximately {age_label} years ago "
            f"in {locus.geographic_origin}, is partially present in your genome."
        )
    elif status == "not_typed":
        note = f" {locus.note}" if locus.note else ""
        return (
            f"The SNP(s) for {name} were not found in your genotyping data, "
            f"so this locus could not be assessed.{note}"
        )
    else:
        return (
            f"You do not carry the derived allele for {name}. "
            f"This adaptive variant, which arose in {locus.geographic_origin}, "
            f"is not detected in your genome."
        )


def _build_introgression_interpretation(locus: EvolutionLocus, status: str) -> str:
    """Build a human-readable interpretation for an archaic introgression locus."""
    name = locus.name
    source = locus.source

    if status == "homozygous":
        return (
            f"You carry two copies of the {source}-derived allele at {name}. "
            f"This archaic variant is strongly present in your genome, suggesting "
            f"this introgressed segment was retained on both chromosomes."
        )
    elif status == "heterozygous":
        return (
            f"You carry one copy of the {source}-derived allele at {name}. "
            f"This indicates you inherited this archaic DNA segment from one parent."
        )
    elif status == "not_typed":
        return (
            f"The SNP(s) for {name} were not found in your genotyping data, "
            f"so this introgression locus could not be assessed."
        )
    else:
        return (
            f"You do not carry the {source}-derived allele at {name}. "
            f"This archaic introgression variant is not detected in your genome."
        )


# ---------------------------------------------------------------------------
# CORE ANALYSIS
# ---------------------------------------------------------------------------

def _analyze_locus(
    locus: EvolutionLocus,
    variants: Dict[str, Tuple[str, str, str]],
) -> LocusResult:
    """Analyze a single evolution locus against user variants."""
    snp_results: List[LocusSNPResult] = []
    max_derived = 0
    any_found = False

    for snp_def in locus.snps:
        rsid_lower = snp_def.rsid.lower()
        if rsid_lower in variants:
            _chrom, _pos, genotype = variants[rsid_lower]
            genotype = genotype.upper().replace("-", "")
            derived = count_allele(genotype, snp_def.derived_allele)
            snp_results.append(LocusSNPResult(
                rsid=snp_def.rsid,
                genotype=genotype,
                derived_copies=derived,
                found=True,
            ))
            max_derived = max(max_derived, derived)
            any_found = True
        else:
            snp_results.append(LocusSNPResult(
                rsid=snp_def.rsid,
                genotype="",
                derived_copies=0,
                found=False,
            ))

    # Determine status
    if not any_found:
        status = "not_typed"
    elif max_derived >= 2:
        status = "homozygous"
    elif max_derived == 1:
        status = "heterozygous"
    else:
        status = "absent"

    # Build interpretation
    if locus.category == "adaptive_sweep":
        interpretation = _build_sweep_interpretation(locus, status)
    else:
        interpretation = _build_introgression_interpretation(locus, status)

    return LocusResult(
        locus=locus,
        snp_results=snp_results,
        status=status,
        interpretation=interpretation,
    )


def analyze_evolution(
    variants: Dict[str, Tuple[str, str, str]],
) -> EvolutionAnalysisResult:
    """
    Analyze genetic variants for evolutionary signals (adaptive sweeps
    and archaic introgression).

    Args:
        variants: Dictionary mapping rsID (lowercase) to (chromosome, position, genotype)

    Returns:
        EvolutionAnalysisResult with sweep and introgression findings
    """
    sweeps: List[LocusResult] = []
    introgression: List[LocusResult] = []

    total_snps = 0
    found_snps = 0

    for locus in ADAPTIVE_SWEEPS:
        result = _analyze_locus(locus, variants)
        sweeps.append(result)
        for sr in result.snp_results:
            total_snps += 1
            if sr.found:
                found_snps += 1

    for locus in ARCHAIC_INTROGRESSION:
        result = _analyze_locus(locus, variants)
        introgression.append(result)
        for sr in result.snp_results:
            total_snps += 1
            if sr.found:
                found_snps += 1

    return EvolutionAnalysisResult(
        sweeps=sweeps,
        introgression=introgression,
        total_snps_checked=total_snps,
        total_snps_found=found_snps,
    )


# ---------------------------------------------------------------------------
# JSON OUTPUT
# ---------------------------------------------------------------------------

def generate_evolution_json(result: EvolutionAnalysisResult) -> dict:
    """
    Generate a JSON-serializable dict for the frontend evolution report.

    Args:
        result: EvolutionAnalysisResult from analyze_evolution()

    Returns:
        Dict matching the evolution report JSON schema
    """
    # Build notable findings
    notable_findings: List[str] = []
    for lr in result.sweeps + result.introgression:
        if lr.status in ("homozygous", "heterozygous"):
            label = "homozygous (2 copies)" if lr.status == "homozygous" else "heterozygous (1 copy)"
            notable_findings.append(f"{lr.locus.name}: {label}")

    # Adaptive sweeps JSON
    adaptive_sweeps_json: List[dict] = []
    for lr in result.sweeps:
        snps_json = []
        for sr in lr.snp_results:
            snps_json.append({
                "rsid": sr.rsid,
                "genotype": sr.genotype if sr.found else None,
                "derivedCopies": sr.derived_copies,
                "found": sr.found,
            })

        adaptive_sweeps_json.append({
            "id": lr.locus.id,
            "name": lr.locus.name,
            "ageYears": lr.locus.age_years,
            "geographicOrigin": lr.locus.geographic_origin,
            "evidence": lr.locus.evidence_level,
            "description": lr.locus.description,
            "phenotypeIcon": lr.locus.phenotype_icon,
            "snps": snps_json,
            "status": lr.status,
            "interpretation": lr.interpretation,
            "note": lr.locus.note if lr.locus.note else None,
        })

    # Archaic introgression JSON
    archaic_introgression_json: List[dict] = []
    for lr in result.introgression:
        snps_json = []
        for sr in lr.snp_results:
            snps_json.append({
                "rsid": sr.rsid,
                "genotype": sr.genotype if sr.found else None,
                "derivedCopies": sr.derived_copies,
                "found": sr.found,
            })

        archaic_introgression_json.append({
            "id": lr.locus.id,
            "name": lr.locus.name,
            "source": lr.locus.source,
            "evidence": lr.locus.evidence_level,
            "description": lr.locus.description,
            "phenotypeIcon": lr.locus.phenotype_icon,
            "geographicOrigin": lr.locus.geographic_origin,
            "snps": snps_json,
            "status": lr.status,
            "interpretation": lr.interpretation,
        })

    return {
        "reportType": "evolution",
        "version": "1.0",
        "summary": {
            "sweepsAnalyzed": len(ADAPTIVE_SWEEPS),
            "introgressionAnalyzed": len(ARCHAIC_INTROGRESSION),
            "snpsUsed": result.total_snps_found,
            "snpsTotal": result.total_snps_checked,
            "notableFindings": notable_findings,
        },
        "adaptiveSweeps": adaptive_sweeps_json,
        "archaicIntrogression": archaic_introgression_json,
    }
