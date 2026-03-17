"""
ROH (Runs of Homozygosity) Analyzer
Detects runs of homozygosity using a PLINK-style sliding window algorithm.
ROH segments reveal ancestral consanguinity and population history.

Algorithm:
1. Group variants by autosome (chr1-22), sort by position
2. Slide a 50-SNP window; flag windows with <= 1 heterozygous SNP
3. Merge consecutive flagged SNPs into ROH segments
4. Filter: >= 1 Mb length, >= 50 SNPs, max 1 Mb gap between consecutive SNPs
5. Compute FROH = total ROH length / autosome genome length

References:
- McQuillan et al. (2008) Am J Hum Genet - ROH and inbreeding
- Pemberton et al. (2012) Am J Hum Genet - population FROH values
- Ceballos et al. (2018) Nat Rev Genet - consanguinity review
- Purcell et al. (2007) Am J Hum Genet - PLINK methodology
"""

from typing import Dict, List, Tuple, Any
from dataclasses import dataclass, field


# ─────────────────────────────────────────────────────────────────────────────
# CONSTANTS
# ─────────────────────────────────────────────────────────────────────────────

AUTOSOME_LENGTHS: Dict[int, int] = {
    1: 248956422, 2: 242193529, 3: 198295559, 4: 190214555, 5: 181538259,
    6: 170805979, 7: 159345973, 8: 145138636, 9: 138394717, 10: 133797422,
    11: 135086622, 12: 133275309, 13: 114364328, 14: 107043718, 15: 101991189,
    16: 90338345, 17: 83257441, 18: 80373285, 19: 58617616, 20: 64444167,
    21: 46709983, 22: 50818468,
}

TOTAL_AUTOSOME_LENGTH = 2_881_033_286  # sum of all autosome lengths (bp)

# Sliding window parameters
WINDOW_SIZE = 50          # SNPs per window
MAX_HET_PER_WINDOW = 1    # max heterozygous SNPs to flag window as homozygous

# ROH segment filters
MIN_ROH_LENGTH_BP = 1_000_000   # 1 Mb minimum
MIN_ROH_SNPS = 50               # minimum SNPs in segment
MAX_GAP_BP = 1_000_000          # max 1 Mb gap between consecutive SNPs

# ROH category boundaries (Mb)
SHORT_MAX_MB = 1.5
MEDIUM_MAX_MB = 5.0

# Population reference FROH values
POPULATION_REFERENCES = [
    ("Continental European (outbred)", 0.0005),
    ("UK / Western European", 0.001),
    ("Finnish / Scandinavian", 0.002),
    ("Sardinian", 0.004),
    ("Ashkenazi Jewish", 0.007),
    ("Middle Eastern isolated", 0.025),
    ("First-cousin parents", 0.0625),
]

# ROH category descriptions
CATEGORY_MEANINGS = {
    "short": "Shared ancestors 10-30 generations ago (~250-750 years)",
    "medium": "Shared ancestors 5-10 generations ago (~125-250 years)",
    "long": "Close parental relatedness (<5 generations / <125 years)",
}


# ─────────────────────────────────────────────────────────────────────────────
# DATA CLASSES
# ─────────────────────────────────────────────────────────────────────────────

@dataclass
class ROHSegment:
    chromosome: int
    start_bp: int
    end_bp: int
    snp_count: int
    het_count: int

    @property
    def length_bp(self) -> int:
        return self.end_bp - self.start_bp

    @property
    def length_mb(self) -> float:
        return self.length_bp / 1_000_000

    @property
    def het_rate(self) -> float:
        return self.het_count / self.snp_count if self.snp_count > 0 else 0.0

    @property
    def category(self) -> str:
        mb = self.length_mb
        if mb < SHORT_MAX_MB:
            return "short"
        elif mb < MEDIUM_MAX_MB:
            return "medium"
        else:
            return "long"


@dataclass
class ROHAnalysisResult:
    segments: List[ROHSegment]
    froh: float
    total_roh_bp: int
    snps_analyzed: int
    snps_per_chromosome: Dict[int, int]


# ─────────────────────────────────────────────────────────────────────────────
# HELPER FUNCTIONS
# ─────────────────────────────────────────────────────────────────────────────

def _is_heterozygous(genotype: str) -> bool:
    """Check if a genotype is heterozygous (two different alleles)."""
    g = genotype.upper().replace("-", "").strip()
    if len(g) != 2:
        return False
    return g[0] != g[1]


def _normalize_chrom(chrom: str) -> int:
    """
    Convert chromosome string to integer (autosomes only).
    Returns 0 for non-autosomal chromosomes.
    """
    c = chrom.upper().replace("CHR", "").strip()
    try:
        n = int(c)
        return n if 1 <= n <= 22 else 0
    except ValueError:
        return 0


def _parse_position(pos: str) -> int:
    """Parse a position string to integer."""
    try:
        return int(pos)
    except (ValueError, TypeError):
        return -1


# ─────────────────────────────────────────────────────────────────────────────
# SLIDING WINDOW ROH DETECTION
# ─────────────────────────────────────────────────────────────────────────────

def _detect_roh_on_chromosome(
    positions: List[int],
    is_het: List[bool],
) -> List[ROHSegment]:
    """
    Detect ROH segments on a single chromosome using PLINK-style sliding windows.

    Args:
        positions: Sorted list of genomic positions (bp)
        is_het: Parallel list indicating heterozygosity at each position

    Returns:
        List of ROHSegment objects passing all filters
    """
    n = len(positions)
    if n < WINDOW_SIZE:
        return []

    # Step 1: Flag each SNP if it appears in at least one homozygous window
    flagged = [False] * n

    # Pre-compute het count for the first window
    het_count = sum(1 for i in range(WINDOW_SIZE) if is_het[i])

    if het_count <= MAX_HET_PER_WINDOW:
        for i in range(WINDOW_SIZE):
            flagged[i] = True

    # Slide the window across the chromosome
    for start in range(1, n - WINDOW_SIZE + 1):
        # Remove the SNP leaving the window
        if is_het[start - 1]:
            het_count -= 1
        # Add the SNP entering the window
        if is_het[start + WINDOW_SIZE - 1]:
            het_count += 1

        if het_count <= MAX_HET_PER_WINDOW:
            for i in range(start, start + WINDOW_SIZE):
                flagged[i] = True

    # Step 2: Merge consecutive flagged SNPs into candidate segments
    segments: List[ROHSegment] = []
    seg_start_idx = None

    for i in range(n):
        if flagged[i]:
            if seg_start_idx is None:
                seg_start_idx = i
        else:
            if seg_start_idx is not None:
                segments.append(_build_segment(
                    positions, is_het, seg_start_idx, i - 1, 0
                ))
                seg_start_idx = None

    # Close final segment if it reaches the end
    if seg_start_idx is not None:
        segments.append(_build_segment(
            positions, is_het, seg_start_idx, n - 1, 0
        ))

    # Step 3: Split segments at large gaps (> MAX_GAP_BP between consecutive SNPs)
    split_segments: List[ROHSegment] = []
    for seg in segments:
        split_segments.extend(_split_at_gaps(positions, is_het, seg))

    # Step 4: Filter by minimum length and minimum SNP count
    filtered = [
        s for s in split_segments
        if s.length_bp >= MIN_ROH_LENGTH_BP and s.snp_count >= MIN_ROH_SNPS
    ]

    return filtered


def _build_segment(
    positions: List[int],
    is_het: List[bool],
    start_idx: int,
    end_idx: int,
    chrom: int,
) -> ROHSegment:
    """Build a ROHSegment from index range."""
    snp_count = end_idx - start_idx + 1
    het_count = sum(1 for i in range(start_idx, end_idx + 1) if is_het[i])
    return ROHSegment(
        chromosome=chrom,
        start_bp=positions[start_idx],
        end_bp=positions[end_idx],
        snp_count=snp_count,
        het_count=het_count,
    )


def _split_at_gaps(
    positions: List[int],
    is_het: List[bool],
    segment: ROHSegment,
) -> List[ROHSegment]:
    """
    Split a segment wherever the gap between consecutive SNPs exceeds MAX_GAP_BP.
    """
    # Find the index range for this segment
    start_idx = None
    end_idx = None
    for i, pos in enumerate(positions):
        if pos == segment.start_bp and start_idx is None:
            start_idx = i
        if pos == segment.end_bp:
            end_idx = i

    if start_idx is None or end_idx is None:
        return [segment]

    result: List[ROHSegment] = []
    sub_start = start_idx

    for i in range(start_idx + 1, end_idx + 1):
        if positions[i] - positions[i - 1] > MAX_GAP_BP:
            # Close current sub-segment
            if i - 1 >= sub_start:
                result.append(_build_segment(
                    positions, is_het, sub_start, i - 1, segment.chromosome
                ))
            sub_start = i

    # Close final sub-segment
    if end_idx >= sub_start:
        result.append(_build_segment(
            positions, is_het, sub_start, end_idx, segment.chromosome
        ))

    return result


# ─────────────────────────────────────────────────────────────────────────────
# CORE ANALYSIS
# ─────────────────────────────────────────────────────────────────────────────

def analyze_roh(
    variants: Dict[str, Tuple[str, str, str]]
) -> ROHAnalysisResult:
    """
    Analyze runs of homozygosity from genome-wide variant data.

    Args:
        variants: Dictionary mapping rsID to (chromosome, position, genotype).
                  Example: {"rs123": ("1", "12345", "AG"), ...}

    Returns:
        ROHAnalysisResult with segments, FROH, and chromosome distribution
    """
    # Step 1: Group variants by autosome and sort by position
    chrom_data: Dict[int, List[Tuple[int, bool]]] = {c: [] for c in range(1, 23)}

    for rsid, (chrom_str, pos_str, genotype) in variants.items():
        chrom = _normalize_chrom(chrom_str)
        if chrom == 0:
            continue

        pos = _parse_position(pos_str)
        if pos < 0:
            continue

        genotype = genotype.upper().replace("-", "").strip()
        if len(genotype) != 2:
            continue

        # Only consider valid biallelic SNPs (two nucleotide characters)
        if not all(c in "ACGT" for c in genotype):
            continue

        het = _is_heterozygous(genotype)
        chrom_data[chrom].append((pos, het))

    # Sort each chromosome by position
    for chrom in chrom_data:
        chrom_data[chrom].sort(key=lambda x: x[0])

    # Track SNP counts per chromosome
    snps_per_chrom: Dict[int, int] = {}
    total_snps = 0

    # Step 2: Detect ROH on each chromosome
    all_segments: List[ROHSegment] = []

    for chrom in range(1, 23):
        data = chrom_data[chrom]
        snps_per_chrom[chrom] = len(data)
        total_snps += len(data)

        if len(data) < WINDOW_SIZE:
            continue

        positions = [d[0] for d in data]
        is_het = [d[1] for d in data]

        segments = _detect_roh_on_chromosome(positions, is_het)

        # Assign chromosome number to segments
        for seg in segments:
            seg.chromosome = chrom

        all_segments.extend(segments)

    # Step 3: Calculate FROH
    total_roh_bp = sum(seg.length_bp for seg in all_segments)
    froh = total_roh_bp / TOTAL_AUTOSOME_LENGTH if TOTAL_AUTOSOME_LENGTH > 0 else 0.0

    return ROHAnalysisResult(
        segments=all_segments,
        froh=froh,
        total_roh_bp=total_roh_bp,
        snps_analyzed=total_snps,
        snps_per_chromosome=snps_per_chrom,
    )


# ─────────────────────────────────────────────────────────────────────────────
# INTERPRETATION
# ─────────────────────────────────────────────────────────────────────────────

def _interpret_froh(froh: float) -> str:
    """Generate a human-readable interpretation of the FROH value."""
    if froh < 0.001:
        return (
            "Your FROH is very low, typical of outbred continental populations. "
            "There is no evidence of recent ancestral consanguinity in your genome."
        )
    elif froh < 0.003:
        return (
            "Your FROH is low, typical of Western and Northern European populations. "
            "This suggests a broadly outbred ancestry with no significant recent consanguinity."
        )
    elif froh < 0.01:
        return (
            "Your FROH indicates moderate ancestral homozygosity, consistent with "
            "populations that have experienced some historical isolation or endogamy."
        )
    elif froh < 0.03:
        return (
            "Your FROH is elevated, consistent with endogamous ancestry or a "
            "historically isolated community. This level of homozygosity suggests "
            "shared ancestors within the last several generations."
        )
    else:
        return (
            "Your FROH is high, suggesting close parental relatedness within recent "
            "generations. This level of homozygosity is above what is typically observed "
            "in most world populations."
        )


# ─────────────────────────────────────────────────────────────────────────────
# JSON OUTPUT
# ─────────────────────────────────────────────────────────────────────────────

def generate_roh_json(result: ROHAnalysisResult) -> dict:
    """
    Generate a JSON-serializable dict for the frontend ROH report.

    Args:
        result: ROHAnalysisResult from analyze_roh()

    Returns:
        Dict matching the ROH report JSON schema
    """
    # Categorize segments
    categories = {"short": [], "medium": [], "long": []}
    for seg in result.segments:
        categories[seg.category].append(seg)

    # Build category summary
    categories_json = {}
    for cat_name in ("short", "medium", "long"):
        segs = categories[cat_name]
        total_mb = sum(s.length_mb for s in segs)
        categories_json[cat_name] = {
            "count": len(segs),
            "totalMb": round(total_mb, 2),
            "meaning": CATEGORY_MEANINGS[cat_name],
        }

    # Population comparison (insert user at the correct sorted position)
    user_froh = result.froh
    pop_entries = [
        {"population": name, "froh": ref_froh, "isUser": False}
        for name, ref_froh in POPULATION_REFERENCES
    ]
    user_entry = {"population": "You", "froh": round(user_froh, 5), "isUser": True}
    pop_entries.append(user_entry)
    pop_entries.sort(key=lambda x: x["froh"])

    # Top segments (sorted by length, descending)
    sorted_segments = sorted(result.segments, key=lambda s: s.length_bp, reverse=True)
    top_segments = []
    for seg in sorted_segments[:10]:  # top 10 longest
        top_segments.append({
            "chromosome": seg.chromosome,
            "startMb": round(seg.start_bp / 1_000_000, 2),
            "endMb": round(seg.end_bp / 1_000_000, 2),
            "lengthMb": round(seg.length_mb, 2),
            "snpCount": seg.snp_count,
            "hetRate": round(seg.het_rate, 4),
        })

    # Chromosome distribution
    chrom_dist = []
    for chrom in range(1, 23):
        chrom_segs = [s for s in result.segments if s.chromosome == chrom]
        if chrom_segs:
            total_mb = sum(s.length_mb for s in chrom_segs)
            chrom_dist.append({
                "chromosome": chrom,
                "segmentCount": len(chrom_segs),
                "totalMb": round(total_mb, 2),
            })

    # Summary statistics
    total_roh_mb = result.total_roh_bp / 1_000_000
    total_segments = len(result.segments)
    longest_mb = sorted_segments[0].length_mb if sorted_segments else 0.0
    mean_mb = total_roh_mb / total_segments if total_segments > 0 else 0.0

    return {
        "reportType": "roh",
        "version": "1.0",
        "summary": {
            "froh": round(result.froh, 5),
            "frohPercentage": round(result.froh * 100, 3),
            "totalRohMb": round(total_roh_mb, 2),
            "totalSegments": total_segments,
            "longestRohMb": round(longest_mb, 2),
            "meanRohMb": round(mean_mb, 2),
            "snpsAnalyzed": result.snps_analyzed,
            "interpretation": _interpret_froh(result.froh),
        },
        "categories": categories_json,
        "populationComparison": pop_entries,
        "topSegments": top_segments,
        "chromosomeDistribution": chrom_dist,
    }
