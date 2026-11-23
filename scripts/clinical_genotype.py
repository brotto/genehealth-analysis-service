"""
Interpret the user's own genotype at each Live Clinical Data position.

`clinvar_live.enrich_with_clinvar()` lists ClinVar-indexed positions found in
the user's file. On its own that list reads like "you carry all of these",
which is how a user misread 14 Lynch-syndrome cards in May 2026. This module
adds, per position, what the user's file actually shows:

    userGenotype    the genotype string from the file ("AG", "--", ...)
    zygosity        homozygous_ref | heterozygous | homozygous_alt |
                    no_call | unknown
    likelyCarrier   True when the file shows at least one non-reference copy,
                    False for homozygous_ref, None when not determinable
    strandFlipped   True when the file reports the opposite DNA strand
    zygosityReason  short machine code explaining an `unknown`

Safety rule: `homozygous_ref` is the only value the UI renders as reassuring,
so it is emitted only when it is established:

  * the reference and variant alleles come from NCBI (dbSNP SPDI on the GRCh38
    plus strand, ClinVar `canonical_spdi` for the classified allele), never
    from guessing that "two identical letters" means reference;
  * the genotype is read both as reported and complemented (opposite strand).
    If both readings fit the alleles but disagree, the answer is `unknown`.
    That is always the case for homozygotes at A/T and C/G sites, where the
    reference on one strand is the variant on the other;
  * if the file shows a real dbSNP allele that is not the one ClinVar
    classified, we do not "fix" it by flipping strands: `unknown`;
  * if ClinVar classifies the reference allele itself (records written as
    "c.388=" with SPDI T:T), two reference copies are not reassuring:
    `unknown`;
  * a reference call also needs a matching ClinVar record to have been
    checked, completely, whenever the card carries any risk-type label;
    otherwise `unknown`.

Assumption, stated plainly: when the ClinVar allele is known, the genotype is
compared against {reference, that allele} only, because consumer chips report
just the two alleles their probe was designed for. So "GG" at a site where
dbSNP also lists a rare C allele is read as reference, even though a
minus-strand file would show a plus-strand rare-C homozygote the same way.
That reading never involves the classified allele (a genotype containing it,
on either strand, is never called reference), which is the property the UI
depends on. Strand itself is read from the alleles: a file is never assumed
to be on the plus strand when the letters fit both strands.
"""

from __future__ import annotations

import logging
from typing import Any, Dict, Iterable, List, Mapping, Optional, Sequence, Tuple

log = logging.getLogger(__name__)

HOMOZYGOUS_REF = "homozygous_ref"
HETEROZYGOUS = "heterozygous"
HOMOZYGOUS_ALT = "homozygous_alt"
NO_CALL = "no_call"
UNKNOWN = "unknown"

ZYGOSITY_VALUES = (HOMOZYGOUS_REF, HETEROZYGOUS, HOMOZYGOUS_ALT, NO_CALL, UNKNOWN)

# Severity tiers (from clinvar_live) for which a reference call must be
# backed by a checked ClinVar record.
_NOTABLE_SEVERITIES = {"high", "moderate"}

# dbSNP significance labels that never call for that check. Any other label
# (including ones clinvar_live's tier table does not know, such as
# "pathogenic-low-penetrance" or "established-risk-allele") does.
_NON_RISK_SIGNIFICANCES = {
    "benign",
    "likely-benign",
    "benign-likely-benign",
    "not-provided",
    "uncertain-significance",
    "no-classifications-from-unflagged-records",
    "protective",
}

_COMPLEMENT = {"A": "T", "T": "A", "C": "G", "G": "C"}
_BASES = frozenset(_COMPLEMENT)
_NO_CALL_CHARS = frozenset("-0N?.")
_NO_CALL_TOKENS = frozenset({"NC", "NOCALL", "NA"})
_SEPARATORS = str.maketrans("", "", "/| \t\r\n")

# Internal marker: a reading with two different non-reference alleles.
_MULTI_ALT = "_multi_alt"


# ── SPDI helpers ───────────────────────────────────────────────────────────


def parse_spdi(spdi: Any) -> Optional[Tuple[str, int, str, str]]:
    """Parse one SPDI string "NC_000006.12:26092912:G:A" → (seq, pos, del, ins)."""
    if not isinstance(spdi, str):
        return None
    parts = spdi.strip().split(":")
    if len(parts) != 4:
        return None
    seq, pos, deleted, inserted = parts
    try:
        pos_i = int(pos)
    except ValueError:
        return None
    if not seq:
        return None
    return seq, pos_i, deleted.strip().upper(), inserted.strip().upper()


def parse_dbsnp_spdi(spdi_field: Any) -> Optional[Dict[str, Any]]:
    """Turn dbSNP's comma-separated SPDI list into reference + alternate alleles.

    "NC_000006.12:26092912:G:A,NC_000006.12:26092912:G:T"
        → {"seq": "NC_000006.12", "pos": 26092912, "ref": "G", "alts": ["A", "T"]}

    Returns None unless every entry is a single-base substitution at the same
    position with the same reference base. Anything else (indels, MNVs,
    disagreeing entries) is not something a genotype letter pair can be
    compared against.
    """
    if not isinstance(spdi_field, str) or not spdi_field.strip():
        return None
    seq: Optional[str] = None
    pos: Optional[int] = None
    ref: Optional[str] = None
    alts: List[str] = []
    for chunk in spdi_field.split(","):
        parsed = parse_spdi(chunk)
        if not parsed:
            return None
        s, p, deleted, inserted = parsed
        if len(deleted) != 1 or len(inserted) != 1:
            return None
        if deleted not in _BASES or inserted not in _BASES or deleted == inserted:
            return None
        if seq is None:
            seq, pos, ref = s, p, deleted
        elif (s, p, deleted) != (seq, pos, ref):
            return None
        if inserted not in alts:
            alts.append(inserted)
    if ref is None or not alts:
        return None
    return {"seq": seq, "pos": pos, "ref": ref, "alts": alts}


_ALLELE_FIELDS = ("refAllele", "altAlleles", "variantAllele")


def apply_clinvar_alleles(
    payload: Dict[str, Any],
    dbsnp: Optional[Mapping[str, Any]],
    clinvar_spdis: Optional[Sequence[str]],
) -> None:
    """Merge the allele ClinVar classified into a card's allele fields, in place.

    `clinvar_spdis` are the canonical SPDIs of the chosen ClinVar record's
    entry for this rsID (normally exactly one).

      * reference-equal record ("c.388=")  → the card describes the reference
        allele, so there is no "variant allele" to name: variantAllele is
        dropped (the reference call is handled by referenceAlleleClassified).
      * clean SNV agreeing with dbSNP      → variantAllele = its ALT (added to
        altAlleles if dbSNP did not list it).
      * anything else (indel, MNV, a different position or reference base, or
        several SPDIs for the same rsID) → the two NCBI sources do not agree on
        what this card is about, so all allele fields are dropped and the card
        is classified `unknown` rather than risk a wrong reading.
    """
    spdis = [s for s in (clinvar_spdis or []) if isinstance(s, str) and s.strip()]
    if not spdis or not dbsnp:
        return
    if len(spdis) == 1 and is_reference_equal_spdi(spdis[0]):
        payload.pop("variantAllele", None)
        return

    parsed = parse_spdi(spdis[0]) if len(spdis) == 1 else None
    consistent = (
        parsed is not None
        and (parsed[0], parsed[1], parsed[2]) == (dbsnp.get("seq"), dbsnp.get("pos"), dbsnp.get("ref"))
        and len(parsed[3]) == 1
        and parsed[3] in _BASES
        and parsed[3] != parsed[2]
    )
    if not consistent:
        for key in _ALLELE_FIELDS:
            payload.pop(key, None)
        return

    alt = parsed[3]
    alts = list(payload.get("altAlleles") or [])
    if alt not in alts:
        alts.append(alt)
    payload["altAlleles"] = alts
    payload["variantAllele"] = alt


def is_reference_equal_spdi(spdi: Any) -> bool:
    """True for ClinVar records that describe the reference allele itself
    (HGVS "c.388=", SPDI "...:T:T")."""
    parsed = parse_spdi(spdi)
    return bool(parsed) and parsed[2] == parsed[3] and parsed[2] != ""


# ── Genotype classification ────────────────────────────────────────────────


def _complement(genotype: str) -> str:
    return "".join(_COMPLEMENT[b] for b in genotype)


def _normalize(genotype: Any) -> str:
    if not isinstance(genotype, str):
        return ""
    return genotype.translate(_SEPARATORS).upper()


def _is_no_call(g: str) -> bool:
    return (not g) or g in _NO_CALL_TOKENS or all(c in _NO_CALL_CHARS for c in g)


def _reading_class(reading: str, ref: str, allowed: frozenset) -> Optional[str]:
    """Classify one strand reading of the genotype, or None if its alleles do
    not fit the allowed set."""
    if not set(reading) <= allowed:
        return None
    n_ref = reading.count(ref)
    if n_ref == len(reading):
        return HOMOZYGOUS_REF
    if len(reading) == 2 and n_ref == 1:
        return HETEROZYGOUS
    if n_ref == 0 and len(set(reading)) == 1:
        return HOMOZYGOUS_ALT
    return _MULTI_ALT


def _result(
    zygosity: str,
    reason: Optional[str] = None,
    strand_flipped: bool = False,
) -> Dict[str, Any]:
    # likelyCarrier is only a boolean when the genotype settles it; for
    # no_call / unknown it is None so nobody reads "False" as "not a carrier".
    if zygosity in (HETEROZYGOUS, HOMOZYGOUS_ALT):
        carrier: Optional[bool] = True
    elif zygosity == HOMOZYGOUS_REF:
        carrier = False
    else:
        carrier = None
    out: Dict[str, Any] = {"zygosity": zygosity, "likelyCarrier": carrier}
    if strand_flipped:
        out["strandFlipped"] = True
    if reason:
        out["zygosityReason"] = reason
    return out


def classify_genotype(
    genotype: Any,
    ref: Optional[str],
    alts: Optional[Sequence[str]],
    variant_allele: Optional[str] = None,
    reference_block: Optional[str] = None,
) -> Dict[str, Any]:
    """Classify a genotype string against NCBI's alleles for that position.

    Args:
        genotype: genotype as it appears in the user's file ("AG", "--", "A").
        ref: reference base on the GRCh38 plus strand (dbSNP SPDI).
        alts: every alternate base dbSNP lists for the position.
        variant_allele: the specific ALT ClinVar classified, when known.
        reference_block: reason code when two reference copies would not be
            reassuring or cannot be confirmed (see module docstring). A
            homozygous-reference result is then reported as `unknown` with
            this code as `zygosityReason`.

    Returns a dict with `zygosity`, `likelyCarrier` and, when relevant,
    `strandFlipped` / `zygosityReason`.
    """
    g = _normalize(genotype)
    if _is_no_call(g):
        return _result(NO_CALL)
    if len(g) not in (1, 2) or any(c not in _BASES for c in g):
        # Indel calls (I/D), multi-base strings, or garbage.
        return _result(UNKNOWN, "not_a_single_base_genotype")

    ref = (ref or "").upper()
    if ref not in _BASES:
        return _result(UNKNOWN, "no_allele_data")
    site_alleles = frozenset({ref} | {a.upper() for a in (alts or []) if a and a.upper() in _BASES})

    variant = (variant_allele or "").upper()
    if variant in _BASES and variant != ref and variant in site_alleles:
        allowed = frozenset({ref, variant})
        plus = _reading_class(g, ref, allowed)
        if plus is None and set(g) <= site_alleles:
            # The file shows a real allele at this site that is not the one
            # ClinVar classified. Reinterpreting it as the opposite strand
            # could turn a non-reference call into a false "reference".
            return _result(UNKNOWN, "other_allele")
    else:
        if len(site_alleles) < 2:
            return _result(UNKNOWN, "no_allele_data")
        allowed = site_alleles
        plus = _reading_class(g, ref, allowed)

    minus = _reading_class(_complement(g), ref, allowed)

    candidates = {c for c in (plus, minus) if c is not None}
    if not candidates:
        return _result(UNKNOWN, "allele_mismatch")
    if len(candidates) > 1:
        return _result(UNKNOWN, "strand_ambiguous")

    zygosity = candidates.pop()
    flipped = plus is None
    if zygosity == _MULTI_ALT:
        return _result(UNKNOWN, "two_different_variant_alleles", flipped)
    if zygosity == HOMOZYGOUS_REF and reference_block:
        return _result(UNKNOWN, reference_block, flipped)
    return _result(zygosity, None, flipped)


# ── Pipeline hook ──────────────────────────────────────────────────────────


def _reference_block(payload: Mapping[str, Any]) -> Optional[str]:
    """Why a homozygous-reference call must not be shown as reassuring."""
    if payload.get("referenceAlleleClassified") is True:
        # ClinVar classifies the reference allele itself (non-benign).
        return "reference_allele_classified"
    if payload.get("referenceCheckIncomplete") is True:
        # Too many ClinVar records to read them all.
        return "reference_not_confirmed"
    if not payload.get("clinvarId"):
        sigs = {str(x).strip().lower() for x in (payload.get("significances") or [])}
        if payload.get("severity") in _NOTABLE_SEVERITIES or (sigs - _NON_RISK_SIGNIFICANCES):
            # Risk-type card with no matching ClinVar record checked: we
            # cannot rule out that the reference allele is the classified one.
            return "reference_not_confirmed"
    return None


def annotate_variant(payload: Dict[str, Any], genotype: Optional[str]) -> None:
    """Add userGenotype / zygosity / likelyCarrier (+ optional fields) in place."""
    for key in ("strandFlipped", "zygosityReason"):
        payload.pop(key, None)

    cleaned = _normalize(genotype)
    payload["userGenotype"] = cleaned or None
    if genotype is None:
        payload.update(_result(NO_CALL, "not_in_file"))
        return

    payload.update(
        classify_genotype(
            cleaned,
            payload.get("refAllele"),
            payload.get("altAlleles"),
            payload.get("variantAllele"),
            reference_block=_reference_block(payload),
        )
    )


def annotate_user_genotypes(
    clinical_live: Dict[str, Dict[str, Any]],
    variants: Mapping[str, Tuple[str, str, str]],
) -> None:
    """Annotate every clinical_live card with the user's genotype, in place.

    `variants` is the parser's {rsid_lower: (chrom, pos, genotype)} dict.
    A failure on one card marks that card `unknown` instead of aborting.
    """
    for rsid, payload in clinical_live.items():
        if not isinstance(payload, dict):
            continue
        try:
            record = variants.get(str(rsid).lower())
            genotype = record[2] if record else None
            annotate_variant(payload, genotype)
        except Exception as e:  # noqa: BLE001
            log.warning("genotype annotation failed for %s: %s", rsid, e)
            payload["zygosity"] = UNKNOWN
            payload["likelyCarrier"] = None
            payload["zygosityReason"] = "annotation_error"
            payload.setdefault("userGenotype", None)


def summarize(clinical_live: Mapping[str, Mapping[str, Any]]) -> Dict[str, int]:
    """Count cards per zygosity (for logs)."""
    counts = {z: 0 for z in ZYGOSITY_VALUES}
    for payload in clinical_live.values():
        z = payload.get("zygosity") if isinstance(payload, Mapping) else None
        if z in counts:
            counts[z] += 1
    return counts


__all__: Iterable[str] = (
    "ZYGOSITY_VALUES",
    "annotate_user_genotypes",
    "annotate_variant",
    "classify_genotype",
    "apply_clinvar_alleles",
    "is_reference_equal_spdi",
    "parse_dbsnp_spdi",
    "parse_spdi",
    "summarize",
)
