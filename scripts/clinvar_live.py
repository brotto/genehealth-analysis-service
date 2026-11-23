"""
Live variant clinical significance enrichment (Biopython feature #2).

Despite the module name, this fetches clinical significance from TWO NCBI
databases, in this order of priority:

  1. **dbSNP** (primary) — rich built-in `clinical_significance` field,
     population frequency hints, function class, gene link. Every modern
     rsID is present; coverage is close to 100%.
  2. **ClinVar** (optional) — richer per-variant interpretation: title,
     condition list, review stars, last-evaluated date, deep-link to the
     variation page. Used as a second pass when dbSNP signals an
     interesting (pathogenic / risk factor / drug response) variant.

Output per rsID:

    {
        "rs429358": {
            "rsid": "rs429358",
            "gene": "APOE",
            "snpClass": "snv",
            "significances": [
                "pathogenic", "risk-factor", "drug-response", ...
            ],
            "topSignificance": "Pathogenic",             # UI-friendly label
            "severity": "high",                           # high|moderate|info
            "functionalClass": "missense_variant",
            "dbsnpUrl": "https://www.ncbi.nlm.nih.gov/snp/rs429358",
            # ClinVar metadata (may be absent if the extra lookup failed)
            "clinvarId": "441268",
            "clinvarTitle": "NM_000041.3(APOE):c.[137T>C;388T>C]",
            "clinvarUrl": "https://www.ncbi.nlm.nih.gov/clinvar/variation/441268/",
            "conditions": ["APOE4(-)-FREIBURG", "Familial hypercholesterolemia"],
            "reviewStatus": "criteria provided, multiple submitters, no conflicts",
            "reviewStars": 2,
            "lastEvaluated": "2024-08-06",
            # Alleles on the GRCh38 plus strand (dbSNP SPDI). variantAllele is
            # the ALT the ClinVar record classifies (or dbSNP's only ALT).
            "refAllele": "T",
            "altAlleles": ["C"],
            "variantAllele": "C",
            # Present (True) only when a ClinVar record classifies the
            # reference allele itself as non-benign ("c.388=").
            "referenceAlleleClassified": True,
            # Present (True) only when the ClinVar search hit its page limit,
            # so not every record could be checked.
            "referenceCheckIncomplete": True
        }
    }

The user's own genotype is added afterwards by
`clinical_genotype.annotate_user_genotypes()` (userGenotype, zygosity, ...).

Any missing field is dropped. Full failure returns an empty dict and logs
a warning — it must never block the analyzer pipeline.
"""

from __future__ import annotations

import logging
import re
from typing import Any, Dict, Iterable, List, Optional, Tuple

from .clinical_genotype import (
    apply_clinvar_alleles,
    is_reference_equal_spdi,
    parse_dbsnp_spdi,
    parse_spdi,
)
from .ncbi_client import NCBIClient, default_client

log = logging.getLogger(__name__)

DBSNP_URL = "https://www.ncbi.nlm.nih.gov/snp/{}"
CLINVAR_URL = "https://www.ncbi.nlm.nih.gov/clinvar/variation/{}/"

# Order matters: first match in the comma-separated `clinical_significance`
# string determines the UI-facing top label + severity tier.
_SIGNIFICANCE_ORDER: List[Tuple[str, str, str]] = [
    # (pattern, UI-friendly label, severity tier)
    ("pathogenic-pathogenic-low-penetrance", "Pathogenic (low-penetrance)", "high"),
    ("pathogenic-established-risk-allele", "Pathogenic (risk allele)", "high"),
    ("likely-pathogenic", "Likely pathogenic", "high"),
    ("pathogenic", "Pathogenic", "high"),
    ("risk-factor", "Risk factor", "high"),
    ("drug-response", "Drug response", "moderate"),
    ("association", "Disease association", "moderate"),
    ("conflicting-interpretations-of-pathogenicity", "Conflicting interpretations", "moderate"),
    ("uncertain-significance", "Uncertain significance", "info"),
    ("likely-benign", "Likely benign", "info"),
    ("benign", "Benign", "info"),
    ("protective", "Protective", "info"),
]

_REVIEW_STARS: Dict[str, int] = {
    "practice guideline": 4,
    "reviewed by expert panel": 3,
    "criteria provided, multiple submitters, no conflicts": 2,
    "criteria provided, conflicting interpretations": 1,
    "criteria provided, conflicting classifications": 1,
    "criteria provided, single submitter": 1,
    "no assertion criteria provided": 0,
    "no assertion provided": 0,
}

# Severities we bother enriching further with a ClinVar deep-dive.
_NOTABLE_SEVERITIES = {"high", "moderate"}

# ClinVar germline descriptions that do NOT make a reference-equal record
# ("c.388=") a reason to doubt a homozygous-reference call. Anything else
# (pathogenic, risk factor, drug response, uncertain, conflicting, empty...)
# does.
_NON_CLASSIFYING_DESCRIPTIONS = {
    "benign",
    "likely benign",
    "benign/likely benign",
    "not provided",
    "no classification provided",
    "no classification for the single variant",
    "no classifications from unflagged records",
}

# How many ClinVar records to read per rsID. Pharmacogenomic sites carry
# 30-40 records; when the search fills the whole page we cannot be sure we
# saw every record, so a reference call is reported as unconfirmed.
_CLINVAR_SEARCH_MAX = 50


def _normalize_rsid(rs: str) -> Optional[str]:
    if not rs:
        return None
    rs = rs.strip().lower()
    if rs.startswith("rs") and rs[2:].isdigit():
        return rs
    return None


def _tier_for(sig_csv: str) -> Tuple[str, str]:
    """Return (ui_label, severity) for a CSV of significance strings.

    Matches on whole tokens only — comparing against the comma-separated
    list rather than substring-matching the whole string. This prevents
    `pathogenic` from falsely matching inside
    `conflicting-interpretations-of-pathogenicity`.
    """
    tokens = {t.strip().lower() for t in (sig_csv or "").split(",") if t.strip()}
    if not tokens:
        return ("", "info")
    for pattern, label, severity in _SIGNIFICANCE_ORDER:
        if pattern in tokens:
            return label, severity
    return ("", "info")


def _get_ci(rec: Dict[str, Any], *keys: str) -> Any:
    """Case-insensitive field fetch — NCBI Entrez.read returns UPPER_CASE
    keys for dbSNP while the JSON API uses lower_case. We normalize."""
    for k in keys:
        if k in rec:
            return rec[k]
        ku = k.upper()
        if ku in rec:
            return rec[ku]
        kl = k.lower()
        if kl in rec:
            return rec[kl]
    return None


def _parse_dbsnp_record(rsid: str, rec: Dict[str, Any]) -> Optional[Dict[str, Any]]:
    """Turn a dbSNP esummary record into the output dict.

    Entrez XML parse surfaces keys in UPPER_CASE (CLINICAL_SIGNIFICANCE,
    GENES, FXN_CLASS, SNP_CLASS, NAME inside genes). The JSON API uses
    snake_case. We handle both via `_get_ci`.
    """
    try:
        sig_csv = str(_get_ci(rec, "clinical_significance") or "")
        significances = [s.strip() for s in sig_csv.split(",") if s.strip()]
        top_label, severity = _tier_for(sig_csv)

        gene = ""
        genes = _get_ci(rec, "genes") or []
        if isinstance(genes, list) and genes:
            first = genes[0]
            if isinstance(first, dict):
                gene = str(_get_ci(first, "name") or "").strip()

        fxn = str(_get_ci(rec, "fxn_class") or "").strip()
        fxn_primary = fxn.split(",")[0] if fxn else ""

        out: Dict[str, Any] = {
            "rsid": rsid,
            "gene": gene,
            "snpClass": str(_get_ci(rec, "snp_class") or "").strip(),
            "significances": significances,
            "topSignificance": top_label,
            "severity": severity,
            "functionalClass": fxn_primary,
            "dbsnpUrl": DBSNP_URL.format(rsid),
        }

        # Reference / alternate alleles on the GRCh38 plus strand. Needed to
        # tell homozygous-reference from homozygous-alternate later.
        alleles = parse_dbsnp_spdi(_get_ci(rec, "spdi"))
        if alleles:
            out["refAllele"] = alleles["ref"]
            out["altAlleles"] = list(alleles["alts"])
            if len(alleles["alts"]) == 1:
                out["variantAllele"] = alleles["alts"][0]
            out["_dbsnp"] = alleles  # internal; removed before returning
        return out
    except (KeyError, TypeError, ValueError) as e:
        log.warning("malformed dbSNP record for %s: %s", rsid, e)
        return None


def _fetch_dbsnp(
    client: NCBIClient,
    rsids: List[str],
) -> Dict[str, Dict[str, Any]]:
    """Batch-fetch dbSNP summaries for up to ~200 rsIDs at once."""
    numeric_ids = [rs[2:] for rs in rsids]  # strip `rs` prefix
    id_to_rs = dict(zip(numeric_ids, rsids))

    try:
        records = client.esummary(db="snp", ids=numeric_ids)
    except RuntimeError as e:
        log.warning("dbSNP esummary failed: %s", e)
        return {}

    out: Dict[str, Dict[str, Any]] = {}
    for rec in records:
        uid = str(_get_ci(rec, "uid") or _get_ci(rec, "Id") or _get_ci(rec, "SNP_ID") or "").strip()
        # SNP_ID may have the numeric form; map back to our rs prefixed form
        rsid = id_to_rs.get(uid)
        if not rsid:
            continue
        parsed = _parse_dbsnp_record(rsid, rec)
        if parsed:
            out[rsid] = parsed
    return out


# ── ClinVar deep-dive (optional second pass) ───────────────────────────────


def _extract_dbsnp_xrefs(record: Dict[str, Any]) -> List[str]:
    ids: List[str] = []
    vs = record.get("variation_set")
    if not isinstance(vs, list):
        return ids
    for v in vs:
        if not isinstance(v, dict):
            continue
        xrefs = v.get("variation_xrefs")
        if not isinstance(xrefs, list):
            continue
        for x in xrefs:
            if isinstance(x, dict) and str(x.get("db_source") or "").lower() == "dbsnp":
                db_id = str(x.get("db_id") or "").strip()
                if db_id:
                    ids.append(db_id)
    return ids


def _extract_entries(record: Dict[str, Any]) -> List[Tuple[List[str], str]]:
    """(dbSNP ids, canonical SPDI) for every variation_set entry.

    Haplotype records ("c.[137T>C;388T>C]") carry one entry per variant, so
    the allele for *our* rsID must come from the entry whose xref matches.
    """
    out: List[Tuple[List[str], str]] = []
    vs = record.get("variation_set")
    if not isinstance(vs, list):
        return out
    for v in vs:
        if not isinstance(v, dict):
            continue
        spdi = str(v.get("canonical_spdi") or "").strip()
        ids: List[str] = []
        xrefs = v.get("variation_xrefs")
        if isinstance(xrefs, list):
            for x in xrefs:
                if isinstance(x, dict) and str(x.get("db_source") or "").lower() == "dbsnp":
                    db_id = str(x.get("db_id") or "").strip()
                    if db_id:
                        ids.append(db_id)
        out.append((ids, spdi))
    return out


def _site_spdis(
    parsed: Dict[str, Any],
    rs_digits: str,
    site: Optional[Dict[str, Any]] = None,
) -> List[str]:
    """SPDIs of the record's entries at our position: matched by dbSNP xref,
    or (when the dbSNP site is known) by sequence + position."""
    found: List[str] = []
    for ids, spdi in parsed.get("_entries") or []:
        if not spdi:
            continue
        at_site = rs_digits in ids
        if not at_site and site:
            loc = parse_spdi(spdi)
            at_site = bool(loc) and (loc[0], loc[1]) == (site.get("seq"), site.get("pos"))
        if at_site and spdi not in found:
            found.append(spdi)
    return found


def _is_reference_only_record(parsed: Dict[str, Any], rs_digits: str, site=None) -> bool:
    """True when the record's entry at our position describes the reference
    allele itself ("c.388=", SPDI T:T) rather than a change."""
    spdis = _site_spdis(parsed, rs_digits, site)
    return bool(spdis) and all(is_reference_equal_spdi(s) for s in spdis)


def _classifies_reference(parsed: Dict[str, Any], rs_digits: str, site=None) -> bool:
    """Does this record give the reference allele at our position a
    non-benign classification?

    Only records where *every* entry is reference-equal count. A haplotype
    record such as CYP2C9 "c.[430C=;1075A>C]" classifies the 1075A>C change;
    the reference base at 430 is context, not what is being classified.
    """
    desc = str(parsed.get("_description") or "").strip().lower()
    if desc in _NON_CLASSIFYING_DESCRIPTIONS:
        return False
    if not _is_reference_only_record(parsed, rs_digits, site):
        return False
    all_spdis = [spdi for _ids, spdi in parsed.get("_entries") or [] if spdi]
    return all(is_reference_equal_spdi(s) for s in all_spdis)


def _parse_clinvar_record(record: Dict[str, Any]) -> Optional[Dict[str, Any]]:
    try:
        cvid = str(record.get("uid") or record.get("Id") or "").strip()
        if not cvid:
            return None
        title = str(record.get("title") or record.get("variation_name") or "").strip()

        gc = record.get("germline_classification") or {}
        if not isinstance(gc, dict):
            gc = {}
        desc = str(gc.get("description") or record.get("clinical_significance") or "").strip()
        review_status = str(gc.get("review_status") or record.get("review_status") or "").strip()
        last_eval = str(gc.get("last_evaluated") or record.get("last_evaluated") or "").strip()
        trait_set = gc.get("trait_set") or record.get("trait_set") or []

        conditions: List[str] = []
        if isinstance(trait_set, list):
            for t in trait_set:
                if isinstance(t, dict):
                    name = str(t.get("trait_name") or "").strip()
                    if name and name.lower() not in ("not provided", "not specified"):
                        conditions.append(name)
        # dedupe preserving order, cap 5
        seen: set[str] = set()
        uniq_conditions: List[str] = []
        for c in conditions:
            if c not in seen:
                seen.add(c)
                uniq_conditions.append(c)
            if len(uniq_conditions) >= 5:
                break

        stars = _REVIEW_STARS.get(review_status.lower().strip(), 0)
        # Normalize "2024/08/06 00:00" → "2024-08-06"
        eval_match = re.match(r"(\d{4})/(\d{2})/(\d{2})", last_eval)
        iso_date = (
            f"{eval_match.group(1)}-{eval_match.group(2)}-{eval_match.group(3)}"
            if eval_match else last_eval
        )

        return {
            "clinvarId": cvid,
            "clinvarTitle": title,
            "clinvarUrl": CLINVAR_URL.format(cvid),
            "conditions": uniq_conditions,
            "reviewStatus": review_status,
            "reviewStars": stars,
            "lastEvaluated": iso_date,
            "_dbsnp_xrefs": _extract_dbsnp_xrefs(record),
            "_entries": _extract_entries(record),
            "_description": desc,
        }
    except (KeyError, TypeError, ValueError) as e:
        log.warning("malformed ClinVar summary: %s", e)
        return None


def _score_clinvar(parsed: Dict[str, Any], rs_digits: str) -> int:
    score = 0
    # Must match the rsID via dbSNP xref; huge bonus.
    if rs_digits in (parsed.get("_dbsnp_xrefs") or []):
        score += 1000
    stars = parsed.get("reviewStars") or 0
    if isinstance(stars, int):
        score += stars * 30
    if parsed.get("conditions"):
        score += 50 + min(len(parsed["conditions"]), 5) * 5
    title = (parsed.get("clinvarTitle") or "")
    if ">" in title and "=" not in title:
        score += 15  # prefer substitution over silent
    return score


def _fetch_clinvar(
    client: NCBIClient,
    rsid: str,
    site: Optional[Dict[str, Any]] = None,
) -> Optional[Dict[str, Any]]:
    """Look up a single rsID in ClinVar, validated via dbSNP xref match.

    `site` is the parsed dbSNP SPDI ({seq, pos, ref, alts}) when known; it lets
    reference-allele records be recognised by position as well as by xref.
    """
    rs_digits = rsid.lower().lstrip("rs")

    ids: List[str] = []
    for term in (rsid, f"{rsid}*"):
        try:
            ids = client.esearch(db="clinvar", term=term, retmax=_CLINVAR_SEARCH_MAX)
        except RuntimeError as e:
            log.warning("ClinVar esearch failed rsid=%s term=%r: %s", rsid, term, e)
            continue
        if ids:
            break
    if not ids:
        return None

    try:
        records = client.esummary(db="clinvar", ids=ids)
    except RuntimeError as e:
        log.warning("ClinVar esummary failed for ids=%s: %s", ids, e)
        return None

    parsed_all: List[Dict[str, Any]] = []
    candidates: List[Dict[str, Any]] = []
    for rec in records:
        parsed = _parse_clinvar_record(rec)
        if not parsed:
            continue
        parsed_all.append(parsed)
        # Only keep records whose dbSNP xref matches our rsID.
        if rs_digits not in (parsed.get("_dbsnp_xrefs") or []):
            continue
        candidates.append(parsed)

    if not candidates:
        return None

    # Show a record that describes a change at this position whenever one
    # exists: a reference-only record ("c.1601=") on the card would make
    # "you don't carry the variant" describe the wrong allele.
    candidates.sort(
        key=lambda p: (not _is_reference_only_record(p, rs_digits, site), _score_clinvar(p, rs_digits)),
        reverse=True,
    )

    # Does any record classify the *reference* allele here as something
    # other than benign? Then two reference copies are not reassuring.
    ref_classified = any(_classifies_reference(p, rs_digits, site) for p in parsed_all)

    best = candidates[0]
    # If even the record shown on the card describes the reference allele
    # here (only reference-only records exist), its classification is about
    # the reference allele, whatever other changes a haplotype includes.
    if _is_reference_only_record(best, rs_digits, site):
        best_desc = str(best.get("_description") or "").strip().lower()
        if best_desc not in _NON_CLASSIFYING_DESCRIPTIONS:
            ref_classified = True
    best["_clinvarSpdis"] = _site_spdis(best, rs_digits, site)
    if ref_classified:
        best["referenceAlleleClassified"] = True
    if len(ids) >= _CLINVAR_SEARCH_MAX:
        best["referenceCheckIncomplete"] = True
    for key in ("_dbsnp_xrefs", "_entries", "_description"):
        best.pop(key, None)
    return best


# ── Public API ─────────────────────────────────────────────────────────────


def enrich_with_clinvar(
    rsids: Iterable[str],
    client: Optional[NCBIClient] = None,
    max_rsids: int = 50,
    include_clinvar_detail: bool = True,
) -> Dict[str, Dict[str, Any]]:
    """
    Look up live clinical significance for up to `max_rsids` rsIDs.

    Returns keyed metadata — caller treats missing keys as "use the static
    snapshot interpretation" (fallback).

    Never raises: any crash per-rsID is logged and skipped.
    """
    cli = client or default_client()

    # Normalize + dedupe
    seen: set[str] = set()
    ordered: List[str] = []
    for raw in rsids:
        norm = _normalize_rsid(raw)
        if norm and norm not in seen:
            seen.add(norm)
            ordered.append(norm)
        if len(ordered) >= max_rsids:
            break
    if not ordered:
        return {}

    # Phase 1: batch dbSNP lookup (one request, up to ~200 rsids).
    out = _fetch_dbsnp(cli, ordered)

    # Phase 2: for "notable" variants, do an extra ClinVar lookup to get
    # conditions + review status + last_evaluated.
    if include_clinvar_detail:
        for rsid, data in out.items():
            severity = data.get("severity")
            if severity not in _NOTABLE_SEVERITIES:
                continue
            try:
                cv = _fetch_clinvar(cli, rsid, data.get("_dbsnp"))
                if cv:
                    clinvar_spdis = cv.pop("_clinvarSpdis", None)
                    data.update(cv)
                    # Prefer the exact allele ClinVar classified over the
                    # dbSNP single-ALT fallback set in phase 1.
                    apply_clinvar_alleles(data, data.get("_dbsnp"), clinvar_spdis)
            except Exception as e:  # noqa: BLE001
                log.warning("ClinVar deep-dive crashed for %s: %s", rsid, e)
                continue

    for data in out.values():
        data.pop("_dbsnp", None)
    return out


def extract_risk_rsids_from_report(report: Dict[str, Any]) -> List[str]:
    """Pull every rsID mentioned in a disease-risk style report."""
    collected: List[str] = []

    def _walk(node: Any) -> None:
        if isinstance(node, dict):
            for k, v in node.items():
                kl = k.lower()
                if kl in ("rsid", "rs_id", "variant", "snp"):
                    if isinstance(v, str):
                        if _normalize_rsid(v):
                            collected.append(v.strip())
                    elif isinstance(v, list):
                        for it in v:
                            if isinstance(it, str) and _normalize_rsid(it):
                                collected.append(it.strip())
                else:
                    _walk(v)
        elif isinstance(node, list):
            for it in node:
                _walk(it)

    _walk(report)

    seen: set[str] = set()
    uniq: List[str] = []
    for rs in collected:
        key = rs.lower()
        if key not in seen:
            seen.add(key)
            uniq.append(rs)
    return uniq
