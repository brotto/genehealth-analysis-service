"""
PubMed evidence enrichment for GeneHealth reports.

Takes a list of genes and/or rsIDs and returns recent peer-reviewed papers for
each, structured so the frontend can render a "Scientific Evidence" section.

Typical call from an analyzer:

    from scripts.pubmed_enrichment import enrich_with_pubmed

    scientific_evidence = enrich_with_pubmed(
        genes=["APOE", "MTHFR"],
        rsids=["rs429358", "rs1801133"],
        per_item=3,
        disease_context="longevity",  # optional — narrows query
    )
    report_json["scientificEvidence"] = scientific_evidence

The return shape is always a list of dicts:

    [
        {
            "key": "APOE",
            "type": "gene",
            "query": "APOE[Gene] AND longevity AND 2022:2026[DP]",
            "references": [
                {"pmid": "12345", "title": "...", "journal": "...",
                 "year": "2024", "authors": ["Doe J", "Smith A"],
                 "doi": "10.xxxx/...", "url": "https://pubmed.ncbi.nlm.nih.gov/12345/"}
            ]
        },
        ...
    ]

Design notes:
  - Queries favor recent papers (default: last 4 years) and sort by pub_date.
  - rsID queries check both the rsID literal and the associated gene if provided.
  - Results cached at the NCBIClient layer (30d TTL by default).
  - If NCBI is unreachable, the function returns whatever it could fetch and
    logs the failure — NEVER raises, because a failed enrichment should not
    block the analysis pipeline.
"""

from __future__ import annotations

import logging
from datetime import datetime
from typing import Any, Dict, Iterable, List, Optional

from .ncbi_client import NCBIClient, default_client

log = logging.getLogger(__name__)

PUBMED_URL = "https://pubmed.ncbi.nlm.nih.gov/{}/"


def _build_gene_query(gene: str, disease_context: Optional[str], years_back: int) -> str:
    year_start = datetime.utcnow().year - years_back
    year_end = datetime.utcnow().year
    base = f'"{gene}"[Gene Name] OR {gene}[Title/Abstract]'
    if disease_context:
        base = f'({base}) AND ("{disease_context}"[Title/Abstract])'
    return f'({base}) AND ({year_start}:{year_end}[DP])'


def _build_rsid_query(rsid: str, gene: Optional[str], years_back: int) -> str:
    year_start = datetime.utcnow().year - years_back
    year_end = datetime.utcnow().year
    if gene:
        base = f'({rsid}[Title/Abstract] OR "{gene}"[Gene Name])'
    else:
        base = f'{rsid}[Title/Abstract]'
    return f'({base}) AND ({year_start}:{year_end}[DP])'


def _parse_summary(record: Dict[str, Any]) -> Optional[Dict[str, Any]]:
    """Convert an Entrez esummary record to a lightweight dict."""
    try:
        pmid = str(record.get("Id") or record.get("uid") or "")
        if not pmid:
            return None

        title = str(record.get("Title") or "").strip().rstrip(".")
        journal = str(record.get("FullJournalName") or record.get("Source") or "").strip()

        pub_date = str(record.get("PubDate") or record.get("EPubDate") or "")
        year = pub_date.split(" ")[0][:4] if pub_date else ""

        raw_authors = record.get("AuthorList") or []
        authors: List[str] = []
        for a in raw_authors[:4]:
            name = str(a).strip()
            if name:
                authors.append(name)
        if len(raw_authors) > 4:
            authors.append("et al.")

        doi = ""
        article_ids = record.get("ArticleIds") or {}
        if isinstance(article_ids, dict):
            doi = str(article_ids.get("doi") or "").strip()

        return {
            "pmid": pmid,
            "title": title,
            "journal": journal,
            "year": year,
            "authors": authors,
            "doi": doi,
            "url": PUBMED_URL.format(pmid),
        }
    except (KeyError, TypeError, ValueError) as e:
        log.warning("malformed esummary record: %s (%s)", record, e)
        return None


def _fetch_for_query(
    client: NCBIClient,
    query: str,
    retmax: int,
) -> List[Dict[str, Any]]:
    try:
        ids = client.esearch(db="pubmed", term=query, retmax=retmax, sort="pub_date")
    except RuntimeError as e:
        log.warning("esearch failed for query=%r: %s", query, e)
        return []

    if not ids:
        return []

    try:
        records = client.esummary(db="pubmed", ids=ids)
    except RuntimeError as e:
        log.warning("esummary failed for ids=%r: %s", ids, e)
        return []

    refs: List[Dict[str, Any]] = []
    for rec in records:
        parsed = _parse_summary(rec)
        if parsed:
            refs.append(parsed)
    return refs


def enrich_with_pubmed(
    genes: Optional[Iterable[str]] = None,
    rsids: Optional[Iterable[str]] = None,
    rsid_gene_map: Optional[Dict[str, str]] = None,
    per_item: int = 3,
    years_back: int = 4,
    disease_context: Optional[str] = None,
    client: Optional[NCBIClient] = None,
) -> List[Dict[str, Any]]:
    """
    Build scientific-evidence payload for a set of genes and/or rsIDs.

    Args:
        genes: list of gene symbols (e.g. ["APOE", "MTHFR"]).
        rsids: list of rsIDs (e.g. ["rs429358"]).
        rsid_gene_map: optional mapping rsID -> gene symbol, used to widen queries.
        per_item: max references per gene/rsID (default 3).
        years_back: how many years into the past to search (default 4).
        disease_context: optional keyword to narrow gene queries
            (e.g. "longevity", "pharmacogenomics").
        client: inject a custom NCBIClient (useful for tests).

    Returns:
        List of {key, type, query, references[]} dicts. Empty list on total failure.
    """
    cli = client or default_client()
    genes = list(genes or [])
    rsids = list(rsids or [])
    rsid_gene_map = rsid_gene_map or {}

    results: List[Dict[str, Any]] = []
    seen_keys: set[str] = set()

    for gene in genes:
        key = gene.upper()
        if key in seen_keys:
            continue
        seen_keys.add(key)

        query = _build_gene_query(gene, disease_context, years_back)
        refs = _fetch_for_query(cli, query, per_item)
        if refs:
            results.append({
                "key": key,
                "type": "gene",
                "query": query,
                "references": refs,
            })

    for rsid in rsids:
        key = rsid.lower()
        if key in seen_keys:
            continue
        seen_keys.add(key)

        gene_hint = rsid_gene_map.get(rsid) or rsid_gene_map.get(key)
        query = _build_rsid_query(rsid, gene_hint, years_back)
        refs = _fetch_for_query(cli, query, per_item)
        if refs:
            results.append({
                "key": key,
                "type": "rsid",
                "query": query,
                "references": refs,
            })

    return results


# -------------------------------------------------------------------- helpers

def extract_genes_and_rsids_from_report(report: Dict[str, Any]) -> Dict[str, List[str]]:
    """
    Walk a report JSON and pull out every gene symbol and rsID mentioned.

    Uses a shallow scan over common keys used across the analyzers:
      - 'gene', 'Gene', 'gene_symbol'
      - 'rsid', 'rsID', 'variant', 'snp'
    """
    genes: List[str] = []
    rsids: List[str] = []

    def _walk(node: Any) -> None:
        if isinstance(node, dict):
            for k, v in node.items():
                kl = k.lower()
                if kl in ("gene", "gene_symbol", "genes"):
                    if isinstance(v, str):
                        genes.append(v)
                    elif isinstance(v, list):
                        genes.extend([str(x) for x in v if x])
                elif kl in ("rsid", "variant", "snp", "rs_id"):
                    if isinstance(v, str) and v.lower().startswith("rs"):
                        rsids.append(v)
                    elif isinstance(v, list):
                        rsids.extend([str(x) for x in v if isinstance(x, str) and x.lower().startswith("rs")])
                else:
                    _walk(v)
        elif isinstance(node, list):
            for item in node:
                _walk(item)

    _walk(report)

    # Deduplicate while preserving order.
    def _dedupe(seq: Iterable[str]) -> List[str]:
        seen: set[str] = set()
        out: List[str] = []
        for x in seq:
            key = x.strip().upper()
            if key and key not in seen:
                seen.add(key)
                out.append(x.strip())
        return out

    return {"genes": _dedupe(genes), "rsids": _dedupe(rsids)}


# -------------------------------------------------------------------- orchestrator

import os as _os


def maybe_enrich_report_inplace(
    report: Dict[str, Any],
    *,
    disease_context: Optional[str] = None,
    gene_cap: int = 12,
    rsid_cap: int = 8,
    per_item: int = 3,
    years_back: int = 4,
    feature_label: str = "pubmed",
) -> None:
    """
    End-to-end helper: gate via env var, extract genes/rsids, enrich, assign
    `scientificEvidence` to the report dict (mutates in place).

    Never raises — any failure is logged and the report is left unchanged
    (no `scientificEvidence` key), so the analyzer pipeline keeps flowing.

    Used by all analyzers that want the PubMed Evidence treatment:
    precision_medicine, longevity_aging, mental_wellbeing, etc.
    """
    gate = _os.environ.get("ENABLE_PUBMED_ENRICHMENT", "true").lower()
    if gate == "false":
        print(f"[{feature_label}] disabled via ENABLE_PUBMED_ENRICHMENT=false", flush=True)
        return

    try:
        extracted = extract_genes_and_rsids_from_report(report)
        genes = extracted["genes"][:gene_cap]
        rsids = extracted["rsids"][:rsid_cap]
        print(
            f"[{feature_label}] extracted {len(genes)} genes + {len(rsids)} rsids "
            f"(context={disease_context!r})",
            flush=True,
        )
        if not genes and not rsids:
            return

        evidence = enrich_with_pubmed(
            genes=genes,
            rsids=rsids,
            per_item=per_item,
            years_back=years_back,
            disease_context=disease_context,
        )
        if evidence:
            report["scientificEvidence"] = evidence
            print(
                f"[{feature_label}] enrichment complete: {len(evidence)} items with refs",
                flush=True,
            )
    except Exception as e:  # noqa: BLE001
        import traceback as _tb
        print(
            f"[{feature_label}] enrichment skipped: {type(e).__name__}: {e}",
            flush=True,
        )
        _tb.print_exc()
