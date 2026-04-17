"""
NCBI Entrez client with rate-limiting, retry/backoff, and on-disk cache.

Wraps Bio.Entrez to make it safe to call from analyzers without worrying about:
  - API key + email setup
  - NCBI rate limits (3 req/s without key, 10 req/s with key)
  - Transient errors (HTTP 429, 500, timeouts)
  - Redundant calls (cached by request hash)

Usage:
    from scripts.ncbi_client import NCBIClient
    client = NCBIClient()
    ids = client.esearch("pubmed", 'APOE[Gene] AND 2024:2026[DP]', retmax=5)
    summaries = client.esummary("pubmed", ids)

The cache lives at NCBI_CACHE_DIR (default: /tmp/ncbi_cache/). In production on
EasyPanel, this is ephemeral per-container; in a future iteration we'll persist
to R2 via the Worker callback path instead.
"""

from __future__ import annotations

import hashlib
import json
import logging
import os
import random
import threading
import time
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional

try:
    from dotenv import load_dotenv  # type: ignore

    load_dotenv()
except Exception:
    # dotenv not installed in all environments — continue with env vars only.
    pass

from Bio import Entrez
from urllib.error import HTTPError, URLError

log = logging.getLogger(__name__)


class NCBIClient:
    """Thread-safe NCBI Entrez wrapper with cache + rate-limit + retry."""

    # Cache TTL defaults (seconds)
    DEFAULT_TTL = 30 * 24 * 3600  # 30 days

    def __init__(
        self,
        api_key: Optional[str] = None,
        email: Optional[str] = None,
        cache_dir: Optional[str] = None,
        ttl_seconds: int = DEFAULT_TTL,
        max_retries: int = 3,
    ) -> None:
        self.api_key = api_key or os.environ.get("NCBI_API_KEY") or ""
        self.email = email or os.environ.get("NCBI_EMAIL") or "contact@genehealth.app"
        self.ttl = ttl_seconds
        self.max_retries = max_retries

        Entrez.email = self.email
        if self.api_key:
            Entrez.api_key = self.api_key
            self._min_interval = 1.0 / 10  # 10 req/s with key
        else:
            self._min_interval = 1.0 / 3  # 3 req/s without key

        self.cache_dir = Path(cache_dir or os.environ.get("NCBI_CACHE_DIR", "/tmp/ncbi_cache"))
        self.cache_dir.mkdir(parents=True, exist_ok=True)

        self._lock = threading.Lock()
        self._last_call = 0.0

    # ---------------------------------------------------------------- cache
    def _cache_key(self, fn: str, **params: Any) -> Path:
        raw = json.dumps({"fn": fn, **params}, sort_keys=True, default=str)
        h = hashlib.sha256(raw.encode()).hexdigest()[:24]
        sub = h[:2]
        (self.cache_dir / sub).mkdir(parents=True, exist_ok=True)
        return self.cache_dir / sub / f"{h}.json"

    def _cache_get(self, path: Path) -> Optional[Any]:
        if not path.exists():
            return None
        try:
            age = time.time() - path.stat().st_mtime
            if age > self.ttl:
                return None
            return json.loads(path.read_text())
        except (json.JSONDecodeError, OSError) as e:
            log.warning("cache miss (corrupt): %s (%s)", path, e)
            return None

    def _cache_put(self, path: Path, value: Any) -> None:
        try:
            path.write_text(json.dumps(value, default=str))
        except OSError as e:
            log.warning("cache write failed: %s (%s)", path, e)

    # ---------------------------------------------------------- rate-limit
    def _throttle(self) -> None:
        """Block until at least min_interval has passed since last call."""
        with self._lock:
            elapsed = time.time() - self._last_call
            wait = self._min_interval - elapsed
            if wait > 0:
                time.sleep(wait)
            self._last_call = time.time()

    # ---------------------------------------------------------------- core
    def _call_with_retry(self, fn, **kwargs):
        """Execute an Entrez call with throttle + retry + exponential backoff."""
        last_err: Optional[Exception] = None
        for attempt in range(self.max_retries):
            self._throttle()
            try:
                handle = fn(**kwargs)
                try:
                    data = Entrez.read(handle)
                finally:
                    handle.close()
                return data
            except (HTTPError, URLError, TimeoutError) as e:
                last_err = e
                sleep = (2 ** attempt) + random.uniform(0, 0.5)
                log.warning(
                    "NCBI call failed (attempt %d/%d): %s — retrying in %.1fs",
                    attempt + 1, self.max_retries, e, sleep,
                )
                time.sleep(sleep)
        raise RuntimeError(f"NCBI call exhausted retries: {last_err}") from last_err

    # -------------------------------------------------------------- public
    def esearch(
        self,
        db: str,
        term: str,
        retmax: int = 5,
        sort: str = "pub_date",
        mindate: Optional[str] = None,
        maxdate: Optional[str] = None,
    ) -> List[str]:
        """Search a database. Returns list of IDs (strings)."""
        cache_key = self._cache_key(
            "esearch", db=db, term=term, retmax=retmax, sort=sort,
            mindate=mindate, maxdate=maxdate,
        )
        cached = self._cache_get(cache_key)
        if cached is not None:
            return cached

        params: Dict[str, Any] = {"db": db, "term": term, "retmax": retmax, "sort": sort}
        if mindate:
            params["mindate"] = mindate
        if maxdate:
            params["maxdate"] = maxdate
            params["datetype"] = "pdat"

        data = self._call_with_retry(Entrez.esearch, **params)
        ids: List[str] = list(data.get("IdList", []))
        self._cache_put(cache_key, ids)
        return ids

    def esummary(self, db: str, ids: Iterable[str]) -> List[Dict[str, Any]]:
        """Fetch summary records for a list of IDs."""
        id_list = list(ids)
        if not id_list:
            return []

        cache_key = self._cache_key("esummary", db=db, ids=sorted(id_list))
        cached = self._cache_get(cache_key)
        if cached is not None:
            return cached

        data = self._call_with_retry(Entrez.esummary, db=db, id=",".join(id_list))
        # Entrez esummary returns a list-like; we normalize to list of dicts.
        out: List[Dict[str, Any]] = []
        for rec in data:
            try:
                out.append(dict(rec))
            except (TypeError, ValueError):
                continue
        self._cache_put(cache_key, out)
        return out

    def elink(
        self,
        dbfrom: str,
        db: str,
        ids: Iterable[str],
        linkname: Optional[str] = None,
    ) -> List[str]:
        """Find related records across databases. Returns list of target IDs."""
        id_list = list(ids)
        if not id_list:
            return []

        cache_key = self._cache_key(
            "elink", dbfrom=dbfrom, db=db, ids=sorted(id_list), linkname=linkname,
        )
        cached = self._cache_get(cache_key)
        if cached is not None:
            return cached

        params: Dict[str, Any] = {"dbfrom": dbfrom, "db": db, "id": ",".join(id_list)}
        if linkname:
            params["linkname"] = linkname
        data = self._call_with_retry(Entrez.elink, **params)

        target_ids: List[str] = []
        try:
            for linkset in data:
                for linksetdb in linkset.get("LinkSetDb", []):
                    for link in linksetdb.get("Link", []):
                        target_ids.append(str(link["Id"]))
        except (KeyError, TypeError):
            pass
        self._cache_put(cache_key, target_ids)
        return target_ids


# Module-level default client for convenience
_default: Optional[NCBIClient] = None


def default_client() -> NCBIClient:
    global _default
    if _default is None:
        _default = NCBIClient()
    return _default
