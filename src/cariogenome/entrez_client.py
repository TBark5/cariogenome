"""Polite, cached access to NCBI Entrez and to the RCSB/UniProt web services.

Every download is written to ``data/raw`` and logged with its access date, so a second
run never touches the network.
"""
from __future__ import annotations

import json
import time
import urllib.request
from datetime import date
from pathlib import Path

from Bio import Entrez

from .config import RAW, load_config

LOG_PATH = RAW / "download_log.json"
_last_request = [0.0]


def _log(key: str, info: dict) -> None:
    """Append an entry to the download log (access date, source, file)."""
    RAW.mkdir(parents=True, exist_ok=True)
    log = json.loads(LOG_PATH.read_text()) if LOG_PATH.exists() else {}
    log[key] = info
    LOG_PATH.write_text(json.dumps(log, indent=2))


def download_log() -> dict:
    """Return the download log, or an empty dict if nothing was downloaded yet."""
    return json.loads(LOG_PATH.read_text()) if LOG_PATH.exists() else {}


def _throttle() -> None:
    """Sleep so consecutive NCBI requests are at least ``request_delay_s`` apart."""
    delay = load_config()["ncbi"]["request_delay_s"]
    wait = _last_request[0] + delay - time.time()
    if wait > 0:
        time.sleep(wait)
    _last_request[0] = time.time()


def _with_retries(func, *args, **kwargs):
    """Call ``func`` with exponential backoff on HTTP 429/5xx and network errors."""
    retries = load_config()["ncbi"]["max_retries"]
    for attempt in range(retries):
        try:
            _throttle()
            return func(*args, **kwargs)
        except Exception as exc:  # urllib raises several unrelated exception types
            if attempt == retries - 1:
                raise
            pause = 2.0 * (attempt + 1)
            print(f"    request failed ({exc}); retrying in {pause:.0f} s")
            time.sleep(pause)
    raise RuntimeError("unreachable")


def fetch_genome(accession: str) -> Path:
    """Download a complete GenBank record (with sequence) unless it is already cached."""
    Entrez.email = load_config()["ncbi"]["email"]
    path = RAW / f"{accession}.gb"
    if path.exists() and path.stat().st_size > 1000:
        return path

    def _get() -> str:
        handle = Entrez.efetch(db="nuccore", id=accession, rettype="gbwithparts", retmode="text")
        text = handle.read()
        handle.close()
        if not text.startswith("LOCUS"):
            raise IOError(f"unexpected response for {accession}")
        return text

    text = _with_retries(_get)
    RAW.mkdir(parents=True, exist_ok=True)
    path.write_text(text)
    _log(accession, {"source": "NCBI nuccore efetch gbwithparts",
                     "accessed": date.today().isoformat(), "file": path.name})
    return path


def fetch_url(url: str, filename: str, source: str) -> Path:
    """Download ``url`` to ``data/raw/<filename>`` once, logging the access date."""
    path = RAW / filename
    if path.exists() and path.stat().st_size > 0:
        return path

    def _get() -> bytes:
        with urllib.request.urlopen(url, timeout=60) as resp:
            return resp.read()

    content = _with_retries(_get)
    RAW.mkdir(parents=True, exist_ok=True)
    path.write_bytes(content)
    _log(filename, {"source": source, "url": url,
                    "accessed": date.today().isoformat(), "file": filename})
    return path


def ncbi_reachable(timeout: float = 10.0) -> bool:
    """True if the NCBI E-utilities endpoint answers."""
    try:
        url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/einfo.fcgi"
        with urllib.request.urlopen(url, timeout=timeout) as resp:
            return resp.status == 200
    except Exception:
        return False
