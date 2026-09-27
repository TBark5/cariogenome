"""Polite, cached access to NCBI Entrez and to the RCSB/UniProt web services.

Every download is written to ``data/raw`` and logged with its access date, so a second
run never touches the network. Setting ``CARIOGENOME_OFFLINE=1`` (``run_all.py --offline``)
forbids all network access: cached files are still used, anything else raises
``OfflineError``.
"""
from __future__ import annotations

import json
import os
import time
import urllib.request
from collections.abc import Callable
from datetime import date
from pathlib import Path
from typing import TypeVar

from Bio import Entrez

from .config import RAW, load_config

LOG_PATH = RAW / "download_log.json"
EMAIL_PLACEHOLDER = "REPLACE_WITH_YOUR_EMAIL"
_last_request = [0.0]
T = TypeVar("T")


class OfflineError(RuntimeError):
    """Raised when a download is needed but offline mode is on."""


def offline() -> bool:
    """True if network access is disabled (``CARIOGENOME_OFFLINE=1``)."""
    return os.environ.get("CARIOGENOME_OFFLINE", "") == "1"


def contact_email() -> str:
    """NCBI contact email: ``NCBI_EMAIL`` environment variable, else ``config.yaml``.

    NCBI requires every E-utilities user to identify themselves. The committed config
    holds a placeholder, so a real address must be supplied before downloading.
    """
    email = os.environ.get("NCBI_EMAIL") or load_config()["ncbi"]["email"]
    if not email or email == EMAIL_PLACEHOLDER or "@" not in email:
        raise ValueError(
            "No NCBI contact email is set. NCBI requires one for downloads. Set the "
            "environment variable NCBI_EMAIL (for example `set NCBI_EMAIL=you@example.org` "
            "on Windows or `export NCBI_EMAIL=you@example.org` on macOS/Linux), or replace "
            f"{EMAIL_PLACEHOLDER!r} in config.yaml. Alternatively run `python run_all.py "
            "--offline` to use the committed sequence cache without downloading."
        )
    return email


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


def _with_retries(func: Callable[[], T]) -> T:
    """Call ``func`` with exponential backoff on HTTP 429/5xx and network errors."""
    if offline():
        raise OfflineError("offline mode: network access is disabled")
    retries = load_config()["ncbi"]["max_retries"]
    for attempt in range(retries):
        try:
            _throttle()
            return func()
        # urllib and http.client raise several unrelated types (URLError, HTTPError,
        # IncompleteRead, timeouts); every failure is reported and retried, then re-raised.
        except Exception as exc:
            if attempt == retries - 1:
                raise
            pause = 2.0 * (attempt + 1)
            print(f"    request failed ({exc}); retrying in {pause:.0f} s")
            time.sleep(pause)
    raise RuntimeError("unreachable")


def fetch_genome(accession: str) -> Path:
    """Download a complete GenBank record (with sequence) unless it is already cached."""
    path = RAW / f"{accession}.gb"
    if path.exists() and path.stat().st_size > 1000:
        return path
    if offline():
        raise OfflineError(f"offline mode: {accession} is not in the local cache")
    Entrez.email = contact_email()

    def _get() -> str:
        handle = Entrez.efetch(db="nuccore", id=accession, rettype="gbwithparts", retmode="text")
        text = handle.read()
        handle.close()
        if not text.startswith("LOCUS"):
            raise OSError(f"unexpected response for {accession}")
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
    """True if the NCBI E-utilities endpoint answers (always False in offline mode)."""
    if offline():
        return False
    try:
        url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/einfo.fcgi"
        with urllib.request.urlopen(url, timeout=timeout) as resp:
            return bool(resp.status == 200)
    except OSError:  # URLError, HTTPError and socket timeouts are all OSError subclasses
        return False
