"""Documentation, figure, provenance and privacy checks.

These tests keep the written claims tied to the files: numbers in README.md and
MORNING_REPORT.md must exist in results/, every figure must be captioned and referenced,
links must resolve, and no tracked file may leak a local path.
"""

from __future__ import annotations

import re
import subprocess
from pathlib import Path

import pandas as pd
import pytest
import yaml

from cariogenome import entrez_client
from cariogenome.captions import CAPTIONS, expected_files

ROOT = Path(__file__).resolve().parents[1]
RES = ROOT / "results"
FIG = ROOT / "figures"
needs_results = pytest.mark.skipif(
    not (RES / "summary_effect_sizes.csv").exists(), reason="run run_all.py first"
)

# Numbers that legitimately appear in the documents without being computed results:
# design constants, thresholds, citation metadata, versions and screenshot geometry.
DOC_CONSTANTS = {
    "0",
    "1",
    "2",
    "3",
    "4",
    "5",
    "6",
    "7",
    "8",
    "9",
    "10",
    "11",
    "12",
    "13",
    "14",
    "15",
    "16",
    "20",
    "22",
    "28",
    "30",
    "36",
    "45",
    "57",
    "60",
    "70",
    "100",
    "300",
    "1000",
    "1250",
    "1500",
    "0.05",
    "0.60",
    "0.75",
    "3.10",
    "3.14",
    "2.5.5",
    "2.1",
    "16S",
    "2011",
    "2012",
    "2013",
    "1997",
    "2006",
    "2026",
    "408",
    "177",
    "518",
    "881",
    "383",
    "226",
    "239",
    "44",
    "21354427",
    "22816041",
    "23228887",
    "9089078",
    "16239014",
}


def _markdown_numbers(text: str) -> list[str]:
    """Numeric tokens outside code blocks and URLs, normalized (typographic minus -> '-')."""
    text = re.sub(r"```.*?```", "", text, flags=re.S)
    text = re.sub(r"\(https?://[^)]*\)|https?://\S+", "", text)
    text = re.sub(r"\]\([^)]*\)", "]", text)  # link targets
    text = re.sub(r"`[^`]*`", "", text)  # inline code (file names, flags)
    text = re.sub(r"\d{4}-\d{2}-\d{2}", "", text)  # ISO dates
    text = text.replace("−", "-")  # noqa: RUF001 (typographic minus in the docs)
    tokens = re.findall(r"(?<![\w.])-?\d+(?:\.\d+)*%?", text)
    return [t.rstrip("%").lstrip("-") for t in tokens]


def _result_numbers() -> set[str]:
    """Every value in results/ tables, formatted the ways the documents may print it."""
    out: set[str] = set()
    for csv in RES.rglob("*.csv"):
        df = pd.read_csv(csv)
        vals = pd.to_numeric(df.stack(), errors="coerce").dropna().abs()
        for v in vals.unique():
            for d in range(5):
                out.add(f"{v:.{d}f}")
                out.add(f"{100 * v:.{d}f}")
    tests = sum(
        p.read_text(encoding="utf-8").count("\ndef test_")
        for p in (ROOT / "tests").glob("test_*.py")
    )
    out.add(str(tests))
    return out


@needs_results
@pytest.mark.parametrize("doc", ["README.md", "MORNING_REPORT.md"])
def test_every_number_in_docs_is_traceable(doc: str) -> None:
    allowed = _result_numbers() | DOC_CONSTANTS
    nums = _markdown_numbers((ROOT / doc).read_text(encoding="utf-8"))
    unknown = sorted({n for n in nums if n not in allowed})
    assert not unknown, f"{doc}: numbers not found in results/ or constants: {unknown}"
    print(f"{doc}: verified {len(nums)} numeric values ({len(set(nums))} distinct)")


def _tracked_markdown() -> list[Path]:
    files = subprocess.run(
        ["git", "ls-files", "*.md"], cwd=ROOT, capture_output=True, text=True, check=True
    ).stdout.split()
    return [ROOT / f for f in files if (ROOT / f).exists()] + [ROOT / "figures" / "CAPTIONS.md"]


def test_captions_one_per_figure_file() -> None:
    lines = (FIG / "CAPTIONS.md").read_text(encoding="utf-8").splitlines()
    captioned = [
        f for line in lines for f in re.findall(r"\*\*(\d\d_[\w]+\.(?:png|html))\*\*", line)
    ]
    assert len(captioned) == len(set(captioned)), "a figure file has more than one caption"
    on_disk = sorted(p.name for p in FIG.iterdir() if p.suffix in (".png", ".html"))
    assert sorted(captioned) == on_disk == sorted(expected_files())
    assert len(CAPTIONS) == len({f.split(".")[0] for f in on_disk})


def test_every_figure_referenced_and_every_reference_exists() -> None:
    texts = {p: p.read_text(encoding="utf-8") for p in _tracked_markdown()}
    texts[ROOT / "app.py"] = (ROOT / "app.py").read_text(encoding="utf-8")
    refs = set()
    for text in texts.values():
        refs |= set(re.findall(r"(\d\d_[A-Za-z0-9_]+\.(?:png|html))", text))
    for r in refs:
        assert (FIG / r).exists(), f"referenced figure missing: {r}"
    app_names = set(
        re.findall(
            r'"(m\d_[A-Za-z0-9_]+|summary_[a-z_]+|validation_[a-z_]+)"', texts[ROOT / "app.py"]
        )
    )
    for f in FIG.iterdir():
        if f.suffix in (".png", ".html"):
            key = f.stem[3:]
            assert f.name in refs or key in app_names, f"figure never referenced: {f.name}"


def _slug(heading: str) -> str:
    s = heading.strip().lower()
    s = re.sub(r"[^\w\- ]", "", s)
    return s.replace(" ", "-")


def test_internal_markdown_links_resolve() -> None:
    for md in _tracked_markdown():
        text = md.read_text(encoding="utf-8")
        text_no_code = re.sub(r"```.*?```", "", text, flags=re.S)
        anchors = {_slug(h) for h in re.findall(r"^#+\s+(.*)$", text, flags=re.M)}
        for target in re.findall(r"\]\(([^)\s]+)\)", text_no_code):
            if target.startswith(("http://", "https://", "mailto:")):
                continue
            path, _, anchor = target.partition("#")
            if path:
                assert (md.parent / path).exists(), f"{md.name}: broken link {target}"
            elif anchor:
                assert anchor in anchors, f"{md.name}: broken anchor #{anchor}"


def test_no_absolute_paths_or_home_directories_in_tracked_files() -> None:
    files = subprocess.run(
        ["git", "ls-files"], cwd=ROOT, capture_output=True, text=True, check=True
    ).stdout.split()
    pattern = re.compile(
        r"[A-Za-z]:[\\/]+(Users|Documents and Settings)[\\/]|/home/\w|/Users/\w"
        r"|AppData[\\/]|(?<![\w/(])~/[A-Za-z]"
    )  # not a regex like /~/g
    offenders = []
    for f in files:
        p = ROOT / f
        if not p.exists() or p.suffix in (".png", ".gz", ".pkl") or p.stat().st_size > 5e6:
            continue
        text = p.read_text(encoding="utf-8", errors="ignore")
        text = re.sub(r"https?://\S+", "", text)  # URLs are not local paths
        if pattern.search(text):
            offenders.append(f)
    assert not offenders, f"absolute/home paths in: {offenders}"


def test_config_keeps_the_email_placeholder() -> None:
    cfg = yaml.safe_load((ROOT / "config.yaml").read_text(encoding="utf-8"))
    assert cfg["ncbi"]["email"] == entrez_client.EMAIL_PLACEHOLDER


def test_placeholder_email_raises_readable_error(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.delenv("NCBI_EMAIL", raising=False)
    with pytest.raises(ValueError, match="NCBI_EMAIL"):
        entrez_client.contact_email()
    monkeypatch.setenv("NCBI_EMAIL", "someone@example.org")
    assert entrez_client.contact_email() == "someone@example.org"


def test_offline_mode_blocks_network(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setenv("CARIOGENOME_OFFLINE", "1")
    assert entrez_client.ncbi_reachable() is False
    with pytest.raises(entrez_client.OfflineError):
        entrez_client.fetch_genome("NC_000000.0")
    with pytest.raises(entrez_client.OfflineError):
        entrez_client.fetch_url("https://example.org/x", "never_cached.bin", "test")


@needs_results
def test_accessions_well_formed_and_consistent() -> None:
    text = (ROOT / "ACCESSIONS.md").read_text(encoding="utf-8")
    genome_rows = re.findall(r"^\| (\w+_[\w-]+) \| \*S\. \w+\* \| (\S+) \|", text, flags=re.M)
    assert len(genome_rows) == 22
    for _, acc in genome_rows:
        assert re.fullmatch(r"(NC|NZ)_[A-Z]{0,2}\d+\.\d+", acc), acc
    cat = pd.read_csv(RES / "m1_catalog.csv")
    rec = re.findall(
        r"^\| (\w+) \| (\w+_[\w-]+) \| \S+ \| [^|]* \| ((?:NC|NZ)_\S+?):(\d+)-(\d+)"
        r"\(([+-])\) \| [\d.]+ \| (yes|no) \|",
        text,
        flags=re.M,
    )
    assert len(rec) == len(cat) == 286
    assert sum(r[6] == "yes" for r in rec) == int(cat["included"].sum()) == 284
    assert {r[2] for r in rec} <= {acc for _, acc in genome_rows}


def test_readme_scorecard_matches_results_discussion() -> None:
    def scorecard(path: Path) -> list[str]:
        lines = path.read_text(encoding="utf-8").splitlines()
        return [ln for ln in lines if re.match(r"^\| P\d", ln)]

    readme, disc = scorecard(ROOT / "README.md"), scorecard(ROOT / "RESULTS_DISCUSSION.md")
    assert readme and readme == disc


def test_markdown_tables_have_consistent_columns() -> None:
    """GitHub silently drops or misrenders table rows whose column count differs."""
    for md in _tracked_markdown():
        block: list[str] = []
        for line in [*md.read_text(encoding="utf-8").splitlines(), ""]:
            if line.startswith("|"):
                block.append(line)
            elif block:
                counts = {row.replace(r"\|", "").count("|") for row in block}  # \| is escaped
                assert len(counts) == 1, f"{md.name}: ragged table starting {block[0][:50]!r}"
                block = []


def test_figures_are_opaque_for_dark_mode() -> None:
    from PIL import Image

    for png in [*FIG.glob("*.png"), *(ROOT / "docs" / "screenshots").glob("*.png")]:
        with Image.open(png) as im:
            alpha = im.convert("RGBA").getchannel("A")
            assert alpha.getextrema()[0] == 255, f"{png.name} has transparent pixels"
