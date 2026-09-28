"""Shared building blocks of the Streamlit dashboard (``app.py`` and ``app_pages/``).

Everything shown is read from ``results/`` and ``figures/`` (produced by ``run_all.py``);
nothing is recomputed on request. Interactive charts use the fixed color semantics of
the static figures: vermillion = virulence-associated, blue = housekeeping control.
"""

from __future__ import annotations

import io
from collections.abc import Sequence
from dataclasses import dataclass
from functools import lru_cache
from typing import Literal
from pathlib import Path

import altair as alt
import numpy as np
import pandas as pd
import streamlit as st
from matplotlib.figure import Figure
from PIL import Image

from .captions import CAPTIONS
from .config import RESULTS, all_genes, gene_class, genomes, load_config
from .m4_figures import draw_tree
from .m4_phylogeny import load_tree
from .plotting import CLASS_COLORS, GREY, W_SINGLE, figure_path

REPO_URL = "https://github.com/TBark5/cariogenome"
ALPHA = 0.05  # BH q-value threshold used throughout the project
CLASS_NAMES = {"virulence": "Virulence-associated", "control": "Housekeeping"}
BadgeColor = Literal["red", "orange", "yellow", "blue", "green", "violet", "gray", "primary"]
CLASS_BADGE_COLORS: dict[str, BadgeColor] = {"virulence": "orange", "control": "blue"}
CSS_PX_PER_INCH = 100  # static figures are 300 dpi; show them at a third of their pixel size

METRIC_LABELS = {
    "gc": "GC content",
    "gc3": "GC3",
    "cai": "Codon adaptation index",
    "gc_skew": "GC skew",
    "mean_pid_Smutans": "Mean aa identity, S. mutans (%)",
    "mean_entropy_Smutans": "Mean per-site entropy, S. mutans",
    "nRF_within_Smutans": "RF distance to reference, S. mutans",
    "n_supported_conflicts": "Supported conflicting splits",
    "n_supported_conflicts_between_species": "Supported species-mixing splits",
    "mean_support": "Mean bootstrap support (%)",
    "omega": "dN/dS (ω)",
    "dN": "dN",
    "dS": "dS",
}


# ---------------------------------------------------------------------------- data


@st.cache_data(max_entries=200, show_spinner=False)
def csv(name: str, index_col: int | None = None) -> pd.DataFrame:
    """Read one results table (cached for every visitor)."""
    return pd.read_csv(RESULTS / name, index_col=index_col)


def results_ready() -> bool:
    return (RESULTS / "m1_catalog.csv").exists()


@st.cache_data(show_spinner=False)
def data_mode() -> str:
    """``REAL`` for NCBI data, ``SYNTHETIC`` for the simulated fallback."""
    return (RESULTS / "data_mode.txt").read_text().strip()


@lru_cache(maxsize=1)
def locus_tags() -> dict[str, str]:
    """UA159 locus tag (or description) of every panel locus, from config.yaml."""
    tags: dict[str, str] = {}
    for group in load_config()["genes"].values():
        tags.update({g: str(t) for g, t in group.items()})
    return tags


def species_of_label() -> dict[str, str]:
    return {g.label: g.species for g in genomes()}


def is_dark() -> bool:
    return st.context.theme.type == "dark"


def accent() -> str:
    """The theme's primary color (for highlights that must work in both modes)."""
    return "#2DD4BF" if is_dark() else "#0F766E"


def muted() -> str:
    return "#2A3B41" if is_dark() else "#E4EAEC"


# ------------------------------------------------------------------ page furniture


@dataclass(frozen=True)
class Prediction:
    code: str
    claim: str
    verdict: str
    color: BadgeColor
    icon: str


_SUPPORTED: tuple[BadgeColor, str] = ("green", ":material/check_circle:")
_PARTLY: tuple[BadgeColor, str] = ("orange", ":material/incomplete_circle:")
_NOT: tuple[BadgeColor, str] = ("gray", ":material/do_not_disturb_on:")
_FALSIFIED: tuple[BadgeColor, str] = ("red", ":material/cancel:")

# Pre-registered predictions (HYPOTHESIS.md) and their verdicts (RESULTS_DISCUSSION.md).
PREDICTIONS = {
    p.code: p
    for p in [
        Prediction("P1", "Virulence genes have higher dN/dS", "Supported", *_SUPPORTED),
        Prediction("P2", "Virulence proteins are less conserved", "Supported", *_SUPPORTED),
        Prediction("P3", "Their compositional signature differs", "Not supported", *_NOT),
        Prediction(
            "P4a",
            "gtf, ftf and spaP are restricted to S. mutans; luxS and controls are universal",
            "Partly supported",
            *_PARTLY,
        ),
        Prediction("P4b", "Virulence gene trees are more discordant", "Not supported", *_NOT),
        Prediction(
            "P5",
            "Catalytic residues are in the top 10% most conserved columns",
            "Narrowly falsified",
            *_FALSIFIED,
        ),
        Prediction(
            "P6", "Conserved residues cluster near the active site", "Supported", *_SUPPORTED
        ),
    ]
}


def prediction_badge(code: str) -> None:
    p = PREDICTIONS[code]
    st.badge(f"{p.code} · {p.verdict.lower()}", icon=p.icon, color=p.color, help=p.claim)


def gene_picker() -> str:
    """The shared gene selector; its value follows the visitor across pages and the URL."""
    gene: str = st.selectbox(
        "Gene",
        all_genes(),
        key="gene",
        bind="query-params",
        persist_state="session",
        help="Drives the gene-level views on the M1 to M5 pages. The choice is kept in the "
        "URL, so a link opens on the same gene.",
    )
    return gene


def page_header(
    title: str,
    icon: str,
    question: str,
    predictions: Sequence[str] = (),
    *,
    gene_select: bool = False,
) -> str | None:
    """Title, research question, verdict badges and (optionally) the gene selector."""
    gene = None
    if gene_select:
        left, right = st.columns([3, 1], vertical_alignment="bottom")
        left.title(title, icon=icon)
        with right:
            gene = gene_picker()
    else:
        st.title(title, icon=icon)
    st.markdown(question)
    if predictions or gene:
        with st.container(horizontal=True, gap="small"):
            for code in predictions:
                prediction_badge(code)
            if gene:
                cls = gene_class(gene)
                st.badge(
                    "rRNA control" if gene == "16S" else CLASS_NAMES[cls],
                    icon=":material/genetics:",
                    color=CLASS_BADGE_COLORS[cls],
                )
                st.badge(f"UA159 {locus_tags()[gene]}", color="gray")
    if data_mode() != "REAL":
        st.warning(
            "These results were produced from SYNTHETIC simulated sequences, not NCBI data.",
            icon=":material/warning:",
        )
    return gene


def section(title: str, icon: str, caption: str | None = None) -> None:
    st.subheader(title, icon=icon)
    if caption:
        st.caption(caption)


# -------------------------------------------------------------------------- tables


def _fmt(x: float) -> str:
    return f"{x:.3g}" if abs(x) < 1000 else f"{x:,.0f}"


def tests_table(df: pd.DataFrame) -> None:
    """Virulence-vs-housekeeping comparisons with effect sizes, CIs and BH q-values."""
    lo = "delta_ci_low" if "delta_ci_low" in df else "ci_low"
    hi = "delta_ci_high" if "delta_ci_high" in df else "ci_high"
    names = df["label"] if "label" in df else df["metric"].map(lambda m: METRIC_LABELS.get(m, m))
    out = pd.DataFrame(
        {
            "Metric": names,
            "Virulence (median)": df["median_virulence"].map(_fmt),
            "Housekeeping (median)": df["median_control"].map(_fmt),
            "Cliff's δ": df["cliffs_delta"],
            "95% CI": [f"[{a:.2f}, {b:.2f}]" for a, b in zip(df[lo], df[hi], strict=True)],
            "q (BH)": df["q_bh"],
            "q < 0.05": df["q_bh"] < ALPHA,
        }
    )
    if "module" in df:
        out.insert(0, "Module", df["module"])
    st.dataframe(
        out,
        hide_index=True,
        column_config={
            "Cliff's δ": st.column_config.NumberColumn(
                format="%+.2f",
                help="Cliff's delta, virulence vs housekeeping: +1 = every virulence gene is "
                "higher, -1 = every virulence gene is lower.",
            ),
            "q (BH)": st.column_config.NumberColumn(
                format="%.4f", help="Benjamini-Hochberg adjusted p-value (Mann-Whitney U)."
            ),
            "q < 0.05": st.column_config.CheckboxColumn(),
        },
    )


# ------------------------------------------------------------------------- figures


@st.cache_data(max_entries=64, show_spinner=False)
def _pixel_width(path: str) -> int:
    with Image.open(path) as im:
        return int(im.width)


def figure(name: str) -> None:
    """A registered static figure (figures/NN_<name>.png) with its caption."""
    path = figure_path(name)
    if not path.exists():
        st.caption(f"Figure {path.name} not found; run `python run_all.py`.")
        return
    width = max(320, _pixel_width(str(path)) * CSS_PX_PER_INCH // 300)
    st.image(str(path), width=width)
    st.caption(f"**{path.name}** · {CAPTIONS[name]}")


def figure_gallery(figures: Sequence[tuple[str, str]], key: str) -> None:
    """Publication figures in lazy tabs: only the open tab is sent to the browser."""
    section(
        "Publication figures",
        ":material/photo_library:",
        "The 300 dpi figures written by run_all.py. Right-click to save.",
    )
    with st.container(border=True):
        tabs = st.tabs([label for _, label in figures], on_change="rerun", key=key)
        for tab, (name, _) in zip(tabs, figures, strict=True):
            if tab.open:
                with tab:
                    figure(name)


@st.cache_data(max_entries=64, show_spinner=False)
def tree_png(gene: str, method: str) -> bytes:
    """A gene tree drawn with the static-figure style (cached; thread-safe Figure API)."""
    tree = load_tree(f"{gene}_{method}")
    n_taxa = len(tree.get_terminals())
    fig = Figure(figsize=(W_SINGLE, 0.32 * n_taxa + 1.5))
    ax = fig.subplots()
    draw_tree(ax, tree, f"{gene} {method.upper()} tree")
    fig.tight_layout()
    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=150, facecolor="white")
    return buf.getvalue()


@st.cache_data(max_entries=4, show_spinner=False)
def viewer_html(path: str) -> str:
    """The self-contained 3D viewer on a white card, so it reads the same in dark mode."""
    html = Path(path).read_text(encoding="utf-8")
    style = (
        "<style>html,body{background:#fff;color:#15232B;margin:0;"
        "font-family:Inter,'Segoe UI',sans-serif}body{padding:12px 16px}</style>"
    )
    return html.replace("<body>", style + "<body>", 1)


# ------------------------------------------------------------- interactive charts


def _class_scale() -> alt.Scale:
    return alt.Scale(
        domain=list(CLASS_NAMES.values()),
        range=[CLASS_COLORS["virulence"], CLASS_COLORS["control"]],
    )


def _with_class(df: pd.DataFrame, gene_col: str = "gene") -> pd.DataFrame:
    return df.assign(Class=df[gene_col].map(lambda g: CLASS_NAMES[gene_class(g)]))


def forest_chart(eff: pd.DataFrame) -> alt.LayerChart:
    """Cliff's delta with 95% CI for every virulence-vs-housekeeping test."""
    df = eff.assign(
        test=eff["module"].str.split().str[0]
        + " · "
        + eff["metric"].map(METRIC_LABELS).fillna(eff["label"]),
        result=np.where(eff["q_bh"] < ALPHA, "q < 0.05", "not significant"),
    )
    color = alt.Color(
        "result:N",
        scale=alt.Scale(domain=["q < 0.05", "not significant"], range=[accent(), GREY]),
        legend=alt.Legend(title=None, orient="top"),
    )
    y = alt.Y("test:N", sort=None, title=None, axis=alt.Axis(labelLimit=260))
    x_scale = alt.Scale(domain=[-1, 1])
    tooltip = [
        alt.Tooltip("test:N", title="Test"),
        alt.Tooltip("cliffs_delta:Q", title="Cliff's δ", format=".2f"),
        alt.Tooltip("ci_low:Q", title="CI low", format=".2f"),
        alt.Tooltip("ci_high:Q", title="CI high", format=".2f"),
        alt.Tooltip("q_bh:Q", title="q (BH)", format=".4f"),
    ]
    zero = alt.Chart(pd.DataFrame({"x": [0]})).mark_rule(strokeDash=[4, 4], color=GREY)
    zero = zero.encode(x=alt.X("x:Q", scale=x_scale))
    base = alt.Chart(df).encode(y=y, color=color, tooltip=tooltip)
    ci = base.mark_rule(strokeWidth=3, opacity=0.55).encode(
        x=alt.X("ci_low:Q", scale=x_scale, title="Cliff's δ  (← lower in virulence · higher →)"),
        x2="ci_high:Q",
    )
    dot = base.mark_circle(size=130, opacity=1).encode(x="cliffs_delta:Q")
    return (zero + ci + dot).properties(height=alt.Step(30))


def gene_interval_chart(
    df: pd.DataFrame,
    value: str,
    low: str,
    high: str,
    selected: str,
    title: str,
    fmt: str = ".3f",
) -> alt.LayerChart:
    """One point + CI per gene, colored by class, with the selected gene highlighted."""
    df = _with_class(df.reset_index()).assign(selected=lambda d: d["gene"] == selected)
    order = [g for g in all_genes() if g in set(df["gene"])]
    x = alt.X("gene:N", sort=order, title=None, axis=alt.Axis(labelAngle=0))
    opacity = alt.condition("datum.selected", alt.value(1.0), alt.value(0.35))
    tooltip = [
        alt.Tooltip("gene:N", title="Gene"),
        alt.Tooltip("Class:N"),
        alt.Tooltip(f"{value}:Q", title=title, format=fmt),
        alt.Tooltip(f"{low}:Q", title="CI low", format=fmt),
        alt.Tooltip(f"{high}:Q", title="CI high", format=fmt),
    ]
    base = alt.Chart(df).encode(
        x=x,
        color=alt.Color("Class:N", scale=_class_scale(), legend=alt.Legend(orient="top")),
        opacity=opacity,
        tooltip=tooltip,
    )
    bars = base.mark_rule(strokeWidth=3).encode(y=alt.Y(f"{low}:Q", title=title), y2=f"{high}:Q")
    dots = base.mark_circle(size=140).encode(y=f"{value}:Q")
    ring = (
        alt.Chart(df[df["selected"]])
        .mark_circle(size=420, filled=False, strokeWidth=2, color=accent())
        .encode(x=x, y=f"{value}:Q")
    )
    return (bars + dots + ring).properties(height=320)


def strain_strip_chart(per: pd.DataFrame, metric: str, selected: str) -> alt.LayerChart:
    """Per-strain values of one composition metric in S. mutans, one column per gene."""
    df = _with_class(per[per["species"] == "S. mutans"].dropna(subset=[metric])).assign(
        selected=lambda d: d["gene"] == selected
    )
    order = [g for g in all_genes() if g in set(df["gene"])]
    title = METRIC_LABELS[metric]
    x = alt.X("gene:N", sort=order, title=None, axis=alt.Axis(labelAngle=0))
    base = alt.Chart(df).encode(
        x=x,
        color=alt.Color("Class:N", scale=_class_scale(), legend=alt.Legend(orient="top")),
        opacity=alt.condition("datum.selected", alt.value(1.0), alt.value(0.4)),
    )
    dots = (
        base.transform_calculate(jitter="(random() - 0.5) * 0.5")
        .mark_circle(size=60)
        .encode(
            y=alt.Y(f"{metric}:Q", title=title, scale=alt.Scale(zero=False)),
            xOffset=alt.XOffset("jitter:Q", scale=alt.Scale(domain=[-1, 1])),
            tooltip=["gene:N", "label:N", alt.Tooltip(f"{metric}:Q", title=title, format=".4f")],
        )
    )
    median = base.mark_tick(size=26, thickness=3).encode(y=f"median({metric}):Q")
    return (dots + median).properties(height=330)


def conservation_chart(cons: pd.DataFrame, gene: str, window: int = 15) -> alt.LayerChart:
    """Per-column conservation along the alignment, with a rolling mean."""
    df = cons[["column", "conservation", "ref_residue"]].copy()
    df["smoothed"] = df["conservation"].rolling(window, center=True, min_periods=1).mean()
    color = CLASS_COLORS[gene_class(gene)]
    x = alt.X(
        "column:Q",
        title="alignment column",
        scale=alt.Scale(domain=[1, int(df["column"].max())], nice=False),
    )
    y_scale = alt.Scale(zero=False, domain=[max(0.0, float(df["conservation"].min())), 1.0])
    raw = (
        alt.Chart(df)
        .mark_area(opacity=0.18, color=color, line=False, clip=True)
        .encode(x=x, y=alt.Y("conservation:Q", title="conservation (1 − H / log₂K)", scale=y_scale))
    )
    smooth = (
        alt.Chart(df).mark_line(strokeWidth=2, color=color, clip=True).encode(x=x, y="smoothed:Q")
    )
    hover = alt.selection_point(fields=["column"], nearest=True, on="pointerover", empty=False)
    points = (
        alt.Chart(df)
        .mark_rule(color=GREY)
        .encode(
            x=x,
            opacity=alt.condition(hover, alt.value(0.6), alt.value(0)),
            tooltip=[
                alt.Tooltip("column:Q", title="Column"),
                alt.Tooltip("ref_residue:N", title="UA159 residue"),
                alt.Tooltip("conservation:Q", format=".3f"),
                alt.Tooltip("smoothed:Q", title=f"{window}-column mean", format=".3f"),
            ],
        )
        .add_params(hover)
    )
    return (raw + smooth + points).properties(height=300)


def identity_heatmap(ident: pd.DataFrame) -> alt.Chart:
    """Pairwise percent identity between strains."""
    order = list(ident.index)
    df = ident.rename_axis("a").reset_index().melt(id_vars="a", var_name="b", value_name="pid")
    lo = float(df["pid"].min())
    return (
        alt.Chart(df)
        .mark_rect(stroke=None)
        .encode(
            x=alt.X(
                "b:N", sort=order, title=None, axis=alt.Axis(labelAngle=-45, labelOverlap=False)
            ),
            y=alt.Y("a:N", sort=order, title=None),
            color=alt.Color(
                "pid:Q",
                title="% identity",
                scale=alt.Scale(scheme="viridis", domain=[lo, 100]),
            ),
            tooltip=[
                alt.Tooltip("a:N", title="Strain"),
                alt.Tooltip("b:N", title="vs"),
                alt.Tooltip("pid:Q", title="% identity", format=".2f"),
            ],
        )
        .properties(height=alt.Step(22 if len(order) > 12 else 34))
    )


def presence_heatmap(mat: pd.DataFrame, selected: str) -> alt.Chart:
    """Gene x genome ortholog calls: included, excluded by QC, or absent."""
    status = {2: "Ortholog, passes QC", 1: "Ortholog, excluded by QC", 0: "No ortholog"}
    species = species_of_label()
    df = (
        mat.rename_axis("gene")
        .reset_index()
        .melt(id_vars="gene", var_name="genome", value_name="code")
        .assign(
            status=lambda d: d["code"].map(status),
            species=lambda d: d["genome"].map(species),
            selected=lambda d: d["gene"] == selected,
        )
    )
    return (
        alt.Chart(df)
        .mark_rect(stroke="white" if not is_dark() else "#0B1417", strokeWidth=1.5)
        .encode(
            x=alt.X(
                "genome:N",
                sort=list(mat.columns),
                title=None,
                axis=alt.Axis(labelAngle=-50, labelOverlap=False),
            ),
            y=alt.Y("gene:N", sort=list(mat.index), title=None),
            color=alt.Color(
                "status:N",
                scale=alt.Scale(domain=list(status.values()), range=[accent(), "#E69F00", muted()]),
                legend=alt.Legend(title=None, orient="top", labelLimit=0),
            ),
            opacity=alt.condition("datum.selected", alt.value(1.0), alt.value(0.55)),
            tooltip=["gene:N", "genome:N", "species:N", "status:N"],
        )
        .properties(height=alt.Step(24))
    )


def window_chart(win: pd.DataFrame, gene: str) -> alt.Chart:
    """pN and pS in sliding windows along the gene."""
    df = win[["mid_codon", "pN", "pS"]].melt("mid_codon", var_name="rate", value_name="p")
    return (
        alt.Chart(df)
        .mark_line(strokeWidth=2, interpolate="monotone")
        .encode(
            x=alt.X(
                "mid_codon:Q",
                title="codon position (window midpoint)",
                scale=alt.Scale(domain=[0, float(win["end_codon"].max())], nice=False),
            ),
            y=alt.Y("p:Q", title="p-distance"),
            color=alt.Color(
                "rate:N",
                scale=alt.Scale(domain=["pN", "pS"], range=[CLASS_COLORS[gene_class(gene)], GREY]),
                legend=alt.Legend(title=None, orient="top"),
            ),
            tooltip=[
                alt.Tooltip("mid_codon:Q", title="Codon"),
                "rate:N",
                alt.Tooltip("p:Q", format=".4f"),
            ],
        )
        .properties(height=280)
    )


def asset(name: str) -> Path:
    return Path(__file__).resolve().parents[2] / "assets" / name
