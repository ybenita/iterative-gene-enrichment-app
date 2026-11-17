"""Interactive drawing of the iGEA network that is tested and sent to the LLM.

The figure is built from the same objects as the statistics and the LLM export: the terms that
passed the q-value filter and the communities, hub genes and bridges computed on them
(``network_modules.analyze_network_modules``). Nothing is re-derived from a DOT file.

- Terms are squares, genes circles (size = degree, i.e. number of terms/collections), bridge genes
  (terms in more than one community) diamonds.
- Colour by community (default) or by collection. Each community is drawn as a shaded area; the
  layout pulls the terms and genes of a community together.
"""
from __future__ import annotations

import math
from typing import Dict, List, Optional

import networkx as nx
import plotly.graph_objects as go

# Categorical palette in fixed order (C1, C2, ...); communities beyond 8 are grey ("other").
COMMUNITY_COLOURS = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#6250d6", "#e34948"]
OTHER = "#9a9a94"
INK = "#3d3d3a"
EDGE = "rgba(120,120,120,0.22)"
MAX_EDGES = 3000


def community_colour(label: str) -> str:
    try:
        k = int(str(label)[1:]) - 1
    except ValueError:
        return OTHER
    return COMMUNITY_COLOURS[k] if 0 <= k < len(COMMUNITY_COLOURS) else OTHER


def _rgba(hex_colour: str, alpha: float) -> str:
    h = hex_colour.lstrip("#")
    return f"rgba({int(h[0:2], 16)},{int(h[2:4], 16)},{int(h[4:6], 16)},{alpha})"


def _short(s: str, n: int = 42) -> str:
    return s if len(s) <= n else s[: n - 1] + "…"


def n_edges(terms: List[Dict]) -> int:
    return sum(len(set(t["genes"])) for t in terms)


def _build(terms, modules):
    term_comm, gene_info = {}, {}
    if modules is not None and len(modules.terms):
        term_comm = {(r.Library, r.Term): r.Community for r in modules.terms.itertuples()}
        gene_info = {r.Gene: r for r in modules.genes.itertuples()}
    deg: Dict[str, int] = {}
    for t in terms:
        for g in set(t["genes"]):
            deg[g] = deg.get(g, 0) + 1
    return term_comm, gene_info, deg


def _layout(terms, term_comm, gene_info, seed: int = 2026):
    """Force-directed layout; edges inside a community pull harder, and an invisible anchor per
    community holds its terms together, so communities appear as groups."""
    G = nx.Graph()
    for i, t in enumerate(terms):
        c = term_comm.get((t["library"], t["term"]))
        G.add_node(("T", i))
        for g in set(t["genes"]):
            gi = gene_info.get(g)
            same = gi is not None and c is not None and gi.Community == c and not gi.Bridge
            G.add_edge(("T", i), ("G", g), weight=1.0 if same else 0.25)
        if c is not None:
            G.add_edge(("T", i), ("C", c), weight=0.6)
    n = G.number_of_nodes()
    pos = nx.spring_layout(G, weight="weight", seed=seed, iterations=300, k=1.6 / math.sqrt(max(n, 1)))
    pos = {k: tuple(v) for k, v in pos.items() if k[0] != "C"}
    # genes linked to exactly the same terms land on one spot: spread them on a small circle
    xs = [p[0] for p in pos.values()] or [0]
    ys = [p[1] for p in pos.values()] or [0]
    step = 0.028 * max(max(xs) - min(xs), max(ys) - min(ys), 1e-9)
    groups: Dict[frozenset, List] = {}
    for node in G.nodes:
        if node[0] == "G":
            groups.setdefault(frozenset(G.neighbors(node)), []).append(node)
    for members in groups.values():
        if len(members) < 2:
            continue
        members.sort()
        cx = sum(pos[m][0] for m in members) / len(members)
        cy = sum(pos[m][1] for m in members) / len(members)
        r = step * max(1.0, len(members) / 4)
        for j, m in enumerate(members):
            a = 2 * math.pi * j / len(members)
            pos[m] = (cx + r * math.cos(a), cy + r * math.sin(a))
    return pos


def _hull(points, pad: float):
    """Padded convex hull (closed polygon) around a set of points."""
    ring = []
    for x, y in points:
        for a in range(16):
            ring.append((x + pad * math.cos(a * math.pi / 8), y + pad * math.sin(a * math.pi / 8)))
    try:
        from scipy.spatial import ConvexHull

        h = ConvexHull(ring)
        poly = [ring[i] for i in h.vertices]
    except Exception:
        poly = ring
    return poly + poly[:1]


def network_figure(terms: List[Dict], modules=None, colour_by: str = "Community", labels: str = "Hub genes",
                   collection_colours: Optional[Dict[str, str]] = None, n_hub_labels: int = 15,
                   height: int = 850) -> go.Figure:
    """Plotly figure of the network.

    :param terms: kept terms, records with library, term, genes and optionally iteration, p, q.
    :param modules: NetworkModules computed on the same terms (None: no communities).
    :param colour_by: "Community" or "Collection".
    :param labels: "Hub genes", "Hub genes and terms" or "None".
    """
    terms = [t for t in terms if list(t.get("genes", []))]
    term_comm, gene_info, deg = _build(terms, modules)
    if not term_comm:
        colour_by = "Collection"
    pos = _layout(terms, term_comm, gene_info)
    collection_colours = collection_colours or {}
    libs = list(dict.fromkeys(t["library"] for t in terms))
    traces: List[go.Scatter] = []
    annotations = []

    # community areas
    xs = [p[0] for p in pos.values()] or [0]
    ys = [p[1] for p in pos.values()] or [0]
    pad = 0.03 * max(max(xs) - min(xs), max(ys) - min(ys), 1e-9)
    if term_comm:
        comms = list(modules.communities.itertuples())
        for row in comms:
            c = row.Community
            pts = [pos[("T", i)] for i, t in enumerate(terms) if term_comm.get((t["library"], t["term"])) == c]
            pts += [pos[("G", g)] for g in deg if g in gene_info and gene_info[g].Community == c
                    and not gene_info[g].Bridge]
            if not pts:
                continue
            poly = _hull(pts, pad)
            col = community_colour(c) if colour_by == "Community" else OTHER
            traces.append(go.Scatter(x=[p[0] for p in poly], y=[p[1] for p in poly], mode="lines", fill="toself",
                                     fillcolor=_rgba(col, 0.10), line=dict(width=1, color=_rgba(col, 0.35)),
                                     hoverinfo="skip", showlegend=False))
            top = max(poly, key=lambda p: p[1])
            annotations.append(dict(x=top[0], y=top[1], text=f"<b>{c}</b>", showarrow=False, yshift=8,
                                    font=dict(size=13, color=INK)))

    # edges (one trace)
    ex, ey = [], []
    for i, t in enumerate(terms):
        x0, y0 = pos[("T", i)]
        for g in set(t["genes"]):
            x1, y1 = pos[("G", g)]
            ex += [x0, x1, None]
            ey += [y0, y1, None]
    traces.append(go.Scatter(x=ex, y=ey, mode="lines", line=dict(width=0.8, color=EDGE), hoverinfo="skip",
                             showlegend=False))

    def term_hover(t):
        c = term_comm.get((t["library"], t["term"]))
        s = f"<b>{t['term']}</b><br>{t['library']}"
        if c:
            s += f"<br>Community {c}"
        if t.get("iteration") not in (None, ""):
            s += f"<br>Iteration {t['iteration']}"
        if t.get("p") is not None:
            s += f"<br>p = {t['p']:.2g}"
        if t.get("q") is not None:
            s += f"<br>q = {t['q']:.2g}"
        return s + f"<br>{len(set(t['genes']))} genes"

    def gene_hover(g):
        gi = gene_info.get(g)
        s = f"<b>{g}</b><br>Degree {deg[g]} (terms)"
        if gi is not None:
            s += f"<br>Community {gi.Community}"
            if gi.Bridge:
                s += f"<br>Bridge gene: {gi.Communities} communities (participation {gi.Participation:.2f})"
        return s

    hubs = [g for g in sorted(deg, key=lambda g: (-deg[g], g)) if deg[g] >= 2][:n_hub_labels]
    label_terms = labels == "Hub genes and terms"
    label_genes = labels != "None"
    gsize = lambda g: min(6 + 3 * deg[g], 30)
    tsize = lambda t: min(10 + 1.6 * math.sqrt(len(set(t["genes"]))) * 2, 30)

    def add_terms(idx, colour, name, group, show):
        if not idx:
            return
        traces.append(go.Scatter(
            x=[pos[("T", i)][0] for i in idx], y=[pos[("T", i)][1] for i in idx],
            mode="markers+text" if label_terms else "markers",
            text=[_short(terms[i]["term"]) for i in idx] if label_terms else None,
            textposition="bottom center", textfont=dict(size=10, color=INK),
            marker=dict(symbol="square", size=[tsize(terms[i]) for i in idx], color=colour,
                        line=dict(width=1.5, color="white")),
            hovertext=[term_hover(terms[i]) for i in idx], hoverinfo="text",
            name=name, legendgroup=group, showlegend=show))

    def add_genes(gs, colours, name, group, show, symbol="circle"):
        if not gs:
            return
        traces.append(go.Scatter(
            x=[pos[("G", g)][0] for g in gs], y=[pos[("G", g)][1] for g in gs],
            mode="markers+text" if label_genes else "markers",
            text=[g if g in hubs else "" for g in gs] if label_genes else None,
            textposition="top center", textfont=dict(size=11, color=INK),
            marker=dict(symbol=symbol, size=[gsize(g) for g in gs], color=colours,
                        line=dict(width=1.5 if symbol == "diamond" else 1, color=INK if symbol == "diamond" else "white")),
            hovertext=[gene_hover(g) for g in gs], hoverinfo="text",
            name=name, legendgroup=group, showlegend=show))

    bridges = [g for g in deg if g in gene_info and gene_info[g].Bridge]
    if colour_by == "Community":
        for _, row in modules.communities.iterrows():
            c = row["Community"]
            col = community_colour(c)
            idx = [i for i, t in enumerate(terms) if term_comm.get((t["library"], t["term"])) == c]
            gs = [g for g in deg if g in gene_info and gene_info[g].Community == c and not gene_info[g].Bridge]
            add_terms(idx, col, f"{c} · {_short(row['Lead term'], 36)} ({row['Terms']} terms, {row['Genes']} genes)", c, True)
            add_genes(gs, col, c, c, False)
        add_genes(bridges, [community_colour(gene_info[g].Community) for g in bridges],
                  f"Bridge genes ({len(bridges)})", "bridge", True, symbol="diamond")
        rest = [g for g in deg if g not in gene_info]
        add_genes(rest, OTHER, "Other genes", "other", False)
    else:
        for k, lib in enumerate(libs):
            col = collection_colours.get(lib) or COMMUNITY_COLOURS[k % len(COMMUNITY_COLOURS)]
            idx = [i for i, t in enumerate(terms) if t["library"] == lib]
            add_terms(idx, col, f"{lib} ({len(idx)})", lib, True)
        plain = [g for g in deg if g not in bridges]
        add_genes(plain, OTHER, "Genes", "genes", True)
        add_genes(bridges, OTHER, f"Bridge genes ({len(bridges)})", "bridge", True, symbol="diamond")

    fig = go.Figure(data=traces)
    fig.update_layout(
        template="plotly_white", height=height, margin=dict(l=10, r=10, t=10, b=10),
        xaxis=dict(visible=False), yaxis=dict(visible=False, scaleanchor="x"),
        legend=dict(title=dict(text="Community" if colour_by == "Community" else "Collection"),
                    font=dict(size=11, color=INK), itemsizing="constant", bgcolor="rgba(255,255,255,0.85)"),
        annotations=annotations, hoverlabel=dict(bgcolor="white", font=dict(color=INK)),
        dragmode="pan",
    )
    return fig
