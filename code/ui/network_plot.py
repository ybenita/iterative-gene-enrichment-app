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


def _wrap(s: str, width: int = 26, max_lines: int = 3) -> str:
    """Wrap a term name onto a few short lines (Plotly <br>)."""
    lines, cur = [], ""
    for w in s.split():
        if cur and len(cur) + 1 + len(w) > width:
            lines.append(cur)
            cur = w
        else:
            cur = f"{cur} {w}".strip()
    lines.append(cur)
    if len(lines) > max_lines:
        lines = lines[:max_lines]
        lines[-1] = lines[-1][: width - 1] + "…"
    return "<br>".join(lines)


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


LABEL_OPTIONS = ["Hub genes and lead terms", "All genes and lead terms", "None"]


def network_figure(terms: List[Dict], modules=None, colour_by: str = "Community",
                   labels: str = "Hub genes and lead terms", collection_colours: Optional[Dict[str, str]] = None,
                   n_hub_labels: int = 15, height: int = 900) -> go.Figure:
    """Plotly figure of the network.

    :param terms: kept terms, records with library, term, genes and optionally iteration, p, q.
    :param modules: NetworkModules computed on the same terms (None: no communities).
    :param colour_by: "Community" or "Collection".
    :param labels: one of LABEL_OPTIONS. Lead term = the community's term with the most genes (the
        "Lead term" of the communities table): named in the community's title and outlined in dark;
        all other term names are shown on hover.

    In community mode every element of a community (area, label, terms, genes including its bridge
    genes, edges) shares one legend group: double-clicking a community in the legend shows it alone.
    """
    terms = [t for t in terms if list(t.get("genes", []))]
    term_comm, gene_info, deg = _build(terms, modules)
    if not term_comm:
        colour_by = "Collection"
    by_comm = colour_by == "Community"
    pos = _layout(terms, term_comm, gene_info)
    collection_colours = collection_colours or {}
    libs = list(dict.fromkeys(t["library"] for t in terms))
    comm_rows = list(modules.communities.to_dict("records")) if term_comm else []
    lead = {r["Community"]: r["Lead term"] for r in comm_rows}
    traces: List[go.Scatter] = []

    def gene_comm(g):
        gi = gene_info.get(g)
        return gi.Community if gi is not None else None

    def term_c(i):
        return term_comm.get((terms[i]["library"], terms[i]["term"]))

    # ---- hover texts
    def term_hover(i):
        t = terms[i]
        c = term_c(i)
        s = f"<b>{t['term']}</b><br>{t['library']}"
        if c:
            s += f"<br>Community {c}" + (" (lead term)" if lead.get(c) == t["term"] else "")
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

    hubs = set([g for g in sorted(deg, key=lambda g: (-deg[g], g)) if deg[g] >= 2][:n_hub_labels])
    label_genes = labels != "None"
    label_lead = labels in ("Hub genes and lead terms", "All genes and lead terms")
    all_genes = labels == "All genes and lead terms"
    gsize = lambda g: min(6 + 3 * deg[g], 30)
    tsize = lambda i: min(10 + 3.2 * math.sqrt(len(set(terms[i]["genes"]))), 30)
    xs = [p[0] for p in pos.values()] or [0]
    ys = [p[1] for p in pos.values()] or [0]
    span = max(max(xs) - min(xs), max(ys) - min(ys), 1e-9)
    pad = 0.03 * span
    px = 950 / (span * 1.2)   # approximate pixels per layout unit (for title sizes)

    # ---- building blocks
    def area(c, colour, group):
        pts = [pos[("T", i)] for i in range(len(terms)) if term_c(i) == c]
        pts += [pos[("G", g)] for g in deg if gene_comm(g) == c and not gene_info[g].Bridge]
        if not pts:
            return
        poly = _hull(pts, pad)
        traces.append(go.Scatter(x=[p[0] for p in poly], y=[p[1] for p in poly], mode="lines", fill="toself",
                                 fillcolor=_rgba(colour, 0.10), line=dict(width=1, color=_rgba(colour, 0.35)),
                                 hoverinfo="skip", showlegend=False, legendgroup=group))
        title = f"<b>{c}</b>" + (f"<br>{_wrap(lead[c], 30)}" if label_lead and lead.get(c) else "")
        x, y, where = _title_spot(poly, title)
        traces.append(go.Scatter(x=[x], y=[y], mode="text", text=[title], textposition=f"{where} center",
                                 textfont=dict(size=13, color=INK), hoverinfo="skip", showlegend=False,
                                 legendgroup=group))
        # invisible point where the title ends, so the axes (also of an isolated community) include it
        n_lines = title.count("<br>") + 1
        dy = (n_lines * 17 + 6) / px
        traces.append(go.Scatter(x=[x], y=[y + dy if where == "top" else y - dy], mode="markers",
                                 marker=dict(size=1, opacity=0), hoverinfo="skip", showlegend=False,
                                 legendgroup=group))

    def _title_spot(poly, title):
        """Above or below the area, whichever covers fewer nodes."""
        lines = title.replace("<b>", "").replace("</b>", "").split("<br>")
        w = max(len(l) for l in lines) * 7.5 / px / 2
        h = (len(lines) * 17 + 6) / px
        best = None
        for where in ("top", "bottom"):
            p = max(poly, key=lambda q: q[1]) if where == "top" else min(poly, key=lambda q: q[1])
            y0, y1 = (p[1], p[1] + h) if where == "top" else (p[1] - h, p[1])
            hits = sum(1 for (qx, qy) in pos.values() if abs(qx - p[0]) <= w and y0 <= qy <= y1)
            if best is None or hits < best[0]:
                best = (hits, p[0], p[1], where)
        return best[1], best[2], best[3]

    def edges(idx, group):
        ex, ey = [], []
        for i in idx:
            x0, y0 = pos[("T", i)]
            for g in set(terms[i]["genes"]):
                x1, y1 = pos[("G", g)]
                ex += [x0, x1, None]
                ey += [y0, y1, None]
        if ex:
            traces.append(go.Scatter(x=ex, y=ey, mode="lines", line=dict(width=0.8, color=EDGE), hoverinfo="skip",
                                     showlegend=False, legendgroup=group))

    def term_marks(idx, colour, name, group, show):
        if not idx:
            return
        traces.append(go.Scatter(
            x=[pos[("T", i)][0] for i in idx], y=[pos[("T", i)][1] for i in idx], mode="markers",
            marker=dict(symbol="square", size=[tsize(i) for i in idx], color=colour,
                        line=dict(width=[2.5 if is_lead(i) else 1.5 for i in idx],
                                  color=[INK if is_lead(i) else "white" for i in idx])),
            hovertext=[term_hover(i) for i in idx], hoverinfo="text", name=name, legendgroup=group, showlegend=show))

    def is_lead(i):
        return label_lead and term_c(i) is not None and lead.get(term_c(i)) == terms[i]["term"]

    def gene_marks(gs, colours, name, group, show, symbol="circle"):
        if not gs:
            return
        diamond = symbol == "diamond"
        traces.append(go.Scatter(
            x=[pos[("G", g)][0] for g in gs], y=[pos[("G", g)][1] for g in gs],
            mode="markers+text" if label_genes else "markers",
            text=[g if (all_genes or g in hubs) else "" for g in gs] if label_genes else None,
            textposition="top center",
            textfont=dict(size=[11 if g in hubs else 9 for g in gs] if all_genes else 11, color=INK),
            marker=dict(symbol=symbol, size=[gsize(g) for g in gs], color=colours,
                        line=dict(width=1.5 if diamond else 1, color=INK if diamond else "white")),
            hovertext=[gene_hover(g) for g in gs], hoverinfo="text", name=name, legendgroup=group, showlegend=show))

    if by_comm:
        for r in comm_rows:
            c = r["Community"]
            col = community_colour(c)
            idx = [i for i in range(len(terms)) if term_c(i) == c]
            area(c, col, c)
            edges(idx, c)
            term_marks(idx, col, f"{c} · {_short(r['Lead term'], 40)} ({r['Terms']} terms, {r['Genes']} genes)", c, True)
            gene_marks([g for g in deg if gene_comm(g) == c and not gene_info[g].Bridge], col, c, c, False)
            gene_marks([g for g in deg if gene_comm(g) == c and gene_info[g].Bridge], col, c, c, False, "diamond")
        rest_i = [i for i in range(len(terms)) if not term_c(i)]
        edges(rest_i, "other")
        term_marks(rest_i, OTHER, "Other terms", "other", True)
        gene_marks([g for g in deg if g not in gene_info], OTHER, "Other genes", "other", False)
    else:
        for r in comm_rows:
            area(r["Community"], OTHER, "areas")
        edges(list(range(len(terms))), "edges")
        for k, lib in enumerate(libs):
            col = collection_colours.get(lib) or COMMUNITY_COLOURS[k % len(COMMUNITY_COLOURS)]
            idx = [i for i in range(len(terms)) if terms[i]["library"] == lib]
            term_marks(idx, col, f"{lib} ({len(idx)})", lib, True)
        bridges = [g for g in deg if g in gene_info and gene_info[g].Bridge]
        gene_marks([g for g in deg if g not in bridges], OTHER, "Genes", "genes", True)
        gene_marks(bridges, OTHER, f"Bridge genes ({len(bridges)})", "bridge", True, symbol="diamond")

    fig = go.Figure(data=traces)
    fig.update_layout(
        template="plotly_white", height=height, margin=dict(l=10, r=10, t=10, b=10),
        xaxis=dict(visible=False), yaxis=dict(visible=False, scaleanchor="x"),
        legend=dict(orientation="h", x=0, xanchor="left", y=-0.01, yanchor="top",
                    title=dict(text="Community: " if by_comm else "Collection: "),
                    font=dict(size=11, color=INK), itemsizing="constant"),
        hoverlabel=dict(bgcolor="white", font=dict(color=INK)),
        dragmode="pan",
    )
    return fig
