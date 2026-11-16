"""Communities and hub genes in an enrichment network (terms connected through shared genes).

The network is bipartite: genes connect only to terms. Communities are found on the term-term
projection (two terms are linked when they share genes, weighted by the number of shared genes)
with the Leiden algorithm (modularity objective), and genes are then assigned to the community of
their terms. In iGEA, terms from the same collection never share genes (genes are removed after each
iteration), so communities form only where different collections agree, and a gene's degree (the
number of terms it connects to) is at most the number of collections.

Per gene:
  degree          number of terms (in iGEA: number of collections placing the gene in the network)
  community       community holding most of the gene's terms (ties: the larger community)
  participation   Guimera & Amaral (2005) participation coefficient, 1 - sum_s (k_s / k)^2;
                  0 = all terms in one community, higher = spread across communities
  bridge          True when the gene's terms fall in more than one community

The analysis is descriptive: it adds no statistical test. Partitions are made reproducible with a
fixed seed; robustness is reported as agreement (adjusted Rand index) of partitions from
``n_seeds`` different seeds with the reported (highest-modularity) partition.
"""
from __future__ import annotations

import random
from collections import Counter, defaultdict
from dataclasses import dataclass, field
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

import pandas as pd

MIN_TERMS_DEFAULT = 6


@dataclass
class NetworkModules:
    terms: pd.DataFrame
    genes: pd.DataFrame
    communities: pd.DataFrame
    summary: Dict[str, object] = field(default_factory=dict)


def _projection(term_genes: Sequence[frozenset]) -> List[Tuple[int, int, int]]:
    by_gene = defaultdict(list)
    for i, gs in enumerate(term_genes):
        for g in gs:
            by_gene[g].append(i)
    w = Counter()
    for idx in by_gene.values():
        for a in range(len(idx)):
            for b in range(a + 1, len(idx)):
                w[(idx[a], idx[b])] += 1
    return [(a, b, c) for (a, b), c in sorted(w.items())]


def _components(n: int, edges) -> List[int]:
    parent = list(range(n))

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    for a, b, _ in edges:
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[max(ra, rb)] = min(ra, rb)
    return [find(i) for i in range(n)]


def _leiden_runs(n: int, edges, seed: int, n_seeds: int):
    """Partitions (list of membership lists) and their modularity, one per seed."""
    try:
        import igraph as ig

        g = ig.Graph(n=n, edges=[(a, b) for a, b, _ in edges])
        g.es["weight"] = [c for _, _, c in edges]
        runs = []
        for s in range(n_seeds):
            ig.set_random_number_generator(random.Random(seed + s))
            part = g.community_leiden(objective_function="modularity", weights="weight", n_iterations=-1)
            runs.append((list(part.membership), g.modularity(part.membership, weights="weight")))
        ig.set_random_number_generator(random)
        return runs, "Leiden (igraph, modularity)"
    except ImportError:  # fall back to networkx Louvain (same objective, no connectivity guarantee)
        import networkx as nx

        g = nx.Graph()
        g.add_nodes_from(range(n))
        g.add_weighted_edges_from(edges)
        runs = []
        for s in range(n_seeds):
            comms = nx.community.louvain_communities(g, weight="weight", seed=seed + s)
            mem = [0] * n
            for k, c in enumerate(comms):
                for v in c:
                    mem[v] = k
            runs.append((mem, nx.community.modularity(g, comms, weight="weight")))
        return runs, "Louvain (networkx, modularity)"


def _ari(a: Sequence[int], b: Sequence[int]) -> float:
    from math import comb

    n = len(a)
    if n < 2:
        return 1.0
    cont = Counter(zip(a, b))
    sa, sb = Counter(a), Counter(b)
    idx = sum(comb(v, 2) for v in cont.values())
    ea = sum(comb(v, 2) for v in sa.values())
    eb = sum(comb(v, 2) for v in sb.values())
    exp = ea * eb / comb(n, 2)
    mx = (ea + eb) / 2
    return 1.0 if mx == exp else (idx - exp) / (mx - exp)


def analyze_network_modules(
    terms: Iterable[Dict[str, object]],
    min_terms: int = MIN_TERMS_DEFAULT,
    seed: int = 2026,
    n_seeds: int = 100,
) -> NetworkModules:
    """Communities, hub genes and bridges for the terms of a network.

    :param terms: records with keys ``library``, ``term`` and ``genes`` (iterable of gene symbols);
        pass the terms that are drawn in the network (after the q-value filter).
    :param min_terms: below this number of terms, communities are the connected components
        (too few terms for a meaningful partition).
    """
    recs = [(str(t["library"]), str(t["term"]), frozenset(t["genes"])) for t in terms]
    recs = [r for r in recs if r[2]]
    n = len(recs)
    edges = _projection([r[2] for r in recs])

    if n == 0:
        empty = pd.DataFrame()
        return NetworkModules(empty, empty, empty, {"n_terms": 0, "method": "none"})

    if n < min_terms:
        membership = _components(n, edges)
        method = f"connected components (fewer than {min_terms} terms)"
        modularity, stability, identical = None, None, None
    else:
        runs, method = _leiden_runs(n, edges, seed, n_seeds)
        best = max(range(len(runs)), key=lambda i: (round(runs[i][1], 12), -i))
        membership, modularity = runs[best]
        aris = [_ari(membership, m) for m, _ in runs]
        stability = sum(aris) / len(aris)
        identical = sum(a > 1 - 1e-12 for a in aris) / len(aris)

    # gene -> term indices
    gene_terms = defaultdict(list)
    for i, (_, _, gs) in enumerate(recs):
        for g in gs:
            gene_terms[g].append(i)

    # order communities by number of genes (then terms, then first term name) and renumber C1, C2, ...
    raw_ids = sorted(set(membership))
    comm_terms = {c: [i for i in range(n) if membership[i] == c] for c in raw_ids}
    comm_genes = {c: set().union(*(recs[i][2] for i in comm_terms[c])) for c in raw_ids}
    order = sorted(raw_ids, key=lambda c: (-len(comm_genes[c]), -len(comm_terms[c]),
                                           min(recs[i][1] for i in comm_terms[c])))
    label = {c: f"C{k + 1}" for k, c in enumerate(order)}
    rank = {c: k for k, c in enumerate(order)}

    gene_rows = []
    for g in sorted(gene_terms):
        idx = gene_terms[g]
        k = len(idx)
        per = Counter(membership[i] for i in idx)
        top = max(per.items(), key=lambda kv: (kv[1], -rank[kv[0]]))[0]
        part = 1.0 - sum((v / k) ** 2 for v in per.values())
        gene_rows.append({
            "Gene": g,
            "Degree": k,
            "Collections": len({recs[i][0] for i in idx}),
            "Community": label[top],
            "Communities": len(per),
            "Participation": round(part, 4),
            "Bridge": len(per) > 1,
            "Terms": "; ".join(sorted(recs[i][1] for i in idx)),
        })
    genes = pd.DataFrame(gene_rows).sort_values(["Degree", "Participation", "Gene"],
                                                 ascending=[False, False, True]).reset_index(drop=True)

    wdeg = Counter()
    for a, b, c in edges:
        wdeg[a] += c
        wdeg[b] += c
    term_rows = []
    for i, (lib, term, gs) in enumerate(recs):
        c = membership[i]
        out = sum(1 for g in gs if any(membership[j] != c for j in gene_terms[g]))
        term_rows.append({
            "Library": lib, "Term": term, "Community": label[c], "Genes": len(gs),
            "Shared genes (weighted degree)": wdeg[i], "Genes bridging to other communities": out,
        })
    terms_df = pd.DataFrame(term_rows).sort_values(
        ["Community", "Genes", "Term"], ascending=[True, False, True],
        key=lambda s: s.map(lambda v: int(v[1:])) if s.name == "Community" else s).reset_index(drop=True)

    comm_rows = []
    for c in order:
        ti = sorted(comm_terms[c], key=lambda i: (-len(recs[i][2]), recs[i][1]))
        cg = genes[genes["Community"] == label[c]]
        hubs = cg.sort_values(["Degree", "Gene"], ascending=[False, True]).head(5)
        comm_rows.append({
            "Community": label[c],
            "Terms": len(ti),
            "Genes": len(comm_genes[c]),
            "Collections": len({recs[i][0] for i in ti}),
            "Lead term": recs[ti[0]][1],
            "Top hub genes (degree)": ", ".join(f"{r.Gene} ({r.Degree})" for r in hubs.itertuples()),
            "Bridge genes": int(cg["Bridge"].sum()),
            "All terms": "; ".join(recs[i][1] for i in ti),
        })
    communities = pd.DataFrame(comm_rows)

    summary = {
        "method": method,
        "n_terms": n,
        "n_genes": len(gene_terms),
        "n_communities": len(order),
        "modularity": None if modularity is None else round(modularity, 4),
        "partition_stability_mean_ARI": None if stability is None else round(stability, 4),
        "partition_identical_fraction": identical,
        "n_seeds": None if modularity is None else n_seeds,
        "seed": seed,
        "n_bridge_genes": int(genes["Bridge"].sum()),
        "max_gene_degree": int(genes["Degree"].max()),
    }
    return NetworkModules(terms_df, genes, communities, summary)


def write_network_modules(mods: NetworkModules, out_dir, prefix: str = "network_modules") -> List[str]:
    """Write ``<prefix>_communities.tsv``, ``_terms.tsv``, ``_genes.tsv`` and ``_summary.tsv``."""
    from pathlib import Path

    out = Path(out_dir)
    paths = []
    for name, df in (("communities", mods.communities), ("terms", mods.terms), ("genes", mods.genes)):
        p = out / f"{prefix}_{name}.tsv"
        df.to_csv(p, sep="\t", index=False)
        paths.append(str(p))
    p = out / f"{prefix}_summary.tsv"
    pd.Series(mods.summary).to_csv(p, sep="\t", header=False)
    paths.append(str(p))
    return paths
