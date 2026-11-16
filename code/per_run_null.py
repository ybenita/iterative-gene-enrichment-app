"""
Per-run statistics for iGEA: term-level empirical q-values and a network-level random-list null.

Two layers, applied identically to the input list and to every random list:

1. Term level (signal-adjusted empirical q-value). iGEA selects terms with the nominal p-value
   threshold. Random lists of the same size, drawn from the same background and analysed with the
   same collections and parameters, show how many terms are selected by chance and at which
   p-values. For a selected term with p-value p in collection c,
        q = E[# chance terms in c with p' <= p] / # selected terms in c with p' <= p,
   made monotone in p (as Benjamini-Hochberg). The expectation first assumes every input gene is
   noise. Genes assigned to terms with q < q_threshold are then treated as signal, the expected
   chance terms are re-estimated from random lists of the remaining number of genes, and this is
   repeated until the set of accepted terms no longer changes (a pi0-style adjustment for the genes
   that real terms have already removed). The number of remaining genes is rounded UP to the next
   point of a size grid, which can only over-estimate chance terms. Terms with q >= q_threshold are
   dropped before the network is built.

2. Network level. The cross-collection network is built from the retained terms and its largest
   connected cluster is measured. Each random list goes through exactly the same steps (its own
   terms are scored against the reference, leaving that list out, then filtered), so the null
   networks are built the same way as the observed one. For each metric:
        z = (observed - mean) / SD,  empirical p = (r + 1) / (n + 1)   (Phipson & Smyth 2010).

Random list i of size s is drawn with numpy.random.default_rng([seed, s, i]) from the sorted
background, so all results are reproducible. q_threshold = 1 keeps every selected term (q-values are
still computed and reported, with the signal adjustment at the default level 0.10); q_threshold = None
skips the term level entirely.
"""
import logging
import math
import os
import time
from multiprocessing import get_context
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

from fast_igea import FISHER, FastCollection
from network_connectivity_benchmark import NetworkConnectivityAnalyzer

logger = logging.getLogger(__name__)

DEFAULT_N_LISTS = 10000
DEFAULT_SEED = 2026
DEFAULT_Q_THRESHOLD = 0.10
DEFAULT_GRID_LISTS = None      # random lists per grid size for the signal adjustment; None = max(200, n_lists // 10)
GRID_STEPS = 20                # grid spacing = ceil(n / GRID_STEPS)
GRID_MIN_FRACTION = 0.25       # smallest grid size = this fraction of n

METRICS = [
    ("cluster_size", "Cluster Size"),
    ("n_genes", "Number of Genes"),
    ("n_terms", "Number of Terms"),
    ("n_edges", "Number of Edges"),
    ("avg_edges_per_gene", "Average Edges per Gene"),
]

# A term record: (collection index, term name, iteration, p-value, genes removed)
Record = Tuple[int, str, int, float, List[str]]


# ----------------------------------------------------------------------------- network metrics
def largest_cluster_metrics(records: List[dict]) -> Dict[str, float]:
    """Largest connected cluster of the gene-term network built from iGEA records (dicts)."""
    analyzer = NetworkConnectivityAnalyzer()
    analyzer.add_igea_results(records)
    clusters = analyzer.get_clusters()
    if not clusters:
        return {"cluster_size": 0, "n_genes": 0, "n_terms": 0, "n_edges": 0,
                "avg_edges_per_gene": 0.0, "n_libraries": 0, "n_terms_total": len(records)}
    c = clusters[0]
    return {"cluster_size": c["size"], "n_genes": c["n_genes"], "n_terms": c["n_terms"],
            "n_edges": c["n_edges"], "avg_edges_per_gene": c["avg_edges_per_gene"],
            "n_libraries": c["n_libraries"], "n_terms_total": len(records)}


def _to_dicts(records: Sequence[Record], names: Sequence[str]) -> List[dict]:
    return [{"Term": t, "Iteration": it, "Library": names[c], "Genes removed for next iteration": g}
            for c, t, it, p, g in records]


# ----------------------------------------------------------------------------- term-level q
def grid_sizes(n: int) -> List[int]:
    step = max(1, math.ceil(n / GRID_STEPS))
    low = max(1, math.ceil(n * GRID_MIN_FRACTION))
    return sorted({n} | set(range(n - step, low - 1, -step)))


class ChanceReference:
    """Pooled p-values of terms selected in random lists, per grid size and collection."""

    def __init__(self, n_collections: int):
        self.n_collections = n_collections
        self.p: Dict[int, List[np.ndarray]] = {}
        self.k: Dict[int, int] = {}

    def add(self, size: int, lists: Sequence[Sequence[Record]]):
        ps = [[] for _ in range(self.n_collections)]
        for recs in lists:
            for c, _, _, p, _ in recs:
                ps[c].append(p)
        self.p[size] = [np.sort(np.array(x)) for x in ps]
        self.k[size] = len(lists)

    def sizes(self) -> List[int]:
        return sorted(self.p)

    def expected(self, c: int, p: float, size: int, own: Optional[np.ndarray] = None) -> float:
        cnt = np.searchsorted(self.p[size][c], p, side="right")
        k = self.k[size]
        if own is not None:                       # leave the scored list itself out
            cnt -= int((own <= p).sum())
            k -= 1
        return cnt / max(k, 1)


def signal_adjusted_q(records: Sequence[Record], n_genes: int, ref: ChanceReference,
                      q_threshold: float, leave_out_self: bool = False,
                      max_rounds: int = 10) -> Tuple[np.ndarray, int]:
    """Signal-adjusted empirical q-values for one list's selected terms.

    Returns (q per record, number of genes treated as noise in the final round).
    """
    if not records:
        return np.array([]), n_genes
    sizes = ref.sizes()
    by_c: Dict[int, List[int]] = {}
    for i, r in enumerate(records):
        by_c.setdefault(r[0], []).append(i)
    pvals = np.array([r[3] for r in records])
    signal: set = set()
    q = np.ones(len(records))
    m0 = n_genes
    for _ in range(max_rounds):
        m0 = n_genes - len(signal)
        size = next((s for s in sizes if s >= m0), sizes[-1])           # round UP (conservative)
        for c, idx in by_c.items():
            op = pvals[idx]
            own = op if (leave_out_self and size == n_genes) else None
            raw = np.array([ref.expected(c, p, size, own) / (op <= p).sum() for p in op])
            order = np.argsort(op, kind="stable")
            qc = np.empty(len(op))
            qc[order] = np.minimum.accumulate(np.clip(raw[order], 0, 1)[::-1])[::-1]
            q[idx] = qc
        new_signal = set().union(*[set(records[i][4]) for i in np.flatnonzero(q < q_threshold)]) \
            if (q < q_threshold).any() else set()
        if new_signal == signal:
            break
        signal = new_signal
    return q, m0


# ----------------------------------------------------------------------------- workers
_W: Dict = {}


def _init_gen(collections, background, seed, params):
    _W.update(collections=collections, background=background, seed=seed, params=params)


def _gen(task: Tuple[int, int]) -> List[Record]:
    size, i = task
    rng = np.random.default_rng([_W["seed"], size, i])
    genes = [_W["background"][j] for j in rng.choice(len(_W["background"]), size, replace=False)]
    out = []
    for c, col in enumerate(_W["collections"]):
        for it, term, p, removed in col.iterate(genes, **_W["params"]):
            out.append((c, term, it, p, removed))
    return out


def _keep(q: float, q_threshold: Optional[float]) -> bool:
    return q_threshold is None or q_threshold >= 1 or q < q_threshold


def _accept_level(q_threshold: Optional[float]) -> float:
    """q level at which genes of accepted terms are treated as signal."""
    return DEFAULT_Q_THRESHOLD if (q_threshold is None or q_threshold >= 1) else q_threshold


def _init_score(names, ref, n_genes, q_threshold):
    _W.update(names=names, ref=ref, n_genes=n_genes, q_threshold=q_threshold)


def _score(item: Tuple[int, List[Record]]) -> Dict[str, float]:
    i, recs = item
    if _W["q_threshold"] is None:
        kept = recs
    else:
        q, _ = signal_adjusted_q(recs, _W["n_genes"], _W["ref"], _accept_level(_W["q_threshold"]), leave_out_self=True)
        kept = [r for r, qq in zip(recs, q) if _keep(qq, _W["q_threshold"])]
    out = largest_cluster_metrics(_to_dicts(kept, _W["names"]))
    out.update(list=i, n_terms_total=len(recs), n_terms_kept=len(kept))
    return out


def _pool_map(fn, items, init, initargs, n_jobs, progress=None, label=""):
    items = list(items)
    if n_jobs == 1:
        init(*initargs)
        res = []
        for j, it in enumerate(items):
            res.append(fn(it))
            if progress and (j + 1) % max(1, len(items) // 20) == 0:
                progress(j + 1, len(items), label)
        return res
    with get_context("spawn").Pool(n_jobs, initializer=init, initargs=initargs) as pool:
        res = []
        for r in pool.imap(fn, items, chunksize=max(1, len(items) // (n_jobs * 20))):
            res.append(r)
            if progress and len(res) % max(1, len(items) // 20) == 0:
                progress(len(res), len(items), label)
        return res


# ----------------------------------------------------------------------------- optional disk cache
CACHE_BLOCK = 500


def _cache_key(background, collections, params, seed) -> str:
    """Identifies the random-list analysis: background, collections (terms and their genes), parameters, seed."""
    import hashlib
    h = hashlib.sha1()
    h.update("\n".join(background).encode())
    for col in collections:
        h.update(col.name.encode())
        h.update("\n".join(col.term_names).encode())
        h.update(np.asarray(col.K).tobytes())
        h.update("\n".join(col.genes).encode())
    h.update(repr(sorted(params.items())).encode())
    h.update(str(seed).encode())
    return h.hexdigest()[:16]


def _generate_cached(tasks, collections, background, seed, params, n_jobs, progress, cache_dir):
    """Generate random-list records in blocks saved to ``cache_dir``; finished blocks are reused, so an
    interrupted run resumes where it stopped. Results are identical with or without the cache."""
    import pickle
    os.makedirs(cache_dir, exist_ok=True)
    key = _cache_key(background, collections, params, seed)
    blocks = [tasks[k:k + CACHE_BLOCK] for k in range(0, len(tasks), CACHE_BLOCK)]
    out, done = [], 0
    pool = None
    try:
        for b in blocks:
            f = os.path.join(cache_dir, f"igea_null_{key}_{b[0][0]}_{b[0][1]}_{len(b)}.pkl")
            if os.path.exists(f):
                with open(f, "rb") as fh:
                    res = pickle.load(fh)
            else:
                if n_jobs == 1:
                    _init_gen(collections, background, seed, params)
                    res = [_gen(t) for t in b]
                else:
                    if pool is None:
                        pool = get_context("spawn").Pool(n_jobs, initializer=_init_gen,
                                                          initargs=(collections, background, seed, params))
                    res = pool.map(_gen, b, chunksize=max(1, len(b) // (n_jobs * 4)))
                with open(f + ".tmp", "wb") as fh:
                    pickle.dump(res, fh, protocol=pickle.HIGHEST_PROTOCOL)
                os.replace(f + ".tmp", f)
            out.extend(res)
            done += len(b)
            if progress:
                progress(done, len(tasks), "random lists")
    finally:
        if pool is not None:
            pool.close()
            pool.join()
    return out


# ----------------------------------------------------------------------------- main entry point
def run_statistics(observed_records: Sequence[Record], libraries, background_genes, n_genes: int,
                   p_value_method: str = FISHER, p_threshold: float = 0.01, min_overlap: int = 3,
                   min_term_size: int = 10, max_term_size: int = 600, max_iterations: Optional[int] = None,
                   n_lists: int = DEFAULT_N_LISTS, seed: int = DEFAULT_SEED,
                   q_threshold: Optional[float] = DEFAULT_Q_THRESHOLD, grid_lists: Optional[int] = DEFAULT_GRID_LISTS,
                   n_jobs: Optional[int] = None, progress_callback=None, cache_dir: Optional[str] = None):
    """Term q-values and network statistics for one analysis.

    Args:
        observed_records: the input list's iGEA terms as (collection index, term, iteration,
            p-value, genes removed); collection index refers to ``libraries``.
        libraries: GeneSetLibrary objects used in the analysis (all of them, including any in
            which no term was selected).
        background_genes: the analysis background (random lists are drawn from it).
        n_genes: validated input size.
        q_threshold: terms with signal-adjusted q >= this are dropped before building networks
            (None = keep all selected terms).
        cache_dir: optional folder for the generated random lists; an interrupted run resumes from it.

    Returns:
        (table, null, observed_terms, observed_metrics, attrs)
    """
    t0 = time.time()
    background = sorted(set(background_genes))
    if n_genes > len(background):
        raise ValueError(f"Input size {n_genes} exceeds background size {len(background)}")
    names = [lib.name for lib in libraries]
    collections = [FastCollection.from_library(lib, background, min_term_size, max_term_size) for lib in libraries]
    params = dict(p_threshold=p_threshold, min_overlap=min_overlap, max_iterations=max_iterations, method=p_value_method)
    n_jobs = n_jobs or os.cpu_count() or 1
    grid_lists = grid_lists or max(200, n_lists // 10)
    tasks = [(n_genes, i) for i in range(n_lists)]
    grid = grid_sizes(n_genes) if q_threshold is not None else [n_genes]
    tasks += [(s, i) for s in grid if s != n_genes for i in range(grid_lists)]
    stage = lambda d, n, lab: progress_callback(d, n, lab) if progress_callback else None
    if cache_dir:
        generated = _generate_cached(tasks, collections, background, seed, params, n_jobs, stage, cache_dir)
    else:
        generated = _pool_map(_gen, tasks, _init_gen, (collections, background, seed, params), n_jobs,
                              stage, "random lists")
    main = generated[:n_lists]
    ref = ChanceReference(len(collections))
    ref.add(n_genes, main)
    pos = n_lists
    for s in grid:
        if s != n_genes:
            ref.add(s, generated[pos:pos + grid_lists])
            pos += grid_lists
    scored = _pool_map(_score, list(enumerate(main)), _init_score, (names, ref, n_genes, q_threshold),
                       n_jobs, stage, "random networks")
    null = pd.DataFrame(scored).sort_values("list").reset_index(drop=True)

    obs = list(observed_records)
    if q_threshold is None:
        q_obs, m0 = np.full(len(obs), np.nan), n_genes
        kept = obs
    else:
        q_obs, m0 = signal_adjusted_q(obs, n_genes, ref, _accept_level(q_threshold))
        kept = [r for r, q in zip(obs, q_obs) if _keep(q, q_threshold)]
    observed = largest_cluster_metrics(_to_dicts(kept, names))
    observed["n_terms_total"], observed["n_terms_kept"] = len(obs), len(kept)
    obs_terms = pd.DataFrame([{"Library": names[c], "Iteration": it, "Term": t, "p-value": p, "q-value": q,
                               "Kept": _keep(q, q_threshold),
                               "Genes": ", ".join(g)} for (c, t, it, p, g), q in zip(obs, q_obs)])
    table = compare_to_null(observed, null)
    attrs = dict(n_genes=n_genes, n_lists=n_lists, seed=seed, n_jobs=n_jobs, background_size=len(background),
                 q_threshold=q_threshold, grid_sizes=grid, grid_lists=grid_lists if q_threshold is not None else 0,
                 noise_genes=m0, expected_chance_terms=round(float(np.mean([len(m) for m in main])), 2),
                 seconds=round(time.time() - t0, 1))
    null.attrs.update(attrs)
    logger.info(f"Per-run statistics: {n_lists} random lists of {n_genes} genes (+{len(tasks) - n_lists} for q) "
                f"in {attrs['seconds']} s ({n_jobs} processes)")
    return table, null, obs_terms, observed, attrs


def generate_null(libraries, background_genes, n_genes: int, n_lists: int = DEFAULT_N_LISTS,
                  seed: int = DEFAULT_SEED, n_jobs: Optional[int] = None, progress_callback=None,
                  q_threshold: Optional[float] = None, **params) -> pd.DataFrame:
    """Network metrics of random lists only (no observed list). Kept for scripts and tests."""
    _, null, _, _, _ = run_statistics([], libraries, background_genes, n_genes, n_lists=n_lists, seed=seed,
                                      n_jobs=n_jobs, q_threshold=q_threshold,
                                      progress_callback=(lambda d, n, lab: progress_callback(d, n)) if progress_callback else None,
                                      **params)
    return null


def compare_to_null(observed: Dict[str, float], null: pd.DataFrame) -> pd.DataFrame:
    """z-score and empirical p-value for each largest-cluster metric."""
    rows = []
    n = len(null)
    for key, label in METRICS:
        x = null[key].to_numpy(dtype=float)
        obs = float(observed[key])
        sd = x.std()
        r = int((x >= obs).sum())
        rows.append({
            "Metric": label, "Observed": obs, "Null Mean": x.mean(), "Null SD": sd, "Null Max": x.max(),
            "Z-score": (obs - x.mean()) / sd if sd > 0 else float("nan"),
            "Random >= Observed": r, "Random Lists": n, "Empirical p": (r + 1) / (n + 1),
        })
    return pd.DataFrame(rows)


def format_report(table: pd.DataFrame, attrs: Dict, libraries: List[str], params: Dict,
                  observed_terms: Optional[pd.DataFrame] = None) -> str:
    """Plain-text statistical report."""
    n = int(table["Random Lists"].iloc[0])
    qt = attrs.get("q_threshold")
    lines = [
        "=" * 100,
        "iGEA STATISTICAL REPORT (per-run random lists)",
        "=" * 100,
        f"Random lists: {n}, each of {attrs['n_genes']} genes (= validated input size), drawn without replacement",
        f"from the analysis background ({attrs['background_size']} genes); seed {attrs['seed']}; {attrs['seconds']} s on {attrs['n_jobs']} process(es).",
        "Each random list was analysed with the same collections and parameters as the input:",
        "  " + ", ".join(f"{k}={v}" for k, v in params.items()),
        "Collections: " + "; ".join(libraries),
        "",
    ]
    if qt is not None and observed_terms is not None:
        kept = int(observed_terms["Kept"].sum())
        lines += [
            "TERM LEVEL (signal-adjusted empirical q-value)",
            f"  A random list of this size selects on average {attrs['expected_chance_terms']} terms by chance.",
            f"  q = expected chance terms with p' <= p / selected terms with p' <= p (per collection, monotone),",
            f"  re-estimated after treating genes of accepted terms as signal ({attrs['noise_genes']} of {attrs['n_genes']} genes",
            f"  treated as noise in the final round; {attrs['grid_lists']} random lists per size step).",
            (f"  Terms kept (q < {qt}): {kept} of {len(observed_terms)}. Random lists are filtered the same way."
             if qt < 1 else f"  No term filter (q threshold {qt}): all {kept} selected terms kept, in the input and the random lists."),
            "",
        ]
    lines += [
        "NETWORK LEVEL (largest connected cluster)",
        "  z = (observed - random mean) / random SD",
        f"  empirical p = (r + 1) / (n + 1), r = random lists with a value >= observed, n = {n}",
        f"  (smallest attainable p = {1/(n+1):.4g}). The metrics describe the same cluster and are not independent tests.",
        "",
        f"  {'Metric':<26}{'Observed':>10}{'Rand mean':>11}{'Rand SD':>9}{'Rand max':>10}{'z':>8}{'r':>6}{'Emp. p':>10}",
        "  " + "-" * 90,
    ]
    for _, r in table.iterrows():
        lines.append(f"  {r['Metric']:<26}{r['Observed']:>10.4g}{r['Null Mean']:>11.3f}{r['Null SD']:>9.3f}"
                     f"{r['Null Max']:>10.4g}{r['Z-score']:>8.2f}{r['Random >= Observed']:>6d}{r['Empirical p']:>10.4g}")
    if observed_terms is not None and len(observed_terms):
        lines += ["", "TERMS", f"  {'Library':<40}{'It':>4}  {'p-value':>10}{'q-value':>9}  Kept  Term"]
        for _, r in observed_terms.sort_values(["Library", "Iteration"]).iterrows():
            lines.append(f"  {r['Library'][:39]:<40}{r['Iteration']:>4}  {r['p-value']:>10.3g}{r['q-value']:>9.3g}  "
                         f"{'yes ' if r['Kept'] else 'no  '}  {r['Term']}")
    return "\n".join(lines) + "\n"
