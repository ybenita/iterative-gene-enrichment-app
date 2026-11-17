"""
Streamlit section: network statistics against random gene lists generated for this run.

The random lists have the validated input size, are drawn from the run's background and are
analysed with the run's collections and parameters (see per_run_null.py). Statistics run automatically
after each new iGEA analysis with the settings chosen under Statistics (random gene lists) (render_statistics_settings),
because the term q-value filter decides which terms enter the network; the network shown and exported
(state.iter_dot) is rebuilt from the retained terms. Changed settings are applied with the Recompute button;
a q-value threshold of 1 keeps all selected terms.
"""
import logging
import time
from typing import Callable, Dict, Optional

import pandas as pd

from per_run_null import (DEFAULT_GRID_LISTS, DEFAULT_Q_THRESHOLD, DEFAULT_SEED, format_report,
                          run_statistics)

logger = logging.getLogger(__name__)


def _run_signature(iter_enrich: Dict) -> tuple:
    first = next(iter(iter_enrich.values()))
    return (tuple(sorted(first.gene_set.genes)), tuple(sorted(iter_enrich.keys())),
            first.background_gene_set.name, first.p_value_method_name, first.p_threshold,
            first.min_overlap, first.min_term_size, first.max_term_size, first.max_iterations,
            tuple((k, len(v.results)) for k, v in sorted(iter_enrich.items())))


def _subset(state) -> Dict:
    """iGEA results for the collections the statistics (and network) refer to: all by default, or
    the subset last chosen for the network."""
    iter_enrich = state.get("iter_enrich") or {}
    chosen = state.get("stats_libraries")
    if chosen:
        sub = {k: v for k, v in iter_enrich.items() if k in set(chosen)}
        if sub:
            return sub
    return iter_enrich


def _compute(st, state, sub: Dict, n_lists: int, q_threshold: float, seed: int,
             n_jobs: Optional[int], grid_lists: Optional[int]):
    """Run the two-layer statistics for the collections in ``sub`` and rebuild their networks."""
    first = next(iter(sub.values()))
    n_genes = first.gene_set.size
    libraries = [it.gene_set_library for it in sub.values()]
    records = []
    for c, it in enumerate(sub.values()):
        for r in it.results:
            records.append((c, r.get("Term_raw", r.get("Term")), r.get("Iteration"), r.get("iteration p-value"),
                            list(r.get("Genes removed for next iteration", []))))
    if not records:
        st.warning("No terms were selected in these collections, so there is nothing to compare.")
        return None
    params = dict(p_value_method=first.p_value_method_name, p_threshold=first.p_threshold,
                  min_overlap=first.min_overlap, min_term_size=first.min_term_size,
                  max_term_size=first.max_term_size, max_iterations=first.max_iterations)
    bar = st.progress(0.0, text="Statistics: analysing random gene lists…")
    t0 = time.time()

    def progress(done, total, label):
        bar.progress(done / total, text=f"Analysing {label}… {done:,}/{total:,} ({time.time() - t0:.0f} s)")
    table, null, terms, observed, attrs = run_statistics(
        records, libraries, first.background_gene_set.genes, n_genes, n_lists=int(n_lists), seed=int(seed),
        q_threshold=float(q_threshold), grid_lists=grid_lists, n_jobs=n_jobs, progress_callback=progress, **params)
    bar.empty()
    shown = dict(params, background=first.background_gene_set.name)
    name = state.get("gene_set_name") or "gene_list"
    report = (f"Gene list: {name}\nValidated input genes: {n_genes}\n"
              f"Largest cluster: {observed['cluster_size']} nodes, {observed['n_genes']} genes, "
              f"{observed['n_terms']} terms, {observed['n_edges']} edges, "
              f"{observed['avg_edges_per_gene']:.4f} edges per gene, {observed['n_libraries']} collections\n\n"
              + format_report(table, attrs, [lib.name for lib in libraries], shown, terms))
    stored = {"signature": _run_signature(sub), "settings": (int(n_lists), float(q_threshold), int(seed)),
              "libraries": list(sub.keys()), "table": table, "report": report,
              "terms": terms, "null": null, "observed": observed, "attrs": dict(attrs)}
    state.per_run_benchmark = stored
    if set(sub) == set(state.get("iter_enrich") or {}):
        state.per_run_full = stored          # all collections: used for the q-values in the iteration tables
    # Rebuild each library's network from the retained terms only
    kept = {(r["Library"], r["Term"]) for _, r in terms.iterrows() if r["Kept"]}
    if "iter_dot" in state:
        for lib_name, it in sub.items():
            res = [r for r in it.results if (lib_name, r.get("Term_raw", r.get("Term"))) in kept]
            state.iter_dot[lib_name] = it.to_dot(res)
    # Communities, hub genes and bridges of the retained network (descriptive)
    try:
        from network_modules import analyze_network_modules
        stored["modules"] = analyze_network_modules(
            [dict(library=lib_name, term=r.get("Term", ""), genes=r.get("Genes removed for next iteration", []))
             for lib_name, it in sub.items() for r in it.results
             if (lib_name, r.get("Term_raw", r.get("Term"))) in kept])
    except Exception:
        stored["modules"] = None
    return stored


def ensure_statistics(st, state, library_names) -> bool:
    """Make sure the statistics refer to exactly ``library_names`` (e.g. the collections chosen for
    the network). Recomputes with the current settings if needed. Returns True if recomputed."""
    iter_enrich = state.get("iter_enrich") or {}
    all_names = list(iter_enrich.keys())
    state.stats_libraries = None if set(library_names) == set(all_names) else list(library_names)
    sub = _subset(state)
    stored = state.get("per_run_benchmark")
    if stored and stored.get("signature") == _run_signature(sub):
        return False
    cfg = state.get("per_run_cfg") or {}
    with st.spinner("Recomputing statistics for the selected collections…"):
        n, q, sd = _settings(state, cfg.get("default_n_lists", 10000))
        _compute(st, state, sub, n, q, sd, cfg.get("n_jobs"), cfg.get("grid_lists"))
    return True


def render_q_threshold(st) -> None:
    """Term q-value threshold, shown next to the iteration p-value threshold."""
    st.number_input("Term q-value threshold (1 = keep all terms)", min_value=0.01, max_value=1.0,
                    value=DEFAULT_Q_THRESHOLD, step=0.01, format="%.2f", key="null_q_threshold",
                    help="Correction of the iteration p-values: the q-value is the empirical false discovery rate "
                         "estimated from the random gene lists. Terms with q below this threshold form the network; "
                         "random lists are filtered the same way.")


def run_statistics_after_analysis(st, state, default_n_lists: int = 10000, n_jobs: Optional[int] = None,
                                  grid_lists: Optional[int] = None):
    """Run the random-list statistics right after the iterations (all collections), with its own progress bar."""
    iter_enrich = state.get("iter_enrich") or {}
    state.per_run_cfg = {"default_n_lists": default_n_lists, "n_jobs": n_jobs, "grid_lists": grid_lists}
    state.stats_libraries = None
    state.per_run_benchmark = state.per_run_full = None
    if not any(it.results for it in iter_enrich.values()):
        return None
    n, q, sd = _settings(state, default_n_lists)
    st.caption(f"Statistics: analysing {n:,} random gene lists (plus lists for the q-value adjustment)…")
    return _compute(st, state, iter_enrich, n, q, sd, n_jobs, grid_lists)


def term_stats(state, library_name: str):
    """{iteration: (q-value, kept)} for one collection, from the all-collection statistics; and the q threshold."""
    stored = state.get("per_run_full")
    if not stored or "terms" not in stored:
        return None, None
    t = stored["terms"]
    t = t[t["Library"] == library_name]
    return ({int(r["Iteration"]): (r["q-value"], bool(r["Kept"])) for _, r in t.iterrows()},
            stored["attrs"].get("q_threshold"))


def render_statistics_settings(st, default_n_lists: int = 10000, max_n_lists: int = 100000) -> None:
    """Random-list settings, shown with the analysis parameters so the first run already uses them."""
    st.number_input("Number of random lists", min_value=100, max_value=max_n_lists, value=default_n_lists,
                    step=100, key="null_n_lists",
                    help="Random gene lists of your input size, drawn from your background and analysed with the "
                         "same collections and parameters. More lists give more precise q-values and p-values "
                         "(smallest attainable p = 1 / (lists + 1)) but take longer.")
    st.number_input("Random seed", min_value=0, value=DEFAULT_SEED, step=1, key="null_seed",
                    help="Fixes the random lists, so results are reproducible.")


def _settings(state, default_n_lists):
    return (int(state.get("null_n_lists", default_n_lists)), float(state.get("null_q_threshold", DEFAULT_Q_THRESHOLD)),
            int(state.get("null_seed", DEFAULT_SEED)))


def render_per_run_benchmark(st, state, download_link: Callable, default_n_lists: int = 10000,
                             max_n_lists: int = 100000, n_jobs: Optional[int] = None, grid_lists: Optional[int] = None,
                             seconds_per_list_hint: Optional[float] = None) -> None:
    """Render the statistics section for the current iGEA results."""
    iter_enrich = state.get("iter_enrich") or {}
    if not iter_enrich:
        st.info("Run iterative enrichment analysis first to compute network statistics.")
        return
    state.per_run_cfg = {"default_n_lists": default_n_lists, "n_jobs": n_jobs, "grid_lists": grid_lists}
    if state.get("stats_libraries") and not set(state.stats_libraries) <= set(iter_enrich):
        state.stats_libraries = None                     # new analysis: back to all collections
    sub = _subset(state)
    first = next(iter(sub.values()))
    stored = state.get("per_run_benchmark")
    if stored and stored.get("signature") != _run_signature(sub):
        state.per_run_benchmark = stored = None

    n_genes = first.gene_set.size
    st.caption(
        f"Your network is compared with networks built from random gene lists of exactly {n_genes} genes "
        f"(your validated input), drawn from your background ({first.background_gene_set.size} genes) and "
        "analysed with the same collections and parameters."
    )
    n_lists, q_threshold, seed = _settings(state, default_n_lists)
    if seconds_per_list_hint:
        st.caption(f"Estimated time: about {max(1, round(n_lists * seconds_per_list_hint / 60))} min "
                   "(longer for large gene lists).")
    changed = bool(stored) and stored.get("settings") not in (None, (n_lists, q_threshold, seed))
    note = st.empty()
    recompute = st.button("Recompute statistics", type="primary" if (changed or not stored) else "secondary",
                          help="Settings are under Statistics (random gene lists) in the Analysis tab.")
    if recompute or not stored:
        stored = _compute(st, state, sub, n_lists, q_threshold, seed, n_jobs, grid_lists)
        if stored and "iter_dot" in state:
            state.network_generated = False
            state.last_merged_dot = ""
        if recompute and stored:
            st.rerun()                        # refresh the q-values in the iteration tables above
    if stored and stored.get("settings") not in (None, (n_lists, q_threshold, seed)):
        u = stored["settings"]
        note.warning(f"The statistics below were computed with {u[0]:,} random lists, q < {u[1]:g} and seed {u[2]}. "
                     f"The settings are now {n_lists:,} lists, q < {q_threshold:g} and seed {seed}: press "
                     "**Recompute statistics** to apply them.")

    if not stored:
        return
    obs, attrs, table = stored["observed"], stored["attrs"], stored["table"]
    qt = attrs.get("q_threshold")
    st.caption(f"Network built from the {obs['n_terms_kept']} of {obs['n_terms_total']} terms with q < {qt:g} "
               "(see the q-values in the iteration tables above)." if qt is not None and qt < 1 else
               f"Network built from all {obs['n_terms_total']} selected terms (no q-value filter).")
    if set(stored.get("libraries", [])) != set(iter_enrich):
        st.warning("Statistics and network refer to the collections chosen for the network: "
                   + ", ".join(stored["libraries"]) + ". Random lists were analysed with the same collections.")
    st.success(f"Largest cluster: {obs['n_genes']} genes, {obs['n_terms']} terms, {obs['n_edges']} edges, "
               f"{obs['n_libraries']} collections. Compared with {attrs['n_lists']:,} random lists "
               f"(seed {attrs['seed']}, {attrs['seconds']} s).")
    show = table[["Metric", "Observed", "Null Mean", "Null SD", "Z-score", "Random >= Observed", "Empirical p"]].copy()
    show.columns = ["Metric", "Your network", "Random mean", "Random SD", "z-score", "Random lists ≥ yours", "Empirical p"]
    st.dataframe(show.round(4), use_container_width=True, hide_index=True)
    mods = stored.get("modules")
    if mods is not None and mods.summary.get("n_terms"):
        sm = mods.summary
        with st.expander(f"Communities and hub genes ({sm['n_communities']} communities)", expanded=True):
            st.caption(
                f"{sm['method']} on the term network (terms linked by shared genes). "
                + (f"Modularity {sm['modularity']}; partition agreement across {sm['n_seeds']} seeds "
                   f"(mean adjusted Rand index) {sm['partition_stability_mean_ARI']}. " if sm.get('modularity') is not None else "")
                + "Gene degree = number of terms the gene connects to (at most one per collection in iGEA). "
                  "Bridge genes connect terms in different communities. Descriptive only: no statistical test.")
            st.dataframe(mods.communities.drop(columns=["All terms"]), use_container_width=True, hide_index=True)
            st.dataframe(mods.genes[["Gene", "Degree", "Community", "Participation", "Bridge", "Terms"]],
                         use_container_width=True, hide_index=True)
            nm = state.get("gene_set_name") or "network"
            st.markdown(f"📄 {download_link(mods.communities.to_csv(sep=chr(9), index=False), nm + '_communities', 'tsv')} · "
                        f"{download_link(mods.terms.rename(columns={'Library': 'Collection'}).to_csv(sep=chr(9), index=False), nm + '_community_terms', 'tsv')} · "
                        f"{download_link(mods.genes.to_csv(sep=chr(9), index=False), nm + '_hub_genes', 'tsv')}",
                        unsafe_allow_html=True)
    with st.expander("📖 How to read these numbers"):
        n = attrs["n_lists"]
        st.markdown(f"""
- **Empirical p** = (r + 1) / (n + 1), where *r* is the number of random lists whose largest cluster
  was at least as large as yours and *n* = {n:,}. The smallest possible value is {1/(n+1):.2g}.
  This makes no assumption about the shape of the random distribution.
- **z-score** = (your value − random mean) / random SD, shown for reference. The random values are
  counts with skewed distributions, so read the empirical p rather than converting z to a probability.
- The five metrics describe the same cluster and are not independent tests.
- Random networks are built exactly like yours: their terms pass the same q-value filter first.
- A small network p means your genes form a more connected cross-collection network than random lists
  of the same size, background and settings.
""")
    name = state.get("gene_set_name") or "statistical_report"
    st.markdown(f"📄 {download_link(stored['report'], name + '_statistical_report', 'txt')} · "
                f"{download_link(stored['null'].to_csv(sep=chr(9), index=False), name + '_random_lists', 'tsv')}",
                unsafe_allow_html=True)


def network_from_state(state):
    """The network that is tested, drawn and sent to the LLM: the terms that passed the q-value filter
    (collections chosen for the network), the communities computed on them and the statistics.
    Returns dict(terms, modules, stats) or None if statistics are not available."""
    stored = state.get("per_run_benchmark")
    if not stored:
        return None
    try:
        from llm_export import network_stats
        sub = _subset(state)
        terms_df = stored["terms"]
        info = {(r["Library"], r["Term"]): (r["q-value"], bool(r["Kept"])) for _, r in terms_df.iterrows()}
        terms = []
        for lib_name, it in sub.items():
            for r in it.results:
                q, kept = info.get((lib_name, r.get("Term_raw", r.get("Term"))), (None, False))
                if kept:
                    terms.append(dict(library=lib_name, term=r.get("Term", ""), iteration=r.get("Iteration"),
                                      p=r.get("iteration p-value"), q=q,
                                      genes=r.get("Genes removed for next iteration", [])))
        n_genes = next(iter(sub.values())).gene_set.size
        stats = network_stats(stored["table"], stored["attrs"], stored["observed"], n_genes)
        return {"terms": terms, "modules": stored.get("modules"), "stats": stats}
    except Exception:
        return None


def llm_export_from_state(state, network=None):
    """(semantic DOT, LLM prompt) for the current network, or None if statistics are not available."""
    network = network or network_from_state(state)
    if not network:
        return None
    try:
        from llm_export import export_for_llm
        return export_for_llm(network["terms"], network["modules"], network["stats"])
    except Exception:
        return None
    try:
        from llm_export import export_for_llm, network_stats
        sub = _subset(state)
        terms_df = stored["terms"]
        info = {(r["Library"], r["Term"]): (r["q-value"], bool(r["Kept"])) for _, r in terms_df.iterrows()}
        terms = []
        for lib_name, it in sub.items():
            for r in it.results:
                q, kept = info.get((lib_name, r.get("Term_raw", r.get("Term"))), (None, False))
                if kept:
                    terms.append(dict(library=lib_name, term=r.get("Term", ""), iteration=r.get("Iteration"),
                                      p=r.get("iteration p-value"), q=q,
                                      genes=r.get("Genes removed for next iteration", [])))
        n_genes = next(iter(sub.values())).gene_set.size
        stats = network_stats(stored["table"], stored["attrs"], stored["observed"], n_genes)
        return export_for_llm(terms, stored.get("modules"), stats)
    except Exception:
        return None
