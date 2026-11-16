"""Network export for LLM interpretation: a semantic DOT file and the prompt that wraps it.

Two DOT files describe the same network:
  * the display DOT (colours, sizes, legend) for drawing, built by ui/dot_utils.merge_iterative_dot;
  * the semantic DOT (this module) for reading by people or language models: no styling, only meaning.

Semantic DOT
  graph attributes   input size, q-value threshold, random lists, network p-values, community method
  subgraph cluster_Cn  one per community (Leiden on the term network, see network_modules.py);
                     holds its terms and the genes that belong only to it
  term nodes         label, library, iteration, p, q, community, genes (= number of genes)
  gene nodes         degree (= number of terms; at most one per collection in iGEA), community;
                     bridge genes (terms in more than one community) sit outside every cluster and
                     carry bridge=true and their participation coefficient
  edges              one statement per term: "term" -- {"GENE1" "GENE2" ...};
Graphviz ignores attributes it does not know, so the file still renders (dot/fdp draw the clusters).
"""
from __future__ import annotations

import re
from typing import Dict, Iterable, List, Optional

LIBRARY_DESCRIPTIONS = {
    "H: Hallmark Gene Sets": "Curated gene sets representing well-defined biological states or processes",
    "C2: BioCarta": "Canonical pathways from the BioCarta database",
    "C2: KEGG MEDICUS": "Metabolic and signaling pathways from KEGG (MEDICUS)",
    "C2: KEGG Legacy": "Metabolic and signaling pathways from KEGG (legacy release)",
    "C2: Pathway Interaction Database": "Curated human signaling pathways",
    "C2: Reactome Pathways": "Expert-curated biological pathways",
    "C2: WikiPathways": "Community-curated biological pathways",
    "C5: Gene Ontology: Biological Process": "Biological processes from Gene Ontology",
    "C5: Gene Ontology: Cellular Component": "Cellular components from Gene Ontology",
    "C5: Gene Ontology: Molecular Function": "Molecular functions from Gene Ontology",
    "C5: Human Phenotype Ontology": "Human phenotypes and diseases",
    "Protein Interaction": "Protein-protein interaction modules",
}

METRIC_LABELS = {
    "Cluster Size": "Nodes (genes + terms) in the largest connected cluster",
    "Number of Genes": "Genes in the largest cluster",
    "Number of Terms": "Terms in the largest cluster",
    "Number of Edges": "Edges in the largest cluster",
    "Average Edges per Gene": "Average edges per gene",
}
SIGNIFICANCE_METRIC = "Cluster Size"


def _sid(s: str) -> str:
    return re.sub(r"_+", "_", re.sub(r"\W+", "_", s)).strip("_")


def _q(v: str) -> str:
    return '"' + str(v).replace('"', '\\"') + '"'


def _fmt(x, digits=2) -> str:
    try:
        x = float(x)
    except (TypeError, ValueError):
        return "NA"
    if x == 0:
        return "0"
    return f"{x:.{digits}g}"


def network_stats(table, attrs: Dict, observed: Dict, n_input_genes: int) -> Dict:
    """Collect the statistics shown to the LLM from ``per_run_null.run_statistics`` outputs."""
    metrics = []
    for _, r in table.iterrows():
        metrics.append(dict(metric=r["Metric"], observed=r["Observed"], random_mean=r["Null Mean"],
                            r=int(r["Random >= Observed"]), p=float(r["Empirical p"])))
    return dict(input_genes=int(n_input_genes), n_lists=int(attrs.get("n_lists")),
                q_threshold=attrs.get("q_threshold"), expected_chance_terms=attrs.get("expected_chance_terms"),
                n_terms_total=observed.get("n_terms_total"), n_terms_kept=observed.get("n_terms_kept"),
                metrics=metrics)


def _num(x) -> bool:
    try:
        return x is not None and float(x) == float(x)
    except (TypeError, ValueError):
        return False


def _p_text(m: Dict, n_lists: int) -> str:
    if m["r"] == 0:
        return f"p < {_fmt(1 / n_lists, 1)} (no random list reached this value)"
    return f"p = {_fmt(m['p'])}"


def _network_p(stats: Optional[Dict]):
    if not stats:
        return None
    for m in stats.get("metrics", []):
        if m["metric"] == SIGNIFICANCE_METRIC:
            return m
    return None


def semantic_dot(terms: Iterable[Dict], modules=None, stats: Optional[Dict] = None,
                 name: str = "enrichment_network") -> str:
    """Semantic DOT for the terms drawn in the network.

    :param terms: records with keys library, term (display name), iteration, p, q (optional), genes.
    :param modules: ``network_modules.NetworkModules`` computed on the same terms (optional).
    :param stats: output of :func:`network_stats` (optional).
    """
    terms = [t for t in terms if list(t.get("genes", []))]
    term_comm, gene_info = {}, {}
    if modules is not None and len(modules.terms):
        term_comm = {(r.Library, r.Term): r.Community for r in modules.terms.itertuples()}
        gene_info = {r.Gene: r for r in modules.genes.itertuples()}
    deg = {}
    for t in terms:
        for g in t["genes"]:
            deg[g] = deg.get(g, 0) + 1

    lines = [f"graph {name} {{"]
    if stats:
        lines.append(f"  // {stats['input_genes']} validated input genes; network = terms passing the term-level filter")
        attrs = [f"input_genes={stats['input_genes']}", f"random_lists={stats['n_lists']}"]
        qt = stats.get("q_threshold")
        attrs.append(f"q_threshold={_q('none' if qt is None or qt >= 1 else qt)}")
        m = _network_p(stats)
        if m:
            attrs.append(f"network_p={_q(('<' + _fmt(1 / stats['n_lists'], 1)) if m['r'] == 0 else _fmt(m['p'], 3))}")
        lines.append(f"  graph [{', '.join(attrs)}];")
    if modules is not None and modules.summary.get("n_terms"):
        s = modules.summary
        extra = (f", modularity={_q(s['modularity'])}" if s.get("modularity") is not None else "")
        lines.append(f"  // communities: {s['method']}{extra}")
        lines.append(f"  graph [community_method={_q(s['method'])}{extra}];")

    def term_node(t):
        a = [f"type=term", f"label={_q(t['term'])}", f"library={_q(t['library'])}"]
        if t.get("iteration") not in (None, ""):
            a.append(f"iteration={int(t['iteration'])}")
        if _num(t.get("p")):
            a.append(f"p={_q(_fmt(t['p']))}")
        if _num(t.get("q")):
            a.append(f"q={_q(_fmt(t['q']))}")
        c = term_comm.get((t["library"], t["term"]))
        if c:
            a.append(f"community={c}")
        a.append(f"genes={len(set(t['genes']))}")
        return f"{_q('term_' + _sid(t['term']))} [{', '.join(a)}];"

    def gene_node(g):
        a = ["type=gene", f"degree={deg[g]}"]
        gi = gene_info.get(g)
        if gi is not None:
            a.append(f"community={gi.Community}")
            if gi.Bridge:
                a += ["bridge=true", f"participation={_q(_fmt(gi.Participation))}"]
        return f"{_q(g)} [{', '.join(a)}];"

    if term_comm and modules is not None:
        for _, row in modules.communities.iterrows():
            c = row["Community"]
            lines.append(f"  subgraph cluster_{c} {{")
            lines.append(f"    label={_q(c)}; lead_term={_q(row['Lead term'])}; terms={row['Terms']}; genes={row['Genes']};")
            for t in terms:
                if term_comm.get((t["library"], t["term"])) == c:
                    lines.append("    " + term_node(t))
            for g in sorted(deg):
                gi = gene_info.get(g)
                if gi is not None and gi.Community == c and not gi.Bridge:
                    lines.append("    " + gene_node(g))
            lines.append("  }")
        bridges = [g for g in sorted(deg) if gene_info.get(g) is not None and gene_info[g].Bridge]
        if bridges:
            lines.append("  // bridge genes: connect terms in different communities")
            lines += ["  " + gene_node(g) for g in bridges]
        rest_t = [t for t in terms if (t["library"], t["term"]) not in term_comm]
        rest_g = [g for g in sorted(deg) if g not in gene_info]
    else:
        rest_t, rest_g = terms, sorted(deg)
    lines += ["  " + term_node(t) for t in rest_t]
    lines += ["  " + gene_node(g) for g in rest_g]
    lines.append("  // edges: each term with its genes")
    for t in terms:
        genes = " ".join(_q(g) for g in sorted(set(t["genes"])))
        lines.append(f"  {_q('term_' + _sid(t['term']))} -- {{{genes}}};")
    lines.append("}")
    return "\n".join(lines)


INTRO = ("Below is an iterative gene set enrichment (iGEA) network computed from a gene list. Iterative means that "
         "for each library source, the top enriched gene set is found, saved and its genes are removed from the gene "
         "list; the remaining genes are then tested again. Thus each gene appears at most once per library source but "
         "can be linked to terms from different libraries. The network of genes and gene sets (= terms) is given in "
         "DOT format, preceded by its statistics and communities.")


def _stats_section(stats: Optional[Dict], has_q: bool = True) -> str:
    if not stats:
        return ("**STATISTICAL CONTEXT:**\nNo statistics were computed for this network. Its terms were not filtered "
                "against random gene lists, so treat every interpretation as tentative.\n")
    qt = stats.get("q_threshold")
    n = stats["n_lists"]
    if qt is None or qt >= 1:
        term = (f"No term-level filter was applied (q-value threshold 1): all {stats['n_terms_total']} terms selected "
                f"by iGEA are in the network, which may include terms that random gene lists of this size also select by chance"
                + (f". Each term carries its q-value (an empirical false discovery rate estimated from {n:,} random gene "
                   f"lists of the same size, background and settings)." if has_q else "."))
    else:
        term = (f"iGEA selected {stats['n_terms_total']} terms; {stats['n_terms_kept']} passed the term-level filter "
                f"(q-value < {qt}, an empirical false discovery rate estimated from {n:,} random gene lists of the same "
                f"size, background and settings) and form this network. A random gene list of this size selects about "
                f"{stats['expected_chance_terms']} terms by chance.")
    rows = "\n".join(f"| {METRIC_LABELS.get(m['metric'], m['metric'])} | {_fmt(m['observed'], 4)} | "
                     f"{_fmt(m['random_mean'], 3)} | {_p_text(m, n)} |" for m in stats["metrics"])
    return f"""**STATISTICAL CONTEXT:**
The input contained {stats['input_genes']} validated genes. {term}

Network-level test: the largest connected cluster of this network compared with the networks built the same way (same filter) from the {n:,} random gene lists. Empirical p = (r + 1) / ({n:,} + 1), where r is the number of random lists with a value at least as large. The metrics describe the same cluster and are not independent tests.

| Metric | This network | Random mean | Empirical p |
|---|---|---|---|
{rows}

How to use the statistics: the network is not distinguishable from what random gene lists produce if the empirical p for the number of nodes in the largest cluster is 0.05 or higher. In that case, say so plainly, keep the rest short and tentative, and do not propose an experimental context. Do not use p-values, q-values or iteration numbers to rank biological importance.
"""


def _community_section(modules) -> str:
    if modules is None or not modules.summary.get("n_terms"):
        return ""
    s = modules.summary
    if s.get("modularity") is not None:
        how = (f"Terms were grouped into communities with {s['method']} on the term network, in which two terms are "
               f"linked when they share genes (modularity {s['modularity']}; the same partition was found with "
               f"{round(100 * s['partition_identical_fraction'])}% of {s['n_seeds']} random seeds).")
    else:
        how = (f"The network has too few terms for community detection; communities are its connected components.")
    rows = "\n".join(f"| {r['Community']} | {r['Terms']} | {r['Genes']} | {r['Lead term']} | "
                     f"{r['Top hub genes (degree)']} | {r['Bridge genes']} |" for _, r in modules.communities.iterrows())
    return f"""**COMMUNITIES:**
{how} Communities are algorithmic groupings, not biological labels; the cluster blocks in the DOT list each community's terms and genes. Gene degree is the number of terms a gene connects to; in iGEA this equals the number of library sources that place the gene in the network. Bridge genes connect terms in different communities.

| Community | Terms | Genes | Lead term (most genes) | Top hub genes (degree) | Bridge genes |
|---|---|---|---|---|---|
{rows}
"""


def _guide_and_task(libs: str, comm_guide: str) -> str:
    """Background on the DOT attributes, then the writing task (narrative interpretation)."""
    return f"""**NETWORK INTERPRETATION GUIDE (background):**
- **Term nodes** (`type=term`): enriched gene sets, with `library` (source collection), `iteration`, `p` (enrichment p-value in its iteration), `q` (term q-value) and `genes` (number of genes linked)
- **Gene nodes** (`type=gene`): input genes, with `degree` (number of terms the gene connects to; in iGEA, the number of library sources that place it in the network){comm_guide}
- **Edges**: each line `"term" -- {{"GENE1" "GENE2" ...}}` links a term to its genes. Edges are gene-set memberships from curated annotation, not physical interactions or directions of change.

**Library Sources:**
{libs}

**YOUR TASK:**
You are an experienced molecular biologist explaining to a colleague what this gene list reveals. Write in clear, natural prose. Explain the biology behind the network and why these processes appear together; do not describe the network itself. Your answer should be a conversation starter: it should make the reader want to dig further, and it ends by asking how the gene list was generated so that the interpretation can be rebuilt around the real experiment.

**How to write it:**
- Ground your claims in the network by naming a few key genes in parentheses, for example (GENE1, GENE2, GENE3). Do not quote term names verbatim or list memberships; describe them in plain words. The reader has the network.
- Use your biological knowledge to explain mechanisms. When a central point depends on knowledge that is not in the network, say so once in plain words ("this is not in the network, but ...").
- Do not repeat general caveats about annotation networks (memberships rather than interactions, no direction of change). Mention a limitation only where it changes the interpretation.
- Point out what is surprising or unexpected.
- About 900-1,200 words.

**Structure:**
1. **Executive summary**: the main findings in 4-6 plain sentences. Start with whether the network is distinguishable from random gene lists, then the overall biological picture and the most important modules and genes.
2. **Communities**: one short section per community (C1, C2, ...), headed by the community ID and a plain biological name you choose. In 3-5 sentences, explain the biology it represents in plain language, which genes anchor it and why they matter, and anything surprising. If a community looks like a mix of unrelated processes, say so; if two communities belong together, say so.
3. **Putting it together**: connected prose that combines the communities into one biological story. Give the key genes a role in that story: hubs that anchor a community, and bridge genes that link communities (what does each link suggest biologically?).
4. **Estimated experimental context**: hypothesize what experiment most likely generated this gene list (for example the type of assay, perturbation, cell type or condition). Reason from the network: which communities and genes point to it, and how confident you are. Name one or two alternatives and what in the network would distinguish them.
5. **Closing**: end with one or two sentences inviting the reader to share how the gene list was actually generated (assay, perturbation, cell type, comparison), so that you can rebuild the interpretation around it: what fits the experiment, what does not, what is already known and what may be new.
"""


def llm_prompt(dot_text: str, terms: Iterable[Dict], modules=None, stats: Optional[Dict] = None) -> str:
    """Full prompt: statistics, communities, the semantic DOT, an interpretation guide and the request."""
    terms = list(terms)
    libs_present = []
    for t in terms:
        if t["library"] not in libs_present:
            libs_present.append(t["library"])
    order = list(LIBRARY_DESCRIPTIONS)
    libs_present.sort(key=lambda l: order.index(l) if l in order else len(order))
    libs = "\n".join(f"- **{l}**: {LIBRARY_DESCRIPTIONS.get(l, 'Gene set collection')}" for l in libs_present)
    has_comm = modules is not None and modules.summary.get("n_terms")
    comm_guide = ("\n- **Communities**: each `subgraph cluster_Cn` block is one community; term and gene nodes carry "
                  "`community`. Bridge genes are listed outside the blocks with `bridge=true` and a `participation` "
                  "value (0 = all terms in one community; higher = spread over communities)." if has_comm else "")
    return f"""{INTRO}

{_stats_section(stats, any(_num(t.get("q")) for t in terms))}
{_community_section(modules)}
**NETWORK (DOT):**
{dot_text}

{_guide_and_task(libs, comm_guide)}"""


FOLLOWUP = """Here is how the gene list was generated: {context}

Rebuild your interpretation of the same network with this experiment in mind. Keep the same writing rules (natural prose, a few key genes in parentheses, outside knowledge flagged in plain words). Structure:
1. **Revised summary**: 4-6 sentences. How well does the network match the experiment, and did your estimate of the experimental context hold?
2. **What fits**: communities and genes that make sense for this experiment, and why. Separate what is well established for this experiment from what is plausible but less documented.
3. **What does not fit**: communities or genes that are unexpected for this experiment. Could they be technical (for example, effects of the assay itself), indirect, or a sign of new biology?
4. **What may be new**: the most interesting candidates for novel biology, with the network evidence for each (which community, which bridge or hub genes) and how confident you are.
5. **Priorities**: the 3 findings you would follow up first, and the experiment that would test each.
"""


def followup_prompt(context: str) -> str:
    """Second turn: the reader reveals how the gene list was generated; the analysis is rebuilt around it."""
    return FOLLOWUP.format(context=context.strip())


def export_for_llm(terms: List[Dict], modules=None, stats: Optional[Dict] = None,
                   name: str = "enrichment_network"):
    """Return (semantic DOT, prompt)."""
    d = semantic_dot(terms, modules, stats, name)
    return d, llm_prompt(d, terms, modules, stats)
