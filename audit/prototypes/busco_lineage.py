"""Taxonomy-walk BUSCO lineage selection (port candidate from fspassemblypipeline).

fsp implementation (subworkflows/local/select_best_assembly_and_qc/main.nf,
~15 lines of Groovy): for each sample walk family -> order -> class -> phylum,
lower-case the rank name, append "_odb12", accept the first name present in a
whitelist file (assets/fungi_busco_lineages.txt, 44 fungal lineages), else use
a generic fallback. The ranks come from the samplesheet; nothing resolves them
from a taxid.

This module is a dependency-free Python equivalent with two additions EGAP
needs that fsp lacks:

* ``select_lineage_pair`` produces the (specific, generic) pair EGAP's BUSCO
  column already expects, so it slots into sample_tsv without schema change.
* The whitelist can be loaded from a file OR built from a live
  ``busco --list-datasets`` / compleasm ``list`` dump, so it stays current.

Not included on purpose: NCBI taxid -> rank resolution. That needs network
(ete3 NCBITaxa or ``datasets summary taxonomy``) and belongs behind an
explicit flag, not inside a unit-tested pure function.
"""
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Set, Tuple

RANK_ORDER: Tuple[str, ...] = ("family", "order", "class", "phylum", "kingdom")


def load_whitelist(source) -> Set[str]:
    """Accept a path (one lineage per line, '#' comments ok) or an iterable
    of names. Names are normalised to lower-case with the _odbNN suffix kept."""
    if isinstance(source, (str, Path)):
        lines = Path(source).read_text().splitlines()
    else:
        lines = list(source)
    out = set()
    for ln in lines:
        ln = ln.strip()
        if not ln or ln.startswith("#"):
            continue
        out.add(ln.lower())
    return out


def _candidate(name: Optional[str], ext: str) -> Optional[str]:
    if not name or not str(name).strip():
        return None
    base = str(name).strip().lower().replace(" ", "_")
    return f"{base}_{ext}" if ext else base


def select_lineage(ranks: Dict[str, Optional[str]], whitelist: Set[str],
                   fallback: str, ext: str = "odb12",
                   rank_order: Sequence[str] = RANK_ORDER) -> Tuple[str, Optional[str]]:
    """Return (lineage, rank_used). rank_used is None when the fallback fired.

    ranks: {"family": "Hymenogastraceae", "order": "Agaricales", ...}
    """
    for rank in rank_order:
        cand = _candidate(ranks.get(rank), ext)
        if cand and cand in whitelist:
            return cand, rank
    return fallback, None


def select_lineage_pair(ranks: Dict[str, Optional[str]], whitelist: Set[str],
                        generic: str, ext: str = "odb12") -> Tuple[str, str]:
    """EGAP runs two lineages. Return (specific, generic) where specific is the
    most specific whitelisted rank that is NOT the generic itself; if nothing
    more specific exists, both slots hold the generic (EGAP tolerates that
    but it wastes a BUSCO run, so callers may dedupe)."""
    specific, _ = select_lineage(ranks, whitelist - {generic.lower()}, generic, ext)
    return specific, generic


def strip_odb(name: str) -> str:
    """'agaricales_odb12' -> 'agaricales'. EGAP's TSV historically stores the
    bare name and lets compleasm resolve the odb version."""
    if "_odb" in name:
        return name.rsplit("_odb", 1)[0]
    return name


def whitelist_from_busco_listing(text: str) -> Set[str]:
    """Parse ``busco --list-datasets`` style output: indented lines like
    '    - agaricales_odb12'. Robust to the tree formatting BUSCO prints."""
    out = set()
    for ln in text.splitlines():
        ln = ln.strip().lstrip("-").strip()
        if ln.endswith(tuple(f"_odb{n}" for n in range(9, 20))):
            out.add(ln.lower())
    return out
