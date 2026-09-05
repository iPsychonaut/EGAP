"""MerquryFK parsing and an explicit assembly-ranking rule (port candidate).

fsp runs FastK on the cleaned R1/R2 reads once per sample, then MerquryFK per
candidate assembly, producing ``<prefix>.qv`` and ``<prefix>.completeness.stats``.
Note what fsp does NOT do: its selection script (bin/select_best_assembly.sh)
ranks by BUSCO complete-single-copy with QUAST auN as the only tie-break. QV
is reported, not used for selection. Matthew's suggestion to put QV in the
EGAP vote goes beyond fsp.

EGAP already carries KMER_COMPLETENESS and QUAL_VAL keys in the per-sample
stats dict (bin/sample_tsv.py) that nothing populates. This module fills them
and replaces the implicit "first metric wins" tie-break in
compare_assemblies.py with a stated lexicographic key.

Tool requirement for real runs: fastk + merquryfk (both on bioconda:
``fastk``, ``merquryfk``). Not needed for these unit tests.
"""
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple


def _read_table(path) -> List[Dict[str, str]]:
    lines = [ln.rstrip("\n") for ln in Path(path).read_text().splitlines() if ln.strip()]
    header = lines[0].split("\t")
    rows = []
    for ln in lines[1:]:
        cells = ln.split("\t")
        rows.append(dict(zip(header, cells)))
    return rows


def parse_qv(path) -> Dict[str, float]:
    """Return {assembly_name: QV} from MerquryFK's .qv (columns: Assembly,
    No Support, Total, Error %, QV). Header-driven, so column order changes
    across MerquryFK versions do not break it."""
    out = {}
    for r in _read_table(path):
        name = r.get("Assembly") or r.get("assembly")
        qv = r.get("QV") or r.get("qv")
        if name is None or qv is None:
            raise ValueError(f"{path}: expected Assembly and QV columns, got {list(r)}")
        out[name] = float(qv)
    return out


def parse_completeness(path) -> Dict[str, float]:
    """Return {assembly_name: % k-mers covered} from .completeness.stats
    (columns: Assembly, Region, Found, Total, % Covered)."""
    out = {}
    for r in _read_table(path):
        name = r.get("Assembly") or r.get("assembly")
        cov = r.get("% Covered") or r.get("%Covered") or r.get("completeness")
        if name is None or cov is None:
            raise ValueError(f"{path}: expected Assembly and % Covered columns, got {list(r)}")
        out[name] = float(cov)
    return out


def qv_from_error_rate(error_fraction: float) -> float:
    """QV = -10 log10(error). Used to sanity-check parsed values."""
    import math
    if error_fraction <= 0:
        return float("inf")
    return -10.0 * math.log10(error_fraction)


# ---------------------------------------------------------------------------
# Ranking
# ---------------------------------------------------------------------------
class AssemblyStats:
    __slots__ = ("name", "busco1", "busco2", "contigs", "n50", "qv", "kmer_completeness")

    def __init__(self, name, busco1=None, busco2=None, contigs=None, n50=None,
                 qv=None, kmer_completeness=None):
        self.name = name
        self.busco1 = busco1
        self.busco2 = busco2
        self.contigs = contigs
        self.n50 = n50
        self.qv = qv
        self.kmer_completeness = kmer_completeness

    def as_dict(self):
        return {k: getattr(self, k) for k in self.__slots__}


def _f(x, default):
    return default if x is None else float(x)


def rank_key(a: AssemblyStats) -> Tuple:
    """Explicit lexicographic priority. Higher is better in every slot.

    1. mean BUSCO completeness over the two lineages (missing -> 0)
    2. MerquryFK QV (missing -> 0, so assemblies without QV never win a tie on it)
    3. k-mer completeness
    4. N50
    5. negative contig count
    6. name (deterministic last resort; alphabetical, documented)
    """
    b = (_f(a.busco1, 0.0) + _f(a.busco2, 0.0)) / 2.0
    return (round(b, 3), round(_f(a.qv, 0.0), 2), round(_f(a.kmer_completeness, 0.0), 2),
            _f(a.n50, 0.0), -_f(a.contigs, float("inf")), a.name)


def select_best(stats: Sequence[AssemblyStats]) -> AssemblyStats:
    if not stats:
        raise ValueError("no assemblies to rank")
    return max(stats, key=rank_key)


def egap_vote(stats: Sequence[AssemblyStats]) -> str:
    """Reference copy of EGAP's current per-metric vote so tests can show
    where the two rules disagree. Mirrors compare_assemblies.py exactly:
    metric order BUSCO_1, BUSCO_2, contigs (min), N50 (max); Counter tie ->
    first metric winner."""
    from collections import Counter
    winners = []
    for idx, attr in enumerate(("busco1", "busco2", "contigs", "n50")):
        cands = [(getattr(s, attr), s.name) for s in stats if getattr(s, attr) is not None]
        if not cands:
            continue
        pick = min(cands, key=lambda x: x[0]) if attr == "contigs" else max(cands, key=lambda x: x[0])
        winners.append(pick[1])
    return Counter(winners).most_common(1)[0][0]
