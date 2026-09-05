"""Tiara + FCS-GX consensus with chimera flagging (port candidate).

What fsp actually does (bin/comparison.R, subworkflows/local/contamination_detection):
  1. FCS-GX taxonomy.rpt rows are collapsed per contig. FCS-GX itself splits
     long sequences into windows and suffixes the id with "~start-end"; the R
     script groups on the id before "~". If the windows disagree on species
     the contig is tagged "possible chimera"; if they disagree on domain,
     "possible cross-domain chimera". The chimera detection is therefore
     entirely FCS-GX's per-window calls; Tiara plays no part in it.
  2. Tiara class_fst_stage is joined per contig.
  3. match = "match" | "mismatch" | "chimera". Contigs < 1 kb are forced to
     "unknown".
  4. Output is a BlobToolKit taxonomy table and a fake hits table. NOTHING IS
     REMOVED. fsp's consensus is a visualisation aid, not a filter.

EGAP's decontaminate_assembly.py removes contigs on Tiara alone.

This module reproduces the fsp collapse/compare logic in Python and then adds
the decision layer EGAP needs. The decision rule is deliberately conservative:

  keep    both classifiers say target domain, or either says unknown/absent
  remove  both classifiers agree on a non-target domain
  flag    classifiers disagree, or chimera -> keep, write to a review table

A "remove on either" rule was rejected: Tiara alone already produces the
>50 % removal warnings in EGAP's own README FAQ, and FCS-GX's fungal
divisions are coarse. Two-classifier agreement is the minimum bar.
"""
import re
from collections import defaultdict
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple

TARGET_DOMAIN = {
    "funga": "eukarya", "flora": "eukarya", "fauna": "eukarya",
    "bacteria": "bacteria", "archaea": "archaea",
}

# FCS-GX "div" strings start with a lineage prefix; fsp maps by prefix.
_DIV_PREFIX = (
    ("prok", "bacteria"), ("bact", "bacteria"), ("arch", "archaea"),
    ("fung", "eukarya"), ("plnt", "eukarya"), ("anml", "eukarya"),
    ("prst", "eukarya"), ("virs", "virus"),
)

UNKNOWN_TAXID = "32644"  # NCBI "unidentified"


def fcs_div_to_domain(div: Optional[str]) -> str:
    if not div:
        return "unknown"
    d = div.strip().lower()
    if d in ("unassigned", "unknown", ""):
        return "unknown"
    for prefix, dom in _DIV_PREFIX:
        if d.startswith(prefix):
            return dom
    return "unknown"


def tiara_to_domain(label: Optional[str]) -> str:
    """Map Tiara classes onto the same domain vocabulary. 'prokarya' is an
    unresolved bacteria/archaea call; 'organelle' is eukaryote-owned."""
    if not label:
        return "unknown"
    l = label.strip().lower()
    return {"eukarya": "eukarya", "organelle": "eukarya", "bacteria": "bacteria",
            "archaea": "archaea", "prokarya": "prokarya"}.get(l, "unknown")


# ---------------------------------------------------------------------------
# Parsers
# ---------------------------------------------------------------------------
def parse_tiara(path) -> Dict[str, str]:
    out = {}
    with open(path) as fh:
        for ln in fh:
            ln = ln.rstrip("\n")
            if not ln or ln.startswith("#"):
                continue
            parts = ln.split("\t")
            if parts[0] == "sequence_id":       # header guard EGAP lacks
                continue
            if len(parts) >= 2:
                out[parts[0].strip()] = parts[1].strip().lower()
    return out


def parse_fcsgx_rpt(path) -> List[Dict[str, str]]:
    """Header-driven parse of taxonomy.rpt. Lines starting '##' are metadata,
    the '#seq-id ...' line is the header."""
    rows = []
    header = None
    with open(path) as fh:
        for ln in fh:
            ln = ln.rstrip("\n")
            if not ln:
                continue
            if ln.startswith("##"):
                continue
            if ln.startswith("#"):
                header = ln[1:].split("\t")
                continue
            if header is None:
                raise ValueError(f"{path}: data before header")
            rows.append(dict(zip(header, ln.split("\t"))))
    return rows


_SEG_RE = re.compile(r"^(.*?)~(\d+)[-.]+(\d+)$")


def split_segment_id(seq_id: str) -> Tuple[str, Optional[Tuple[int, int]]]:
    m = _SEG_RE.match(seq_id)
    if not m:
        return seq_id, None
    return m.group(1), (int(m.group(2)), int(m.group(3)))


def collapse_fcsgx(rows: Iterable[Dict[str, str]]) -> Dict[str, dict]:
    """Per-contig summary mirroring comparison.R Step 3B."""
    groups = defaultdict(list)
    for r in rows:
        base, seg = split_segment_id(r["seq-id"])
        groups[base].append((r, seg))
    out = {}
    for base, items in groups.items():
        domains, species, taxids = set(), set(), []
        seq_len = 0
        for r, seg in items:
            seq_len += int(float(r.get("seq-len", 0) or 0))
            dom = fcs_div_to_domain(r.get("div"))
            sp = (r.get("tax-name-1") or r.get("best-tax-name") or "").strip()
            if dom != "unknown":
                domains.add(dom)
            if sp and sp.lower() != "unknown":
                species.add(sp)
                taxids.append(r.get("tax-id-1") or r.get("best-taxid"))
        if len(domains) > 1:
            call, dom = "possible cross-domain chimera", "mixed"
        elif len(species) > 1:
            call, dom = "possible chimera", next(iter(domains)) if domains else "unknown"
        elif len(species) == 1:
            call, dom = next(iter(species)), next(iter(domains)) if domains else "unknown"
        else:
            call, dom = "unknown", "unknown"
        out[base] = {
            "seq_len": seq_len, "domain": dom, "species": call,
            "chimera": call.startswith("possible"),
            "species_list": sorted(species),
            "taxid": UNKNOWN_TAXID if call in ("unknown",) or call.startswith("possible")
                     else (taxids[0] if taxids else UNKNOWN_TAXID),
            "n_segments": len(items),
        }
    return out


# ---------------------------------------------------------------------------
# Consensus + decision
# ---------------------------------------------------------------------------
def compare(fcs: Dict[str, dict], tiara: Dict[str, str], min_len: int = 1000) -> Dict[str, dict]:
    """comparison.R Steps 5-6: match/mismatch/chimera per contig, small
    contigs -> unknown. Contigs known only to Tiara are included with FCS
    domain 'unknown' (the R left_join would drop them; EGAP must not)."""
    ids = set(fcs) | set(tiara)
    out = {}
    for cid in sorted(ids):
        f = fcs.get(cid, {"seq_len": None, "domain": "unknown", "species": "unknown",
                          "chimera": False, "species_list": [], "taxid": UNKNOWN_TAXID,
                          "n_segments": 0})
        t_dom = tiara_to_domain(tiara.get(cid))
        if f["chimera"]:
            match = "chimera"
        elif t_dom == f["domain"]:
            match = "match"
        else:
            match = "mismatch"
        seq_len = f["seq_len"]
        if seq_len is not None and seq_len < min_len:
            tag = "unknown"
        elif match == "match":
            tag = f["species"].replace(" ", "_")
        elif match == "chimera":
            tag = "possible_chimera"
        else:
            tag = f"{t_dom}_{f['domain']}"
        out[cid] = {"tiara_domain": t_dom, "fcs_domain": f["domain"], "match": match,
                    "blob_tag": tag, "seq_len": seq_len, "chimera": f["chimera"],
                    "species": f["species"], "taxid": f["taxid"]}
    return out


def decide(kingdom: str, tiara_domain: str, fcs_domain: str, chimera: bool) -> str:
    """Return 'keep' | 'remove' | 'flag'. See module docstring for the rule."""
    target = TARGET_DOMAIN.get((kingdom or "").strip().lower(), "eukarya")
    if chimera:
        return "flag"
    t_non = tiara_domain not in ("unknown", target)
    f_non = fcs_domain not in ("unknown", target)
    # Tiara 'prokarya' agrees with FCS bacteria/archaea for a eukaryote target
    if target == "eukarya" and tiara_domain == "prokarya" and fcs_domain in ("bacteria", "archaea"):
        return "remove"
    if t_non and f_non and tiara_domain == fcs_domain:
        return "remove"
    if t_non or f_non:
        return "flag"
    return "keep"


def partition(kingdom: str, comparison: Dict[str, dict]) -> Dict[str, List[str]]:
    buckets = {"keep": [], "remove": [], "flag": []}
    for cid, c in comparison.items():
        buckets[decide(kingdom, c["tiara_domain"], c["fcs_domain"], c["chimera"])].append(cid)
    return buckets


def write_review_table(path, comparison: Dict[str, dict], decisions: Dict[str, List[str]]):
    inv = {cid: d for d, cids in decisions.items() for cid in cids}
    with open(path, "w") as fh:
        fh.write("seq_id\tseq_len\ttiara_domain\tfcs_domain\tmatch\tspecies\tdecision\n")
        for cid, c in comparison.items():
            fh.write(f"{cid}\t{c['seq_len']}\t{c['tiara_domain']}\t{c['fcs_domain']}\t"
                     f"{c['match']}\t{c['species']}\t{inv[cid]}\n")
    return Path(path)
