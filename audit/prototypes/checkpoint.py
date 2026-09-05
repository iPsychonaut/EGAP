"""Declared-outputs checkpoint markers (Matthew's "checkpoint integrity as code").

Current EGAP pattern (e.g. decontaminate_assembly.py):
    if os.path.exists(done_marker) and os.path.exists(final_out): SKIP
    if os.path.exists(tiara_out): SKIP   # no size check
The marker is free text and vouches for nothing. A run killed mid-write
leaves a truncated tiara_out that is trusted forever on --resume.

Replacement contract:
    write_done_marker(marker, outputs)  records path, size, sha256 (optional),
                                        mtime for every declared output
    step_can_skip(marker, outputs)      True only if the marker parses, every
                                        declared output still exists, is
                                        non-empty, and matches recorded size
                                        (and hash when recorded)
Anything else -> False, and the caller re-runs the step. Cheap, no deps,
and the JSON marker doubles as provenance the HTML reporter can show.
"""
import hashlib
import json
import os
import time
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple

MARKER_VERSION = 1


def _sha256(path: Path, chunk: int = 1 << 20) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(chunk), b""):
            h.update(block)
    return h.hexdigest()


def write_done_marker(marker_path, outputs: Iterable, step: str = "",
                      hash_files: bool = False, extra: Optional[dict] = None,
                      min_size: int = 1) -> Path:
    """Refuses to write if any declared output is missing or below min_size.
    That refusal is the whole point: a marker must never precede its outputs."""
    records = []
    for o in outputs:
        p = Path(o)
        if not p.exists():
            raise FileNotFoundError(f"declared output missing, marker not written: {p}")
        size = p.stat().st_size
        if size < min_size:
            raise ValueError(f"declared output is empty ({size} B), marker not written: {p}")
        rec = {"path": str(p.resolve()), "size": size, "mtime": p.stat().st_mtime}
        if hash_files:
            rec["sha256"] = _sha256(p)
        records.append(rec)
    doc = {"version": MARKER_VERSION, "step": step, "written": time.time(),
           "outputs": records, "extra": extra or {}}
    marker_path = Path(marker_path)
    tmp = marker_path.with_suffix(marker_path.suffix + ".tmp")
    tmp.write_text(json.dumps(doc, indent=2))
    os.replace(tmp, marker_path)   # atomic on POSIX
    return marker_path


def read_done_marker(marker_path) -> Optional[dict]:
    p = Path(marker_path)
    if not p.exists():
        return None
    try:
        doc = json.loads(p.read_text())
    except (json.JSONDecodeError, UnicodeDecodeError):
        return None   # legacy free-text marker or corruption -> not trusted
    if not isinstance(doc, dict) or doc.get("version") != MARKER_VERSION:
        return None
    return doc


def verify_outputs(doc: dict, expected: Optional[Iterable] = None,
                   check_hash: bool = True) -> Tuple[bool, List[str]]:
    """Return (ok, reasons). ``expected`` lets the caller assert that the set of
    declared outputs matches what this code version expects, so a marker from
    an older EGAP that declared fewer files does not satisfy a newer step."""
    reasons = []
    declared = {r["path"]: r for r in doc.get("outputs", [])}
    if expected is not None:
        want = {str(Path(e).resolve()) for e in expected}
        missing_decl = want - set(declared)
        if missing_decl:
            reasons.append(f"marker does not declare: {sorted(missing_decl)}")
    for path, rec in declared.items():
        p = Path(path)
        if not p.exists():
            reasons.append(f"missing: {path}")
            continue
        size = p.stat().st_size
        if size == 0:
            reasons.append(f"empty: {path}")
            continue
        if size != rec.get("size"):
            reasons.append(f"size changed: {path} ({rec.get('size')} -> {size})")
            continue
        if check_hash and "sha256" in rec and _sha256(p) != rec["sha256"]:
            reasons.append(f"hash changed: {path}")
    return (not reasons), reasons


def step_can_skip(marker_path, expected: Optional[Iterable] = None,
                  check_hash: bool = True) -> Tuple[bool, List[str]]:
    doc = read_done_marker(marker_path)
    if doc is None:
        return False, ["no valid marker"]
    return verify_outputs(doc, expected, check_hash)
