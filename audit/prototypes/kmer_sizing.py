"""Data-driven k-mer sizing (port candidate from fspassemblypipeline).

fsp has three strategies run as separate assembly branches:
  manual      the assembler default k list
  kmergenie   k parsed from KmerGenie's HTML "<h2>Predicted best k: N</h2>"
  read_length k = int(median_read_length * 2 / 3), forced odd (GETSEQKITK)
The predicted k is validated (odd, 15 <= k <= assembler max, not already in
the default list) and appended to the default list.

EGAP today: assemble_spades.py hard-codes ["21","33","55","77","99"] and
removes the matching K* directories on cleanup, so a k-list change touches
both places.

Everything below is pure Python and needs no tools. The FASTQ reader is a
streaming median so it works on 100M-read files without loading them.
"""
import gzip
import re
import statistics
from pathlib import Path
from typing import Iterable, List, Optional, Sequence

SPADES_DEFAULT_K: List[int] = [21, 33, 55, 77]      # fsp's SPAdes default
EGAP_SPADES_K: List[int] = [21, 33, 55, 77, 99]     # EGAP's current hard-coded list
SPADES_MAX_K = 127
MEGAHIT_MAX_K = 141
MIN_K = 15

_KMERGENIE_RE = re.compile(r"Predicted best k:\s*(\d+)")


def read_lengths(fastq_path, limit: Optional[int] = None) -> Iterable[int]:
    """Yield sequence lengths from a FASTQ (.gz ok). ``limit`` caps reads."""
    p = Path(fastq_path)
    opener = gzip.open if p.suffix == ".gz" else open
    n = 0
    with opener(p, "rt") as fh:
        while True:
            header = fh.readline()
            if not header:
                return
            seq = fh.readline().rstrip("\n")
            fh.readline()
            fh.readline()
            yield len(seq)
            n += 1
            if limit and n >= limit:
                return


def median_read_length(fastq_path, limit: Optional[int] = 200_000) -> int:
    lens = list(read_lengths(fastq_path, limit))
    if not lens:
        raise ValueError(f"no reads in {fastq_path}")
    return int(statistics.median(lens))


def force_odd(k: int) -> int:
    """fsp's awk: (k % 2 == 0) ? k - 1 : k."""
    return k - 1 if k % 2 == 0 else k


def heuristic_k(median_len: int) -> int:
    """k = 2/3 of median read length, odd."""
    return force_odd(int(median_len * 2 / 3))


def parse_kmergenie_html(text: str) -> Optional[int]:
    m = _KMERGENIE_RE.search(text)
    return int(m.group(1)) if m else None


def validate_k(k: Optional[int], max_k: int = SPADES_MAX_K, min_k: int = MIN_K,
               read_len: Optional[int] = None) -> bool:
    """Odd, inside [min_k, max_k], and strictly below read length when known
    (SPAdes refuses k >= read length)."""
    if k is None or k % 2 == 0 or k < min_k or k > max_k:
        return False
    if read_len is not None and k >= read_len:
        return False
    return True


def build_k_list(default: Sequence[int], predicted: Optional[int],
                 max_k: int = SPADES_MAX_K, read_len: Optional[int] = None) -> List[int]:
    """fsp rule: append predicted k to the defaults if valid and absent, then
    sort. Also drops any default k that is not below read_len, which is the
    failure EGAP would hit with a 99-mer on reads cropped under 100 bp."""
    ks = [k for k in default if validate_k(k, max_k, read_len=read_len)]
    if validate_k(predicted, max_k, read_len=read_len) and predicted not in ks:
        ks.append(predicted)
    return sorted(set(ks))


def spades_k_arg(ks: Sequence[int]) -> str:
    return ",".join(str(k) for k in ks)


def spades_k_dirs(ks: Sequence[int]) -> List[str]:
    """The per-kmer working directories assemble_spades.py deletes on cleanup
    must be derived from the same list, not a second hard-coded copy."""
    return [f"K{k}" for k in ks]
