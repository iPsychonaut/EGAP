"""Compare two A/B output trees from tests/ab: each assembler's assembly.

usage: python compare_ab.py <out/variantA> <out/variantB>

For every <sample>_<assembler>.fasta[.gz] in variant A, finds the same file in
B and reports whether the sequence content is identical (gzip is undone, so a
compressed and an uncompressed copy compare equal), then the BUSCO / N50 /
contig-count lines of both <sample>_<assembler>_stats.txt files.
"""
import gzip
import hashlib
import re
import sys
from pathlib import Path

PAT = re.compile(r"_(spades|masurca|flye|hifiasm|final_curated)\.fasta(\.gz)?$")


def content_sha(p):
    opener = gzip.open if p.suffix == ".gz" else open
    h = hashlib.sha256()
    with opener(p, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def find(root, rel_no_gz):
    for cand in (root / rel_no_gz, root / (str(rel_no_gz) + ".gz")):
        if cand.exists():
            return cand
    return None


def stats(fasta):
    name = PAT.sub(r"_\1_stats.txt", fasta.name)
    for d in (fasta.parent, fasta.parent.parent):
        if (d / name).exists():
            keep = ("First BUSCO", "Second BUSCO", "N50", "Contig Count")
            return [l.strip() for l in (d / name).read_text().splitlines() if l.startswith(keep)]
    return []


def main(a, b):
    a, b = Path(a), Path(b)
    same = diff = missing = 0
    for fa in sorted(p for p in a.rglob("*") if PAT.search(p.name)):
        rel = fa.relative_to(a)
        rel_no_gz = rel.with_suffix("") if rel.suffix == ".gz" else rel
        fb = find(b, rel_no_gz)
        if fb is None:
            missing += 1
            print(f"MISSING IN B  {rel_no_gz}")
            continue
        ident = content_sha(fa) == content_sha(fb)
        same += ident
        diff += not ident
        print(f"{'IDENTICAL' if ident else 'DIFFERENT'}  {rel_no_gz}")
        sa, sb = stats(fa), stats(fb)
        if sa != sb:
            for la, lb in zip(sa, sb):
                print(f"    A: {la:<42} B: {lb}")
    print(f"\n{same} identical, {diff} different, {missing} missing in B")


if __name__ == "__main__":
    main(*sys.argv[1:3])
