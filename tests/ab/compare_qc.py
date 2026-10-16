"""Side-by-side QC columns of each sample's final stats CSV in two A/B output trees.

usage: python compare_qc.py <out/variantA> <out/variantB>
"""
import csv
import sys
from pathlib import Path

COLS = ("FIRST_COMPLEASM_C", "SECOND_COMPLEASM_C", "FIRST_BUSCO_C", "SECOND_BUSCO_C",
        "FIRST_BUSCO_D", "FIRST_BUSCO_M", "ASSEMBLY_CONTIGS", "ASSEMBLY_N50", "GENOME_SIZE")


def row(p):
    with open(p, newline="") as f:
        return next(csv.DictReader(f))


def main(a, b):
    a, b = Path(a), Path(b)
    for fa in sorted(a.rglob("*_final_EGAP_assembly_final_stats.csv")):
        fb = b / fa.relative_to(a)
        if not fb.exists():
            continue
        ra, rb = row(fa), row(fb)
        print(ra["SAMPLE_ID"])
        for c in COLS:
            va, vb = ra.get(c, ""), rb.get(c, "")
            print(f"    {c:<20} {va:>12} {vb:>12}{'' if va == vb else '   <- differs'}")


if __name__ == "__main__":
    main(*sys.argv[1:3])
