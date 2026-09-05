"""Deterministic mock-data generators for the EGAP test suite.

Nothing here calls an external tool. Every generator is seeded so the same
inputs reproduce byte-identical files on GitHub Actions and on a laptop.

Formats produced:
    genome FASTA          random ACGT contigs, optional planted contaminant
    paired FASTQ          fixed-length or variable-length reads sampled from a genome
    Tiara classification  sequence_id <TAB> class_fst_stage <TAB> class_snd_stage ...
    FCS-GX taxonomy.rpt   header-driven layout; long contigs split into "~start-end" segments
    BUSCO short_summary   the text block EGAP/fsp both grep for
    QUAST report.tsv      transposed metric rows, one column per assembly
    MerquryFK .qv / .completeness.stats
"""
import gzip
import random
from pathlib import Path

BASES = "ACGT"


def random_seq(rng, length, gc=0.45):
    """Random nucleotide string with approximate GC fraction."""
    weights = [(1 - gc) / 2, gc / 2, gc / 2, (1 - gc) / 2]
    return "".join(rng.choices(BASES, weights=weights, k=length))


def write_fasta(path, records, width=80):
    """records: iterable of (id, seq). Writes plain or gzipped by suffix."""
    path = Path(path)
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "wt") as fh:
        for rid, seq in records:
            fh.write(f">{rid}\n")
            for i in range(0, len(seq), width):
                fh.write(seq[i:i + width] + "\n")
    return path


def make_genome(rng, n_contigs=5, min_len=2000, max_len=20000, gc=0.45,
                contaminants=0, contam_gc=0.65, contam_len=3000):
    """Return list of (id, seq). Contaminant contigs are high-GC and named 'contam_N'."""
    recs = []
    for i in range(1, n_contigs + 1):
        recs.append((f"contig_{i}", random_seq(rng, rng.randint(min_len, max_len), gc)))
    for j in range(1, contaminants + 1):
        recs.append((f"contam_{j}", random_seq(rng, contam_len, contam_gc)))
    return recs


def _revcomp(seq):
    return seq.translate(str.maketrans("ACGT", "TGCA"))[::-1]


def write_paired_fastq(rng, genome_records, r1_path, r2_path, n_pairs=500,
                       read_len=150, insert=350, variable=False, qual_char="I"):
    """Sample read pairs uniformly from the genome. ``variable`` draws
    read lengths from [read_len*0.5, read_len] so median-length code paths
    get exercised."""
    seqs = [s for _, s in genome_records if len(s) > insert + 10]
    r1_path, r2_path = Path(r1_path), Path(r2_path)
    o1 = gzip.open if r1_path.suffix == ".gz" else open
    o2 = gzip.open if r2_path.suffix == ".gz" else open
    with o1(r1_path, "wt") as f1, o2(r2_path, "wt") as f2:
        for n in range(n_pairs):
            s = rng.choice(seqs)
            start = rng.randint(0, len(s) - insert)
            frag = s[start:start + insert]
            L = rng.randint(read_len // 2, read_len) if variable else read_len
            r1 = frag[:L]
            r2 = _revcomp(frag)[:L]
            f1.write(f"@read_{n}/1\n{r1}\n+\n{qual_char * len(r1)}\n")
            f2.write(f"@read_{n}/2\n{r2}\n+\n{qual_char * len(r2)}\n")
    return r1_path, r2_path


def write_tiara(path, labels):
    """labels: dict seq_id -> class (eukarya/bacteria/archaea/organelle/prokarya/unknown)."""
    with open(path, "w") as fh:
        fh.write("sequence_id\tclass_fst_stage\tclass_snd_stage\n")
        for sid, lab in labels.items():
            snd = "mitochondrion" if lab == "organelle" else "n/a"
            fh.write(f"{sid}\t{lab}\t{snd}\n")
    return Path(path)


FCSGX_HEADER = ["seq-id", "seq-len", "masked-len", "best-taxid", "best-tax-name",
                "div", "cvg-by-div", "cvg-by-tax", "score",
                "tax-id-1", "tax-name-1", "div-1", "cvg-by-div-1", "cvg-by-tax-1", "score-1"]


def write_fcsgx_rpt(path, rows):
    """rows: list of dicts with keys seq_id, seq_len, taxid, tax_name, div.
    A row may carry 'segment': (start, end) to emulate FCS-GX splitting a
    long contig into windows, which is what the chimera logic keys on."""
    with open(path, "w") as fh:
        fh.write("##[[\"FCS genome report\", 2, 1], {\"git-rev\": \"mock\", \"db\": {\"build-date\": \"2024-01-01\"}}]\n")
        fh.write("#" + "\t".join(FCSGX_HEADER) + "\n")
        for r in rows:
            sid = r["seq_id"]
            if "segment" in r:
                sid = f"{sid}~{r['segment'][0]}-{r['segment'][1]}"
            vals = [sid, str(r["seq_len"]), "0", str(r["taxid"]), r["tax_name"],
                    r["div"], "90", "85", "1000",
                    str(r["taxid"]), r["tax_name"], r["div"], "90", "85", "1000"]
            fh.write("\t".join(vals) + "\n")
    return Path(path)


def write_busco_summary(path, lineage, complete=95.0, single=90.0, dup=5.0,
                        frag=2.0, missing=3.0, n=758):
    body = f"""# BUSCO version is: 5.8.2
# The lineage dataset is: {lineage} (Creation date: 2024-01-08, number of genomes: 100, number of BUSCOs: {n})
# Summarized benchmarking in BUSCO notation for file mock.fasta
# BUSCO was run in mode: euk_genome_min
# Gene predictor used: miniprot

	***** Results: *****

	C:{complete:.1f}%[S:{single:.1f}%,D:{dup:.1f}%],F:{frag:.1f}%,M:{missing:.1f}%,n:{n}
	{int(round(complete * n / 100))}	Complete BUSCOs (C)
	{int(round(single * n / 100))}	Complete and single-copy BUSCOs (S)
	{int(round(dup * n / 100))}	Complete and duplicated BUSCOs (D)
	{int(round(frag * n / 100))}	Fragmented BUSCOs (F)
	{int(round(missing * n / 100))}	Missing BUSCOs (M)
	{n}	Total BUSCO groups searched
"""
    Path(path).write_text(body)
    return Path(path)


def write_quast_tsv(path, assemblies):
    """assemblies: dict name -> dict of metric -> value (N50, auN, # contigs, Total length, L50)."""
    names = list(assemblies)
    metrics = ["# contigs", "Total length", "N50", "L50", "auN", "GC (%)"]
    with open(path, "w") as fh:
        fh.write("Assembly\t" + "\t".join(names) + "\n")
        for m in metrics:
            fh.write(m + "\t" + "\t".join(str(assemblies[n].get(m, "-")) for n in names) + "\n")
    return Path(path)


def write_merquryfk(prefix, qv, completeness, error_pct=None, kmers_total=1_000_000):
    """Write <prefix>.qv and <prefix>.completeness.stats in MerquryFK layout."""
    prefix = Path(prefix)
    if error_pct is None:
        error_pct = 100 * (10 ** (-qv / 10))
    no_support = int(kmers_total * error_pct / 100)
    (prefix.parent / (prefix.name + ".qv")).write_text(
        "Assembly\tNo Support\tTotal\tError %\tQV\n"
        f"{prefix.name}\t{no_support}\t{kmers_total}\t{error_pct:.4f}\t{qv:.1f}\n")
    found = int(kmers_total * completeness / 100)
    (prefix.parent / (prefix.name + ".completeness.stats")).write_text(
        "Assembly\tRegion\tFound\tTotal\t% Covered\n"
        f"{prefix.name}\tall\t{found}\t{kmers_total}\t{completeness:.2f}\n")
    return prefix


def write_kmergenie_html(path, best_k):
    Path(path).write_text(
        "<html><body><h1>KmerGenie report</h1>\n"
        f"<h2>Predicted best k: {best_k}</h2>\n"
        "<p>Predicted assembly size: 12345678</p></body></html>\n")
    return Path(path)


def seeded(seed=42):
    return random.Random(seed)
