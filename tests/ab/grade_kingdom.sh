#!/bin/bash
# Grade each kingdom benchmark sample's final assembly against its reference.
# usage: grade_kingdom.sh <out dir (out/kingdom)> <reference dir> [image]
# Reference files (from NCBI): yeast.fasta GCF_000146045.2, mjannaschii.fasta
# GCF_000091665.1, arabidopsis.fasta GCA_028009825.2 (Col-CC),
# celegans.fasta GCA_052911875.1 (CGC1 T2T).
set -uo pipefail
OUT=$(cd "$1" && pwd) REF=$(cd "$2" && pwd) IMG=${3:-localhost:5050/egap:pd3}
w() { cygpath -w "$1" 2>/dev/null || echo "$1"; }
ref_for() {
  case $1 in
    Saccharomyces*) echo yeast ;; Methanocaldococcus*) echo mjannaschii ;;
    Arabidopsis*) echo arabidopsis ;; Caenorhabditis*) echo celegans ;;
  esac
}
for summary in "$OUT"/*_ab_summary.txt; do
  [ -e "$summary" ] || continue
  sample=$(basename "$summary" _ab_summary.txt)
  asm=$(find "$OUT" -path "*/$sample/${sample}_final_EGAP_assembly.fasta*" | head -1)
  stats=$(find "$OUT" -path "*/$sample/${sample}_final_EGAP_assembly_final_stats.csv" | head -1)
  ref=$(ref_for "$sample")
  [ -n "$asm" ] && [ -n "$ref" ] || { echo "$sample: no final assembly or reference"; continue; }
  large=""; [ "$ref" = arabidopsis ] || [ "$ref" = celegans ] && large="--large"
  rel=${asm#"$OUT"/}
  MSYS_NO_PATHCONV=1 docker run --rm --cpus 14 -v "$(w "$OUT")":/o:ro -v "$(w "$REF")":/r:ro --entrypoint bash "$IMG" -c "
    export PATH=/EGAP_env/bin:\$PATH
    a=/o/$rel; case \$a in *.gz) gzip -dc \$a > /tmp/a.fasta; a=/tmp/a.fasta;; esac
    quast.py -r /r/$ref.fasta -t 14 $large --min-contig 1000 -o /tmp/q \$a > /dev/null 2>&1
    grep -E '^(# contigs|Total length|N50|# misassemblies|Genome fraction|# mismatches per 100 kbp|# indels per 100 kbp|Largest alignment)\b' /tmp/q/report.txt | sed 's/  */ /g'
  " > "$OUT/${sample}_grade.txt"
  python3 - "$sample" "$summary" "$stats" "$OUT/${sample}_grade.txt" <<'PY'
import csv, sys
sample, summary, stats, grade = sys.argv[1:5]
s = dict(kv.split("=") for kv in open(summary).read().split()[1:])
r = next(csv.DictReader(open(stats))) if stats else {}
q = {l.rsplit(" ", 1)[0]: l.rsplit(" ", 1)[1] for l in open(grade).read().splitlines() if l.strip()}
print(f"{sample}: wall {int(s['wall_s'])/3600:.1f} h, peak {int(s['peak_mem_bytes'])/2**30:.1f} GiB, "
      f"BUSCO {r.get('FIRST_BUSCO_C','?')}/{r.get('SECOND_BUSCO_C','?')}, "
      f"{q.get('# contigs','?')} contigs, N50 {q.get('N50','?')}, genome fraction {q.get('Genome fraction (%)','?')}%, "
      f"misassemblies {q.get('# misassemblies','?')}, mismatches {q.get('# mismatches per 100 kbp','?')} + "
      f"indels {q.get('# indels per 100 kbp','?')} /100 kbp")
PY
done
