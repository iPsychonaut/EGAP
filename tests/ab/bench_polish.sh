#!/bin/bash
# Polishing speed/memory/accuracy benchmark (speed & memory A/B, PY312 plan).
# Each step runs in its own egap:py312-bench container, capped like an A/B pod
# (7 CPUs, 24 GB), on the same pre-polish assembly and reads.
#
# usage: bench_polish.sh <work dir with in/{asm.fasta,ont.fastq,r1.fastq,r2.fastq}> [ref.fasta]
# Writes <work dir>/results.tsv and, when a reference is given, quast_*.txt per assembly.
set -uo pipefail
W=$(cd "$1" && pwd); REF=${2:-}
IMG=localhost:5050/egap:py312-bench
STEP=$(cd "$(dirname "$0")" && pwd)/bench_step.sh
w() { cygpath -w "$1" 2>/dev/null || echo "$1"; }
run() {  # run <step> <in name> <out name>
  MSYS_NO_PATHCONV=1 docker run --rm --cpus 7 --memory 24g \
    -v "$(w "$W/in")":/in:ro -v "$(w "$W")":/w -v "$(w "$STEP")":/bench_step.sh:ro \
    ${REF:+-v "$(w "$(dirname "$REF")")":/ref:ro} \
    --entrypoint bash "$IMG" /bench_step.sh "$1" "/w/$2" "/w/$3" 7 | tee -a "$W/results.tsv" | grep BENCH
}
[ -n "$REF" ] && [ "$(basename "$REF")" != ref.fasta ] && { echo "name the reference ref.fasta" >&2; exit 2; }
cp "$W/in/asm.fasta" "$W/prepolish.fasta"
run racon2     prepolish.fasta        racon2.fasta
run medaka     prepolish.fasta        medaka.fasta
for lr in racon2 medaka; do
  run pilon      $lr.fasta              ${lr}_pilon.fasta
  run polypolish $lr.fasta              ${lr}_polypolish.fasta
  run pypolca    ${lr}_polypolish.fasta ${lr}_polypolish_pypolca.fasta
done
if [ -n "$REF" ]; then
  for a in prepolish racon2 medaka racon2_pilon racon2_polypolish_pypolca medaka_pilon medaka_polypolish_pypolca; do
    run quast $a.fasta quast_$a.txt > /dev/null
  done
fi
