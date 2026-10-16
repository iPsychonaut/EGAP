#!/bin/bash
# One polishing/eval step, run inside egap:py312-bench by bench_polish.sh, one
# container per step so the cgroup's memory.peak is that step's peak.
#
# usage: bench_step.sh <step> <in.fasta> <out.fasta> <threads>
#   steps: racon2 | medaka | pilon | polypolish | pypolca | quast
# Reads are read from /in: ont.fastq, r1.fastq, r2.fastq; reference /ref/ref.fasta.
# Prints: BENCH step=<step> in=<in> wall_s=<n> peak_mem_bytes=<n> rc=<n>
set -uo pipefail
step=$1 in=$2 out=$3 t=$4
export PATH=/EGAP_env/bin:$PATH
work=$(mktemp -d /tmp/bench.XXXX)
SECONDS=0
case $step in
  racon2)  # EGAP's long-read polish: two rounds of minimap2 + racon
    cur=$in
    for i in 1 2; do
      minimap2 -t "$t" -x map-ont "$cur" /in/ont.fastq > "$work/r$i.paf" &&
      racon -t "$t" /in/ont.fastq "$work/r$i.paf" "$cur" > "$work/r$i.fasta" || break
      cur=$work/r$i.fasta
    done
    cp "$cur" "$out" ;;
  medaka)  # Medaka 2 bacterial model, as recommended for ONT R10 bacterial assemblies.
    # Named explicitly: --bacteria auto-detects from read headers, which SRA strips.
    PATH=/bench_env/bin:$PATH medaka_consensus -i /in/ont.fastq -d "$in" -o "$work/m" -t "$t" \
        -m r1041_e82_400bps_bacterial_methylation > "$work/medaka.log" 2>&1 &&
    cp "$work/m/consensus.fasta" "$out" ;;
  pilon)   # EGAP's short-read polish: bwa-mem2 + samtools sort + Pilon at 80% of 24 GB
    cp "$in" "$work/a.fasta" && bwa-mem2 index "$work/a.fasta" > /dev/null 2>&1 &&
    bwa-mem2 mem -t "$t" "$work/a.fasta" /in/r1.fastq /in/r2.fastq 2>/dev/null | samtools sort -@ "$t" -o "$work/a.bam" - &&
    samtools index "$work/a.bam" &&
    pilon -Xmx19g --genome "$work/a.fasta" --frags "$work/a.bam" --output p --outdir "$work" \
          --changes --chunksize 5000000 --fix indels,local,snps > "$work/pilon.log" 2>&1 &&
    cp "$work/p.fasta" "$out" ;;
  polypolish)  # Polypolish: every alignment per read (-a), each read file aligned separately
    cp "$in" "$work/a.fasta" && bwa index "$work/a.fasta" 2>/dev/null &&
    bwa mem -t "$t" -a "$work/a.fasta" /in/r1.fastq > "$work/1.sam" 2>/dev/null &&
    bwa mem -t "$t" -a "$work/a.fasta" /in/r2.fastq > "$work/2.sam" 2>/dev/null &&
    /bench_env/bin/polypolish filter --in1 "$work/1.sam" --in2 "$work/2.sam" --out1 "$work/f1.sam" --out2 "$work/f2.sam" 2>/dev/null &&
    /bench_env/bin/polypolish polish "$work/a.fasta" "$work/f1.sam" "$work/f2.sam" > "$out" 2>/dev/null ;;
  pypolca)
    /bench_env/bin/pypolca run -a "$in" -1 /in/r1.fastq -2 /in/r2.fastq -t "$t" -o "$work/pp" --careful > "$work/pypolca.log" 2>&1 &&
    cp "$work/pp/pypolca_corrected.fasta" "$out" ;;
  quast)   # accuracy against the same-strain reference
    quast.py -r /ref/ref.fasta -t "$t" --min-contig 1000 -o "$work/q" "$in" > /dev/null 2>&1 &&
    grep -E "^(# mismatches per 100 kbp|# indels per 100 kbp|# misassemblies|Genome fraction|Total length|# contigs)\b" "$work/q/report.txt" | sed 's/  */ /g' > "$out" ;;
  *) echo "unknown step $step" >&2; exit 2 ;;
esac
rc=$?
echo "BENCH step=$step in=$(basename "$in") wall_s=$SECONDS peak_mem_bytes=$(cat /sys/fs/cgroup/memory.peak 2>/dev/null || echo NA) rc=$rc"
rm -rf "$work"
exit $rc
