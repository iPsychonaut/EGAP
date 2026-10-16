#!/bin/bash
# Fetch the kingdom benchmark read subsets (depths from the assembler audit:
# nanopore 35x, HiFi 30x, Illumina 40x eukaryotes / 50x prokaryotes).
# usage: kingdom_fetch.sh <reads dir>     (needs curl, seqtk, pigz on PATH)
set -uo pipefail
D=$1
S=$(cd "$(dirname "$0")" && pwd)/kingdom_subsample.sh
# .sra files can reach ~60 GB: keep them under the reads dir (a Windows folder
# when run from WSL), whose space is freed on delete, unlike WSL's disk image.
export TMPDIR=$D/.tmp
sub() {
  if [ -e "$D/$4.fastq.gz" ] || [ -e "$D/$4_1.fastq.gz" ]; then echo "have $4"; return; fi
  mkdir -p "$(dirname "$D/$4")"
  bash "$S" "$1" "$2" "$3" "$D/$4" > "$D/$4.log" 2>&1 && echo "ok   $4" || echo "FAIL $4 (see $D/$4.log)"
}
mkdir -p "$D" "$TMPDIR"
# genome sizes: S. cerevisiae 12.07 Mb, M. jannaschii 1.74 Mb,
# A. thaliana Col-CC 142.5 Mb (assembly), C. elegans CGC1 106.3 Mb
sub SRR18210286 12071326  30 yeast/hifi &
sub SRR17374240 12071326  35 yeast/ont &
# Short reads: SRR17374239 (BGISEQ, PRJNA792930) holds only mate 2 of each pair,
# so use ERR5526889 instead: NextSeq paired, WT of a YPH (S288C-derived) strain.
# Grading against S288C will count its few strain differences as errors.
sub ERR5526889  12071326  40 yeast/illumina &
sub ERR15459741 1739933   35 mjannaschii/ont &
wait
sub SRR14728885 142500000 30 arabidopsis/hifi &
sub SRR14739179 142500000 35 arabidopsis/ont &
sub SRR14715607 142500000 40 arabidopsis/illumina &
wait
sub DRR619744   106300000 30 celegans/hifi &
sub DRR656639   106300000 35 celegans/ont &
sub SRR7594466  106300000 40 celegans/illumina &
wait
