#!/bin/bash
# Assembler-stage speed/memory benchmark (speed & memory A/B, PY312 plan).
# Runs one EGAP stage script (assemble_spades.py / assemble_masurca.py) on a
# fresh copy of a baseline sample's inputs, in its own egap:py312-bench
# container (7 CPU, 24 GB), with one A/B switch set. Wall time and cgroup peak
# memory cover the stage including its QC (compleasm + QUAST), as in a real run.
#
# usage: bench_stage.sh <baseline out dir> <bench dir> <tsv> <species_id> <sample_id> <stage> <case> [VAR=value]
#   stage: spades | masurca      case: a name for this run, e.g. frugal1
set -uo pipefail
BASE=$1 BENCH=$2 TSV=$3 SPECIES=$4 SAMPLE=$5 STAGE=$6 CASE=$7 ENVSET=${8:-}
IMG=localhost:5050/egap:py312-bench
w() { cygpath -w "$1" 2>/dev/null || echo "$1"; }
D=$BENCH/$SAMPLE/$STAGE-$CASE
rm -rf "${D:?}" && mkdir -p "$D/$SPECIES/Illumina" "$D/$SPECIES/$SAMPLE"
cp "$BASE/$SPECIES/Illumina/"*_dedup.fastq "$D/$SPECIES/Illumina/"
if ls "$BASE/$SPECIES/ONT/"*_highest_mean_qual_long_reads.fastq > /dev/null 2>&1; then
  mkdir -p "$D/$SPECIES/ONT" && cp "$BASE/$SPECIES/ONT/"*_highest_mean_qual_long_reads.fastq "$D/$SPECIES/ONT/"
fi
[ -d "$BASE/$SPECIES/RefSeq" ] && mkdir -p "$D/$SPECIES/RefSeq" && cp "$BASE/$SPECIES/RefSeq/"*_RefSeq.fasta* "$D/$SPECIES/RefSeq/"
awk -F'\t' -v s="$SAMPLE" 'NR==1 || $2==s' "$TSV" > "$D/sample.tsv"
MSYS_NO_PATHCONV=1 docker run --rm --cpus 7 --memory 24g -v "$(w "$D")":/w ${ENVSET:+-e "$ENVSET"} \
  --entrypoint bash "$IMG" -c "
    export PATH=/EGAP_env/bin:\$PATH CONDA_PREFIX=/EGAP_env
    ${PRE:-}
    cd /w && SECONDS=0
    /EGAP_env/bin/python /EGAP_env/bin/assemble_$STAGE.py $SAMPLE /w/sample.tsv /w 7 24 > /w/stage.log 2>&1
    rc=\$?
    echo \"BENCH sample=$SAMPLE stage=$STAGE case=$CASE env=${ENVSET:-none} wall_s=\$SECONDS peak_mem_bytes=\$(cat /sys/fs/cgroup/memory.peak) rc=\$rc\"
  " | tee -a "$BENCH/results.tsv"
# keep the assembly, its stats and the log; drop reads and work dirs (KEEP=1 keeps them)
find "$D" -name "*_${STAGE}.fasta" -exec cp {} "$D/" \;
find "$D/$SPECIES/$SAMPLE" -maxdepth 1 -name "*_${STAGE}_stats.txt" -exec cp {} "$D/" \;
[ -n "${KEEP:-}" ] || rm -rf "${D:?}/${SPECIES:?}"
