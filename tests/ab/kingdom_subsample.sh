#!/bin/bash
# Fetch one SRA run and keep a random subset sized to a target depth. The .sra
# file comes from NCBI (prefetch; ~12 MB/s here, while EBI's FASTQ mirror ran at
# ~60 KB/s), is converted with fasterq-dump straight into the sampler, and is
# deleted afterwards, so no full FASTQ is ever written. Paired runs are sampled
# a spot (both mates) at a time, which keeps mates paired.
#
# usage: kingdom_subsample.sh <run accession> <genome size bp> <target depth> <out prefix>
#   writes <out prefix>.fastq.gz (single) or <out prefix>_1/_2.fastq.gz (paired)
# needs: curl, prefetch, fasterq-dump, seqtk, pigz on PATH; TMPDIR with room for the .sra
set -euo pipefail
acc=$1 gsize=$2 depth=$3 out=$4
SEED=11
# Base count and file layout from ENA's API, which intermittently answers with
# an HTML error page; retry until the answer is "<accession>\t<bases>\t<urls>".
info=""
for try in $(seq 1 12); do
  info=$(curl -s "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=$acc&result=read_run&fields=base_count,fastq_ftp&format=tsv" | tail -n 1)
  [[ $info =~ ^$acc$'\t'[0-9]+$'\t'ftp ]] && break
  sleep 15
done
[[ $info =~ ^$acc$'\t'[0-9]+$'\t'ftp ]] || { echo "$acc: no usable ENA filereport after 12 tries" >&2; exit 1; }
bases=$(cut -f2 <<< "$info")
frac=$(awk -v g="$gsize" -v d="$depth" -v b="$bases" 'BEGIN{f=g*d/b; if (f>1) f=1; printf "%.5f", f}')
mkdir -p "$(dirname "$out")"

tmp=$(mktemp -d "${TMPDIR:-/tmp}/kingdom.XXXXXX")
# keep the converter's log next to the output if anything fails
trap 'rc=$?; [ $rc -eq 0 ] || cat "$tmp"/*.log > "$out.failed.log" 2>/dev/null; rm -rf "${tmp:?}"' EXIT
prefetch --max-size 500g -O "$tmp" "$acc" > "$tmp/prefetch.log" 2>&1
sra=$(find "$tmp/$acc" -maxdepth 1 -name "$acc.sra*" | head -1)
[ -n "$sra" ] || { echo "$acc: prefetch produced no .sra file" >&2; cat "$tmp/prefetch.log" >&2; exit 1; }

# Paired or not, from the data: ENA's file list and layout field disagree for
# some runs (SRR17374239 is listed as one file but is paired). With
# --split-spot, a paired spot is two records, <acc>.<n>/1 then <acc>.<n>/2.
paired=$(fasterq-dump --split-spot --stdout -e 2 -t "$tmp" "$sra" 2>/dev/null |
  awk 'NR==1{a=$1} NR==5{b=$1; exit} END{print (a ~ /\/1$/ && b ~ /\/2$/ && substr(a,1,length(a)-2) == substr(b,1,length(b)-2)) ? "yes" : "no"}' || true)
echo "$acc: $bases bases, keeping fraction $frac for ${depth}x of $gsize bp ($([ "$paired" = yes ] && echo paired || echo single))"

if [ "$paired" = yes ]; then
  fasterq-dump --split-spot --stdout -e 4 -t "$tmp" "$sra" 2> "$tmp/fqd.log" |
    awk -v f="$frac" -v seed="$SEED" -v o="$out" '
      BEGIN { srand(seed); c1 = "pigz -p 2 > \"" o "_1.fastq.gz.part\""; c2 = "pigz -p 2 > \"" o "_2.fastq.gz.part\"" }
      { l[(NR - 1) % 8] = $0 }
      NR % 8 == 0 && rand() < f {
        print l[0] "\n" l[1] "\n" l[2] "\n" l[3] | c1
        print l[4] "\n" l[5] "\n" l[6] "\n" l[7] | c2
      }
      END { close(c1); close(c2) }'
  for i in 1 2; do mv "${out}_$i.fastq.gz.part" "${out}_$i.fastq.gz"; done
else
  # fastq-dump streams with no scratch space; fasterq-dump --stdout failed on the
  # 57-67 Gb nanopore runs.
  fastq-dump --stdout "$sra" 2> "$tmp/fqd.log" | seqtk sample -s$SEED - "$frac" | pigz -p 2 > "$out.fastq.gz.part"
  mv "$out.fastq.gz.part" "$out.fastq.gz"
fi
echo "$acc: done"
