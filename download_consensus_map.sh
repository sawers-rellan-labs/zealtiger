#!/usr/bin/env bash
# Download the MaizeGDB consensus genetic-map files (chr1-10) — the input for
# Stage 1 (map compilation). Run ONCE from the project root; it populates
#   data/maizegdb_consensus_map/consensus_chr01.txt … consensus_chr10.txt
#
# Each file is a MaizeGDB "map_text" export (consensus genetic map with the
# per-assembly physical anchors, incl. the B73 NAM v5 columns we parse).
# Map ids are 1203637 (chr1) … 1203646 (chr10).
set -euo pipefail

outdir="data/maizegdb_consensus_map"
mkdir -p "$outdir"

for i in {1..10}; do
  id_suffix=$((36 + i))                                  # 37..46  -> id 1203637..1203646
  filename=$(printf "consensus_chr%02d.txt" "$i")
  curl -s -o "$outdir/$filename" \
    "https://www.maizegdb.org/map_text?id=12036${id_suffix}"
done

echo "Downloaded $(ls "$outdir"/consensus_chr*.txt | wc -l | tr -d ' ') files to $outdir/"
