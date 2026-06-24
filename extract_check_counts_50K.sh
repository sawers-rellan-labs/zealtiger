#!/usr/bin/env bash
# Extract SNP50K allele counts for the CHECK samples (B73 + Purple) that
# make_rtiger_50K_input.sh deliberately skips (they have donor==NA).
#
# Source: /Volumes/BZea/bzeaseq/50K/allelic_counts50K.tsv (2.5 GB, sample-grouped,
#   cols: SAMPLE CONTIG POSITION REF_COUNT ALT_COUNT REF_NUCLEOTIDE ALT_NUCLEOTIDE).
# One streaming pass; writes RTIGER 6-col (Chr Pos RefBase RefCount AltBase AltCount)
# to data/rtiger_50K/counts_checks/<group>/<sample>.tsv and an expDesign_checks.csv.
set -euo pipefail

IN=${IN:-/Volumes/BZea/bzeaseq/50K/allelic_counts50K.tsv}
LIST=${LIST:-/tmp/skim_checks.tsv}          # sample<TAB>group
OUT=${OUT:-data/rtiger_50K}
DEST="$OUT/counts_checks"

[ -r "$IN" ]   || { echo "ERROR: cannot read $IN" >&2; exit 1; }
[ -r "$LIST" ] || { echo "ERROR: cannot read $LIST" >&2; exit 1; }

rm -rf "$DEST"; mkdir -p "$DEST"
echo "files,name,group" > "$OUT/expDesign_checks.csv"

awk -F'\t' -v DEST="$DEST" -v ED="$OUT/expDesign_checks.csv" '
  FNR==NR { grp[$1]=$2; want[$1]=1; next }   # pass 1: wanted sample -> group
  FNR==1  { next }                            # pass 2 header
  {
    s=$1
    if (!(s in want)) next
    if (s != cur) {
      if (cur != "") { close(curfile); done[cur]=1 }
      if (s in done) { print "ERROR: not grouped by sample ("s")"; exit 3 }
      cur=s; g=grp[s]; dir=DEST"/"g
      system("mkdir -p \"" dir "\"")
      curfile=dir"/"s".tsv"
      print curfile","s","g >> ED
      nsamp++
    }
    ref=$6; rc=$4; alt=$7; ac=$5
    if (alt=="N") alt=(ref=="A" ? "C" : "A")
    print $2"\t"$3"\t"ref"\t"rc"\t"alt"\t"ac >> curfile
  }
  END { printf("Wrote %d check sample files.\n", nsamp+0) > "/dev/stderr" }
' "$LIST" "$IN"

echo "Per-sample marker counts (should be ~52K each):"
for f in "$DEST"/*/*.tsv; do printf "  %s\t%s\n" "$(basename "$f" .tsv)" "$(wc -l < "$f")"; done
echo "Catalogue: $OUT/expDesign_checks.csv"