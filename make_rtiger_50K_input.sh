#!/usr/bin/env bash
# Build per-donor RTIGER input from the real SNP50K allele counts (Source B).
#
# Source B = /Volumes/.../bzeaseq/50K/allelic_counts50K.tsv : one concatenated
#   file, columns  SAMPLE CONTIG POSITION REF_COUNT ALT_COUNT REF_NUCLEOTIDE
#   ALT_NUCLEOTIDE, grouped by sample. ~52K array sites/sample, B73 NAM v5.
#
# We reformat each sample to RTIGER's 6-col, NO-header, tab-separated layout:
#   Chr  Position  RefBase  RefCount  AltBase  AltCount
# (uncovered sites kept as 0/0; ALT base "N" -> a non-ref base, cosmetic only —
#  RTIGER uses counts, not bases). Files are partitioned by donor accession:
#   data/rtiger_50K/counts/<donor_id>/<sample>.tsv
# and catalogued in data/rtiger_50K/expDesign_all.csv (files,name,donor).
#
# The sample -> donor_id join comes from bzeaseq/sample_metatada.csv (project
# "bzea", donor_id != NA). Samples not in that map (checks/other projects) are
# skipped and counted.
#
# Assumes Source B is grouped by sample (the pipeline concatenated per-sample
# outputs). The splitter keeps ONE file handle open and aborts loudly if a
# sample reappears after its block closed — so a wrong assumption fails safe
# rather than silently corrupting output.
#
# Usage:  bash make_rtiger_50K_input.sh
set -euo pipefail

SHARE=/Volumes/rsstu/users/r/rrellan/BZea/bzeaseq
IN=${IN:-$SHARE/50K/allelic_counts50K.tsv}
META=${META:-$SHARE/sample_metatada.csv}
OUT=${OUT:-$(cd "$(dirname "$0")" && pwd)/data/rtiger_50K}

[ -r "$IN" ]   || { echo "ERROR: cannot read counts file: $IN"   >&2; exit 1; }
[ -r "$META" ] || { echo "ERROR: cannot read metadata file: $META" >&2; exit 1; }

echo "Source counts : $IN"
echo "Sample->donor : $META"
echo "Output root   : $OUT"

# Fresh start (>> appends, so clear any previous split).
rm -rf "$OUT/counts" "$OUT/expDesign_all.csv"
mkdir -p "$OUT/counts"

# B73 NAM v5 chromosome lengths (from the BAM/reference @SQ header used to make
# the counts). RTIGER seqlengths must cover the max marker position per chr.
cat > "$OUT/seqlengths.csv" <<'EOF'
chr_label,len
chr1,308452471
chr2,243675191
chr3,238017767
chr4,250330460
chr5,226353449
chr6,181357234
chr7,185808916
chr8,182411202
chr9,163004744
chr10,152435371
EOF

echo "files,name,donor" > "$OUT/expDesign_all.csv"

# One streaming pass. -F'[,\t]' splits BOTH files: metadata is comma-delimited
# (8 fields), counts is tab-delimited (7 fields) — no field uses the other's
# delimiter, so this is unambiguous.
awk -F'[,\t]' -v OUT="$OUT" '
  # ---- pass 1: metadata (first file) -> sample->donor map ----
  FNR==NR {
    if (FNR==1) next                       # header
    proj=$2; samp=$3; donor=$7
    if (proj=="bzea" && donor!="NA" && donor!="") map[samp]=donor
    next
  }
  # ---- pass 2: allele counts ----
  FNR==1 { next }                          # header: SAMPLE CONTIG ...
  {
    s=$1
    if (!(s in map)) { if (!(s in skipped)) { skipped[s]=1; nskip++ }; next }
    if (s != cur) {
      if (cur != "") { close(curfile); done[cur]=1 }
      if (s in done) {
        print "ERROR: counts not grouped by sample ("s" reappears). Aborting." > "/dev/stderr"
        exit 3
      }
      cur=s
      donor=map[s]
      dir=OUT"/counts/"donor
      system("mkdir -p \"" dir "\"")
      curfile=dir"/"s".tsv"
      print curfile","s","donor >> (OUT"/expDesign_all.csv")
      nsamp++; perdonor[donor]++
    }
    # reformat: Chr Pos RefBase RefCount AltBase AltCount
    ref=$6; rc=$4; alt=$7; ac=$5
    if (alt=="N") alt=(ref=="A" ? "C" : "A")   # ensure a valid, non-ref base
    print $2"\t"$3"\t"ref"\t"rc"\t"alt"\t"ac >> curfile
  }
  END {
    if (cur != "") close(curfile)
    ndon=0; for (d in perdonor) ndon++
    printf("\nWrote %d sample files across %d donor accessions.\n", nsamp, ndon) > "/dev/stderr"
    printf("Skipped %d samples not in the bzea donor map (checks/other projects).\n", nskip) > "/dev/stderr"
  }
' "$META" "$IN"

echo
echo "Donor group sizes (samples per accession):"
tail -n +2 "$OUT/expDesign_all.csv" | cut -d, -f3 | sort | uniq -c | sort -rn | head -20
echo "..."
echo "Total accessions: $(tail -n +2 "$OUT/expDesign_all.csv" | cut -d, -f3 | sort -u | wc -l | tr -d ' ')"
echo "Total sample files: $(tail -n +2 "$OUT/expDesign_all.csv" | wc -l | tr -d ' ')"
echo "Output: $OUT/{counts/<donor>/<sample>.tsv, expDesign_all.csv, seqlengths.csv}"
