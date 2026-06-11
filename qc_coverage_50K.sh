#!/usr/bin/env bash
# Per-sample coverage QC for the SNP50K RTIGER input.
#
# RTIGER requires every chromosome of every sample to carry at least 2*rigidity
# COVERED (REF+ALT>0) markers; a failed/empty library (e.g. one chromosome with
# <6 covered sites at r=3) aborts the entire joint fit it sits in. This scans
# every per-sample count file and records coverage so the fitter can drop the
# failures.
#
# Output: agent/coverage_qc.csv  (sample,donor,n_chr_cov,min_chr_cov,total_cov)
# "pass" is decided downstream against the chosen rigidity (min_chr_cov >= 2r
# AND n_chr_cov == 10).
#
# Usage:  bash qc_coverage_50K.sh
set -euo pipefail
cd "$(dirname "$0")"

ED=data/rtiger_50K/expDesign_all.csv
OUT=agent/coverage_qc.csv
echo "sample,donor,n_chr_cov,min_chr_cov,total_cov" > "$OUT"

tail -n +2 "$ED" | while IFS=, read -r file name donor; do
  awk -F'\t' -v name="$name" -v donor="$donor" '
    ($4+$6)>0 { c[$1]++ }
    END{
      m=1e9; tot=0; nchr=0
      for (k in c){ nchr++; tot+=c[k]; if(c[k]<m) m=c[k] }
      if (nchr==0) m=0
      print name","donor","nchr","m","tot
    }' "$file" >> "$OUT"
done

echo "Wrote $OUT  ($(($(wc -l < "$OUT")-1)) samples)"
echo
echo "=== failed libraries at r=3 (min_chr_cov < 6 or n_chr_cov < 10) ==="
awk -F, 'NR>1 && ($3<10 || $4<6){print "  ",$0}' "$OUT"
n_fail=$(awk -F, 'NR>1 && ($3<10 || $4<6){c++} END{print c+0}' "$OUT")
echo "Total failed: $n_fail / $(($(wc -l < "$OUT")-1))"
