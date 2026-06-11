#!/usr/bin/env bash
# Sample-size scaling sweep: fit a series of REAL per-donor accession groups of
# log-spaced size, each as its own joint RTIGER fit under /usr/bin/time -l, to
# trace peak RSS + walltime vs group n. Confirms the per-donor regime (max n=100)
# is laptop-tractable. Results -> agent/sweep_results.csv.
#
# Usage:  bash sweep_rtiger_50K.sh "Zx.0190 Zx.0040 Zx.0490 Zv.0490 Zl.0040 Zd.0030 Zx.0580"
set -uo pipefail
cd "$(dirname "$0")"

# Resolve JULIA_HOME once from the juliaup launcher (version-independent) so each
# fit reuses it instead of re-detecting.
if [ -z "${JULIA_HOME:-}" ] || [ ! -x "${JULIA_HOME}/julia" ]; then
  export JULIA_HOME="$("$HOME/.juliaup/bin/julia" --startup-file=no -e 'print(Sys.BINDIR)')"
fi
echo "JULIA_HOME = $JULIA_HOME"

DONORS=${1:-"Zx.0190 Zx.0040 Zx.0490 Zv.0490 Zl.0040 Zd.0030 Zx.0580"}
OUT=agent/sweep_results.csv
LOGDIR=agent/tmp/sweep_logs
mkdir -p "$LOGDIR"
echo "donor,n,peak_rss_gb,fit_seconds,wall_seconds,status" > "$OUT"

for d in $DONORS; do
  n=$(awk -F, -v D="$d" '$3==D{c++} END{print c+0}' data/rtiger_50K/expDesign_all.csv)
  log="$LOGDIR/${d}.log"
  echo ">>> fitting $d (n=$n) ..."
  ws=$(date +%s)
  # macOS /usr/bin/time -l prints "maximum resident set size" in BYTES to stderr.
  /usr/bin/time -l Rscript fit_rtiger_by_donor.R "$d" 3 0 > "$log" 2>&1
  rc=$?
  we=$(date +%s)
  wall=$((we - ws))
  rss_bytes=$(grep -E "maximum resident set size" "$log" | awk '{print $1}' | tail -1)
  rss_gb=$(awk -v b="${rss_bytes:-0}" 'BEGIN{printf "%.2f", b/1024/1024/1024}')
  fit_s=$(grep -oE "done in [0-9.]+s" "$log" | grep -oE "[0-9.]+" | tail -1)
  [ -z "$fit_s" ] && fit_s="NA"
  status=$([ "$rc" -eq 0 ] && echo ok || echo "exit$rc")
  echo "    n=$n  peak_rss=${rss_gb}GB  fit=${fit_s}s  wall=${wall}s  [$status]"
  echo "$d,$n,$rss_gb,$fit_s,$wall,$status" >> "$OUT"
done

echo
echo "=== sweep complete -> $OUT ==="
column -t -s, "$OUT"
