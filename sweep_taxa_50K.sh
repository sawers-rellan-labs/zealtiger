#!/usr/bin/env bash
# Resource/scaling sweep over the 5 teosinte-taxon RTIGER fits (Zh<Zl<Zd<Zv<Zx), each as
# its own joint fit at r=8 under /usr/bin/time -l, to trace peak RSS + walltime vs group
# size n. Smallest-first so the cheap points land before the Zx (~575) memory stress.
# Also produces the per-taxon calls. Results -> agent/taxa_sweep.csv.
#
# Usage:  bash sweep_taxa_50K.sh
set -uo pipefail
cd "$(dirname "$0")"

if [ -z "${JULIA_HOME:-}" ] || [ ! -x "${JULIA_HOME}/julia" ]; then
  export JULIA_HOME="$("$HOME/.juliaup/bin/julia" --startup-file=no -e 'print(Sys.BINDIR)')"
fi
export JULIA_NUM_THREADS=4     # so RTIGER threads=4 actually engages (set before startup)
echo "JULIA_HOME = $JULIA_HOME | JULIA_NUM_THREADS = $JULIA_NUM_THREADS"

TAXA=${1:-"Zh Zl Zd Zv Zx"}
OUT=agent/taxa_sweep.csv
LOGDIR=agent/tmp/taxa_logs
mkdir -p "$LOGDIR"
echo "taxon,n,peak_rss_gb,fit_seconds,wall_seconds,status" > "$OUT"

for tx in $TAXA; do
  log="$LOGDIR/${tx}.log"
  echo ">>> fitting $tx ..."
  ws=$(date +%s)
  /usr/bin/time -l Rscript fit_rtiger_by_taxa.R "$tx" 8 0 1 4 > "$log" 2>&1
  rc=$?
  we=$(date +%s); wall=$((we - ws))
  n=$(grep -oE "n=[0-9]+" "$log" | head -1 | grep -oE "[0-9]+")
  rss_bytes=$(grep -E "maximum resident set size" "$log" | awk '{print $1}' | tail -1)
  rss_gb=$(awk -v b="${rss_bytes:-0}" 'BEGIN{printf "%.2f", b/1024/1024/1024}')
  fit_s=$(grep -oE "done in [0-9.]+s" "$log" | grep -oE "[0-9.]+" | tail -1)
  [ -z "$fit_s" ] && fit_s="NA"
  status=$([ "$rc" -eq 0 ] && echo ok || echo "exit$rc")
  echo "    n=${n:-?}  peak_rss=${rss_gb}GB  fit=${fit_s}s  wall=${wall}s  [$status]"
  echo "${tx},${n:-NA},$rss_gb,$fit_s,$wall,$status" >> "$OUT"
done

echo; echo "=== taxa sweep complete -> $OUT ==="
column -t -s, "$OUT"
