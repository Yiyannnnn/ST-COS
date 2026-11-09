#!/usr/bin/env bash
set -euo pipefail

# Full run for computational efficiency test (R methods only)
# - Runs comprehensive_efficiency_test_fixed.R with full gene set
# - Uses internal RSS sampler for memory to avoid peakRAM dependency

cd "$(dirname "$0")"
export GENE_COUNTS=50,100,500,1000,3000
export USE_PEAKRAM=0

echo "[full] Running comprehensive_efficiency_test_fixed.R with GENE_COUNTS=${GENE_COUNTS} (USE_PEAKRAM=${USE_PEAKRAM})"
Rscript comprehensive_efficiency_test_fixed.R

echo "[full] Verifying outputs..."
ok=1
for f in \
  fixed_comprehensive_results.csv \
  fixed_method_summary_stats.csv \
  fixed_mem_time_summary.csv \
  fixed_performance_by_genes.csv
do
  if [[ -s "$f" ]]; then
    echo "[full] OK: $f"
  else
    echo "[full] MISSING or EMPTY: $f" >&2
    ok=0
  fi
done

if [[ "$ok" -ne 1 ]]; then
  echo "[full] One or more expected outputs missing."
  exit 1
fi

echo "[full] Preview summaries:"
echo "--- fixed_mem_time_summary.csv ---"
head -n 20 fixed_mem_time_summary.csv || true
echo "--- fixed_method_summary_stats.csv ---"
head -n 20 fixed_method_summary_stats.csv || true

echo "[full] Done."

