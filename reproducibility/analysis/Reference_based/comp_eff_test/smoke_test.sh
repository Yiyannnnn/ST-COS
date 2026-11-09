#!/usr/bin/env bash
set -euo pipefail

# Smoke test for computational efficiency scripts (R methods only)
# - Runs comprehensive_efficiency_test_fixed.R with 50 and 100 genes
# - Verifies expected output files exist and are non-empty

cd "$(dirname "$0")"
export GENE_COUNTS=50,100
# Prefer internal RSS sampler (avoid peakRAM) for smoke test
export USE_PEAKRAM=0

echo "[smoke] Running comprehensive_efficiency_test_fixed.R with GENE_COUNTS=${GENE_COUNTS}"
Rscript comprehensive_efficiency_test_fixed.R

echo "[smoke] Verifying outputs..."
ok=1
for f in \
  fixed_comprehensive_results.csv \
  fixed_method_summary_stats.csv \
  fixed_mem_time_summary.csv
do
  if [[ -s "$f" ]]; then
    echo "[smoke] OK: $f"
  else
    echo "[smoke] MISSING or EMPTY: $f" >&2
    ok=0
  fi
done

if [[ "$ok" -ne 1 ]]; then
  echo "[smoke] One or more expected outputs missing."
  exit 1
fi

echo "[smoke] Preview summaries:"
echo "--- fixed_mem_time_summary.csv ---"
head -n 10 fixed_mem_time_summary.csv || true
echo "--- fixed_method_summary_stats.csv ---"
head -n 10 fixed_method_summary_stats.csv || true

echo "[smoke] Done."
