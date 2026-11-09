#!/usr/bin/env bash
set -euo pipefail

slides=(04 09 14 19 24)

run_script="./Reference_based/merfish_run_methods_by_slide.R"
plot_script="./Reference_based/merfish_matrix_comparison_by_slide.R"

out_dir="./Reference_based"
combined_pdf="${out_dir}/MERFISH_all_slides_matrix_comparison.pdf"

pdfs=()
for slide in "${slides[@]}"; do
  echo "[INFO] Running methods for slide ${slide}"
  Rscript "${run_script}" "${slide}"
  echo "[INFO] Plotting matrix comparison for slide ${slide}"
  Rscript "${plot_script}" "${slide}"
  echo "[INFO] Done slide ${slide}"
  echo
  
  # Store expected output path for later
  pdfs+=("${out_dir}/MERFISH_X0_${slide}_matrix_comparison.pdf")
done

if command -v pdfunite >/dev/null 2>&1; then
  echo "[INFO] Combining per-slide PDFs into ${combined_pdf}"
  pdfunite "${pdfs[@]}" "${combined_pdf}"
  echo "[INFO] Combined PDF written to ${combined_pdf}"
else
  echo "[WARN] pdfunite not found; skipping combined PDF."
  echo "[WARN] Per-slide PDFs are in ${out_dir}"
fi
