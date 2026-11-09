#!/bin/bash

# Batch script to run hdWGCNA analysis for all scenarios
# S1-S6 with suffixes _1 to _5 (30 total analyses)

# Set up environment
module load r
export R_LIBS_USER="${R_LIBS_USER:-$HOME/R/library}"

# Change to working directory
cd "${STCOS_ANALYSIS_ROOT:-$(pwd)}"

# Create log directory
mkdir -p STSimulator/coexpression/logs

# Function to run hdWGCNA analysis with error handling
run_hdwgcna() {
    local S=$1
    local suffix=$2
    local log_file="STSimulator/coexpression/logs/hdWGCNA_${S}${suffix}.log"
    
    echo "Starting hdWGCNA analysis for ${S}${suffix}..."
    echo "Log file: $log_file"
    
    # Run the analysis and capture both stdout and stderr
    if Rscript STSimulator/coexpression/hdWGCNA.R "$S" "$suffix" > "$log_file" 2>&1; then
        echo "✓ Completed successfully: ${S}${suffix}"
        echo "COMPLETED: ${S}${suffix}" >> "$log_file"
    else
        echo "✗ Failed: ${S}${suffix}"
        echo "FAILED: ${S}${suffix}" >> "$log_file"
    fi
}

# Array of scenarios and suffixes
scenarios=("S1" "S2" "S3" "S4" "S5" "S6")
suffixes=("_1" "_2" "_3" "_4" "_5")

echo "Starting batch hdWGCNA analysis for all scenarios..."
echo "Total analyses to run: $((${#scenarios[@]} * ${#suffixes[@]}))"
echo "Start time: $(date)"

# Initialize counters
total_analyses=0
completed_analyses=0
failed_analyses=0

# Run analysis for each combination
for S in "${scenarios[@]}"; do
    for suffix in "${suffixes[@]}"; do
        total_analyses=$((total_analyses + 1))
        
        echo ""
        echo "=== Analysis $total_analyses/$((${#scenarios[@]} * ${#suffixes[@]})): ${S}${suffix} ==="
        
        # Check if input files exist before running
        gene_file="./Data/${S}_gene/simulated/${S}_gene_expression${suffix}.csv"
        coord_file="./Data/${S}_gene/simulated/${S}_coord${suffix}.csv"
        layer_file="./Data/${S}_gene/simulated/${S}_layer${suffix}.csv"
        
        if [[ -f "$gene_file" && -f "$coord_file" && -f "$layer_file" ]]; then
            run_hdwgcna "$S" "$suffix"
            
            # Check if analysis completed successfully by looking for the success message
            if grep -q "COMPLETED" "STSimulator/coexpression/logs/hdWGCNA_${S}${suffix}.log"; then
                completed_analyses=$((completed_analyses + 1))
            else
                failed_analyses=$((failed_analyses + 1))
            fi
        else
            echo "✗ Skipping ${S}${suffix}: Input files not found"
            echo "Missing files:"
            [[ ! -f "$gene_file" ]] && echo "  - $gene_file"
            [[ ! -f "$coord_file" ]] && echo "  - $coord_file"
            [[ ! -f "$layer_file" ]] && echo "  - $layer_file"
            failed_analyses=$((failed_analyses + 1))
        fi
    done
done

echo ""
echo "=== Batch Analysis Summary ==="
echo "End time: $(date)"
echo "Total analyses: $total_analyses"
echo "Completed successfully: $completed_analyses"
echo "Failed: $failed_analyses"
echo "Success rate: $(( completed_analyses * 100 / total_analyses ))%"

echo ""
echo "All results saved to: ./coexpression/hdWGCNA_Agent/"
echo "Logs saved to: ./coexpression/logs/"

# Create a summary report
summary_file="STSimulator/coexpression/hdWGCNA_batch_summary.txt"
echo "Batch hdWGCNA Analysis Summary" > "$summary_file"
echo "Generated on: $(date)" >> "$summary_file"
echo "" >> "$summary_file"
echo "Total analyses: $total_analyses" >> "$summary_file"
echo "Completed successfully: $completed_analyses" >> "$summary_file"
echo "Failed: $failed_analyses" >> "$summary_file"
echo "Success rate: $(( completed_analyses * 100 / total_analyses ))%" >> "$summary_file"
echo "" >> "$summary_file"

echo "Results summary saved to: $summary_file"