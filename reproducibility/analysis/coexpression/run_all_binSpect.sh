#!/bin/bash

# Batch script to run binSpect analysis for all scenarios and replicates

echo "Starting binSpect batch analysis..."
echo "Generated on: $(date)"

# Create output directory
mkdir -p ./coexpression/binSpect_Agent

# Initialize counters
total_analyses=0
completed_successfully=0
failed=0

# Define scenarios and replicates
scenarios=("S1" "S2" "S3" "S4" "S5" "S6")
replicates=(1 2 3 4 5)

# Run binSpect for each scenario and replicate
for scenario in "${scenarios[@]}"; do
    echo "Processing scenario: $scenario"
    
    for replicate in "${replicates[@]}"; do
        echo "  Processing replicate: $replicate"
        
        total_analyses=$((total_analyses + 1))
        
        # Run binSpect analysis
        Rscript ./coexpression/binSpect.R "$scenario" "_$replicate"
        
        # Check if analysis was successful by looking for the summary file
        summary_file="./coexpression/binSpect_Agent/${scenario}_${replicate}_summary_stats.csv"
        
        if [[ -f "$summary_file" ]]; then
            # Check if processing was successful
            status=$(grep "Processing_status" "$summary_file" | cut -d',' -f2 | tr -d '"')
            if [[ "$status" == "Success" ]]; then
                completed_successfully=$((completed_successfully + 1))
                echo "    ✓ Completed successfully"
            else
                failed=$((failed + 1))
                echo "    ✗ Failed"
            fi
        else
            failed=$((failed + 1))
            echo "    ✗ Failed (no output file)"
        fi
    done
done

# Calculate success rate
success_rate=$((completed_successfully * 100 / total_analyses))

# Generate summary report
summary_file="./coexpression/binSpect_batch_summary.txt"
cat > "$summary_file" << EOF
Batch binSpect Analysis Summary
Generated on: $(date)

Total analyses: $total_analyses
Completed successfully: $completed_successfully
Failed: $failed
Success rate: ${success_rate}%
EOF

echo ""
echo "=== BATCH ANALYSIS SUMMARY ==="
cat "$summary_file"
echo ""
echo "Detailed results saved to: $summary_file"