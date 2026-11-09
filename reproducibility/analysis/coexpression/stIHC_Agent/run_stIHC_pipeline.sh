#!/bin/bash

# stIHC Pipeline Batch Execution Script
# Runs stIHC analysis on all available datasets (S1-S6 with suffixes 1-5)

echo "=== stIHC Pipeline Execution Started ==="
echo "Date: $(date)"
echo "Output directory: ./coexpression/stIHC_Agent"
echo ""

# Create log directory
mkdir -p ./coexpression/stIHC_Agent/logs

# Summary variables
total_jobs=0
completed_jobs=0
failed_jobs=0
start_time=$(date +%s)

# Progress tracking file
progress_file="./coexpression/stIHC_Agent/progress.txt"
echo "=== stIHC Pipeline Progress ===" > $progress_file
echo "Started: $(date)" >> $progress_file
echo "" >> $progress_file

# Run stIHC for all combinations
for S in {1..6}; do
    for suffix in {1..5}; do
        total_jobs=$((total_jobs + 1))
        
        echo "[$total_jobs/30] Processing S${S}_${suffix}..."
        
        # Log files
        log_file="./coexpression/stIHC_Agent/logs/S${S}_${suffix}.log"
        error_file="./coexpression/stIHC_Agent/logs/S${S}_${suffix}.err"
        
        # Run stIHC analysis
        if Rscript ./coexpression/stIHC.R $S $suffix > $log_file 2> $error_file; then
            completed_jobs=$((completed_jobs + 1))
            echo "  ✓ S${S}_${suffix} completed successfully"
            echo "✓ S${S}_${suffix}: SUCCESS" >> $progress_file
        else
            failed_jobs=$((failed_jobs + 1))
            echo "  ✗ S${S}_${suffix} failed (check logs)"
            echo "✗ S${S}_${suffix}: FAILED" >> $progress_file
        fi
    done
done

# Calculate execution time
end_time=$(date +%s)
execution_time=$((end_time - start_time))

echo ""
echo "=== Pipeline Execution Summary ==="
echo "Total jobs: $total_jobs"
echo "Completed: $completed_jobs"
echo "Failed: $failed_jobs"
echo "Success rate: $(echo "scale=2; $completed_jobs * 100 / $total_jobs" | bc)%"
echo "Total execution time: $execution_time seconds"

# Update progress file
echo "" >> $progress_file
echo "=== Final Summary ===" >> $progress_file
echo "Completed: $(date)" >> $progress_file
echo "Total jobs: $total_jobs" >> $progress_file
echo "Completed: $completed_jobs" >> $progress_file
echo "Failed: $failed_jobs" >> $progress_file
echo "Success rate: $(echo "scale=2; $completed_jobs * 100 / $total_jobs" | bc)%" >> $progress_file
echo "Total execution time: $execution_time seconds" >> $progress_file

# List generated files
echo ""
echo "=== Generated Files ==="
ls -la ./coexpression/stIHC_Agent/*.RDS 2>/dev/null | wc -l | xargs echo "RDS files generated:"
ls -la ./coexpression/stIHC_Agent/*.RDS 2>/dev/null || echo "No RDS files found"

echo ""
echo "=== stIHC Pipeline Execution Completed ==="