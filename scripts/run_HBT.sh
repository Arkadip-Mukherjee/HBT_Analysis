#!/bin/bash

# Set number of loops
num_loops=5
successful_runs=0

# Set up directories
mkdir -p run_HBT_outputs
#mkdir -p input_data
echo ""
echo "Temporary directory run_HBT_outputs created"

# Copy required files to current directory
cp config/config.ini .  # Essential configuration
cp input_data/particle_list_set_0.bin .  # Input data

# Verify files exist
if [ ! -f "config.ini" ]; then
    echo "ERROR: config.ini not found!" >&2
    exit 1
fi

if [ ! -f "particle_list_set_0.bin" ]; then
    echo "ERROR: particle_list_set_0.bin not found!" >&2
    exit 1
fi

# Compile ALL required files
g++ -o run CQinv_Calculation.cpp main.cpp -fopenmp -O3

echo "$num_loops runs started..."
echo ""

# Run the program num_loops times
for ((i = 1; i <= num_loops; i++)); do

    # Run the program in the background with nohup
    nohup ./run > "run_HBT_outputs/run_${i}.log" 2>&1 &
    wait  # Wait for the background process to finish

    # Check and move the output file
    if [ -f "Cqinv_kT_0.data" ]; then
        mv Cqinv_kT_0.data Cqinv_kT_0_${i}.data
        mv Cqinv_kT_0_${i}.data run_HBT_outputs/
        successful_runs=$((successful_runs + 1))
        echo "Run $i completed successfully"
    else
        echo "ERROR: Cqinv_kT_0.data not found after run $i" >&2
        # Don't decrement counter, just report error
    fi
done

echo "Completed $successful_runs successful runs out of $num_loops"

# Process results and calculate averages
if [ "$successful_runs" -gt 0 ]; then
    # Create a list of all data files
    data_files=(run_HBT_outputs/Cqinv_kT_0_*.data)

    # Use awk to process all files at once
    awk -v total_runs="$successful_runs" '
    BEGIN {
        # Initialize arrays
        for (row=1; row<=16; row++) {
            sum1[row] = 0
            sum2[row] = 0
        }
    }
    {
        # Get the current row number
        row = FNR

        # Only process the first 16 rows
        if (row <= 16) {
            sum1[row] += $1
            sum2[row] += $2
        }
    }
    END {
        # Calculate and print averages
        for (row=1; row<=16; row++) {
            avg1 = sum1[row] / total_runs
            avg2 = sum2[row] / total_runs
            printf "%f %f\n", avg1, avg2
        }
    }' "${data_files[@]}" > Avg_Cqinv.data
    echo "Averages calculated. Results saved to Avg_Cqinv.data"
else
    echo "No successful runs. Creating empty Avg_Cqinv.data"
fi

# Clean up - remove the temporary output directory
rm -rf run_HBT_outputs
echo "The newly created temporary directory run_HBT_outputs has been removed"

# Cleanup copied files
rm config.ini particle_list_set_0.bin
