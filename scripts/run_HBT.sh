#!/bin/bash

num_loops=5
successful_runs=0

mkdir -p run_HBT_outputs

echo ""
echo "Temporary directory run_HBT_outputs created"
echo "Copying required files to current directory"
echo ""

cp config/config.ini .
cp input_data/particle_list_set_0.bin .

if [ ! -f "config.ini" ]; then
    echo "ERROR: config.ini not found!" >&2
    exit 1
fi

if [ ! -f "particle_list_set_0.bin" ]; then
    echo "ERROR: particle_list_set_0.bin not found!" >&2
    exit 1
fi

echo "All necessary files copied successfully"

g++ -o run CQ_Calculation.cpp main.cpp -fopenmp -O3

echo "$num_loops runs started..."
echo ""

for ((i = 1; i <= num_loops; i++)); do

    nohup ./run > "run_HBT_outputs/run_${i}.log" 2>&1 &
    wait 

    if [ -f "Cqinv_kT_0.data" ]; then
        mv Cqinv_kT_0.data Cqinv_kT_0_${i}.data
        mv Cqinv_kT_0_${i}.data run_HBT_outputs/
        successful_runs=$((successful_runs + 1))
        echo "Run $i completed successfully"
    else
        echo "ERROR: Cqinv_kT_0.data not found after run $i" >&2
    fi
done

echo "Completed $successful_runs successful runs out of $num_loops"

if [ "$successful_runs" -gt 0 ]; then
    data_files=(run_HBT_outputs/Cqinv_kT_0_*.data)

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
    echo "No successful runs :("
fi

rm -rf run_HBT_outputs
echo "The newly created temporary directory run_HBT_outputs has been removed"

rm -rf config.ini particle_list_set_0.bin
