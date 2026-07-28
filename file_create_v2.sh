#!/bin/bash

read -p "Output filename: " output_file
read -p "How many jobs? " njobs

gfortran -O3 old_v1.f

if [ $? -ne 0 ]; then
    echo "Error: Compilation failed."
    exit 1
fi

echo "Compilation successful. Executing $njobs runs in parallel..."

pids=()

# Launch jobs in parallel
for ((i=0; i<njobs; i++)); do
    ./a.exe "$i" &
    pids+=($!)
done

# Wait for all jobs to finish
for pid in "${pids[@]}"; do
    wait "$pid"
done

echo "All runs completed."

files=()

# Build list of output files
for ((i=0; i<njobs; i++)); do
    files+=("config${i}_total.txt")
done

# Combine outputs
paste "${files[@]}" > "$output_file"

# Cleanup
rm "${files[@]}"
