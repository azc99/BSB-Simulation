#!/bin/bash
set -euo pipefail

read -r -p "Simulation code: " sim_code
read -r -p "Output filename: " output_file
read -r -p "How many jobs? " njobs
read -r -p "OpenMP threads per job [default: 4]: " threads_per_job

threads_per_job=${threads_per_job:-4}
total_cpus=$(nproc)
required_cpus=$((njobs * threads_per_job))

if (( njobs < 1 || threads_per_job < 1 )); then
    echo "Error: jobs and threads must be positive integers."
    exit 1
fi

if (( required_cpus > total_cpus )); then
    echo "Error: requested ${required_cpus} CPU threads,"
    echo "but only ${total_cpus} logical CPUs are available."
    exit 1
fi

echo "Compiling ${sim_code}..."

gfortran \
    -O3 \
    -march=native \
    -fopenmp \
    -ffixed-line-length-none \
    "$sim_code" \
    -o neighbor_modes

echo "Compilation successful."
echo "Running ${njobs} jobs with ${threads_per_job} threads each."
echo "Total CPU threads used: ${required_cpus}/${total_cpus}"

export OMP_NUM_THREADS="$threads_per_job"
export OMP_DYNAMIC=FALSE
export OMP_PROC_BIND=close
export OMP_PLACES=cores


exe="$(realpath ./neighbor_modes)"
pids=()
run_dirs=()

# Start high-resolution simulation timer
start_ns=$(date +%s%N)

for ((i=0; i<njobs; i++)); do
    run_dir="run_${i}"
    mkdir -p "$run_dir"
    run_dirs+=("$run_dir")

    first_cpu=$((i * threads_per_job))
    last_cpu=$((first_cpu + threads_per_job - 1))

    echo "Starting config ${i} on CPUs ${first_cpu}-${last_cpu}"

    (
        cd "$run_dir"

        # write(6,*) and error output go into this log file
        taskset -c "${first_cpu}-${last_cpu}" \
            "$exe" "$i" > "config${i}.log" 2>&1
    ) &

    pids+=("$!")
done

failed=0

for ((i=0; i<njobs; i++)); do
    if ! wait "${pids[$i]}"; then
        echo "Error: config ${i} failed."
        echo "See run_${i}/config${i}.log"
        failed=1
    fi
done

# Stop timer immediately after simulations finish
end_ns=$(date +%s%N)

elapsed_ms=$(((end_ns - start_ns) / 1000000))

minutes=$((elapsed_ms / 60000))
seconds=$(((elapsed_ms % 60000) / 1000))
milliseconds=$((elapsed_ms % 1000))

if (( failed != 0 )); then
    printf 'Simulation runtime: %02d:%02d.%03d\n' \
        "$minutes" "$seconds" "$milliseconds"

    echo "One or more simulations failed."
    exit 1
fi

echo "All runs completed."

files=()

for ((i=0; i<njobs; i++)); do
    file="run_${i}/config${i}_total.txt"

    if [[ ! -f "$file" ]]; then
        echo "Error: expected output file not found: $file"
        exit 1
    fi

    files+=("$file")
done

paste "${files[@]}" > "$output_file"

echo "Combined output written to: $output_file"

printf 'Simulation runtime: %02d:%02d.%03d\n' \
    "$minutes" "$seconds" "$milliseconds"

read -r -p "Delete individual run directories? [y/N]: " cleanup

if [[ "$cleanup" =~ ^[Yy]$ ]]; then
    rm -rf "${run_dirs[@]}"
    echo "Individual run directories removed."
fi
