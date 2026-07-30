#!/bin/bash
set -euo pipefail

read -r -p "Simulation code: " sim_code
read -r -p "Output filename: " output_file
read -r -p "How many configurations? [1-7]: " njobs
read -r -p "Maximum simultaneous GPU jobs [default: 1]: " max_parallel
read -r -p "GPU number [default: 0]: " gpu_id

max_parallel=${max_parallel:-1}
gpu_id=${gpu_id:-0}

# Validate integer inputs
if ! [[ "$njobs" =~ ^[0-9]+$ ]] ||
   ! [[ "$max_parallel" =~ ^[0-9]+$ ]] ||
   ! [[ "$gpu_id" =~ ^[0-9]+$ ]]; then
    echo "Error: job counts and GPU number must be nonnegative integers."
    exit 1
fi

if (( njobs < 1 || njobs > 7 )); then
    echo "Error: number of configurations must be between 1 and 7."
    exit 1
fi

if (( max_parallel < 1 )); then
    echo "Error: simultaneous GPU jobs must be at least 1."
    exit 1
fi

if (( max_parallel > njobs )); then
    max_parallel=$njobs
fi

if [[ ! -f "$sim_code" ]]; then
    echo "Error: simulation source file not found: $sim_code"
    exit 1
fi

if ! command -v nvfortran >/dev/null 2>&1; then
    echo "Error: nvfortran was not found in PATH."
    echo "Check with: nvfortran -V"
    exit 1
fi

if ! command -v nvidia-smi >/dev/null 2>&1; then
    echo "Warning: nvidia-smi was not found."
    echo "The compiler may work, but GPU availability could not be checked."
else
    gpu_count=$(nvidia-smi --query-gpu=name \
        --format=csv,noheader 2>/dev/null | wc -l)

    if (( gpu_id >= gpu_count )); then
        echo "Error: GPU ${gpu_id} does not exist."
        echo "Detected ${gpu_count} GPU(s), numbered 0 through $((gpu_count - 1))."
        exit 1
    fi

    gpu_name=$(nvidia-smi \
        --query-gpu=name \
        --format=csv,noheader \
        --id="$gpu_id")

    echo "Using GPU ${gpu_id}: ${gpu_name}"
fi

exe_name="neighbor_modes_gpu"

echo "Compiling ${sim_code}..."

nvfortran \
    -O3 \
    -mp=gpu \
    -gpu=cc89 \
    -Minfo=mp \
    -Mextend \
    "$sim_code" \
    -o "$exe_name"

echo "Compilation successful."
echo "Running ${njobs} configurations."
echo "Maximum simultaneous GPU jobs: ${max_parallel}"
echo "GPU number: ${gpu_id}"

exe="$(realpath "./${exe_name}")"

run_dirs=()

for ((i=0; i<njobs; i++)); do
    run_dir="run_${i}"
    mkdir -p "$run_dir"
    run_dirs+=("$run_dir")
done

# GPU target regions must execute on the GPU.
export OMP_TARGET_OFFLOAD=MANDATORY

# Your timestep loop is serial on the host. Each process only needs one
# regular CPU thread to launch and coordinate GPU kernels.
export OMP_NUM_THREADS=1
export OMP_DYNAMIC=FALSE

# Start high-resolution simulation timer after compilation.
start_ns=$(date +%s%N)

failed=0

# Run configurations in batches so no more than max_parallel jobs
# compete for the GPU at once.
for ((batch_start=0; batch_start<njobs; batch_start+=max_parallel)); do
    batch_end=$((batch_start + max_parallel))

    if (( batch_end > njobs )); then
        batch_end=$njobs
    fi

    pids=()
    configs=()

    for ((i=batch_start; i<batch_end; i++)); do
        run_dir="run_${i}"

        echo "Starting config ${i} on GPU ${gpu_id}"

        (
            cd "$run_dir"

            CUDA_VISIBLE_DEVICES="$gpu_id" \
                "$exe" "$i" \
                > "config${i}.log" 2>&1
        ) &

        pids+=("$!")
        configs+=("$i")
    done

    # Wait for the current batch before starting the next batch.
    for ((j=0; j<${#pids[@]}; j++)); do
        config="${configs[$j]}"

        if ! wait "${pids[$j]}"; then
            echo "Error: config ${config} failed."
            echo "See run_${config}/config${config}.log"
            failed=1
        else
            echo "Config ${config} completed."
        fi
    done
done

# Stop timer immediately after all simulations finish.
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
