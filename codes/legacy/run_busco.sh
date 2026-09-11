#!/usr/bin/env bash

# Function to check if a utility is installed
check_utility() {
    local utility=$1
    if command -v "$utility" &> /dev/null; then
        echo "$utility is installed and available."
    else
        echo "$utility is not installed. Please install it before proceeding."
        exit 1
    fi
}

# Function to check if a directory exists and optionally create it
check_directory() {
    local dir=$1
    local create_if_missing=${2:-F}  # Default value is F (False)

    if [ -d "$dir" ]; then
        echo "Directory '$dir' exists."
    else
        if [ "$create_if_missing" == "T" ]; then
            mkdir -p "$dir"
            echo "Directory '$dir' did not exist, so it was created."
        else
            echo "Directory '$dir' does not exist."
            exit 1
        fi
    fi
}

# Function to log messages
log_message() {
    local message=$1
    echo "$(date '+%Y-%m-%d %H:%M:%S'): $message" | tee -a "$log_file"
}

# Function to parse and validate parameters
parse_params() {
    if [[ "$#" -lt 4 ]]; then
        echo "Usage: $0 <gnm_dir> <out_dir> <lineage> <type> [threads]"
        exit 1
    fi

    gnm_dir=$1
    out_dir=$2
    lineage=$3
    tipo=$4
    threads=${5:-1}  # Default to 1 thread if not specified

    if ! [[ "$threads" =~ ^[0-9]+$ ]]; then
        echo "Argument for threads must be a positive integer."
        exit 2
    fi
}

# Function to run BUSCO on a given FASTA file
run_busco() {
    local fasta=$1
    local nm=$(basename "$fasta" | cut -d '.' -f 1)
    local tar_option=${2:-false}

    log_message "Running BUSCO on $fasta"
    if ! busco \
        --offline \
        --in "$fasta" \
        --cpu "$threads" \
        --lineage_dataset "$lineage" \
        --mode "$tipo" \
        --out "$nm" \
        --out_path "$out_dir" \
        --force > "$out_dir/$nm"_std.out 2>&1; then
        log_message "BUSCO analysis failed for $fasta"
        exit 1
    fi

    if [ "$tar_option" = true ]; then
        tar -czf "$out_dir/$nm".tar.gz "$out_dir/$nm"
        rm -rf "$out_dir/$nm"
    fi

    log_message "Completed BUSCO for $fasta"
}

# Main script
main() {
    parse_params "$@"
    check_directory "$gnm_dir"
    cd "$gnm_dir" || { log_message "Failed to change directory to $gnm_dir"; exit 1; }

    check_utility "busco"
    check_directory "$out_dir" "T"
    check_directory "$lineage"

    log_dir="./logs"
    check_directory "$log_dir" "T"
    log_file="$log_dir/busco.log"

    local fasta_files=()
    for file in ./*.{faa,fa,fna}; do
        if [ -f "$file" ]; then
            fasta_files+=("$file")
        fi
    done

    for fasta in "${fasta_files[@]}"; do
        run_busco "$fasta"
    done

    log_message "All BUSCO analyses completed."

    if ! busco_summary.py "$out_dir" --outfile "$out_dir/busco_summary.csv"; then
        log_message "BUSCO summary generation failed."
        exit 1
    fi
    log_message "BUSCO summary table generated at $out_dir/busco_summary.csv"

    if [ -d "$./busco_downloads" ]; then
        rm -rf "./busco_downloads"
        log_message "Removed unnecessary busco_downloads directory."
    fi

    mv "$out_dir"/*.out "$log_dir/"
    log_message "Moved log files to $log_dir."
}

main "$@"
