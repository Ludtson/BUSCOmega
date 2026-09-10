#!/usr/bin/env bash

log_file="script_log_$(date +'%Y%m%d_%H%M%S').log"
max_concurrent_processes=${6:-4}  # Set the maximum number of concurrent processes with a default value of 4

log() {
    echo "$(date +'%Y-%m-%d %H:%M:%S') - $1" | tee -a "$log_file"
}

check_command() {
    if ! command -v "$1" &> /dev/null; then
        log "Error: $1 could not be found in your PATH. Please make sure it's installed and in the PATH."
        exit 1
    fi
}

validate_arguments() {
    if [[ -z "$seq_dir" || -z "$n" || -z "$TREEFILE" || -z "$CONTROL_TYPE" ]]; then
        log "Error: Missing required arguments."
        display_usage
        exit 1
    fi

    if [[ ! -d "$seq_dir" ]]; then
        log "Error: Directory $seq_dir does not exist."
        exit 1
    fi

    if [[ "$CONTROL_TYPE" == "reduced_all_sites" && -z "$NS_SITES" ]]; then
        log "Error: --ns_sites must be specified for reduced_all_sites control type."
        exit 1
    fi
}

display_usage() {
    script_name=$(basename "$0")
    echo "Usage: $script_name seq_dir n TREEFILE CONTROL_TYPE [NS_SITES] [max_concurrent_processes]"
    echo
    echo "seq_dir                 Directory containing .pml files"
    echo "n                       Number of groups to divide the files into"
    echo "TREEFILE                Newick format tree file of species"
    echo "CONTROL_TYPE            Type of control file to generate (zero_model, all_sites, reduced_all_sites, free_model)"
    echo "NS_SITES                (Optional) NSsites values for reduced model (e.g., '0,1,2' or '0,7,8'). Required if type is 'reduced_all_sites'"
    echo "max_concurrent_processes (Optional) Maximum number of concurrent processes (default: 4)"
    echo
}

create_control_file() {
    local output_dir=$1
    local output_file=$2
    local TREEFILE=$3
    local CONTROL_TYPE=$4
    local NS_SITES=$5
    local num_genes=$6

    local cmd="generate_codeml_control.py -t $CONTROL_TYPE -g $num_genes -o $output_dir/$output_dir.ctrl $output_file $TREEFILE"

    if [[ "$CONTROL_TYPE" == "reduced_all_sites" && -n "$NS_SITES" ]]; then
        cmd+=" --ns_sites $NS_SITES"
    fi

    log "Running: $cmd"
    eval "$cmd"

    if [[ $? -ne 0 ]]; then
        log "Error: generate_codeml_control.py script failed."
        exit 1
    fi
}

run_codeml() {
    local output_dir=$1
    local ctrl_file=$2
    local base_ctrl_file=$(basename "$ctrl_file")

    log "Running CODEML for $ctrl_file"
    codeml "$ctrl_file" > "$output_dir/${base_ctrl_file%.ctrl}_$(date +'%Y%m%d_%H%M%S').log" &
    manage_concurrency
}

manage_concurrency() {
    while [[ $(jobs -r -p | wc -l) -ge $max_concurrent_processes ]]; do
        wait -n
    done
}

main() {
    if [[ "$1" == "-h" || "$1" == "--help" ]]; then
        display_usage
        exit 0
    fi

    seq_dir=$1
    n=$2
    TREEFILE=$3
    CONTROL_TYPE=$4
    NS_SITES=$5

    validate_arguments
    check_command "generate_codeml_control.py"
    check_command "codeml"

    cd "$seq_dir"
    files=(*.pml)
    num_files=${#files[@]}
    files_per_group=$((num_files / n))
    remainder=$((num_files % n))
    log "Total number of .pml files: ${#files[@]}"

    start_index=0
    stop_index=$((files_per_group))

    for i in $(seq 1 "$n"); do
        if [[ ! $i -eq 1 ]]; then
            start_index="$stop_index"
        fi
        stop_index=$((start_index + files_per_group))
        if [[ $i -eq $n ]]; then
            files_per_group=$((files_per_group + remainder))
        fi
        output_file="files_$((start_index + 1))_to_$((stop_index)).pml"
        output_dir="files_$((start_index + 1))_to_$((stop_index))"
        ctd_files=("${files[@]:start_index:files_per_group}")
        mkdir -p "$output_dir"
        for file in "${ctd_files[@]}"; do
            cat "$file" >> "$output_dir/$output_file"
        done
        log "Concatenated ${#ctd_files[@]} files into $output_file"
        create_control_file "$output_dir" "$output_file" "$TREEFILE" "$CONTROL_TYPE" "$NS_SITES" "${#ctd_files[@]}"
        run_codeml "$output_dir" "$output_dir/$output_dir.ctrl"
    done

    wait
    log "All CODEML runs are complete."
}

main "$@"
