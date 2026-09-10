#!/usr/bin/env bash

# Function to check if necessary software is available
check_software() {
  local missing=0

  for cmd in mafft perl; do
    if ! command -v "$cmd" &> /dev/null; then
      echo "Error: $cmd is not installed. Please install $cmd to proceed." | tee -a "$log_file"
      if [[ "$cmd" == "mafft" ]]; then
        echo "You can install it by running: sudo apt install mafft" | tee -a "$log_file"
      elif [[ "$cmd" == "perl" ]]; then
        echo "You can install it by running: sudo apt install perl" | tee -a "$log_file"
      fi
      missing=1
    fi
  done

  if [[ $missing -eq 1 ]]; then
    exit 1
  fi
}

# Function to print usage information
print_usage() {
  echo "Usage: $(basename "$0") <directory> [threads]"
  echo "Run MAFFT alignments on all protein FASTA files (*.faa) in the specified directory."
  echo
  echo "Arguments:"
  echo "  <directory>  The directory containing the input files."
  echo "  [threads]    (Optional) The number of threads to use. Defaults to single-threaded mode if not provided."
  echo
  echo "Example:"
  echo "  $(basename "$0") /path/to/files 4"
  exit 1
}

# Function to run MAFFT alignment
run_alignment() {
  local file=$1
  local threads=$2
  local out=$3
  local log_file=$4

  if [[ -n "$threads" ]]; then
    echo "Running <mafft> with $threads cores" | tee -a "$log_file"
    mafft --inputorder --auto --thread "$threads" "$file" > "$out" 2>>"$log_file"
  else
    echo "Running <mafft> in single-threaded mode" | tee -a "$log_file"
    mafft --inputorder --auto "$file" > "$out" 2>>"$log_file"
  fi
}

# Main script
main() {
  if [[ $# -lt 1 || $# -gt 2 ]]; then
    print_usage
  fi

  local bg_dir=$1
  local threads=$2
  local datetime
  datetime=$(date +"%Y%m%d_%H%M%S")
  local log_file="$bg_dir/run_mafft_$datetime.log"

  if [[ ! -d "$bg_dir" ]]; then
    echo "$bg_dir is not a directory or does not exist!" | tee -a "$log_file"
    exit 1
  fi

  if [[ "$#" -ge 2 ]]; then
    if [[ ! "$2" =~ ^[0-9]+$ ]]; then
      echo "Argument #2 is the number of threads; $2 is not a number" | tee -a "$log_file"
      exit 2
    fi
  fi

  cd "$bg_dir" || { echo "Failed to change directory to $bg_dir" | tee -a "$log_file"; exit 1; }

  # Check if there are no .faa files and give a friendly message
  if [[ "$(find . -maxdepth 1 -type f -name '*.faa')" == "" ]]; then
    echo "No protein FASTA files (*.faa) found in $bg_dir. Please ensure your files have the correct .faa extension." | tee -a "$log_file"
    exit 1
  fi

  local m_a_d=../m_aln_dir
  mkdir -p "$m_a_d" || { echo "Failed to create directory $m_a_d" | tee -a "$log_file"; exit 1; }

  # Check for necessary software
  check_software

  find . -maxdepth 1 -type f -name "*.faa" -print0 | while IFS= read -r -d '' file; do
    local filename
    filename=$(basename "$file")
    local bg_name="${filename%.*}"
    local out="$m_a_d/$bg_name.aln"
    run_alignment "$file" "$threads" "$out" "$log_file"
  done
}

main "$@"

  # find . -maxdepth 1 -type f \( -name "*.fna" -o -name "*.faa" \) -print0 | while IFS= read -r -d '' file; do

  