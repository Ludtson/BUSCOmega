#!/usr/bin/env bash

# Function to check if necessary software is available
check_software() {
  local missing=0

  for cmd in pal2nal.pl; do
    if ! command -v "$cmd" &> /dev/null; then
      echo "Error: $cmd is not installed. Please install $cmd to proceed." | tee -a "$log_file"
      if [[ "$cmd" == "pal2nal.pl" ]]; then
        echo "You can install it by following the instructions at: http://www.bork.embl.de/pal2nal/" | tee -a "$log_file"
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
  echo "Usage: $(basename "$0") <prt_aln_dir> <cds_dir> [-e <extension>]"
  echo "Run PAL2NAL alignments on all CDS files (*.fna) in a dir using corresponding protein alignment files (*.aln by default) in another dir."
  echo
  echo "Arguments:"
  echo "  <prt_aln_dir>  The directory containing the protein alignment files."
  echo "  <cds_dir>      The directory containing the corresponding CDS files."
  echo "  -e <extension> (Optional) The input file extension. Defaults to .aln if not provided."
  echo
  echo "Example:"
  echo "  $(basename "$0") /path/to/protein_alignments /path/to/cds_files -e .aln"
  exit 1
}

# Function to run PAL2NAL alignment
run_pal2nal_alignment() {
  local prt_aln=$1
  local cds_file=$2
  local out=$3
  local log_file=$4

  echo "Running CDS alignment for BUSCO gene: $bg_name" | tee -a "$log_file"
  pal2nal.pl "$prt_aln" "$cds_file" -output paml -nogap -nomismatch > "$out" 2>>"$log_file"
}

# Main script
main() {
  local ext=".aln"  # Default extension

  while getopts ":e:" opt; do
    case ${opt} in
      e )
        ext=$OPTARG
        ;;
      \? )
        echo "Invalid option: -$OPTARG" 1>&2
        print_usage
        ;;
      : )
        echo "Invalid option: -$OPTARG requires an argument" 1>&2
        print_usage
        ;;
    esac
  done
  shift $((OPTIND -1))

  if [[ $# -lt 2 ]]; then
    echo "#ERROR!!! 2 positional arguments required"
    print_usage
  fi

  local prt_aln_dir=$1
  local cds_dir=$2
  local datetime
  datetime=$(date +"%Y%m%d_%H%M%S")
  local log_file="$prt_aln_dir/run_pal2nal_$datetime.log"

  if [[ ! -d "$prt_aln_dir" ]]; then
    echo "$prt_aln_dir is not a directory or does not exist!" | tee -a "$log_file"
    exit 1
  fi

  if [[ ! -d "$cds_dir" ]]; then
    echo "$cds_dir is not a directory or does not exist!" | tee -a "$log_file"
    exit 1
  fi

  cd "$prt_aln_dir" || { echo "Failed to change directory to $prt_aln_dir" | tee -a "$log_file"; exit 1; }

  # Check if there are no .ext files and give a friendly message
  if [[ "$(find . -maxdepth 1 -type f -name "*$ext")" == "" ]]; then
    echo "No protein alignment files (*$ext) found in $prt_aln_dir. Please ensure your files have the correct $ext extension." | tee -a "$log_file"
    exit 1
  fi

  local p_a_d=../p_aln_dir
  mkdir -p "$p_a_d" || { echo "Failed to create directory $p_a_d"; tee -a "$log_file"; exit 1; }

  # Check for necessary software
  check_software

  find . -maxdepth 1 -type f -name "*$ext" -print0 | while IFS= read -r -d '' prt_aln; do
    local prt_aln_fname
    prt_aln_fname=$(basename "$prt_aln")
    local bg_name="${prt_aln_fname%.*}"
    local out="$p_a_d"/"$bg_name"_"$datetime".pml
    for ext in fa fna faa; do
      if [ -f "$cds_dir/$bg_name.$ext" ]; then
        local cds_file="$cds_dir/$bg_name.$ext"
        run_pal2nal_alignment "$prt_aln" "$cds_file" "$out" "$log_file"
      fi
    done
  done
}

main "$@"
