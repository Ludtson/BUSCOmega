#!/usr/bin/env bash

set -e  # Exit immediately if a command exits with a non-zero status.
set -u  # Treat unset variables as an error and exit immediately.

# Function to check if necessary software is available
check_software() {
  local missing=0

  for cmd in unzip gunzip; do
    if ! command -v "$cmd" &> /dev/null; then
      echo "Error: $cmd is not installed. Please install $cmd to proceed." | tee -a "$log_file"
      if [[ "$cmd" == "unzip" ]]; then
        echo "You can install it by running: sudo apt install unzip" | tee -a "$log_file"
      elif [[ "$cmd" == "gunzip" ]]; then
        echo "You can install it by running: sudo apt install gzip" | tee -a "$log_file"
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
  echo "Usage: $(basename "$0") <zipfile> [delete]"
  echo "Sort and process Phytozome files from a zip archive."
  echo
  echo "Arguments:"
  echo "  <zipfile>  The zip file containing Phytozome data."
  echo "  [delete]   (Optional) Delete the zip file after processing. Defaults to 'keep'."
  echo
  echo "Example:"
  echo "  $(basename "$0") data.zip delete"
  exit 1
}

# Function to clean genome name
clean_genome_name() {
  echo "$1" | sed -E 's/ [A-Za-z0-9]+\.[0-9]+$//'
}

# Main script
main() {
  if [[ $# -lt 1 ]]; then
    print_usage
  fi

  local zipfile=$1
  local delete_zip=${2:-keep}
  local datetime
  datetime=$(date +"%Y%m%d_%H%M%S")
  local zipname
  zipname=$(basename "$zipfile" .zip)
  local working_dir="${zipname}_${datetime}"
  local summary_csv="${working_dir}.csv"
  local log_file="${zipname}_process_${datetime}.log"

  # Check for necessary software
  check_software

  # Create working directory and unzip file
  mkdir -p "$working_dir"
  sudo unzip -X "$zipfile" -d "$working_dir" | tee -a "$log_file"

  # Change to working directory
  cd "$working_dir" || { echo "Failed to change directory to $working_dir" | tee -a "$log_file"; exit 1; }

  # Find the CSV file
  local csvfile
  csvfile=$(find . -name "*.csv")
  if [[ -z "$csvfile" ]]; then
    echo "CSV file not found" >&2 | tee -a "$log_file"
    exit 1
  fi

  # Create directories
  mkdir -p gnm cds gff prt

  # Initialize summary CSV
  echo "shortname,ID,genome,species,gnm,prt,cds,gff" > "$summary_csv"

  # Initialize an associative array to keep track of processed files
  declare -A summary

  # Read CSV file, skip header
  tail -n +2 "$csvfile" | while IFS=, read -r filename file_id jgi_grouping_id directory_path shortname genome_name md5_checksum file_size; do
    # Remove quotes and non-printable characters
    filename=$(echo "$filename" | tr -d '"' | tr -cd '[:print:]')
    directory_path=$(echo "$directory_path" | tr -d '"' | tr -cd '[:print:]')
    shortname=$(echo "$shortname" | tr -d '"' | tr -cd '[:print:]')
    genome_name=$(echo "$genome_name" | tr -d '"' | tr -cd '[:print:]')

    # Clean the genome name
    local clean_name
    clean_name=$(clean_genome_name "$genome_name")

    # Construct the full file path
    local full_path="./$directory_path/$filename"

    # Initialize the summary entry if it doesn't exist
    if [[ -z "${summary[$shortname]:-}" ]]; then
      summary[$shortname]="$jgi_grouping_id,$genome_name,$clean_name,no,no,no,no"
    fi

    # Determine file type and destination directory
    local dest_dir
    local new_ext
    local file_type
    case $filename in
      *.fa.gz)
        if [[ $filename == *cds.fa.gz ]]; then
          dest_dir="cds"
          new_ext="fna"
          file_type="cds"
        elif [[ $filename == *protein.fa.gz ]]; then
          dest_dir="prt"
          new_ext="faa"
          file_type="prt"
        else
          dest_dir="gnm"
          new_ext="fa"
          file_type="gnm"
        fi
        ;;
      *.gff*)
        dest_dir="gff"
        new_ext="gff"
        file_type="gff"
        ;;
      *)
        echo "Unknown file type: $filename" >&2 | tee -a "$log_file"
        continue
        ;;
    esac

    if [[ -f "$full_path" ]]; then
      local new_name="${shortname}.${new_ext}"
      gunzip -c "$full_path" > "$dest_dir/$new_name" | tee -a "$log_file"

      # Update summary for the current shortname
      IFS=, read -r ID genome species gnm prt cds gff <<< "${summary[$shortname]}"
      case $file_type in
        gnm) gnm="yes" ;;
        prt) prt="yes" ;;
        cds) cds="yes" ;;
        gff) gff="yes" ;;
      esac
      summary[$shortname]="$ID,$genome,$species,$gnm,$prt,$cds,$gff"

      # Check if all file types are found for the species and write to CSV
      if [[ $gnm == "yes" && $prt == "yes" && $cds == "yes" && $gff == "yes" ]]; then
        echo "$shortname,${summary[$shortname]}" >> "$summary_csv"
        unset summary[$shortname]  # Clear the entry after writing
      fi
    else
      echo "File not found: $full_path" >&2 | tee -a "$log_file"
    fi

  done

  # Check if CSV was written correctly
  if [[ ! -s "$summary_csv" ]]; then
    echo "Summary CSV is empty or not written correctly" >&2 | tee -a "$log_file"
    exit 1
  fi

  # Delete the Phytozome directory
  rm -rf "Phytozome" | tee -a "$log_file"

  # Delete the zip file if specified
  if [[ "$delete_zip" == "delete" ]]; then
    rm -f "../$zipfile" | tee -a "$log_file"
  fi

  echo "Unzipping completed." | tee -a "$log_file"
  echo "Directories created: gnm, cds, gff, prt" | tee -a "$log_file"
  echo "Files moved and renamed." | tee -a "$log_file"
  echo "Summary CSV generated: $summary_csv" | tee -a "$log_file"
  echo "ZIP file deleted: $delete_zip" | tee -a "$log_file"
}

main "$@"
