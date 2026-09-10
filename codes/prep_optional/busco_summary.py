#!/usr/bin/env python3

import argparse
import logging
import os
import re
import sys

# Configure logging
def configure_logging(log_level):
    logging.basicConfig(level=log_level, format='%(asctime)s - %(levelname)s - %(message)s')

def parse_busco_output(filename):
    """
    Parse a BUSCO output file for species name and BUSCO scores.

    Args:
        filename (str): The path to the BUSCO output file.

    Returns:
        tuple: species name and a list of BUSCO scores.
    """
    data = []
    species = ""
    busco_score_pattern = r"^\t[0-9]+\t"
    species_pattern = re.compile(r"^# Summari[sz]ed benchmarking in BUSCO notation for file .*/(.*)\.faa")

    try:
        with open(filename, 'r') as f:
            lines = f.readlines()
            for line in lines:
                species_match = species_pattern.match(line)
                if species_match:
                    species = species_match.group(1)
                if re.search(busco_score_pattern, line):
                    col = line.strip().split("\t")
                    if len(col) > 0:
                        data.append(col[0])
    except FileNotFoundError:
        logging.error(f"File not found: {filename}")
    except Exception as e:
        logging.error(f"Error reading file {filename}: {e}")
    return species, data

def busco_dir(dir_name):
    """
    Process the given directory for a BUSCO summary file.

    Args:
        dir_name (str): The directory to search for a BUSCO summary file.

    Returns:
        tuple: species name and a list of BUSCO scores.
    """
    for file in os.listdir(dir_name):
        if re.search(r"^short_summary", file) and file.endswith('.txt'):
            return parse_busco_output(os.path.join(dir_name, file))
    logging.warning(f"No short_summary text file found in directory: {dir_name}")
    return "", []

def busco_table(folder):
    """
    Generate a summary table of BUSCO results for all directories in the folder.

    Args:
        folder (str): The folder containing directories of BUSCO results.

    Returns:
        list: A list of dictionaries containing BUSCO results for each species.
    """
    busco_data = []
    for d in os.listdir(folder):
        dir_path = os.path.join(folder, d)
        if os.path.isdir(dir_path):  # Only process directories
            species_name, busco_scores = busco_dir(dir_path)
            if busco_scores and len(busco_scores) >= 6:
                entry = {
                    "Species": species_name,
                    "CompleteBUSCOs": busco_scores[0],
                    "SingleBUSCOs": busco_scores[1],
                    "DuplicatedBUSCOs": busco_scores[2],
                    "FragmentedBUSCOs": busco_scores[3],
                    "MissingBUSCOs": busco_scores[4],
                    "TotalBUSCOs": busco_scores[5]
                }
                busco_data.append(entry)
    return busco_data

def main():
    """
    Main function to execute the script.
    """
    parser = argparse.ArgumentParser(description="Summarize BUSCO results.")
    parser.add_argument("folder", help="Folder containing directories of BUSCO results.")
    parser.add_argument("--outfile", "-o", help="Output file name", default="busco_table_all.csv")
    parser.add_argument("--log", "-l", help="Set the logging level", default="INFO", choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"])
    args = parser.parse_args()

    # Set logging level
    configure_logging(getattr(logging, args.log.upper(), None))

    try:
        busco_summary = busco_table(args.folder)

        with open(args.outfile, 'w') as f:
            f.write("Species,CompleteBUSCOs,SingleBUSCOs,DuplicatedBUSCOs,FragmentedBUSCOs,MissingBUSCOs,TotalBUSCOs\n")
            for row in busco_summary:
                f.write(','.join([row["Species"], row["CompleteBUSCOs"], row["SingleBUSCOs"], row["DuplicatedBUSCOs"],
                                  row["FragmentedBUSCOs"], row["MissingBUSCOs"], row["TotalBUSCOs"]]) + '\n')
        logging.info(f"BUSCO summary saved to {args.outfile}")
    except Exception as e:
        logging.critical(f"An unexpected error occurred: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()
