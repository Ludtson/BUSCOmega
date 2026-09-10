#!/usr/bin/env python3
import os
import argparse
import logging

# Set up logging to file
logging.basicConfig(level=logging.DEBUG,
                    format='%(levelname)s: %(message)s',
                    filename='busco_seq_extractor.log',
                    filemode='w')

# Define a custom logging level to console (only warnings and errors)
console = logging.StreamHandler()
console.setLevel(logging.WARNING)
formatter = logging.Formatter('%(levelname)s: %(message)s')
console.setFormatter(formatter)
logging.getLogger('').addHandler(console)

# Define a function to log messages to both file and console
def log_to_file_and_console(logger, message, level=logging.INFO):
    logger.log(level, message)
    if level == logging.INFO:
        print(message)

def fetch_sequences(fasta_file):
    """
    Fetches sequences from a FASTA file.

    Parameters:
    fasta_file (str): Path to the FASTA file.

    Returns:
    dict: A dictionary of sequences with IDs as keys.
    """
    sequences = {}
    current_id = ""
    sequence = ""
    try:
        with open(fasta_file, 'r') as f:
            for line in f:
                if line.startswith('>'):
                    if current_id:
                        sequences[current_id] = sequence
                    current_id = line[1:]
                    sequence = ""
                else:
                    sequence += line.strip()
            if current_id:
                sequences[current_id] = sequence
    except FileNotFoundError:
        logging.error(f"FASTA file not found: {fasta_file}")
    return sequences

def get_seq(seq_dict, hdr):
    """
    Fetches the sequence from the dictionary where the key matches or contains the given header.

    Parameters:
    seq_dict (dict): Dictionary containing sequences.
    hdr (str): The header to match against the keys in the dictionary.

    Returns:
    str or None: The sequence if a match is found, otherwise None.
    """
    for key, seq in seq_dict.items():
        if key in hdr or hdr in key:
            return seq
    return None

def process_line(line, header):
    """
    Processes a line from the CSV file and maps it to header columns.

    Parameters:
    line (str): A line from the CSV file.
    header (list): The header columns of the CSV file.

    Returns:
    tuple: (BUSCO ID, dictionary of species data)
    """
    columns = line.strip().split(',')
    if len(columns) != len(header):
        logging.error(f"Mismatch in number of columns. Line: {line.strip()}")
        return None, None

    busco_id = columns[0]
    species_data = {header[i]: columns[i] for i in range(1, len(header))}
    return busco_id, species_data

def write_sequences(directory, busco_id, sequences, file_ext):
    """
    Writes the sequences to a file in the specified directory.

    Parameters:
    directory (str): Output directory.
    busco_id (str): BUSCO ID.
    sequences (dict): Dictionary of sequences to write.
    file_ext (str): File extension of the output file.
    """
    file_path = os.path.join(directory, f'{busco_id}.{file_ext}')
    with open(file_path, 'w') as f:
        for header, sequence in sequences.items():
            f.write(f'>{header}\n{sequence}\n')

def main(csv_file_path, cds_fasta_dir, prt_fasta_dir):
    """
    Main function to process the CSV file and extract sequences.

    Parameters:
    csv_file_path (str): Path to the CSV file summarizing BUSCO gene analysis.
    cds_fasta_dir (str): Directory containing species CDS FASTA files.
    prt_fasta_dir (str): Directory containing species protein FASTA files.
    """
    output_dir = 'busco_seqs'
    cds_dir = os.path.join(output_dir, 'busco_cds')
    prt_dir = os.path.join(output_dir, 'busco_prt')
    os.makedirs(cds_dir, exist_ok=True)
    os.makedirs(prt_dir, exist_ok=True)

    with open(csv_file_path, 'r') as csvfile:
        header = csvfile.readline().strip().split(',')
        lines = csvfile.readlines()

    cds_files = {species: os.path.join(cds_fasta_dir, f'{species}.fna') for species in header[1:]}
    prt_files = {species: os.path.join(prt_fasta_dir, f'{species}.faa') for species in header[1:]}

    cds_sequences = {species: fetch_sequences(fasta_file) for species, fasta_file in cds_files.items()}
    prt_sequences = {species: fetch_sequences(fasta_file) for species, fasta_file in prt_files.items()}

    total_busco_ids = 0
    busco_seqs_created = 0

    log_to_file_and_console(logging, "Fetching sequences...")
    for line in lines:
        busco_id, species_data = process_line(line, header)
        if not busco_id or not species_data:
            continue

        if all(value != 'NA' for value in species_data.values()):
            total_busco_ids += 1
            log_to_file_and_console(logging, f"Processing BUSCO ID: {busco_id}")
            valid_species = True

            for species, ids in species_data.items():
                if ids:
                    species_ids = set(ids.split(';'))
                    cds_seqs = cds_sequences[species]
                    prt_seqs = prt_sequences[species]

                    found_cds = found_prt = False
                    for species_id in species_ids:
                        exact_match_cds = get_seq(cds_seqs, species_id)
                        exact_match_prt = get_seq(prt_seqs, species_id)

                        if exact_match_cds and exact_match_prt:
                            found_cds = found_prt = True
                            break
                        elif exact_match_cds:
                            found_cds = True
                        elif exact_match_prt:
                            found_prt = True

                    if not found_cds or not found_prt:
                        logging.warning(f"    {species}:- No matching sequences found for both CDS and PRT. Skipping this BUSCO ID for all species.")
                        valid_species = False
                        break

            if valid_species:
                cds_to_write = {}
                prt_to_write = {}

                for species, ids in species_data.items():
                    logging.info(f"    {species}:")
                    if ids:
                        species_ids = set(ids.split(';'))
                        cds_seqs = cds_sequences[species]
                        prt_seqs = prt_sequences[species]

                        for species_id in species_ids:
                            exact_match_cds = get_seq(cds_seqs, species_id)
                            exact_match_prt = get_seq(prt_seqs, species_id)

                            if exact_match_cds and exact_match_prt:
                                logging.info(f"        {species_id}:- CDS_len={len(exact_match_cds)};PRT_len={len(exact_match_prt)}")
                                cds_to_write[f'{species} {species_id}'] = exact_match_cds
                                prt_to_write[f'{species} {species_id}'] = exact_match_prt
                                busco_seqs_created += 2
                            elif exact_match_cds:
                                logging.warning(f"========================{species_id}:- CDS found but PRT missing, skipping...")
                            elif exact_match_prt:
                                logging.warning(f"========================{species_id}:- PRT found but CDS missing, skipping...")

                if cds_to_write:
                    write_sequences(cds_dir, busco_id, cds_to_write, 'fna')
                if prt_to_write:
                    write_sequences(prt_dir, busco_id, prt_to_write, 'faa')

    log_to_file_and_console(logging, f"Total number of BUSCO IDs processed: {total_busco_ids}")
    log_to_file_and_console(logging, f"Total number of BUSCO sequences created: {busco_seqs_created}")

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Extracts corresponding BUSCO genes sequences for species with same lineage.db.")
    parser.add_argument('csv_file_path', type=str, help='Path to CSV summarizing BUSCO gene analysis. See merge_busco_table.py')
    parser.add_argument('cds_fasta_dir', type=str, help='Directory containing species CDS FASTA files')
    parser.add_argument('prt_fasta_dir', type=str, help='Directory containing species protein FASTA files')
    args = parser.parse_args()
    main(args.csv_file_path, args.cds_fasta_dir, args.prt_fasta_dir)
