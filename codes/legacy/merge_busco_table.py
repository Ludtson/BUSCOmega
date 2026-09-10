#!/usr/bin/env python3

from pathlib import Path
import csv
from datetime import datetime
import argparse

def parse_full_table(file_path):
    busco_data = []
    with open(file_path, 'r') as file:
        for line in file:
            if not line.startswith('#'):
                columns = line.strip().split('\t')
                if len(columns) >= 3:  # Ensure there are at least 3 columns
                    busco_data.append(columns)
                else:
                    # If the row has insufficient columns, treat it as Missing/Fragmented
                    busco_data.append([columns[0], 'NA', 'NA'])
    return busco_data

def merge_busco_tables(busco_dir, output_file=None):
    busco_path = Path(busco_dir).resolve()
    species_paths = [p for p in busco_path.iterdir() if p.is_dir()]

    all_busco_data = {}

    for species_path in species_paths:
        species_name = species_path.name
        full_table_path = species_path.glob('**/full_table.tsv')

        for file_path in full_table_path:
            busco_data = parse_full_table(file_path)

            for row in busco_data:
                busco_id = row[0]
                status = row[1]
                sequence = row[2]

                if busco_id not in all_busco_data:
                    all_busco_data[busco_id] = {}

                if species_name not in all_busco_data[busco_id]:
                    all_busco_data[busco_id][species_name] = []

                if status in ['Complete', 'Duplicated']:
                    if all_busco_data[busco_id][species_name] == ['NA']:
                        all_busco_data[busco_id][species_name] = [sequence]
                    else:
                        all_busco_data[busco_id][species_name].append(sequence)
                elif status in ['Missing', 'Fragmented'] or not all_busco_data[busco_id][species_name]:
                    all_busco_data[busco_id][species_name] = ['NA']

    # Prepare the final merged data
    merged_data = {'busco_id': []}

    for species_name in species_paths:
        merged_data[species_name.name] = []

    for busco_id, species_data in all_busco_data.items():
        merged_data['busco_id'].append(busco_id)

        for species_name in species_paths:
            sequences = species_data.get(species_name.name, ['NA'])
            merged_data[species_name.name].append("; ".join(sequences))

    # Determine the output file name
    if output_file is None:
        dt_string = datetime.now().strftime("%Y%m%d_%H%M%S")
        output_file = busco_path / f"{busco_path.name}_{dt_string}.csv"
    else:
        output_file = Path(output_file)

    # Write the merged data to a CSV file
    headers = ['busco_id'] + [species_path.name for species_path in species_paths]

    with open(output_file, 'w', newline='') as csvfile:
        csvwriter = csv.writer(csvfile, delimiter=',')
        csvwriter.writerow(headers)

        for i in range(len(merged_data['busco_id'])):
            row = [merged_data['busco_id'][i]] + [merged_data[species_path.name][i] for species_path in species_paths]
            csvwriter.writerow(row)

    print(f"Merged BUSCO table written to {output_file}")

def main():
    parser = argparse.ArgumentParser(description='Merge BUSCO full_table.tsv files.')
    parser.add_argument('busco_dir', help='Path to the BUSCO directory')
    parser.add_argument('output_file', nargs='?', default=None, help='Output CSV file name (optional)')

    args = parser.parse_args()
    merge_busco_tables(args.busco_dir, args.output_file)

if __name__ == "__main__":
    main()
