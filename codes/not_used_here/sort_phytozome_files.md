README: sort_phytozome_files.sh
Description
The sort_phytozome_files.sh script is designed to organize and process files from a Phytozome download. It extracts the files from the zip archive, sorts them into appropriate directories based on their type (genome, protein, cds, gff), renames them for consistency, and generates a summary CSV file.

Usage
bash
./sort_phytozome_files.sh <zipfile> [delete]
Parameters
<zipfile>: Path to the Phytozome zip file that contains the files and a CSV file with metadata.

[delete]: Optional parameter. If set to delete, the script will remove the zip file after extraction. Default behavior is to keep the zip file.

Prerequisites
Ensure you have bash and gunzip installed on your system.

The zip file should contain the following:

A directory named phytozome with all necessary files.

A CSV file listing the metadata of the files.

Steps
Create a Working Directory: The script creates a working directory named phytozome_<datetime>, where <datetime> is the current date and time in the format YYYYMMDD_HHMMSS.

Unzip the Archive: The script unzips the contents of the provided zip file into the working directory.

Process the CSV File: The script finds the CSV file in the working directory and reads it line by line, skipping the header.

Organize Files: Based on the file type (genome, protein, cds, gff), the script moves the files into their respective directories (gnm, cds, prt, gff) and renames them for consistency.

Generate Summary CSV: The script creates a summary CSV file named summary.csv with the following columns:

shortname: Short name of the organism.

ID: File ID from the CSV.

genome: Genome name from the CSV.

gnm, prt, cds, gff: Indicate if the respective files were found and processed (yes or no).

Optional Cleanup: If the delete parameter is specified, the script removes the original zip file after processing.

Example
bash
./sort_phytozome_files.sh phytozome_download.zip delete
This example runs the script to process the phytozome_download.zip file and deletes the zip file after processing.

Output
Working Directory: phytozome_<datetime>

Organized Directories: gnm, cds, gff, prt

Summary CSV: summary.csv

Notes
Ensure the CSV file and the Phytozome directory are correctly structured as expected.

If any files are not found during processing, a message will be printed, and the summary CSV will be updated accordingly.

