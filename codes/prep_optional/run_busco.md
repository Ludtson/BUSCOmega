Documentation for BUSCO Analysis Bash Script
Overview
The run_busco.sh script automates the process of running BUSCO (Benchmarking Universal Single-Copy Orthologs) on a set of FASTA files and generates a summary table of the results. It includes utility checks, directory validations, and parallel processing to enhance performance.

Dependencies
BUSCO: The main tool for benchmarking single-copy orthologs.
GNU Parallel: For parallel processing of multiple FASTA files.
Python: To run the summary table generation script.
curl: To download BUSCO lineage data.
tar: To extract downloaded lineage data.
Installation
Install BUSCO:
conda install -c bioconda busco

Install GNU Parallel:
sudo apt-get install parallel

Install Python 3:
sudo apt-get install python3

Download BUSCO Lineage Data:
curl https://busco-data.ezlab.org/v5/data/lineages/viridiplantae_odb10.2024-01-08.tar.gz --output ./busco/viridiplantae_odb10.2024-01-08.tar.gz
tar -xvf viridiplantae_odb10.2024-01-08.tar.gz -C /scratch/user/ludtson/busco_downloads/lineages
Note: The lineage data URL and file name will vary depending on the specific lineage you are working with. Adjust the URL and file name accordingly.

Script Usage
Save the script as run_busco.sh.
Ensure You Have the busco_summary.py
Run the Bash Script:
run_busco.sh /path/to/fasta/files /path/to/output /path/to/lineage genome 4

This command will:
Check for required utilities and directories.
Run BUSCO on each FASTA file in parallel.
Generate a summary table of the BUSCO results and save it as busco_summary.csv in the output directory.

Parameters for run_busco.sh Script
Directory with FASTA Files:
    Description: The path to the directory containing the FASTA files you want to analyze with BUSCO.
    Example: /path/to/fasta/files
Desired Output Directory:
    Description: The path to the directory where you want to save the BUSCO output results.
    Example: /path/to/output
Lineage Directory:
    Description: The path to the directory containing the BUSCO lineage data. This data is necessary for BUSCO to perform its analysis.
    Example: /path/to/lineage
    Note: You can download the lineage data using the following commands:
        url https://busco-data.ezlab.org/v5/data/lineages/viridiplantae_odb10.2024-01-08.tar.gz --output ./busco/viridiplantae_odb10.2024-01-08.tar.gz
        tar -xvf viridiplantae_odb10.2024-01-08.tar.gz -C /scratch/user/ludtson/busco_downloads/lineages
FASTA Data Type:
    Description: The type of data in the FASTA files. This can be genome, protein, or transcript.
    Example: genome
Number of Cores to Use (default is 1):
    Description: The number of CPU cores to use for running BUSCO. If not specified, the script will use 1 core by default.
    Example: 4

Example Command
bash run_busco.sh /path/to/fasta/files /path/to/output /path/to/lineage genome 4

This command will:
    Check for required utilities and directories.
    Run BUSCO on each FASTA file in parallel using 4 cores.
    Generate a summary table of the BUSCO results and save it as busco_summary.csv in the output directory.