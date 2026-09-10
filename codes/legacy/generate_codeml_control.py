#!/usr/bin/env python3

import argparse
from datetime import datetime
from pathlib import Path

def now_stamp():
    return datetime.now().strftime("%Y%m%d%H%M%S")

def create_codeml_control_file(seqfile, treefile, outfile, ndata=1, model=0, ns_sites='0', fix_kappa=0, cleandata=0):
    sites = ' '.join(ns_sites.split(','))
    ctrl = f"""
          seqfile = {seqfile} * sequence data filename
         treefile = {treefile} * tree structure file name
          outfile = {outfile}.cdml * main result file name
     
            noisy = 3   * 0,1,2,3,9: how much rubbish on the screen
          verbose = 0   * 1: detailed output, 0: concise output 
          runmode = 0   * 0: user tree; 1: semi-automatic; 2: automatic
                        * 3: StepwiseAddition; (4,5):PerturbationNNI; ...
    
          seqtype = 1   * 1:codons; 2:AAs; 3:codons-->AAs
        CodonFreq = 2   * 0:1/61 each, 1:F1X4, 2:F3X4, 3:codon table
            clock = 0   * 0:no clock, 1:clock; 2:local clock
    
            ndata = {ndata}   * number of gene alignments to be analysed
            model = {model}   * models for codons: 0:one, 1:b, 2:2 or more dN/dS
                        * ratios for branches
          NSsites = {sites}   * 0:one w;1:neutral;2:selection;
                        * 3:discrete;4:freqs;
                        * 5:gamma;6:2gamma;7:beta;8:beta&w;9:beta&gamma;
                        
            icode = 0   * 0:universal code; 1:mammalian mt; 2-10:see below
 
        fix_omega = 0   * 1: omega or omega_1 fixed, 0: estimate
            omega = .4  * initial or fixed omega for codons
            
        cleandata = {cleandata}   * remove sites with ambiguity data (1:yes, 0:no)?
    """
    return ctrl

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('seq_file', help="INPUT:  aligned_cds_sequences.pml", type=str)
    parser.add_argument('tree_file', help="INPUT:  newick format tree file of species", type=str)
    parser.add_argument('-t', '--type', choices=['zero_model', 'all_sites', 'reduced_all_sites', 'free_model'], required=True, help="Type of control file to generate")
    parser.add_argument('-g', '--genes', default=1, help="INPUT:  # of genes to be analyzed. Default: 1", type=int)
    parser.add_argument('--ns_sites', default='0', help="NSsites values for reduced model (e.g., '0,1,2' or '0,7,8'). Required if type is 'reduced_all_sites'.", type=str)
    parser.add_argument('-o', '--output_file', help="OUTPUT:  Path to save the generated control file", type=str)

    arguments = parser.parse_args()

    if not Path(arguments.seq_file).exists():
        print(f"Sequence file <{arguments.seq_file}> does not exist")
        exit(1)
    if not Path(arguments.tree_file).exists():
        print(f"Tree file <{arguments.tree_file}> does not exist")
        exit(1)

    SEQFILE = arguments.seq_file
    TREEFILE = arguments.tree_file
    OUTFILE = f"mlc_{now_stamp()}"
    NDATA = arguments.genes

    if arguments.type == 'zero_model':
        MODEL = 0
        NSSTIES = '0'
        CLEANDATA = 0

    elif arguments.type == 'all_sites':
        MODEL = 0
        NSSTIES = '0,1,2,7,8'
        CLEANDATA = 0

    elif arguments.type == 'reduced_all_sites':
        MODEL = 0
        NSSTIES = arguments.ns_sites
        CLEANDATA = 0

    elif arguments.type == 'free_model':
        MODEL = 1
        NSSTIES = '0'
        CLEANDATA = 1

    control_file = create_codeml_control_file(SEQFILE, TREEFILE, OUTFILE, NDATA, MODEL, NSSTIES, cleandata=CLEANDATA)

    if arguments.output_file:
        with open(arguments.output_file, 'w') as f:
            f.write(control_file)
        print(f"Control file saved to {arguments.output_file}")
    else:
        print(control_file)

if __name__ == "__main__":
    main()
