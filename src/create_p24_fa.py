#! bin/env/ python3

import sys
import os

# I'm sure there's a FASTA reader that could make this faster - but hopefully text parsing on 120k sequences isn't too bad...

#Make a dictionary with the FASTA file, then make a smaller dictionary with just the CA
def pull_p24_seqs(file, aa_range):
    all_seqs = {}
    for line in open(file, "r"):
        if line.startswith(">"):
            seq_name = line.strip()[1:]
            all_seqs[seq_name] = ""
        else:
            all_seqs[seq_name] += line.strip() 
    p24_seqs = {}
    for key, value in all_seqs.items():
        p24_seqs[key] = value[aa_range].strip("-") ## since they come prealligned as a whole to Hxb2 and we are excising a small region, we need to remove the "gap" notations and rerun our alignment - we will do a MUSCLE alignment to "refresh" our alignemnts as an unsupervised MSA

    return p24_seqs

def write_fasta(dict, new_file_name):
    fs = open(new_file_name, "w")
    for key, value in dict.items():
        fs.write(f'>{key}\n{value}\n')
    fs.close()

if __name__ == "__main__":
    ## code will be run like: python create_p24_fa.py old_file new_file
    whole_file = sys.argv[1]
    p24_range = slice(132, 366) ## LANL uses HXB2 coordinates for annotations - this is pulled directly from that - any gaps get condensed and our MSA gets rerun on the new fasta made
    all_seqs = pull_p24_seqs(whole_file, p24_range)
    write_fasta(all_seqs, sys.argv[2])
