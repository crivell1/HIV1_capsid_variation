#!/usr/bin/env python3

import sys
from io import StringIO
from Bio import AlignIO
from pymsaviz import MsaViz, get_msa_testdata

## This is adapted from the documentation of pymsaviz https://www.biostars.org/p/9545724/
def msa_viz(file_path):
    align_list = []
    with open(file_path) as f:
        for line in f:
            if line.startswith("<"):
                if align_list:
                    break  # Stop when the second alignment begins
                continue
            align_list.append(line)
    alignment = AlignIO.read(StringIO("".join(align_list)), "fasta")
    mv = MsaViz(alignment, wrap_length=250, show_count=True)
    mv.savefig("trial1_msaviz.png")


if __name__ == "__main__":
    msa_viz(sys.argv[1])
