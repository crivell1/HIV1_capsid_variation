#!/bin/bash
##SBATCH --job-name=alignment_viz
##SBATCH --output=/home/crivell1/HIV_qb_proj/alignment_viz.out
##SBATCH --error=/home/crivell1/HIV_qb_proj/alignment_viz.err
##SBATCH --time=08:00:00
##SBATCH --nodes=1
##SBATCH --ntasks-per-node=1
##SBATCH --cpus-per-task=1
##SBATCH --mem=8G

set -euo pipefail

echo "Job ID: $SLURM_JOB_ID"
echo "Array Task ID: $SLURM_JOB_NAME"
echo "Started on: $(hostname)"
echo "Start time: $(date)"
echo "Working directory: $(pwd)"

# Load required modules
ml python/3.12.3

#
DATA_DIR = "/home/crivell1/HIV_qb_proj/data"

#change file names before running this script
python create_p24_fa.py $DATA_DIR/gag.fasta $DATA_DIR/p24.fasta # this will create the p24 fasta from the gag fasta file

muscle -align p24.fasta -stratified -output align_all_sequence.efa # this will align the p24 sequences and output the alignment in EFA format
muscle -maxcc align_all_sequence.efa -output best_sequence.fasta # this will output the best sequence in FASTA format

## This will need to be changed to take in arguments for the output files, so that it can be run from the command line
python MSA_viz.py align_all_sequence.efa $DATA_DIR/msa_visualization.png # this will visualize the multiple sequence alignment and output the results in a PNG file

## the count_p24_frequency.py needs to be changed to take in arguments for the input and output files, so that it can be run from the command line
python count_p24_frequency.py $DATA_DIR/p24.fasta $DATA_DIR/p24_frequency.csv # this will count the frequency of each p24 sequence and output the results in a CSV file
python entropy_calc.py $DATA_DIR/p24.fasta $DATA_DIR/p24_entropy.csv # this will calculate the entropy of each p24 sequence and output the results in a CSV file
python entropy_viz.py $DATA_DIR/p24_entropy.csv $DATA_DIR/p24_entropy_visualization.png # this will visualize the entropy of each p24 sequence and output the results in a PNG file

python count_p24_frequency_by_variant.py $DATA_DIR/p24.fasta $DATA_DIR/p24_frequency_by_variant.csv # this will count the frequency of each p24 sequence and output the results in a CSV file
python entropy_calc_by_variant.py $DATA_DIR/p24.fasta $DATA_DIR/p24_entropy_by_variant.csv # this will calculate the entropy of each p24 sequence and output the results in a CSV file
python entropy_heatmap.py $DATA_DIR/p24_entropy_by_variant.csv $DATA_DIR/p24_entropy_visualization_by_variant.png # this will visualize the entropy of each p24 sequence and output the results in a PNG file
