#!/bin/bash
#SBATCH --job-name=alignment_viz
#SBATCH --output=/home/crivell1/HIV_qb_proj/alignment_viz.out
#SBATCH --error=/home/crivell1/HIV_qb_proj/alignment_viz.err
#SBATCH --time=08:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4 
#SBATCH --mem=8G

#I'm putting 4 CPUs because the MUSCLE alignment tool can use multiple threads for faster processing

set -euo pipefail

export MPLBACKEND=Agg # Google said to use this to avoid the "RuntimeError: Invalid DISPLAY variable" error when running matplotlib on HPC

echo "Job ID: $SLURM_JOB_ID"
echo "Array Task ID: $SLURM_JOB_NAME"
echo "Started on: $(hostname)"
echo "Start time: $(date)"
echo "Working directory: $(pwd)"

# Load required environment modules
module anaconda
conda activate hiv_p24
module load muscle/5.1

# Set up project directories
PROJECT_DIR="/home/crivell1/HIV_qb_proj"
DATA_DIR="$PROJECT_DIR/data"
SRC_DIR="$PROJECT_DIR/src"

#change file names before running this script
python $SRC_DIR/create_p24_fa.py $DATA_DIR/gag.fasta $DATA_DIR/p24.fasta # this will create the p24 fasta from the gag fasta file

muscle -align $DATA_DIR/p24.fasta -stratified -output $DATA_DIR/align_all_sequence.efa -threads "$SLURM_CPUS_PER_TASK" # this will align the p24 sequences and output the alignment in EFA format
muscle -maxcc $DATA_DIR/align_all_sequence.efa -output $DATA_DIR/best_alignment.fasta # this will output the best sequence in FASTA format

python $SRC_DIR/count_p24_frequency.py $DATA_DIR/best_alignment.fasta $DATA_DIR/p24_frequency.csv # this will count the frequency of each p24 sequence and output the results in a CSV file
python $SRC_DIR/entropy_calc.py $DATA_DIR/p24_frequency.csv $DATA_DIR/p24_entropy.csv # this will calculate the entropy of each p24 sequence and output the results in a CSV file
python $SRC_DIR/entropy_viz.py $DATA_DIR/p24_entropy.csv $DATA_DIR/p24_entropy_visualization.png # this will visualize the entropy of each p24 sequence amino acid and output the results in a PNG file

python $SRC_DIR/count_p24_frequency_by_variant.py $DATA_DIR/best_alignment.fasta $DATA_DIR/p24_frequency_by_variant.csv # this will count the frequency of each p24 sequence and output the results in a CSV file
python $SRC_DIR/entropy_calc_by_variant.py $DATA_DIR/p24_frequency_by_variant.csv $DATA_DIR/p24_entropy_by_variant.csv # this will calculate the entropy of each p24 sequence and output the results in a CSV file
python $SRC_DIR/entropy_heatmap.py $DATA_DIR/p24_entropy_by_variant.csv $DATA_DIR/p24_entropy_visualization_by_variant.png # this will visualize the entropy of each p24 sequence amino acid and output the results in a PNG file

# Full MSA visualization might be huge...
python $SRC_DIR/MSA_viz.py $DATA_DIR/best_alignment.fasta $DATA_DIR/msa_visualization.png # this will visualize the multiple sequence alignment and output the results in a PNG file
