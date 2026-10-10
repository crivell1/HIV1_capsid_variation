# Save as src/count_p24_frequency.py
import csv
from collections import Counter
from io import StringIO
from Bio import SeqIO

input_file = "align100_p24.efa"
output_file = "p24_frequencies.csv"
amino_acids = set("ACDEFGHIKLMNPQRSTVWY")

# Extract the first alignment from the EFA file.
fasta_lines = []
with open(input_file) as file:
    for line in file:
        if line.startswith("<"):
            if fasta_lines:
                break
            continue
        fasta_lines.append(line)

sequences = list(SeqIO.parse(StringIO("".join(fasta_lines)), "fasta"))
if not sequences:
    raise ValueError("No sequences found.")

lengths = {len(record.seq) for record in sequences}
if len(lengths) != 1:
    raise ValueError("Aligned sequences have different lengths.")

with open(output_file, "w", newline="") as file:
    writer = csv.writer(file)
    writer.writerow(["alignment_column", "amino_acid", "count", "valid_sequences", "frequency_percent"])

    for column in range(len(sequences[0].seq)):
        counts = Counter(
            str(record.seq[column]).upper()
            for record in sequences
            if str(record.seq[column]).upper() in amino_acids
        )
        total = sum(counts.values())

        for amino_acid, count in sorted(counts.items()):
            writer.writerow([column + 1, amino_acid, count, total, round(100 * count / total, 2)])

print(f"Counted {len(sequences)} sequences. Wrote {output_file}")
