import csv
from collections import Counter
from io import StringIO
from Bio import SeqIO
from collections import defaultdict, Counter

input_file = "align100_p24.efa"
output_file = "p24_frequencies_by_variant.csv"
amino_acids = set("ACDEFGHIKLMNPQRSTVWY")

# Extract the first alignment from the EFA file.
fasta_lines = []
with open(input_file) as file:
    for line in file:
        if line.startswith("<"):
            if fasta_lines: # "if the list contains any lines, break" - only keeps first alignment
                break
            continue
        fasta_lines.append(line)

sequences = list(SeqIO.parse(StringIO("".join(fasta_lines)), "fasta")) # Parse the first alignment as FASTA sequences
if not sequences:
    raise ValueError("No sequences found.")

lengths = {len(record.seq) for record in sequences}
if len(lengths) != 1:
    raise ValueError("Aligned sequences have different lengths.")

# Group sequences using the first field of each FASTA header.
variants = defaultdict(list) #makes each variant a key in a dictionary
for record in sequences:
    variant = record.id.split(".")[0]
    variants[variant].append(record)

with open(output_file, "w", newline="") as file: #making the new csv file with variant info
    writer = csv.writer(file)
    writer.writerow([
        "variant", "alignment_column", "amino_acid",
        "count", "valid_sequences", "frequency_percent",
        "group_sequences"
    ])
    for variant, records in sorted(variants.items()):
            for column in range(len(records[0].seq)):
                counts = Counter(
                    str(record.seq[column]).upper()
                    for record in records
                    if str(record.seq[column]).upper() in amino_acids
                )
                total = sum(counts.values()) ## Fancy way to make the dictionary counter similar to our bootcamp project - just using Counter function

                # Preserve positions with no valid amino acids as missing data.
                if total == 0:
                    writer.writerow([
                        variant, column + 1, "", 0, 0, "", len(records)
                    ])
                    continue

                for amino_acid, count in sorted(counts.items()):
                    writer.writerow([
                        variant, column + 1, amino_acid,
                        count, total, 100 * count / total, len(records) # Frequency percent and group size
                    ])

print(f"Counted {len(sequences)} sequences. Wrote {output_file}")