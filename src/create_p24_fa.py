import sys
import os

def read_fasta(file):
    """Read FASTA file into dictionary."""
    all_seqs = {}
    seq_name = None
    for line in open(file, "r"):
        line = line.strip()
        if line.startswith(">"):
            seq_name = line[1:]
            all_seqs[seq_name] = ""
        elif seq_name is not None:
            all_seqs[seq_name] += line
    return all_seqs

def find_p24_boundaries_from_hxb2(all_seqs, p24_start_aa=133, p24_end_aa=363):
    """
    Find alignment columns corresponding to p24 using HXB2 
    amino acid numbering.
    """
    # Find HXB2
    hxb2_key = None
    for key in all_seqs:
        if 'HXB2' in key.upper() or 'B.FR.1983' in key:
            hxb2_key = key
            break
    
    if hxb2_key is None:
        print("ERROR: HXB2 not found!")
        print("First 5 sequence names:")
        for i, key in enumerate(all_seqs.keys()):
            print(f"  {key}")
            if i >= 4:
                break
        sys.exit(1)
    
    hxb2_seq = all_seqs[hxb2_key]
    aa_count = 0
    start_col = None
    end_col = None
    
    for col_idx, char in enumerate(hxb2_seq):
        if char not in ['-', '.']:
            aa_count += 1
            if aa_count == p24_start_aa and start_col is None:
                start_col = col_idx
            if aa_count == p24_end_aa:
                end_col = col_idx + 1
                break
    
    print(f"Found HXB2: {hxb2_key}")
    print(f"Alignment columns for p24: {start_col} to {end_col}")
    
    # Verify by showing extracted HXB2 p24
    extracted = hxb2_seq[start_col:end_col].replace('-', '')
    print(f"HXB2 p24 length (no gaps): {len(extracted)}")
    print(f"HXB2 p24 starts with: {extracted[:20]}")
    print(f"Expected start: PIVQNLQGQMVHQAISPR...")
    
    return start_col, end_col

def extract_p24(all_seqs, start_col, end_col, remove_gaps=True):
    """Extract p24 region from all sequences."""
    p24_seqs = {}
    skipped = 0
    
    for key, value in all_seqs.items():
        extracted = value[start_col:end_col]
        
        if remove_gaps:
            extracted = extracted.replace('-', '').replace('.', '')
        
        # Skip empty or very short sequences
        if len(extracted) < 200:  # p24 should be ~231 aa
            skipped += 1
            continue
            
        p24_seqs[key] = extracted
    
    print(f"\nExtracted {len(p24_seqs)} sequences")
    print(f"Skipped {skipped} (too short after extraction)")
    
    return p24_seqs

def write_fasta(seq_dict, new_file_name):
    """Write dictionary to FASTA file."""
    with open(new_file_name, "w") as fs:
        for key, value in seq_dict.items():
            fs.write(f'>{key}\n{value}\n')

if __name__ == "__main__":
    whole_file = sys.argv[1]
    output_file = sys.argv[2]
    
    # Step 1: Read all sequences
    print("Reading FASTA file...")
    all_seqs = read_fasta(whole_file)
    print(f"Total sequences: {len(all_seqs)}")
    
    # Step 2: Find p24 boundaries using HXB2
    print("\nFinding p24 boundaries...")
    start_col, end_col = find_p24_boundaries_from_hxb2(all_seqs)
    
    # Step 3: Extract p24 from all sequences
    print("\nExtracting p24...")
    p24_seqs = extract_p24(all_seqs, start_col, end_col, remove_gaps=True)
    
    # Step 4: Write output
    write_fasta(p24_seqs, output_file)
    print(f"\nOutput written to: {output_file}")

