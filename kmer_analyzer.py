import sys

def validate_sequence(sequence, k):
    """
    Validate a DNA sequence for k-mer analysis.

    A sequence is considered valid if:
    - Its length is at least k
    - It contains only valid DNA nucleotides (A, C, G, T)

    Parameters:
        sequence (str): DNA sequence to validate
        k (int): Length of k-mers

    Returns:
        bool: True if the sequence is valid, False otherwise
    """
    if len(sequence) < k:
        return False
    #edit this section to only accept capital letters ACGT
    for nucleotide in sequence:
        if nucleotide not in "ACGT":
            return False

    return True

def update_kmer_count(kmer_data, kmer, next_char):
    """
    Update k-mer frequency data with a new observation.

    Tracks:
    - Total count of each k-mer
    - Frequency of characters that follow each k-mer

    Parameters:
        kmer_data (dict): Dictionary storing k-mer counts and context
        kmer (str): The k-mer substring
        next_char (str): The character that follows the k-mer

    Returns:
        dict: Updated kmer_data dictionary
    """

    # change starting count to 0 so doesn't double count first time kmer is seen
    if kmer not in kmer_data:
        kmer_data[kmer] = {'count': 0, 'next_chars': {}}
    
    kmer_data[kmer]['count'] += 1
    
    if next_char not in kmer_data[kmer]['next_chars']:
        kmer_data[kmer]['next_chars'][next_char] = 0
    kmer_data[kmer]['next_chars'][next_char] += 1

    return kmer_data

def count_kmers_with_context(sequence, k):
    """
    Generate k-mer counts and their subsequent character frequencies.

    Iterates through the sequence to extract overlapping k-mers and
    records how often each k-mer appears and what character follows it.
    The final k-mer is included even if it has no following character.

    Parameters:
        sequence (str): DNA sequence
        k (int): Length of k-mers

    Returns:
        dict: Dictionary containing k-mer counts and next character frequencies
    """

    kmer_data = {}

    # chnage to k +1 to include overlapping k-mers and last k-mer
    for i in range(len(sequence) - k + 1):
        kmer = sequence[i:i+k]

        # if there IS a next character
        if i + k < len(sequence):
            next_char = sequence[i + k]
            kmer_data = update_kmer_count(kmer_data, kmer, next_char)

        # if this is the last k-mer (no next character)
        else:
            if kmer not in kmer_data:
                kmer_data[kmer] = {'count': 0, 'next_chars': {}}
            kmer_data[kmer]['count'] += 1

    return kmer_data


def write_results_to_file(kmer_data, output_filename):
    """
    Write k-mer results to a file in sorted order.

    Each line of the output file contains:
    - A k-mer
    - The frequencies of characters that follow it (formatted as char:count)

    Parameters:
        kmer_data (dict): Dictionary of k-mer counts and next character frequencies
        output_filename (str): Path to the output file

    Returns:
        None
    """
    sorted_kmers = sorted(kmer_data.keys())
    
    with open(output_filename, 'w') as f:
        for kmer in sorted_kmers:
            next_chars = kmer_data[kmer]['next_chars']
            
            next_char_str = " ".join(
                f"{char}:{freq}" 
                for char, freq in sorted(next_chars.items())
            )
            
            f.write(f"{kmer} {next_char_str}\n")


def main():
    """
    Main driver function for k-mer analysis.

    Reads DNA sequences from an input file, validates them, computes k-mer
    counts and context, and writes the aggregated results to an output file.

    Command-line arguments:
        sys.argv[1]: Input file containing DNA sequences (one per line)
        sys.argv[2]: Integer value of k (k-mer length)
        sys.argv[3]: Output file path

    Returns:
        None
    """
    sequence_file = sys.argv[1]
    k = int(sys.argv[2])
    output_file = sys.argv[3]
    
    print(f"Reading sequences from {sequence_file}...")
    
    kmer_data = {}

    with open(sequence_file, 'r') as f:
        for sequence in f:
            sequence = sequence.strip().upper()
            
            # validate and give reason for why sequence is incorrect 
            if not validate_sequence(sequence, k):
                if len(sequence) < k:
                    print(f"  Warning: Sequence length is shorter than k (k={k}), skipping.")
                else:
                    print(f"  Warning: Invalid characters in sequence, skipping.")
                continue
            
            new_data = count_kmers_with_context(sequence, k)
            
            # merge results
            for kmer in new_data:
                if kmer not in kmer_data:
                    kmer_data[kmer] = new_data[kmer]
                else:
                    kmer_data[kmer]['count'] += new_data[kmer]['count']
                    
                    for char, freq in new_data[kmer]['next_chars'].items():
                        kmer_data[kmer]['next_chars'][char] = (
                            kmer_data[kmer]['next_chars'].get(char, 0) + freq
                        )

    # write ONCE after loop
    write_results_to_file(kmer_data, output_file)


if __name__ == "__main__":
    main()
            
  
