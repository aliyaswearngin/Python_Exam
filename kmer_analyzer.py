import sys

def validate_sequence(sequence, k):
    if len(sequence) < k:
        return False
#edit this section to only accept capital letters ACGT
    for nucleotide in sequence:
        if nucleotide not in "ACGT":
            return False

    return True

def update_kmer_count(kmer_data, kmer, next_char):
    # change starting count to zero so doesn't double count first time kmer appears
    if kmer not in kmer_data:
        kmer_data[kmer] = {'count': 0, 'next_chars': {}}
    
    kmer_data[kmer]['count'] += 1
    
    if next_char not in kmer_data[kmer]['next_chars']:
        kmer_data[kmer]['next_chars'][next_char] = 0
    kmer_data[kmer]['next_chars'][next_char] += 1

    return kmer_data

def count_kmers_with_context(sequence, k):
    kmer_data = {}
    #change to k + 1 so range includes overlapping and last kmer
    for i in range(len(sequence) - k + 1):
        kmer = sequence[i:i+k]
        # add line to say get next_ char if i + k is in range of sequence  
        if i + k < len(sequence):
          next_char = sequence[i+k]
          kmer_data = update_kmer_count(kmer_data, kmer, next_char)
        #if not in range, then last k mer should stil be counted
        else:
            # get last k-mer even if no next character
            if kmer not in kmer_data:
                kmer_data[kmer] = {'count': 0, 'next_chars': {}}
            kmer_data[kmer]['count'] += 1
    
    return kmer_data


def write_results_to_file(kmer_data, output_filename):
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
    sequence_file = sys.argv[1]
    k = int(sys.argv[2])
    output_file = sys.argv[3]
    
    print(f"Reading sequences from {sequence_file}...")

    with open(sequence_file, 'r') as f:
        for sequence in f:
            sequence = sequence.strip()

            if not validate_sequence(sequence, k):
                print(f"  Warning: Skipping sequence")
                continue
            
            kmer_data = count_kmers_with_context(sequence, k) 
            
            write_results_to_file(kmer_data, output_file)

if __name__ == '__main__':
    main()
