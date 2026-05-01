## K-mer Genome Sequence Analyzer

The script `kmer_analyzer.py` is created to analyze a DNA sequence and give you information about the kmers and the next character in the sequence.

To run `kmer_analyzer.py` in the shell terminal use:
  
  `python kmer_analyzer.py <input_file> <k> <output_file>`

Example:
  `python kmer_analyzer.py dna.txt 2 output.txt`

#### Input format:
Each line in the input file should contain a DNA sequence.

Valid characters: A, C, G, T (uppercase only)

Example:
ACTGTCGATA

#### Output format:
The output file contains each k-mer followed by the nucleotide and frequency of the next character.

Format:
  KMER A:count C:count G:count T:count

Example:
  GAC T:1
  CTA G:2

#### Testing:
This project was tested using pytest. Tests for each function in `kmer_analyzer.py` can be found in the script `test_kmer_analyzer.py`.

#### AI Statement: AI was used to help generate test ideas and code for each function and code on how to fix errors exposed by tests. AI was also used to write the docstrings for each function.
