import pytest

from kmer_analyzer import validate_sequence, update_kmer_count, count_kmers_with_context, write_results_to_file, main

#############################
## validate_sequence tests ##
#############################

#test basics
def test_validate_sequences_basic():
  result = validate_sequence("ATGC", 2)
  assert result is True

#test edge cases-too long 
def test_validate_sequence_k_too_large():
    result = validate_sequence("AT", 5)
    assert result is False

#empty sequence
def test_validate_empty_sequence():
    result = validate_sequence("", 2) 
    assert result is False

#test for non numbers
def test_validate_rejects_numbers():
    result = validate_sequence("AT1GC", 2)
    assert result is False

#test non DNA letters outside of kmer length
#failed - need to fix function to only recognize letters ATCG
def test_validate_rejects_non_dna_letters():
    result = validate_sequence("ANGC", 3) 
    assert result is False
