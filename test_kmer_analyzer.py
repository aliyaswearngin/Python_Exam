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
    
#test for non letters
def test_validate_rejects_numbers():
    result = validate_sequence("AT1GC", 2)
    assert result is False

#test non DNA letters 
def test_validate_rejects_non_dna_letters():
    result = validate_sequence("ANGC", 3) 
    assert result is False

#test lowercase letters
def test_validate_rejects_lower_case_letters():
    result = validate_sequence("aggc", 3) 
    assert result is False
    
#############################
## update_kmer_count tests ##
#############################
#how many times kmer appears
#what character comes after kmer in the sequence
# creates a new entry for new k-mer found-# double counts the first time-need to fix this

# test what first insertion of kmer returns 
def test_update_kmer_first_insert():
    data = {}
    result = update_kmer_count(data, "AT", "G")

    assert result["AT"]["count"] == 1
    assert result["AT"]["next_chars"]["G"] == 1

#use same kmer twice, but different next character
def test_update_kmer_multiple_updates():
    data = {}

    update_kmer_count(data, "AT", "G")
    update_kmer_count(data, "AT", "C")

    assert data["AT"]["count"] == 2
    assert data["AT"]["next_chars"]["G"] == 1
    assert data["AT"]["next_chars"]["C"] == 1

#same kmer and same next character
def test_update_kmer_repeated_next_char():
    data = {}

    update_kmer_count(data, "AT", "G")
    update_kmer_count(data, "AT", "G")

    assert data["AT"]["count"] == 2
    assert data["AT"]["next_chars"]["G"] == 2

# test multiple different kmers 
def test_update_multiple_kmers():
    data = {}

    update_kmer_count(data, "AT", "G")
    update_kmer_count(data, "TG", "A")

    assert "AT" in data
    assert "TG" in data
    assert data["AT"]["count"] == 1
    assert data["TG"]["count"] == 1


