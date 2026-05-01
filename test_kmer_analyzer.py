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

####################################
## count_kmers_with_context tests ##
####################################

#basic function-pass
def test_basic_kmers():
    result = count_kmers_with_context("ATGTGA", 2)

    assert result["AT"]["count"] == 1
    assert result["TG"]["count"] == 2
    assert result["GT"]["count"] == 1

# Overlapping k-mers - failed
def test_overlapping_kmers():
    result = count_kmers_with_context("ATAT", 2)
    
    assert result["AT"]["count"] == 2
    assert result["TA"]["count"] == 1

#Repeated patterns- failed
def test_repeated_pattern():
    result = count_kmers_with_context("AAAA", 2)

    assert result["AA"]["count"] == 3  

#Next-character tracking correctness
def test_next_char_tracking():
    result = count_kmers_with_context("ATGTGA", 2)

    assert result["TG"]["count"] == 2
    assert result["TG"]["next_chars"]["T"] == 1
    assert result["TG"]["next_chars"]["A"] == 1
    
# Edge case: sequence too short
def test_short_sequence():
    result = count_kmers_with_context("A", 2)

    assert result == {}

#test is recognize last kmer if no next character - failed   
def test_last_kmer_handling():
    result = count_kmers_with_context("ATG", 2)

    assert "TG" in result 
