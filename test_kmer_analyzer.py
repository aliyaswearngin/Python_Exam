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

#basic function
def test_basic_kmers():
    result = count_kmers_with_context("ATGTGA", 2)

    assert result["AT"]["count"] == 1
    assert result["TG"]["count"] == 2
    assert result["GT"]["count"] == 1

# Overlapping k-mers 
def test_overlapping_kmers():
    result = count_kmers_with_context("ATAT", 2)
    
    assert result["AT"]["count"] == 2
    assert result["TA"]["count"] == 1

#Repeated patterns
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

#test is recognize last kmer if no next character    
def test_last_kmer_handling():
    result = count_kmers_with_context("ATG", 2)

    assert "TG" in result 
    
#################################
## write_results_to_file tests ##
#################################
#Takes your processed k-mer dictionary and writes it to a file


#test basic funtion
def test_write_results_basic(tmp_path):
    data = {
        "AT": {"count": 2, "next_chars": {"G": 1, "C": 1}}
    }

    file = tmp_path / "output.txt"

    write_results_to_file(data, file)

    content = file.read_text()

    assert "AT C:1 G:1" in content

#test how sorts multiple inputs 
def test_write_results_sorted(tmp_path):
    data = {
        "TG": {"count": 1, "next_chars": {"A": 1}},
        "AT": {"count": 1, "next_chars": {"G": 1}}
    }

    file = tmp_path / "output.txt"
    write_results_to_file(data, file)

    content = file.read_text().splitlines()

    assert content[0].startswith("AT")
    assert content[1].startswith("TG")

# test how handles multiple next characters  
def test_write_multiple_next_chars(tmp_path):
    data = {
        "AT": {"count": 2, "next_chars": {"G": 2, "C": 1}}
    }

    file = tmp_path / "output.txt"
    write_results_to_file(data, file)

    content = file.read_text()

    assert "G:2" in content
    assert "C:1" in content
    
################
## main tests ##
################

#test how writes files for multiple sequences 
import sys
def test_main_multiple_sequences(tmp_path, monkeypatch):
    # create input file with TWO sequences
    input_file = tmp_path / "input.txt"
    input_file.write_text("ATG\nTGA\n")

    output_file = tmp_path / "output.txt"

    # fake command-line arguments
    monkeypatch.setattr(sys, "argv", [
        "script.py",
        str(input_file),
        "2",
        str(output_file)
    ])

    # run main
    main()

    content = output_file.read_text()

    # Should contain results from BOTH sequences
    assert "AT" in content
    assert "TG" in content

#verify correct warning when kmer is longer than sequence 
def test_warning_when_k_larger_than_sequence(tmp_path, monkeypatch, capfd):
    import sys

    input_file = tmp_path / "input.txt"
    input_file.write_text("ATG\n")

    output_file = tmp_path / "output.txt"

    monkeypatch.setattr(sys, "argv", [
        "script.py",
        str(input_file),
        "10",   # k is larger than sequence
        str(output_file)
    ])

    main()

    out, err = capfd.readouterr()

    assert "Sequence length is shorter than k" in out
