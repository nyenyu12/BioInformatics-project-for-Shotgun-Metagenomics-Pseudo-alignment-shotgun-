import subprocess
import pytest
import os
import json
import gzip
import pickle
from pathlib import Path


# Helper functions
def run_command(command):
    """Runs a subprocess command and returns output."""
    result = subprocess.run(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    return result.stdout, result.stderr, result.returncode


@pytest.fixture
def test_data_dir(tmp_path):
    """Creates a temporary directory with test data."""
    genome_file = tmp_path / "test_genome.fa"
    fastq_file = tmp_path / "test_reads.fq"
    reference_file = tmp_path / "test_reference.kdb"
    align_file = tmp_path / "test_alignment.aln"

    # Create a sample genome file
    genome_data = ">Genome1\nAGCTAGCTAGCTAGCTAGCT\n>Genome2\nTGCATGCATGCATGCATGCA\n"
    genome_file.write_text(genome_data)

    # Create a sample reads file
    reads_data = "@Read1\nAGCTAGCT\n+\nIIIIIIII\n@Read2\nTGCATGCA\n+\nIIIIIIII\n"
    fastq_file.write_text(reads_data)

    return {
        "genome_file": genome_file,
        "fastq_file": fastq_file,
        "reference_file": reference_file,
        "align_file": align_file,
    }


def test_create_reference(test_data_dir):
    """Tests creating a k-mer reference database."""
    command = [
        "python3", "main.py",
        "-t", "reference",
        "-g", str(test_data_dir["genome_file"]),
        "-k", "4",
        "-r", str(test_data_dir["reference_file"])
    ]
    stdout, stderr, returncode = run_command(command)

    assert returncode == 0, f"Error creating reference: {stderr}"
    assert test_data_dir["reference_file"].exists(), "Reference file was not created."

    # Check that the reference file is a valid pickle
    with gzip.open(test_data_dir["reference_file"], "rb") as f:
        kmer_reference = pickle.load(f)
    assert kmer_reference is not None, "Failed to load reference from file."


def test_dump_reference(test_data_dir):
    """Tests dumping a reference database to JSON."""
    # Ensure reference is created first
    test_create_reference(test_data_dir)

    command = [
        "python3", "main.py",
        "-t", "dumpref",
        "-r", str(test_data_dir["reference_file"])
    ]
    stdout, stderr, returncode = run_command(command)

    assert returncode == 0, f"Error dumping reference: {stderr}"
    
    # Check that output is valid JSON
    try:
        json_output = json.loads(stdout)
        assert "Summary" in json_output, "Reference dump output does not contain expected summary."
    except json.JSONDecodeError:
        pytest.fail("Reference dump output is not valid JSON.")


def test_align_reads(test_data_dir):
    """Tests aligning reads using a pre-built reference."""
    # Ensure reference is created first
    test_create_reference(test_data_dir)

    command = [
        "python3", "main.py",
        "-t", "align",
        "-r", str(test_data_dir["reference_file"]),
        "--reads", str(test_data_dir["fastq_file"]),
        "-a", str(test_data_dir["align_file"])
    ]
    stdout, stderr, returncode = run_command(command)

    assert returncode == 0, f"Error aligning reads: {stderr}"
    assert test_data_dir["align_file"].exists(), "Alignment file was not created."

    # Check that alignment file is a valid pickle
    with gzip.open(test_data_dir["align_file"], "rb") as f:
        pseudo_alignment = pickle.load(f)
    assert pseudo_alignment is not None, "Failed to load alignment from file."


def test_dump_alignment(test_data_dir):
    """Tests dumping an alignment file."""
    # Ensure alignment is created first
    test_align_reads(test_data_dir)

    command = [
        "python3", "main.py",
        "-t", "dumpalign",
        "-a", str(test_data_dir["align_file"])
    ]
    stdout, stderr, returncode = run_command(command)

    assert returncode == 0, f"Error dumping alignment: {stderr}"
    
    # Check that output is valid JSON
    try:
        json_output = json.loads(stdout)
        assert "Statistics" in json_output, "Alignment dump output does not contain expected statistics."
    except json.JSONDecodeError:
        pytest.fail("Alignment dump output is not valid JSON.")


def test_build_reference_align_and_dump(test_data_dir):
    """Tests building a reference, aligning reads, and dumping results in one step."""
    command = [
        "python3", "main.py",
        "-t", "dumpalign",
        "-g", str(test_data_dir["genome_file"]),
        "-k", "4",
        "--reads", str(test_data_dir["fastq_file"])
    ]
    stdout, stderr, returncode = run_command(command)

    assert returncode == 0, f"Error in build-reference-align-dump: {stderr}"
    
    # Check that output is valid JSON
    try:
        json_output = json.loads(stdout)
        assert "Statistics" in json_output, "Build-reference-align-dump output does not contain expected statistics."
    except json.JSONDecodeError:
        pytest.fail("Build-reference-align-dump output is not valid JSON.")


def test_invalid_file_handling(test_data_dir):
    """Tests handling of missing or invalid files."""
    invalid_file = test_data_dir["genome_file"].parent / "invalid_file.fa"

    # Test missing genome file
    command = [
        "python3", "main.py",
        "-t", "reference",
        "-g", str(invalid_file),
        "-k", "4",
        "-r", str(test_data_dir["reference_file"])
    ]
    stdout, stderr, returncode = run_command(command)
    assert returncode != 0, "Reference creation should fail with missing genome file."
    assert "does not exist" in stderr, "Error message should mention missing file."

    # Test invalid extension handling
    invalid_ext_file = test_data_dir["genome_file"].parent / "invalid_file.txt"
    invalid_ext_file.touch()

    command = [
        "python3", "main.py",
        "-t", "reference",
        "-g", str(invalid_ext_file),
        "-k", "4",
        "-r", str(test_data_dir["reference_file"])
    ]
    stdout, stderr, returncode = run_command(command)
    assert returncode != 0, "Reference creation should fail with invalid extension."
    assert "Invalid file extension" in stderr, "Error message should indicate invalid extension."

    # Test missing reference file for alignment
    command = [
        "python3", "main.py",
        "-t", "align",
        "-r", str(invalid_file),
        "--reads", str(test_data_dir["fastq_file"]),
        "-a", str(test_data_dir["align_file"])
    ]
    stdout, stderr, returncode = run_command(command)
    assert returncode != 0, "Alignment should fail with missing reference file."
    assert "does not exist" in stderr, "Error message should mention missing reference file."

    # Test missing reads file
    command = [
        "python3", "main.py",
        "-t", "align",
        "-r", str(test_data_dir["reference_file"]),
        "--reads", str(invalid_file),
        "-a", str(test_data_dir["align_file"])
    ]
    stdout, stderr, returncode = run_command(command)
    assert returncode != 0, "Alignment should fail with missing reads file."
    assert "does not exist" in stderr, "Error message should mention missing reads file."


def test_invalid_task():
    """Tests invalid task handling."""
    command = [
        "python3", "main.py",
        "-t", "invalidtask"
    ]
    stdout, stderr, returncode = run_command(command)
    assert returncode != 0, "Invalid task should return an error."
    assert "Unsupported task" in stderr, "Error message should mention unsupported task."

