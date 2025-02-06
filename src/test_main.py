"""
@file test_main.py
@brief Unit tests for main.py functionalities.
@details This test suite verifies reference creation, dumping, alignment, file validation,
         argument parsing, and EXTSIM filtering via the command-line interface.
"""

import subprocess
import pytest
import os
import json
import gzip
import pickle
from pathlib import Path
from typing import List, Tuple, Dict, Any
from unittest.mock import patch

from kmer import KmerReference
from records import FASTARecordContainer
from main import (
    validate_file_readable,
    validate_file_writable,
    parse_arguments,
)

## ===========================================================
## Helper Functions
## ===========================================================

def run_command(command: List[str]) -> Tuple[str, str, int]:
    """
    @brief Runs a subprocess command and returns its output.
    @param command List of command arguments.
    @return A tuple containing (stdout, stderr, return code).
    """
    result = subprocess.run(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    return result.stdout, result.stderr, result.returncode


## ===========================================================
## Fixtures
## ===========================================================

@pytest.fixture
def test_data_dir(tmp_path: Path) -> Dict[str, Path]:
    """
    @brief Creates a temporary directory with test data files.
    @param tmp_path A temporary directory provided by pytest.
    @return Dictionary with paths for genome, FASTQ, reference, and alignment files.
    """
    genome_file: Path = tmp_path / "test_genome.fa"
    fastq_file: Path = tmp_path / "test_reads.fq"
    reference_file: Path = tmp_path / "test_reference.kdb"
    align_file: Path = tmp_path / "test_alignment.aln"

    # Create a sample genome file.
    genome_data: str = ">Genome1\nAGCTAGCTAGCTAGCTAGCT\n>Genome2\nTGCATGCATGCATGCATGCA\n"
    genome_file.write_text(genome_data)

    # Create a sample reads file.
    reads_data: str = "@Read1\nAGCTAGCT\n+\nIIIIIIII\n@Read2\nTGCATGCA\n+\nIIIIIIII\n"
    fastq_file.write_text(reads_data)

    return {
        "genome_file": genome_file,
        "fastq_file": fastq_file,
        "reference_file": reference_file,
        "align_file": align_file,
    }


## ===========================================================
## Tests: Reference and Alignment Commands
## ===========================================================

def test_create_reference(test_data_dir: Dict[str, Path]) -> None:
    """
    @brief Test creating a k-mer reference database.
    """
    stdout, stderr, returncode = run_command([
        "python3", "main.py",
        "-t", "reference",
        "-g", str(test_data_dir["genome_file"]),
        "-k", "4",
        "-r", str(test_data_dir["reference_file"])
    ])
    assert returncode == 0, f"Error creating reference: {stderr}"
    assert test_data_dir["reference_file"].exists(), "Reference file was not created."

    # Check that the reference file is a valid pickle.
    with gzip.open(test_data_dir["reference_file"], "rb") as f:
        kmer_reference = pickle.load(f)
    assert kmer_reference is not None, "Failed to load reference from file."


def test_dump_reference(test_data_dir: Dict[str, Path]) -> None:
    """
    @brief Test dumping a reference database to JSON.
    """
    # Ensure reference is created first.
    test_create_reference(test_data_dir)

    stdout, stderr, returncode = run_command([
        "python3", "main.py",
        "-t", "dumpref",
        "-r", str(test_data_dir["reference_file"])
    ])
    assert returncode == 0, f"Error dumping reference: {stderr}"
    try:
        json_output: Dict[str, Any] = json.loads(stdout)
        assert "Summary" in json_output, "Reference dump output does not contain expected summary."
    except json.JSONDecodeError:
        pytest.fail("Reference dump output is not valid JSON.")


def test_align_reads(test_data_dir: Dict[str, Path]) -> None:
    """
    @brief Test aligning reads using a pre-built reference.
    """
    # Ensure reference is created first.
    test_create_reference(test_data_dir)
    stdout, stderr, returncode = run_command([
        "python3", "main.py",
        "-t", "align",
        "-r", str(test_data_dir["reference_file"]),
        "--reads", str(test_data_dir["fastq_file"]),
        "-a", str(test_data_dir["align_file"])
    ])
    assert returncode == 0, f"Error aligning reads: {stderr}"
    assert test_data_dir["align_file"].exists(), "Alignment file was not created."
    with gzip.open(test_data_dir["align_file"], "rb") as f:
        pseudo_alignment = pickle.load(f)
    assert pseudo_alignment is not None, "Failed to load alignment from file."


def test_dump_alignment(test_data_dir: Dict[str, Path]) -> None:
    """
    @brief Test dumping an alignment file.
    """
    # Ensure alignment is created first.
    test_align_reads(test_data_dir)
    stdout, stderr, returncode = run_command([
        "python3", "main.py",
        "-t", "dumpalign",
        "-a", str(test_data_dir["align_file"])
    ])
    assert returncode == 0, f"Error dumping alignment: {stderr}"
    try:
        json_output: Dict[str, Any] = json.loads(stdout)
        assert "Statistics" in json_output, "Alignment dump output does not contain expected statistics."
    except json.JSONDecodeError:
        pytest.fail("Alignment dump output is not valid JSON.")


def test_build_reference_align_and_dump(test_data_dir: Dict[str, Path]) -> None:
    """
    @brief Test building a reference, aligning reads, and dumping results in one step.
    """
    stdout, stderr, returncode = run_command([
        "python3", "main.py",
        "-t", "dumpalign",
        "-g", str(test_data_dir["genome_file"]),
        "-k", "4",
        "--reads", str(test_data_dir["fastq_file"])
    ])
    assert returncode == 0, f"Error in build-reference-align-dump: {stderr}"
    try:
        json_output: Dict[str, Any] = json.loads(stdout)
        assert "Statistics" in json_output, "Build-reference-align-dump output does not contain expected statistics."
    except json.JSONDecodeError:
        pytest.fail("Build-reference-align-dump output is not valid JSON.")


## ===========================================================
## Tests: Invalid File Handling
## ===========================================================

def test_invalid_file_handling(test_data_dir: Dict[str, Path]) -> None:
    """
    @brief Test handling of missing or invalid files.
    """
    invalid_file: Path = test_data_dir["genome_file"].parent / "invalid_file.fa"
    stdout, stderr, returncode = run_command([
        "python3", "main.py",
        "-t", "reference",
        "-g", str(invalid_file),
        "-k", "4",
        "-r", str(test_data_dir["reference_file"])
    ])
    assert returncode != 0, "Reference creation should fail with missing genome file."
    assert "does not exist" in stderr, "Error message should mention missing file."

    invalid_ext_file: Path = test_data_dir["genome_file"].parent / "invalid_file.txt"
    invalid_ext_file.touch()
    stdout, stderr, returncode = run_command([
        "python3", "main.py",
        "-t", "reference",
        "-g", str(invalid_ext_file),
        "-k", "4",
        "-r", str(test_data_dir["reference_file"])
    ])
    assert returncode != 0, "Reference creation should fail with invalid extension."
    assert "Invalid file extension" in stderr, "Error message should indicate invalid extension."

    stdout, stderr, returncode = run_command([
        "python3", "main.py",
        "-t", "align",
        "-r", str(invalid_file),
        "--reads", str(test_data_dir["fastq_file"]),
        "-a", str(test_data_dir["align_file"])
    ])
    assert returncode != 0, "Alignment should fail with missing reference file."
    assert "does not exist" in stderr, "Error message should mention missing reference file."

    stdout, stderr, returncode = run_command([
        "python3", "main.py",
        "-t", "align",
        "-r", str(test_data_dir["reference_file"]),
        "--reads", str(invalid_file),
        "-a", str(test_data_dir["align_file"])
    ])
    assert returncode != 0, "Alignment should fail with missing reads file."
    assert "does not exist" in stderr, "Error message should mention missing reads file."


## ===========================================================
## Tests: Invalid Task Handling
## ===========================================================

def test_invalid_task() -> None:
    """
    @brief Test that an invalid task returns an error.
    """
    stdout, stderr, returncode = run_command([
        "python3", "main.py",
        "-t", "invalidtask"
    ])
    assert returncode != 0, "Invalid task should return an error."
    assert "Unsupported task" in stderr, "Error message should mention unsupported task."


## ===========================================================
## Tests: EXTSIM (Extension 4.4)
## ===========================================================

def test_extsim_filter_similar_genomes(tmp_path: Path) -> None:
    """
    @brief Test EXTSIM filtering of highly similar genomes.
    """
    fasta_content: str = (
        ">GenomeA\nAGCTAGCTAGCT\n"
        ">GenomeB\nAGCTAGCTAGCT\n"  # Identical to GenomeA; should be filtered.
        ">GenomeC\nTGCATGCATGCA\n"   # Distinct; should be kept.
    )
    fasta_path: Path = tmp_path / "test_extsim.fa"
    fasta_path.write_text(fasta_content)
    
    stdout, stderr, returncode = run_command([
        "python3", "main.py",
        "-t", "reference",
        "-g", str(fasta_path),
        "-k", "4",
        "-r", str(tmp_path / "extsim_reference.kdb"),
        "--filter-similar",
        "--similarity-threshold", "0.95"
    ])
    assert returncode == 0, f"Error creating reference with EXTSIM: {stderr}"
    
    reference_file: Path = tmp_path / "extsim_reference.kdb"
    with gzip.open(reference_file, "rb") as f:
        kmer_reference: Any = pickle.load(f)
    assert hasattr(kmer_reference, "similarity_info"), "EXTSIM filtering info missing."
    sim_info: Dict[str, Any] = kmer_reference.similarity_info
    kept_ids = {genome.identifier for genome in kmer_reference.genomes}
    assert "GenomeC" in kept_ids, "GenomeC should be kept."
    if "GenomeA" in kept_ids:
        assert "GenomeB" not in kept_ids, "GenomeB should be filtered out."
    else:
        assert "GenomeB" in kept_ids, "GenomeA should be filtered out."

def test_no_filter_similar_genomes():
    """
    @brief Test that when filtering is disabled, all genomes are retained.
    """
    data = (
        ">GenomeA\nAGCTAGCTAGCT\n"
        ">GenomeB\nAGCTAGCTAGCT\n"
    )
    container = FASTARecordContainer()
    container.parse_records(data)
    kref = KmerReference(4, container, filter_similar=False)
    assert not hasattr(kref, "similarity_info")
    kept_ids = {genome.identifier for genome in kref.genomes}
    assert "GenomeA" in kept_ids and "GenomeB" in kept_ids

def test_similarity_info_in_summary():
    """
    @brief Test that the JSON summary includes a 'Similarity' section when filtering is enabled.
    """
    data = (
        ">GenomeA\nAGCTAGCTAGCT\n"
        ">GenomeB\nAGCTAGCTAGCT\n"
        ">GenomeC\nTGCATGCATGCA\n"
    )
    container = FASTARecordContainer()
    container.parse_records(data)
    kref = KmerReference(4, container, filter_similar=True, similarity_threshold=0.95)
    summary = kref.get_summary()
    assert "Similarity" in summary
    sim_info = summary["Similarity"]
    assert "GenomeA" in sim_info
    assert "GenomeB" in sim_info
    assert "GenomeC" in sim_info

def test_single_genome_no_filter():
    """
    @brief Test that with a single genome, filtering does nothing.
    """
    data = ">GenomeA\nAGCTAGCTAGCT\n"
    container = FASTARecordContainer()
    container.parse_records(data)
    kref = KmerReference(4, container, filter_similar=True, similarity_threshold=0.95)
    kept_ids = {genome.identifier for genome in kref.genomes}
    assert kept_ids == {"GenomeA"}
    summary = kref.get_summary()
    assert "Similarity" in summary
    assert summary["Similarity"]["GenomeA"]["kept"] == "yes"
    assert summary["Similarity"]["GenomeA"]["similar_to"] == "NA"

## ===========================================================
## Tests: Additional Main.py Functions
## ===========================================================

def test_validate_file_readable(tmp_path: Path) -> None:
    """
    @brief Test that validate_file_readable correctly detects a missing file.
    """
    test_file: Path = tmp_path / "missing_file.fq"
    with pytest.raises(SystemExit):
        validate_file_readable(str(test_file), "FASTQ reads")

def test_validate_file_writable(tmp_path: Path) -> None:
    """
    @brief Test that validate_file_writable correctly detects unwritable locations.
    """
    unwritable_dir: Path = tmp_path / "readonly_dir"
    unwritable_dir.mkdir()
    os.chmod(unwritable_dir, 0o400)  # Make directory read-only.
    unwritable_file: Path = unwritable_dir / "output.aln"
    with pytest.raises(SystemExit):
        validate_file_writable(str(unwritable_file), "Alignment output")
    os.chmod(unwritable_dir, 0o700)  # Reset permissions.

@patch("sys.exit")
def test_parse_arguments_default_values(mock_exit: Any) -> None:
    """
    @brief Test that parse_arguments sets default values when optional arguments are missing.
    """
    test_args: List[str] = ["-t", "align", "--reads", "reads.fq", "-a", "output.aln"]
    with patch("sys.argv", ["main.py"] + test_args):
        args = parse_arguments()
    assert args.min_read_quality is None
    assert args.min_kmer_quality is None
    assert args.max_genomes is None
