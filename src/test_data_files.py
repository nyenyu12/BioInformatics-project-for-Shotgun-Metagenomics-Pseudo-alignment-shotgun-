"""
@file test_data_file.py
@brief Unit tests for data file loading and dumping functionality.
@details Tests for FASTAFile and FASTAQFile classes, including gzip compressed files,
         invalid extensions, empty files, and file dumping.
"""

import gzip
import pickle
import pytest
from pathlib import Path
from data_file import FASTAFile, FASTAQFile, NoRecordsInDataFile, InvalidExtensionError


## ===========================================================
## FASTAFile Tests
## ===========================================================

## @brief Test loading a regular FASTA file.
def test_fasta_file_loading(tmp_path: Path) -> None:
    """
    @brief Test that a FASTA file is loaded correctly.
    @details Writes sample FASTA content to a temporary file and verifies that two records are loaded.
    """
    fasta_content: str = ">Genome1\nAGCTAGCTAGCTA\n>Genome2\nTGCATGCATGCA\n"
    fasta_path: Path = tmp_path / "test.fa"
    fasta_path.write_text(fasta_content)
    
    fasta_file: FASTAFile = FASTAFile(str(fasta_path))
    assert len(list(fasta_file.container)) == 2


## @brief Test loading a gzip-compressed FASTA file.
def test_fasta_gz_file_loading(tmp_path: Path) -> None:
    """
    @brief Test that a gzipped FASTA file is loaded correctly.
    @details Writes sample FASTA content to a gzip file and verifies that two records are loaded.
    """
    fasta_content: str = ">Genome1\nAGCTAGCTAGCTA\n>Genome2\nTGCATGCATGCA\n"
    fasta_gz_path: Path = tmp_path / "test.fa.gz"
    
    with gzip.open(fasta_gz_path, 'wt', encoding='utf-8') as f:
        f.write(fasta_content)
    
    fasta_file: FASTAFile = FASTAFile(str(fasta_gz_path))
    assert len(list(fasta_file.container)) == 2


## @brief Test that a file with an invalid extension for FASTAFile raises an error.
def test_fasta_file_invalid_extension(tmp_path: Path) -> None:
    """
    @brief Test that FASTAFile raises InvalidExtensionError when loaded from a file with an invalid extension.
    """
    invalid_path: Path = tmp_path / "test.txt"
    invalid_path.write_text(">Genome1\nAGCTAGCTAGCTA\n")
    
    with pytest.raises(InvalidExtensionError):
        FASTAFile(str(invalid_path))


## @brief Test that an empty FASTA file raises NoRecordsInDataFile.
def test_no_records_in_fasta_file(tmp_path: Path) -> None:
    """
    @brief Test that parsing an empty FASTA file raises NoRecordsInDataFile.
    """
    empty_fasta_path: Path = tmp_path / "empty.fa"
    empty_fasta_path.write_text("")
    
    with pytest.raises(NoRecordsInDataFile):
        FASTAFile(str(empty_fasta_path))


## @brief Test dumping a FASTAFile container.
def test_fasta_file_dump(tmp_path: Path) -> None:
    """
    @brief Test that a FASTAFile container can be dumped and reloaded.
    @details Dumps the parsed FASTA container to a pickle file and reloads it, verifying two records are present.
    """
    fasta_content: str = ">Genome1\nAGCTAGCTAGCTA\n>Genome2\nTGCATGCATGCA\n"
    fasta_path: Path = tmp_path / "test.fa"
    dump_path: Path = tmp_path / "dump.pkl"
    fasta_path.write_text(fasta_content)
    
    fasta_file: FASTAFile = FASTAFile(str(fasta_path))
    fasta_file.dump(str(dump_path))
    
    with open(dump_path, 'rb') as f:
        data = pickle.load(f)
    
    assert len(list(data)) == 2


## ===========================================================
## FASTAQFile Tests
## ===========================================================

## @brief Test loading a regular FASTQ file.
def test_fastq_file_loading(tmp_path: Path) -> None:
    """
    @brief Test that a FASTQ file is loaded correctly.
    @details Writes sample FASTQ content to a temporary file and verifies that two records are loaded.
    """
    fastq_content: str = "@Read1\nAGCTAGCT\n+\nIIIIIIII\n@Read2\nTGCATGCA\n+\nIIIIIIII\n"
    fastq_path: Path = tmp_path / "test.fq"
    fastq_path.write_text(fastq_content)
    
    fastq_file: FASTAQFile = FASTAQFile(str(fastq_path))
    assert len(list(fastq_file.container)) == 2


## @brief Test loading a gzip-compressed FASTQ file.
def test_fastq_gz_file_loading(tmp_path: Path) -> None:
    """
    @brief Test that a gzipped FASTQ file is loaded correctly.
    @details Writes sample FASTQ content to a gzip file and verifies that two records are loaded.
    """
    fastq_content: str = "@Read1\nAGCTAGCT\n+\nIIIIIIII\n@Read2\nTGCATGCA\n+\nIIIIIIII\n"
    fastq_gz_path: Path = tmp_path / "test.fq.gz"
    
    with gzip.open(fastq_gz_path, 'wt', encoding='utf-8') as f:
        f.write(fastq_content)
    
    fastq_file: FASTAQFile = FASTAQFile(str(fastq_gz_path))
    assert len(list(fastq_file.container)) == 2


## @brief Test that a file with an invalid extension for FASTAQFile raises an error.
def test_fastq_file_invalid_extension(tmp_path: Path) -> None:
    """
    @brief Test that FASTAQFile raises InvalidExtensionError when loaded from a file with an invalid extension.
    """
    invalid_path: Path = tmp_path / "test.txt"
    invalid_path.write_text("@Read1\nAGCTAGCT\n+\nIIIIIIII\n")
    
    with pytest.raises(InvalidExtensionError):
        FASTAQFile(str(invalid_path))


## @brief Test that an empty FASTQ file raises NoRecordsInDataFile.
def test_no_records_in_fastq_file(tmp_path: Path) -> None:
    """
    @brief Test that parsing an empty FASTQ file raises NoRecordsInDataFile.
    """
    empty_fastq_path: Path = tmp_path / "empty.fq"
    empty_fastq_path.write_text("")
    
    with pytest.raises(NoRecordsInDataFile):
        FASTAQFile(str(empty_fastq_path))
