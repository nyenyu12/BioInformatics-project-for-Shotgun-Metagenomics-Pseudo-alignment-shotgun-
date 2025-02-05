import pytest
import gzip
import pickle
from data_file import FASTAFile, FASTAQFile, NoRecordsInDataFile, InvalidExtensionError

def test_fasta_file_loading(tmp_path):
    fasta_content = ">Genome1\nAGCTAGCTAGCTA\n>Genome2\nTGCATGCATGCA\n"
    fasta_path = tmp_path / "test.fa"
    fasta_path.write_text(fasta_content)
    
    fasta_file = FASTAFile(str(fasta_path))
    assert len(list(fasta_file.container)) == 2

def test_fasta_gz_file_loading(tmp_path):
    fasta_content = ">Genome1\nAGCTAGCTAGCTA\n>Genome2\nTGCATGCATGCA\n"
    fasta_gz_path = tmp_path / "test.fa.gz"
    
    with gzip.open(fasta_gz_path, 'wt', encoding='utf-8') as f:
        f.write(fasta_content)
    
    fasta_file = FASTAFile(str(fasta_gz_path))
    assert len(list(fasta_file.container)) == 2

def test_fasta_file_invalid_extension(tmp_path):
    invalid_path = tmp_path / "test.txt"
    invalid_path.write_text(">Genome1\nAGCTAGCTAGCTA\n")
    
    with pytest.raises(InvalidExtensionError):
        FASTAFile(str(invalid_path))

def test_fastq_file_loading(tmp_path):
    fastq_content = "@Read1\nAGCTAGCT\n+\nIIIIIIII\n@Read2\nTGCATGCA\n+\nIIIIIIII\n"
    fastq_path = tmp_path / "test.fq"
    fastq_path.write_text(fastq_content)
    
    fastq_file = FASTAQFile(str(fastq_path))
    assert len(list(fastq_file.container)) == 2

def test_fastq_gz_file_loading(tmp_path):
    fastq_content = "@Read1\nAGCTAGCT\n+\nIIIIIIII\n@Read2\nTGCATGCA\n+\nIIIIIIII\n"
    fastq_gz_path = tmp_path / "test.fq.gz"
    
    with gzip.open(fastq_gz_path, 'wt', encoding='utf-8') as f:
        f.write(fastq_content)
    
    fastq_file = FASTAQFile(str(fastq_gz_path))
    assert len(list(fastq_file.container)) == 2

def test_fastq_file_invalid_extension(tmp_path):
    invalid_path = tmp_path / "test.txt"
    invalid_path.write_text("@Read1\nAGCTAGCT\n+\nIIIIIIII\n")
    
    with pytest.raises(InvalidExtensionError):
        FASTAQFile(str(invalid_path))

def test_no_records_in_fasta_file(tmp_path):
    empty_fasta_path = tmp_path / "empty.fa"
    empty_fasta_path.write_text("")
    
    with pytest.raises(NoRecordsInDataFile):
        FASTAFile(str(empty_fasta_path))

def test_no_records_in_fastq_file(tmp_path):
    empty_fastq_path = tmp_path / "empty.fq"
    empty_fastq_path.write_text("")
    
    with pytest.raises(NoRecordsInDataFile):
        FASTAQFile(str(empty_fastq_path))

def test_fasta_file_dump(tmp_path):
    fasta_content = ">Genome1\nAGCTAGCTAGCTA\n>Genome2\nTGCATGCATGCA\n"
    fasta_path = tmp_path / "test.fa"
    dump_path = tmp_path / "dump.pkl"
    fasta_path.write_text(fasta_content)
    
    fasta_file = FASTAFile(str(fasta_path))
    fasta_file.dump(str(dump_path))
    
    with open(dump_path, 'rb') as f:
        data = pickle.load(f)
    
    assert len(list(data)) == 2
