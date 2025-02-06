import pytest
import random
import pickle
import gzip
from records import FASTARecordContainer, FASTAQRecordContainer
from kmer import (
    Read,
    ReadMappingType,
    KmerReference,
    extract_kmers_from_genome,
    ReadMapping,
    AddingExistingRead,
    PseudoAlignment
)


@pytest.fixture
def kmer_reference(sample_fasta_container):
    return KmerReference(3, sample_fasta_container)


@pytest.fixture
def sample_reads():
    data = (
        "@Read1\nAGCTAGCT\n+\nIIIIIIII\n"
        "@Read2\nTGCATGCA\n+\nIIIIIIII\n"
        "@Read3\nGGGGGGGG\n+\nIIIIIIII\n"  # No matching k-mers
        "@Read4\nAGCTTCGC\n+\nIIIIIIII\n"  # Partial match in multiple genomes
        "@Read5\nTGCATGCTAGCTA\n+\nIITT55IIIIII9\n"  # Long read covering k-mers in multiple genomes
        "@Read6\nCCGGAAGCTTGCATGCA\n+\nIIIIIIG$IIIIIIIII\n"  # Long read that matches Genome4
        "@Read7\nAGCGTAGCTAGCTAGCT\n+\nIIIIIII&IIIIIIIII\n"  # Fully matching Genome1
    )
    container = FASTAQRecordContainer()
    container.parse_records(data)
    reads = [Read(read) for read in container]
    for r in reads:
        print(r)
    return reads


@pytest.fixture
def sample_fasta_container():
    data = (
        ">Genome1\nAGCTAGCTAGCTAGCTAGCT\n"
        ">Genome2\nTGCATGCATGCATGCATGCA\n"
        ">Genome3\nAGCTTGCATGCAGCTAGCTA\n"  # Contains k-mers appearing in multiple genomes
        ">Genome4\nCCGGAAGCTTGCATGCAGCTA\n"  # Another genome with overlapping k-mers
    )
    container = FASTARecordContainer()
    container.parse_records(data)
    return container


def test_extract_kmers_from_genome():
    genome = "AGCTAGCTAGCT"
    k = 3
    expected_kmers = [
        (0, "AGC"),
        (1, "GCT"),
        (2, "CTA"),
        (3, "TAG"),
        (4, "AGC"),
        (5, "GCT"),
        (6, "CTA"),
        (7, "TAG"),
        (8, "AGC"),
        (9, "GCT"),
    ]
    assert list(extract_kmers_from_genome(k, genome)) == expected_kmers


def test_kmer_reference_building(kmer_reference):
    assert kmer_reference.get_kmer_references("AGC")
    assert kmer_reference.get_kmer_references("TGC")
    assert not kmer_reference.get_kmer_references(
        "GGG"
    )  # k-mer not present in reference genomes
    assert kmer_reference.get_kmer_references(
        "GCT"
    )  # k-mer appearing in multiple genomes
    assert kmer_reference.get_kmer_references("CCG")  # k-mer exclusive to Genome4


# Case 1: Unmapped Read
@pytest.fixture
def fasta_container_unmapped():
    data = (
        ">Genome1\nAACCGGTTAACC\n"  # No matching k-mers
        ">Genome2\nGGTTCCAAGGTT\n"  # No matching k-mers
    )
    container = FASTARecordContainer()
    container.parse_records(data)
    return container


@pytest.fixture
def kmer_reference_unmapped(fasta_container_unmapped):
    return KmerReference(4, fasta_container_unmapped)


@pytest.fixture
def read_unmapped():
    data = "@Read1\nTAGGCAT\n+\nIIIIIII\n"
    container = FASTAQRecordContainer()
    container.parse_records(data)
    return Read(list(container)[0])


def test_unmapped_read(kmer_reference_unmapped, read_unmapped):
    assert (
        read_unmapped.pseudo_align(kmer_reference_unmapped) == ReadMappingType.UNMAPPED
    )


# Case 2: Uniquely Mapped Read
@pytest.fixture
def fasta_container_unique():
    data = ">Genome1\nATGGCTATGCTA\n" ">Genome2\nCTATGGCAGGCA\n"
    container = FASTARecordContainer()
    container.parse_records(data)
    return container


@pytest.fixture
def kmer_reference_unique(fasta_container_unique):
    return KmerReference(4, fasta_container_unique)


@pytest.fixture
def read_unique():
    data = "@Read2\nATGGCTAT\n+\nIIIIIIII\n"
    container = FASTAQRecordContainer()
    container.parse_records(data)
    return Read(list(container)[0])


def test_uniquely_mapped_read(kmer_reference_unique, read_unique):
    assert (
        read_unique.pseudo_align(kmer_reference_unique)
        == ReadMappingType.UNIQUELY_MAPPED
    )


# Case 3: Ambiguous Mapping (With Multiple Genomes)
@pytest.fixture
def fasta_container_ambiguous():
    data = (
        ">Genome1\nATCGACGGTCGTTA\n"
        ">Genome2\nCGATGATCAGTACGA\n"
        ">Genome3\nATCCACCTAACGTACGGT\n"
        ">Genome4\nCTAGGGACTGCACTA\n"
    )
    container = FASTARecordContainer()
    container.parse_records(data)
    return container


@pytest.fixture
def kmer_reference_ambiguous(fasta_container_ambiguous):
    return KmerReference(4, fasta_container_ambiguous)


@pytest.fixture
def read_ambiguous():
    data = "@Read3\nATCGATCCTAG\n+\nIIIIIIIIIII\n"
    container = FASTAQRecordContainer()
    container.parse_records(data)
    return Read(list(container)[0])


def test_ambiguously_mapped_read(kmer_reference_ambiguous, read_ambiguous):
    result = read_ambiguous.pseudo_align(kmer_reference_ambiguous)
    print("Ambiguous Read Alignment Result:", result)
    assert result == ReadMappingType.AMBIGUOUSLY_MAPPED
    assert len(read_ambiguous.mapping.genomes_mapped_to) == 4


# Case 4: Initially Unique but Changed to Ambiguous
@pytest.fixture
def fasta_container_initially_unique():
    data = (
        ">Genome1\nATGCCTTTTCGGGG\n"
        ">Genome2\nGCCGTTTTCGGGGCTA\n"
        ">Genome3\nCCGG\n"
        ">Genome4\nAAAAAAAAGGGCT\n"
        ">Genome5\nTTTTTTTTGCTAA\n"
    )
    container = FASTARecordContainer()
    container.parse_records(data)
    return container


@pytest.fixture
def kmer_reference_initially_unique(fasta_container_initially_unique):
    return KmerReference(4, fasta_container_initially_unique)


@pytest.fixture
def read_initially_unique():
    data = "@Read4\nATGCCGGGGCTAA\n+\nIIIIIIIIIIIII\n"
    container = FASTAQRecordContainer()
    container.parse_records(data)
    return Read(list(container)[0])


def test_initially_unique_but_changes_to_ambiguous_tests_default_p(
    kmer_reference_initially_unique, read_initially_unique
):
    result = read_initially_unique.pseudo_align(kmer_reference_initially_unique)
    print("Initially Unique But Changed to Ambiguous Read Alignment Result:", result)
    assert result == ReadMappingType.AMBIGUOUSLY_MAPPED


def test_initially_unique_but_changes_to_ambiguous_tests_big_p(
    kmer_reference_initially_unique, read_initially_unique
):
    result = read_initially_unique.pseudo_align(kmer_reference_initially_unique, p=5)
    print("Initially Unique But Changed to Ambiguous Read Alignment Result:", result)
    assert result == ReadMappingType.UNIQUELY_MAPPED


@pytest.fixture
def big_read():
    data = "@BigRead\nAGCTAGCTAGAGGTCCTAATCCTAGCTAGCTAGCTAGCTAGCTAGCTGGTCATCAAAACCTTT\n+\nIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII\n"
    container = FASTAQRecordContainer()
    container.parse_records(data)
    return Read(list(container)[0])


@pytest.fixture
def synthetic_kmer_reference(big_read):
    k = 31
    extracted_kmers = list(extract_kmers_from_genome(k, big_read._Read__raw_read))

    # Randomly distribute k-mers across genomes, allowing duplicates
    num_genomes = 4
    genomes = {f">Genome{i+1}": [] for i in range(num_genomes)}

    for _, kmer in extracted_kmers:
        selected_genomes = random.sample(
            list(genomes.keys()), k=random.randint(1, num_genomes)
        )
        for genome in selected_genomes:
            genomes[genome].append(kmer)

    fasta_data = "".join(f"{name}\n{"NN".join(seq)}\n" for name, seq in genomes.items())

    for name, sequence in genomes.items():
        print(f"{name}: {sequence}")

    container = FASTARecordContainer()
    container.parse_records(fasta_data)
    return KmerReference(k, container)


def test_kmer_reference_counts(synthetic_kmer_reference):
    kmer_counts = {}
    unspecific_count = 0

    for kmer in synthetic_kmer_reference.kmers.values():
        for genome, positions in kmer.items():
            kmer_counts[genome] = kmer_counts.get(genome, 0) + len(positions)

    print("Specific k-mer counts per genome:", kmer_counts)
    print("Unspecific k-mer count:", unspecific_count)

    assert sum(kmer_counts.values()) > 0  # Ensure k-mers are assigned


@pytest.mark.parametrize("execution_number", range(10))
def test_pseudo_align_random_big_kmers(
    capsys, big_read, synthetic_kmer_reference, execution_number
):
    """
    This is meant to almost implement the algorithm a second time to verify the main implementation.
    In addition, because the spererator for the k-mers is 'NN', then this tests that as well :)
    """
    m = 1
    p = 1
    specific_counts = {}
    total_counts = {}

    for k in synthetic_kmer_reference.kmers.keys():
        print(k)
    # Notice all k-mers are in the Read since the genomes are generated from k-mers in the read
    for kmer, genome_map in synthetic_kmer_reference.kmers.items():
        mapped_genomes = list(genome_map.keys())
        if len(mapped_genomes) == 1:
            specific_counts[mapped_genomes[0]] = specific_counts.get(
                mapped_genomes[0], 0
            ) + len(genome_map[mapped_genomes[0]])
        for genome, positions in genome_map.items():
            total_counts[genome] = total_counts.get(genome, 0) + len(positions)

    print("Specific k-mer counts per genome:", specific_counts)
    print("Total k-mer counts per genome:", total_counts)
    sorted_genomes = sorted(specific_counts, key=specific_counts.get, reverse=True)
    max_genome = sorted_genomes[0]
    second_max_genome = sorted_genomes[1] if len(sorted_genomes) > 1 else None

    specific_max_count = specific_counts[max_genome]
    specific_second_max_count = (
        specific_counts[second_max_genome] if second_max_genome else 0
    )
    m_difference = specific_max_count - specific_second_max_count

    alignment_result = big_read.pseudo_align(
        synthetic_kmer_reference, p=p, m=1, debug=True
    )
    captured = capsys.readouterr()

    max_total_kmers = max(total_counts.values())
    map_total_kmers = total_counts[max_genome]
    p_difference = max_total_kmers - map_total_kmers

    print(captured.out)
    print("Specific MaxCount:", specific_max_count)
    print("Second Specific Max Count:", specific_second_max_count)
    print("Difference:", m_difference)
    print("MaxCount", max_total_kmers)
    print("MapCount", map_total_kmers)
    print("Alignment Result:", alignment_result)

    # Because of the test is built there should be no case of ReadMappingType.UNMAPPED
    if (m_difference >= m) and (not p_difference > p):
        assert alignment_result == ReadMappingType.UNIQUELY_MAPPED
    elif (m_difference >= m) and (p_difference > p):
        assert alignment_result == ReadMappingType.AMBIGUOUSLY_MAPPED
    else:
        assert alignment_result == ReadMappingType.AMBIGUOUSLY_MAPPED


# Testing PseudoAlignment
# Mock KmerReference to avoid unnecessary dependencies
class MockKmerReference:
    def get_kmer_references(self, kmer):
        return {}


@pytest.fixture
def mock_kmer_reference():
    return MockKmerReference()


@pytest.fixture
def pseudo_alignment(mock_kmer_reference):
    return PseudoAlignment(mock_kmer_reference)


@pytest.fixture
def create_read():
    """Creates a Read object from FASTAQRecordContainer"""
    def _create_read(identifier, sequence, quality, genomes_mapped_to=None, mapping_type="UNMAPPED"):
        fastaq_data = f"@{identifier}\n{sequence}\n+\n{quality}\n"
        container = FASTAQRecordContainer()
        container.parse_records(fastaq_data)
        read = Read(list(container)[0])

        if genomes_mapped_to or genomes_mapped_to == []:
            read.mapping = ReadMapping(ReadMappingType[mapping_type], genomes_mapped_to)

        return read

    return _create_read


@pytest.fixture
def create_record():
    """Creates a Record object using FASTARecordContainer"""
    def _create_record(genome_id, genome_sequence="AGCTAGCTAG"):
        fasta_data = f">{genome_id}\n{genome_sequence}\n"
        container = FASTARecordContainer()
        container.parse_records(fasta_data)
        return list(container)[0]  # Return the first (and only) record

    return _create_record


def test_add_read_for_pseudo_alignment(pseudo_alignment, create_read, create_record):
    genome1 = create_record("Genome1")
    read1 = create_read("read1", "AGCTAGCT", "IIIIIIII", [genome1], "UNIQUELY_MAPPED")

    pseudo_alignment.add_read(read1)

    assert "read1" in pseudo_alignment.reads
    assert pseudo_alignment.reads["read1"]["mapping_type"] == ReadMappingType.UNIQUELY_MAPPED
    assert pseudo_alignment.reads["read1"]["genomes_mapped_to"] == [genome1.identifier]


def test_add_duplicate_read_raises_exception_for_pseudo_alignment(pseudo_alignment, create_read, create_record):
    genome1 = create_record("Genome1")
    read1 = create_read("read1", "AGCTAGCT", "IIIIIIII", [genome1], "UNIQUELY_MAPPED")

    pseudo_alignment.add_read(read1)

    with pytest.raises(AddingExistingRead):
        pseudo_alignment.add_read(read1)


def test_get_summary_for_pseudo_alignment(pseudo_alignment, create_read, create_record):
    genome1 = create_record("Genome1")
    genome2 = create_record("Genome2")

    read1 = create_read("read1", "AGCTAGCT", "IIIIIIII", [genome1], "UNIQUELY_MAPPED")
    read2 = create_read("read2", "TGCATGCA", "IIIIIIII", [genome1, genome2], "AMBIGUOUSLY_MAPPED")
    read3 = create_read("read3", "GGGGGGGG", "IIIIIIII", [], "UNMAPPED")

    pseudo_alignment.add_read(read1)
    pseudo_alignment.add_read(read2)
    pseudo_alignment.add_read(read3)

    summary = pseudo_alignment.get_summary()

    assert summary["Statistics"]["unique_mapped_reads"] == 1
    assert summary["Statistics"]["ambiguous_mapped_reads"] == 1
    assert summary["Statistics"]["unmapped_reads"] == 1

    assert summary["Summary"]["Genome1"]["unique_reads"] == 1
    assert summary["Summary"]["Genome1"]["ambiguous_reads"] == 1
    assert summary["Summary"]["Genome2"]["unique_reads"] == 0
    assert summary["Summary"]["Genome2"]["ambiguous_reads"] == 1


def test_save_and_load_for_pseudo_alignment(pseudo_alignment, create_read, create_record, tmp_path):
    genome1 = create_record("Genome1")
    read1 = create_read("read1", "AGCTAGCT", "IIIIIIII", [genome1], "UNIQUELY_MAPPED")

    pseudo_alignment.add_read(read1)

    align_file = tmp_path / "test.aln"
    pseudo_alignment.save(str(align_file))

    with gzip.open(str(align_file), "rb") as f:
        loaded_alignment = pickle.load(f)

    assert "read1" in loaded_alignment.reads
    assert loaded_alignment.reads["read1"]["mapping_type"] == ReadMappingType.UNIQUELY_MAPPED
    assert loaded_alignment.reads["read1"]["genomes_mapped_to"] == [genome1.identifier]


@pytest.fixture
def sample_fasta_container():
    """Create a mock FASTA file with multiple genomes for testing."""
    data = (
        ">Genome1\nAGCTAGCTAGCTAGCTAGCT\n"
        ">Genome2\nTGCATGCATGCATGCATGCA\n"
        ">Genome3\nAGCTTGCATGCAGCTAGCTA\n"
        ">Genome4\nCCGGAAGCTTGCATGCAGCTA\n"
    )
    container = FASTARecordContainer()
    container.parse_records(data)
    return container

@pytest.fixture
def kmer_reference(sample_fasta_container):
    return KmerReference(3, sample_fasta_container)

@pytest.fixture
def sample_fastq_container():
    """Create a mock FASTQ file with quality scores."""
    data = (
        "@Read1\nAGCTAGCT\n+\nIIIIIIII\n"  # High-quality read (ASCII 'I' -> high score)
        "@Read2\nTGCATGCA\n+\n!!!!!!!!\n"  # Low-quality read (ASCII '!' -> very low)
        "@Read3\nGGGGGGGG\n+\n!!IIIIII\n"  # Another high-quality read
    )
    container = FASTAQRecordContainer()
    container.parse_records(data)
    return container

@pytest.fixture
def sample_reads(sample_fastq_container):
    return [Read(record) for record in sample_fastq_container]

@pytest.fixture
def pseudo_alignment(kmer_reference):
    return PseudoAlignment(kmer_reference)


### **Test Read-Level Quality Filtering**
def test_min_read_quality_filter(pseudo_alignment, sample_fastq_container):
    """Ensure that reads below `min_read_quality` are ignored."""
    pseudo_alignment.align_reads_from_container(sample_fastq_container, min_read_quality=40)
    
    summary = pseudo_alignment.get_summary()
    print (summary)
    assert summary["Statistics"]["filtered_quality_reads"] == 1  # One read is filtered
    assert summary["Statistics"]["unmapped_reads"] >= 1  # At least one unmapped read


### **Test K-mer-Level Quality Filtering**
def test_min_kmer_quality_filter(pseudo_alignment, sample_fastq_container):
    """Ensure that k-mers below `min_kmer_quality` are ignored."""
    pseudo_alignment.align_reads_from_container(sample_fastq_container, min_kmer_quality=60)
    
    summary = pseudo_alignment.get_summary()
    print (summary)
    assert summary["Statistics"]["filtered_quality_kmers"] > 0  # Some k-mers should be filtered


### **Test Highly Redundant K-mer Filtering**
def test_max_genomes_filter(pseudo_alignment, sample_fastq_container):
    """Ensure that k-mers appearing in too many genomes are ignored."""
    pseudo_alignment.align_reads_from_container(sample_fastq_container, max_genomes=2)
    
    summary = pseudo_alignment.get_summary()
    assert summary["Statistics"]["filtered_hr_kmers"] > 0  # Some k-mers should be removed


### **Test Combined Filtering (Read + K-mer + Max Genomes)**
def test_combined_filters(pseudo_alignment, sample_fastq_container):
    """Ensure the system correctly applies all three filters together."""
    pseudo_alignment.align_reads_from_container(
        sample_fastq_container, min_read_quality=40, min_kmer_quality=50, max_genomes=2
    )
    
    summary = pseudo_alignment.get_summary()
    assert summary["Statistics"]["filtered_quality_reads"] == 1  # Read2 should be filtered
    assert summary["Statistics"]["filtered_quality_kmers"] == 1  # Some k-mers should be removed
    assert summary["Statistics"]["filtered_hr_kmers"] == 5  # Some k-mers should be removed


### **Test Pseudo-Alignment with Filters Applied**
def test_pseudo_alignment_with_quality_and_max_genomes(pseudo_alignment, sample_fastq_container):
    """Ensure alignment still works correctly when filtering is applied."""
    pseudo_alignment.align_reads_from_container(
        sample_fastq_container, min_read_quality=30, min_kmer_quality=30, max_genomes=3
    )
    
    summary = pseudo_alignment.get_summary()
    assert summary["Statistics"]["unique_mapped_reads"] == 0 
    assert summary["Statistics"]["ambiguous_mapped_reads"] == 2 
    assert summary["Statistics"]["unmapped_reads"] == 1 
    assert summary["Statistics"]["filtered_quality_reads"] == 0
    assert summary["Statistics"]["filtered_quality_kmers"] == 0 
    assert summary["Statistics"]["filtered_hr_kmers"] == 0
