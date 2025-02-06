"""
@file test_kmer.py
@brief Test suite for KmerReference, Read, and PseudoAlignment functionalities.
"""

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

# Global constant for parametrized tests.
TIMES_TEST_TEST_PSEUDO_ALIGN_RANDOM_BIG_KMERS = 10

## ===========================================================
## Fixtures
## ===========================================================

# --- FASTA and FASTQ Fixtures ---

@pytest.fixture
def sample_fasta_container():
    """
    @brief Provides a sample FASTA container.
    @details Contains four genomes: Genome1, Genome2, Genome3, and Genome4.
    """
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
    """
    @brief Creates a KmerReference instance with k=3.
    """
    return KmerReference(3, sample_fasta_container)


@pytest.fixture
def sample_fastq_container():
    """
    @brief Provides a sample FASTQ container.
    @details Contains three reads with different quality scores.
    """
    data = (
        "@Read1\nAGCTAGCT\n+\nIIIIIIII\n"  # High quality
        "@Read2\nTGCATGCA\n+\n!!!!!!!!\n"  # Low quality
        "@Read3\nGGGGGGGG\n+\n!!IIIIII\n"   # Mixed
    )
    container = FASTAQRecordContainer()
    container.parse_records(data)
    return container


@pytest.fixture
def sample_reads(sample_fastq_container):
    """
    @brief Creates a list of Read objects from sample_fastq_container.
    """
    return [Read(record) for record in sample_fastq_container]


# --- Fixtures for Specific Read Mapping Cases ---

@pytest.fixture
def fasta_container_unmapped():
    """
    @brief Provides a FASTA container with genomes that have no matching k-mers.
    """
    data = (
        ">Genome1\nAACCGGTTAACC\n"
        ">Genome2\nGGTTCCAAGGTT\n"
    )
    container = FASTARecordContainer()
    container.parse_records(data)
    return container


@pytest.fixture
def kmer_reference_unmapped(fasta_container_unmapped):
    """
    @brief Creates a KmerReference (k=4) from fasta_container_unmapped.
    """
    return KmerReference(4, fasta_container_unmapped)


@pytest.fixture
def read_unmapped():
    """
    @brief Provides a Read expected to be unmapped.
    """
    data = "@Read1\nTAGGCAT\n+\nIIIIIII\n"
    container = FASTAQRecordContainer()
    container.parse_records(data)
    return Read(list(container)[0])


@pytest.fixture
def fasta_container_unique():
    """
    @brief Provides a FASTA container with two genomes that have distinct k-mers.
    """
    data = ">Genome1\nATGGCTATGCTA\n>Genome2\nCTATGGCAGGCA\n"
    container = FASTARecordContainer()
    container.parse_records(data)
    return container


@pytest.fixture
def kmer_reference_unique(fasta_container_unique):
    """
    @brief Creates a KmerReference (k=4) from fasta_container_unique.
    """
    return KmerReference(4, fasta_container_unique)


@pytest.fixture
def read_unique():
    """
    @brief Provides a Read expected to uniquely map.
    """
    data = "@Read2\nATGGCTAT\n+\nIIIIIIII\n"
    container = FASTAQRecordContainer()
    container.parse_records(data)
    return Read(list(container)[0])


@pytest.fixture
def fasta_container_ambiguous():
    """
    @brief Provides a FASTA container with four genomes to induce ambiguous mapping.
    """
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
    """
    @brief Creates a KmerReference (k=4) from fasta_container_ambiguous.
    """
    return KmerReference(4, fasta_container_ambiguous)


@pytest.fixture
def read_ambiguous():
    """
    @brief Provides a Read expected to map ambiguously.
    """
    data = "@Read3\nATCGATCCTAG\n+\nIIIIIIIIIII\n"
    container = FASTAQRecordContainer()
    container.parse_records(data)
    return Read(list(container)[0])


@pytest.fixture
def fasta_container_initially_unique():
    """
    @brief Provides a FASTA container for a case where a read is initially unique.
    """
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
    """
    @brief Creates a KmerReference (k=4) from fasta_container_initially_unique.
    """
    return KmerReference(4, fasta_container_initially_unique)


@pytest.fixture
def read_initially_unique():
    """
    @brief Provides a Read that is initially unique but may become ambiguous.
    """
    data = "@Read4\nATGCCGGGGCTAA\n+\nIIIIIIIIIIIII\n"
    container = FASTAQRecordContainer()
    container.parse_records(data)
    return Read(list(container)[0])


@pytest.fixture
def big_read():
    """
    @brief Provides a large Read used for synthetic k-mer tests.
    """
    data = "@BigRead\nAGCTAGCTAGAGGTCCTAATCCTAGCTAGCTAGCTAGCTAGCTAGCTGGTCATCAAAACCTTT\n+\n" \
           "IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII\n"
    container = FASTAQRecordContainer()
    container.parse_records(data)
    return Read(list(container)[0])


@pytest.fixture
def synthetic_kmer_reference(big_read):
    """
    @brief Creates a synthetic KmerReference from a big read.
    @details Extracts 31-mers from big_read, randomly distributes them among 4 genomes,
             and constructs a KmerReference.
    """
    k = 31
    extracted_kmers = list(extract_kmers_from_genome(k, big_read._Read__raw_read))
    num_genomes = 4
    genomes = {f">Genome{i+1}": [] for i in range(num_genomes)}
    for _, kmer in extracted_kmers:
        selected_genomes = random.sample(list(genomes.keys()), k=random.randint(1, num_genomes))
        for genome in selected_genomes:
            genomes[genome].append(kmer)
    fasta_data = "".join(f"{name}\n{'NN'.join(seq)}\n" for name, seq in genomes.items())
    for name, sequence in genomes.items():
        print(f"{name}: {sequence}")
    container = FASTARecordContainer()
    container.parse_records(fasta_data)
    return KmerReference(k, container)


@pytest.fixture
def create_read():
    """
    @brief Returns a factory function to create Read objects from FASTQ data.
    """
    def _create_read(identifier, sequence, quality, genomes_mapped_to=None, mapping_type="UNMAPPED"):
        fastaq_data = f"@{identifier}\n{sequence}\n+\n{quality}\n"
        container = FASTAQRecordContainer()
        container.parse_records(fastaq_data)
        read = Read(list(container)[0])
        if genomes_mapped_to is not None:
            read.mapping = ReadMapping(ReadMappingType[mapping_type], genomes_mapped_to)
        return read
    return _create_read


@pytest.fixture
def create_record():
    """
    @brief Returns a factory function to create FASTA Record objects.
    """
    def _create_record(genome_id, genome_sequence="AGCTAGCTAG"):
        fasta_data = f">{genome_id}\n{genome_sequence}\n"
        container = FASTARecordContainer()
        container.parse_records(fasta_data)
        return list(container)[0]
    return _create_record


@pytest.fixture
def pseudo_alignment(kmer_reference):
    """
    @brief Creates a PseudoAlignment instance using a real KmerReference.
    """
    return PseudoAlignment(kmer_reference)


## ===========================================================
## Basic Function Tests
## ===========================================================

def test_extract_kmers_from_genome():
    """
    @brief Verify extraction of k-mers from a genome string.
    @details Checks that both positions and sequences are as expected.
    """
    genome = "AGCTAGCTAGCT"
    k = 3
    expected_kmers = [
        (0, "AGC"), (1, "GCT"), (2, "CTA"), (3, "TAG"),
        (4, "AGC"), (5, "GCT"), (6, "CTA"), (7, "TAG"),
        (8, "AGC"), (9, "GCT"),
    ]
    assert list(extract_kmers_from_genome(k, genome)) == expected_kmers


def test_kmer_reference_building(kmer_reference):
    """
    @brief Verify that the k-mer reference is built correctly.
    @details Checks for expected k-mer presence and absence.
    """
    assert kmer_reference.get_kmer_references("AGC")
    assert kmer_reference.get_kmer_references("TGC")
    assert not kmer_reference.get_kmer_references("GGG")
    assert kmer_reference.get_kmer_references("GCT")
    assert kmer_reference.get_kmer_references("CCG")


## ===========================================================
## Read Mapping Tests
## ===========================================================

def test_unmapped_read(kmer_reference_unmapped, read_unmapped):
    """
    @brief Verify that a read with no matching k-mers is unmapped.
    """
    assert read_unmapped.pseudo_align(kmer_reference_unmapped) == ReadMappingType.UNMAPPED


def test_uniquely_mapped_read(kmer_reference_unique, read_unique):
    """
    @brief Verify that a read maps uniquely.
    """
    assert read_unique.pseudo_align(kmer_reference_unique) == ReadMappingType.UNIQUELY_MAPPED


def test_ambiguously_mapped_read(kmer_reference_ambiguous, read_ambiguous):
    """
    @brief Verify that a read maps ambiguously.
    """
    result = read_ambiguous.pseudo_align(kmer_reference_ambiguous)
    print("Ambiguous Read Alignment Result:", result)
    assert result == ReadMappingType.AMBIGUOUSLY_MAPPED
    assert len(read_ambiguous.mapping.genomes_mapped_to) == 4


def test_initially_unique_but_changes_to_ambiguous_tests_default_p(kmer_reference_initially_unique, read_initially_unique):
    """
    @brief Verify that a read initially mapping uniquely becomes ambiguous with default parameters.
    """
    result = read_initially_unique.pseudo_align(kmer_reference_initially_unique)
    print("Initially Unique But Changed to Ambiguous Result:", result)
    assert result == ReadMappingType.AMBIGUOUSLY_MAPPED


def test_initially_unique_but_changes_to_ambiguous_tests_big_p(kmer_reference_initially_unique, read_initially_unique):
    """
    @brief Verify that a read remains uniquely mapped when using a high p-threshold.
    """
    result = read_initially_unique.pseudo_align(kmer_reference_initially_unique, p=5)
    print("Initially Unique With Big p-threshold Result:", result)
    assert result == ReadMappingType.UNIQUELY_MAPPED


@pytest.mark.parametrize("execution_number", range(TIMES_TEST_TEST_PSEUDO_ALIGN_RANDOM_BIG_KMERS))
def test_pseudo_align_random_big_kmers(capsys, big_read, synthetic_kmer_reference, execution_number):
    """
    @brief Test pseudo-alignment using synthetic k-mers.
    @details Recomputes specific and total k-mer counts and compares with the alignment result.
    """
    m = 1
    p = 1
    specific_counts = {}
    total_counts = {}

    # Debug print of all k-mer keys.
    for k in synthetic_kmer_reference.kmers.keys():
        print(k)
    # Compute specific and total k-mer counts.
    for kmer, genome_map in synthetic_kmer_reference.kmers.items():
        mapped_genomes = list(genome_map.keys())
        if len(mapped_genomes) == 1:
            specific_counts[mapped_genomes[0]] = specific_counts.get(mapped_genomes[0], 0) + len(genome_map[mapped_genomes[0]])
        for genome, positions in genome_map.items():
            total_counts[genome] = total_counts.get(genome, 0) + len(positions)

    print("Specific k-mer counts per genome:", specific_counts)
    print("Total k-mer counts per genome:", total_counts)

    sorted_genomes = sorted(specific_counts, key=specific_counts.get, reverse=True)
    if not sorted_genomes:
        alignment_result = big_read.pseudo_align(synthetic_kmer_reference, p=p, m=m, debug=True)
        captured = capsys.readouterr()
        print(captured.out)
        assert alignment_result == ReadMappingType.AMBIGUOUSLY_MAPPED
        return

    max_genome = sorted_genomes[0]
    second_max_genome = sorted_genomes[1] if len(sorted_genomes) > 1 else None
    specific_max_count = specific_counts[max_genome]
    specific_second_max_count = specific_counts[second_max_genome] if second_max_genome else 0
    m_difference = specific_max_count - specific_second_max_count

    alignment_result = big_read.pseudo_align(synthetic_kmer_reference, p=p, m=m, debug=True)
    captured = capsys.readouterr()

    max_total_kmers = max(total_counts.values())
    map_total_kmers = total_counts[max_genome]
    p_difference = max_total_kmers - map_total_kmers

    print(captured.out)
    print("Specific MaxCount:", specific_max_count)
    print("Second Specific Max Count:", specific_second_max_count)
    print("Difference:", m_difference)
    print("MaxCount:", max_total_kmers)
    print("MapCount:", map_total_kmers)
    print("Alignment Result:", alignment_result)

    if (m_difference >= m) and (not p_difference > p):
        assert alignment_result == ReadMappingType.UNIQUELY_MAPPED
    else:
        assert alignment_result == ReadMappingType.AMBIGUOUSLY_MAPPED


## ===========================================================
## PseudoAlignment Tests
## ===========================================================

def test_add_read_for_pseudo_alignment(pseudo_alignment, create_read, create_record):
    """
    @brief Test adding a read to PseudoAlignment.
    """
    genome1 = create_record("Genome1")
    read1 = create_read("read1", "AGCTAGCT", "IIIIIIII", [genome1], "UNIQUELY_MAPPED")
    pseudo_alignment.add_read(read1)
    assert "read1" in pseudo_alignment.reads
    assert pseudo_alignment.reads["read1"]["mapping_type"] == ReadMappingType.UNIQUELY_MAPPED
    assert pseudo_alignment.reads["read1"]["genomes_mapped_to"] == [genome1.identifier]


def test_add_duplicate_read_raises_exception_for_pseudo_alignment(pseudo_alignment, create_read, create_record):
    """
    @brief Test that adding a duplicate read raises an exception.
    """
    genome1 = create_record("Genome1")
    read1 = create_read("read1", "AGCTAGCT", "IIIIIIII", [genome1], "UNIQUELY_MAPPED")
    pseudo_alignment.add_read(read1)
    with pytest.raises(AddingExistingRead):
        pseudo_alignment.add_read(read1)


def test_get_summary_for_pseudo_alignment(pseudo_alignment, create_read, create_record):
    """
    @brief Test that PseudoAlignment summary returns correct statistics.
    """
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
    """
    @brief Test saving and loading of PseudoAlignment.
    """
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


## ===========================================================
## Quality Filtering Tests Extension (4.1)
## ===========================================================

def test_min_read_quality_filter(pseudo_alignment, sample_fastq_container):
    """
    @brief Test that reads below the minimum read quality are filtered.
    """
    pseudo_alignment.align_reads_from_container(sample_fastq_container, min_read_quality=40)
    summary = pseudo_alignment.get_summary()
    print(summary)
    assert summary["Statistics"]["filtered_quality_reads"] == 1
    assert summary["Statistics"]["unmapped_reads"] >= 1


def test_min_kmer_quality_filter(pseudo_alignment, sample_fastq_container):
    """
    @brief Test that k-mers below the minimum k-mer quality are filtered.
    """
    pseudo_alignment.align_reads_from_container(sample_fastq_container, min_kmer_quality=60)
    summary = pseudo_alignment.get_summary()
    print(summary)
    assert summary["Statistics"]["filtered_quality_kmers"] > 0


def test_max_genomes_filter(pseudo_alignment, sample_fastq_container):
    """
    @brief Test that k-mers mapping to too many genomes are filtered.
    """
    pseudo_alignment.align_reads_from_container(sample_fastq_container, max_genomes=2)
    summary = pseudo_alignment.get_summary()
    assert summary["Statistics"]["filtered_hr_kmers"] > 0


def test_combined_filters(pseudo_alignment, sample_fastq_container):
    """
    @brief Test combined filtering: read quality, k-mer quality, and max genomes.
    """
    pseudo_alignment.align_reads_from_container(sample_fastq_container, min_read_quality=40, min_kmer_quality=50, max_genomes=2)
    summary = pseudo_alignment.get_summary()
    assert summary["Statistics"]["filtered_quality_reads"] == 1
    assert summary["Statistics"]["filtered_quality_kmers"] == 1
    assert summary["Statistics"]["filtered_hr_kmers"] == 5


def test_pseudo_alignment_with_quality_and_max_genomes(pseudo_alignment, sample_fastq_container):
    """
    @brief Test pseudo-alignment when both quality and genome filters are applied.
    """
    pseudo_alignment.align_reads_from_container(sample_fastq_container, min_read_quality=30, min_kmer_quality=30, max_genomes=3)
    summary = pseudo_alignment.get_summary()
    assert summary["Statistics"]["unique_mapped_reads"] == 0 
    assert summary["Statistics"]["ambiguous_mapped_reads"] == 2 
    assert summary["Statistics"]["unmapped_reads"] == 1 
    assert summary["Statistics"]["filtered_quality_reads"] == 0
    assert summary["Statistics"]["filtered_quality_kmers"] == 0 
    assert summary["Statistics"]["filtered_hr_kmers"] == 0


## ===========================================================
## KmerReference Filtering Tests (Extension 4.4)
## ===========================================================

def test_filter_similar_genomes():
    """
    @brief Test that highly similar genomes are filtered out when filtering is enabled.
    """
    data = (
        ">GenomeA\nAGCTAGCTAGCT\n"
        ">GenomeB\nAGCTAGCTAGCT\n"  # Identical to GenomeA; should be filtered.
        ">GenomeC\nTGCATGCATGCA\n"   # Different; should be kept.
    )
    container = FASTARecordContainer()
    container.parse_records(data)
    kref = KmerReference(4, container, filter_similar=True, similarity_threshold=0.95)
    assert hasattr(kref, "similarity_info")
    sim_info = kref.similarity_info
    kept_ids = {genome.identifier for genome in kref.genomes}
    assert "GenomeC" in kept_ids
    if "GenomeA" in kept_ids:
        assert "GenomeB" not in kept_ids
        kept_genome = "GenomeA"
        filtered_genome = "GenomeB"
    else:
        assert "GenomeB" in kept_ids
        assert "GenomeA" not in kept_ids
        kept_genome = "GenomeB"
        filtered_genome = "GenomeA"
    assert sim_info[filtered_genome]["kept"] == "no"
    assert sim_info[filtered_genome]["similar_to"] == kept_genome
    for genome_id in kept_ids:
        assert sim_info[genome_id]["kept"] == "yes"
        assert sim_info[genome_id]["similar_to"] == "NA"


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
