import pytest
import random
from collections import namedtuple
from records import FASTARecordContainer, FASTAQRecordContainer
from k_mer import Read, ReadMappingType, KmerReference, extract_k_mers_from_genome


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


def test_extract_k_mers_from_genome():
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
    assert list(extract_k_mers_from_genome(k, genome)) == expected_kmers


def test_kmer_reference_building(kmer_reference):
    assert kmer_reference.get_k_mer_references("AGC")
    assert kmer_reference.get_k_mer_references("TGC")
    assert not kmer_reference.get_k_mer_references(
        "GGG"
    )  # k-mer not present in reference genomes
    assert kmer_reference.get_k_mer_references(
        "GCT"
    )  # k-mer appearing in multiple genomes
    assert kmer_reference.get_k_mer_references("CCG")  # k-mer exclusive to Genome4


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
    extracted_kmers = list(extract_k_mers_from_genome(k, big_read._Read__raw_read))
    
    # Randomly distribute k-mers across genomes, allowing duplicates
    num_genomes = 4
    genomes = {f">Genome{i+1}": [] for i in range(num_genomes)}
    
    for _, k_mer in extracted_kmers:
        selected_genomes = random.sample(list(genomes.keys()), k=random.randint(1, num_genomes))
        for genome in selected_genomes:
            genomes[genome].append(k_mer)
    
    fasta_data = "".join(f"{name}\n{"NN".join(seq)}\n" for name, seq in genomes.items())
    
    for name, sequence in genomes.items():
        print(f"{name}: {sequence}")
    
    container = FASTARecordContainer()
    container.parse_records(fasta_data)
    return KmerReference(k, container)


def test_kmer_reference_counts(synthetic_kmer_reference):
    kmer_counts = {}
    unspecific_count = 0

    for kmer in synthetic_kmer_reference.k_mers.values():
        for genome, positions in kmer.items():
            kmer_counts[genome] = kmer_counts.get(genome, 0) + len(positions)

    print("Specific k-mer counts per genome:", kmer_counts)
    print("Unspecific k-mer count:", unspecific_count)

    assert sum(kmer_counts.values()) > 0  # Ensure k-mers are assigned

@pytest.mark.parametrize('execution_number', range(10))
def test_pseudo_align_random_big_kmers(capsys, big_read, synthetic_kmer_reference, execution_number):
    """
    This is meant to almost implement the algorithm a second time to verify the main implementation.
    In addition, because the spererator for the k-mers is 'NN', then this tests that as well :)
    """
    m = 1
    p = 1
    specific_counts = {}
    total_counts = {}

    for k in synthetic_kmer_reference.k_mers.keys():
        print (k)
    # Notice all k-mers are in the Read since the genomes are generated from k-mers in the read
    for kmer, genome_map in synthetic_kmer_reference.k_mers.items():
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
