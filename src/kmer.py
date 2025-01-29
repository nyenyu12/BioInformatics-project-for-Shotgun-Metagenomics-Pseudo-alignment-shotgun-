import constants
from collections import namedtuple
from enum import Enum


class NotValidatingUniqueMapping(Exception):
    def __init__(self, message):
        super().__init__(message)


class ReadMappingType(Enum):
    UNMAPPED = 1
    UNIQUELY_MAPPED = 2
    AMBIGUOUSLY_MAPPED = 3


class KmerSpecifity(Enum):
    SPECIFIC = 1
    UNSPECIFIC = 2


ReadKmer = namedtuple("ReadKmer", ["specifity", "references"])
ReadMapping = namedtuple("ReadMappinng", ["type", "genome_mapped_to"])


def extract_k_max_value_keys_from_dict(d, k):
    """
    Extracts the two keys with the maximum values from the dictionary.

    Args:
        d (dict): A dictionary with keys and values.
        k: Amount of keys to find
    Returns:
        list: A list containing the two keys with the highest values.
              If there are fewer than two items, it returns all keys.
    """
    if not isinstance(d, dict):
        raise ValueError("Input must be a dictionary.")
    if len(d) == 0:
        return []

    sorted_keys = sorted(d, key=d.get, reverse=True)

    return sorted_keys[:k]


def extract_k_mers_from_genome(k: int, genome: str):
    if k > len(genome) or k <= 0:
        return iter([])

    for i in range(len(genome) - k + 1):
        yield (i, genome[i : i + k])


class KmerReference(object):

    # TODO: not basic implementation, wtf is it
    def __init__(self, k, fasta_record_container):
        self.k_mer_len = k
        self.k_mers = {}

        for fasta_record in fasta_record_container:
            genome_sequence = fasta_record["genome"]
            for k_mer_pos, k_mer in extract_k_mers_from_genome(k, genome_sequence):
                if constants.NULL_NUCLEOTIDES_CHAR not in k_mer:
                    if fasta_record not in self.k_mers:
                        self.k_mers[fasta_record] = {}

                    if k_mer not in self.k_mers[fasta_record]:
                        self.k_mers[fasta_record][k_mer] = set()

                    self.k_mers[fasta_record][k_mer].add(k_mer_pos)

    def __getitem__(self, fasta_record):
        return self.k_mers.get(fasta_record, None)

    def get_k_mer_references(self, k_mer):
        return {
            record: positions[k_mer]
            for record, positions in self.k_mers.items()
            if k_mer in positions
        }

class Read:
    def __init__(self, fastaq_record):
        self.mapping = ReadMapping(ReadMappingType.UNMAPPED, None)
        self.k_mers = {}
        self.__raw_read = fastaq_record["sequence"]
        self.__genomes_map_count = None

    def extract_kmer_references(self, k_mer_reference):
        for _, k_mer in extract_k_mers_from_genome(
            k_mer_reference.k_mer_len, self.__raw_read
        ):
            single_k_mer_reference = k_mer_reference.get_k_mer_references(k_mer)
            if single_k_mer_reference:
                k_mer_specifity = (
                    KmerSpecifity.SPECIFIC
                    if len(single_k_mer_reference) == 1
                    else KmerSpecifity.UNSPECIFIC
                )
                self.k_mers[k_mer] = ReadKmer(
                    specifity=k_mer_specifity, references=single_k_mer_reference
                )

    def generate_genome_counts(self, map_count=False):
        genomes_and_count = {}
        for k_mer in self.k_mers.values():
            if not map_count or k_mer.specifity == KmerSpecifity.SPECIFIC:
                for genome in k_mer.references.keys():
                    genomes_and_count[genome] = genomes_and_count.get(genome, 0) + 1
        return genomes_and_count

    def try_to_align_specific(self, m):
        self.__genomes_map_count = self.generate_genome_counts(map_count=True)
        if len(self.__genomes_map_count) == 1:
            self.mapping = ReadMapping(
                ReadMappingType.UNIQUELY_MAPPED,
                list(self.__genomes_map_count.keys())[0],
            )
            return True
        if len(self.__genomes_map_count) > 1:
            sorted_genomes = sorted(
                self.__genomes_map_count, key=self.__genomes_map_count.get, reverse=True
            )
            if (
                self.__genomes_map_count[sorted_genomes[0]]
                >= self.__genomes_map_count[sorted_genomes[1]] + m
            ):
                self.mapping = ReadMapping(
                    ReadMappingType.UNIQUELY_MAPPED, sorted_genomes[0]
                )
                return True
        self.mapping = ReadMapping(
            ReadMappingType.AMBIGUOUSLY_MAPPED, list(self.__genomes_map_count.keys())
        )
        return False

    def pseudo_align(self, k_mer_reference, p=1, m=1):
        self.extract_kmer_references(k_mer_reference)
        if not self.k_mers:
            return ReadMappingType.UNMAPPED
        if self.try_to_align_specific(m):
            return self.mapping.type
        return ReadMappingType.AMBIGUOUSLY_MAPPED
