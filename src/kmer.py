import constants
from collections import namedtuple
from enum import Enum
from data_files import FASTARecordContainer, FASTAFile


class NotValidatingUniqueMapping(Exception):
    def __init__(self, message, errors):
        super().__init__(message)


class ReadMappingType(Enum):
    UNMAPPED = 1
    UNIQUELY_MAPPED = 2
    AMBIGUOUSLY_MAPPED = 3


class KmerSpecifity(Enum):
    SPECIFIC = 1
    UNSPECIFIC = 2


ReadKmer = namedtuple("ReadedKmer", ["specifity", "references"])
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
        return

    for i in range(len(genome) - k + 1):
        yield (i, genome[i : i + k])


class Kmer(object):

    def __init__(self):
        pass


class KmerReference(object):

    def __init__(self, k, fasta_record_container):
        self.__k_mers = {}

        for fasta_record in fasta_record_container:
            for k_mer_pos, k_mer in extract_k_mers_from_genome(fasta_record.genome):
                # TODO: not basic implementation, wtf is it
                if not constants.NULL_NUCLEOTIDES_CHAR in k_mer:
                    if not k_mer in self.__k_mers:
                        self.__k_mers[k_mer] = {}

                    single_k_mer_references = self.__k_mers[k_mer]
                    if not fasta_record.genome in single_k_mer_references.keys():
                        single_k_mer_references[fasta_record.genome] = set(k_mer_pos)

                    single_k_mer_references[fasta_record.genome].add(k_mer_pos)

    def __getitem__(self, k_mer):
        if k_mer in self.__k_mers:
            return self.__k_mers[k_mer]
        return None


class Read(object):

    def __init__(self, fastaq_record):
        self.mapping = ReadMapping(ReadMappingType.UNMAPPED, None)
        self.__k_mers = {}
        self.__raw_read = fastaq_record.raw_read
        self.__genomes_map_count = None
        self.__genomes_max_count = None

    def extract_kmer_references(self, k_mer_reference):
        for _, k_mer in extract_k_mers_from_genome(self.__raw_read):
            single_k_mer_reference = k_mer_reference[k_mer]

            if single_k_mer_reference != None:
                k_mer_specifity = KmerSpecifity.SPECIFIC
                if len(single_k_mer_reference.keys()) > 1:
                    k_mer_specifity = KmerSpecifity.UNSPECIFIC

                self.__k_mers[k_mer] = ReadKmer(
                    specifity=k_mer_specifity, references=single_k_mer_reference
                )

    def generate_genome_counts(self, map_count=False):
        genomes_and_count = {}

        for k_mer in self.__k_mers.values():
            if map_count and k_mer.specifity == KmerSpecifity.SPECIFIC:
                for genome in k_mer.references.keys():
                    if genome in genomes_and_count:
                        genomes_and_count[genome] += 1
                    else:
                        genomes_and_count[genome] = 1

        return genomes_and_count

    def try_to_align_specific(self, m):
        """
        Returns True if mapped uniquely, False otherwise
        """
        self.__genomes_map_count = self.generate_genome_counts(map_count=True)

        if len(self.__genomes_map_count.keys()) == 1:
            self.mapping = ReadMapping(
                ReadMappingType.UNIQUELY_MAPPED,
                self.__genomes_map_count.keys()[0],
            )
            return True

        if len(self.__genomes_map_count.keys()) > 1:
            most_referred_genomes = extract_k_max_value_keys_from_dict(
                self.__genomes_map_count
            )

            if (
                self.__genomes_map_count[most_referred_genomes[0]]
                >= self.__genomes_map_count[most_referred_genomes[1]] + m
            ):
                self.mapping = ReadMapping(
                    ReadMappingType.UNIQUELY_MAPPED,
                    most_referred_genomes[0],
                )
                return True

        for genome in self.__genomes_map_count.keys():
            self.__genomes_map_count[genome] += 1

        self.mapping = ReadMapping(
            ReadMappingType.UNIQUELY_MAPPED, self.__genomes_map_count
        )
        return False

    def validate_unique_mappings(self, p):
        if self.mapping.type != ReadMappingType.UNIQUELY_MAPPED:
            raise NotValidatingUniqueMapping

        self.__genome_max_counts = self.generate_genome_counts()
        sorted_max_count_genomes = extract_k_max_value_keys_from_dict(
            self.__genome_max_counts, len(self.__genome_max_counts)
        )

        if self.__genome_max_counts[sorted_max_count_genomes[0]] > (
            self.__genomes_map_count[self.mapping.genome_mapped_to] + p
        ):
            self.mapping.type = ReadMappingType.AMBIGUOUSLY_MAPPED
            for genome in sorted_max_count_genomes:
                if (
                    self.__genome_max_counts[genome]
                    >= self.__genomes_map_count[self.mapping.genome_mapped_to]
                ):
                    self.__genome_max_counts[genome] += 1

    def pseudo_align(self, k_mer_reference, p=1, m=1):
        self.extract_kmer_references(self, k_mer_reference)

        if len(self.__k_mers.keys()) == 0:
            return ReadMappingType.UNMAPPED

        if self.try_to_align_specific(m):
            self.validate_unique_mappings(p)
