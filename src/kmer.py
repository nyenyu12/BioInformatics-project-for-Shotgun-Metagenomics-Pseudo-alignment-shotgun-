import constants
import json
import gzip
import pickle
from collections import namedtuple, defaultdict
from enum import Enum
from records import FASTAQRecordContainer, FASTARecordContainer, Record
from typing import List, Dict, Iterator, Set, Tuple, Union


IGNORE_AMBIGUOUS_THRESHOLD = 0
M_THRESHOLD = 0


class NotValidatingUniqueMapping(Exception):
    def __init__(self, message):
        super().__init__(message)


class AddingExistingRead(Exception):
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
ReadMapping = namedtuple("ReadMapping", ["type", "genomes_mapped_to"])


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


def extract_kmers_from_genome(k: int, genome: str) -> Iterator[Tuple[int, str]]:
    if k > len(genome) or k <= 0:
        return iter([])

    for i in range(len(genome) - k + 1):
        yield (i, genome[i : i + k])


class KmerReference(object):
    def __init__(self, k: int, fasta_record_container: FASTARecordContainer):
        self.kmer_len: int = k
        self.kmers: Dict[str, Dict[Record, Set[int]]] = {}

        for fasta_record in fasta_record_container:
            genome_sequence = fasta_record["genome"]
            genome_name = (
                fasta_record  # Assuming fasta_record is uniquely identifying the genome
            )

            for kmer_pos, kmer in extract_kmers_from_genome(k, genome_sequence):
                if constants.NULL_NUCLEOTIDES_CHAR not in kmer:
                    if kmer not in self.kmers:
                        self.kmers[kmer] = {}

                    if genome_name not in self.kmers[kmer]:
                        self.kmers[kmer][genome_name] = set()

                    self.kmers[kmer][genome_name].add(kmer_pos)

    def __getitem__(self, kmer):
        return self.kmers.get(kmer, None)

    def get_kmer_references(self, kmer: str) -> Dict[Record, Set[int]]:
        return self.kmers.get(kmer, {})

    def get_summary(self) -> Dict[str, Dict[str, Union[Dict[str, List[int]], Dict[str, int]]]]:
        kmer_details = {
            kmer: {
                genome["description"]: sorted(positions)
                for genome, positions in genomes.items()
            }
            for kmer, genomes in self.kmers.items()
        }

        genome_summary: defaultdict[str, Dict[str, int]] = defaultdict(
            lambda: {"total_bases": 0, "unique_kmers": 0, "multi_mapping_kmers": 0}
        )
        genome_kmer_count = defaultdict(set)

        for kmer, genomes in self.kmers.items():
            for genome, positions in genomes.items():
                description = genome["description"]
                genome_summary[description]["total_bases"] = len(genome["genome"])
                genome_kmer_count[description].add(kmer)

        for genome in genome_kmer_count:
            unique_kmers = sum(
                1 for kmer in genome_kmer_count[genome] if len(self.kmers[kmer]) == 1
            )
            multi_mapping_kmers = len(genome_kmer_count[genome]) - unique_kmers
            genome_summary[genome]["unique_kmers"] = unique_kmers
            genome_summary[genome]["multi_mapping_kmers"] = multi_mapping_kmers

        return {"Kmers": kmer_details, "Summary": dict(genome_summary)}


class Read:
    def __init__(self, fastaq_record: Record):
        self.identifier = fastaq_record["identifier"]
        self.mapping = ReadMapping(ReadMappingType.UNMAPPED, None)
        self.kmers: Dict[str, ReadKmer] = {}
        self.__raw_read = fastaq_record["sequence"]
        self.__genomes_map_count: Dict[Record, int] = None

    def __str__(self):
        rows = [f"Mapping: {self.mapping}"]
        for kmer in self.kmers.keys():
            rows.append(f"k-mer: {kmer}")
            rows.append(f"specifity: {self.kmers[kmer].specifity}")
            rows.append(f"Genome References:")
            for references in self.kmers[kmer]:
                rows.append(f"\t{references}")
        return "\n".join(rows)

    def __repr__(self):
        return self.__str__()

    def extract_kmer_references(self, kmer_reference: KmerReference):
        for _, kmer in extract_kmers_from_genome(
            kmer_reference.kmer_len, self.__raw_read
        ):
            single_kmer_reference = kmer_reference.get_kmer_references(kmer)
            if single_kmer_reference:
                kmer_specifity = (
                    KmerSpecifity.SPECIFIC
                    if len(single_kmer_reference) == 1
                    else KmerSpecifity.UNSPECIFIC
                )
                self.kmers[kmer] = ReadKmer(
                    specifity=kmer_specifity, references=single_kmer_reference
                )

    def generate_genome_counts(self, map_count: bool=False) -> Dict[Record, int]:
        genomes_and_count: Dict[Record, int] = {}
        for kmer in self.kmers.values():
            if not map_count or kmer.specifity == KmerSpecifity.SPECIFIC:
                for genome in kmer.references.keys():
                    genomes_and_count[genome] = genomes_and_count.get(genome, 0) + 1
        return genomes_and_count

    def try_to_align_specific(self, m: int):
        if m < 0:
            raise ValueError("m must be non-negative")

        self.__genomes_map_count = self.generate_genome_counts(map_count=True)
        if len(self.__genomes_map_count) == 1:
            self.mapping = ReadMapping(
                ReadMappingType.UNIQUELY_MAPPED,
                [list(self.__genomes_map_count.keys())[0]],
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
                    ReadMappingType.UNIQUELY_MAPPED, [sorted_genomes[0]]
                )
                return True
        self.mapping = ReadMapping(
            ReadMappingType.AMBIGUOUSLY_MAPPED, list(self.__genomes_map_count.keys())
        )
        return False

    def validate_unique_mappings(self, p: int):
        if (
            self.mapping.type != ReadMappingType.UNIQUELY_MAPPED
            or p < IGNORE_AMBIGUOUS_THRESHOLD
        ):
            return  # Skip validation if not uniquely mapped

        genome_total_counts = self.generate_genome_counts(map_count=False)
        mapped_genome = self.mapping.genomes_mapped_to[0]
        max_total_kmers = max(genome_total_counts.values(), default=0)
        mapped_genome_kmers = genome_total_counts.get(mapped_genome, 0)

        if max_total_kmers - mapped_genome_kmers > p:
            ambiguous_genomes = [mapped_genome]
            for genome, count in genome_total_counts.items():
                if count >= mapped_genome_kmers:
                    ambiguous_genomes.append(genome)
            self.mapping = ReadMapping(
                ReadMappingType.AMBIGUOUSLY_MAPPED, ambiguous_genomes
            )

    def pseudo_align(self, kmer_reference: KmerReference, p: int=1, m: int=1, debug: bool=False) -> ReadMappingType:
        if not (
            isinstance(kmer_reference, KmerReference)
            and isinstance(p, int)
            and isinstance(m, int)
            and isinstance(debug, bool)
        ):
            raise TypeError(
                f"Invalid types given to pseudo align: {type(kmer_reference)}, {type(p)}, {type(m)}, {type(debug)}"
            )
        if m < M_THRESHOLD:
            raise ValueError(f"m must be bigger than or equal to {M_THRESHOLD}")

        self.extract_kmer_references(kmer_reference)
        if not self.kmers:
            return ReadMappingType.UNMAPPED
        if self.try_to_align_specific(m):
            if debug:
                print(
                    f"[DEBUG pseudo_align]: After try_to_align_specific self.mapping: {self.mapping.type}"
                )
            self.validate_unique_mappings(p)
            return self.mapping.type
        elif debug:
            print(
                f"[DEBUG pseudo_align]: After try_to_align_specific self.mapping: {self.mapping.type}"
            )
        return ReadMappingType.AMBIGUOUSLY_MAPPED


class PseudoAlignment:
    def __init__(self, kmer_reference: KmerReference):
        self.kmer_reference = kmer_reference
        self.reads: Dict[str, Dict[str, Union[ReadMappingType, List[str]]]] = {}

    def add_read(self, read: Read):
        if read.identifier in self.reads:
            raise AddingExistingRead(
                f"There already exists a read with identifier: {read.identifier}"
            )

        self.reads[read.identifier] = {
            "mapping_type": read.mapping.type, 
            "genomes_mapped_to": [genome.identifier for genome in read.mapping.genomes_mapped_to],
        }

    def add_read_from_read_record(self, read_record: Record, m: int=1, p: int=1):
        read = Read(read_record)
        read.pseudo_align(self.kmer_reference, m=m, p=p)
        self.add_read(read)

    def align_reads_from_container(self, reads_container: FASTAQRecordContainer, m: int=1, p: int=1):
        for read_record in reads_container:
            self.add_read_from_read_record(read_record, m=1, p=1)

    def save(self, align_file: str):
        with gzip.open(align_file, "wb") as f:
            pickle.dump(self, f)

    def get_summary(self) -> Dict[str, Dict[str, Union[int, Dict[str, int]]]]:
        summary = {
            "unique_mapped_reads": 0,
            "ambiguous_mapped_reads": 0,
            "unmapped_reads": 0,
        }
        genome_mapping: Dict[str, Dict[str, int]] = {}

        for _, read_alignment_details in self.reads.items():
            mapping_type = read_alignment_details["mapping_type"]
            genomes_mapped_to = read_alignment_details["genomes_mapped_to"]
            
            if mapping_type != ReadMappingType.UNMAPPED:
                read_type_string: str = None
                
                if mapping_type == ReadMappingType.UNIQUELY_MAPPED:
                    summary["unique_mapped_reads"] += 1
                    read_type_string = "unique_reads"
                    
                elif mapping_type == ReadMappingType.AMBIGUOUSLY_MAPPED:
                    summary["ambiguous_mapped_reads"] += 1
                    read_type_string = "ambiguous_reads"

                for genome in genomes_mapped_to:
                    genome_mapping.setdefault(
                        genome, {"unique_reads": 0, "ambiguous_reads": 0}
                    )
                    genome_mapping[genome][read_type_string] += 1
            else:
                summary["unmapped_reads"] += 1

        return {"Statistics": summary, "Summary": genome_mapping}

    def __repr__(self):
        return json.dumps(self.get_summary(), indent=4)
