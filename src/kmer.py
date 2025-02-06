"""
@file kmer.py
@brief Core bioinformatics functionality for k-mer reference construction and pseudo-alignment.
"""

import constants
import json
import gzip
import pickle
from collections import namedtuple, defaultdict
from enum import Enum
from records import FASTAQRecordContainer, FASTARecordContainer, Record
from typing import List, Dict, Iterator, Set, Tuple, Union, Optional

"""Global thresholds"""
IGNORE_AMBIGUOUS_THRESHOLD = 0
M_THRESHOLD = 0

##===========================================================
## Custom Exceptions
##===========================================================

class NotValidatingUniqueMapping(Exception):
    """
    @brief Exception raised when unique mapping validation fails.
    """
    def __init__(self, message: str) -> None:
        super().__init__(message)

class AddingExistingRead(Exception):
    """
    @brief Exception raised when adding a read that already exists.
    """
    def __init__(self, message: str) -> None:
        super().__init__(message)

##===========================================================
## Enumerations
##===========================================================

class ReadMappingType(Enum):
    """
    @brief Represents the type of mapping for a read.
    """
    UNMAPPED = 1
    UNIQUELY_MAPPED = 2
    AMBIGUOUSLY_MAPPED = 3

class KmerSpecifity(Enum):
    """
    @brief Indicates whether a k-mer is specific to one genome or unspecific.
    """
    SPECIFIC = 1
    UNSPECIFIC = 2

##===========================================================
## Named Tuples
##===========================================================

ReadKmer = namedtuple("ReadKmer", ["specifity", "references"])
"""@brief A tuple representing a read's k-mer with its specificity and references."""

ReadMapping = namedtuple("ReadMapping", ["type", "genomes_mapped_to"])
"""@brief A tuple representing a read's mapping result."""

"""===========================================================
Helper Functions
============================================================"""

def extract_k_max_value_keys_from_dict(d: Dict[str, int], k: int) -> List[str]:
    """
    @brief Extracts up to k keys with the maximum values from a dictionary.
    @param d Dictionary with keys and values.
    @param k Number of keys to extract.
    @return List containing up to k keys with the highest values.
    """
    if not isinstance(d, dict):
        raise ValueError("Input must be a dictionary.")
    if len(d) == 0:
        return []
    sorted_keys = sorted(d, key=d.get, reverse=True)
    return sorted_keys[:k]

def extract_kmers_from_genome(k: int, genome: str) -> Iterator[Tuple[int, str]]:
    """
    @brief Extracts all k-mers from a given genome sequence.
    @param k Length of the k-mer.
    @param genome Genome sequence as a string.
    @return Iterator of tuples (position, k-mer string).
    """
    if k > len(genome) or k <= 0:
        return iter([])
    for i in range(len(genome) - k + 1):
        yield (i, genome[i: i + k])

def reverse_complement(seq: str) -> str:
    """
    @brief Returns the reverse complement of a nucleotide sequence.
    @param seq Input nucleotide sequence.
    @return Reverse complement of the input sequence.
    """
    complement = str.maketrans("ACGT", "TGCA")
    return seq.translate(complement)[::-1]

##===========================================================
## Class: KmerReference
##===========================================================

class KmerReference(object):
    """
    @brief Constructs a k-mer reference database from FASTA records.
    """
    def __init__(self,
                 k: int,
                 fasta_record_container: FASTARecordContainer,
                 filter_similar: bool = False,
                 similarity_threshold: float = 0.95) -> None:
        """
        @brief Initializes the k-mer reference database.
        @param k Length of the k-mer.
        @param fasta_record_container Container of genome records (FASTARecordContainer).
        @param filter_similar If True, applies similarity filtering to remove highly similar genomes.
        @param similarity_threshold Threshold (between 0 and 1) for considering genomes similar.
        """
        if filter_similar and not (0 <= similarity_threshold <= 1):
            raise ValueError("similarity_threshold must be between 0 and 1")
        fasta_records = list(fasta_record_container)
        self.genomes = fasta_records
        self.kmer_len: int = k
        self.kmers: Dict[str, Dict[Record, Set[int]]] = {}
        self._build_kmer_mapping(fasta_records, k)
        if filter_similar:
            self._filter_similar_genomes(similarity_threshold)

    def _build_kmer_mapping(self, fasta_records: List[Record], k: int) -> None:
        """
        @brief Builds the internal k-mer mapping from FASTA records.
        @param fasta_records List of genome records.
        @param k Length of the k-mer.
        """
        for fasta_record in fasta_records:
            genome_sequence = fasta_record["genome"]
            # Using the record itself as a unique key.
            genome_record = fasta_record
            for kmer_pos, kmer in extract_kmers_from_genome(k, genome_sequence):
                if constants.NULL_NUCLEOTIDES_CHAR not in kmer:
                    if kmer not in self.kmers:
                        self.kmers[kmer] = {}
                    if genome_record not in self.kmers[kmer]:
                        self.kmers[kmer][genome_record] = set()
                    self.kmers[kmer][genome_record].add(kmer_pos)

    def _compute_genome_stats(self) -> Tuple[Dict[str, Dict[str, Union[int, float]]], Dict[str, Set[str]]]:
        """
        @brief Computes per-genome statistics and maps genome IDs to their k-mers.
        @return A tuple (genome_stats, genome_to_kmers).
                genome_stats: Dictionary mapping genome_id to its statistics.
                genome_to_kmers: Dictionary mapping genome_id to the set of k-mers found.
        """
        genome_to_kmers: Dict[str, Set[str]] = {}
        for kmer, mapping in self.kmers.items():
            for genome_record in mapping.keys():
                genome_id = genome_record.identifier
                genome_to_kmers.setdefault(genome_id, set()).add(kmer)
        genome_stats: Dict[str, Dict[str, Union[int, float]]] = {}
        for order, genome in enumerate(self.genomes):
            genome_id = genome.identifier
            kmers_set = genome_to_kmers.get(genome_id, set())
            total_kmers = len(kmers_set)
            unique_kmers = sum(1 for kmer in kmers_set if len(self.kmers[kmer]) == 1)
            genome_length = len(genome["genome"])
            genome_stats[genome_id] = {
                "unique_kmers": unique_kmers,
                "total_kmers": total_kmers,
                "genome_length": genome_length,
                "order": order,
            }
        return genome_stats, genome_to_kmers

    def _sort_genomes_for_filtering(self, genome_stats: Dict[str, Dict[str, Union[int, float]]]) -> List[Tuple[str, Dict[str, Union[int, float]]]]:
        """
        @brief Sorts genomes based on uniqueness, total k-mers, genome length, and input order.
        @param genome_stats Dictionary of genome statistics.
        @return A sorted list of tuples (genome_id, stats).
        """
        return sorted(genome_stats.items(),
                      key=lambda x: (x[1]["unique_kmers"], x[1]["total_kmers"], x[1]["genome_length"], x[1]["order"]))

    def _apply_greedy_filter(self,
                               sorted_genomes: List[Tuple[str, Dict[str, Union[int, float]]]],
                               genome_to_kmers: Dict[str, Set[str]],
                               similarity_threshold: float) -> Tuple[Set[str], Dict[str, Dict[str, Union[str, int, float]]]]:
        """
        @brief Applies a greedy filtering process to decide which genomes to keep.
        @param sorted_genomes Sorted list of (genome_id, stats).
        @param genome_to_kmers Mapping from genome_id to its set of k-mers.
        @param similarity_threshold Threshold for filtering.
        @return A tuple (kept_ids, similarity_info).
        """
        kept: List[Tuple[str, Set[str]]] = []
        similarity_info: Dict[str, Dict[str, Union[str, int, float]]] = {}
        for genome_id, stats in sorted_genomes:
            current_kmers = genome_to_kmers.get(genome_id, set())
            filtered_flag = False
            for kept_genome_id, kept_kmers in kept:
                min_count = min(len(current_kmers), len(kept_kmers))
                sim_score = (len(current_kmers.intersection(kept_kmers)) / min_count) if min_count > 0 else 0
                if sim_score > similarity_threshold:
                    similarity_info[genome_id] = {
                        "kept": "no",
                        "unique_kmers": stats["unique_kmers"],
                        "total_kmers": stats["total_kmers"],
                        "genome_length": stats["genome_length"],
                        "similar_to": kept_genome_id,
                        "similarity_score": round(sim_score, 2),
                    }
                    filtered_flag = True
                    break
            if not filtered_flag:
                similarity_info[genome_id] = {
                    "kept": "yes",
                    "unique_kmers": stats["unique_kmers"],
                    "total_kmers": stats["total_kmers"],
                    "genome_length": stats["genome_length"],
                    "similar_to": "NA",
                    "similarity_score": "NA",
                }
                kept.append((genome_id, current_kmers))
        kept_ids = {genome_id for genome_id, info in similarity_info.items() if info["kept"] == "yes"}
        return kept_ids, similarity_info

    def _remove_filtered_genomes_from_kmers(self, kept_ids: Set[str]) -> None:
        """
        @brief Removes k-mer mappings for genomes not in the kept_ids set.
        @param kept_ids Set of genome IDs that are kept.
        """
        for kmer in list(self.kmers.keys()):
            mapping = self.kmers[kmer]
            for genome_record in list(mapping.keys()):
                if genome_record.identifier not in kept_ids:
                    del mapping[genome_record]
            if not mapping:
                del self.kmers[kmer]

    def _update_genomes_list(self, kept_ids: Set[str]) -> None:
        """
        @brief Updates the list of genomes to include only those in kept_ids.
        @param kept_ids Set of genome IDs that are kept.
        """
        self.genomes = [genome for genome in self.genomes if genome.identifier in kept_ids]

    def _filter_similar_genomes(self, similarity_threshold: float) -> None:
        """
        @brief Applies the overall filtering process to remove highly similar genomes.
        @param similarity_threshold Threshold for filtering.
        @details Updates self.kmers, self.genomes, and stores filtering details in self.similarity_info.
        """
        genome_stats, genome_to_kmers = self._compute_genome_stats()
        sorted_genomes = self._sort_genomes_for_filtering(genome_stats)
        kept_ids, similarity_info = self._apply_greedy_filter(sorted_genomes, genome_to_kmers, similarity_threshold)
        self._remove_filtered_genomes_from_kmers(kept_ids)
        self._update_genomes_list(kept_ids)
        self.similarity_info = similarity_info

    def save(self, ref_file: str) -> None:
        """
        @brief Saves the KmerReference to a gzipped pickle file.
        @param ref_file Output file path.
        """
        with gzip.open(ref_file, "wb") as f:
            pickle.dump(self, f)

    @classmethod
    def load(cls, ref_file: str) -> "KmerReference":
        """
        @brief Loads a KmerReference from a gzipped pickle file.
        @param ref_file Input file path.
        @return The loaded KmerReference instance.
        """
        with gzip.open(ref_file, "rb") as f:
            ref = pickle.load(f)
        return ref

    def __getitem__(self, kmer: str) -> Optional[Dict[Record, Set[int]]]:
        """
        @brief Gets the k-mer mapping for a given k-mer.
        @param kmer The k-mer string.
        @return Mapping of genome records to positions or None if k-mer not found.
        """
        return self.kmers.get(kmer, None)

    def get_kmer_references(self, kmer: str) -> Dict[Record, Set[int]]:
        """
        @brief Returns the genome references for a given k-mer.
        @param kmer The k-mer string.
        @return A dictionary mapping genome records to a set of positions.
        """
        return self.kmers.get(kmer, {})

    def get_summary(self) -> Dict[str, Dict[str, Union[Dict[str, List[int]], Dict[str, int], Dict[str, Dict]]]]:
        """
        @brief Returns a summary of the k-mer reference database.
        @return A dictionary containing k-mer details, genome summary, and filtering similarity info if available.
        """
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
            unique_kmers = sum(1 for kmer in genome_kmer_count[genome] if len(self.kmers[kmer]) == 1)
            multi_mapping_kmers = len(genome_kmer_count[genome]) - unique_kmers
            genome_summary[genome]["unique_kmers"] = unique_kmers
            genome_summary[genome]["multi_mapping_kmers"] = multi_mapping_kmers
        summary_dict = {"Kmers": kmer_details, "Summary": dict(genome_summary)}
        if hasattr(self, "similarity_info"):
            summary_dict["Similarity"] = self.similarity_info
        return summary_dict

    def get_kmer_and_reverse_references(self, kmer: str) -> Dict[Record, Set[int]]:
        """
        @brief Returns genome references for a k-mer and its reverse complement.
        @param kmer The k-mer string.
        @return A dictionary mapping genome records to positions (merged for direct and reverse complement).
        """
        result: Dict[Record, Set[int]] = {}
        direct_refs = self.get_kmer_references(kmer)
        if direct_refs:
            for genome, positions in direct_refs.items():
                result[genome] = set(positions)
        rev_kmer = reverse_complement(kmer)
        if rev_kmer != kmer:
            rev_refs = self.get_kmer_references(rev_kmer)
            if rev_refs:
                for genome, positions in rev_refs.items():
                    if genome in result:
                        result[genome].update(positions)
                    else:
                        result[genome] = set(positions)
        return result

##===========================================================
## Class: Read
##===========================================================

class Read:
    """
    @brief Represents a sequencing read.
    """
    def __init__(self, fastaq_record: Record) -> None:
        """
        @brief Initializes a Read from a FASTAQ record.
        @param fastaq_record A record containing identifier, sequence, and quality.
        """
        self.identifier = fastaq_record.identifier
        self.mapping = ReadMapping(ReadMappingType.UNMAPPED, [])
        self.kmers: Dict[str, ReadKmer] = {}
        self.__raw_read: str = fastaq_record["sequence"]
        self.__quality_scores: str = fastaq_record["quality_sequence"]
        self.num_quality_filtered_kmers: int = 0
        self.num_redundant_kmers: int = 0
        self.__genomes_map_count: Optional[Dict[Record, int]] = None

    def __str__(self) -> str:
        """
        @brief Returns a string representation of the read.
        """
        rows = [f"Mapping: {self.mapping}"]
        for kmer in self.kmers.keys():
            rows.append(f"k-mer: {kmer}")
            rows.append(f"specifity: {self.kmers[kmer].specifity}")
            rows.append("Genome References:")
            for references in self.kmers[kmer]:
                rows.append(f"\t{references}")
        return "\n".join(rows)

    def __repr__(self) -> str:
        """
        @brief Returns a representation of the read.
        """
        return self.__str__()

    def mean_quality(self) -> float:
        """
        @brief Calculates the mean quality of the read.
        @return Mean quality score.
        """
        return sum(map(ord, self.__quality_scores)) / len(self.__quality_scores)

    def kmer_quality(self, start: int, k: int) -> float:
        """
        @brief Calculates the mean quality of a k-mer starting at a given position.
        @param start Starting index of the k-mer.
        @param k Length of the k-mer.
        @return Mean quality score of the k-mer.
        """
        return sum(map(ord, self.__quality_scores[start: start + k])) / k

    def extract_kmer_references(self, kmer_reference: KmerReference,
                                min_kmer_quality: Optional[int] = None,
                                max_genomes: Optional[int] = None) -> None:
        """
        @brief Extracts k-mer references from the read using a KmerReference.
        @param kmer_reference The KmerReference instance.
        @param min_kmer_quality Optional minimum quality threshold for a k-mer.
        @param max_genomes Optional maximum allowed genome mappings for a k-mer.
        """
        for start, kmer in extract_kmers_from_genome(kmer_reference.kmer_len, self.__raw_read):
            if min_kmer_quality is not None and self.kmer_quality(start, kmer_reference.kmer_len) < min_kmer_quality:
                self.num_quality_filtered_kmers += 1
                continue
            single_kmer_reference = kmer_reference.get_kmer_references(kmer)
            if single_kmer_reference:
                if max_genomes is not None and len(single_kmer_reference) > max_genomes:
                    self.num_redundant_kmers += 1
                    continue
                kmer_specifity = (KmerSpecifity.SPECIFIC if len(single_kmer_reference) == 1
                                  else KmerSpecifity.UNSPECIFIC)
                self.kmers[kmer] = ReadKmer(specifity=kmer_specifity, references=single_kmer_reference)

    def generate_genome_counts(self, map_count: bool = False) -> Dict[Record, int]:
        """
        @brief Generates counts of genome mappings from the read's k-mers.
        @param map_count If True, only counts specific mappings.
        @return Dictionary mapping genome records to counts.
        """
        genomes_and_count: Dict[Record, int] = {}
        for kmer in self.kmers.values():
            if not map_count or kmer.specifity == KmerSpecifity.SPECIFIC:
                for genome in kmer.references.keys():
                    genomes_and_count[genome] = genomes_and_count.get(genome, 0) + 1
        return genomes_and_count

    def try_to_align_specific(self, m: int) -> bool:
        """
        @brief Attempts to align the read using its specific k-mer counts.
        @param m Minimum difference required for unique mapping.
        @return True if uniquely mapped, False otherwise.
        """
        if m < 0:
            raise ValueError("m must be non-negative")
        self.__genomes_map_count = self.generate_genome_counts(map_count=True)
        if len(self.__genomes_map_count) == 1:
            self.mapping = ReadMapping(ReadMappingType.UNIQUELY_MAPPED, [list(self.__genomes_map_count.keys())[0]])
            return True
        if len(self.__genomes_map_count) > 1:
            sorted_genomes = sorted(self.__genomes_map_count, key=self.__genomes_map_count.get, reverse=True)
            if self.__genomes_map_count[sorted_genomes[0]] >= self.__genomes_map_count[sorted_genomes[1]] + m:
                self.mapping = ReadMapping(ReadMappingType.UNIQUELY_MAPPED, [sorted_genomes[0]])
                return True
        self.mapping = ReadMapping(ReadMappingType.AMBIGUOUSLY_MAPPED, list(self.__genomes_map_count.keys()))
        return False

    def validate_unique_mappings(self, p: int) -> None:
        """
        @brief Validates the unique mapping by comparing total k-mer counts.
        @param p Allowable difference threshold.
        """
        if self.mapping.type != ReadMappingType.UNIQUELY_MAPPED or p < IGNORE_AMBIGUOUS_THRESHOLD:
            return
        genome_total_counts = self.generate_genome_counts(map_count=False)
        mapped_genome = self.mapping.genomes_mapped_to[0]
        max_total_kmers = max(genome_total_counts.values(), default=0)
        mapped_genome_kmers = genome_total_counts.get(mapped_genome, 0)
        if max_total_kmers - mapped_genome_kmers > p:
            ambiguous_genomes = [mapped_genome]
            for genome, count in genome_total_counts.items():
                if count >= mapped_genome_kmers:
                    ambiguous_genomes.append(genome)
            self.mapping = ReadMapping(ReadMappingType.AMBIGUOUSLY_MAPPED, ambiguous_genomes)

    def pseudo_align(self,
                     kmer_reference: KmerReference,
                     p: int = 1,
                     m: int = 1,
                     min_read_quality: Optional[int] = None,
                     min_kmer_quality: Optional[int] = None,
                     max_genomes: Optional[int] = None,
                     debug: bool = False) -> ReadMappingType:
        """
        @brief Performs pseudo-alignment applying quality and redundancy filters.
        @param kmer_reference The KmerReference instance.
        @param p Ambiguity threshold.
        @param m Unique mapping threshold.
        @param min_read_quality Optional minimum read quality threshold.
        @param min_kmer_quality Optional minimum k-mer quality threshold.
        @param max_genomes Optional maximum allowed genome mappings for a k-mer.
        @param debug If True, prints debug information.
        @return The resulting ReadMappingType.
        """
        if not (isinstance(kmer_reference, KmerReference)
                and isinstance(p, int)
                and isinstance(m, int)
                and (min_read_quality is None or isinstance(min_read_quality, int))
                and (min_kmer_quality is None or isinstance(min_kmer_quality, int))
                and (max_genomes is None or isinstance(max_genomes, int))
                and isinstance(debug, bool)):
            raise TypeError(f"Invalid types given to pseudo align: {type(kmer_reference)}, {type(p)}, {type(m)}, {type(debug)}")
        if m < M_THRESHOLD:
            raise ValueError(f"m must be bigger than or equal to {M_THRESHOLD}")

        if min_read_quality is not None and self.mean_quality() < min_read_quality:
            return ReadMappingType.UNMAPPED

        self.extract_kmer_references(kmer_reference, min_kmer_quality, max_genomes)
        if not self.kmers:
            return ReadMappingType.UNMAPPED

        if self.try_to_align_specific(m):
            if debug:
                print(f"[DEBUG pseudo_align]: After try_to_align_specific self.mapping: {self.mapping.type}")
            self.validate_unique_mappings(p)
            return self.mapping.type
        elif debug:
            print(f"[DEBUG pseudo_align]: After try_to_align_specific self.mapping: {self.mapping.type}, mapped to: {self.mapping}")
        return ReadMappingType.AMBIGUOUSLY_MAPPED

##===========================================================
## Class: PseudoAlignment
##===========================================================

class PseudoAlignment:
    """
    @brief Manages pseudo-alignment of reads using a k-mer reference.
    """
    def __init__(self, kmer_reference: KmerReference) -> None:
        """
        @brief Initializes a PseudoAlignment with a given k-mer reference.
        @param kmer_reference The KmerReference instance.
        """
        self.kmer_reference: KmerReference = kmer_reference
        self.reads: Dict[str, Dict[str, Union[ReadMappingType, List[str]]]] = {}
        self.filtered_quality_reads: int = 0
        self.filtered_quality_kmers: int = 0
        self.filtered_hr_kmers: int = 0

        self.filter_read_quality_flag: bool = False
        self.filter_kmer_quality_flag: bool = False
        self.filter_max_genomes_flag: bool = False

    def add_read(self, read: Read) -> None:
        """
        @brief Adds an aligned read to the alignment.
        @param read The Read object to add.
        """
        if read.identifier in self.reads:
            raise AddingExistingRead(f"There already exists a read with identifier: {read.identifier}")
        self.reads[read.identifier] = {
            "mapping_type": read.mapping.type,
            "genomes_mapped_to": [genome.identifier for genome in read.mapping.genomes_mapped_to],
        }

    def add_read_from_read_record(self,
                                  read_record: Record,
                                  m: int = 1,
                                  p: int = 1,
                                  min_read_quality: Optional[int] = None,
                                  min_kmer_quality: Optional[int] = None,
                                  max_genomes: Optional[int] = None) -> None:
        """
        @brief Aligns a read from a FASTAQ record while applying filters.
        @param read_record The FASTAQ record.
        @param m Unique mapping threshold.
        @param p Ambiguity threshold.
        @param min_read_quality Optional minimum read quality threshold.
        @param min_kmer_quality Optional minimum k-mer quality threshold.
        @param max_genomes Optional maximum allowed genome mappings for a k-mer.
        """
        if min_read_quality is not None:
            self.filter_read_quality_flag = True
        if min_kmer_quality is not None:
            self.filter_kmer_quality_flag = True
        if max_genomes is not None:
            self.filter_max_genomes_flag = True

        read = Read(read_record)
        if min_read_quality is not None and read.mean_quality() < min_read_quality:
            self.filtered_quality_reads += 1
            return
        read.pseudo_align(self.kmer_reference, m=m, p=p,
                          min_read_quality=min_read_quality,
                          min_kmer_quality=min_kmer_quality,
                          max_genomes=max_genomes)
        if min_kmer_quality is not None:
            self.filtered_quality_kmers += read.num_quality_filtered_kmers
        if max_genomes is not None:
            self.filtered_hr_kmers += read.num_redundant_kmers
        self.add_read(read)

    def align_reads_from_container(self,
                                   reads_container: FASTAQRecordContainer,
                                   m: int = 1,
                                   p: int = 1,
                                   min_read_quality: Optional[int] = None,
                                   min_kmer_quality: Optional[int] = None,
                                   max_genomes: Optional[int] = None) -> None:
        """
        @brief Aligns all reads from a FASTAQ container.
        @param reads_container The FASTAQRecordContainer.
        @param m Unique mapping threshold.
        @param p Ambiguity threshold.
        @param min_read_quality Optional minimum read quality threshold.
        @param min_kmer_quality Optional minimum k-mer quality threshold.
        @param max_genomes Optional maximum allowed genome mappings for a k-mer.
        """
        for read_record in reads_container:
            self.add_read_from_read_record(read_record, m=m, p=p,
                                           min_read_quality=min_read_quality,
                                           min_kmer_quality=min_kmer_quality,
                                           max_genomes=max_genomes)

    def get_summary(self) -> Dict[str, Dict[str, Union[int, Dict[str, int]]]]:
        """
        @brief Returns a summary of the alignment.
        @return Dictionary with alignment statistics and genome mapping summary.
        """
        summary: Dict[str, Union[int, Dict[str, int]]] = {
            "unique_mapped_reads": 0,
            "ambiguous_mapped_reads": 0,
            "unmapped_reads": 0,
        }
        if self.filter_read_quality_flag:
            summary["filtered_quality_reads"] = self.filtered_quality_reads
        if self.filter_kmer_quality_flag:
            summary["filtered_quality_kmers"] = self.filtered_quality_kmers
        if self.filter_max_genomes_flag:
            summary["filtered_hr_kmers"] = self.filtered_hr_kmers

        genome_mapping: Dict[str, Dict[str, int]] = {}
        for _, read_alignment_details in self.reads.items():
            mapping_type = read_alignment_details["mapping_type"]
            genomes_mapped_to = read_alignment_details["genomes_mapped_to"]
            if mapping_type != ReadMappingType.UNMAPPED:
                read_type_string: Optional[str] = None
                if mapping_type == ReadMappingType.UNIQUELY_MAPPED:
                    summary["unique_mapped_reads"] += 1
                    read_type_string = "unique_reads"
                elif mapping_type == ReadMappingType.AMBIGUOUSLY_MAPPED:
                    summary["ambiguous_mapped_reads"] += 1
                    read_type_string = "ambiguous_reads"
                for genome in genomes_mapped_to:
                    genome_mapping.setdefault(genome, {"unique_reads": 0, "ambiguous_reads": 0})
                    genome_mapping[genome][read_type_string] += 1
            else:
                summary["unmapped_reads"] += 1
        return {"Statistics": summary, "Summary": genome_mapping}

    def save(self, align_file: str) -> None:
        """
        @brief Saves the alignment to a gzipped pickle file.
        @param align_file Output file path.
        """
        with gzip.open(align_file, "wb") as f:
            pickle.dump(self, f)

    def __repr__(self) -> str:
        """
        @brief Returns a JSON string representation of the alignment summary.
        """
        return json.dumps(self.get_summary(), indent=4)

    @classmethod
    def load(cls, align_file: str) -> "PseudoAlignment":
        """
        @brief Loads a PseudoAlignment object from a gzipped pickle file.
        @param align_file Input file path.
        @return The loaded PseudoAlignment instance.
        """
        with gzip.open(align_file, "rb") as f:
            alignment = pickle.load(f)
        return alignment

    def export_summary_to_json(self, json_file: str) -> None:
        """
        @brief Exports the alignment summary to a JSON file.
        @param json_file Output JSON file path.
        """
        summary = self.get_summary()
        with open(json_file, "w") as f:
            json.dump(summary, f, indent=4)

    def get_reads_by_mapping_type(self, mapping_type: ReadMappingType) -> List[str]:
        """
        @brief Retrieves a list of read identifiers with the specified mapping type.
        @param mapping_type The desired ReadMappingType.
        @return List of read identifiers.
        """
        return [read_id for read_id, details in self.reads.items() if details["mapping_type"] == mapping_type]
