#!/usr/bin/env python3
"""
@file main.py
@brief Command-line interface for the Biosequence project.
@details Supports building a k-mer reference database, performing pseudo-alignment,
         and dumping alignment/reference summaries.
"""

import os
import argparse
import pickle
import gzip
import sys
import json
from typing import Optional, List

from data_file import FASTAFile, FASTAQFile, InvalidExtensionError, NoRecordsInDataFile
from kmer import (
    KmerReference,
    PseudoAlignment,
    NotValidatingUniqueMapping,
    AddingExistingRead,
)


## ===========================================================
## File Validation Functions
## ===========================================================

## @brief Checks whether a file exists and is readable.
#  @param filepath Path to the file.
#  @param description Description of the file (used for error messages).
#  @return None. Exits the program if the file is not readable.
def validate_file_readable(filepath: str, description: str) -> None:
    if not os.path.isfile(filepath):
        sys.exit(f"Error: {description} file '{filepath}' does not exist or is not a file.")
    if not os.access(filepath, os.R_OK):
        sys.exit(f"Error: {description} file '{filepath}' is not readable.")


## @brief Checks whether a file (or its directory) is writable.
#  @param filepath Path to the file.
#  @param description Description of the file (used for error messages).
#  @return None. Exits the program if the file is not writable.
def validate_file_writable(filepath: str, description: str) -> None:
    dir_path: str = os.path.dirname(filepath) or "."
    if os.path.exists(filepath) and not os.access(filepath, os.W_OK):
        sys.exit(f"Error: {description} file '{filepath}' is not writable.")
    if not os.path.exists(filepath) and not os.access(dir_path, os.W_OK):
        sys.exit(f"Error: Directory '{dir_path}' is not writable to create {description} file '{filepath}'.")


## ===========================================================
## Argument Parsing
## ===========================================================

## @brief Parses command-line arguments.
#  @param args Optional list of arguments (for testing).
#  @return An argparse.Namespace object containing parsed arguments.
def parse_arguments(args: Optional[List[str]] = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(prog="Biosequence project")
    parser.add_argument("-t", "--task", required=True, help="Task to execute")
    parser.add_argument("-g", "--genomefile", help="Genome FASTA file (multiple records)")
    parser.add_argument("-k", "--kmer-size", type=int, help="Length of k-mers")
    parser.add_argument("-r", "--referencefile", help="KDB file (input/output)")
    parser.add_argument("-a", "--alignfile", help="aln file. Can be either input or name for output file")
    parser.add_argument("--reads", help="FASTQ reads file")
    parser.add_argument("-m", "--unique-threshold", help="unique k-mer threshold", default=1, type=int)
    parser.add_argument("-p", "--ambiguous-threshold", help="ambiguous k-mer threshold", default=1, type=int)
    parser.add_argument("--reverse-complement", action="store_true")
    parser.add_argument("--min-read-quality", type=int, default=None)
    parser.add_argument("--min-kmer-quality", type=int, default=None)
    parser.add_argument("--filter-similar", action="store_true")
    parser.add_argument("--similarity-threshold", type=float, default=0.95)
    return parser.parse_args(args)


## ===========================================================
## Reference Creation Functions
## ===========================================================

## @brief Creates a k-mer reference from a FASTA file.
#  @param fasta_file Path to the FASTA file.
#  @param kmer_size Length of the k-mers.
#  @param filter_similar Flag to enable filtering of similar genomes.
#  @param similarity_threshold Threshold for genome similarity.
#  @return A KmerReference instance.
def create_reference(fasta_file: str, kmer_size: int,
                     filter_similar: bool = False, similarity_threshold: float = 0.95) -> KmerReference:
    fasta_container = FASTAFile(fasta_file).container
    kmer_reference = KmerReference(kmer_size, fasta_container,
                                   filter_similar=filter_similar,
                                   similarity_threshold=similarity_threshold)
    return kmer_reference


## @brief Creates a k-mer reference and saves it to a file.
#  @param fasta_file Path to the FASTA file.
#  @param kmer_size Length of the k-mers.
#  @param reference_file Output file path.
#  @param filter_similar Flag to enable filtering of similar genomes.
#  @param similarity_threshold Threshold for genome similarity.
#  @return None.
def create_reference_and_save_it(fasta_file: str, kmer_size: int, reference_file: str,
                                   filter_similar: bool = False, similarity_threshold: float = 0.95) -> None:
    kmer_reference = create_reference(fasta_file, kmer_size, filter_similar, similarity_threshold)
    with gzip.open(reference_file, "wb") as f:
        pickle.dump(kmer_reference, f)


## @brief Dumps the k-mer reference summary to stdout.
#  @param kmer_reference A KmerReference instance.
#  @return None.
def dump_reference(kmer_reference: KmerReference) -> None:
    print(json.dumps(kmer_reference.get_summary(), indent=4))


## @brief Loads a k-mer reference from a file and dumps its summary.
#  @param reference_file Input file path.
#  @return None.
def dump_reference_file(reference_file: str) -> None:
    with gzip.open(reference_file, "rb") as f:
        kmer_reference = pickle.load(f)
    dump_reference(kmer_reference)


## @brief Builds a k-mer reference from a FASTA file and dumps its summary.
#  @param fasta_file Path to the FASTA file.
#  @param kmer_size Length of the k-mers.
#  @param filter_similar Flag to enable filtering of similar genomes.
#  @param similarity_threshold Threshold for genome similarity.
#  @return None.
def build_reference_and_dump_from_file(fasta_file: str, kmer_size: int,
                                         filter_similar: bool = False, similarity_threshold: float = 0.95) -> None:
    kmer_reference = create_reference(fasta_file, kmer_size, filter_similar, similarity_threshold)
    dump_reference(kmer_reference)


## ===========================================================
## Alignment Creation Functions
## ===========================================================

## @brief Creates a pseudo-alignment from a k-mer reference and a FASTQ file.
#  @param kmer_reference A KmerReference instance.
#  @param reads_file Path to the FASTQ file.
#  @param p Ambiguity threshold.
#  @param m Unique mapping threshold.
#  @param min_read_quality Minimum read quality threshold.
#  @param min_kmer_quality Minimum k-mer quality threshold.
#  @param max_genomes Maximum allowed genome mappings for a k-mer.
#  @return A PseudoAlignment instance.
def create_alignment_from_reference(kmer_reference: KmerReference, reads_file: str,
                                      p: int, m: int, min_read_quality: int,
                                      min_kmer_quality: int, max_genomes: int) -> PseudoAlignment:
    reads_container = FASTAQFile(reads_file).container
    pseudo_alignment = PseudoAlignment(kmer_reference)
    pseudo_alignment.align_reads_from_container(reads_container, p, m,
                                                min_read_quality, min_kmer_quality, max_genomes)
    return pseudo_alignment


## @brief Creates a pseudo-alignment from a k-mer reference and saves it to a file.
#  @param kmer_reference A KmerReference instance.
#  @param reads_file Path to the FASTQ file.
#  @param align_file Output file path.
#  @param p Ambiguity threshold.
#  @param m Unique mapping threshold.
#  @param min_read_quality Minimum read quality threshold.
#  @param min_kmer_quality Minimum k-mer quality threshold.
#  @param max_genomes Maximum allowed genome mappings for a k-mer.
#  @return None.
def create_alignment_file_from_reference(kmer_reference: KmerReference, reads_file: str,
                                           align_file: str, p: int, m: int,
                                           min_read_quality: int, min_kmer_quality: int,
                                           max_genomes: int) -> None:
    pseudo_alignment = create_alignment_from_reference(kmer_reference, reads_file, p, m,
                                                       min_read_quality, min_kmer_quality, max_genomes)
    pseudo_alignment.save(align_file)


## @brief Loads a k-mer reference from file and creates an alignment file.
#  @param reference_file Path to the k-mer reference file.
#  @param reads_file Path to the FASTQ reads file.
#  @param align_file Output alignment file path.
#  @param p Ambiguity threshold.
#  @param m Unique mapping threshold.
#  @param min_read_quality Minimum read quality threshold.
#  @param min_kmer_quality Minimum k-mer quality threshold.
#  @param max_genomes Maximum allowed genome mappings for a k-mer.
#  @return None.
def create_alignment_from_reference_file(reference_file: str, reads_file: str,
                                           align_file: str, p: int, m: int,
                                           min_read_quality: int, min_kmer_quality: int,
                                           max_genomes: int) -> None:
    with gzip.open(reference_file, "rb") as f:
        kmer_reference = pickle.load(f)
    create_alignment_file_from_reference(kmer_reference, reads_file, align_file, p, m,
                                           min_read_quality, min_kmer_quality, max_genomes)


## @brief Builds a k-mer reference from a FASTA file and creates an alignment file.
#  @param fasta_file Path to the FASTA file.
#  @param kmer_size Length of the k-mers.
#  @param reads_file Path to the FASTQ reads file.
#  @param align_file Output alignment file path.
#  @param m Unique mapping threshold.
#  @param p Ambiguity threshold.
#  @param min_read_quality Minimum read quality threshold.
#  @param min_kmer_quality Minimum k-mer quality threshold.
#  @param max_genomes Maximum allowed genome mappings for a k-mer.
#  @param filter_similar Flag to enable filtering of similar genomes.
#  @param similarity_threshold Threshold for genome similarity.
#  @return None.
def build_reference_and_create_alignment_file(fasta_file: str, kmer_size: int,
                                                reads_file: str, align_file: str,
                                                m: int, p: int, min_read_quality: int,
                                                min_kmer_quality: int, max_genomes: int,
                                                filter_similar: bool = False,
                                                similarity_threshold: float = 0.95) -> None:
    kmer_reference = create_reference(fasta_file, kmer_size, filter_similar, similarity_threshold)
    create_alignment_file_from_reference(kmer_reference, reads_file, align_file, m, p,
                                           min_read_quality, min_kmer_quality, max_genomes)


## @brief Loads an alignment file and dumps its summary.
#  @param align_file Path to the alignment file.
#  @return None.
def dump_alignment_file(align_file: str) -> None:
    with gzip.open(align_file, "rb") as f:
        pseudo_alignment = pickle.load(f)
    print(json.dumps(pseudo_alignment.get_summary(), indent=4))


## @brief Loads a k-mer reference from file, creates an alignment, and dumps its summary.
#  @param reference_file Path to the k-mer reference file.
#  @param reads_file Path to the FASTQ reads file.
#  @param m Ambiguity threshold.
#  @param p Unique mapping threshold.
#  @param min_read_quality Minimum read quality threshold.
#  @param min_kmer_quality Minimum k-mer quality threshold.
#  @param max_genomes Maximum allowed genome mappings for a k-mer.
#  @return None.
def dump_alignment_from_reference(reference_file: str, reads_file: str,
                                    m: int, p: int, min_read_quality: int,
                                    min_kmer_quality: int, max_genomes: int) -> None:
    with gzip.open(reference_file, "rb") as f:
        kmer_reference = pickle.load(f)
    pseudo_alignment = create_alignment_from_reference(kmer_reference, reads_file, p, m,
                                                       min_read_quality, min_kmer_quality, max_genomes)
    print(json.dumps(pseudo_alignment.get_summary(), indent=4))


## @brief Builds a k-mer reference, creates an alignment, and dumps its summary.
#  @param fasta_file Path to the FASTA file.
#  @param kmer_size Length of the k-mers.
#  @param reads_file Path to the FASTQ reads file.
#  @param m Unique mapping threshold.
#  @param p Ambiguity threshold.
#  @param min_read_quality Minimum read quality threshold.
#  @param min_kmer_quality Minimum k-mer quality threshold.
#  @param max_genomes Maximum allowed genome mappings for a k-mer.
#  @param filter_similar Flag to enable filtering of similar genomes.
#  @param similarity_threshold Threshold for genome similarity.
#  @return None.
def build_reference_align_and_dump(fasta_file: str, kmer_size: int, reads_file: str,
                                     m: int, p: int, min_read_quality: int,
                                     min_kmer_quality: int, max_genomes: int,
                                     filter_similar: bool = False, similarity_threshold: float = 0.95) -> None:
    kmer_reference = create_reference(fasta_file, kmer_size, filter_similar, similarity_threshold)
    pseudo_alignment = create_alignment_from_reference(kmer_reference, reads_file, p, m,
                                                       min_read_quality, min_kmer_quality, max_genomes)
    print(json.dumps(pseudo_alignment.get_summary(), indent=4))


## ===========================================================
## Main Entry Point
## ===========================================================

## @brief Main function that dispatches tasks based on command-line arguments.
#  @return None.
def main() -> None:
    args = parse_arguments()

    try:
        if args.task == "reference":
            if not args.genomefile or not args.kmer_size or not args.referencefile:
                sys.exit("Error: -g (genome file), -k (kmer size), and -r (reference file) are required for this task.")
            validate_file_readable(args.genomefile, "Genome FASTA")
            validate_file_writable(args.referencefile, "Reference database output")
            create_reference_and_save_it(args.genomefile, args.kmer_size, args.referencefile,
                                           args.filter_similar, args.similarity_threshold)
        elif args.task == "dumpref":
            if args.referencefile:
                validate_file_readable(args.referencefile, "Reference database")
                dump_reference_file(args.referencefile)
            elif args.genomefile and args.kmer_size:
                validate_file_readable(args.genomefile, "Genome FASTA")
                build_reference_and_dump_from_file(args.genomefile, args.kmer_size,
                                                     args.filter_similar, args.similarity_threshold)
            else:
                sys.exit("Error: Either -r (reference file) or -g (genome file) and -k (kmer size) are required for this task.")
        elif args.task == "align":
            validate_file_readable(args.reads, "FASTQ reads")
            validate_file_writable(args.alignfile, "Alignment output")
            if args.referencefile and args.reads and args.alignfile:
                validate_file_readable(args.referencefile, "Reference database")
                create_alignment_from_reference_file(args.referencefile, args.reads, args.alignfile,
                                                       args.unique_threshold, args.ambiguous_threshold,
                                                       args.min_read_quality, args.min_kmer_quality,
                                                       args.max_genomes)
            elif args.genomefile and args.kmer_size and args.reads and args.alignfile:
                validate_file_readable(args.genomefile, "Genome FASTA")
                build_reference_and_create_alignment_file(args.genomefile, args.kmer_size, args.reads,
                                                            args.alignfile, args.unique_threshold,
                                                            args.ambiguous_threshold, args.min_read_quality,
                                                            args.min_kmer_quality, args.max_genomes,
                                                            args.filter_similar, args.similarity_threshold)
            else:
                sys.exit("Error: Either provide -r (reference file) or -g (genome file) and -k (kmer size) along with --reads and -a (alignment output file).")
        elif args.task == "dumpalign":
            if args.referencefile and args.reads:
                validate_file_readable(args.reads, "FASTQ reads")
                dump_alignment_from_reference(args.referencefile, args.reads, args.unique_threshold,
                                              args.ambiguous_threshold, args.min_read_quality,
                                              args.min_kmer_quality, args.max_genomes)
            if args.genomefile and args.kmer_size and args.reads:
                validate_file_readable(args.reads, "FASTQ reads")
                validate_file_readable(args.genomefile, "Genome FASTA")
                build_reference_align_and_dump(args.genomefile, args.kmer_size, args.reads,
                                               args.unique_threshold, args.ambiguous_threshold,
                                               args.min_read_quality, args.min_kmer_quality, args.max_genomes,
                                               args.filter_similar, args.similarity_threshold)
            elif args.alignfile:
                validate_file_writable(args.alignfile, "Alignment output")
                dump_alignment_file(args.alignfile)
            else:
                sys.exit("Error: Provide either -g (genome file) and -k (kmer size) with --reads (FASTQ file), or -a (alignment output file).")
        else:
            sys.exit("Error: Unsupported task.")
    except (InvalidExtensionError, NoRecordsInDataFile, NotValidatingUniqueMapping, AddingExistingRead) as err:
        sys.exit(err)


if __name__ == "__main__":
    main()
