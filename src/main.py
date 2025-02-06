#!/usr/bin/env python3
"""
@file main.py
@brief Command-line interface for the Biosequence project.
@details Supports building a k-mer reference database, performing pseudo-alignment,
         and dumping alignment/reference summaries.
"""

import os
import argparse
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

from constants import DEFAULT_AMBIGUOUS_THRESHOLD, DEFAULT_SIMILARITY_THRESHOLD, DEFAULT_UNIQUE_THRESHOLD

## ===========================================================
## File Validation Functions
## ===========================================================

def validate_file_readable(filepath: str, description: str) -> None:
    """
    @brief Checks whether a file exists and is readable.
    @param filepath Path to the file.
    @param description Description of the file (used for error messages).
    @return None. Exits the program if the file is not readable.
    """
    if not os.path.isfile(filepath):
        sys.exit(f"Error: {description} file '{filepath}' does not exist or is not a file.")
    if not os.access(filepath, os.R_OK):
        sys.exit(f"Error: {description} file '{filepath}' is not readable.")


def validate_file_writable(filepath: str, description: str) -> None:
    """
    @brief Checks whether a file (or its directory) is writable.
    @param filepath Path to the file.
    @param description Description of the file (used for error messages).
    @return None. Exits the program if the file is not writable.
    """
    dir_path: str = os.path.dirname(filepath) or "."
    if os.path.exists(filepath) and not os.access(filepath, os.W_OK):
        sys.exit(f"Error: {description} file '{filepath}' is not writable.")
    if not os.path.exists(filepath) and not os.access(dir_path, os.W_OK):
        sys.exit(f"Error: Directory '{dir_path}' is not writable to create {description} file '{filepath}'.")


## ===========================================================
## Argument Parsing
## ===========================================================

def parse_arguments(args: Optional[List[str]] = None) -> argparse.Namespace:
    """
    @brief Parses command-line arguments.
    @param args Optional list of arguments (for testing).
    @return An argparse.Namespace object containing parsed arguments.
    """
    parser = argparse.ArgumentParser(prog="Biosequence project")
    parser.add_argument("-t", "--task", required=True, help="Task to execute")
    parser.add_argument("-g", "--genomefile", help="Genome FASTA file (multiple records)")
    parser.add_argument("-k", "--kmer-size", type=int, help="Length of k-mers")
    parser.add_argument("-r", "--referencefile", help="KDB file (input/output)")
    parser.add_argument("-a", "--alignfile", help="aln file. Can be either input or name for output file")
    parser.add_argument("--reads", help="FASTQ reads file")
    parser.add_argument("-m", "--unique-threshold", help="unique k-mer threshold", type=int)
    parser.add_argument("-p", "--ambiguous-threshold", help="ambiguous k-mer threshold", type=int)
    parser.add_argument("--reverse-complement", action="store_true")
    parser.add_argument("--min-read-quality", type=int, default=None)
    parser.add_argument("--min-kmer-quality", type=int, default=None)
    parser.add_argument("--max-genomes", type=int, default=None)
    parser.add_argument("--filter-similar", action="store_true")
    parser.add_argument("--similarity-threshold", type=float)
    return parser.parse_args(args)


## ===========================================================
## Reference Creation Functions
## ===========================================================

def create_reference(fasta_file: str, kmer_size: int,
                     filter_similar: bool = False, similarity_threshold: float = 0.95) -> KmerReference:
    """
    @brief Creates a k-mer reference from a FASTA file.
    @param fasta_file Path to the FASTA file.
    @param kmer_size Length of the k-mers.
    @param filter_similar Flag to enable filtering of similar genomes.
    @param similarity_threshold Threshold for genome similarity.
    @return A KmerReference instance.
    """
    fasta_container = FASTAFile(fasta_file).container
    kmer_reference = KmerReference(kmer_size, fasta_container,
                                   filter_similar=filter_similar,
                                   similarity_threshold=similarity_threshold)
    return kmer_reference


def create_reference_and_save_it(fasta_file: str, kmer_size: int, reference_file: str,
                                   filter_similar: bool = False, similarity_threshold: float = 0.95) -> None:
    """
    @brief Creates a k-mer reference and saves it to a file.
    @param fasta_file Path to the FASTA file.
    @param kmer_size Length of the k-mers.
    @param reference_file Output file path.
    @param filter_similar Flag to enable filtering of similar genomes.
    @param similarity_threshold Threshold for genome similarity.
    @return None.
    """
    kmer_reference = create_reference(fasta_file, kmer_size, filter_similar, similarity_threshold)
    kmer_reference.save(reference_file)


def dump_reference(kmer_reference: KmerReference) -> None:
    """
    @brief Dumps the k-mer reference summary to stdout.
    @param kmer_reference A KmerReference instance.
    @return None.
    """
    print(json.dumps(kmer_reference.get_summary(), indent=4))


def dump_reference_file(reference_file: str) -> None:
    """
    @brief Loads a k-mer reference from a file and dumps its summary.
    @param reference_file Input file path.
    @return None.
    """
    try:
        kmer_reference = KmerReference.load(reference_file)
    except gzip.BadGzipFile:
        sys.exit("Error: Incorrect format of input file.")
    dump_reference(kmer_reference)


def build_reference_and_dump_from_file(fasta_file: str, kmer_size: int,
                                         filter_similar: bool = False, similarity_threshold: float = 0.95) -> None:
    """
    @brief Builds a k-mer reference from a FASTA file and dumps its summary.
    @param fasta_file Path to the FASTA file.
    @param kmer_size Length of the k-mers.
    @param filter_similar Flag to enable filtering of similar genomes.
    @param similarity_threshold Threshold for genome similarity.
    @return None.
    """
    kmer_reference = create_reference(fasta_file, kmer_size, filter_similar, similarity_threshold)
    dump_reference(kmer_reference)


## ===========================================================
## Alignment Creation Functions
## ===========================================================

def create_alignment_from_reference(kmer_reference: KmerReference, reads_file: str,
                                      p: int, m: int, min_read_quality: int,
                                      min_kmer_quality: int, max_genomes: int) -> PseudoAlignment:
    """
    @brief Creates a pseudo-alignment from a k-mer reference and a FASTQ file.
    @param kmer_reference A KmerReference instance.
    @param reads_file Path to the FASTQ file.
    @param p Ambiguity threshold.
    @param m Unique mapping threshold.
    @param min_read_quality Minimum read quality threshold.
    @param min_kmer_quality Minimum k-mer quality threshold.
    @param max_genomes Maximum allowed genome mappings for a k-mer.
    @return A PseudoAlignment instance.
    """
    reads_container = FASTAQFile(reads_file).container
    pseudo_alignment = PseudoAlignment(kmer_reference)
    pseudo_alignment.align_reads_from_container(reads_container, p, m,
                                                min_read_quality, min_kmer_quality, max_genomes)
    return pseudo_alignment


def create_alignment_file_from_reference(kmer_reference: KmerReference, reads_file: str,
                                           align_file: str, p: int, m: int,
                                           min_read_quality: int, min_kmer_quality: int,
                                           max_genomes: int) -> None:
    """
    @brief Creates a pseudo-alignment from a k-mer reference and saves it to a file.
    @param kmer_reference A KmerReference instance.
    @param reads_file Path to the FASTQ file.
    @param align_file Output file path.
    @param p Ambiguity threshold.
    @param m Unique mapping threshold.
    @param min_read_quality Minimum read quality threshold.
    @param min_kmer_quality Minimum k-mer quality threshold.
    @param max_genomes Maximum allowed genome mappings for a k-mer.
    @return None.
    """
    pseudo_alignment = create_alignment_from_reference(kmer_reference, reads_file, p, m,
                                                       min_read_quality, min_kmer_quality, max_genomes)
    pseudo_alignment.save(align_file)


def create_alignment_from_reference_file(reference_file: str, reads_file: str,
                                           align_file: str, p: int, m: int,
                                           min_read_quality: int, min_kmer_quality: int,
                                           max_genomes: int) -> None:
    """
    @brief Loads a k-mer reference from file and creates an alignment file.
    @param reference_file Path to the k-mer reference file.
    @param reads_file Path to the FASTQ reads file.
    @param align_file Output alignment file path.
    @param p Ambiguity threshold.
    @param m Unique mapping threshold.
    @param min_read_quality Minimum read quality threshold.
    @param min_kmer_quality Minimum k-mer quality threshold.
    @param max_genomes Maximum allowed genome mappings for a k-mer.
    @return None.
    """
    try:
        kmer_reference = KmerReference.load(reference_file)
    except gzip.BadGzipFile:
        sys.exit("Error: Incorrect format of input file.")
    create_alignment_file_from_reference(kmer_reference, reads_file, align_file, p, m,
                                           min_read_quality, min_kmer_quality, max_genomes)


def build_reference_and_create_alignment_file(fasta_file: str, kmer_size: int,
                                                reads_file: str, align_file: str,
                                                m: int, p: int, min_read_quality: int,
                                                min_kmer_quality: int, max_genomes: int,
                                                filter_similar: bool = False,
                                                similarity_threshold: float = 0.95) -> None:
    """
    @brief Builds a k-mer reference from a FASTA file and creates an alignment file.
    @param fasta_file Path to the FASTA file.
    @param kmer_size Length of the k-mers.
    @param reads_file Path to the FASTQ reads file.
    @param align_file Output alignment file path.
    @param m Unique mapping threshold.
    @param p Ambiguity threshold.
    @param min_read_quality Minimum read quality threshold.
    @param min_kmer_quality Minimum k-mer quality threshold.
    @param max_genomes Maximum allowed genome mappings for a k-mer.
    @param filter_similar Flag to enable filtering of similar genomes.
    @param similarity_threshold Threshold for genome similarity.
    @return None.
    """
    kmer_reference = create_reference(fasta_file, kmer_size, filter_similar, similarity_threshold)
    create_alignment_file_from_reference(kmer_reference, reads_file, align_file, p, m,
                                           min_read_quality, min_kmer_quality, max_genomes)


def dump_alignment_file(align_file: str) -> None:
    """
    @brief Loads an alignment file and dumps its summary.
    @param align_file Path to the alignment file.
    @return None.
    """
    try:
        pseudo_alignment = PseudoAlignment.load(align_file)
    except gzip.BadGzipFile:
        sys.exit("Error: Incorrect format of input file.")
    print(json.dumps(pseudo_alignment.get_summary(), indent=4))


def dump_alignment_from_reference(reference_file: str, reads_file: str,
                                    m: int, p: int, min_read_quality: int,
                                    min_kmer_quality: int, max_genomes: int) -> None:
    """
    @brief Loads a k-mer reference from file, creates an alignment, and dumps its summary.
    @param reference_file Path to the k-mer reference file.
    @param reads_file Path to the FASTQ reads file.
    @param m Ambiguity threshold.
    @param p Unique mapping threshold.
    @param min_read_quality Minimum read quality threshold.
    @param min_kmer_quality Minimum k-mer quality threshold.
    @param max_genomes Maximum allowed genome mappings for a k-mer.
    @return None.
    """
    try:
        kmer_reference = KmerReference.load(reference_file)
    except gzip.BadGzipFile:
        sys.exit("Error: Incorrect format of input file.")
    pseudo_alignment = create_alignment_from_reference(kmer_reference, reads_file, p, m,
                                                       min_read_quality, min_kmer_quality, max_genomes)
    print(json.dumps(pseudo_alignment.get_summary(), indent=4))


def build_reference_align_and_dump(fasta_file: str, kmer_size: int, reads_file: str,
                                     m: int, p: int, min_read_quality: int,
                                     min_kmer_quality: int, max_genomes: int,
                                     filter_similar: bool = False, similarity_threshold: float = 0.95) -> None:
    """
    @brief Builds a k-mer reference, creates an alignment, and dumps its summary.
    @param fasta_file Path to the FASTA file.
    @param kmer_size Length of the k-mers.
    @param reads_file Path to the FASTQ reads file.
    @param m Unique mapping threshold.
    @param p Ambiguity threshold.
    @param min_read_quality Minimum read quality threshold.
    @param min_kmer_quality Minimum k-mer quality threshold.
    @param max_genomes Maximum allowed genome mappings for a k-mer.
    @param filter_similar Flag to enable filtering of similar genomes.
    @param similarity_threshold Threshold for genome similarity.
    @return None.
    """
    kmer_reference = create_reference(fasta_file, kmer_size, filter_similar, similarity_threshold)
    pseudo_alignment = create_alignment_from_reference(kmer_reference, reads_file, p, m,
                                                       min_read_quality, min_kmer_quality, max_genomes)
    print(json.dumps(pseudo_alignment.get_summary(), indent=4))


## ===========================================================
## Main Entry Point
## ===========================================================

def main() -> None:
    args = parse_arguments()

    # Validate flag combinations based on the task.
    if args.task == "reference":
        if args.reads or args.alignfile or args.unique_threshold or args.ambiguous_threshold or args.min_read_quality or args.min_kmer_quality or args.max_genomes:
            sys.exit("Error: For task 'reference', only -g, -k, -r, --filter-similar, and --similarity-threshold are allowed.")
    elif args.task == "dumpref":
        if args.reads or args.alignfile or args.unique_threshold or args.ambiguous_threshold or args.min_read_quality or args.min_kmer_quality or args.max_genomes:
            sys.exit ((args.reads, args.alignfile, args.unique_threshold, args.ambiguous_threshold, args.min_read_quality, args.min_kmer_quality, args.max_genomes))
            sys.exit("Error: For task 'dumpref', only -r or (-g and -k) with --filter-similar and --similarity-threshold are allowed.")
    elif args.task == "align":
        if not ((args.referencefile and args.reads and args.alignfile) or (args.genomefile and args.kmer_size and args.reads and args.alignfile)):
            sys.exit("Error: For task 'align', provide either -r (reference file) or -g and -k (genome file and kmer size) along with --reads and -a.")
    elif args.task == "dumpalign":
        if not ((args.referencefile and args.reads) or (args.genomefile and args.kmer_size and args.reads) or args.alignfile):
            sys.exit("Error: For task 'dumpalign', provide either -r and --reads, or -g, -k, and --reads, or -a.")
    else:
        sys.exit("Error: Unsupported task.")

    # Default values for flags that are not None set here for the previous section to work properly
    args.unique_threshold = DEFAULT_UNIQUE_THRESHOLD
    args.ambiguous_threshold = DEFAULT_AMBIGUOUS_THRESHOLD
    args.similarity_threshold =  DEFAULT_SIMILARITY_THRESHOLD
    try:
        if args.task == "reference":
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
                # Use load and save methods from the classes.
                from data_file import FASTAFile
                kmer_ref = create_reference(args.genomefile, args.kmer_size,
                                            args.filter_similar, args.similarity_threshold)
                kmer_ref.save(args.referencefile)  # Save the reference using its method.
                create_alignment_from_reference_file(args.referencefile, args.reads, args.alignfile,
                                                       args.unique_threshold, args.ambiguous_threshold,
                                                       args.min_read_quality, args.min_kmer_quality,
                                                       args.max_genomes)
        elif args.task == "dumpalign":
            if args.referencefile and args.reads:
                validate_file_readable(args.reads, "FASTQ reads")
                dump_alignment_from_reference(args.referencefile, args.reads,
                                              args.unique_threshold, args.ambiguous_threshold,
                                              args.min_read_quality, args.min_kmer_quality,
                                              args.max_genomes)
            elif args.genomefile and args.kmer_size and args.reads:
                validate_file_readable(args.reads, "FASTQ reads")
                validate_file_readable(args.genomefile, "Genome FASTA")
                build_reference_align_and_dump(args.genomefile, args.kmer_size, args.reads,
                                               args.unique_threshold, args.ambiguous_threshold,
                                               args.min_read_quality, args.min_kmer_quality,
                                               args.max_genomes, args.filter_similar,
                                               args.similarity_threshold)
            elif args.alignfile:
                validate_file_writable(args.alignfile, "Alignment output")
                dump_alignment_file(args.alignfile)
            else:
                sys.exit("Error: Provide either -g and -k with --reads, or -r with --reads, or -a.")
        else:
            sys.exit("Error: Unsupported task.")
    except gzip.BadGzipFile:
        sys.exit("Error: Incorrect format of input file.")
    except (InvalidExtensionError, NoRecordsInDataFile, NotValidatingUniqueMapping, AddingExistingRead) as err:
        sys.exit(err)


if __name__ == "__main__":
    main()
