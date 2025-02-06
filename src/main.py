#!/usr/bin/env python3
import os
import argparse
import pickle
import gzip
import sys
import json

from data_file import FASTAFile, FASTAQFile, InvalidExtensionError, NoRecordsInDataFile
from kmer import (
    KmerReference,
    PseudoAlignment,
    NotValidatingUniqueMapping,
    AddingExistingRead,
)


def validate_file_readable(filepath: str, description: str):
    """Check if a file exists and is readable."""
    if not os.path.isfile(filepath):
        sys.exit(
            f"Error: {description} file '{filepath}' does not exist or is not a file."
        )
    if not os.access(filepath, os.R_OK):
        sys.exit(f"Error: {description} file '{filepath}' is not readable.")


def validate_file_writable(filepath: str, description: str):
    """Check if a file is writable (or if the directory is writable)."""
    dir_path = os.path.dirname(filepath) or "."
    if os.path.exists(filepath) and not os.access(filepath, os.W_OK):
        sys.exit(f"Error: {description} file '{filepath}' is not writable.")
    if not os.path.exists(filepath) and not os.access(dir_path, os.W_OK):
        sys.exit(
            f"Error: Directory '{dir_path}' is not writable to create {description} file '{filepath}'."
        )


def parse_arguments(args=None):
    parser = argparse.ArgumentParser(prog="Biosequence project")

    parser.add_argument("-t", "--task", required=True, help="Task to execute")
    parser.add_argument(
        "-g", "--genomefile", help="Genome FASTA file (multiple records)"
    )
    parser.add_argument("-k", "--kmer-size", type=int, help="Length of k-mers")
    parser.add_argument("-r", "--referencefile", help="KDB file (input/output)")

    # Task specific arguments
    parser.add_argument(
        "-a",
        "--alignfile",
        help="aln file. Can be either input or name for output file",
    )
    parser.add_argument("--reads", help="fastq reads file")
    parser.add_argument(
        "-m", "--unique-threshold", help="unique k-mer threshold", default=1, type=int
    )
    parser.add_argument(
        "-p",
        "--ambiguous-threshold",
        help="ambiguous k-mer threshold",
        default=1,
        type=int,
    )
    parser.add_argument("--reverse-complement", action="store_true")
    parser.add_argument("--min-read-quality", type=int, default=None)
    parser.add_argument("--min-kmer-quality", type=int, default=None)
    parser.add_argument("--max-genomes", type=int, default=None)
    parser.add_argument("--genomes")
    parser.add_argument("--coverage")
    parser.add_argument("--window-size", type=int, default=100)
    parser.add_argument("--min-coverage", type=int, default=1)
    parser.add_argument("--full-coverage", action="store_true")
    parser.add_argument("--filter-similar", action="store_true")
    parser.add_argument("--similarity-threshold", type=float, default=0.95)

    return parser.parse_args(args)


def create_reference(fasta_file, kmer_size):
    fasta_container = FASTAFile(fasta_file).container
    kmer_reference = KmerReference(kmer_size, fasta_container)
    return kmer_reference


def create_reference_and_save_it(fasta_file, kmer_size, reference_file):
    kmer_reference = create_reference(fasta_file, kmer_size)
    with gzip.open(reference_file, "wb") as f:
        pickle.dump(kmer_reference, f)


def dump_reference(kmer_reference):
    print(json.dumps(kmer_reference.get_summary(), indent=4))


def dump_reference_file(reference_file):
    with gzip.open(reference_file, "rb") as f:
        kmer_reference = pickle.load(f)
    dump_reference(kmer_reference)


def build_reference_and_dump_from_file(fasta_file, kmer_size):
    kmer_reference = create_reference(fasta_file, kmer_size)
    dump_reference(kmer_reference)


def create_alignment_from_reference(
    kmer_reference, reads_file, p, m, min_read_quality, min_kmer_quality, max_genomes
):
    reads_container = FASTAQFile(reads_file).container
    pseudo_alignment = PseudoAlignment(kmer_reference)
    pseudo_alignment.align_reads_from_container(
        reads_container, p, m, min_read_quality, min_kmer_quality, max_genomes
    )
    return pseudo_alignment


def create_alignment_file_from_reference(
    kmer_reference,
    reads_file,
    align_file,
    p,
    m,
    min_read_quality,
    min_kmer_quality,
    max_genomes
):
    pseudo_alignment = create_alignment_from_reference(
        kmer_reference,
        reads_file,
        p,
        m,
        min_read_quality,
        min_kmer_quality,
        max_genomes
    )
    pseudo_alignment.save(align_file)


def create_alignment_from_reference_file(
    reference_file,
    reads_file,
    align_file,
    p,
    m,
    min_read_quality,
    min_kmer_quality,
    max_genomes
):
    with gzip.open(reference_file, "rb") as f:
        kmer_reference = pickle.load(f)
    create_alignment_file_from_reference(
        kmer_reference,
        reads_file,
        align_file,
        p,
        m,
        min_read_quality,
        min_kmer_quality,
        max_genomes
    )


def build_reference_and_create_alignment_file(
    fasta_file,
    kmer_size,
    reads_file,
    align_file,
    m,
    p,
    min_read_quality,
    min_kmer_quality,
    max_genomes
):
    kmer_reference = create_reference(fasta_file, kmer_size)
    create_alignment_file_from_reference(
        kmer_reference,
        reads_file,
        align_file,
        m,
        p,
        min_read_quality,
        min_kmer_quality,
        max_genomes
    )


def dump_alignment_file(align_file):
    with gzip.open(align_file, "rb") as f:
        pseudo_alignment = pickle.load(f)
    print(json.dumps(pseudo_alignment.get_summary(), indent=4))


def dump_alignment_from_reference(
    reference_file, reads_file, m, p, min_read_quality, min_kmer_quality, max_genomes
):
    with gzip.open(reference_file, "rb") as f:
        kmer_reference = pickle.load(f)

    pseudo_alignment = create_alignment_from_reference(
        kmer_reference,
        reads_file,
        p,
        m,
        min_read_quality,
        min_kmer_quality,
        max_genomes
    )
    print(json.dumps(pseudo_alignment.get_summary(), indent=4))


def build_reference_align_and_dump(
    fasta_file,
    kmer_size,
    reads_file,
    m,
    p,
    min_read_quality,
    min_kmer_quality,
    max_genomes
):
    kmer_reference = create_reference(fasta_file, kmer_size)
    pseudo_alignment = create_alignment_from_reference(
        kmer_reference,
        reads_file,
        p,
        m,
        min_read_quality,
        min_kmer_quality,
        max_genomes
    )
    print(json.dumps(pseudo_alignment.get_summary(), indent=4))


def main():
    args = parse_arguments()

    try:
        if args.task == "reference":
            if not args.genomefile or not args.kmer_size or not args.referencefile:
                sys.exit(
                    "Error: -g (genome file), -k (kmer size), and -r (reference file) are required for this task."
                )
            validate_file_readable(args.genomefile, "Genome FASTA")
            validate_file_writable(args.referencefile, "Reference database output")
            create_reference_and_save_it(
                args.genomefile, args.kmer_size, args.referencefile
            )
        elif args.task == "dumpref":
            if args.referencefile:
                validate_file_readable(args.referencefile, "Reference database")
                dump_reference_file(args.referencefile)

            elif args.genomefile and args.kmer_size:
                validate_file_readable(args.genomefile, "Genome FASTA")
                build_reference_and_dump_from_file(args.genomefile, args.kmer_size)
            else:
                sys.exit(
                    "Error: Either -r (reference file) or -g (genome file) and -k (kmer size) are required for this task."
                )
        elif args.task == "align":
            validate_file_readable(args.reads, "FASTQ reads")
            validate_file_writable(args.alignfile, "Alignment output")

            if args.referencefile and args.reads and args.alignfile:
                validate_file_readable(args.referencefile, "Reference database")
                create_alignment_from_reference_file(
                    args.referencefile,
                    args.reads,
                    args.alignfile,
                    args.unique_threshold,
                    args.ambiguous_threshold,
                    args.min_read_quality,
                    args.min_kmer_quality,
                    args.max_genomes
                )
            elif args.genomefile and args.kmer_size and args.reads and args.alignfile:
                validate_file_readable(args.genomefile, "Genome FASTA")
                build_reference_and_create_alignment_file(
                    args.genomefile,
                    args.kmer_size,
                    args.reads,
                    args.alignfile,
                    args.unique_threshold,
                    args.ambiguous_threshold,
                    args.min_read_quality,
                    args.min_kmer_quality,
                    args.max_genomes
                )
            else:
                sys.exit(
                    "Error: Either provide -r (reference file) or -g (genome file) and -k (kmer size) along with --reads and -a (alignment output file)."
                )
        elif args.task == "dumpalign":
            if args.referencefile and args.reads:
                validate_file_readable(args.reads, "FASTQ reads")

                dump_alignment_from_reference(
                    args.referencefile,
                    args.reads,
                    args.unique_threshold,
                    args.ambiguous_threshold,
                    args.min_read_quality,
                    args.min_kmer_quality,
                    args.max_genomes
                )
            if args.genomefile and args.kmer_size and args.reads:
                validate_file_readable(args.reads, "FASTQ reads")
                validate_file_readable(args.genomefile, "Genome FASTA")
                build_reference_align_and_dump(
                    args.genomefile,
                    args.kmer_size,
                    args.reads,
                    args.unique_threshold,
                    args.ambiguous_threshold,
                    args.min_read_quality,
                    args.min_kmer_quality,
                    args.max_genomes
                )
            elif args.alignfile:
                validate_file_writable(args.alignfile, "Alignment output")
                dump_alignment_file(args.alignfile)
            else:
                sys.exit(
                    "Error: Provide either -g (genome file) and -k (kmer size) with --reads (FASTQ file), or -a (alignment output file)."
                )
        else:
            sys.exit("Error: Unsupported task.")
    except (
        InvalidExtensionError,
        NoRecordsInDataFile,
        NotValidatingUniqueMapping,
        AddingExistingRead,
    ) as err:
        sys.exit(err)


if __name__ == "__main__":
    main()
