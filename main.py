#!/usr/bin/env python3

import argparse

def readargs(args=None):
    parser = argparse.ArgumentParser(
                    prog='Biosequence project',
    )
    #General arguments
    parser.add_argument('-t', '--task',
                        help="task",
                        required=True
                        )
    parser.add_argument('-g', '--genomefile',
                        help="Genome fasta file (multiple records)",
                        )
    parser.add_argument('-r', '--referencefile',
                        help="kdb file. Can be either input or name for output file",
                        )
    parser.add_argument('-k', '--kmer-size',
                        help="length of kmers",
                        )
    #Task specific arguments
    #align
    parser.add_argument('-a', '--alignfile',
                        help="aln file. Can be either input or name for output file",
                        )
    parser.add_argument('--reads',
                        help="fastq reads file",
                        )
    parser.add_argument('-m', '--unique-threshold',
                        help="unique k-mer threshold",
                        )
    parser.add_argument('-p', '--ambiguous-threhold',
                        help="ambiguous k-mer threshold",
                        )
                        
    #align+rc
    parser.add_argument('--reverse-complement',
                        action='store_true',
                        )
    #align+quality
    parser.add_argument('--min-read-quality',
                        type=int,
                        )
    parser.add_argument('--min-kmer-quality',
                        type=int,
                        )
    parser.add_argument('--max-genomes',
                        type=int,
                        )
    #coverage
    parser.add_argument('--genomes',
                        )
    parser.add_argument('--coverage', # fixed regions-->coverage
                        )
    parser.add_argument('--window-size',
                        type=int,
                        default=100,
                        )
    parser.add_argument('--min-coverage',
                        type=int,
                        default=1,
                        )
    parser.add_argument('--full-coverage',
                        action='store_true',
                        )
    #similarity
    parser.add_argument('--filter-similar',
                        action='store_true',
                        )
    parser.add_argument('--similarity-threshold',
                        type=float,
                        default=0.95,
                        )
    
    return parser.parse_args(args)

def main():
    args=readargs()
    ...

if __name__=="__main__":
    main()
