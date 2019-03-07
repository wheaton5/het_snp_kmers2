#!/usr/bin/env python

import argparse

parser = argparse.ArgumentParser(description="filter het kmer set on coverage thresholds")

parser.add_argument("input",help="input kmer set file as output by het_snp_kmers")

parser.add_argument("--min", help = "min coverage",default=0)
parser.add_argument("--max", help = "max coverage",default=100000000)
parser.add_argument("--max_offsite", help = "max value of 3rd+4th middle base kmers",default=100000)
parser.add_argument("--max_sum", help = "max sum of coverage",default = 1000000)
parser.add_argument("--min_sum", help = "min sum of coverage",default =0)

args = parser.parse_args()

with open(args.input) as kmers:
    for line in kmers:
        tokens = line.strip().split()
        c1 = int(tokens[1])
        c2 = int(tokens[3])
        ct = c1+c2
        ce = int(tokens[4])
        if c1 >= args.min and c2 >= args.min and c1 <= args.max and c2 <= args.max and ct <= args.max_sum and ct >= args.min_sum and ce <= args.max_offsite:
            print(line.strip())
