#!/usr/bin/env python3
import argparse
import re
from Bio import SeqIO

parser = argparse.ArgumentParser(description="Filter sequences with too many ambiguous bases (N/IUPAC)")
parser.add_argument("--infile", required=True, help="Input FASTA (duplicates already restored)")
parser.add_argument("--outfile", required=True, help="Output FASTA after filtering")
parser.add_argument("--excluded", help="Optional: write list of excluded sequence IDs")
parser.add_argument("--maxN", type=int, default=6, help="Maximum allowed Ns per sequence")
parser.add_argument("--maxIUPAC", type=int, default=16, help="Maximum allowed ambiguous IUPAC codes per sequence")
args = parser.parse_args()

# Regex for ambiguities
iupac_pattern = re.compile("[YRSWKMBDHVN]")
n_pattern = re.compile("N")

def passes_filter(seq: str) -> bool:
    n_count = len(n_pattern.findall(seq))
    iupac_count = len(iupac_pattern.findall(seq))
    return n_count <= args.maxN and iupac_count <= args.maxIUPAC

excluded_ids = []

with open(args.outfile, "w") as fout:
    for record in SeqIO.parse(args.infile, "fasta"):
        seq = str(record.seq).upper()
        if passes_filter(seq):
            fout.write(f">{record.id}\n{seq}\n")
        else:
            excluded_ids.append(record.id)

# Optionally write excluded sequences
if args.excluded:
    with open(args.excluded, "w") as ex:
        for seqid in excluded_ids:
            ex.write(f"{seqid}\n")
