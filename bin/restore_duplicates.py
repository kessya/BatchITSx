#!/usr/bin/env python3
import argparse
from Bio import SeqIO

parser = argparse.ArgumentParser(description="Restore duplicates from vsearch dereplication")
parser.add_argument("--infile", required=True, help="Dereplicated fasta (combined)")
parser.add_argument("--uc", required=True, help="vsearch .uc file")
parser.add_argument("--outfile", required=True, help="Output full fasta with duplicates restored")
args = parser.parse_args()

# Parse uc file into mapping {seed_id: original_ids}
mapping = {}
with open(args.uc, "r") as f:
    for line in f:
        if line.startswith("#") or not line.strip():
            continue
        fields = line.strip().split("\t")
        record_type = fields[0]
        query_label = fields[8]  # sequence id
        target_label = fields[9] if len(fields) > 9 else None

        if record_type == "S":  # seed is its own representative
            rep = query_label
            mapping.setdefault(rep, []).append(rep)
        elif record_type == "H":  # duplicate, maps to representative
            rep = target_label
            dup = query_label
            mapping.setdefault(rep, []).append(dup)

# Write full fasta with duplicates restored
with open(args.outfile, "w") as out:
    with open(args.infile, "r") as handle:
        for record in SeqIO.parse(handle, "fasta"):
            rep_id = record.id
            if rep_id in mapping:
                for orig_id in mapping[rep_id]:
                    out.write(f">{orig_id}\n{str(record.seq)}\n")
