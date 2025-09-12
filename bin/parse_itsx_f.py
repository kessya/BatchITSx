#!/usr/bin/env python3
import argparse
import csv
import logging
import re
from pathlib import Path
from Bio import SeqIO

parser = argparse.ArgumentParser(description="Script to format ITSx output (fungi)")
parser.add_argument("--its2", required=True, help="ITS2 fasta file from ITSx")
parser.add_argument("--positions", required=True, help="ITSx positions.txt file")
parser.add_argument("--full", required=True, help="Full ITS fasta file from ITSx")
parser.add_argument("--output", required=True, help="Output fasta file")
parser.add_argument("--region", choices=["its2", "itsfull"], required=True, help="Target region")
args = parser.parse_args()

# input files
its2_file = Path(args.its2)
pos_file = Path(args.positions)
full_file = Path(args.full)
outfile = Path(args.output)

# logging setup
log_file = outfile.with_suffix(".log")
logging.basicConfig(
    filename=log_file,
    filemode="a",
    format="%(asctime)s - %(levelname)s - %(message)s",
    level="INFO",
)

positions_dict = {}
length_dict = {}
full_seq_dict = {}
no_coverage_count = 0
full_counter = 0
its2_counter = 0
ex_f_ct = 0
len_limit = 140 if args.region == "itsfull" else 100

# parse positions.txt
with open(pos_file) as pos:
    dataReader_pos = csv.reader(pos, delimiter="\t")
    for row in dataReader_pos:
        if row[0] == "--END--":
            continue

        if args.region == "itsfull":
            if (
                row[3] != "ITS1: Not found"
                and row[5] not in ["ITS2: Not found", "ITS2: No start", "ITS2: No end"]
            ):
                chim_match = re.search("Chimeric!", row[7])
                part_58S_match = re.search("Broken or partial sequence, only partial 5.8S!", row[7])
                no_58S_match = re.search("Broken or partial sequence, no 5.8S!", row[7])
                too_long_match = re.search("ITS region too long!", row[7])
                if not (chim_match or part_58S_match or no_58S_match or too_long_match):
                    positions_dict[row[0]] = 1
                    length_dict[row[0]] = int(row[1].split(" ")[0])
                else:
                    logging.info(f"Excluded {row[0]}: Chimeric or broken sequence according to ITSx.")
            else:
                ex_f_ct += 1
                logging.info(f"Excluded {row[0]}: ITS1 or ITS2 sequence not detected.")
        elif args.region == "its2":
            if row[5] not in ["ITS2: Not found", "ITS2: No start", "ITS2: No end"]:
                chim_match = re.search("Chimeric!", row[7])
                too_long_match = re.search("ITS region too long!", row[7])
                if not (chim_match or too_long_match):
                    positions_dict[row[0]] = 1
                    length_dict[row[0]] = int(row[1].split(" ")[0])
                else:
                    logging.info(f"Excluded {row[0]}: Chimeric or broken sequence according to ITSx.")
            else:
                ex_f_ct += 1
                logging.info(f"Excluded {row[0]}: ITS2 sequence not detected.")

print("Positions dict size:", len(positions_dict))

# load full sequences
with open(full_file, "r") as handle:
    for record in SeqIO.parse(handle, "fasta"):
        record_id = record.id
        full_seq_dict[record_id] = str(record.seq)

# filter sequences
with open(outfile, "w") as o, open(its2_file, "r") as handle:
    for record in SeqIO.parse(handle, "fasta"):
        record_id = record.id
        if record_id in positions_dict:
            if length_dict[record_id] >= len_limit:
                if record_id in full_seq_dict and len(full_seq_dict[record_id]) >= len_limit:
                    o.write(f">{record_id}\n{full_seq_dict[record_id]}\n")
                    full_counter += 1
                else:
                    if len(str(record.seq)) >= len_limit:
                        o.write(f">{record_id}\n{record.seq}\n")
                        its2_counter += 1
            else:
                logging.info(f"Excluded {record_id}: Sequence too short.")
        else:
            no_coverage_count += 1

# summary log
logging.info(f"No coverage for {no_coverage_count} sequences.")
logging.info(f"No. of full seqs: {full_counter}; No of ITS2 seqs: {its2_counter}")
logging.info(f"No of seqs excluded (fungi, no ITS region detected): {ex_f_ct}")
