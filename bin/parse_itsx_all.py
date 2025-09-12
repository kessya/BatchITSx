#!/usr/bin/env python3
import argparse
import csv
import logging
import re
from pathlib import Path
from Bio import SeqIO

parser = argparse.ArgumentParser(description="Script to format ITSx output (non-fungi sequences)")
parser.add_argument("--its2_fungi", required=True, help="ITS2 fasta file (fungi run)")
parser.add_argument("--its2_other", required=True, help="ITS2 fasta file (other run)")
parser.add_argument("--positions_fungi", required=True, help="ITSx positions.txt (fungi run)")
parser.add_argument("--positions_other", required=True, help="ITSx positions.txt (other run)")
parser.add_argument("--full", required=True, help="Full ITS fasta file (other run)")
parser.add_argument("--no_detect_fungi", required=False, help="ITSx no_detections.txt (fungi run)")
parser.add_argument("--no_detect_other", required=False, help="ITSx no_detections.txt (other run)")
parser.add_argument("--output", required=True, help="Output fasta file")
parser.add_argument("--region", choices=["its2", "itsfull"], required=True, help="Target region")
args = parser.parse_args()

outfile = Path(args.output)
log_file = outfile.with_suffix(".log")

logging.basicConfig(
    filename=log_file,
    filemode="a",
    format="%(asctime)s - %(levelname)s - %(message)s",
    level="INFO",
)

fungi_dict = {}
positions_dict = {}
length_dict = {}
new_positions_dict = {}
new_length_dict = {}
full_seq_dict = {}
no_coverage_count = 0
full_counter = 0
its2_counter = 0
len_limit = 140 if args.region == "itsfull" else 100
ex_o_ct = 0

# --------------------
# Load fungi ITS2 set
# --------------------
with open(args.its2_fungi, "r") as handle:
    for record in SeqIO.parse(handle, "fasta"):
        fungi_dict[record.id] = 1

# --------------------
# Parse fungi positions.txt
# --------------------
with open(args.positions_fungi) as pos:
    dataReader = csv.reader(pos, delimiter="\t")
    for row in dataReader:
        if row[0] == "--END--":
            continue
        if args.region == "itsfull":
            if (
                row[3] != "ITS1: Not found"
                and row[5] not in ["ITS2: Not found", "ITS2: No start", "ITS2: No end"]
            ):
                chim_match = re.search("Chimeric!", row[7])
                too_long_match = re.search("ITS region too long!", row[7])
                part_58S_match = re.search("Broken or partial sequence, only partial 5.8S!", row[7])
                no_58S_match = re.search("Broken or partial sequence, no 5.8S!", row[7])
                if not (chim_match or too_long_match or part_58S_match or no_58S_match):
                    positions_dict[row[0]] = 1
                    length_dict[row[0]] = int(row[1].split(" ")[0])
        elif args.region == "its2":
            if row[5] not in ["ITS2: Not found", "ITS2: No start", "ITS2: No end"]:
                chim_match = re.search("Chimeric!", row[7])
                too_long_match = re.search("ITS region too long!", row[7])
                if not (chim_match or too_long_match):
                    positions_dict[row[0]] = 1
                    length_dict[row[0]] = int(row[1].split(" ")[0])

print("Fungi pos: " + str(len(positions_dict)))

# --------------------
# Parse other positions.txt
# --------------------
with open(args.positions_other) as pos:
    dataReader = csv.reader(pos, delimiter="\t")
    for row in dataReader:
        if row[0] == "--END--":
            continue
        if args.region == "itsfull":
            if (
                row[3] != "ITS1: Not found"
                and row[5] not in ["ITS2: Not found", "ITS2: No start", "ITS2: No end"]
            ):
                chim_match = re.search("Chimeric!", row[7])
                too_long_match = re.search("ITS region too long!", row[7])
                part_58S_match = re.search("Broken or partial sequence, only partial 5.8S!", row[7])
                no_58S_match = re.search("Broken or partial sequence, no 5.8S!", row[7])
                if not (chim_match or too_long_match or part_58S_match or no_58S_match):
                    new_positions_dict[row[0]] = 1
                    new_length_dict[row[0]] = int(row[1].split(" ")[0])
                else:
                    logging.info(f"Excluded {row[0]}: Chimeric or broken sequence according to ITSx.")
            else:
                ex_o_ct += 1
                logging.info(f"Excluded {row[0]}: ITS1 or ITS2 sequence not detected.")
        elif args.region == "its2":
            if row[5] not in ["ITS2: Not found", "ITS2: No start", "ITS2: No end"]:
                chim_match = re.search("Chimeric!", row[7])
                too_long_match = re.search("ITS region too long!", row[7])
                if not (chim_match or too_long_match):
                    new_positions_dict[row[0]] = 1
                    new_length_dict[row[0]] = int(row[1].split(" ")[0])
                else:
                    logging.info(f"Excluded {row[0]}: Chimeric or broken sequence according to ITSx.")
            else:
                ex_o_ct += 1
                logging.info(f"Excluded {row[0]}: ITS2 sequence not detected.")

print("New pos: " + str(len(new_positions_dict)))

# --------------------
# Load full sequences (from "other" run)
# --------------------
with open(args.full, "r") as handle:
    for record in SeqIO.parse(handle, "fasta"):
        full_seq_dict[record.id] = str(record.seq)

print("Full seqs: " + str(len(full_seq_dict)))

# --------------------
# Filter sequences (non-fungal only)
# --------------------
with open(outfile, "w") as o, open(args.its2_other, "r") as handle:
    for record in SeqIO.parse(handle, "fasta"):
        record_id = record.id
        if record_id not in fungi_dict:
            if record_id in new_positions_dict:
                if new_length_dict[record_id] >= len_limit:
                    if record_id in full_seq_dict and len(full_seq_dict[record_id]) >= len_limit:
                        o.write(f">{record_id}\n{full_seq_dict[record_id]}\n")
                        full_counter += 1
                    elif len(str(record.seq)) >= len_limit:
                        o.write(f">{record_id}\n{record.seq}\n")
                        its2_counter += 1
                else:
                    logging.info(f"Excluded {record_id}: Sequence too short.")
            else:
                no_coverage_count += 1
                logging.info(f"Excluded {record_id}: ITS1 or ITS2 sequence not detected.")

# --------------------
# Handle no_detections.txt files
# --------------------
no_detect_f_set = set()
if args.no_detect_fungi and Path(args.no_detect_fungi).is_file():
    with open(args.no_detect_fungi, "r") as det:
        dataReader = csv.reader(det, delimiter="\t")
        for row in dataReader:
            no_detect_f_set.add(row[0])

if args.no_detect_other and Path(args.no_detect_other).is_file():
    with open(args.no_detect_other, "r") as det:
        dataReader = csv.reader(det, delimiter="\t")
        for row in dataReader:
            if row[0] in no_detect_f_set:
                logging.info(f"Excluded {row[0]}: Sequence not detected as ITS by ITSx.")

# --------------------
# Final logging
# --------------------
logging.info(f"No coverage for {no_coverage_count} sequences.")
logging.info(f"No. of full seqs: {full_counter}; No of ITS2 seqs: {its2_counter}")
logging.info(f"No of seqs (other, excluded no ITS detected): {ex_o_ct}")
