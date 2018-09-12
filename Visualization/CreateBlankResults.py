#!/usr/bin/env python

# -----------------------------------------------------------------------------------------------------------
# Created by: Lee Bergstrand
# Description: A simple program that takes a FASTA file query and makes a csv of blank BLAST results.
#
# Requirements: - This script requires the Biopython module: http://biopython.org/wiki/Download
#     
# Usage: CreateBlankResults.py <sequence_files.txt> <query.faa>
# Example: CreateBlankResults.py sequences_files_to_replace.txt query_proteins.faa
# TODO - replace the original with the NEW below.
# NEW! Example: CreateBlankResults.py blast_output_file.csv query_proteins.faa > new_blast_output_file.csv
#       Will make a blank file if the input blast file is blank, and will otherwise return the original file.
# -----------------------------------------------------------------------------------------------------------
# ===========================================================================================================

# Imports:
import sys
import os

from Bio import SeqIO
from shutil import copyfile


# ===========================================================================================================
# Functions:

# 1: Checks if the proper number of arguments are passed gives instructions on proper use.
def argsCheck(numArgs):
    if len(sys.argv) < numArgs or len(sys.argv) > numArgs:
        print(
            "Takes a nucleotide FASTA file and returns the exact same FASTA file with a reverse complemented sequence.")
        print("By Lee Bergstrand\n")
        print("Usage: " + sys.argv[0] + "  <blast_output_file.csv> <query_proteins.faa> > <new_blast_output_file.csv>")
        exit(1)  # Aborts program. (exit(1) indicates that an error occurred)

# Get user-provided variables
blast_output_file = sys.argv[1]
query_proteins = sys.argv[2]

# ===========================================================================================================
# Functions:

def check_if_input_CSV_is_empty(input_csv):
    # Returns True is the CSV is empty and False if not

    if os.stat(input_csv).st_size == 0:
        return True
    elif os.stat(input_csv).st_size > 0:
        return False
    else
        print("ERROR: the input file appears to neither be empty nor non-empty. Exiting...")
        exit(1)


# Stores file one for input checking.
print(">> Opening FASTA file...")


def file_extension_check(input_filename):
    # Throws a warning if file extension does not appear to be appropriate

    if not blast_output_file.endswith(".csv"):
        print("[Warning] " + blast_output_file + " may not be a CSV file!")

    if not query_proteins.endswith(".faa"):
        print("[Warning] " + query_proteins + " may not be an amino acid FASTA file!")

# Reads BLAST output file and stores it as a string.
# TODO - should not be a string
try:
    with open(blast_output_file, "rU") as newFile:
        blast_output_file = newFile.read()
        newFile.close()
except IOError:
    print("Failed to open " + blast_output_file)
    exit(1)

# Splits string into a list. Each element is a single line from the string.
filesToReplaceList = filesToReplace.splitlines()

FakeResults = []
try:
    handle = open(query_proteins, "rU")
    SeqRecords = SeqIO.parse(handle, "fasta")
    for record in SeqRecords:
        FakeResults.append(
            record.id + ",sseqid,0,evalue,qcovhsp,bitscore")  # qseqid sseqid pident evalue qcovhsp bitscore
    handle.close()
except IOError:
    print("Failed to open " + blast_output_file + " or " + query_proteins)
    exit(1)

FakeResultsOut = "\n".join(FakeResults)

print(">> Writing Fake Results...")

for x in filesToReplaceList:
    replacementFileName = x + ".csv"

    try:
        writeFile = open(replacementFileName, "w")
        writeFile.write(FakeResultsOut)
        writeFile.close()
    except IOError:
        print("Failed to create " + replacementFileName)
        exit(1)

print(">> Done.")

def main(args):
    # ===========================================================================================================
    # Main program code:

    # House keeping...
    argsCheck(3)  # Checks if the number of arguments are correct.

    # TODO - make input_csv and output_csv variables (global? Then capitalize)

    # If the CSV file is not empty, then just keep the file as is.
    if check_if_input_CSV_is_empty(input_csv) == False:
        copyfile(input_csv, output_csv)

    # If it is empty, then make a fake file.
    elif check_if_input_CSV_is_empty(input_csv) == True:
        # TODO - keep moving the above code into functions and then get it into this block. Then add proper argument parser below...
        # Doesn't have to be pretty.
        

