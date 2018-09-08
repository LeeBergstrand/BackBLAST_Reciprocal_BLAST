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

# 1: Checks if in proper number of arguments are passed gives instructions on proper use.
def argsCheck(numArgs):
    if len(sys.argv) < numArgs or len(sys.argv) > numArgs:
        print(
            "Takes a nucleotide FASTA file and returns the exact same FASTA file with a reverse complemented sequence.")
        print("By Lee Bergstrand\n")
        print("Usage: " + sys.argv[0] + "  <sequence_files.txt> <query.faa>")
        print("Examples: " + sys.argv[0] + " sequences_files_to_replace.txt query_proteins.faa\n")
        exit(1)  # Aborts program. (exit(1) indicates that an error occurred)


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
filesToReplace = sys.argv[1]
geneList = sys.argv[2]

# File extension check
if not filesToReplace.endswith(".txt"):
    print("[Warning] " + geneList + " may not be a text file!")

# File extension check
if not geneList.endswith(".faa"):
    print("[Warning] " + geneList + " may not be a FASTA file!")

# Reads sequence file list and stores it as a string.
try:
    with open(filesToReplace, "rU") as newFile:
        filesToReplace = newFile.read()
        newFile.close()
except IOError:
    print("Failed to open " + filesToReplace)
    exit(1)

# Splits string into a list. Each element is a single line from the string.
filesToReplaceList = filesToReplace.splitlines()

FakeResults = []
try:
    handle = open(geneList, "rU")
    SeqRecords = SeqIO.parse(handle, "fasta")
    for record in SeqRecords:
        FakeResults.append(
            record.id + ",sseqid,0,evalue,qcovhsp,bitscore")  # qseqid sseqid pident evalue qcovhsp bitscore
    handle.close()
except IOError:
    print("Failed to open " + filesToReplace + " or " + geneList)
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
        #

