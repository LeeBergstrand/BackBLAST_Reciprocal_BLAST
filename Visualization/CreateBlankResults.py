#!/usr/bin/env python

# -----------------------------------------------------------------------------------------------------------
# Created by: Lee Bergstrand
# Description: A simple program that takes a FASTA file query and makes a csv of blank BLAST results.
#
# Requirements: - This script requires the Biopython module: http://biopython.org/wiki/Download
#
# Usage: CreateBlankResults.py blast_output_file.csv query_proteins.faa new_blast_output_file.csv
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

# Checks if the proper number of arguments are passed gives instructions on proper use.
def argsCheck(numArgs):
    if len(sys.argv) < numArgs or len(sys.argv) > numArgs:
        print(
            "Takes a nucleotide FASTA file and returns the exact same FASTA file with a reverse complemented sequence.")
        print("By Lee Bergstrand\n")
        print("Usage: " + sys.argv[0] + "  <blast_output_file.csv> <query_proteins.faa> <new_blast_output_file.csv>")
        exit(1)  # Aborts program. (exit(1) indicates that an error occurred)


def check_if_input_CSV_is_empty(input_csv):
    # Returns True is the CSV is empty and False if not

    if os.stat(input_csv).st_size == 0:
        return True

    elif os.stat(input_csv).st_size > 0:
        return False

    else
        print("ERROR: the input file appears to neither be empty nor non-empty. Exiting...")
        exit(1)


def file_extension_check(filename, extension):
    # Throws a warning if file extension does not appear to be appropriate

    if extension == "csv":
        if not filename.endswith(".csv"):
            print("[Warning] " + filename + " may not be a CSV file!")

    elif extension == "faa":
        if not filename.endswith(".faa"):
            print("[Warning] " + filename + " may not be an amino acid FASTA file!")

    else:
        print("[ERROR]: extension must be 'csv' or 'faa'. Exiting...")
        exit(1)


def generate_blank_results(query_proteins):

    blank_results_list = []

    try:
        handle = open(query_proteins, "rU")
        fasta_entries = SeqIO.parse(handle, "fasta")
        for entry in fasta_entries:
            blank_results_list.append(
                record.id + ",sseqid,0,evalue,qcovhsp,bitscore")  # qseqid sseqid pident evalue qcovhsp bitscore
        handle.close()

    except IOError:
        print("Failed to open " + query_proteins)
        exit(1)

    blank_results = "\n".join(blank_results_list)

    return(blank_results)


def write_blank_results(blank_results, new_blast_file):

    try:
        writeFile = open(new_blast_file, "w")
        writeFile.write(blank_results)
        writeFile.close()

    except IOError:
        print("Failed to create " + new_blast_file)
        exit(1)


def main(args):
    # ===========================================================================================================
    # Main program code:

    # House keeping...
    argsCheck(3)  # Checks if the number of arguments are correct.

    # Get user-provided variables
    original_blast_file = sys.argv[1]
    query_proteins = sys.argv[2]
    new_blast_file = sys.argv[3]

    # If the CSV file is not empty, then just keep the file as is.
    if check_if_input_CSV_is_empty(original_blast_file) == False:
        copyfile(original_blast_file, new_blast_file)

    # If it is empty, then make a fake file.
    elif check_if_input_CSV_is_empty(original_blast_file) == True:

        # Check file extensions
        file_extension_check(original_blast_file, "csv")
        file_extension_check(query_proteins, "faa")
        file_extension_check(new_blast_file, "csv")

        # Stores file one for input checking.
        print(">> Generating blank BLAST table based on FAA file provided")
        blank_results = generate_blank_results(query_proteins)

        print(">> Writing to CSV file '" + new_blast_file + "'")
        write_blank_results(blank_results, new_blast_file)

        print(">> Done")

    else:
        print("Something went wrong. Exiting...")
        exit(1)

