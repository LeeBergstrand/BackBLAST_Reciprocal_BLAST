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
import datetime

from Bio import SeqIO
from shutil import copyfile


# ===========================================================================================================
# Functions:

# Checks if the proper number of arguments are passed gives instructions on proper use.
def argsCheck(numArgs):
    if len(sys.argv) < numArgs or len(sys.argv) > numArgs:
        print(
            "Takes a nucleotide FASTA file and returns the exact same FASTA file with a reverse complemented sequence.")
        print("Copyright Lee H. Bergstrand and Jackson M. Tsuji, 2018\n")
        print("Usage: " + sys.argv[0] + "  <blast_output_file.csv> <query_proteins.faa> <new_blast_output_file.csv>")
        exit(1)  # Aborts program. (exit(1) indicates that an error occurred)


# Gives a timestamp (string) for help with logging
def ts():
    current_time = datetime.datetime.utcnow()
    time_formatted = current_time.strftime("%a %b %d %H:%M:%S UTC %Y")

    timestamp = "[ " + time_formatted + " ]: "

    return(timestamp)


# Checks whether or not the provided input_csv is an empty file. Returns logical.
def check_if_input_CSV_is_empty(input_csv):
    # Returns True is the CSV is empty and False if not

    if os.stat(input_csv).st_size == 0:
        return True

    elif os.stat(input_csv).st_size > 0:
        return False

    else:
        print(ts() + "ERROR: the input file appears to neither be empty nor non-empty. Exiting...")
        exit(1)


# Checks whether or not the provided file has the expected extension. Throws warning if not. No return.
def file_extension_check(filename, extension):
    # Throws a warning if file extension does not appear to be appropriate

    if extension == "csv":
        if not filename.endswith(".csv"):
            print(ts() + "WARNING: '" + filename + "' may not be a CSV file!")

    elif extension == "faa":
        if not filename.endswith(".faa"):
            print(ts() + "WARNING: '" + filename + "' may not be an amino acid FASTA file!")

    else:
        print(ts() + "ERROR: extension must be 'csv' or 'faa'. Exiting...")
        exit(1)


# Uses input query_proteins FAA file to build a blank BLAST CSV table. Returns the table.
def generate_blank_results(query_proteins):

    blank_results_list = []

    try:
        handle = open(query_proteins, "rU")
        fasta_entries = SeqIO.parse(handle, "fasta")
        for entry in fasta_entries:
            blank_results_list.append(
                entry.id + ",NA,NA,NA,NA,NA")  # qseqid sseqid pident evalue qcovhsp bitscore
        handle.close()

    except IOError:
        print(ts() + "Failed to open " + query_proteins)
        exit(1)

    blank_results = "\n".join(blank_results_list)

    return(blank_results)


# Writes a blank_results table to a file with name new_blast_file. No return.
def write_blank_results(blank_results, new_blast_file):

    try:
        writeFile = open(new_blast_file, "w")
        writeFile.write(blank_results)
        writeFile.close()

    except IOError:
        print(ts() + "Failed to create " + new_blast_file)
        exit(1)


def main(args):
    # Housekeeping
    argsCheck(4)  # Checks if the number of arguments are correct.

    # Get user-provided variables
    original_blast_file = sys.argv[1]
    query_proteins = sys.argv[2]
    new_blast_file = sys.argv[3]

    # Print settings
    print(ts() + "Running "+ sys.argv[0])
    print(ts() + "Settings: Input BLAST CSV table: " + original_blast_file)
    print(ts() + "Settings: Query protein FastA file: " + query_proteins)
    print(ts() + "Settings: Output BLAST CSV table: " + new_blast_file)

    # Check file extensions
    file_extension_check(original_blast_file, "csv")
    file_extension_check(query_proteins, "faa")
    file_extension_check(new_blast_file, "csv")

    # If the CSV file is not empty, then just keep the file as is.
    if check_if_input_CSV_is_empty(original_blast_file) == False:

        print(ts() + "Provided BLAST file has content - no need to replace. Copying to output file.")
        copyfile(original_blast_file, new_blast_file)

    # If it is empty, then make a fake file.
    elif check_if_input_CSV_is_empty(original_blast_file) == True:

        # Stores file one for input checking.
        print(ts() + "Generating blank BLAST table based on FAA file provided")
        blank_results = generate_blank_results(query_proteins)

        print(ts() + "Writing to CSV file '" + new_blast_file + "'")
        write_blank_results(blank_results, new_blast_file)

    else:
        print(ts() + "Something went wrong. Exiting...")
        exit(1)

    print(ts() + "Done")

main(sys.argv)
