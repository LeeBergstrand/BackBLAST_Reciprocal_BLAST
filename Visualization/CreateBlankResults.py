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
import logging
import time

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
        sys.exit(1)


# Checks whether or not the provided input_csv is an empty file. Returns logical.
def check_if_input_csv_is_empty(input_csv):
    # Returns True is the CSV is empty and False if not

    if os.stat(input_csv).st_size == 0:
        return True

    else:
        return False



# Checks whether or not the provided file has the expected extension. Throws warning if not. No return.
def file_extension_check(filename, extension):
    # Throws a warning if file extension does not appear to be appropriate

    if not filename.endswith('.' + extension):
        logger.warning("'" + filename + "' may not be a '" + extension + "' file!")


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
        logger.error("Failed to open " + query_proteins)
        sys.exit(1)

    blank_results = "\n".join(blank_results_list)

    return(blank_results)


# Writes a blank_results table to a file with name new_blast_file. No return.
def write_blank_results(blank_results, new_blast_file):

    try:
        writeFile = open(new_blast_file, "w")
        writeFile.write(blank_results)
        writeFile.close()

    except IOError:
        logger.error("Failed to create " + new_blast_file)
        sys.exit(1)


def main(args):
    # Set up the logger
    logging.basicConfig(level=logging.DEBUG, format="[ %(asctime)s UTC ]: %(levelname)s: %(module)s: %(message)s")
    logging.Formatter.converter = time.gmtime
    logger = logging.getLogger(__name__)  # Not sure if this is needed

    # Check if the number of arguments are correct.
    argsCheck(4)

    # Get user-provided variables
    original_blast_file = sys.argv[1]
    query_proteins = sys.argv[2]
    new_blast_file = sys.argv[3]

    # Report the input arguments
    logger.debug("Running " + os.path.basename(sys.argv[0]))
    logger.debug("Settings: Input BLAST CSV table: " + original_blast_file)
    logger.debug("Settings: Query protein FastA file: " + query_proteins)
    logger.debug("Settings: Output BLAST CSV table: " + new_blast_file)

    # Check file extensions
    file_extension_check(original_blast_file, "csv")
    file_extension_check(query_proteins, "faa")
    file_extension_check(new_blast_file, "csv")

    # If the CSV file is not empty, then just keep the file as is.
    if check_if_input_csv_is_empty(original_blast_file) == False:

        logger.debug("Provided BLAST file has content - no need to replace. Copying to output file.")
        copyfile(original_blast_file, new_blast_file)

    # If it is empty, then make a fake file.
    elif check_if_input_csv_is_empty(original_blast_file) == True:

        # Stores file one for input checking.
        logger.debug("Generating blank BLAST table based on FAA file provided")
        blank_results = generate_blank_results(query_proteins)

        logger.debug("Writing to CSV file '" + new_blast_file + "'")
        write_blank_results(blank_results, new_blast_file)

    else:
        logger.error("Something went wrong. Exiting...")
        sys.exit(1)

    logger.debug("Done")

main(sys.argv)
