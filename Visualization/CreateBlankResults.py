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
import argparse
import logging
import time

from Bio import SeqIO
from shutil import copyfile

# Set up the logger
logging.basicConfig(level=logging.DEBUG, format="[ %(asctime)s UTC ]: %(levelname)s: %(module)s: %(message)s")
logging.Formatter.converter = time.gmtime
logger = logging.getLogger(__name__)

# ===========================================================================================================
# Functions:

# Checks whether or not the provided input_csv is an empty file. Returns logical.
def check_if_input_csv_is_empty(input_csv):
    if os.stat(input_csv).st_size == 0: # Returns True is the CSV is empty and False if not
        return True
    else:
        return False


# Checks whether or not the provided file has the expected extension. Throws warning if not. No return.
def file_extension_check(filename, extension):
    if not filename.endswith('.' + extension):
        logger.warning("'" + filename + "' may not be a '" + extension + "' file!")


# Uses input query_proteins FAA file to build a blank BLAST CSV table. Returns the table.
def generate_blank_results(query_proteins):
    blank_results_list = []

    with open(query_proteins, "rU") as fasta_file:
        fasta_entries = SeqIO.parse(fasta_file, "fasta")

        for entry in fasta_entries:
            blank_results_list.append(
                entry.id + ",NA,NA,NA,NA,NA")  # qseqid sseqid pident evalue qcovhsp bitscore

    blank_results = "\n".join(blank_results_list)
    return(blank_results)


# Writes a blank_results table to a file with name new_blast_file. No return.
def write_blank_results(blank_results, new_blast_file):
    with open(new_blast_file, "w") as writeFile:
        writeFile.write(blank_results)


def main(args):
    # Get user-provided variables
    original_blast_file = args.blast_results
    query_proteins = args.query_proteome
    new_blast_file = args.output_file

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
        logger.debug("Generating blank BLAST table based on FAA file provided")
        blank_results = generate_blank_results(query_proteins) # Stores file one for input checking.

        logger.debug("Writing to CSV file '" + new_blast_file + "'")
        write_blank_results(blank_results, new_blast_file)
    else:
        logger.error("Something went wrong. Exiting...")
        sys.exit(1)

    logger.debug("Done")


if __name__ == '__main__':
    """Command Line Interface Options"""

    parser = argparse.ArgumentParser(description = "A simple utility for working with BLAST results in batch. "
                                                   "Creates a BLAST results template based on the query_proteome "
                                                   "if the input blast file is blank. "
                                                   "Returns the original file if not blank. "
                                                   "Copyright Lee H. Bergstrand and Jackson M. Tsuji, 2018.")
    parser.add_argument('-i', '--blast_results', metavar='BLAST_IN', required=True,
                        help='''The path to CSV-format BLAST results (to be checked by this script if empty or not).''')

    parser.add_argument('-q', '--query_proteome', metavar='FASTA', required=True,
                        help='''The path to a protein FASTA file used as the original BLAST query.''')

    parser.add_argument('-o', '--output_file', metavar='BLAST_OUT', required=True,
                        help='''The path to write the output CSV-format BLAST results to.''')

    cli_args = parser.parse_args()
    main(cli_args)
