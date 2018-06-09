#!/usr/bin/env python

# ----------------------------------------------------------------------------------------
# Created by: Lee Bergstrand 
# Description: A simple program that takes a FASTA file protein list and makes a csv of fake BLAST results.
#
# Requirements: - This script requires the Biopython module: http://biopython.org/wiki/Download
#     
# Usage: getGenbankSeqs.py <sequences.txt> <geneList.fna>
# Example: getGenbankSeqs.py mySeqs.txt geneList.fna
# ----------------------------------------------------------------------------------------
# ===========================================================================================================

# Imports:
import sys

from Bio import SeqIO


# ===========================================================================================================
# Functions:

# 1: Checks if in proper number of arguments are passed gives instructions on proper use.
def argsCheck(numArgs):
    if len(sys.argv) < numArgs or len(sys.argv) > numArgs:
        print(
            "Takes a nucleotide FASTA file and returns the exact same FASTA file with a reverse complemented sequence.")
        print("By Lee Bergstrand\n")
        print("Usage: " + sys.argv[0] + " <sequences.txt> <geneList.fna>")
        print("Examples: " + sys.argv[0] + " mySeq.txt geneList.fna\n")
        exit(1)  # Aborts program. (exit(1) indicates that an error occurred)


# ===========================================================================================================
# Main program code:

# House keeping...
argsCheck(3)  # Checks if the number of arguments are correct.

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

# Reads sequence file list and stores it as a string object. Safely closes file.try:
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
