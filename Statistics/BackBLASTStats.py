#!/usr/bin/env python 
# Created by: Lee Bergstrand 
# Description: Counts up BLAST hits per subject.
#
# Usage: BackBLASTStats.py <input.csv>
# Example: BackBLASTStats.py AUUJ0000000000.csv
# ----------------------------------------------------------------------------------------

# Imports
import csv
import sys

# ========================================================================================

# If in proper number of arguments are passed gives instructions on proper use.
if len(sys.argv) < 3 or len(sys.argv) > 3:
    print("Filters Blast Hits")
    print("By Lee Bergstrand\n")
    print("Filters Blast Hits to ...\n")
    print("Usage: " + sys.argv[0] + " <input.csv>")
    print("Examples: " + sys.argv[0] + ' AUUJ0000000000.csv')
    exit(1)  # Aborts program. (exit(1) indicates that an error occurred)

# Stores stores argument data.
print(">> Opening CSV Blast Result file...")
# ClusterBoundaries
BLASTResults = sys.argv[1]
outFile = BLASTResults + ".out"

# File extension check
if not BLASTResults.endswith(".csv"):
    print("[Warning] " + BLASTResults + " may not be a CSV file!")
else:
    print(">> Good file extension.")

# Opens CSV file for reading.
try:
    readFile = open(BLASTResults, "r")
    reader = csv.reader(
        readFile)  # opens file with csv module which takes into account varying csv formats and parses correctly
    print(">> Good CSV file.")

    BlastHits = []

    print(">> Reading BLAST hits...")
    # Modifies changes elements in column according to regex. New modified rows to new file.
    for row in reader:
        BlastHits.append(row)
    readFile.close()

    QueryCounter = {}

    for hit in BlastHits:
        QueryProtein = hit[0]
        if QueryProtein in QueryCounter:
            QueryCounter[QueryProtein] += 1
        else:
            QueryCounter.update({QueryProtein: 1})

    for x in QueryCounter:
        print(x + ":", QueryCounter[x])

except IOError:
    print("Failed to open " + BLASTResults)
    exit(1)
