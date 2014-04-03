#!/usr/bin/env python 
# Created by: Lee Bergstrand 
# Descript: Counts up BLAST hits per subject.
#
# Usage: QueryProtienRepeats.py <input.csv> 
# Example: QueryProtienRepeats.py AUUJ0000000000.csv
#----------------------------------------------------------------------------------------
# Imports 
import csv
import sys
#========================================================================================

# If in proper number of arguments are passed gives instructions on proper use.
if len(sys.argv) < 2 or len(sys.argv) > 2:
	print "Filtres Blast Hits"
	print "By Lee Bergstrand\n"
	print "Filtres Blast Hits to ...\n"
	print "Usage: " + sys.argv[0] + " <input.csv>"
	print "Examples: " + sys.argv[0] + ' AUUJ0000000000.csv'
	exit(1) # Aborts program. (exit(1) indicates that an error occured)

# Stores stores argument data.
print ">> Opening CSV Blast Result file..."
inFile  = sys.argv[1]
outFile = inFile + ".out"

# File extension check
if not inFile.endswith(".csv"):
	print "[Warning] " + inFile + " may not be a CSV file!"
else:
	print ">> Good file extention."	
	
# Opens CSV file for reading.
try:
	readFile = open(inFile, "r")
	reader  = csv.reader(readFile) # opens file with csv module which takes into account verying csv formats and parses correctly
	print ">> Good CSV file."
except IOError:
	print "Failed to open " + inFile
	exit(1)

BlastHits = []

print ">> Reading BLAST hits..."
# Modifies changes elements in column according to regex. New modified rows to new file.
for row in reader:
	BlastHits.append(row)
readFile.close()

QueryCounter = {}

for hit in BlastHits:
	QueryProtien = hit[0]
	if QueryProtien in QueryCounter:
		QueryCounter[QueryProtien] += 1
	else:
		QueryCounter.update({QueryProtien:1})

for x in QueryCounter:
	print x + ":", QueryCounter[x]