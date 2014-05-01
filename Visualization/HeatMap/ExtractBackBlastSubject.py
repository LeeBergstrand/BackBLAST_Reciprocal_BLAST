#!/usr/bin/env python 
# Created by: Lee Bergstrand 
# Descript: A simple script that extracts the subject BLAST hits for a given subject protien from the csv results from BackBLAST.py.
#
# Usage: ExtractBackBlastSubject.py <myInputList.txt> <myInputBLASTCSV.csv>
# Example: ExtractBackBlastSubject.py myInputList.txt myInputBLASTCSV.csv 
#----------------------------------------------------------------------------------------

import csv
import sys
import re

# If in proper number of arguments are passed gives instructions on proper use.
if len(sys.argv) < 3 or len(sys.argv) > 3:
	print "BLAST Subject Extractor"
	print "By Lee Bergstrand\n"
	print "A simple script that extracts the subject BLAST hits for a given subject protien from the csv results from BackBLAST.py.\n"
	print "Usage: " + sys.argv[0] + "<myInputList.txt> <myInputBLASTCSV.csv>"
	print "Examples: " + sys.argv[0] + 'myInputList.txt myInputBLASTCSV.csv'
	exit(1) # Aborts program. (exit(1) indicates that an error occured)

# Store argument data.
inListFile  = sys.argv[1]
inCSVFile   = sys.argv[2]
outListFile = sys.argv[1] + ".out"

print ">> Reading accession list from " + inListFile + "..." # File extension check
if not inListFile.endswith(".txt"):
	print "[Warning] " + inListFile + "May not be a txt file!"
else:
	print ">> Good file extention."	

print ">> Reading back blast results from " + inCSVFile + "..." # File extension check
if not inCSVFile.endswith(".csv"):
	print "[Warning] " + inCSVFile + "May not be a CSV file!"
else:
	print ">> Good file extention."	
	
# Reads query sequence file list and stores it as a string object. Safely closes file.try:
try:	
	with open(inListFile,"r") as newFile:
		sequences = newFile.read()
		seqList = sequences.splitlines() # Splits string into a list. Each element is a single line from the string.
		newFile.close()
except IOError:
	print "Failed to open " + inListFile
	exit(1)

# Opens BLAST CSV output for reading.
try:
	csvFile = open(inCSVFile, "r")
	BLASTreader  = csv.reader(csvFile) # opens file with csv module which takes into account verying csv formats and parses correctly
	print ">> Good CSV file."
except IOError:
	print "Failed to open " + inCSVFile
	exit(1)

accessionsToAdd = []	
	
for row in BLASTreader:
	if row[0] in seqList:
		if float(row[2]) >= 30:
			if (row[1] + "\n") not in accessionsToAdd:
				accessionsToAdd.append(row[1] + "\n")
				
try:
	print ">> Output file created."
	print ">> Writing Data..."
	writeFile = open(outListFile, "w") 	 
	writeFile.writelines(accessionsToAdd)	
	writeFile.close()
except IOError:
	print "Failed to create " + outListFile
	exit(1)