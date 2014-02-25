#!/usr/bin/env python 
# Created by: Lee Bergstrand
# Descript: A Bio-Python program that takes a list of query proteins and uses local BLASTp to search
#           for highly similer proteins within a local blast database (usally a local db of a target 
#           proteome). The program then BLASTps backward from the found subject protein to the proteome 
#           for which the original query protein if found in order to confirm gene orthology. 
#             
# Requirements: - This program requires the Biopython module: http://biopython.org/wiki/Download
#               - This script requires BLAST+ 2.2.9 or later.
#               - All operations are done with protien sequences.
#               - All query proteins should be from sequenced genomes in order to facilitate backwards BLAST. 
#               - MakeBlastDB must be used to create BLASTp databases for both query and subject proteomes.
#               - BLAST databases require the FASTA file they were made from to be in the same directory.
#  
# Usage: BackBLAST.py <queryGeneList.faa> <queryProteomes.csv> <subject1.faa> 
# Example: BackBLAST.py queryGeneList.faa queryProteomes.csv AUUJ00000000.faa
#----------------------------------------------------------------------------------------
#===========================================================================================================
#Imports:
import sys
import csv
import subprocess
from Bio import SeqIO
#===========================================================================================================
# Functions:

# 1: Checks if in proper number of arguments are passed gives instructions on proper use.
def argsCheck():
	if len(sys.argv) < 4:
		print "Orthologous Gene Finder"
		print "By Lee Bergstrand\n"
		print "Please refer to source code for documentation\n"
		print "Usage: " + sys.argv[0] + " <queryProteomes.csv> <queryGeneList.faa> <subject1.faa>\n"
		print "Examples:" + sys.argv[0] + " queryProteomes.csv queryGeneList.faa AUUJ0000000.faa"
		exit(1) # Aborts program. (exit(1) indicates that an error occured)
#-------------------------------------------------------------------------------------------------
# 2: Runs BLAST, can either be sent a fasta formatted string or a file ...
def runBLAST(query, BLASTDBFile):
	
	query = query.strip()
	
	# If the string starts with > it is a fasta file and should be sent to blast via stdin.
	if query.startswith(">"):
		BlastArgs = "blastp -db "+ BLASTDBFile + " -num_threads 16 -outfmt \"10 qseqid sseqid pident evalue qcovhsp score\""
		blastp = subprocess.Popen(BlastArgs, shell=True, stdin=subprocess.PIPE,stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
		BLASTOut, err = blastp.communicate(query)
	
	# Else the string is a file name and should be submitted as a file argument.
	else:
		BLASTOut = subprocess.check_output(["blastp", "-db", BLASTDBFile, "-query", query, "-evalue", "1e-40", "-num_threads", "16", "-outfmt", "10 qseqid sseqid pident evalue qcovhsp score"]) # Runs BLASTp and save output to a string. Blastp is set to output csv which can be parsed.
	return BLASTOut
#-------------------------------------------------------------------------------------------------
# 3: Filters HSPs by Percent Identity...
def filtreBLASTCSV(BLASTOut, minIdent):
	
	BLASTCSVOut = BLASTOut.splitlines(True) # Converts raw BLAST csv output into list of csv rows.
	BLASTreader = csv.reader(BLASTCSVOut) # Reads BLAST csv rows as a csv.

	BLASTCSVOutFiltred = [] # Note should simply delete unwanted HSPs from curent list rather than making new list. 
					        # Rather than making a new one.
	for HSP in BLASTreader:
		if HSP[2] >= minIdent: # Filtres by minimum ident.
			# Converts each HSP parameter that should be a number to a number.
			HSP[2] = float(HSP[2]) 
			HSP[3] = float(HSP[3])
			HSP[4] = int(HSP[4])
			HSP[5] = int(HSP[5]) 
			BLASTCSVOutFiltred.append(HSP) # Appends to output array.
	
	return BLASTCSVOutFiltred
#-------------------------------------------------------------------------------------------------
# 4: Finds Top Scoring Hit For Each Query Protien In BLAST result... Could be more elegantly done...
def getTopHits(BLASTCSVOut): 
	
	topHits = []
	topHits.append(BLASTCSVOut[0])	
	currentQuery = BLASTCSVOut[0]
	
	# If ties occur include these tied in the top hit list. (Ties should have the same score)
	for x in range(1, len(BLASTCSVOut)):
		if BLASTCSVOut[x][5] == currentQuery[5]:			topHits.append(BLASTCSVOut[x])
		else:
			break # Break out of loop if hit has lower score than top hit. # Should remove break as it is bad voodoo...

	return topHits
#-------------------------------------------------------------------------------------------------
# 4: Returns Accession of proteome for which the query protein is found
def getQueryProteome(proteomesCSV, queryProtein):
	
	for row in proteomesCSV:
		if row[1].strip() == queryProtein.strip(): 
			proteomeAccession = row[0]
			break # Should remove break as it is bad voodoo...
	return proteomeAccession
#===========================================================================================================
# Main program code:

# House keeping...
argsCheck() # Checks if the number of arguments are correct.

queryFile = sys.argv[1]
queryProteomesFile = sys.argv[2]

# File extension check
if not queryFile.endswith(".faa"):
	print "[Warning] " + queryFile + " may not be a amino acid fasta file!"
# File extension check
if not queryProteomesFile.endswith(".csv"):
	print "[Warning] " + queryProteomesFile + " may not be a csv file!"
	
BLASTDBFile = sys.argv[3]

print "Forward Blasting to subject proteome..."
BLASTOut = runBLAST(queryFile, BLASTDBFile) # Forward BLASTs from query protien to subject proteome
BLASTOut = filtreBLASTCSV(BLASTOut, 30) # Filtres BLAST results by PIdnet.

# Attemps to open csv contain a csv that links each query protien to its respective proteome. 
try:
	QueryProteomes = [] # Stores this CSV in memory for later use 
	inFile = open(queryProteomesFile, "r")
	reader = csv.reader(inFile) # opens file with csv module which takes into account verying csv formats and parses correctly
	for row in reader:
		QueryProteomes.append(row)
except IOError:
	print "Failed to open " + queryProteomesFile
	exit(1)

BackBlastOutput = []

print "Back-Blasting hits to query proteome..."
# For each top Hit...
for hit in BLASTOut:
	subjectProtein = hit[1]
	queryProtein = hit[0]
	
	# Takes the csv list of query proteomes and a queryProtien and find what query proteome the query protien is part of.    
	CurrentQueryProteome = getQueryProteome(QueryProteomes, queryProtein) + ".faa" 
	
	# Opens query proteome and extracts fasta formated sequence of the query protien.
	# This will be used query protien fasta will be used as a BLAST query.
	try:
		handle = open(BLASTDBFile, "rU")
		for record in SeqIO.parse(handle, "fasta") :
		    if record.id == subjectProtein: 
				subjectProtienFASTA = record.format("fasta") 
		handle.close()
	except IOError:
		print "Failed to open " + BLASTDBFile
		exit(1)
	
	BackBlastOut = runBLAST(subjectProtienFASTA, CurrentQueryProteome) #Backwards BLASTs from subject protien hit to query proteome.
	BackBlastOut = filtreBLASTCSV(BackBlastOut, 30) # Filtres BLAST results by PIdnet.
	BackHits = getTopHits(BackBlastOut) # Gets top hits from the BackBlast.	
	match = False # By default set match to false
	for BackHit in BackHits:
		if BackHit[1] == queryProtein:
			match = True
			
	if match == True:
		BackBlastOutput.append(hit)

		
OutFile = BLASTDBFile.rstrip(".faa") + ".csv"

# Attempts to write reciprocal BLAST output to file.
try:
	writeFile = open(OutFile, "w") 	
	writer = csv.writer(writeFile)
	print ">> Output file created."
	print ">> Writing Data..."
except IOError:
	print "Failed to create " + outFile
	exit(1)
	
for row in BackBlastOutput:
	writer.writerow(row)

writeFile.close()
print "done"	