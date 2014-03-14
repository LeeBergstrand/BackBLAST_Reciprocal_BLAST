#!/usr/bin/env python 
# Created by: Lee Bergstrand
# Descript: A Bio-Python program that takes a list of query proteins and uses local BLASTp to search
#           for highly similer proteins within a local blast database (usally a local db of a target 
#           proteome). The program then BLASTps backward from the found subject proteins to the proteome 
#           for which the original query protein is found in order to confirm gene orthology. 
#             
# Requirements: - This program requires the Biopython module: http://biopython.org/wiki/Download
#               - This script requires BLAST+ 2.2.9 or later.
#               - All operations are done with protien sequences.
#               - All query proteins should be from sequenced genomes in order to facilitate backwards BLAST. 
#               - MakeBlastDB must be used to create BLASTp databases for both query and subject proteomes.
#               - BLAST databases require the FASTA file they were made from be in the same directory.
#  
# Usage: BackBLAST.py <queryGeneList.faa> <queryProteomes.csv> <subject1.faa> 
# Example: BackBLAST.py queryGeneList.faa queryProteomes.csv AUUJ00000000.faa
#----------------------------------------------------------------------------------------
#===========================================================================================================
#Imports & Setup:
import sys
import csv
import subprocess
from Bio import SeqIO
from Graph import Vertex
from Graph import Graph
from multiprocessing import cpu_count

processors = cpu_count() # Gets number of processor cores for BLAST.

# Dev Imports:
import time # For profiling purposes.
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
	# Runs BLASTp and saves the output to a string. Blastp is set to output a csv which can be parsed by Pythons CSV module.
	BLASTOut = subprocess.check_output(["blastp", "-db", BLASTDBFile, "-query", query, "-evalue", "1e-30", "-num_threads", str(processors), "-outfmt", "10 qseqid sseqid pident evalue qcovhsp score"]) 
	return BLASTOut
#-------------------------------------------------------------------------------------------------
# 3: Filters HSPs by Percent Identity...
def filtreBLASTCSV(BLASTOut):
	
	minIdent = 25
	
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
# 4: Returns a list of accession of the proteomes for which the query proteins is are found.
def GetQueryProteomeAccessions(queryProteomesFile):
	# Reads sequence file list and stores it as a string object. Safely closes file.try:
	try:	
		with open(queryProteomesFile,"r") as newFile:
			proteomeAccessions = newFile.read()
			newFile.close()
	except IOError:
		print "Failed to open " + queryProteomesFile
		exit(1)

	QueryProteomeAccessions = proteomeAccessions.splitlines() # Splits string into a list. Each element is a single line from the string.	
	return QueryProteomeAccessions
#----------------------------------------------------------------------
# 5: Creates a python dictionary (hash table) that contains the the fasta for each protien in the proteome.
def createProteomeHash(ProteomeFile):
	ProteomeHash = dict() 
	try:
		handle = open(ProteomeFile, "rU")
		proteome = SeqIO.parse(handle, "fasta")
		for record in proteome:
			ProteomeHash.update({ record.id : record.format("fasta") })
		handle.close()
	except IOError:
		print "Failed to open " + ProteomeFile
		exit(1)
		
	return ProteomeHash
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
if not queryProteomesFile.endswith(".txt"):
	print "[Warning] " + queryProteomesFile + " may not be a txt file!"
	
BLASTDBFile = sys.argv[3]
print "Opening " + BLASTDBFile + "..."

BLASTGraph = Graph() # Creates graph to map BLAST hits.

print ">> Forward Blasting to subject proteome..."
BLASTForward = runBLAST(queryFile, BLASTDBFile) # Forward BLASTs from query protiens to subject proteome
BLASTForward = filtreBLASTCSV(BLASTForward) # Filtres BLAST results by PIdnet.

if len(BLASTForward) == 0:
	print ">> No Forward hits in subject proteome were found."
	print ">> Exiting."
	exit(1) # Aborts program. (exit(1) indicates that an error occured)

SubjectProteomeHash = createProteomeHash(BLASTDBFile) # Creates python dictionary contianing every protien in the subject Proteome.
BackBlastQueryFASTAs = []

print ">> Creating Back-Blasting Query from found subject protiens..."
# For each top Hit...
for hit in BLASTForward:
	subjectProtein = hit[1]
	queryProtein = hit[0]
	subjectProtienFASTA = SubjectProteomeHash.get(subjectProtein) # Extracts subjectProtien from python dictionary.
	subjectProtienFASTA.strip()
	BackBlastQueryFASTAs.append(subjectProtienFASTA) # Addes to master protien list.
	
CompleteBackBlastQuery = "".join(BackBlastQueryFASTAs)

try:
	writeFile = open("tempQuery.faa", "w") 
	writeFile.write(CompleteBackBlastQuery) 
	writeFile.close()
except IOError:
	print "Failed to create " + "tempQuery.faa"
	exit(1)

BackBLASTResults = []
print ">> Back-Blasting hits to query proteomes..."
for proteome in GetQueryProteomeAccessions(queryProteomesFile):
	proteomeFile = proteome + ".faa"
	BackBLASTResults.append(runBLAST("tempQuery.faa", proteomeFile))

BLASTBackward = "".join(BackBLASTResults)
BLASTBackward = filtreBLASTCSV(BLASTBackward) # Filtres BLAST results by PIdnet.

print ">> Creating Graph..."
for hit in BLASTForward:
	BLASTGraph.addEdge(hit[0],hit[1],hit[5])
for hit in BLASTBackward:
	BLASTGraph.addEdge(hit[0],hit[1],hit[5])
	
BackBlastOutput = list(BLASTForward)

print ">> Checking if forward hit subjects have better reciprocal hits than query."
for hit in BLASTForward: 
	queryProtien = BLASTGraph.getVertex(hit[0])
	subjectProtien = BLASTGraph.getVertex(hit[1])
	
	topBackHitScore = 0
	# Find the top score of the best reciprical BLAST hit.
	for backHit in subjectProtien.getConnections():
		backHitScore = subjectProtien.getWeight(backHit) # The edge weight between the subject and its reciprocal BLAST hit is the BLAST score.
		if backHitScore >= topBackHitScore:
			topBackHitScore = backHitScore
	
	# Check if the query is the best reciprical BLAST hit for the subject.
	deleteHit = False
	if queryProtien in subjectProtien.getConnections():
		BackHitToQueryScore = subjectProtien.getWeight(queryProtien) # The edge weight between the subject and the query is the reciprocal BLAST score.
		if BackHitToQueryScore < topBackHitScore: 
			deleteHit = True  # If the query is not the best reciprocal BLAST hit simply delete it from the BackBlastOutput.
	else:
		deleteHit = True # If the query is not a reciprocal BLAST hit simply delete it from the BackBlastOutput.
	
	if deleteHit == True:
		del BackBlastOutput[BackBlastOutput.index(hit)] # Delete the forward BLAST hit from BackBlastOutput.
print ">> Done"

OutFile = BLASTDBFile.rstrip(".faa") + ".csv" 

#Attempts to write reciprocal BLAST output to file.
try:
	writeFile = open(OutFile, "w") 	
	writer = csv.writer(writeFile)
	print ">> Output file created."
	print ">> Writing Data..."
	for row in BackBlastOutput:
		writer.writerow(row)
except IOError:
	print "Failed to create " + outFile
	exit(1)
print "done"	