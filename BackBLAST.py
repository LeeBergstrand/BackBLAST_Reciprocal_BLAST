#!/usr/bin/env python 
# Created by: Lee Bergstrand
# Descript: A Biopython program that takes a list of query proteins and uses local BLASTp to search 
#	    for highly similer proteins within a local blast database (usally a local db of a target 
#	    proteome). The program then BLASTps backwards from the found subject proteins to the query 
#	    proteome to confirm gene orthology.
#             
# Requirements: - This program requires the Biopython module: http://biopython.org/wiki/Download
#               - This script requires BLAST+ 2.2.9 or later.
#               - All operations are done with protien sequences.
#               - All query proteins should be from sequenced genomes in order to facilitate backwards BLAST. 
#               - MakeAABlastDB must be used to create BLASTn databases for both query and subject proteomes.
#               - BLAST databases require that the FASTA file they were made from remain in the same directory.
#  
# Usage: BackBLAST.py <queryGeneList.faa> <queryBLASTDB.faa> <subjectBLASTDB.faa> 
# Example: BackBLAST.py queryGeneList.faa AL123456.3.faa AUUJ00000000.faa
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
def argsCheck(argsCount):
	if len(sys.argv) < argsCount:
		print "Orthologous Gene Finder"
		print "By Lee Bergstrand\n"
		print "Please refer to source code for documentation\n"
		print "Usage: " + sys.argv[0] + " <queryGeneList.faa> <queryBLASTDB.faa> <subjectBLASTDB.faa> \n"
		print "Examples:" + sys.argv[0] + " queryGeneList.faa AL123456.3.faa AUUJ00000000.faa"
		sys.exit(1) # Aborts program. (exit(1) indicates that an error occured)
#-------------------------------------------------------------------------------------------------
# 2: Runs BLAST, can either be sent a fasta formatted string or a file ...
def runBLAST(query, BLASTDBFile):
	# Runs BLASTp and saves the output to a string. Blastp is set to output a csv which can be parsed by Pythons CSV module.
	BLASTOut = subprocess.check_output(["blastp", "-db", BLASTDBFile, "-query", query, "-evalue", "1e-25", "-num_threads", str(processors), "-outfmt", "10 qseqid sseqid pident evalue qcovhsp bitscore"]) 
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
			HSP[4] = float(HSP[4])
			HSP[5] = float(HSP[5]) 
			BLASTCSVOutFiltred.append(HSP) # Appends to output array.
	
	return BLASTCSVOutFiltred
#-------------------------------------------------------------------------------------------------
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
		sys.exit(1)
		
	return ProteomeHash
#===========================================================================================================
# Main program code:
# House keeping...
argsCheck(4) # Checks if the number of arguments are correct.

queryFile   = sys.argv[1]
queryBLASTDBFile = sys.argv[2]
subjectBLASTDBFile = sys.argv[3]

print "Opening " + subjectBLASTDBFile + "..."

# File extension checks
if not queryFile.endswith(".faa"):
	print "[Warning] " + queryFile + " may not be a amino acid fasta file!"
if not queryBLASTDBFile.endswith(".faa"):
	print "[Warning] " + queryBLASTDBFile + " may not be a amino acid fasta file!"
if not subjectBLASTDBFile.endswith(".faa"):
	print "[Warning] " + subjectBLASTDBFile + " may not be a amino acid fasta file!"

OutFile = subjectBLASTDBFile.rstrip(".faa") + ".csv" 

BLASTGraph = Graph() # Creates graph to map BLAST hits.

print ">> Forward Blasting to subject proteome..."
BLASTForward = runBLAST(queryFile, subjectBLASTDBFile) # Forward BLASTs from query protiens to subject proteome
BLASTForward = filtreBLASTCSV(BLASTForward) # Filtres BLAST results by PIdnet.

if len(BLASTForward) == 0:
	print ">> No Forward hits in subject proteome were found."
	# Writes empty file for easier data processing.
	try:
		writeFile = open(OutFile, "w") 	
		writer = csv.writer(writeFile)
	except IOError:
		print ">> Failed to create " + outFile
		sys.exit(1)
	print ">> Exiting.\n\n"
	sys.exit(0) # Aborts program. (exit(0) indicates that no error occured)

SubjectProteomeHash = createProteomeHash(subjectBLASTDBFile) # Creates python dictionary contianing every protien in the subject Proteome.
BackBlastQueryFASTAs = []

print ">> Creating Back-Blasting Query from found subject protiens..."
# For each top Hit...
for hit in BLASTForward:
	subjectProtein = hit[1]
	queryProtein = hit[0]
	subjectProtienFASTA = SubjectProteomeHash.get(subjectProtein) # Extracts subjectProtien from python dictionary.
	subjectProtienFASTA.strip()
	BackBlastQueryFASTAs.append(subjectProtienFASTA) # Addes current subject to overall protien list.
	
CompleteBackBlastQuery = "".join(BackBlastQueryFASTAs)

# Attempt to write a temporary fasta file for the reverse BLAST to use. 
try:
	writeFile = open("tempQuery.faa", "w") 
	writeFile.write(CompleteBackBlastQuery) 
	writeFile.close()
except IOError:
	print "Failed to create tempQuery.faa"
	sys.exit(1)

print ">> Blasting backwards from subject genome to query genome."
# Run backwards BLAST towards query proteome.
BLASTBackward = runBLAST("tempQuery.faa", queryBLASTDBFile)
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

#Attempts to write reciprocal BLAST output to file.
try:
	writeFile = open(OutFile, "w") 	
	writer = csv.writer(writeFile)
	print ">> Output file created."
	print ">> Writing Data..."
	for row in BackBlastOutput:
		writer.writerow(row)
except IOError:
	print ">> Failed to create " + outFile
	sys.exit(1)
print ">> Done\n"	
