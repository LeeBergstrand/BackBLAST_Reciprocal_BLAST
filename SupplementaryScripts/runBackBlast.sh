#!/usr/bin/env bash
# A simple script for the batch creation of blast databases from a directory with fasta files.

timestamp=$(date +"%T")

for fasta;
do
	echo Running Back-BLAST on ${fasta}
	python -u ../BackBLAST.py ../ExampleData/ExampleQueries/MycobacteriumCholesterolCluster.faa ../ExampleData/ExampleQueryProteomes/AL123456.3.faa ${fasta} 2>&1 | tee -a "BackLog$timestamp.txt"
	date +"%T" >> "BackLog$timestamp.txt" 
done
echo All files BLASTed.
exit 0

