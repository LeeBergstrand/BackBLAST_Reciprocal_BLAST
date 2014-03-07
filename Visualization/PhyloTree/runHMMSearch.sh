#!/bin/bash
# A simple script for the batch creation of blast databases from a directory with fasta files.

for fasta;
do
	echo Running hmmsearch on $fasta
	echo $fasta
	sed "s/\.\/TestData\///" $fasta
	#hmmsearch --tblout "hmmResults$fasta.tsv" -A "hmmAlign$fasta" --acc SSU.hmm $fasta | tee "hmmHumanResults$fasta.txt"
done
echo All files searched.
exit 0

