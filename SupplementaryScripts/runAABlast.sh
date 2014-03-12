#!/bin/bash
# A simple script for the batch creation of blast databases from a directory with fasta files.

for fasta;
do
	echo Running blast on $fasta
	blastp -db $fasta -query SteriodQueryV3.faa -out "$fasta.csv" -evalue 1e-40 -num_threads 16 -outfmt "10 qseqid sseqid pident evalue qcovhsp score" 
done
echo All files BLASTed.
exit 0

