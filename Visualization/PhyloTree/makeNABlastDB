#!/bin/bash
# A simple script for the batch creation of blast databases from a directory with fasta files.

for fasta;
do
	echo Makeing blast database for $fasta
	makeblastdb -in $fasta -input_type fasta -dbtype nucl -parse_seqids;
done
echo All databases created.
exit 0
