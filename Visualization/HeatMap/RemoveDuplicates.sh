#!/bin/bash
# A simple script for cleans duplicate sequences  batch creation of blast databases from a directory with fasta files.

for fasta;
do
	echo Removing duplicate hits from on $fasta
	sort -k 1,1 -t , -u $fasta > "$fasta.out"
done
echo Done.
exit 0


