#!/bin/bash
# A simple script that removes multiple hits from BackBlast results.

for fasta;
do
	echo Removing duplicate hits from on ${fasta}
	sort -k 1,1 -t , -u ${fasta} > "$fasta.out"
done
echo Done.
exit 0


