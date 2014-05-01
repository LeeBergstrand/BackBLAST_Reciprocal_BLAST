#!/bin/bash
# A simple script that extracts BLAST hit subjects for a certain set of query accessions from from BackBlast results.

for file;
do
	echo Etracting subjects from $file
	python /Users/lee/Dropbox/RandD/Repositories/BackBLAST-Gene-Cluster-Finder/Visualization/HeatMap/ExtractBackBlastSubject.py $file /Users/lee/Dropbox/RandD/Repositories/BackBLAST-Gene-Cluster-Finder/Visualization/HeatMap/TestData/SecondRun/CompleteResults.csv
	grep "_" "$file.out" > "$file.prod"
	grep -v "_" "$file.out" > "$file.ncbi"
done
echo Done.
exit 0


