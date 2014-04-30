#!/bin/bash
# A simple script that removes multiple hits from BackBlast results.

for file;
do
	echo Etracting subjects from $file
	python ExtractBackBlastSubject.py $file /Users/lee/Dropbox/RandD/Repositories/BackBLAST-Gene-Cluster-Finder/Visualization/HeatMap/TestData/SecondRun/CompleteResults.csv
	grep "_" "$file.out" > "$file.prod"
	grep -v "_" "$file.out" > "$file.ncbi"
done
echo Done.
exit 0


