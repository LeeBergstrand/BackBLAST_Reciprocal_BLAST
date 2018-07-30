#!/usr/bin/env bash
# A simple script that removes multiple hits from BackBlast results.

for csv;
do
	echo Removing duplicate hits from on ${csv}
	sort -k 1,1 -t , -u ${csv} > "$csv.out"
done
echo Done.
exit 0
