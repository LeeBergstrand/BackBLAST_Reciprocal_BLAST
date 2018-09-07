#!/usr/bin/env bash
set -euo pipefail
# Copyright Lee H. Bergstrand and Jackson M. Tsuji, 2018
# A simple script that removes multiple hits from BackBLAST results.

script_name=${0##*/}

# If input field is empty, give error message and end script
if [ $# == 0 ]; then
    printf "${script_name}: A simple script that removes multiple hits from BackBLAST results.\n"
    printf "Usage: ${script_name} input.csv > output.csv\n\n"
    exit 1
fi

# Get input from user
input=$1

# Note: (>&2 ... ) prints to STDERR
(>&2 echo "[$(date -u)] Removing duplicate hits from ${input}")

# Prints output to STDOUT
sort -k 1,1 -t , -u ${input}

(>&2 echo "[$(date -u)] ${script_name}: done.")

exit 0
