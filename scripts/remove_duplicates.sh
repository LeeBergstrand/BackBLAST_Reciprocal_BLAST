#!/usr/bin/env bash
set -euo pipefail
# remove_duplicates.sh
# Copyright Lee H. Bergstrand and Jackson M. Tsuji, 2025
# A simple script that removes multiple hits from BackBLAST results
# Part of the BackBLAST pipeline

SCRIPT_NAME=${0##*/}

# Print help statement
if [[ $# -ne 1 ]]; then
  echo "Incorrect number of arguments provided. Please run '-h' or '--help' to see help. Exiting..." >&2
  exit 1
elif [[ $1 = "-h" ]] || [[ $1 = "--help" ]]; then
  printf "${SCRIPT_NAME}: A simple script that removes multiple hits from BackBLAST results.\n"
  printf "Copyright Lee H. Bergstrand and Jackson M. Tsuji, 2025\n\n"
  printf "Usage: ${SCRIPT_NAME} input.csv > output.csv\n\n"
  printf "Note: log information is printed to STDERR.\n\n"
  exit 0
fi

# Print stack trace to log
set -x

# Get input from user
input=$1

echo "[ $(date -u) ]: Removing duplicate hits from ${input}" >&2

# Prints output to STDOUT
sort -k 1,1 -t , -u ${input}

echo "[ $(date -u) ]: ${SCRIPT_NAME}: done." >&2
