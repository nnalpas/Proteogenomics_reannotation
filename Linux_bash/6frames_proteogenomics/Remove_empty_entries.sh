#!/bin/bash

# Script to clean up fasta file (from nucleotide translation) by removing empty entries

# Time scripts starts
echo "$0"
echo "Start $(date +"%T %d-%m-%Y")."

# Check for required command line arguments
if [ $# -lt 1 ]; then
	echo "Usage: $0 <FastaFile>"
	exit 1
fi

# Define variables
SEQ="${@:1}"

# Use awk to check consecutive line in fasta files for header symbol, and remove header line if not followed by actual sequence
for files in ${SEQ[@]}; do
	out=`echo $files | perl -p -e 's/^(.+)(\\..+?)$/$1_FIXED$2/'`
	awk '/^>/,/^>/ {lastheader = $0; next}
	{ if (lastheader != "") { print lastheader; lastheader = ""; }
	print }' $files > $out
done


