#!/bin/bash

# Script to make a blast database from fasta file

# Time scripts starts
echo "$0"
echo "Start $(date +"%T %d-%m-%Y")."

# Check for required command line arguments
if [ $# -lt 3 ]; then
	echo "Usage: $0 <FastaFile>"
	exit 1
fi

# Define variables
INPUTTYPE=$1
DBTYPE=$2
FASTA="${@:3}"
echo "Working directory: $WKDIR"

# Loop through all submitted fasta
for fasta in ${FASTA[@]}; do
	title=`basename $fasta | perl -p -e 's/.fasta$//'`
	makeblastdb -in ${fasta} -input_type ${INPUTTYPE} -title ${title} -parse_seqids -dbtype ${DBTYPE}
done

# Time scripts ends
echo "Completed $(date +"%T %d-%m-%Y")."
