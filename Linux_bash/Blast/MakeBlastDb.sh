#!/bin/bash

# Script to make a blast database from fasta file

# Time scripts starts
echo "$0"
echo "Start $(date +"%T %d-%m-%Y")."

# Check for required command line arguments
if [ $# -lt 1 ]; then
	echo "Usage: $0 <FastaFile>"
	exit 1
fi

# Define variables
FASTA="${@:1}"
echo "Working directory: $WKDIR"

# Loop throough all submitted fasta
for fasta in ${FASTA[@]}; do
	title=`basename $fasta | perl -p -e 's/.fasta$//'`
	makeblastdb -in ${fasta} -input_type 'fasta' -title ${title} -parse_seqids -dbtype 'prot'
done

# Time scripts ends
echo "Completed $(date +"%T %d-%m-%Y")."
