#!/bin/bash

# Script to make a blast database from fasta file

# Time scripts starts
echo "$0"
echo "Start $(date +"%T %d-%m-%Y")."

# Check for required command line arguments
if [ $# -lt 4 ]; then
	echo "Usage: $0 <InputType> <DbType> <TaxId> <FastaFiles>"
	exit 1
fi

# Define variables
INPUTTYPE=$1
DBTYPE=$2
TAXID=$3
FASTA="${@:4}"

# Loop through all submitted fasta
for fasta in ${FASTA[@]}; do
	title=`basename $fasta | perl -p -e 's/.fasta$//'`
	makeblastdb -in ${fasta} -input_type ${INPUTTYPE} -title ${title} -parse_seqids -dbtype ${DBTYPE} -taxid ${TAXID}
done

# Time scripts ends
echo "Completed $(date +"%T %d-%m-%Y")."
