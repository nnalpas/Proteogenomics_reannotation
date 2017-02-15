#!/bin/bash

# Script to get the ORF from sequence file

# Time scripts starts
echo "$0"
echo "Start $(date +"%T %d-%m-%Y")."

# Check for required command line arguments
if [ $# -lt 6 ]; then
	echo "Usage: $0 <FastaFile>"
	exit 1
fi

# Define variables
WKDIR=$1
echo "Working directory: $WKDIR"
FIND=$2
TABLE=$3
MINSIZE=$4
CIRCULAR=$5
SEQ="${@:6}"


# Create the output folder
if [ ! -d ${WKDIR} ] ; then
	mkdir ${WKDIR}
fi

# Loop through all submitted fasta
for files in ${SEQ[@]}; do

	outfile=`basename $files`
	$HOME/Software/EMBOSS/bin/getorf -sequence ${files} -outseq ${WKDIR}/Find${FIND}_${outfile} -find=${FIND} -table=${TABLE} -minsize=${MINSIZE} -circular=${CIRCULAR}
	
	# Remove all entries that are empty
	${HOME}/bin/Remove_empty_entries.sh ${WKDIR}/Find${FIND}_${outfile}
	
done

# Time scripts ends
echo "Completed $(date +"%T %d-%m-%Y")."
