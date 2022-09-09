#!/bin/bash

# Script to automate change of the different parameters and input files

# Time scripts starts
echo "$0"
echo "Start $(date +"%T %d-%m-%Y")."

# Check for required command line arguments
if [ $# -ne 3 ]; then
	echo "Usage: $0 <MultiRunFile>"
	exit 1
fi

# Define variables
WKDIR="$1"
MULTI="$2"
PARAM="$3"

# Read and parse the multi genome file
while IFS= read -r line; do

	IFS=$'\t' read -r -a array <<< "$line"

	# Create new folders and copy working directory data
	mkdir -p "$WKDIR/${array[0]}/Genome"; mkdir -p "$WKDIR/${array[0]}/MaxQuant"; mkdir -p "$WKDIR/${array[0]}/Phenodata"
	cp "$WKDIR/Genome/${array[5]}" "$WKDIR/${array[0]}/Genome/"; cp "$PARAM" "$WKDIR/${array[0]}/Phenodata/"

	# Edit the parameter file
	NEWPARAM=`basename $PARAM`
	sed -ie "s/PROJECT_REPLACE/${array[0]}/g" "$WKDIR/${array[0]}/Phenodata/$NEWPARAM"
	sed -ie "s/GENOME_REPLACE/${array[5]}/g" "$WKDIR/${array[0]}/Phenodata/$NEWPARAM"

done < ${MULTI}

# Time scripts ends
echo "Completed $(date +"%T %d-%m-%Y")."
