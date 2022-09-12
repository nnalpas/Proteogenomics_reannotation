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
	NEWPARAM=`basename $PARAM | sed "s/.txt/_${array[0]}.txt/"`
	mkdir -p "$WKDIR/${array[0]}/Genome"; mkdir -p "$WKDIR/${array[0]}/MaxQuant"; mkdir -p "$WKDIR/${array[0]}/Phenodata"
	cp "$WKDIR/Genome/${array[5]}" "$WKDIR/${array[0]}/Genome/"; cp "$PARAM" "$WKDIR/${array[0]}/Phenodata/$NEWPARAM"
	find $WKDIR/Genome/ -name "*FIXED.fasta" -exec cp '{}' $WKDIR/${array[0]}/Genome/ \;

	# Edit the parameter file
	sed -i -e "s/PROJECT_REPLACE/${array[0]}/g" "$WKDIR/${array[0]}/Phenodata/$NEWPARAM"
	sed -i -e "s/GENOME_REPLACE/${array[5]}/g" "$WKDIR/${array[0]}/Phenodata/$NEWPARAM"

	# Source the variable from the parameter file
	ProjDir="$WKDIR/${array[0]}"
	source "$WKDIR/${array[0]}/Phenodata/$NEWPARAM"
	echo "Genome: $GENOME"
	echo "Proteome: $UNIREFPROT"

	# Edit the mqpar file
	cp "$WKDIR/MaxQuant/${array[4]}" "$WKDIR/${array[0]}/MaxQuant/"; cp "$WKDIR/MaxQuant/mqpar_posix.xml" "$WKDIR/${array[0]}/MaxQuant/";
	SIXFRAMEPROT=`basename $GENOME | perl -p -e "s%^%${WKDIR}\\/${array[0]}\\/Nuc_translation\\/Find0_%" | perl -p -e "s/\\.fasta/_FIXED.fasta/"`
	echo "ORF: $SIXFRAMEPROT"
	sed -i -e "s/ORF_REPLACE/$SIXFRAMEPROT/g" "$WKDIR/${array[0]}/MaxQuant/mqpar_posix.xml"
	sed -i -e "s/REF_REPLACE/$UNIREFPROT/g" "$WKDIR/${array[0]}/MaxQuant/mqpar_posix.xml"
	sed -i -e "s/RAW_REPLACE/$WKDIR\\/${array[0]}\\/MaxQuant\\/${array[4]}/g" "$WKDIR/${array[0]}/MaxQuant/mqpar_posix.xml"
	sed -i -e "s/PROJECT_REPLACE/${array[0]}/g" "$WKDIR/${array[0]}/MaxQuant/mqpar_posix.xml"

	#qsub -d ${HOME}/ws/tu_kxmna01-Proteogenomics-0 -m abe -M nicolas.nalpas@ifiz.uni-tuebingen.de -j eo -V -v SCRIPT_FLAGS="$WKDIR/${array[0]}/Phenodata/$NEWPARAM" -q short -l nodes=1:ppn=1,walltime=24:00:00 -N maxquant ${HOME}/bin/SixFrame_proteogenomics.sh
	
done < ${MULTI}

# Time scripts ends
echo "Completed $(date +"%T %d-%m-%Y")."
