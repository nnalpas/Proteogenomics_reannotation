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
	#mkdir -p "$WKDIR/${array[0]}/Genome"; mkdir -p "$WKDIR/${array[0]}/MQ_6frame"; mkdir -p "$WKDIR/${array[0]}/Phenodata"
	#cp "$WKDIR/Genome/${array[5]}" "$WKDIR/${array[0]}/Genome/"
	find $WKDIR/Genome/ -name "*FIXED.fasta" -exec cp '{}' $WKDIR/${array[0]}/Genome/ \;

	# Edit the parameter file
	NEWPARAM=`basename $PARAM | sed "s/.txt/_${array[0]}.txt/"`
	cp -f "$PARAM" "$WKDIR/${array[0]}/Phenodata/$NEWPARAM"
	sed -i -e "s/PROJECT_REPLACE/${array[0]}/g" "$WKDIR/${array[0]}/Phenodata/$NEWPARAM"
	sed -i -e "s/GENOME_REPLACE/${array[5]}/g" "$WKDIR/${array[0]}/Phenodata/$NEWPARAM"

	# Source the variable from the parameter file
	ProjDir="$WKDIR/${array[0]}"
	source "$WKDIR/${array[0]}/Phenodata/$NEWPARAM"

	# Edit the mqpar file
	#cp "$WKDIR/MQ_6frame/${array[4]}" "$WKDIR/${array[0]}/MQ_6frame/"; cp "$WKDIR/MQ_6frame/mqpar_posix.xml" "$WKDIR/${array[0]}/MQ_6frame/";
	SIXFRAMEPROT=`basename $GENOME | perl -p -e "s%^%${WKDIR}\\/${array[0]}\\/Nuc_translation\\/Find0_%" | perl -p -e "s/\\.fasta/_FIXED.fasta/"`
	#sed -i -e "s%ORF_REPLACE%$SIXFRAMEPROT%g" "$WKDIR/${array[0]}/MQ_6frame/mqpar_posix.xml"
	#sed -i -e "s%REF_REPLACE%$UNIREFPROT%g" "$WKDIR/${array[0]}/MQ_6frame/mqpar_posix.xml"
	#sed -i -e "s%RAW_REPLACE%$WKDIR\\/${array[0]}\\/MQ_6frame\\/${array[4]}%g" "$WKDIR/${array[0]}/MQ_6frame/mqpar_posix.xml"
	#sed -i -e "s%PROJECT_REPLACE%${array[0]}%g" "$WKDIR/${array[0]}/MQ_6frame/mqpar_posix.xml"

	# Edit the serial blast file
	cp -f "$WKDIR/Phenodata/Blast_iterations.txt" "$WKDIR/${array[0]}/Phenodata/Blast_iterations.txt"
	sed -i -e "s%REF_REPLACE%$UNIREFPROT%g" "$WKDIR/${array[0]}/Phenodata/Blast_iterations.txt"
	sed -i -e "s%GENOME_REPLACE%$GENOME%g" "$WKDIR/${array[0]}/Phenodata/Blast_iterations.txt"
	sed -i -e "s%CDS_REPLACE%$UNIREFGENE%g" "$WKDIR/${array[0]}/Phenodata/Blast_iterations.txt"

	qsub -d ${HOME}/ws/tu_kxmna01-Proteogenomics-0 -m abe -M nicolas.nalpas@ifiz.uni-tuebingen.de -j eo -V -v SCRIPT_FLAGS="$WKDIR/${array[0]}/Phenodata/$NEWPARAM" -q short -l nodes=1:ppn=14,walltime=48:00:00 -N blast ${HOME}/bin/SixFrame_proteogenomics.sh
	
done < ${MULTI}

# Time scripts ends
echo "Completed $(date +"%T %d-%m-%Y")."
