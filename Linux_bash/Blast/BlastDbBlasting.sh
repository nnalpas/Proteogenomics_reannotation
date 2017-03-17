#!/bin/bash

# Script to blast two blast databases against each other

# Time scripts starts
echo "$0"
echo "Start $(date +"%T %d-%m-%Y")."

# Check for required command line arguments
if [ "$#" -eq 9 ]; then
	echo "Usage: $0 <WorkDir> <DbType> <Entry> <Task> <InputDb> <Db> <Eval> <NumAlign> <Threads> <OutBlast>"
	exit 1
fi

# Define variables
WKDIR=$1
echo "Working directory: $WKDIR"
DBTYPE=$2
ENTRY=$3
TASK=$4
INPUTSEQ=$5
DB=$6
EVAL=$7
NUMALIGN=$8
THREADS=$9; shift
OUTBLAST=$9

# Create the output folder
if [ ! -d ${WKDIR} ] ; then
	mkdir ${WKDIR}
fi

if [ ${ENTRY} == 'all' ]; then

	# All entries retrieval and blasting of retrieved entries against another database
	blastdbcmd -db ${INPUTSEQ} \
	  -dbtype ${DBTYPE} \
	  -entry ${ENTRY} | \
	  eval "${TASK} \
	  -query - \
	  -task ${TASK} \
	  -db ${DB} \
	  -out ${WKDIR}/${OUTBLAST} \
	  -evalue ${EVAL} \
	  -num_alignments ${NUMALIGN} \
	  -num_threads ${THREADS} \
	  -outfmt '6 qseqid sseqid pident nident mismatch length gapopen qstart qend sstart send evalue bitscore score'"

else

	# Specific entries retrieval and blasting of retrieved entries against another database
	blastdbcmd -db ${INPUTSEQ} \
	  -dbtype ${DBTYPE} \
	  -entry_batch ${ENTRY} | \
	  eval "${TASK} \
	  -query - \
	  -task ${TASK} \
	  -db ${DB} \
	  -out ${WKDIR}/${OUTBLAST} \
	  -evalue ${EVAL} \
	  -num_alignments ${NUMALIGN} \
	  -num_threads ${THREADS} \
	  -outfmt '6 qseqid sseqid pident nident mismatch length gapopen qstart qend sstart send evalue bitscore score'"

fi

# Time scripts ends
echo "Completed $(date +"%T %d-%m-%Y")."
