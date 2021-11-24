#!/bin/bash

# Script to make a blast database from fasta file

# Time scripts starts
echo "$0"
echo "Start $(date +"%T %d-%m-%Y")."

# Function holding the usage
display_usage() { 
	echo "
		Usage: $0 [Options] Fastafiles"
	echo "
		Options:
        	-i	[str]	The input file format.
			-y	[str]	The type of the data ('nucl' or 'prot')
			-t	[int]	The taxon ID for the database(s) (if unspecified, it must be the name of the input file)
        	-h	[]	To display the help.
	"
}

# Parse user parameters
while getopts "i:y:t:h" opt; do
	case $opt in
		i)
			INPUTTYPE=$OPTARG
			echo "Working directory: $WKDIR"
			shift $((OPTIND-1)); OPTIND=1
			;;
        y)
            DBTYPE=$OPTARG
            shift $((OPTIND-1)); OPTIND=1
			;;
		t)
			TAXID=$OPTARG
			shift $((OPTIND-1)); OPTIND=1
			;;
		h)
			display_usage
			shift $((OPTIND-1)); OPTIND=1
			exit 0
			;;
		\?)
			display_usage
			exit 1
			;;
		:)
			display_usage
			exit 1
			;;
	esac
done

# Also get target folders from parameter
INPUT=${@:1}

# Check for mandatory parameters
if [ -z "${INPUTTYPE}" ]; then
	echo "The input file format must be specified." >&2
	display_usage
	exit 1
fi
if [ -z "${DBTYPE}" ]; then
	echo "The database type must be specified." >&2
	display_usage
	exit 1
fi
if [ -z "${TAXID}" ]; then
	echo "The taxon ID is not specified. Will try to parse from file name." >&2
fi
if [ -z "${INPUT}" ]; then
	echo "At least one path must be specified." >&2
	display_usage
	exit 1
fi

# Loop through all submitted fasta
ls -1 ${INPUTS[@]} | 
while read file; do
	title=`basename "$file" | perl -p -e 's/.(fasta|faa)$//'`
	if [ -z "${TAXID}" ]; then
		TAXID=${title}
	fi
	makeblastdb -in "${file}" -input_type ${INPUTTYPE} -title ${title} -parse_seqids -dbtype ${DBTYPE} -taxid ${TAXID}
done

# Time scripts ends
echo "Completed $(date +"%T %d-%m-%Y")."
