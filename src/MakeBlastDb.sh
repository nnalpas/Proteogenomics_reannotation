#!/bin/bash

# Script to make a blast database from fasta file

# Time scripts starts
echo "$0"
echo "Start $(date +"%T %d-%m-%Y")."

# Load environment module depending on server
if [[ `hostname` == *"core-login"* ]] || [[ `hostname` == *"ifb"* ]] || [[ `hostname` == *"cpu-node"* ]]; then
	module load blast/2.12.0
elif [[ `hostname` == *"binac"* ]]; then
	module load bio/blastplus/2.11.0
else
	echo "Unknown server: $(hostname), do not know what module to load" >&2
	exit 1
fi

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
INPUTS=${@:1}

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
if [ -z "${INPUTS}" ]; then
	echo "At least one path must be specified." >&2
	display_usage
	exit 1
fi

# Loop through all submitted fasta
ls -1 ${INPUTS[@]} | 
while read file; do
	if [ ! -f ${file}.pdb ]; then
		title=`basename "$file" | perl -p -e 's/.(fasta|faa)$//'`
		add_param=""
		if [ -z "${TAXID}" ]; then
			add_param+="-taxid ${title} "
		else
			add_param+="-taxid ${TAXID} "
		fi
		makeblastdb -in "${file}" -input_type ${INPUTTYPE} -title ${title} -parse_seqids -dbtype ${DBTYPE} ${add_param}
	fi
done

# Time scripts ends
echo "Completed $(date +"%T %d-%m-%Y")."
