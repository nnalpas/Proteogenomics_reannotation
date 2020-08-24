#!/bin/bash

# Script to find best blast for multiple input files

# Time scripts starts
echo "$0"
echo "Start $(date +"%T %d-%m-%Y")."

# Function holding the usage
display_usage() { 
	echo "
		Usage: $0 [Options] <Blast files>"
	echo "
		Options:
        	-o	[str]	The output directory.
			-f	[str]	One or more blast fields filtering in R language (e.g. evalue <= 0.0001 & score >= 100)
			-m	[str]	What to do for multi-hits entries (i.e. 'keep' (default), 'remove', 'uniquify')
        	-h	[]	To display the help.
	"
}

# Define default values
o_default="`pwd`/Blast"
f_default="''"
m_default="keep"

# Parse user parameters
while getopts "o:f:m:h" opt; do
	case $opt in
		o)
			WKDIR=$OPTARG
			echo "Working directory: $WKDIR"
			shift $((OPTIND-1)); OPTIND=1
			;;
		f)
            FILTER=$OPTARG
            shift $((OPTIND-1)); OPTIND=1
			;;
		m)
            MULTI=$OPTARG
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

# Get all blast result files
TARGETFILES=${@:1}

# Setting default values ${VARIABLE=DEFAULT_VALUE}
# This command will set VARIABLE to DEFAULT_VALUE if it is currently 
# undefined, then return VARIABLE
: ${WKDIR=$o_default}
: ${FILTER=$f_default}
: ${MULTI=$m_default}

# Check for mandatory parameters
if [ -z "${TARGETFILES}" ]; then
	echo "At least one blast file must be specified." >&2
	display_usage
	exit 1
fi

# Create the output folder
if [ ! -d ${WKDIR} ] ; then
	mkdir ${WKDIR}
fi

# Loop through all blast files
for file in `ls ${TARGETFILES}`; do
	Best_blast.R \
		-i ${file} \
		-f ${FILTER} \
		-m ${MULTI} \
		-o ${WKDIR} 2>&1
done

# Time scripts ends
echo "Completed $(date +"%T %d-%m-%Y")."
