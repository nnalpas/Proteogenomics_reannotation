#!/bin/bash

# Script to find reciprocal best blast for multiple input files

# Time scripts starts
echo "$0"
echo "$@"
echo "Start $(date +"%T %d-%m-%Y")."

# Function holding the usage
display_usage() { 
	echo "
		Usage: $0 [Options] -b <Blast files> -r <Reciprocal blast folder>"
	echo "
		Options:
			-o	[str]	The output directory.
			-b	[str]	One or more blast results files
			-r	[str]	One folder containing the reciprocal blast results
			-f	[str]	One or more blast fields filtering in R language (e.g. evalue <= 0.0001 & score >= 100)
			-m	[str]	What to do for multi-hits entries (i.e. 'keep' (default), 'remove', 'uniquify')
        	-h	[]	To display the help.
	"
}

# Define default values
o_default="`pwd`/ReciprocalBlast"
m_default="keep"

# Parse user parameters
while getopts "o:b:r:f:m:h" opt; do
	case $opt in
		o)
			WKDIR=$OPTARG
			echo "Working directory: $WKDIR"
			shift $((OPTIND-1)); OPTIND=1
			;;
		b)
            BLASTS=$OPTARG
            shift $((OPTIND-1)); OPTIND=1
			;;
		r)
            RECIPROCALS=$OPTARG
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

# Setting default values ${VARIABLE=DEFAULT_VALUE}
# This command will set VARIABLE to DEFAULT_VALUE if it is currently 
# undefined, then return VARIABLE
: ${WKDIR=$o_default}
: ${MULTI=$m_default}

# Check for mandatory parameters
if [ -z "${BLASTS}" ]; then
	echo "At least one blast file must be specified." >&2
	display_usage
	exit 1
fi
if [ -z "${RECIPROCALS}" ]; then
	echo "A reciprocal blast folder must be specified." >&2
	display_usage
	exit 1
fi
FILTER_ARG=""
if [ -n "${FILTER}" ]; then
	FILTER_ARG="-f ${FILTER}"
fi

# Create the output folder
if [ ! -d ${WKDIR} ] ; then
	mkdir ${WKDIR}
fi

# Loop through all blast files
for file in `ls ${BLASTS}`; do
	rec_pat=`basename $file`
	rec_file=`find ${RECIPROCALS} -name "${rec_pat}*"`
	Reciprocal_best_blast.R \
		-b ${file} \
		-r ${rec_file} ${FILTER_ARG} \
		-m ${MULTI} \
		-o ${WKDIR} 2>&1
done

# Time scripts ends
echo "Completed $(date +"%T %d-%m-%Y")."
