#!/bin/bash

# Script to determine location of potential signal peptide via Signalp.

# Time scripts starts
echo "$0"
echo "Start $(date +"%T %d-%m-%Y")."

# Function holding the usage
display_usage() { 
	echo "
		Usage: $0 [Options] <Input FASTA 1> <...> <FASTA N>"
	echo "
		Options:
			-o	[str]	The output directory.
			-x	[str]	The organism of interest; i.e. Archaea: 'arch', Gram-positive: 'gram+', Gram-negative: 'gram-' or Eukarya: 'euk'
			-t	[int]	The number of threads.
			-h	[]	To display the help.
	"
}

# Define default values
o_default="`pwd`/Signalp"
x_default="euk"
t_default=1

# Parse user parameters
while getopts "o:x:t:h" opt; do
	case $opt in
		o)
			WKDIR=$OPTARG
			echo "Working directory: $WKDIR"
			shift $((OPTIND-1)); OPTIND=1
			;;
		x)
			ORGANISM=$OPTARG
			shift $((OPTIND-1)); OPTIND=1
			;;
		t)
			THREADS=$OPTARG
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

# Also get target sequence location files from parameter
INPUTS=${@:1}

# Setting default values ${VARIABLE=DEFAULT_VALUE}
# This command will set VARIABLE to DEFAULT_VALUE if it is currently 
# undefined, then return VARIABLE
: ${WKDIR=$o_default}
: ${ORGANISM=$x_default}
: ${THREADS=$t_default}

# Check for mandatory parameters
if [ -z "${INPUTS}" ]; then
	echo "At least one input fasta file must be specified." >&2
	display_usage
	exit 1
fi

# Create the output folder
if [ ! -d ${WKDIR}/TMP ] ; then
	mkdir -p ${WKDIR}/TMP
fi

# Loop through all fasta files
ls -1 ${INPUTS[@]} | 
while read file; do
    
	outfile=`basename "$file" | sed -E 's/(.+)\..+?/\1/'`

	signalp \
		-fasta "${file}" \
		-format "long" \
		-gff3 \
		-mature \
		-org ${ORGANISM} \
		-prefix "${WKDIR}/${outfile}"
	
done

# Time scripts ends
echo "Completed $(date +"%T %d-%m-%Y")."
