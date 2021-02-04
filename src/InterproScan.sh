#!/bin/bash

# Script to provides functional analysis of proteins by classifying them into families and predicting domains and important sites via InterProScan.

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
			-s	[str]	The type of the input sequences (dna/rna (n) or protein (p) [default]).
			-t	[int]	The number of threads.
			-h	[]	To display the help.
	"
}

# Define default values
o_default="`pwd`/InterproScan"
s_default="p"
t_default=1

# Parse user parameters
while getopts "o:s:t:h" opt; do
	case $opt in
		o)
			WKDIR=$OPTARG
			echo "Working directory: $WKDIR"
			shift $((OPTIND-1)); OPTIND=1
			;;
		s)
			SEQTYPE=$OPTARG
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
: ${SEQTYPE=$s_default}
: ${THREADS=$t_default}

# Check for mandatory parameters
if [ -z "${INPUTS}" ]; then
	echo "At least one input folder (containing proteinGroups.txt or grange.RDS) must be specified." >&2
	display_usage
	exit 1
fi

# Create the output folder
if [ ! -d ${WKDIR}/TMP ] ; then
	mkdir -p ${WKDIR}/TMP
fi

# Loop through all proteinGroups or sites files
ls -1 ${INPUTS[@]} | 
while read file; do
    
	outfile=`basename "$file" | sed -E 's/(.+)\..+?/\1/'`
	
	interproscan.sh \
		--output-file-base "${WKDIR}/${outfile}" \
		--input "${file}" \
		--goterms \
		--pathways \
		--seqtype ${SEQTYPE} \
		--tempdir "${WKDIR}/TMP" \
		--cpu ${THREADS}

	sed -iE "1 i\Protein accession\tSequence MD5 digest\tSequence length\tAnalysis\tSignature accession\tSignature description\tStart location\tStop location\tScore\tStatus\tDate\tInterPro accession\tInterPro description\tGO annotations\tPathways annotations" ${outfile}.tsv
	
done

# Time scripts ends
echo "Completed $(date +"%T %d-%m-%Y")."
