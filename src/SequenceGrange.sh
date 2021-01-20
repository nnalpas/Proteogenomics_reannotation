#!/bin/bash

# Script to generate Grange from a list of sequence location files

# Time scripts starts
echo "$0"
echo "Start $(date +"%T %d-%m-%Y")."

# Function holding the usage
display_usage() { 
	echo "
		Usage: $0 [Options] <Input folder containing file(s) '*_location.RDS'> <...>"
	echo "
		Options:
        	-o	[str]	The GRange output directory.
			-c	[str]	The coordinates output directory.
			-n	[str]	The sequence novelty reason file.
			-r	[str]	The reference protein coordinate file.
			-f	[str]	The ORF protein coordinate file.
			-g	[str]	The genome as single Fasta file.
			-b	[str]	The BSgenome R package name.
        	-h	[]	To display the help.
	"
}

# Define default values
o_default="`pwd`/GRanges"
c_default="`pwd`/Coordinates"

# Parse user parameters
while getopts "o:c:n:r:f:g:b:h" opt; do
	case $opt in
		o)
			WKDIR=$OPTARG
			echo "Working directory: $WKDIR"
			shift $((OPTIND-1)); OPTIND=1
			;;
		c)
			INTERDIR=$OPTARG
			echo "Intermediary directory: $INTERDIR"
			shift $((OPTIND-1)); OPTIND=1
			;;
		n)
            NOVELTY=$OPTARG
            shift $((OPTIND-1)); OPTIND=1
			;;
		r)
            REF=$OPTARG
            shift $((OPTIND-1)); OPTIND=1
			;;
		f)
            ORF=$OPTARG
            shift $((OPTIND-1)); OPTIND=1
			;;
		g)
            GENOME=$OPTARG
            shift $((OPTIND-1)); OPTIND=1
			;;
		b)
            BSGENOME=$OPTARG
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
: ${INTERDIR=$c_default}

# Check for mandatory parameters
if [ -z "${NOVELTY}" ]; then
	echo "The sequence novelty reason file must be specified." >&2
	display_usage
	exit 1
fi
if [ -z "${REF}" ]; then
	echo "The reference protein coordinate file must be specified." >&2
	display_usage
	exit 1
fi
if [ -z "${ORF}" ]; then
	echo "The ORF protein coordinate file must be specified." >&2
	display_usage
	exit 1
fi
if [ -z "${GENOME}" ]; then
	echo "The genome as single Fasta file must be specified." >&2
	display_usage
	exit 1
fi
if [ -z "${BSGENOME}" ]; then
	echo "The BSgenome R package name must be specified." >&2
	display_usage
	exit 1
fi

# Create the output folder
if [ ! -d ${WKDIR} ] ; then
	mkdir ${WKDIR}
fi
if [ ! -d ${INTERDIR} ] ; then
	mkdir ${INTERDIR}
fi

# Loop through all location files
find ${INPUTS[@]} -type f -name "*_location.RDS" -print0 | 
while IFS= read -r -d '' file; do
    
	echo "Genomic_position_from_within_proteins.R \
		-p ${file} \
		-r ${NOVELTY} \
		-k ${REF} \
		-n ${ORF} \
		-o ${INTERDIR}"
	
	coord_file=`basename "$file" | sed "s/_location.RDS/_coordinates.txt/"`
	
	echo "GRanges_generation.R \
		-c ${INTERDIR}/${coord_file} \
		-g ${GENOME} \
		-b ${BSGENOME} \
		-o ${WKDIR}"
	
done

# Time scripts ends
echo "Completed $(date +"%T %d-%m-%Y")."
