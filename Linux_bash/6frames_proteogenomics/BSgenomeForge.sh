#!/bin/bash

# Script to generate a BSgenome package

# Time scripts starts
echo "$0"
echo "Start $(date +"%T %d-%m-%Y")."

# Function holding the usage
display_usage() { 
	echo "
		Usage: $0 [Options] <Genome fasta>"
	echo "
		Options:
        	-o	[str]	The output directory.
			-s	[str]	The BSgenome seed file.
        	-h	[]	To display the help.
	"
}

# Parse user parameters
while getopts "o:s:h" opt; do
	case $opt in
		o)
			WKDIR=$OPTARG
			echo "Working directory: $WKDIR"
			shift $((OPTIND-1)); OPTIND=1
			;;
		s)
            SEED=$OPTARG
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

# Also get target genome sequence files from parameter
FASTAS=${@:1}

# Check for mandatory parameters
if [ -z "${WKDIR}" ]; then
	echo "The output directory must be specified." >&2
	display_usage
	exit 1
fi
if [ -z "${FASTAS}" ]; then
	echo "At least one genome fasta file must be specified." >&2
	display_usage
	exit 1
fi
if [ -z "${SEED}" ]; then
	echo "The BSgenome seed file must be specified." >&2
	display_usage
	exit 1
fi

# Create the output folder
if [ ! -d ${WKDIR} ] ; then
	mkdir ${WKDIR}
fi

# Copy seed file and all fasta files into the output directory
cp ${SEED} ${WKDIR}
for file in `ls ${TARGETFILES}`; do
	cp ${file} ${WKDIR}
done

# Append the genome sequence directory to the seed file
if [[ `grep "seqs_srcdir" ${WKDIR}/basename ${SEED}` ]]; then
	echo "Value already present for seqs_srcdir at ${WKDIR}/`basename ${SEED}`!"
	exit 1
fi
echo "seqs_srcdir: ${WKDIR}" >> ${WKDIR}/`basename ${SEED}`

# Use R script to generate all configuration for BSgenome package
${PBS_O_HOME}/bin/BSgenome_forging.R -s ${WKDIR}/`basename ${SEED}` -f ${WKDIR} -o ${WKDIR}

# Build the source package
R CMD build ${WKDIR}

# Check the package
R CMD check ${WKDIR}/*.tar.gz

# Install the package
R CMD INSTALL ${WKDIR}/*.tar.gz

# Time scripts ends
echo "Completed $(date +"%T %d-%m-%Y")."
