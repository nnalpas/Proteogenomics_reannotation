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
for file in `ls ${FASTAS}`; do
	nmb_entries=`grep -c "^>" $file`
	if [[ "$nmb_entries" -eq 1 ]]; then
		name=`grep "^>" $FASTAS | sed -E "s/^>([^ ]*).*/\1/"`
		ext=`basename ${file} | sed -E "s/.*(\.fa(sta)?(\.gz)?)$/\1/"`
		cp ${file} ${WKDIR}/${name}${ext}
	else
		echo "More than one fasta header at: ${file}!"
		exit 1
	fi
done

# Append the genome sequence directory to the seed file
NEWSEED=`basename ${SEED}`
if [[ `grep "seqs_srcdir" ${WKDIR}/${NEWSEED}` ]]; then
	echo "Value already present for seqs_srcdir at ${WKDIR}/`basename ${SEED}`!"
	exit 1
fi
echo "seqs_srcdir: ${WKDIR}" >> ${WKDIR}/${NEWSEED}

# Use R script to generate all configuration for BSgenome package
${PBS_O_HOME}/bin/BSgenome_forging.R -s ${WKDIR}/${NEWSEED} -f ${WKDIR} -o ${WKDIR}

# Build the source package
R CMD build ${WKDIR}

# Check for presence of a package tar.gz file
for pkg in `ls *.tar.gz`; do
	
	# Check the package
	R CMD check ${WKDIR}/${pkg}
	
	# Install the package
	R CMD INSTALL ${WKDIR}/${pkg}
	
done

# Time scripts ends
echo "Completed $(date +"%T %d-%m-%Y")."
