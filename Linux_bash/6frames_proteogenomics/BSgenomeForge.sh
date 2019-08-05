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
			-x	[]	To split a fasta file containing multiple sequences.
        	-h	[]	To display the help.
	"
}

# Define default values
o_default="`pwd`/BSgenome"
x_default=false

# Parse user parameters
while getopts "o:s:xh" opt; do
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
		x)
            SPLIT=true
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

# Setting default values ${VARIABLE=DEFAULT_VALUE}
# This command will set VARIABLE to DEFAULT_VALUE if it is currently 
# undefined, then return VARIABLE
: ${WKDIR=$o_default}
: ${SPLIT=$x_default}

# Check for mandatory parameters
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
CURDIR=`pwd`
cp ${SEED} ${WKDIR}
if $SPLIT; then
	echo "More than one fasta header at: ${file}! Attempting to split!"
	cd $WKDIR
	faidx -x $file
	cd $CURDIR
else
	cp -t ${WKDIR} `ls ${FASTAS}`
fi

# Append the fasta sequence name to the seed file
NEWSEED=`basename ${SEED}`
if [[ `grep "seqnames" ${WKDIR}/${NEWSEED}` ]]; then
	echo "Value already present for seqnames at ${WKDIR}/${NEWSEED}!"
	exit 1
fi
seqnames=""
for file in `ls ${WKDIR}/*.fa`; do
	seqs=`basename ${file} | sed -E "s/(.*)\.fa(sta)?(\.gz)?$/\1/"`
	seqnames="$seqnames '${seqs}'"
done
echo "seqnames: c(${seqnames})" >> ${WKDIR}/${NEWSEED}

# Append the genome sequence directory to the seed file
if [[ `grep "seqs_srcdir" ${WKDIR}/${NEWSEED}` ]]; then
	echo "Value already present for seqs_srcdir at ${WKDIR}/${NEWSEED}!"
	exit 1
fi
echo "seqs_srcdir: ${WKDIR}" >> ${WKDIR}/${NEWSEED}

# Use R script to generate all configuration for BSgenome package
BSgenome_forging.R -s ${WKDIR}/${NEWSEED} -f ${WKDIR} -o ${WKDIR}

# Build the source package
#pkg=`grep "^Package: " ${WKDIR}/${NEWSEED} | sed -E "s/^Package: (.+)$/\1/"`
#cd ${WKDIR}
#R CMD build ${WKDIR}/${pkg}
#cd ${PBS_O_WORKDIR}

#exit 1

# Check for presence of a package tar.gz file
for tarball in `ls *.tar.gz`; do
	
	# Check the package
	R CMD check ${WKDIR}/${tarball}
	
	# Install the package
	R CMD INSTALL ${WKDIR}/${tarball}
	
done

# Time scripts ends
echo "Completed $(date +"%T %d-%m-%Y")."
