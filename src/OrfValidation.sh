#!/bin/bash

# Script to compute and interprete the novel ORF identified from at least one MaxQuant processing

# Time scripts starts
echo "$0"
echo "Start $(date +"%T %d-%m-%Y")."

# Function holding the usage
display_usage() { 
	echo "
		Usage: $0 [Options] <Input file containing a table of MaxQuant processing path>"
	echo "
		Options:
            -o	[str]	The validation output directory.
            -r	[str]	The reference protein coordinate file.
            -n	[str]	The novel ORF protein coordinate file.
            -g	[str]	The genome as single Fasta file.
            -b	[str]	The BSgenome R package name.
            -w	[str]	The path to the reference best reciprocal blast file.
            -x	[str]	The path to the folder containing all best reciprocal blast files.
            -y	[str]	The path to the reference protein coordinates file.
            -z	[str]	The path to the novel ORF protein coordinates file.
            -t  [int]   The number of threads to use.
            -h	[]	To display the help.
	"
}

# Define default values
o_default="`pwd`/OrfValidation"
t_default=1

# Parse user parameters
while getopts "o:r:n:g:b:w:x:y:z:t:h" opt; do
	case $opt in
		o)
			WKDIR=$OPTARG
			echo "Working directory: $WKDIR"
			shift $((OPTIND-1)); OPTIND=1
			;;
		r)
			REFFASTA=$OPTARG
			shift $((OPTIND-1)); OPTIND=1
			;;
        n)
			ORFFASTA=$OPTARG
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
		w)
            REFRECIP=$OPTARG
            shift $((OPTIND-1)); OPTIND=1
			;;
        x)
            RECIP=$OPTARG
            shift $((OPTIND-1)); OPTIND=1
			;;
		y)
            REFCOORD=$OPTARG
            shift $((OPTIND-1)); OPTIND=1
			;;
		z)
            ORFCOORD=$OPTARG
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
: ${THREADS=$t_default}

# Check for mandatory parameters
if [ -z "${REFFASTA}" ]; then
	echo "The reference sequence fasta file must be specified." >&2
	display_usage
	exit 1
fi
if [ -z "${ORFFASTA}" ]; then
	echo "The ORF sequence fasta file must be specified." >&2
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
if [ -z "${REFRECIP}" ]; then
	echo "The reference reciprocal blast file must be specified." >&2
	display_usage
	exit 1
fi
if [ -z "${RECIP}" ]; then
	echo "The reciprocal blast folder must be specified." >&2
	display_usage
	exit 1
fi
if [ -z "${REFCOORD}" ]; then
	echo "The reference coordinates file must be specified." >&2
	display_usage
	exit 1
fi
if [ -z "${ORFCOORD}" ]; then
	echo "The ORF coordinates file must be specified." >&2
	display_usage
	exit 1
fi

# Loop through all location files
while IFS=$'\t' read -r prefix path; do

    if [ ! -d ${WKDIR}/${prefix} ] ; then
	    mkdir -p ${WKDIR}/${prefix}
    fi

	GetNovelEntries.R \
		-o ${WKDIR}/${prefix} \
		-m ${path}/combined/txt \
		-r ${REFFASTA} \
		-n ${ORFFASTA} \
		-t ${THREADS}
	Novelty_discovery_peptides.R \
		-e ${WKDIR}/${prefix}/Group_evidence.RDS \
		-f ${REFFASTA} \
		-r ${REFRECIP} \
		-u ${RECIP} \
		-p ${WKDIR}/${prefix}/Peptides_location.RDS \
		-t ${THREADS} \
		-o ${WKDIR}/${prefix}
	${PBS_O_HOME}/bin/SequenceGrange.sh \
		-o ${WKDIR}/${prefix} \
		-c ${WKDIR}/${prefix} \
		-n ${WKDIR}/${prefix}/Sequence_novelty_reason.RDS \
		-r ${REFCOORD} \
		-f ${ORFCOORD} \
		-g ${GENOME} \
		-b ${BSGENOME} \
		${WKDIR}/${prefix}
   
done < "$INPUTS"

# Time scripts ends
echo "Completed $(date +"%T %d-%m-%Y")."
