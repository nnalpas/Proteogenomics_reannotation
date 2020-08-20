#!/bin/bash

# Script to blast two blast databases against each other using Diamond

# Time scripts starts
echo "$0"
echo "$@"
echo "Start $(date +"%T %d-%m-%Y")."

# Function holding the usage
display_usage() { 
	echo "
		Usage: $0 [Options] -a <Task> -q <Query> -d <Database>  -b <Basename>"
	echo "
		Options:
			-o	[str]	The output directory.
			-l	[str]	Which entry to use, either 'all' or a file containing entries
			-a	[str]	What blast to perform ('blastx' or 'blastp')
			-q	[str]	The query sequences
			-d	[str]	The database (subject)
			-e	[float]	The e-value threshold for blast to be reported
			-n	[int]	The maximum number of alignment to report
			-x	[str]	Additional blast parameters.
			-y	[str]	Additional makedb parameters.
			-b	[str]	The name of the output file
			-t	[int]	Number of threads.
			-h	[]		To display the help.
	"
}

# Define default values
o_default="`pwd`/Blast"
l_default="all"
e_default=0.1
n_default=100
t_default=1
x_default=""
y_default=""

# Parse user parameters
while getopts "o:l:a:q:d:e:n:x:y:b:t:h" opt; do
	case $opt in
		o)
			WKDIR=$OPTARG
			echo "Working directory: $WKDIR"
			shift $((OPTIND-1)); OPTIND=1
			;;
		l)
			ENTRY=$OPTARG
			shift $((OPTIND-1)); OPTIND=1
			;;
		a)
			TASK=$OPTARG
			shift $((OPTIND-1)); OPTIND=1
			;;
		q)
			QUERY=$OPTARG
			shift $((OPTIND-1)); OPTIND=1
			;;
		d)
			DATABASE=$OPTARG
			shift $((OPTIND-1)); OPTIND=1
			;;
		e)
			EVAL=$OPTARG
			shift $((OPTIND-1)); OPTIND=1
			;;
		n)
			NUMALIGN=$OPTARG
			shift $((OPTIND-1)); OPTIND=1
			;;
		b)
			BASENAME=$OPTARG
			shift $((OPTIND-1)); OPTIND=1
			;;
		x)
			while [[ -n "${!OPTIND}" ]] && [[ ${!OPTIND} != ';' ]]; do
				BLAST_ADD+=("${!OPTIND}")
				let OPTIND++
			done
			let OPTIND++
			#BLAST_ADD=$OPTARG
			#shift $((OPTIND-1)); OPTIND=1
			;;
		y)
			while [[ -n "${!OPTIND}" ]] && [[ "${!OPTIND}" != ';' ]]; do
				MAKEDB_ADD+=("${!OPTIND}")
				let OPTIND++
			done
			let OPTIND++
			#MAKEDB_ADD=$OPTARG
			#shift $((OPTIND-1)); OPTIND=1
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

# Setting default values ${VARIABLE=DEFAULT_VALUE}
# This command will set VARIABLE to DEFAULT_VALUE if it is currently 
# undefined, then return VARIABLE
: ${WKDIR=$o_default}
: ${ENTRY=$l_default}
: ${EVAL=$e_default}
: ${NUMALIGN=$n_default}
: ${THREADS=$t_default}
: ${BLAST_ADD=$x_default}
: ${MAKEDB_ADD=$y_default}

if [ -z "${TASK}" ]; then
	echo "The blast type to perform must be specified." >&2
	display_usage
	exit 1
fi
if [ -z "${QUERY}" ]; then
	echo "The query sequences must be specified." >&2
	display_usage
	exit 1
fi
if [ -z "${DATABASE}" ]; then
	echo "The database (subject) must be specified." >&2
	display_usage
	exit 1
fi
if [ -z "${BASENAME}" ]; then
	echo "The output file basename must be specified." >&2
	display_usage
	exit 1
fi

# Create the output folder
if [ ! -d ${WKDIR} ] ; then
	mkdir ${WKDIR}
fi

# Make a database for the query if it does not already exists
if [ ! -e ${DATABASE}.dmnd ] ; then
	echo "diamond makedb --in ${DATABASE} \
		--db ${DATABASE}.dmnd ${MAKEDB_ADD[@]}"
fi

# Check whether queries should be filtered
entry_retrieval="${ENTRY}"
if [ ${ENTRY} != 'all' ]; then
	entry_retrieval="-entry ${ENTRY}"
fi

# All entries retrieval and blasting of retrieved entries against another database
echo "blastdbcmd -db ${QUERY} \
       -dbtype ${DBTYPE} \
       ${entry_retrieval} | \
       eval diamond ${TASK} \
       --query - \
       --db ${DATABASE} \
       --out ${WKDIR}/${BASENAME} \
       --evalue ${EVAL} \
       --max-target-seqs ${NUMALIGN} \
       --threads ${THREADS} ${BLAST_ADD[@]} \
       --outfmt '6 qseqid sseqid pident nident mismatch gaps length gapopen qstart qend qlen qframe qseq qstrand sstart send slen sframe sseq staxids sscinames evalue bitscore score'"

# Time scripts ends
echo "Completed $(date +"%T %d-%m-%Y")."
