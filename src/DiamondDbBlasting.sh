#!/bin/bash

# Script to blast two blast databases against each other using Diamond

# Time scripts starts
echo "$0"
echo "Start $(date +"%T %d-%m-%Y")."

# Load environment module
# Load environment module depending on server
if [[ `hostname` == *"core-login"* ]] || [[ `hostname` == *"ifb"* ]]; then
	module load blast/2.12.0
	module load diamond/2.0.9
elif [[ `hostname` == *"binac"* ]]; then
	module load bio/blastplus/2.11.0
	module load diamond/2.0.13
else
	echo "Unknown server: $(hostname), do not know what module to load" >&2
	exit 1
fi

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
			BLAST_ADD=("$OPTARG")
			while [[ -n "${!OPTIND}" ]] && [[ ${!OPTIND} != ';' ]]; do
				BLAST_ADD+=("${!OPTIND}")
				let OPTIND++
			done
			let OPTIND++
			;;
		y)
			MAKEDB_ADD=("$OPTARG")
			while [[ -n "${!OPTIND}" ]] && [[ "${!OPTIND}" != ';' ]]; do
				MAKEDB_ADD+=("${!OPTIND}")
				let OPTIND++
			done
			let OPTIND++
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
dir=$( dirname ${DATABASE} )
base=$( basename ${DATABASE} )
if [[ -e ${DATABASE}.dmnd ]]; then
	echo "Provided database is already formatted for Diamond!"
elif [[ -e ${dir}/DiamondDB/${base}.dmnd ]]; then
	echo "Update provided database for Diamond children folder!"
	DATABASE="${dir}/DiamondDB/${base}"
else
	mkdir ${dir}/DiamondDB
	for (( i=0; i<(${#MAKEDB_ADD[@]}-1); i++ )); do
		j=$(($i+1))
		#if [[ "${MAKEDB_ADD[$i]}" == '--taxonmap' ]] && [[ ! "${MAKEDB_ADD[$j]}" =~ "prot.accession2taxid.gz$" ]]; then
		if [[ "${MAKEDB_ADD[$i]}" == '--taxonmap' ]] && [[ ! "${MAKEDB_ADD[$j]}" =~ "prot.accession2taxid.FULL.gz" ]]; then
			#grep -E "^>" ${DATABASE} | sed 's/^>//' | sed 's/ .*//' | awk -v taxid="${MAKEDB_ADD[$j]}" '{ print $1, "\t", $1, "\t", taxid, "\t", 0 }' > ${DATABASE}.prot.accession2taxid
			echo -e "accession.version\ttaxid" > ${dir}/DiamondDB/${base}.prot.accession2taxid.FULL
			grep -E "^>" ${DATABASE} | sed 's/^>//' | sed 's/ .*//' | awk -v taxid="${MAKEDB_ADD[$j]}" '{ print $1, "\t", taxid }' >> ${dir}/DiamondDB/${base}.prot.accession2taxid.FULL
			#gzip ${DATABASE}.prot.accession2taxid
			gzip ${dir}/DiamondDB/${base}.prot.accession2taxid.FULL
			#MAKEDB_ADD[$j]="${DATABASE}.prot.accession2taxid.gz"
			MAKEDB_ADD[$j]="${dir}/DiamondDB/${base}.prot.accession2taxid.FULL.gz"
		else
			eval MAKEDB_ADD[$j]=${MAKEDB_ADD[$j]}
		fi
	done
	diamond makedb --in ${DATABASE} \
		--db ${dir}/DiamondDB/${base}.dmnd ${MAKEDB_ADD[@]}
	DATABASE="${dir}/DiamondDB/${base}"
fi

# Check whether queries should be filtered
entry_retrieval="-entry_batch ${ENTRY}"
if [ ${ENTRY} == 'all' ]; then
	entry_retrieval="-entry ${ENTRY}"
fi

# All entries retrieval and blasting of retrieved entries against another database
blastdbcmd -db ${QUERY} \
       ${entry_retrieval} | \
       eval "diamond ${TASK} \
       --db ${DATABASE} \
       --out ${WKDIR}/${BASENAME} \
       --evalue ${EVAL} \
       --max-target-seqs ${NUMALIGN} \
       --threads ${THREADS} ${BLAST_ADD[@]} \
       --outfmt 6 qseqid sseqid pident nident mismatch gaps length gapopen qstart qend qlen qframe qseq qstrand sstart send slen sseq staxids sscinames evalue bitscore score"

# Include column names at start of file
sed -i '1s/^/qseqid\tsseqid\tpident\tnident\tmismatch\tgaps\tlength\tgapopen\tqstart\tqend\tqlen\tqframe\tqseq\tqstrand\tsstart\tsend\tslen\tsseq\tstaxid\tssciname\tevalue\tbitscore\tscore\n/' ${WKDIR}/${BASENAME}

# Time scripts ends
echo "Completed $(date +"%T %d-%m-%Y")."
