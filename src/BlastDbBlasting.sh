#!/bin/bash

# Script to blast two blast databases against each other using Blast

# Time scripts starts
echo "$0"
echo "Start $(date +"%T %d-%m-%Y")."

# Load environment module depending on server
if [[ `hostname` == *"core-login"* ]] || [[ `hostname` == *"ifb"* ]] || [[ `hostname` == *"cpu-node"* ]]; then
	module load blast/2.12.0
elif [[ `hostname` == *"binac"* ]]; then
	module load bio/blastplus/2.11.0
else
	echo "Unknown server: $(hostname), do not know what module to load" >&2
	exit 1
fi

# Function holding the usage
display_usage() { 
	echo "
		Usage: $0 [Options] -a <Task> -q <Query> -d <Database> -b <Basename>"
	echo "
		Options:
        	-o	[str]	The output directory.
			-l	[str]	Which entry to use, either 'all' or a file containing entries
			-a	[str]	What blast to perform ('blastn' or 'blastp'...)
			-q	[str]	The query sequences (as database)
			-d	[str]	The database (subject)
			-e	[float]	The e-value threshold for blast to be reported
			-n	[int]	The maximum number of alignment to report
			-x	[str]	Additional blast parameters.
			-b	[str]	The name of the output file
			-t	[int]	Number of threads.
        	-h	[]	To display the help.
	"
}

# Define default values
o_default="`pwd`/Blast"
l_default="all"
e_default=0.1
n_default=100
t_default=1
x_default=""

# Parse user parameters
while getopts "o:l:a:q:d:e:n:x:b:t:h" opt; do
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
			BLAST_ADD=$OPTARG
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
			;;
		:)
			display_usage
			exit 1
			;;
	esac
done

shift $((OPTIND-1)); OPTIND=1
BLAST_ADD=$@

# Setting default values ${VARIABLE=DEFAULT_VALUE}
# This command will set VARIABLE to DEFAULT_VALUE if it is currently 
# undefined, then return VARIABLE
: ${WKDIR=$o_default}
: ${ENTRY=$l_default}
: ${EVAL=$e_default}
: ${NUMALIGN=$n_default}
: ${THREADS=$t_default}
: ${BLAST_ADD=$x_default}

# Check for mandatory parameters
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

# Check whether subject should be filtered
entry_retrieval="-entry_batch ${ENTRY}"
if [ ${ENTRY} == 'all' ]; then
	entry_retrieval="-entry ${ENTRY}"
fi

# All entries retrieval and blasting of retrieved entries against another database
echo "${TASK} \
       -query - \
       -task ${TASK} \
       -db ${DATABASE} \
       -out ${WKDIR}/${BASENAME} \
       -evalue ${EVAL} \
       -num_alignments ${NUMALIGN} \
       -num_threads ${THREADS} ${BLAST_ADD} \
       -outfmt '6 qseqid sseqid pident nident mismatch gaps length gapopen qstart qend qlen qframe qseq sstart send slen sframe sseq staxid ssciname sstrand evalue bitscore score'"

blastdbcmd -db ${QUERY} \
       ${entry_retrieval} | \
       eval "${TASK} \
       -query - \
       -task ${TASK} \
       -db ${DATABASE} \
       -out ${WKDIR}/${BASENAME} \
       -evalue ${EVAL} \
       -num_alignments ${NUMALIGN} \
       -num_threads ${THREADS} ${BLAST_ADD} \
       -outfmt '6 qseqid sseqid pident nident mismatch gaps length gapopen qstart qend qlen qframe qseq sstart send slen sframe sseq staxid ssciname sstrand evalue bitscore score'"

# Include column names at start of file
sed -i '1s/^/qseqid\tsseqid\tpident\tnident\tmismatch\tgaps\tlength\tgapopen\tqstart\tqend\tqlen\tqframe\tqseq\tsstart\tsend\tslen\tsframe\tsseq\tstaxid\tssciname\tsstrand\tevalue\tbitscore\tscore\n/' ${WKDIR}/${BASENAME}

# Time scripts ends
echo "Completed $(date +"%T %d-%m-%Y")."
