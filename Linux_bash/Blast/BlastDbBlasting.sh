#!/bin/bash

# Script to blast two blast databases against each other

# Time scripts starts
echo "$0"
echo "Start $(date +"%T %d-%m-%Y")."

# Function holding the usage
display_usage() { 
	echo "
		Usage: $0 [Options]"
	echo "
		Options:
        	-o	[str]	The output directory.
			-y	[str]	The type of the data ('nucl' or 'prot')
			-l	[str]	Which entry to use, either 'all' or a file containing entries
			-a	[str]	What blast to perform ('blastn' or 'blastp'...)
			-s	[str]	The subject database
			-q	[str]	The query database
			-e	[float]	The e-value threshold for blast to be reported
			-n	[int]	The maximum number of alignment to report
        	-t	[int]	Number of threads. Only applicable to FastQC.
			-b	[str]	The name of the output file
        	-h	[]	To display the help.
	"
}

# Define default values
o_default="`pwd`/Blast"
l_default="all"
e_default=0.1
n_default=100
t_default=1

# Parse user parameters
while getopts "o:y:l:a:s:q:e:n:b:t:h" opt; do
	case $opt in
		o)
			WKDIR=$OPTARG
			echo "Working directory: $WKDIR"
			shift $((OPTIND-1)); OPTIND=1
			;;
        y)
            DBTYPE=$OPTARG
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
		s)
            SUBJECT=$OPTARG
            shift $((OPTIND-1)); OPTIND=1
			;;
		q)
            QUERY=$OPTARG
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

# Check for mandatory parameters
if [ -z "${DBTYPE}" ]; then
	echo "The database type must be specified." >&2
	display_usage
	exit 1
fi
if [ -z "${TASK}" ]; then
	echo "The blast type to perform must be specified." >&2
	display_usage
	exit 1
fi
if [ -z "${SUBJECT}" ]; then
	echo "The subject database must be specified." >&2
	display_usage
	exit 1
fi
if [ -z "${QUERY}" ]; then
	echo "The query database must be specified." >&2
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

if [ ${ENTRY} == 'all' ]; then

	# All entries retrieval and blasting of retrieved entries against another database
	blastdbcmd -db ${SUBJECT} \
	  -dbtype ${DBTYPE} \
	  -entry ${ENTRY} | \
	  eval "${TASK} \
	  -query - \
	  -task ${TASK} \
	  -db ${QUERY} \
	  -out ${WKDIR}/${BASENAME} \
	  -evalue ${EVAL} \
	  -num_alignments ${NUMALIGN} \
	  -num_threads ${THREADS} \
	  -outfmt '6 qseqid sseqid pident nident mismatch length gapopen qstart qend sstart send evalue bitscore score'"

else

	# Specific entries retrieval and blasting of retrieved entries against another database
	blastdbcmd -db ${SUBJECT} \
	  -dbtype ${DBTYPE} \
	  -entry_batch ${ENTRY} | \
	  eval "${TASK} \
	  -query - \
	  -task ${TASK} \
	  -db ${QUERY} \
	  -out ${WKDIR}/${BASENAME} \
	  -evalue ${EVAL} \
	  -num_alignments ${NUMALIGN} \
	  -num_threads ${THREADS} \
	  -outfmt '6 qseqid sseqid pident nident mismatch length gapopen qstart qend sstart send evalue bitscore score'"

fi

# Time scripts ends
echo "Completed $(date +"%T %d-%m-%Y")."
