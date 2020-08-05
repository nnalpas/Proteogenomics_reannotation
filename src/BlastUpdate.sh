#!/bin/bash

# Download pre-formatted BLAST database from NCBI

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
        	-t	[int]	Number of threads.
        	-h	[]	To display the help.
	"
}

# Define default values
o_default="`pwd`/BlastDB"
t_default=1

# Parse user parameters
while getopts "o:t:h" opt; do
	case $opt in
		o)
			WKDIR=$OPTARG
			echo "Working directory: $WKDIR"
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

# Get all database anme to download
DB=${@:1}

# Setting default values ${VARIABLE=DEFAULT_VALUE}
# This command will set VARIABLE to DEFAULT_VALUE if it is currently 
# undefined, then return VARIABLE
: ${WKDIR=$o_default}
: ${THREADS=$t_default}

# Check for mandatory parameters
if [ -z "${DB}" ]; then
	echo "The database must be specified. If you are not sure try 'update_blastdb.pl --showall'" >&2
	display_usage
	exit 1
fi

# Create the output folder
if [ ! -d ${WKDIR} ] ; then
	mkdir ${WKDIR}
fi
cd $WKDIR

# Retrieve all databases
update_blastdb.pl \
	--decompress \
	--timeout 240 \
	--num_threads $THREADS \
	$DB

# Time scripts ends
echo "Completed $(date +"%T %d-%m-%Y")."


