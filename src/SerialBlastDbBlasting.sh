#!/bin/bash

# Script to initiate serial blast based on parameters stored on file

# Time scripts starts
echo "$0"
echo "Start $(date +"%T %d-%m-%Y")."

# Function holding the usage
display_usage() { 
	echo "
		Usage: $0 [Options] -x <Cross map file>"
	echo "
		Options:
        	-o	[str]	The output directory.
			-s	[str]	The blasting software too use (i.e. 'Blast' or 'Diamond')
			-x	[str]	A 9 columns tab-separated file containing a list of blast to perform (<Task> <Query> <Database> <Basename> <Subject filter> <E value> <Number of alignment> <Blasting additional parameters> <Make DB additional parameters>)
			-t	[int]	Number of threads.
        	-h	[]	To display the help.
	"
}

# Define default values
o_default="`pwd`/Blast"
t_default=1

# Parse user parameters
while getopts "o:s:x:t:h" opt; do
	case $opt in
		o)
			WKDIR=$OPTARG
			echo "Working directory: $WKDIR"
			shift $((OPTIND-1)); OPTIND=1
			;;
		s)
			SOFTWARE=$OPTARG
			shift $((OPTIND-1)); OPTIND=1
			;;
		x)
			CROSSMAP=$OPTARG
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
: ${THREADS=$t_default}

# Check for mandatory parameters
if [ -z "${SOFTWARE}" ]; then
	echo "The blasting software must be specified." >&2
	display_usage
	exit 1
fi
if [ -z "${CROSSMAP}" ]; then
	echo "The blast crossmap file must be specified." >&2
	display_usage
	exit 1
fi

# Read the crossmap file containing blast instructions
while IFS= read -r line; do
	
	IFS=$'\t' read -r -a array <<< "$line"
	
	# Generate the complete blast command for current entry
	BLASTCMD="-o ${WKDIR} "
	if [[ ! -z "${array[0]}" ]]; then BLASTCMD+="-a ${array[0]} "; fi;
	if [[ ! -z "${array[1]}" ]]; then eval QUERY=${array[1]}; BLASTCMD+="-q $QUERY "; fi;
	if [[ ! -z "${array[2]}" ]]; then eval DATABASE=${array[2]}; BLASTCMD+="-d $DATABASE "; fi;
	if [[ ! -z "${array[3]}" ]]; then BLASTCMD+="-b ${array[3]} "; fi;
	if [[ ! -z "${array[4]}" ]]; then eval LIST=${array[4]}; BLASTCMD+="-l $LIST "; fi;
	if [[ ! -z "${array[5]}" ]]; then BLASTCMD+="-e ${array[5]} "; fi;
	if [[ ! -z "${array[6]}" ]]; then BLASTCMD+="-n ${array[6]} "; fi;
	if [[ ! -z "${array[7]}" ]]; then eval BLAST_ADD=${array[7]}; BLASTCMD+="-x \"$BLAST_ADD\" "; fi;
	if [[ ! -z "${array[8]}" ]]; then eval DB_ADD=${array[8]}; BLASTCMD+="-y \"${DB_ADD}\" "; fi;
	
	# Run the blast for each iteration
	if [[   ]]; then
		echo "${PBS_O_HOME}/bin/BlastDbBlasting.sh $BLASTCMD -t ${THREADS}" 2>&1
	elif
		echo "${PBS_O_HOME}/bin/DiamondDbBlasting.sh $BLASTCMD -t ${THREADS}" 2>&1
	else
		echo "The blasting software must be either 'Blast' or 'Diamond'." >&2
		exit 1
	fi
	
done < ${CROSSMAP}

# Time scripts ends
echo "Completed $(date +"%T %d-%m-%Y")."
