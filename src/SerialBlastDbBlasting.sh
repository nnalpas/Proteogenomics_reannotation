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
			-x	[str]	A 9 columns tab-separated file containing a list of blast to perform (<DB type> <Task> <Subject> <Query> <Basename> <Subject filter> <E value> <Number of alignment> <Additional parameters>)
			-t	[int]	Number of threads.
        	-h	[]	To display the help.
	"
}

# Define default values
o_default="`pwd`/Blast"
t_default=1

# Parse user parameters
while getopts "o:x:t:h" opt; do
	case $opt in
		o)
			WKDIR=$OPTARG
			echo "Working directory: $WKDIR"
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
	if [[ ! -z "${array[0]}" ]]; then BLASTCMD+="-y ${array[0]} "; fi;
	if [[ ! -z "${array[1]}" ]]; then BLASTCMD+="-a ${array[1]} "; fi;
	if [[ ! -z "${array[2]}" ]]; then BLASTCMD+="-s ${array[2]} "; fi;
	if [[ ! -z "${array[3]}" ]]; then BLASTCMD+="-q ${array[3]} "; fi;
	if [[ ! -z "${array[4]}" ]]; then BLASTCMD+="-b ${array[4]} "; fi;
	if [[ ! -z "${array[5]}" ]]; then BLASTCMD+="-l ${array[5]} "; fi;
	if [[ ! -z "${array[6]}" ]]; then BLASTCMD+="-e ${array[6]} "; fi;
	if [[ ! -z "${array[7]}" ]]; then BLASTCMD+="-n ${array[7]} "; fi;
	if [[ ! -z "${array[8]}" ]]; then BLASTCMD+="-m ${array[8]} "; fi;
	
	# Run the blast for each iteration
	echo "${PBS_O_HOME}/bin/BlastDbBlasting.sh ${BLASTCMD} -t ${THREADS}" >&2
	
done < ${CROSSMAP}

# Time scripts ends
echo "Completed $(date +"%T %d-%m-%Y")."
