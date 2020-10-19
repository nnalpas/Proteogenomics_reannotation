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
	if [[ ! -z "${array[7]}" ]]; then BLASTCMD+="-x ${array[7]} ; "; fi;
	if [[ ! -z "${array[8]}" ]]; then BLASTCMD+="-y ${array[8]} ; "; fi;
	
	# Check software compatibility with selected task
	if [[ "${array[0]}" != "blastp" ]] && [[ "${array[0]}" != "blastx" ]]; then
		echo "The selected task: ${array[0]} is compatible only with Blast software." >&2
		SOFTWARE="Blast"
	fi
	
	# Run the blast for each iteration
	#if [[ "$SOFTWARE" == "Blast" ]]; then
	#	${PBS_O_HOME}/bin/BlastDbBlasting.sh $BLASTCMD -t ${THREADS} 2>&1
	#elif [[ "$SOFTWARE" == "Diamond" ]]; then
	#	${PBS_O_HOME}/bin/DiamondDbBlasting.sh $BLASTCMD -t ${THREADS} 2>&1
	#else
	#	echo "The blasting software must be either 'Blast' or 'Diamond'." >&2
	#	exit 1
	#fi
	
	# Retrieve blast ids from blast results and the corresponding sequence header
	echo -e "qseqid\tDescription\tTaxon\tTaxonID" > ${WKDIR}/${array[3]}_header
	cut -f 1 ${WKDIR}/${array[3]} > ${WKDIR}/${array[3]}_qseqid.txt
	blastdbcmd -db ${QUERY} -outfmt "%a;%t;%S;%T" -target_only -entry_batch ${WKDIR}/${array[3]}_qseqid.txt | sed 's/;;/\t/g' >> ${WKDIR}/${array[3]}_header
	rm ${WKDIR}/${array[3]}_qseqid.txt
	
	# If the qseqid is identical (after retrieving the header) then concatenate the new columns with the reciprocal results
	diff_count=`diff <(awk '{print $1}' ${WKDIR}/${array[3]}) <(awk '{print $1}' ${WKDIR}/${array[3]}_header) | wc -l`
	if (( "$diff_count" == 0 )); then
		echo "Including annotation into: ${WKDIR}/${array[3]}_annot"
		cut -f 2-4 ${WKDIR}/${array[3]}_header > ${WKDIR}/${array[3]}_header.tmp
		paste -d '\t' ${WKDIR}/${array[3]} ${WKDIR}/${array[3]}_header.tmp > ${WKDIR}/${array[3]}_annot
		rm ${WKDIR}/${array[3]}_header.tmp
	else
		echo "Error while retrieving annotation!"
	fi
		
done < ${CROSSMAP}

# Time scripts ends
echo "Completed $(date +"%T %d-%m-%Y")."
