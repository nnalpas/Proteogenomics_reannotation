#!/bin/bash

# Script to generate expression and phenodata files from MaxQuant processing folder

# Time scripts starts
echo "$0"
echo "Start $(date +"%T %d-%m-%Y")."

# Function holding the usage
display_usage() { 
	echo "
		Usage: $0 [Options] <Input folder 1 (containing proteinGroups.txt or grange.RDS)> <...> <folder N>"
	echo "
		Options:
        	-o	[str]	The output directory.
        	-h	[]	To display the help.
	"
}

# Define default values
o_default="`pwd`/SummarizedExp"

# Parse user parameters
while getopts "o:h" opt; do
	case $opt in
		o)
			WKDIR=$OPTARG
			echo "Working directory: $WKDIR"
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

# Check for mandatory parameters
if [ -z "${INPUTS}" ]; then
	echo "At least one input folder (containing proteinGroups.txt or grange.RDS) must be specified." >&2
	display_usage
	exit 1
fi

# Create the output folder
if [ ! -d ${WKDIR} ] ; then
	mkdir ${WKDIR}
fi

# Loop through all proteinGroups or sites files
find ${INPUTS[@]} -type f \( -name "*proteinGroups.txt" -o -name "*Sites_grange.RDS" \) -print0 | 
while IFS= read -r -d '' file; do
    
	ext=`basename "$file" | sed -E 's/.+\.//'`
	isgrange="FALSE"
	
	if [[ "${ext}" == "RDS" ]]; then
		isgrange="TRUE"
	fi
	
	SummarizedExp_MQ_preparation.R \
		-i "${file}" \
		-g "${isgrange}" \
		-o ${WKDIR}
	
done

# Time scripts ends
echo "Completed $(date +"%T %d-%m-%Y")."
