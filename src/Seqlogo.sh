#!/bin/bash

# Script to generate a sequence logo from alignment file

# Time scripts starts
echo "$0"
echo "Start $(date +"%T %d-%m-%Y")."

# Function holding the usage
display_usage() { 
	echo "
		Usage: $0 [Options] <AlignmentFiles>"
	echo "
		Options:
		-o	[str]	The output directory.
		-h	[]		To display the help.
	" 
}

# Define default values
o_default="`pwd`/seqlogo"

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

# Also get target folders from parameter
TARGETFILES=${@:1}

# Setting default values ${VARIABLE=DEFAULT_VALUE}
# This command will set VARIABLE to DEFAULT_VALUE if it is currently 
# undefined, then return VARIABLE
: ${WKDIR=$o_default}

# Check for mandatory parameters
if [ -z "${TARGETFILES}" ]; then
	echo "At least one alignment file must be specified." >&2
	display_usage
	exit 1
fi

# Create the output folder
if [ ! -d ${WKDIR} ] ; then
	mkdir ${WKDIR}
fi

# Loop through all alignment files provided from command line
for file in $TARGETFILES; do
    
    # Define the basename of the ouptut file
    outfile=`basename $file | sed "s/\\..*//"`
    
    # Run Weblogo to get sequence alignment logo
    kweblogo -seqall $file \
        -goutfile ${WKDIR}/${outfile}.pdf \
        -format 'pdf'
    
done

# Time scripts ends
echo "Completed $(date +"%T %d-%m-%Y")."
