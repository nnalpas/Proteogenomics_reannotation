#!/bin/bash

# Genomic location of novel ORF



##############################
# Parameters and logs set-up #
##############################

# Check whether a parameter file was provided
if [ $# -ne 1 ]; then
	echo "Usage: $0 <ParametersPath>"
	exit 1
fi

# Load all parameters
Param=$1
source $Param

# Create project directory
ProjDir=${MainDir}/${ProjectName}

# Create the log folder
DateStart=$(date -I)
LogDir=${ProjDir}/Log/${DateStart}
mkdir -p $LogDir

# Copy the parameter file to log
cp $Param ${LogDir}/Parameters.txt



############
# Get ORFs #
############

# Check whether to get the ORFs from genomic sequence
if [ $GetOrf == 1 ]; then

    ${HOME}/bin/GetOrf.sh ${ProjDir}/Nuc_translation ${FIND} ${TABLE} ${MINSIZE} ${CIRCULAR} ${SEQUENCE} > ${LogDir}/GetOrf.log 2>&1
    
fi


