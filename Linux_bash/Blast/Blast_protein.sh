#!/bin/bash

# Blast pipeline for protein sequence



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



##############################
# Make custom Blast database #
##############################

# Check whether to create a blast database
if [ $MakeBlastDb == 1 ]; then

    ${HOME}/bin/MakeBlastDb.sh ${RefProteome} ${frameProteome} > ${LogDir}/MakeBlastDb.log 2>&1
    
fi


