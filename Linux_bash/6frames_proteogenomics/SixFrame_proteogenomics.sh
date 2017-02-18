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

    ${HOME}/bin/GetOrf.sh ${ProjDir}/Nuc_translation 0 ${TABLE} ${MINSIZE} ${CIRCULAR} ${GENOME} > ${LogDir}/GetOrf.log 2>&1
    ${HOME}/bin/GetOrf.sh ${ProjDir}/Nuc_translation 2 ${TABLE} ${MINSIZE} ${CIRCULAR} ${GENOME} >> ${LogDir}/GetOrf.log 2>&1
    
fi



##############################
# Make custom Blast database #
##############################

# Define the file names for the ORF of the genomic sequence
SIXFRAMEPROT=`echo $GENOME | perl -p -e 's/^(.*\\/)(.*)\\.fasta/$1Find0_$2_FIXED.fasta/'`
SIXFRAMEGENE=`echo $GENOME | perl -p -e 's/^(.*\\/)(.*)\\.fasta/$1Find2_$2_FIXED.fasta/'`

# Check whether to create a blast database for proteome data
if [ $MakeBlastDbProt == 1 ]; then
    
    ${HOME}/bin/MakeBlastDb.sh ${InputType} "prot" ${TaxId} ${UNIREFPROT} ${SIXFRAMEPROT} > ${LogDir}/MakeBlastDb.log 2>&1
    
fi



# Check whether to create a blast database for genome data
if [ $MakeBlastDbNuc == 1 ]; then
    
    ${HOME}/bin/MakeBlastDb.sh ${InputType} "nucl" ${TaxId} ${UNIREFGENE} ${SIXFRAMEGENE} >> ${LogDir}/MakeBlastDb.log 2>&1
    
fi



###################################
# Blast entries against databases #
###################################

# Check whether to blast the ORFs against taxon reference protein and RNA and against all uniprot and ncbi
if [ $BlastDbBlasting == 1 ]; then
    
    ${HOME}/bin/BlastDbBlasting.sh ${ProjDir}/Blast "prot" "all" "blastp" ${SIXFRAMEPROT} ${UNIREFPROT} ${Eval} ${NumAlign} ${MultiThreads} "ORFprot_vs_Refprot" > ${LogDir}/BlastDbBlasting.log 2>&1
    ${HOME}/bin/BlastDbBlasting.sh ${ProjDir}/Blast "prot" "all" "blastp" ${SIXFRAMEPROT} ${ALLUNIPROT} ${Eval} ${NumAlign} ${MultiThreads} "ORFprot_vs_Uniprot" >> ${LogDir}/BlastDbBlasting.log 2>&1
    ${HOME}/bin/BlastDbBlasting.sh ${ProjDir}/Blast "prot" "all" "blastp" ${SIXFRAMEPROT} ${ALLNCBIPROT} ${Eval} ${NumAlign} ${MultiThreads} "ORFprot_vs_NCBIprot" >> ${LogDir}/BlastDbBlasting.log 2>&1
    ${HOME}/bin/BlastDbBlasting.sh ${ProjDir}/Blast "prot" "all" "blastp" ${SIXFRAMEGENE} ${UNIREFGENE} ${Eval} ${NumAlign} ${MultiThreads} "ORFnucl_vs_Refrna" >> ${LogDir}/BlastDbBlasting.log 2>&1
    ${HOME}/bin/BlastDbBlasting.sh ${ProjDir}/Blast "prot" "all" "blastp" ${SIXFRAMEGENE} ${ALLNCBIRNA} ${Eval} ${NumAlign} ${MultiThreads} "ORFnucl_vs_NCBIrna" >> ${LogDir}/BlastDbBlasting.log 2>&1
    
fi


