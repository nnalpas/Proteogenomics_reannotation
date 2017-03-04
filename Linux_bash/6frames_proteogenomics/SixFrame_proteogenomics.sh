#!/bin/bash

# Genomic location of novel ORF



##############################
# Parameters and logs set-up #
##############################

# Check if $SCRIPT_FLAGS is set
if [ -n "${SCRIPT_FLAGS}" ] ; then

	# But if positional parameters are already present
	# we are going to ignore $SCRIPT_FLAGS
	if [ -z "${*}" ] ; then
		set -- ${SCRIPT_FLAGS}
	fi

fi

# Check whether a parameter file was provided
#if [ $# -ne 1 ]; then
#	echo "Usage: $0 <ParametersPath>"
#	exit 1
#fi

# Load all parameters
source ${SCRIPT_FLAGS}
export OMP_NUM_THREADS=${MOAB_PROCCOUNT}

# Go to workspace
cd ${PBS_O_WORKDIR}

# Create project directory
ProjDir=${PBS_O_INITDIR}/${ProjectName}

# Create the log folder
DateStart=$(date -I)
LogDir=${ProjDir}/Log/${DateStart}
mkdir -p $LogDir

# Copy the parameter file to log
cp ${SCRIPT_FLAGS} ${LogDir}/Parameters_${DateStart}.txt



############
# Get ORFs #
############

# Check whether to get the ORFs from genomic sequence
if [ $GetOrf == 1 ]; then

    ${PBS_O_HOME}/bin/GetOrf.sh ${ProjDir}/Nuc_translation 0 ${TABLE} ${MINSIZE} ${CIRCULAR} ${GENOME} > ${LogDir}/GetOrf.log 2>&1
    ${PBS_O_HOME}/bin/GetOrf.sh ${ProjDir}/Nuc_translation 2 ${TABLE} ${MINSIZE} ${CIRCULAR} ${GENOME} >> ${LogDir}/GetOrf.log 2>&1
    
fi

# Define the file names for the ORF of the genomic sequence
SIXFRAMEPROT=`echo $GENOME | perl -p -e 's/^(.*\\/)(.*)\\.fasta/$1Find0_$2_FIXED.fasta/'`
SIXFRAMEGENE=`echo $GENOME | perl -p -e 's/^(.*\\/)(.*)\\.fasta/$1Find2_$2_FIXED.fasta/'`



##################
# Get novel ORFs #
##################

# Check whether to get the novel ORFs from MaxQuant seach results
if [ $GetNovelEntries == 1 ]; then

    ${PBS_O_HOME}/bin/GetNovelEntries.R ${ProjDir}/MaxQuant_txt "." ${UNIREFPROT} ${SIXFRAMEPROT} > ${LogDir}/GetNovelEntries.log 2>&1

fi



##############################
# Make custom Blast database #
##############################

# Check whether to create a blast database for proteome data
if [ $MakeBlastDbProt == 1 ]; then
    
    ${PBS_O_HOME}/bin/MakeBlastDb.sh ${InputType} "prot" ${TaxId} ${UNIREFPROT} ${SIXFRAMEPROT} > ${LogDir}/MakeBlastDb.log 2>&1
    
fi



# Check whether to create a blast database for genome data
if [ $MakeBlastDbNuc == 1 ]; then
    
    ${PBS_O_HOME}/bin/MakeBlastDb.sh ${InputType} "nucl" ${TaxId} ${UNIREFGENE} ${SIXFRAMEGENE} >> ${LogDir}/MakeBlastDb.log 2>&1
    
fi



###################################
# Blast entries against databases #
###################################

# Check whether to blast the ORFs against taxon reference protein and RNA and against all uniprot and ncbi
if [ $BlastDbBlasting == 1 ]; then
    
    ${PBS_O_HOME}/bin/BlastDbBlasting.sh ${ProjDir}/Blast "prot" "all" "blastp" ${SIXFRAMEPROT} ${UNIREFPROT} ${Eval} ${NumAlign} ${OMP_NUM_THREADS} "ORFprot_vs_Refprot" > ${LogDir}/BlastDbBlasting.log 2>&1
    ${PBS_O_HOME}/bin/BlastDbBlasting.sh ${ProjDir}/Blast "prot" ${ProjDir}/MaxQuant_txt/Novel_ORF.txt "blastp" ${SIXFRAMEPROT} ${ALLUNIPROT} ${Eval} ${NumAlign} ${OMP_NUM_THREADS} "ORFprot_vs_Uniprot" >> ${LogDir}/BlastDbBlasting.log 2>&1
    ${PBS_O_HOME}/bin/BlastDbBlasting.sh ${ProjDir}/Blast "prot" ${ProjDir}/MaxQuant_txt/Novel_ORF.txt "blastp" ${SIXFRAMEPROT} ${ALLNCBIPROT} ${Eval} ${NumAlign} ${OMP_NUM_THREADS} "ORFprot_vs_NCBIprot" >> ${LogDir}/BlastDbBlasting.log 2>&1
    ${PBS_O_HOME}/bin/BlastDbBlasting.sh ${ProjDir}/Blast "nucl" "all" "blastn" ${SIXFRAMEGENE} ${UNIREFGENE} ${Eval} ${NumAlign} ${OMP_NUM_THREADS} "ORFnucl_vs_Refrna" >> ${LogDir}/BlastDbBlasting.log 2>&1
    ${PBS_O_HOME}/bin/BlastDbBlasting.sh ${ProjDir}/Blast "nucl" ${ProjDir}/MaxQuant_txt/Novel_ORF.txt "blastn" ${SIXFRAMEGENE} ${ALLNCBIRNA} ${Eval} ${NumAlign} ${OMP_NUM_THREADS} "ORFnucl_vs_NCBIrna" >> ${LogDir}/BlastDbBlasting.log 2>&1
    
fi


