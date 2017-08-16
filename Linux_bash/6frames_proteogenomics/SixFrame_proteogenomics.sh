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

# Go to workspace
cd ${PBS_O_WORKDIR}

# Load the required modules
module load ${PBS_O_HOME}/modulefiles/blast+/2.6.0
module load math/R/3.2.3-mkl-11.3

# Create project directory
ProjDir=${PBS_O_INITDIR}/${ProjectName}

# Create the log folder
DateStart=$(date -I)
LogDir=${ProjDir}/Log/${DateStart}
mkdir -p $LogDir

# Copy the parameter file to log
cp ${SCRIPT_FLAGS} ${LogDir}/Parameters_${DateStart}.txt

# Get the number of threads to use
THREADS=`wc -l ${PBS_NODEFILE} | grep -o "^[0-9]*"`



############
# Get ORFs #
############

# Check whether to get the ORFs from genomic sequence
if [ $GetOrf == 1 ]; then

    ${PBS_O_HOME}/bin/GetOrf.sh ${ProjDir}/Nuc_translation 0 ${TABLE} ${MINSIZE} ${CIRCULAR} ${GENOME} > ${LogDir}/GetOrf.log 2>&1
    ${PBS_O_HOME}/bin/GetOrf.sh ${ProjDir}/Nuc_translation 2 ${TABLE} ${MINSIZE} ${CIRCULAR} ${GENOME} >> ${LogDir}/GetOrf.log 2>&1
    
fi

# Define the file names for the ORF of the genomic sequence
#SIXFRAMEPROT=`echo $GENOME | perl -p -e 's/^(.*\\/)(.*)\\.fasta/$1Find0_$2_FIXED.fasta/'`
#SIXFRAMEGENE=`echo $GENOME | perl -p -e 's/^(.*\\/)(.*)\\.fasta/$1Find2_$2_FIXED.fasta/'`



##################
# Get novel ORFs #
##################

# Check whether to get the novel ORFs from MaxQuant seach results
if [ $GetNovelEntries == 1 ]; then

    ${PBS_O_HOME}/bin/GetNovelEntries.R -o ${ProjDir}/Novel_res -m ${ProjDir}/MaxQuant_txt -r ${UNIREFPROT} -n ${SIXFRAMEPROT} > ${LogDir}/GetNovelEntries.log 2>&1

fi



##############################
# Make custom Blast database #
##############################

# Check whether to create a blast database for proteome data
if [ $MakeBlastDbProt == 1 ]; then
    
    ${PBS_O_HOME}/bin/MakeBlastDb.sh -i ${InputType} -y "prot" -t ${TaxId} ${UNIREFPROT} ${SIXFRAMEPROT} > ${LogDir}/MakeBlastDb.log 2>&1
    
fi



# Check whether to create a blast database for genome data
if [ $MakeBlastDbNuc == 1 ]; then
    
    ${PBS_O_HOME}/bin/MakeBlastDb.sh -i ${InputType} -y "nucl" -t ${TaxId} ${UNIREFGENE} ${SIXFRAMEGENE} >> ${LogDir}/MakeBlastDb.log 2>&1
    
fi



###################################
# Blast entries against databases #
###################################

# Check whether to blast the ORFs against taxon reference protein and RNA and against all uniprot and ncbi
if [ $BlastDbBlasting == 1 ]; then
    
    ${PBS_O_HOME}/bin/BlastDbBlasting.sh -o ${ProjDir}/Blast -y "prot" -l "all" -a "blastp" -s ${SIXFRAMEPROT} -q ${UNIREFPROT} -e ${Eval} -n ${NumAlign} -t ${THREADS} -b "ORFprot_vs_Refprot" > ${LogDir}/BlastDbBlasting.log 2>&1
    ${PBS_O_HOME}/bin/BlastDbBlasting.sh -o ${ProjDir}/Blast -y "prot" -l ${ProjDir}/Novel_res/*_Novel_ORF.txt -a "blastp" -s ${SIXFRAMEPROT} -q ${ALLUNIPROT} -e ${Eval} -n ${NumAlign} -t ${THREADS} -b "ORFprot_vs_Uniprot" >> ${LogDir}/BlastDbBlasting.log 2>&1
    ${PBS_O_HOME}/bin/BlastDbBlasting.sh -o ${ProjDir}/Blast -y "prot" -l ${ProjDir}/Novel_res/*_Novel_ORF.txt -a "blastp" -s ${SIXFRAMEPROT} -q ${ALLNCBIPROT} -e ${Eval} -n ${NumAlign} -t ${THREADS} -b "ORFprot_vs_NCBIprot" >> ${LogDir}/BlastDbBlasting.log 2>&1
    #${PBS_O_HOME}/bin/BlastDbBlasting.sh -o ${ProjDir}/Blast -y "nucl" -l "all" -a "blastn" -s ${SIXFRAMEGENE} -q ${UNIREFGENE} -e ${Eval} -n ${NumAlign} -t ${THREADS} -b "ORFnucl_vs_Refrna" >> ${LogDir}/BlastDbBlasting.log 2>&1
    #${PBS_O_HOME}/bin/BlastDbBlasting.sh -o ${ProjDir}/Blast -y "nucl" -l ${ProjDir}/Novel_res/Novel_ORF.txt -a "blastn" -s ${SIXFRAMEGENE} -q ${ALLNCBIRNA} -e ${Eval} -n ${NumAlign} -t ${THREADS} -b "ORFnucl_vs_NCBIrna" >> ${LogDir}/BlastDbBlasting.log 2>&1
    
fi



##################################
# Best blast hits identification #
##################################

# Check whether to identify the best blast hits from all results of previous step
if [ $BestBlast == 1 ]; then

    ${PBS_O_HOME}/bin/Best_blast.R -i ${ProjDir}/Blast/ORFprot_vs_Refprot -o ${ProjDir}/Blast > ${LogDir}/BestBlast.log 2>&1
    ${PBS_O_HOME}/bin/Best_blast.R -i ${ProjDir}/Blast/ORFprot_vs_Uniprot -o ${ProjDir}/Blast >> ${LogDir}/BestBlast.log 2>&1
    ${PBS_O_HOME}/bin/Best_blast.R -i ${ProjDir}/Blast/ORFprot_vs_NCBIprot -o ${ProjDir}/Blast >> ${LogDir}/BestBlast.log 2>&1
    #${PBS_O_HOME}/bin/Best_blast.R -i ${ProjDir}/Blast/ORFnucl_vs_Refrna -o ${ProjDir}/Blast/ORFnucl_vs_Refrna_besthit.txt >> ${LogDir}/BestBlast.log 2>&1
    #${PBS_O_HOME}/bin/Best_blast.R -i ${ProjDir}/Blast/ORFnucl_vs_NCBIrna -o ${ProjDir}/Blast/ORFnucl_vs_NCBIrna_besthit.txt >> ${LogDir}/BestBlast.log 2>&1

fi



#########################
# Reciprocal best blast #
#########################

# Check whether to perform the reciprocal best blast against previously used databases
if [ $ReciprocalBlast == 1 ]; then
    
    ${PBS_O_HOME}/bin/BlastDbBlasting.sh -o ${ProjDir}/ReciprocalBlast -y "prot" -l ${ProjDir}/Blast/Reciprocal_id_ORFprot_vs_Refprot -a "blastp" -s ${UNIREFPROT} -q ${SIXFRAMEPROT} -e ${Eval} -n ${NumAlign} -t ${THREADS} -b "Refprot_vs_ORFprot" > ${LogDir}/ReciprocalBlast.log 2>&1
    ${PBS_O_HOME}/bin/BlastDbBlasting.sh -o ${ProjDir}/ReciprocalBlast -y "prot" -l ${ProjDir}/Blast/Reciprocal_id_ORFprot_vs_Uniprot -a "blastp" -s ${ALLUNIPROT} -q ${SIXFRAMEPROT} -e ${Eval} -n ${NumAlign} -t ${THREADS} -b "Uniprot_vs_ORFprot" >> ${LogDir}/ReciprocalBlast.log 2>&1
    ${PBS_O_HOME}/bin/BlastDbBlasting.sh -o ${ProjDir}/ReciprocalBlast -y "prot" -l ${ProjDir}/Blast/Reciprocal_id_ORFprot_vs_NCBIprot -a "blastp" -s ${ALLNCBIPROT} -q ${SIXFRAMEPROT} -e ${Eval} -n ${NumAlign} -t ${THREADS} -b "NCBIprot_vs_ORFprot" >> ${LogDir}/ReciprocalBlast.log 2>&1
    #${PBS_O_HOME}/bin/BlastDbBlasting.sh -o ${ProjDir}/Blast -y "nucl" -l ${ProjDir}/Blast/ORFnucl_vs_Refrna_besthit.txt -a "blastn" -s ${UNIREFGENE} -q ${SIXFRAMEGENE} -e ${Eval} -n ${NumAlign} -t ${THREADS} -b "Refrna_vs_ORFnucl" >> ${LogDir}/ReciprocalBlast.log 2>&1
    #${PBS_O_HOME}/bin/BlastDbBlasting.sh -o ${ProjDir}/Blast -y "nucl" -l ${ProjDir}/Blast/ORFnucl_vs_NCBIrna_besthit.txt -a "blastn" -s ${ALLNCBIRNA} -q ${SIXFRAMEGENE} -e ${Eval} -n ${NumAlign} -t ${THREADS} -b "NCBIrna_vs_ORFnucl" >> ${LogDir}/ReciprocalBlast.log 2>&1

fi



######################################
# Compile reciprocal best blast data #
######################################

# Check whether to perform reciprocal best blast hit data processing
if [ $ReciprocalBestBlast == 1 ]; then

    ${PBS_O_HOME}/bin/Reciprocal_best_blast.R -b ${ProjDir}/Blast/ORFprot_vs_Refprot -r ${ProjDir}/ReciprocalBlast/Refprot_vs_ORFprot -o ${ProjDir}/ReciprocalBlast > ${LogDir}/ReciprocalBestBlast.log 2>&1
    ${PBS_O_HOME}/bin/Reciprocal_best_blast.R -b ${ProjDir}/Blast/ORFprot_vs_Uniprot -r ${ProjDir}/ReciprocalBlast/Uniprot_vs_ORFprot -o ${ProjDir}/ReciprocalBlast >> ${LogDir}/ReciprocalBestBlast.log 2>&1
    ${PBS_O_HOME}/bin/Reciprocal_best_blast.R -b ${ProjDir}/Blast/ORFprot_vs_NCBIprot -r ${ProjDir}/ReciprocalBlast/NCBIprot_vs_ORFprot -o ${ProjDir}/ReciprocalBlast >> ${LogDir}/ReciprocalBestBlast.log 2>&1
    #${PBS_O_HOME}/bin/Reciprocal_best_blast.R -b ${ProjDir}/Blast/ORFnucl_vs_Refrna -r ${ProjDir}/Blast/Refrna_vs_ORFnucl -o ${ProjDir}/Blast/ORFnucl_vs_Refrna_reciprocbesthit.txt >> ${LogDir}/ReciprocalBestBlast.log 2>&1
    #${PBS_O_HOME}/bin/Reciprocal_best_blast.R -b ${ProjDir}/Blast/ORFnucl_vs_NCBIrna -r ${ProjDir}/Blast/NCBIrna_vs_ORFnucl -o ${ProjDir}/Blast/ORFnucl_vs_NCBIrna_reciprocbesthit.txt >> ${LogDir}/ReciprocalBestBlast.log 2>&1

fi



##################################
# Protein coordinate computation #
##################################

# Check whether to perform protein coordinate identification
if [ $ProteinCoordinate == 1 ]; then

    ${PBS_O_HOME}/bin/BlastDbBlasting.sh               -o ${ProjDir}/ProtPosition -y "prot" -l "all" -a "blastp" -s ${UNIREFPROT} -q ${SIXFRAMEPROT} -e ${Eval} -n ${NumAlign} -t ${THREADS} -b "AllRefprot_vs_ORFprot" > ${LogDir}/ProteinCoordinate.log 2>&1
    ${PBS_O_HOME}/bin/BestBlast.sh       
    ${PBS_O_HOME}/bin/ORF_coordinates.R                -f ${SIXFRAMEPROT} -o ${ProjDir}/ProtPosition/orf_coordinates.txt >> ${LogDir}/ProteinCoordinate.log 2>&1
    ${PBS_O_HOME}/bin/ORF_coordinates_transfer.R       -f ${UNIREFPROT} -b ${ProjDir}/ProtPosition/AllRefprot_vs_ORFprot -g ${ProjDir}/ProtPosition/orf_coordinates.txt -o ${ProjDir}/ProtPosition/refprot_coordinates.txt >> ${LogDir}/ProteinCoordinate.log 2>&1

fi


