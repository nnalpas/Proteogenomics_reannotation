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
    
    ${PBS_O_HOME}/bin/BlastDbBlasting.sh ${ProjDir}/Blast "prot" "all" "blastp" ${SIXFRAMEPROT} ${UNIREFPROT} ${Eval} ${NumAlign} ${THREADS} "ORFprot_vs_Refprot" > ${LogDir}/BlastDbBlasting.log 2>&1
    ${PBS_O_HOME}/bin/BlastDbBlasting.sh ${ProjDir}/Blast "prot" ${ProjDir}/MaxQuant_txt/Novel_ORF.txt "blastp" ${SIXFRAMEPROT} ${ALLUNIPROT} ${Eval} ${NumAlign} ${THREADS} "ORFprot_vs_Uniprot" >> ${LogDir}/BlastDbBlasting.log 2>&1
    ${PBS_O_HOME}/bin/BlastDbBlasting.sh ${ProjDir}/Blast "prot" ${ProjDir}/MaxQuant_txt/Novel_ORF.txt "blastp" ${SIXFRAMEPROT} ${ALLNCBIPROT} ${Eval} ${NumAlign} ${THREADS} "ORFprot_vs_NCBIprot" >> ${LogDir}/BlastDbBlasting.log 2>&1
    ${PBS_O_HOME}/bin/BlastDbBlasting.sh ${ProjDir}/Blast "nucl" "all" "blastn" ${SIXFRAMEGENE} ${UNIREFGENE} ${Eval} ${NumAlign} ${THREADS} "ORFnucl_vs_Refrna" >> ${LogDir}/BlastDbBlasting.log 2>&1
    ${PBS_O_HOME}/bin/BlastDbBlasting.sh ${ProjDir}/Blast "nucl" ${ProjDir}/MaxQuant_txt/Novel_ORF.txt "blastn" ${SIXFRAMEGENE} ${ALLNCBIRNA} ${Eval} ${NumAlign} ${THREADS} "ORFnucl_vs_NCBIrna" >> ${LogDir}/BlastDbBlasting.log 2>&1
    
fi



##################################
# Best blast hits identification #
##################################

# Check whether to identify the best blast hits from all results of previous step
if [ $BestBlast ]; then

    ${PBS_O_HOME}/bin/Best_blast.R -i ${ProjDir}/Blast/ORFprot_vs_Refprot -o ${ProjDir}/Blast/ORFprot_vs_Refprot_besthit.txt > ${LogDir}/BestBlast.log 2>&1
    ${PBS_O_HOME}/bin/Best_blast.R -i ${ProjDir}/Blast/ORFprot_vs_Uniprot -o ${ProjDir}/Blast/ORFprot_vs_Uniprot_besthit.txt >> ${LogDir}/BestBlast.log 2>&1
    ${PBS_O_HOME}/bin/Best_blast.R -i ${ProjDir}/Blast/ORFprot_vs_NCBIprot -o ${ProjDir}/Blast/ORFprot_vs_NCBIprot_besthit.txt >> ${LogDir}/BestBlast.log 2>&1
    ${PBS_O_HOME}/bin/Best_blast.R -i ${ProjDir}/Blast/ORFnucl_vs_Refrna -o ${ProjDir}/Blast/ORFnucl_vs_Refrna_besthit.txt >> ${LogDir}/BestBlast.log 2>&1
    ${PBS_O_HOME}/bin/Best_blast.R -i ${ProjDir}/Blast/ORFnucl_vs_NCBIrna -o ${ProjDir}/Blast/ORFnucl_vs_NCBIrna_besthit.txt >> ${LogDir}/BestBlast.log 2>&1

fi



#########################
# Reciprocal best blast #
#########################

# Check whether to perform the reciprocal best blast against previously used databases
if [ $ReciprocalBlast ]; then
    
    ${PBS_O_HOME}/bin/BlastDbBlasting.sh ${ProjDir}/Blast "prot" ${ProjDir}/Blast/ORFprot_vs_Refprot_besthit.txt "blastp" ${UNIREFPROT} ${SIXFRAMEPROT} ${Eval} ${NumAlign} ${THREADS} "Refprot_vs_ORFprot" > ${LogDir}/ReciprocalBlast.log 2>&1
    ${PBS_O_HOME}/bin/BlastDbBlasting.sh ${ProjDir}/Blast "prot" ${ProjDir}/Blast/ORFprot_vs_Uniprot_besthit.txt "blastp" ${ALLUNIPROT} ${SIXFRAMEPROT} ${Eval} ${NumAlign} ${THREADS} "Uniprot_vs_ORFprot" >> ${LogDir}/ReciprocalBlast.log 2>&1
    ${PBS_O_HOME}/bin/BlastDbBlasting.sh ${ProjDir}/Blast "prot" ${ProjDir}/Blast/ORFprot_vs_NCBIprot_besthit.txt "blastp" ${ALLNCBIPROT} ${SIXFRAMEPROT} ${Eval} ${NumAlign} ${THREADS} "NCBIprot_vs_ORFprot" >> ${LogDir}/ReciprocalBlast.log 2>&1
    ${PBS_O_HOME}/bin/BlastDbBlasting.sh ${ProjDir}/Blast "nucl" ${ProjDir}/Blast/ORFnucl_vs_Refrna_besthit.txt "blastn" ${UNIREFGENE} ${SIXFRAMEGENE} ${Eval} ${NumAlign} ${THREADS} "Refrna_vs_ORFnucl" >> ${LogDir}/ReciprocalBlast.log 2>&1
    ${PBS_O_HOME}/bin/BlastDbBlasting.sh ${ProjDir}/Blast "nucl" ${ProjDir}/Blast/ORFnucl_vs_NCBIrna_besthit.txt "blastn" ${ALLNCBIRNA} ${SIXFRAMEGENE} ${Eval} ${NumAlign} ${THREADS} "NCBIrna_vs_ORFnucl" >> ${LogDir}/ReciprocalBlast.log 2>&1

fi



######################################
# Compile reciprocal best blast data #
######################################

# Check whether to perform reciprocal best blast hit data processing
if [ $ReciprocalBestBlast ]; then

    ${PBS_O_HOME}/bin/Reciprocal_best_blast.R -b ${ProjDir}/Blast/ORFprot_vs_Refprot -r ${ProjDir}/Blast/Refprot_vs_ORFprot -o ${ProjDir}/Blast/ORFprot_vs_Refprot_reciprocbesthit.txt > ${LogDir}/ReciprocalBestBlast.log 2>&1
    ${PBS_O_HOME}/bin/Reciprocal_best_blast.R -b ${ProjDir}/Blast/ORFprot_vs_Uniprot -r ${ProjDir}/Blast/Uniprot_vs_ORFprot -o ${ProjDir}/Blast/ORFprot_vs_Uniprot_reciprocbesthit.txt >> ${LogDir}/ReciprocalBestBlast.log 2>&1
    ${PBS_O_HOME}/bin/Reciprocal_best_blast.R -b ${ProjDir}/Blast/ORFprot_vs_NCBIprot -r ${ProjDir}/Blast/NCBIprot_vs_ORFprot -o ${ProjDir}/Blast/ORFprot_vs_NCBIprot_reciprocbesthit.txt >> ${LogDir}/ReciprocalBestBlast.log 2>&1
    ${PBS_O_HOME}/bin/Reciprocal_best_blast.R -b ${ProjDir}/Blast/ORFnucl_vs_Refrna -r ${ProjDir}/Blast/Refrna_vs_ORFnucl -o ${ProjDir}/Blast/ORFnucl_vs_Refrna_reciprocbesthit.txt >> ${LogDir}/ReciprocalBestBlast.log 2>&1
    ${PBS_O_HOME}/bin/Reciprocal_best_blast.R -b ${ProjDir}/Blast/ORFnucl_vs_NCBIrna -r ${ProjDir}/Blast/NCBIrna_vs_ORFnucl -o ${ProjDir}/Blast/ORFnucl_vs_NCBIrna_reciprocbesthit.txt >> ${LogDir}/ReciprocalBestBlast.log 2>&1

fi



##################################
# Protein coordinate computation #
##################################

# Check whether to perform protein coordinate identification
if [ $ProteinCoordinate ]; then

    ${PBS_O_HOME}/bin/BlastDbBlasting.sh ${ProjDir}/ProtPosition "prot" "all" "blastp" ${UNIREFPROT} ${SIXFRAMEPROT} ${Eval} ${NumAlign} ${THREADS} "AllRefprot_vs_ORFprot" > ${LogDir}/ProteinCoordinate.log 2>&1
    ${PBS_O_HOME}/bin/ORF_coordinates.R -f ${SIXFRAMEPROT} -o ${ProjDir}/ProtPosition/orf_coordinates.txt >> ${LogDir}/ProteinCoordinate.log 2>&1
    ${PBS_O_HOME}/bin/Genomic_position_from_blast.R -f ${UNIREFPROT} -b ${ProjDir}/Blast/AllRefprot_vs_ORFprot -g ${ProjDir}/ProtPosition/orf_coordinates.txt -o ${ProjDir}/ProtPosition/refprot_coordinates.txt >> ${LogDir}/ProteinCoordinate.log 2>&1

fi


