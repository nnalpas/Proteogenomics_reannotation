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
module load curl/7.57.0
module load zlib/1.2.11
module load blast+/2.6.0
module load math/R/3.2.3-mkl-11.3
module load clustal_omega/1.2.4
module load emboss/6.6.0

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
	#SIXFRAMEPROT=`echo $GENOME | perl -p -e 's/^(.*\\/)(.*)\\.fasta/$1Find0_$2_FIXED.fasta/'`
	#SIXFRAMEGENE=`echo $GENOME | perl -p -e 's/^(.*\\/)(.*)\\.fasta/$1Find2_$2_FIXED.fasta/'`

fi



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
    
    ${PBS_O_HOME}/bin/MakeBlastDb.sh -i ${InputType} -y "prot" -t ${TaxId} ${UNIREFPROT} ${SIXFRAMEPROT} ${COORDSIXFRAMEPROT} > ${LogDir}/MakeBlastDb.log 2>&1
    
fi



# Check whether to create a blast database for genome data
if [ $MakeBlastDbNuc == 1 ]; then
    
	${PBS_O_HOME}/bin/MakeBlastDb.sh -i ${InputType} -y "nucl" -t ${TaxId} ${UNIREFGENE} ${SIXFRAMEGENE} ${GENOME} >> ${LogDir}/MakeBlastDb.log 2>&1

fi



###################################
# Blast entries against databases #
###################################

# Check whether to blast the ORFs against taxon reference protein and RNA and against all uniprot and ncbi
if [ $BlastDbBlasting == 1 ]; then
    
    ${PBS_O_HOME}/bin/BlastDbBlasting.sh -o ${ProjDir}/Blast -y "prot" -l "all" -a "blastp" -s ${SIXFRAMEPROT} -q ${UNIREFPROT} -e ${Eval} -n ${NumAlign} -t ${THREADS} -b "ORFprot_vs_Refprot" > ${LogDir}/BlastDbBlasting.log 2>&1
    ${PBS_O_HOME}/bin/BlastDbBlasting.sh -o ${ProjDir}/Blast -y "prot" -l ${ProjDir}/Novel_res/Novel_ORF.txt -a "blastp" -s ${SIXFRAMEPROT} -q ${ALLUNIPROT} -e ${Eval} -n ${NumAlign} -t ${THREADS} -b "ORFprot_vs_Uniprot" >> ${LogDir}/BlastDbBlasting.log 2>&1
    ${PBS_O_HOME}/bin/BlastDbBlasting.sh -o ${ProjDir}/Blast -y "prot" -l ${ProjDir}/Novel_res/Novel_ORF.txt -a "blastp" -s ${SIXFRAMEPROT} -q ${ALLNCBIPROT} -e ${Eval} -n ${NumAlign} -t ${THREADS} -b "ORFprot_vs_NCBIprot" >> ${LogDir}/BlastDbBlasting.log 2>&1
    #${PBS_O_HOME}/bin/BlastDbBlasting.sh -o ${ProjDir}/Blast -y "nucl" -l "all" -a "blastn" -s ${SIXFRAMEGENE} -q ${UNIREFGENE} -e ${Eval} -n ${NumAlign} -t ${THREADS} -b "ORFnucl_vs_Refrna" >> ${LogDir}/BlastDbBlasting.log 2>&1
    #${PBS_O_HOME}/bin/BlastDbBlasting.sh -o ${ProjDir}/Blast -y "nucl" -l ${ProjDir}/Novel_res/Novel_ORF.txt -a "blastn" -s ${SIXFRAMEGENE} -q ${ALLNCBIRNA} -e ${Eval} -n ${NumAlign} -t ${THREADS} -b "ORFnucl_vs_NCBIrna" >> ${LogDir}/BlastDbBlasting.log 2>&1
    ${PBS_O_HOME}/bin/BlastDbBlasting.sh -o ${ProjDir}/Blast -y "prot" -l "all" -a "tblastn" -s ${UNIREFPROT} -q ${GENOME} -e ${Eval} -n ${NumAlign} -m "${TblastnParam}" -t 1 -b "Refprot_vs_Genome" >> ${LogDir}/BlastDbBlasting.log 2>&1
    ${PBS_O_HOME}/bin/BlastDbBlasting.sh -o ${ProjDir}/Blast -y "prot" -l "all" -a "tblastn" -s ${SIXFRAMEPROT} -q ${GENOME} -e ${Eval} -n ${NumAlign} -m "${TblastnParam}" -t 1 -b "ORFprot_vs_Genome" >> ${LogDir}/BlastDbBlasting.log 2>&1
    
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

    ${PBS_O_HOME}/bin/Genomic_position_from_blast.R -f ${UNIREFPROT} -b ${ProjDir}/Blast/Refprot_vs_Genome -g ${GENOME} -o ${ProjDir}/ProtPosition/Ref_prot_coordinates.txt > ${LogDir}/ProteinCoordinate.log 2>&1
	${PBS_O_HOME}/bin/Genomic_position_from_blast.R -f ${SIXFRAMEPROT} -b ${ProjDir}/Blast/ORFprot_vs_Genome -g ${GENOME} -o ${ProjDir}/ProtPosition/Orf_prot_coordinates.txt >> ${LogDir}/ProteinCoordinate.log 2>&1
	
fi



######################
# Protein annotation #
######################

# Check whether to perform protein annotation
if [ $ProteinAnnotation == 1 ]; then

    ${PBS_O_HOME}/bin/Fasta_annotation.R -f ${UNIREFPROT} -t ${TaxId} -c "ENTRY-NAME,GENES,PROTEIN-NAMES,SUBCELLULAR-LOCATIONS,FAMILIES" -k "UNIPROTKB" -o ${ProjDir}/ProtAnnotation/Ref_prot_annotations.txt > ${LogDir}/ProteinAnnotation.log 2>&1
	${PBS_O_HOME}/bin/Transfer_annotation.R -c ${ProjDir}/ReciprocalBlast/Best_Reciproc_Blast_cross-map_Refprot_vs_ORFprot -a ${ProjDir}/ProtAnnotation/Ref_prot_annotations.txt -k "UNIPROTKB" -r "sseqid" -o ${ProjDir}/ProtAnnotation/Orf_prot_annotations.txt >> ${LogDir}/ProteinAnnotation.log 2>&1
	
fi



#######################
# BSgenome generation #
#######################

# Check whether to perform BSgenome package installation
if [ $BSgenomeForge == 1 ]; then

    ${PBS_O_HOME}/bin/BSgenomeForge.sh -o ${ProjDir}/BSgenome -s ${SEED} ${GENOME} > ${LogDir}/BSgenomeForge.log 2>&1
	
fi



##########################
# GRange object creation #
##########################

# Check whether to perform GRange object creation
if [ $GRangeCreate == 1 ]; then
	
	${PBS_O_HOME}/bin/GRanges_generation.R -g ${GENOME} -n ${GenomeName} -t ${Circular} -o ${ProjDir}/GRanges/Genome_grange.RDS > ${LogDir}/GRangesGeneration.log 2>&1
    ${PBS_O_HOME}/bin/GRanges_generation.R -c ${ProjDir}/ProtPosition/Ref_prot_coordinates.txt -g ${GENOME} -n ${GenomeName} -t ${Circular} -a ${ProjDir}/ProtAnnotation/Ref_prot_annotations.txt -o ${ProjDir}/GRanges/Ref_prot_grange.RDS >> ${LogDir}/GRangesGeneration.log 2>&1
	${PBS_O_HOME}/bin/GRanges_generation.R -c ${ProjDir}/ProtPosition/Orf_prot_coordinates.txt -g ${GENOME} -n ${GenomeName} -t ${Circular} -a ${ProjDir}/ProtAnnotation/Orf_prot_annotations.txt -o ${ProjDir}/GRanges/Orf_prot_grange.RDS >> ${LogDir}/GRangesGeneration.log 2>&1
	
fi



##################################
# Peptide coordinate computation #
##################################

# Check whether to perform peptide coordinate identification
if [ $PeptideCoordinate == 1 ]; then

    ${PBS_O_HOME}/bin/Genomic_position_for_peptides.R -p ${ProjDir}/Novel_res/Peptides_location.RDS -k ${ProjDir}/ProtPosition/Ref_prot_coordinates.txt -n ${ProjDir}/ProtPosition/Orf_prot_coordinates.txt -o ${ProjDir}/PeptPosition/Pept_seq_coordinates.txt > ${LogDir}/PeptideCoordinate.log 2>&1
	${PBS_O_HOME}/bin/GRanges_generation.R -c ${ProjDir}/PeptPosition/Pept_seq_coordinates.txt -g ${GENOME} -n ${GenomeName} -t ${Circular} -o ${ProjDir}/GRanges/Pept_seq_grange.RDS >> ${LogDir}/PeptideCoordinate.log 2>&1
	
fi



#################################
# Operon coordinate computation #
#################################

# Check whether to perform operon coordinate identification
if [ $OperonCoordinate == 1 ]; then

    ${PBS_O_HOME}/bin/Operon_reformatting.R -i ${OPERON} -o ${ProjDir}/ProtPosition/Operon_coordinates.txt > ${LogDir}/OperonCoordinate.log 2>&1
	${PBS_O_HOME}/bin/GRanges_generation.R -c ${ProjDir}/ProtPosition/Operon_coordinates.txt -g ${GENOME} -n ${GenomeName} -t ${Circular} -o ${ProjDir}/GRanges/Operon_grange.RDS >> ${LogDir}/OperonCoordinate.log 2>&1
	
fi



################################
# Sanger validation coordinate #
################################

# Check whether to compute coordinate of Sanger sequence
if [ $SangerCoordinate == 1 ]; then
    
    ${PBS_O_HOME}/bin/MakeBlastDb.sh -i ${InputType} -y "nucl" -t ${TaxId} ${SANGER} > ${LogDir}/SangerCoordinate.log 2>&1
    ${PBS_O_HOME}/bin/BlastDbBlasting.sh -o ${ProjDir}/Sanger_validation -y "nucl" -l "all" -a "blastn" -s ${SANGER} -q ${GENOME} -e ${Eval} -n ${NumAlign} -t ${THREADS} -b "Sanger_vs_Genome" >> ${LogDir}/SangerCoordinate.log 2>&1
	${PBS_O_HOME}/bin/Genomic_position_from_blast.R -f ${SANGER} -b ${ProjDir}/Sanger_validation/Sanger_vs_Genome -g ${GENOME} -o ${ProjDir}/Sanger_validation/Sanger_seq_coordinates.txt >> ${LogDir}/SangerCoordinate.log 2>&1
	${PBS_O_HOME}/bin/GRanges_generation.R -c ${ProjDir}/Sanger_validation/Sanger_seq_coordinates.txt -g ${GENOME} -n ${GenomeName} -t ${Circular} -o ${ProjDir}/GRanges/Sanger_seq_grange.RDS >> ${LogDir}/SangerCoordinate.log 2>&1
	
fi



##############################
# Novelty reason explanation #
##############################

# Check whether to perform the novelty reason explanation
if [ $NoveltyReason == 1 ]; then
	
	${PBS_O_HOME}/bin/Novelty_discovery_peptides.R -e ${ProjDir}/Novel_res/Sequence_group_mapping.RDS -f ${UNIREFPROT} -r ${ProjDir}/ReciprocalBlast/Best_Reciproc_Blast_Refprot_vs_ORFprot -u ${ProjDir}/ReciprocalBlast/Best_Reciproc_Blast_Uniprot_vs_ORFprot -n ${ProjDir}/ReciprocalBlast/Best_Reciproc_Blast_NCBIprot_vs_ORFprot -p ${ProjDir}/Novel_res/Peptides_location.RDS -t ${THREADS} -o ${ProjDir}/NoveltyExplain > ${LogDir}/NoveltyReason.log 2>&1
	${PBS_O_HOME}/bin/Novelty_discovery_ORFs.R -i ${ProjDir}/NoveltyExplain/Sequence_novelty_reason.RDS -r ${ProjDir}/GRanges/Ref_prot_grange.RDS -n ${ProjDir}/GRanges/Orf_prot_grange.RDS -p ${ProjDir}/GRanges/Operon_grange.RDS -e ${ProjDir}/Novel_res/Peptides_location.RDS -d ${ProjDir}/GRanges/Pept_seq_grange.RDS -s ${ProjDir}/GRanges/Sanger_seq_grange.RDS -g ${ProjDir}/GRanges/Genome_grange.RDS -a "${AddRBS}" -t ${THREADS} -o ${ProjDir}/NoveltyExplain >> ${LogDir}/NoveltyReason.log 2>&1

fi



##############################
# Clustal multiple alignment #
##############################

# Check whether to perform multiple pair-wise alignments for RBS analysis
if [ $ClustalAlign == 1 ]; then

    ${PBS_O_HOME}/bin/ClustalAlignment.sh -o ${ProjDir}/RBS_clustalo -t ${THREADS} ${RBSLOW}.fasta ${RBSHIGH}.fasta > ${LogDir}/ClustalAlign.log 2>&1
    ${PBS_O_HOME}/bin/Seqlogo.sh -o ${ProjDir}/RBS_clustalo ${RBSLOW}_alignment_full.clustal ${RBSHIGH}_alignment_full.clustal >> ${LogDir}/ClustalAlign.log 2>&1
    
fi


