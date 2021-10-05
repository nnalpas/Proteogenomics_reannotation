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
	
	# Get the number of threads to use
	THREADS=`wc -l ${PBS_NODEFILE} | grep -o "^[0-9]*"`

elif [ $# -eq 3 ]; then
	
	SCRIPT_FLAGS=$1
	PBS_O_WORKDIR=$2
	THREADS=$3
	PBS_O_HOME=$HOME

else
	
	echo "Usage: $0 <ParametersPath> <WorkDir> <Threads>"
	exit 1

fi

# Load all parameters
source ${SCRIPT_FLAGS}

# Go to workspace
cd ${PBS_O_WORKDIR}

# Load the required modules
#module load curl/7.57.0
#module load zlib/1.2.11
module load blast+/2.10.1
#module load math/R/3.5.2-mkl-2018 # the prerequisite GNU is incompatible with blast+ (must find solution)
#module load clustal_omega/1.2.4
module load emboss/6.6.0
module load devel/perl/5.26
module load diamond/2.0.2
module load interproscan/5.48-83.0
module load signalp/5.0b

# Create project directory
ProjDir=${PBS_O_INITDIR}/${ProjectName}

# Create the log folder
DateStart=$(date +%F_%H-%M)
LogDir=${ProjDir}/Log/${DateStart}
mkdir -p $LogDir

# Copy the parameter file to log
cp ${SCRIPT_FLAGS} ${LogDir}/Parameters_${DateStart}.txt



############
# Get ORFs #
############

# Check whether to get the ORFs from genomic sequence
if [ $GetOrf == 1 ]; then

    ${PBS_O_HOME}/bin/GetOrf.sh \
		${ProjDir}/Nuc_translation \
		0 \
		${TABLE} \
		${MINSIZE} \
		${CIRCULAR} \
		${GENOME} > ${LogDir}/GetOrf.log 2>&1
    ${PBS_O_HOME}/bin/GetOrf.sh \
		${ProjDir}/Nuc_translation \
		2 \
		${TABLE} \
		${MINSIZE} \
		${CIRCULAR} \
		${GENOME} >> ${LogDir}/GetOrf.log 2>&1

fi
export SIXFRAMEPROT=`basename $GENOME | perl -p -e "s%^%${ProjDir}\\/Nuc_translation\\/Find0_%" | perl -p -e "s/\\.fasta/_FIXED.fasta/"`
export SIXFRAMEGENE=`basename $GENOME | perl -p -e "s%^%${ProjDir}\\/Nuc_translation\\/Find2_%" | perl -p -e "s/\\.fasta/_FIXED.fasta/"`



#######################
# MaxQuant processing #
#######################

# Check whether to perform the maxquant processing
if [ $MaxquantProc == 1 ]; then
	
	singularity run $mq_version $mqpar > ${LogDir}/MaxquantProc.log 2>&1
	
fi



##################
# Get novel ORFs #
##################

# Check whether to get the novel ORFs from MaxQuant search results
if [ $GetNovelEntries == 1 ]; then

    ${PBS_O_HOME}/bin/GetNovelEntries.R \
		-o ${ProjDir}/Novel_res \
		-m ${ProjDir}/MaxQuant/combined/txt \
		-r ${UNIREFPROT} \
		-n ${SIXFRAMEPROT} \
		-t ${THREADS} > ${LogDir}/GetNovelEntries.log 2>&1

fi



################################
# Download NCBI Blast database #
################################

# Check whether to get the NCBI blast databases
if [ $GetNCBIBlastDb == 1 ]; then
    
    ${PBS_O_HOME}/bin/BlastUpdate.sh \
		-o ${PBS_O_INITDIR}/BlastDB \
		-t ${THREADS} \
		"${BlastDbs}" > ${LogDir}/GetNCBIBlastDb.log 2>&1
    
fi



##############################
# Make custom Blast database #
##############################

# Check whether to create a blast database for proteome data
if [ $MakeBlastDbProt == 1 ]; then
    
    ${PBS_O_HOME}/bin/MakeBlastDb.sh \
		-i ${InputType} \
		-y "prot" \
		-t ${TaxId} \
		${OTHERPROT[@]} > ${LogDir}/MakeBlastDb.log 2>&1
    
fi



# Check whether to create a blast database for genome data
if [ $MakeBlastDbNuc == 1 ]; then
    
	${PBS_O_HOME}/bin/MakeBlastDb.sh \
		-i ${InputType} \
		-y "nucl" \
		-t ${TaxId} \
		${UNIREFGENE} \
		${SIXFRAMEGENE} \
		${GENOME} ${OTHERGENE} >> ${LogDir}/MakeBlastDb.log 2>&1

fi



###################################
# Blast entries against databases #
###################################

# Check whether to blast the ORFs against taxon reference protein and RNA and against all uniprot and ncbi
if [ $DbBlasting == 1 ]; then
	
	${PBS_O_HOME}/bin/SerialDbBlasting.sh \
		-o ${ProjDir}/Blast \
		-s ${BlastSoftware} \
		-x ${SerialBlastMap} \
		-t ${THREADS} > ${LogDir}/DbBlasting.log 2>&1
    
fi



##################################
# Best blast hits identification #
##################################

# Check whether to identify the best blast hits from all results of previous step
if [ $BestBlast == 1 ]; then

    ${PBS_O_HOME}/bin/BestBlasts.sh \
		-o ${ProjDir}/Blast ${ProjDir}/Blast/ORFnew_vs_*_annot > ${LogDir}/BestBlast.log 2>&1
    
fi



#########################
# Reciprocal best blast #
#########################

# Check whether to perform the reciprocal best blast against previously used databases
if [ $ReciprocalBlast == 1 ]; then
	
	${PBS_O_HOME}/bin/SerialDbBlasting.sh \
		-o ${ProjDir}/ReciprocalBlast \
		-s ${BlastSoftware} \
		-x ${SerialRecipBlastMap} \
		-t ${THREADS} > ${LogDir}/ReciprocalBlast.log 2>&1
    
fi



######################################
# Compile reciprocal best blast data #
######################################

# Check whether to perform reciprocal best blast hit data processing
if [ $ReciprocalBestBlast == 1 ]; then

    ${PBS_O_HOME}/bin/ReciprocalBestBlasts.sh \
		-o ${ProjDir}/ReciprocalBlast \
		-r ${ProjDir}/ReciprocalBlast \
		${ProjDir}/Blast/ORFnew_vs_*_annot > ${LogDir}/ReciprocalBestBlast.log 2>&1
    
fi



##################################
# Protein coordinate computation #
##################################

# Check whether to perform protein coordinate identification
if [ $ProteinCoordinate == 1 ]; then

    ${PBS_O_HOME}/bin/Genomic_position_from_blast.R \
		-f ${UNIREFPROT} \
		-b ${ProjDir}/Blast/Refprot_vs_Genome \
		-g ${GENOME} \
		-o ${ProjDir}/ProtPosition/Ref_prot_coordinates.txt > ${LogDir}/ProteinCoordinate.log 2>&1
	${PBS_O_HOME}/bin/ORF_coordinates.R \
		-f ${SIXFRAMEPROT} \
		-o ${ProjDir}/ProtPosition/Orf_prot_coordinates.txt >> ${LogDir}/ProteinCoordinate.log 2>&1
	
fi



######################
# Protein annotation #
######################

# Check whether to perform protein annotation
if [ $ProteinAnnotation == 1 ]; then

    ${PBS_O_HOME}/bin/Fasta_annotation.R \
		-f ${UNIREFPROT} \
		-t ${TaxId} \
		-c ${AnnotColumn} \
		-k ${AnnotKey} \
		-s ${AnnotSeparator} \
		-o ${ProjDir}/ProtAnnotation/Ref_prot_annotations.txt > ${LogDir}/ProteinAnnotation.log 2>&1
	${PBS_O_HOME}/bin/Transfer_annotation.R \
		-c ${ProjDir}/ReciprocalBlast/Best_Reciproc_Blast_cross-map_ORFprot_vs_Refprot_recip_annot \
		-a ${ProjDir}/ProtAnnotation/Ref_prot_annotations.txt \
		-k ${AnnotKey} \
		-r "sseqid" \
		-o ${ProjDir}/ProtAnnotation/Orf_prot_annotations.txt >> ${LogDir}/ProteinAnnotation.log 2>&1
	
fi



#######################
# BSgenome generation #
#######################

# Check whether to perform BSgenome package installation
if [ $BSgenomeForge == 1 ]; then

    ${PBS_O_HOME}/bin/BSgenomeForge.sh \
		-o ${ProjDir}/BSgenome \
		-s ${SEED} \
		-c ${Circular} \
		-x ${GENOME} > ${LogDir}/BSgenomeForge.log 2>&1
	
fi



##########################
# GRange object creation #
##########################

# Check whether to perform GRange object creation
if [ $GRangeCreate == 1 ]; then
	
	${PBS_O_HOME}/bin/GRanges_generation.R \
		-g ${GENOME} \
		-b ${PKGNAME} \
		-o ${ProjDir}/GRanges > ${LogDir}/GRangesGeneration.log 2>&1
    ${PBS_O_HOME}/bin/GRanges_generation.R \
		-c ${ProjDir}/ProtPosition/Ref_prot_coordinates.txt \
		-g ${GENOME} \
		-b ${PKGNAME} \
		-a ${ProjDir}/ProtAnnotation/Ref_prot_annotations.txt \
		-o ${ProjDir}/GRanges >> ${LogDir}/GRangesGeneration.log 2>&1
	${PBS_O_HOME}/bin/GRanges_generation.R \
		-c ${ProjDir}/ProtPosition/Orf_prot_coordinates.txt \
		-g ${GENOME} \
		-b ${PKGNAME} \
		-a ${ProjDir}/ProtAnnotation/Orf_prot_annotations.txt \
		-o ${ProjDir}/GRanges >> ${LogDir}/GRangesGeneration.log 2>&1
	
fi



######################################
# Peptide novelty reason explanation #
######################################

# Check whether to perform the peptide novelty reason explanation
if [ $PepNoveltyReason == 1 ]; then
	
	${PBS_O_HOME}/bin/Novelty_discovery_peptides.R \
		-e ${ProjDir}/Novel_res/Group_evidence.RDS \
		-f ${UNIREFPROT} \
		-r ${ProjDir}/ReciprocalBlast/Best_Reciproc_Blast_ORFprot_vs_Refprot_recip_annot \
		-u ${ProjDir}/ReciprocalBlast \
		-p ${ProjDir}/Novel_res/Peptides_location.RDS \
		-t ${THREADS} \
		-o ${ProjDir}/NoveltyExplain > ${LogDir}/PepNoveltyReason.log 2>&1

fi



##################################
# Peptide coordinate computation #
##################################

# Check whether to perform peptide coordinate identification
if [ $SequencesCoordinate == 1 ]; then

	${PBS_O_HOME}/bin/SequenceGrange.sh \
		-o ${ProjDir}/GRanges \
		-c ${ProjDir}/PeptPosition \
		-n ${ProjDir}/NoveltyExplain/Sequence_novelty_reason.RDS \
		-r ${ProjDir}/ProtPosition/Ref_prot_coordinates.txt \
		-f ${ProjDir}/ProtPosition/Orf_prot_coordinates.txt \
		-g ${GENOME} \
		-b ${PKGNAME} \
		${ProjDir}/Novel_res >> ${LogDir}/SequencesCoordinate.log 2>&1
	
fi



#################################
# Operon coordinate computation #
#################################

# Check whether to perform operon coordinate identification
if [ $OperonCoordinate == 1 ]; then

    ${PBS_O_HOME}/bin/Operon_reformatting.R \
		-i ${OPERON} \
		-o ${ProjDir}/ProtPosition/Operon_coordinates.txt > ${LogDir}/OperonCoordinate.log 2>&1
	${PBS_O_HOME}/bin/GRanges_generation.R \
		-c ${ProjDir}/ProtPosition/Operon_coordinates.txt \
		-g ${GENOME} \
		-b ${PKGNAME} \
		-o ${ProjDir}/GRanges >> ${LogDir}/OperonCoordinate.log 2>&1
	
fi



################################
# Sanger validation coordinate #
################################

# Check whether to compute coordinate of Sanger sequence
if [ $SangerCoordinate == 1 ]; then
    
    ${PBS_O_HOME}/bin/MakeBlastDb.sh \
		-i ${InputType} \
		-y "nucl" \
		-t ${TaxId} \
		${SANGER} > ${LogDir}/SangerCoordinate.log 2>&1
    ${PBS_O_HOME}/bin/BlastDbBlasting.sh \
		-o ${ProjDir}/Sanger_validation \
		-y "nucl" \
		-l "all" \
		-a "blastn" \
		-s ${SANGER} \
		-q ${GENOME} \
		-e ${Eval} \
		-n ${NumAlign} \
		-t ${THREADS} \
		-b "Sanger_vs_Genome" >> ${LogDir}/SangerCoordinate.log 2>&1
	${PBS_O_HOME}/bin/Genomic_position_from_blast.R
		-f ${SANGER} \
		-b ${ProjDir}/Sanger_validation/Sanger_vs_Genome \
		-g ${GENOME} \
		-o ${ProjDir}/Sanger_validation/Sanger_seq_coordinates.txt >> ${LogDir}/SangerCoordinate.log 2>&1
	${PBS_O_HOME}/bin/GRanges_generation.R \
		-c ${ProjDir}/Sanger_validation/Sanger_seq_coordinates.txt \
		-g ${GENOME} \
		-n ${GenomeName} \
		-t ${Circular} \
		-o ${ProjDir}/GRanges >> ${LogDir}/SangerCoordinate.log 2>&1
	
fi



##################################
# ORF novelty reason explanation #
##################################

# Check whether to perform the ORF novelty reason explanation
if [ $OrfNoveltyReason == 1 ]; then
	
	${PBS_O_HOME}/bin/Novelty_discovery_ORFs.R \
		-i ${ProjDir}/NoveltyExplain/Sequence_novelty_reason.RDS \
		-r ${ProjDir}/GRanges/Ref_prot_grange.RDS \
		-n ${ProjDir}/GRanges/Orf_prot_grange.RDS \
		-p ${ProjDir}/GRanges/Operon_grange.RDS \
		-e ${ProjDir}/Novel_res/Peptides_location.RDS \
		-d ${ProjDir}/GRanges/Peptides_grange.RDS \
		-s ${ProjDir}/GRanges/Sanger_seq_grange.RDS \
		-g ${ProjDir}/GRanges/Genome_grange.RDS \
		-a "${AddRBS}" \
		-b ${PKGNAME} \
		-t ${THREADS} \
		-o ${ProjDir}/NoveltyExplain > ${LogDir}/OrfNoveltyReason.log 2>&1

fi



##############################
# Clustal multiple alignment #
##############################

# Check whether to perform multiple pair-wise alignments for RBS analysis
if [ $ClustalAlign == 1 ]; then

    ${PBS_O_HOME}/bin/ClustalAlignment.sh \
		-o ${ProjDir}/RBS_clustalo \
		-t ${THREADS} \
		${RBSLOW}.fasta ${RBSHIGH}.fasta > ${LogDir}/ClustalAlign.log 2>&1
    ${PBS_O_HOME}/bin/Seqlogo.sh \
		-o ${ProjDir}/RBS_clustalo \
		${RBSLOW}_alignment_full.clustal \
		${RBSHIGH}_alignment_full.clustal >> ${LogDir}/ClustalAlign.log 2>&1
    
fi



#########################
# Conservation analysis #
#########################

# Check whether to perform the conservation blast
if [ $ConservationBlast == 1 ] ; then
	
	${PBS_O_HOME}/bin/BlastDbBlasting.sh \
		-o ${ProjDir}/Conservation \
		-y "nucl" \
		-l ${ProjDir}/Novel_res/Novel_ORF.txt \
		-a "blastn" \
		-s ${SIXFRAMEGENE} \
		-q ${NCBInt} \
		-e ${Eval} \
		-n ${NumAlign} \
		-m "${BlastnParam}" \
		-t 1 \
		-b "ORFnucl_vs_NCBInt" >> ${LogDir}/ConservationBlast.log 2>&1
	
fi



##################################
# MaxQuant processing validation #
##################################

# Check whether to perform the maxquant processing to validate novel ORFs using public data
if [ $MaxquantValidation == 1 ]; then
	
	singularity run $mq_version $mqpar_valid > ${LogDir}/MaxquantValidation.log 2>&1
	
	find ${ProjDir}/ORF_validation/ -type f -name "proteinGroups.txt" -print0 | 
	while IFS= read -r -d '' file; do
		prefix=`dirname "$file" | xargs dirname | xargs dirname | xargs basename`
		outdir=`dirname "$file"`
		outname=`basename "$file" | sed "s/^/${prefix}_/"`
		mv "${file}" "${outdir}/${outname}"
	done
	
fi



###################################
# SummarizedExperiment generation #
###################################

# Check whether to generate the summarizedExperiment for all datasets
if [ $SummarizedExp == 1 ]; then
	
	SummarizedExpPreparation.sh \
		-o ${ProjDir}/SummarizedExp \
		${ProjDir}/ORF_validation ${ProjDir}/MaxQuant ${ProjDir}/GRanges > ${LogDir}/SummarizedExp.log 2>&1
	SummarizedExp_generation.R \
		-i ${ProjDir}/SummarizedExp \
		-o ${ProjDir}/SummarizedExp >> ${LogDir}/SummarizedExp.log 2>&1
	
fi



###########################
# InterProScan annotation #
###########################

# Check whether to get annotation via InterProScan for all fasta files
if [ $InterproScan == 1 ]; then
	
	InterproScan.sh \
		-o ${ProjDir}/InterPro \
		-t ${THREADS} \
		${UNIREFPROT} ${OTHERPROT[@]} ${ProjDir}/Novel_res/*.fasta > ${LogDir}/InterproScan.log 2>&1
	
fi



######################
# Signalp prediction #
######################

# Check whether to predict signal sequence via signalp for all fasta files
if [ $SignalpPrediction == 1 ]; then
	
	SignalpCheck.sh \
		-o ${ProjDir}/Signalp \
		-x ${Organism} \
		-t ${THREADS} \
		${UNIREFPROT} ${OTHERPROT[@]} ${ProjDir}/Novel_res/*.fasta > ${LogDir}/SignalpPrediction.log 2>&1
	
fi


