#!/bin/bash


# Script that creates a blast database from Uniprot complete entries

# Go to workspace
cd ${PBS_O_WORKDIR}

# Download Uniprot entries (swissprot and trembl)
#wget "ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_*.fasta.gz"

# To get a single database make sure to compile entries
#gunzip -c uniprot_*.fasta.gz >> uniprot_sprot_trembl.fasta
 
# Create the blast database
${PBS_O_HOME}/bin/MakeBlastDb.sh "fasta" "prot" 131567 "${PBS_O_INITDIR}/uniprot_sprot_trembl.fasta" >> ./MakeBlastDb_{PBS_JOBID}.log


