#!/bin/bash


# Script that creates a blast database from Uniprot complete entries

# Go to workspace
WKDIR="~/work/BlastDB"
mkdir ${WKDIR}

# Download Uniprot entries (swissprot and trembl)
wget "ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_*.fasta.gz"

# To get a single database make sure to compile entries
gunzip -c uniprot_*.fasta.gz >> uniprot_sprot_trembl.fasta
 
# Create the blast database
MakeBlastDb.sh "fasta" "prot" 131567 "${WKDIR}/uniprot_sprot_trembl.fasta" >> ${WKDIR}/MakeBlastDb.log 2>&1


