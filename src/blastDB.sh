#!/bin/bash

# Go to workspace
cd ${PBS_O_WORKDIR}

# Load the required modules
module load blast+/2.6.0

cd ${PBS_O_INITDIR}/BlastDB

# Download NCBI reference protein
wget "ftp://ftp.ncbi.nlm.nih.gov/blast/db/refseq_protein.*.tar.gz"
for file in `ls refseq_protein.*.tar.gz`; do
	tar -xzvf  $file
done


# Download uniprot databases to use in blast
wget "ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_trembl.fasta.gz"
gunzip ./*.gz
cat *.fasta > uniprot_sprot_trembl.fasta

makeblastdb -in ./uniprot_sprot_trembl.fasta -input_type "fasta" -title "uniprot_sprot_trembl" -parse_seqids -dbtype "prot"

