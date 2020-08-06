

# Extract CDS amino acid sequences
from Bio import SeqIO
gbk_filename = "H:\\Cyanobacteria\\genbank\\pSYSM.genbank"
faa_filename = "H:\\Cyanobacteria\\genbank\\pSYSM_cds_aa.fasta"
input_handle  = open(gbk_filename, "r")
output_handle = open(faa_filename, "w")

for seq_record in SeqIO.parse(input_handle, "genbank") :
    print("Dealing with GenBank record %s" % seq_record.id)
    for seq_feature in seq_record.features :
        if seq_feature.type=="CDS" :
            if 'translation' in seq_feature.qualifiers :
                output_handle.write(">%s | %s | %s\n%s\n" % (
                    seq_feature.qualifiers['gene'][0],
                    #seq_feature.qualifiers['locus_tag'][0],
                    #seq_feature.qualifiers['label'][0],
                    seq_record.name,
                    #seq_feature.qualifiers['product'][0],
                    seq_feature.location,
                    seq_feature.qualifiers['translation'][0]))
            else :
                print("Translation qualifier not found %s" % seq_feature.qualifiers['locus_tag'][0])

output_handle.close()
input_handle.close()


# Extract genomic sequence
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_dna
gbk_filename = "H:\\Cyanobacteria\\genbank\\pSYSX.genbank"
faa_filename = "H:\\Cyanobacteria\\genbank\\pSYSX.fasta"
input_handle  = open(gbk_filename, "r")

for seq_record in SeqIO.parse(input_handle, "genbank") :
    print("Dealing with GenBank record %s" % seq_record.id)
    SeqIO.write(seq_record, faa_filename, "fasta")

input_handle.close()


# Extract gene nucleotide sequences
from Bio import SeqIO
gbk_filename = "H:\\data\\Synechocystis_6frame\\genbank\\pSYSX.genbank"
faa_filename = "H:\\data\\Synechocystis_6frame\\genbank\\pSYSX_gene_nucl.fasta"
input_handle  = open(gbk_filename, "r")
output_handle = open(faa_filename, "w")

for seq_record in SeqIO.parse(input_handle, "genbank") :
    print("Dealing with GenBank record %s" % seq_record.id)
    for seq_feature in seq_record.features :
        if seq_feature.type=="gene" :
            output_handle.write(">%s | %s | %s\n%s\n" % (
                seq_feature.qualifiers['gene'][0],
                #seq_feature.qualifiers['locus_tag'][0],
                seq_record.name,
                seq_feature.location,
                seq_feature.location.extract(seq_record).seq))

output_handle.close()
input_handle.close()


