


# Peptide novelty
opt <- list(
    evidence = "C:\\Users\\kxmna01\\Dropbox\\Home_work_sync\\Work\\Srim_6frame\\Novel_res\\Group_evidence.RDS",
    reference_fasta = "C:\\Users\\kxmna01\\Dropbox\\Home_work_sync\\Work\\Srim_6frame\\Genome\\Streptomyces_rimosus_allStrains_2019-02-13.fasta",
    reciprocal_blast_ref = "C:\\Users\\kxmna01\\Dropbox\\Home_work_sync\\Work\\Srim_6frame\\ReciprocalBlast\\Best_Reciproc_Blast_Refprot_vs_ORFprot_annot",
    reciprocal_blast_uniprot = "C:\\Users\\kxmna01\\Dropbox\\Home_work_sync\\Work\\Srim_6frame\\ReciprocalBlast\\Best_Reciproc_Blast_Uniprot_vs_ORFprot_annot",
    reciprocal_blast_ncbi = "C:\\Users\\kxmna01\\Dropbox\\Home_work_sync\\Work\\Srim_6frame\\ReciprocalBlast\\Best_Reciproc_Blast_NCBIprot_vs_ORFprot_annot",
    peptide_location = "C:\\Users\\kxmna01\\Dropbox\\Home_work_sync\\Work\\Srim_6frame\\Novel_res\\Peptides_location.RDS",
    threads = 2,
    output = "NoveltyExplain"
)

# ORF novelty
opt <- list(
    novel_reason = "C:\\Users\\kxmna01\\Desktop\\NoveltyExplain\\Sequence_novelty_reason.RDS",
    ref_grange = "C:\\Users\\kxmna01\\Desktop/GRanges/Ref_prot_grange.RDS",
    orf_grange = "C:\\Users\\kxmna01\\Desktop/GRanges/Orf_prot_grange.RDS",
    operon_grange = character(0),
    pep_pos = "C:\\Users\\kxmna01\\Dropbox\\Home_work_sync\\Work\\Srim_6frame\\Novel_res\\Peptides_location.RDS",
    pep_grange = "C:\\Users\\kxmna01\\Dropbox\\Home_work_sync\\Work\\Srim_6frame\\GRanges\\Pept_seq_grange.RDS",
    sanger_grange = character(0),
    genome_grange = "C:\\Users\\kxmna01\\Dropbox\\Home_work_sync\\Work\\Srim_6frame\\GRanges\\Genome_grange.RDS",
    add_rbs = "",
    pep_class = c('class 1', 'class 2'),
    threads = 2,
    output = "NoveltyExplain"
)

# ref grange
opt <- list(
    coordinates = "C:\\Users\\kxmna01\\Dropbox\\Home_work_sync\\Work\\Srim_6frame\\ProtPosition\\Ref_prot_coordinates.txt",
    genome = "C:\\Users\\kxmna01\\Dropbox\\Home_work_sync\\Work\\Srim_6frame\\Genome\\Srim_G7_assembly.fasta",
    geno_name = "062019",
    circular = FALSE,
    annotations = "C:\\Users\\kxmna01\\Desktop\\ProtAnnotation\\Ref_prot_annotations.txt",
    output = "./GRanges/Ref_prot_grange.RDS"
)

# orf grange
opt <- list(
    coordinates = "C:\\Users\\kxmna01\\Dropbox\\Home_work_sync\\Work\\Srim_6frame\\ProtPosition\\Orf_prot_coordinates.txt",
    genome = "C:\\Users\\kxmna01\\Dropbox\\Home_work_sync\\Work\\Srim_6frame\\Genome\\Srim_G7_assembly.fasta",
    geno_name = "062019",
    circular = FALSE,
    annotations = "C:\\Users\\kxmna01\\Desktop\\ProtAnnotation\\Orf_prot_annotations.txt",
    output = "./GRanges/Orf_prot_grange.RDS"
)

# strathclyde grange
opt <- list(
    coordinates = "C:\\Users\\kxmna01\\Dropbox\\Home_work_sync\\Work\\Srim_6frame\\ProtPosition\\Strath_prot_coordinates.txt",
    genome = "C:\\Users\\kxmna01\\Dropbox\\Home_work_sync\\Work\\Srim_6frame\\Genome\\Srim_G7_assembly.fasta",
    geno_name = "062019",
    circular = FALSE,
    annotations = "",
    output = "C:\\Users\\kxmna01\\Dropbox\\Home_work_sync\\Work\\Srim_6frame/GRanges/Strath_prot_grange.RDS"
)


# Peptide novelty
opt <- list(
    evidence = "C:\\Users\\kxmna01\\Dropbox\\Home_work_sync\\Work\\Srim_6frame\\Novel_res\\Group_evidence.RDS",
    reference_fasta = "C:\\Users\\kxmna01\\Dropbox\\Home_work_sync\\Work\\Srim_6frame\\Genome\\Streptomyces_rimosus_allStrains_2019-02-13.fasta",
    reciprocal_blast_ref = "C:\\Users\\kxmna01\\Dropbox\\Home_work_sync\\Work\\Srim_6frame\\ReciprocalBlast\\Best_Reciproc_Blast_Refprot_vs_ORFprot_annot",
    reciprocal_blast_uniprot = "C:\\Users\\kxmna01\\Dropbox\\Home_work_sync\\Work\\Srim_6frame\\ReciprocalBlast\\Best_Reciproc_Blast_Uniprot_vs_ORFprot_annot",
    reciprocal_blast_ncbi = "C:\\Users\\kxmna01\\Dropbox\\Home_work_sync\\Work\\Srim_6frame\\ReciprocalBlast\\Best_Reciproc_Blast_NCBIprot_vs_ORFprot_annot",
    reciprocal_blast_strath = "C:\\Users\\kxmna01\\Dropbox\\Home_work_sync\\Work\\Srim_6frame\\ReciprocalBlast\\Best_Reciproc_Blast_Strathprot_vs_ORFprot_annot",
    peptide_location = "C:\\Users\\kxmna01\\Dropbox\\Home_work_sync\\Work\\Srim_6frame\\Novel_res\\Peptides_location.RDS",
    threads = 2,
    output = "NoveltyExplain"
)
# Import the reciprocal blast best hits between novel ORFs and
# all uniprot proteins
reciprocal_blast_strath <- opt$reciprocal_blast_strath %>%
    as.character(.) %>%
    read.table(
        file = ., header = TRUE, sep = "\t", quote = "",
        as.is = TRUE, comment.char = "") %>%
    dplyr::mutate(
        ., staxid_blast = as.integer(staxid_blast),
        staxid_reciproc = as.integer(staxid_reciproc))
# Compile all reciprocal best hits
reciprocal_blast_all <- dplyr::bind_rows(
    data.frame(DB = "Reference", reciprocal_blast_ref, stringsAsFactors = F),
    data.frame(DB = "UniProt", reciprocal_blast_uniprot, stringsAsFactors = F),
    data.frame(DB = "NCBI", reciprocal_blast_ncbi, stringsAsFactors = F),
    data.frame(DB = "Strathclyde", reciprocal_blast_strath, stringsAsFactors = F))

# ORF novelty
opt <- list(
    novel_reason = "C:\\Users\\kxmna01\\Dropbox\\Home_work_sync\\Work\\Srim_6frame\\NoveltyExplain\\Sequence_novelty_reason.RDS",
    ref_grange = "C:\\Users\\kxmna01\\Dropbox\\Home_work_sync\\Work\\Srim_6frame\\GRanges/Ref_prot_grange.RDS",
    strath_grange = "C:\\Users\\kxmna01\\Dropbox\\Home_work_sync\\Work\\Srim_6frame\\GRanges/Strath_prot_grange.RDS",
    orf_grange = "C:\\Users\\kxmna01\\Dropbox\\Home_work_sync\\Work\\Srim_6frame\\GRanges/Orf_prot_grange.RDS",
    operon_grange = character(0),
    pep_pos = "C:\\Users\\kxmna01\\Dropbox\\Home_work_sync\\Work\\Srim_6frame\\Novel_res\\Peptides_location.RDS",
    pep_grange = "C:\\Users\\kxmna01\\Dropbox\\Home_work_sync\\Work\\Srim_6frame\\GRanges\\Pept_seq_grange.RDS",
    sanger_grange = character(0),
    genome_grange = "C:\\Users\\kxmna01\\Dropbox\\Home_work_sync\\Work\\Srim_6frame\\GRanges\\Genome_grange.RDS",
    add_rbs = "",
    pep_class = c('class 1', 'class 2'),
    threads = 2,
    output = "NoveltyExplain"
)
# Import the reference genomic ranges file
strath_grange <- opt$strath_grange %>%
    as.character(.) %>%
    readRDS(file = .)


