


library(magrittr)

my_fasta_f <- "C:/Users/kxmna01/Desktop/CP023688_protein.fasta"

my_fasta <- seqinr::read.fasta(
    file = my_fasta_f, seqtype = "AA", as.string = TRUE, whole.header = TRUE)

my_data <- names(my_fasta) %>%
    data.frame(Header = ., stringsAsFactors = FALSE)

my_data_format <- my_data %>%
    tidyr::separate(
        data = ., col = "Header", into = c("ID", "REST"),
        sep = " ", extra = "merge") %>%
    dplyr::mutate(
        ., ID = sub(".+(_.+?)$", "CP023688\\1", ID),
        locus_tag = ifelse(
            grepl("locus_tag", REST),
            sub(".*\\[locus_tag=(.+?)\\].*", "\\1", REST),
            ""),
        gene_id = ifelse(
            grepl("GeneID", REST),
            sub(".*GeneID:(.+?)\\].*", "\\1", REST),
            ""),
        protein = ifelse(
            grepl("protein=", REST),
            sub(".*\\[protein=(.+?)\\].*", "\\1", REST),
            ""),
        protein_id = ifelse(
            grepl("protein_id", REST),
            sub(".*\\[protein_id=(.+?)\\].*", "\\1", REST),
            ""),
        pseudo = ifelse(
            grepl("pseudo=", REST),
            sub(".*\\[pseudo=(.+?)\\].*", "\\1", REST),
            "false"))
my_data_format$Header <- paste(
    my_data_format$ID, my_data_format$locus_tag,
    my_data_format$gene_id, my_data_format$protein,
    my_data_format$protein_id, my_data_format$pseudo,
    sep = " | ")

seqinr::write.fasta(
    sequences = my_fasta,
    names = my_data_format$Header,
    file.out = "CP023688_protein_FIXED.fasta",
    open = "w", nbchar = 60, as.string = TRUE)


