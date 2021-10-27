


library(magrittr)

my_fasta_f <- "H:/data/Srim_6frame_3rd_analysis/Genome/CP048261_CP048262_prot_sequence_FIXED.fasta"

my_fasta <- seqinr::read.fasta(
    file = my_fasta_f, seqtype = "AA", as.string = TRUE, whole.header = TRUE)

# Make sure there are no NA
if (any(is.na(my_fasta))) {
    stop("Cannot have NA in fasta sequences!")
}

# Check for presence of stop sign within sequence (sequence will be removed)
stop_filter <- grepl(".\\*.", my_fasta)
if (any(stop_filter)) {
    
    warning(paste0(
        sum(stop_filter),
        " entries will be deleted due to '*' in the middle of sequence!"))
    my_fasta <- my_fasta[!stop_filter]
    
}

my_data <- names(my_fasta) %>%
    data.frame(Header = ., stringsAsFactors = FALSE)

my_data_format <- my_data %>%
    tidyr::separate(
        data = ., col = "Header", into = c("ID", "REST"),
        sep = " ", extra = "merge") %>%
    dplyr::mutate(
        ., ID = sub(".+(_.+?)$", "AL939104\\1", ID),
        locus_tag = ifelse(
            grepl("locus_tag", REST),
            sub(".*\\[locus_tag=(.+?)\\].*", "\\1", REST),
            ""),
        #gene_id = ifelse(
        #    grepl("GeneID", REST),
        #    sub(".*GeneID:(.+?)\\].*", "\\1", REST),
        #    ""),
        gene_id = ifelse(
            grepl("gene", REST),
            sub(".*\\[gene=(.+?)\\].*", "\\1", REST),
            ""),
        protein = ifelse(
            grepl("protein=", REST),
            sub(".*\\[protein=(.+?)\\].*", "\\1", REST),
            ""),
        protein_id = ifelse(
            grepl("protein_id", REST),
            sub(".*\\[protein_id=(.+?)\\].*", "\\1", REST),
            ""),
        #protein_id = ifelse(
        #    grepl("CP.+_prot", REST),
        #    sub(".*\\|CP.+_(prot.+?) .*", "\\1", REST),
        #    ""),
        pseudo = ifelse(
            grepl("pseudo=", REST),
            sub(".*\\[pseudo=(.+?)\\].*", "\\1", REST),
            "false"))
my_data_format$Header <- paste(
    my_data_format$ID, my_data_format$locus_tag,
    my_data_format$gene_id, my_data_format$protein,
    my_data_format$protein_id, my_data_format$pseudo,
    sep = " | ")

# Check that IDs are no more than 15 characters long or provide warning
if (any(nchar(my_data_format$ID) > 14)) {
    stop("Some ID are too long (> 14 char.), it's a problem for Diamond!")
}

seqinr::write.fasta(
    sequences = my_fasta,
    names = my_data_format$Header,
    file.out = sub("\\.fasta", "_FIXED.fasta", my_fasta_f),
    open = "w", nbchar = 60, as.string = TRUE)


