


library(seqinr)

my_fasta_file <- "C:/Users/kxmna01/Dropbox/Home_work_sync/Work/Srim_6frame/G7 Sequence Strathclyde/Strathclyde_protein.fasta"

my_fasta <- seqinr::read.fasta(
    file = my_fasta_file, seqtype = "AA",
    as.string = TRUE, whole.header = TRUE)

my_fasta_info <- names(my_fasta) %>%
    data.frame(header = .) %>%
    cSplit(
        indt = ., splitCols = "header", direction = "wide",
        sep = " | ", fixed = TRUE, drop = FALSE) %>%
    dplyr::mutate(., header_5 = gsub("\\[|\\)", "", header_5)) %>%
    cSplit(
        indt = ., splitCols = "header_5", direction = "wide",
        sep = "(:|\\]\\()", fixed = FALSE, drop = FALSE) %>%
    set_colnames(c(
        "header", "id", "gene_name", "source", "protein_name", "location",
        "start", "end", "strand"))

my_fasta_format <- my_fasta_info %>%
    dplyr::filter(., strand == "+") %>%
    dplyr::select(
        ., id, start, end, strand,
        `GENE-NAME` = gene_name, `PROTEIN-NAME` = protein_name)

my_fasta_format <- my_fasta_info %>%
    dplyr::filter(., strand == "-") %>%
    dplyr::select(
        ., id, start = end, end = start, strand,
        `GENE-NAME` = gene_name, `PROTEIN-NAME` = protein_name) %>%
    dplyr::bind_rows(my_fasta_format, .) %>%
    dplyr::mutate(., chromosome = ifelse(grepl("^p", id), "plasmid1", "chr1"))

write.table(
    x = my_fasta_format, file = "Strath_prot_coordinates.txt",
    append = FALSE, quote = FALSE, sep = "\t",
    row.names = FALSE, col.names = TRUE)


