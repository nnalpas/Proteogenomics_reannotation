


rm(list = ls())

library(magrittr)

my_path <- "/mnt/storage/kxmna01/data/Synechocystis_6frame/Codon_usage_other_species"

my_fasta_f <- list.files(
    path = my_path, pattern = "*_rna.fna$",
    full.names = TRUE, recursive = TRUE) %>%
    set_names(sub("_rna.fna", "", basename(.)))
my_prot_f <- list.files(
    path = my_path, pattern = "*_top_prot.txt",
    full.names = TRUE, recursive = TRUE) %>%
    set_names(sub("_top_prot.txt", "", basename(.)))

my_res <- lapply(names(my_fasta_f), function(x) {
#for (x in names(my_fasta_f)) {
    
    my_nucl_char <- seqinr::read.fasta(
        file = my_fasta_f[[x]], seqtype = "DNA",
        as.string = FALSE, whole.header = TRUE)
    
    my_high_exp <- data.table::fread(
        input = my_prot_f[[x]], sep = "\t", quote = "", header = FALSE,
        stringsAsFactors = FALSE, colClasses = "character", col.names = "ID")
    
    my_nucl_char_high <- my_nucl_char[
        grepl(paste0(
            "(", paste0(my_high_exp$ID, collapse = "|"),
            ")(\\.|,|;| |\\]|$)"), names(my_nucl_char))]
    my_w <- lapply(my_nucl_char_high, function(y) {
        seqinr::uco(
            seq = y, frame = 0, index = "eff")
    }) %>%
        plyr::ldply(., data.table::data.table) %>%
        set_colnames(c("qseqid", "CODON", "Count")) %>%
        dplyr::full_join(x = ., y = seqinr::SEQINR.UTIL$CODON.AA) %>%
        dplyr::group_by(., AA, L, CODON) %>%
        dplyr::summarise(., Sum = sum(Count, na.rm = TRUE)) %>%
        dplyr::group_by(., AA, L) %>%
        dplyr::mutate(., w = Sum/max(Sum)) %>%
        dplyr::ungroup(.) %>%
        dplyr::arrange(., CODON) %>%
        dplyr::select(., -Sum)
    
    my_w
    
}) %>%
    set_names(names(my_fasta_f)) %>%
    plyr::ldply(., dplyr::bind_rows, .id = "Species") %>%
    tidyr::pivot_wider(data = ., names_from = "Species", values_from = "w")

data.table::fwrite(
    x = my_res, file = "Codon_preference_other_species.txt",
    append = FALSE, quote = FALSE, sep = "\t",
    row.names = FALSE, col.names = TRUE)


