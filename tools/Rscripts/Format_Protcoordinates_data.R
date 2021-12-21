


library(magrittr)

my_files <- c(
    "H:/data/Synechocystis_6frame/ProtPosition/Ref_prot_coordinates.txt",
    "H:/data/Synechocystis_6frame/ProtPosition/Orf_prot_coordinates.txt")

my_loc <- lapply(my_files, function(x) {
    data.table::fread(
        input = x, sep = "\t", quote = "", header = TRUE,
        stringsAsFactors = FALSE)
}) %>%
    plyr::ldply(., data.table::data.table, .id = NULL)

data.table::fwrite(
    x = my_loc,
    file = "H:/data/Synechocystis_6frame/ProtPosition/Ref_and_ORF_prot_coordinates.txt",
    append = FALSE, quote = FALSE, sep = "\t",
    row.names = FALSE, col.names = TRUE)


