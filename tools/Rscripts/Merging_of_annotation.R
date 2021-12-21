


library(magrittr)

my_files <- c(
    SCyCode = "H:/data/Synechocystis_6frame/Custom_annotation/2021-12-21_annotation_session.txt",
    UniProt = "H:/data/Synechocystis_6frame/UniprotAnnotation/2021-12-21_UniProt_annotation_2_perseus.txt",
    EggNOG = "H:/data/Synechocystis_6frame/EggnogMapper/2021-12-21_EggNOG_annotation_2_perseus.txt"
)

my_annot <- lapply(my_files, function(x) {
    data.table::fread(
        input = x, sep = "\t", quote = "", header = TRUE,
        stringsAsFactors = FALSE, colClasses = "character")
})

my_annot_final <- dplyr::full_join(
    x = my_annot$SCyCode, y = my_annot$UniProt,
    by = "#query_name", suffix = c(".SCyCode", ".UniProt")) %>%
    dplyr::full_join(
        x = ., y = my_annot$EggNOG,
        by = "#query_name", suffix = c("", ".EggNOG"))

data.table::fwrite(
    x = my_annot_final,
    file = "H:/data/Synechocystis_6frame/Custom_annotation/2021-12-21_Custom_Uniprot_Eggnog_annotations.txt",
    append = FALSE, quote = FALSE, sep = "\t",
    row.names = FALSE, col.names = TRUE)


