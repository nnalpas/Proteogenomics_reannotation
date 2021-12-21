


library(magrittr)

my_ps_f <- "H:/data/Synechocystis_6frame/EggnogMapper/EggNOG_annotation_2_perseus_PS.txt"
my_ps <- data.table::fread(
    input = my_ps_f, sep = "\t", quote = "", header = TRUE,
    stringsAsFactors = FALSE)

my_ps_format <- my_ps %>%
    dplyr::select(
        ., `#query_name`,
        Custom_classification = `Classification according to annotation`) %>%
    tidyr::separate_rows(data = ., Custom_classification, sep = "/") %>%
    dplyr::mutate(
        ., Custom_classification = sub("^ +", "", Custom_classification) %>%
            sub(" +$", "", .)) %>%
    dplyr::group_by(., `#query_name`) %>%
    dplyr::summarise_all(~paste0(., collapse = ";"))

data.table::fwrite(
    x = my_ps_format,
    file = "H:/data/Synechocystis_6frame/EggnogMapper/EggNOG_annotation_2_perseus_PS_parsed.txt",
    append = FALSE, quote = FALSE, sep = "\t",
    row.names = FALSE, col.names = TRUE)



