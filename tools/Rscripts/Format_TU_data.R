


my_tu_f <- "H:/data/Synechocystis_6frame/Kopf_2014_TU/Transcriptional_units_expr.txt"
my_tu <- data.table::fread(
    input = my_tu_f, sep = "\t", header = TRUE,
    stringsAsFactors = FALSE)

my_tu_final <- my_tu %>%
    dplyr::select(., `TU ID`, `TU type`, `Sense tags`, `Antisense TUs`) %>%
    tidyr::separate_rows(data = ., `Sense tags`, sep = ", ") %>%
    tidyr::separate_rows(data = ., `Antisense TUs`, sep = ", ") %>%
    dplyr::filter(., !is.na(`Sense tags`)) %>%
    dplyr::group_by(., `Sense tags`) %>%
    dplyr::summarise_all(~paste0(., collapse = ";")) %>%
    dplyr::ungroup(.)

data.table::fwrite(
    x = my_tu_final,
    file = "H:/data/Synechocystis_6frame/Kopf_2014_TU/Transcriptional_units_per_gene.txt",
    append = FALSE, quote = FALSE, sep = "\t",
    row.names = FALSE, col.names = TRUE)



