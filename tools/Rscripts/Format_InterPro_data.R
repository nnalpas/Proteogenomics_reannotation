


library(magrittr)

my_interproscan_f <- "H:/data/Srim_6frame_3rd_analysis/InterPro/CP048261_CP048262_prot_sequence_FIXED3.tsv"

my_interproscan <- data.table::fread(
    input = my_interproscan_f, sep = "\t", quote = "", header = TRUE,
    stringsAsFactors = FALSE, fill = TRUE)

my_interproscan_format <- my_interproscan %>%
    dplyr::select(
        ., `Protein accession`, Analysis, `Signature accession`,
        `Signature description`, `Start location`, `Stop location`, Score)

my_interproscan_format <- my_interproscan %>%
    dplyr::select(
        ., `Protein accession`, `Signature accession` = `InterPro accession`,
        `Signature description` = `InterPro description`,
        `Start location`, `Stop location`, Score) %>%
    dplyr::mutate(., Analysis = "InterPro") %>%
    dplyr::bind_rows(my_interproscan_format, .) %>%
    dplyr::filter(., `Signature accession` != "-") %>%
    dplyr::mutate_all(~sub("^-$", NA, .)) %>%
    dplyr::mutate(
        ., `Start location` = as.integer(`Start location`),
        `Stop location` = as.integer(`Stop location`),
        Score = as.numeric(Score))

my_interproscan_filter <- my_interproscan_format %>%
    dplyr::filter(., Score <= 0.0001) %>%
    dplyr::mutate_all(~as.character(.)) %>%
    tidyr::pivot_longer(
        data = ., cols = c(-`Protein accession`, -`Analysis`),
        names_to = "Param", values_to = "Value") %>%
    tidyr::unite(data = ., col = "Name", Analysis, Param, sep = " ") %>%
    dplyr::mutate(., Value = gsub (";", " |", Value)) %>%
    dplyr::group_by(., `Protein accession`, Name) %>%
    dplyr::summarise(., Value = paste0(Value, collapse = ";")) %>%
    tidyr::pivot_wider(
        data = ., names_from = "Name", values_from = "Value") %>%
    dplyr::ungroup(.)

data.table::fwrite(
    x = my_interproscan_filter,
    file = sub(".txt", "parsed.txt", my_interproscan_f),
    append = FALSE, quote = FALSE, sep = "\t",
    row.names = FALSE, col.names = TRUE)


