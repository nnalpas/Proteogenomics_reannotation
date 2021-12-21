


library(magrittr)

my_prosite_f <- "D:/Local_databases/X-databases/Prosite/2021-12/prosite.dat"
my_prosite <- data.table::fread(
    input = my_prosite_f, sep = "\t", quote = "", header = FALSE,
    stringsAsFactors = FALSE, col.names = "Header")

my_prosite_fields <- c(
    ID = "ID",
    CC = "Comment",
    AC = "Accession",
    DT = "Dates",
    DE = "Description",
    PA = "Pattern",
    PR = "Rule",
    DO = "Documentation",
    NR = "Notes",
    DR = "Matches",
    `3D` = "3D",
    MA = "Matrix",
    PP = "Exception"
)

my_prosite_format <- my_prosite %>%
    tidyr::separate(
        data = ., col = "Header", into = c("Name", "Value"),
        sep = "   ", extra = "merge") %>%
    dplyr::filter(., Name %in% c("ID", "AC", "DT", "DE", "PA", "//"))

my_prosite_format %<>%
    dplyr::mutate(., ID = dplyr::case_when(
        Name == "ID" ~ Value,
        TRUE ~ NA_character_)) %>%
    tidyr::fill(data = ., ID, .direction = "down") %>%
    dplyr::filter(., !Name %in% c("ID", "//")) %>%
    dplyr::group_by(., ID, Name) %>%
    dplyr::summarise_all(~paste0(., collapse = ";")) %>%
    tidyr::pivot_wider(data = ., names_from = "Name", values_from = "Value") %>%
    dplyr::mutate(
        ., ID = gsub(";", " -", ID), AC = gsub(";", "", AC))

data.table::fwrite(
    x = my_prosite_format,
    file = "D:/Local_databases/X-databases/Prosite/2021-12/prosite_parsed.txt",
    append = FALSE, quote = FALSE, sep = "\t",
    row.names = FALSE, col.names = TRUE)


