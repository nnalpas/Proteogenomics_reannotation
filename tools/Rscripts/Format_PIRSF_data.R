


library(magrittr)

my_pirsf_f <- "D:/Local_databases/X-databases/PIRSF/pirsfinfo_format.dat"
pirsf_col <- c(
    "PIRSF_number",
    "curation_status",
    "name",
    "parent"
)
my_pirsf <- data.table::fread(
    input = my_pirsf_f, sep = "\t", quote = "", header = FALSE,
    stringsAsFactors = FALSE, col.names = "Header")

my_pirsf_format <- my_pirsf %>%
    tidyr::separate(
        data = ., col = "Header", into = c("PIRSF_number", "Rest"),
        sep = " ", extra = "merge") %>%
    tidyr::separate(
        data = ., col = "Rest", into = c("curation_status", "Header"),
        sep = "\\)", extra = "merge") %>%
    tidyr::separate(
        data = ., col = "Header", into = c("name", "parent"),
        sep = "\\[Parent=", extra = "warn") %>%
    dplyr::mutate_all(~sub("^(>| |\\()", "", .) %>% sub("( |\\])$", "", .))

data.table::fwrite(
    x = my_pirsf_format,
    file = "D:/Local_databases/X-databases/PIRSF/pirsfinfo_parsed.txt",
    append = FALSE, quote = FALSE, sep = "\t",
    row.names = FALSE, col.names = TRUE)


