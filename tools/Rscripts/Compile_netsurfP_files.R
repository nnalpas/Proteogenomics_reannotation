


rm(list = ls())

library(magrittr)
library(ggplot2)

my_netsurf_csv_f <- "C:/Users/nalpanic/SynologyDrive/Work/Abaumannii_trimeth/Netsurfp/Trimethylated/Trimethylated_proteins.csv"

my_netsurf_csv <- readLines(con = my_netsurf_csv_f)

my_netsurf_csv %<>%
    sub(">(p?1?ABYAL[0-9]+)_.+_(,[A-Z],[0-9])", "\\1\\2", .)

writeLines(text = my_netsurf_csv, con = "C:/Users/nalpanic/SynologyDrive/Work/Abaumannii_trimeth/Netsurfp/Trimethylated/Trimethylated_proteins_cleaned.csv")


my_netsurf_csv_f <- "C:/Users/nalpanic/SynologyDrive/Work/Abaumannii_trimeth/Netsurfp/Trimethylated/Trimethylated_proteins_cleaned.csv"

my_netsurf_txt_f <- "C:/Users/nalpanic/SynologyDrive/Work/Abaumannii_trimeth/Netsurfp/Trimethylated/Trimethylated_proteins.netsurfp.txt"

my_netsurf_csv <- data.table::fread(
    input = my_netsurf_csv_f, sep = ",",
    header = TRUE, stringsAsFactors = FALSE)

my_netsurf_txt <- data.table::fread(
    input = my_netsurf_txt_f, sep = "\t",
    header = TRUE, stringsAsFactors = FALSE) %>%
    dplyr::mutate(., Proteins = sub("(p?1?ABYAL[0-9]+)_.+", "\\1", Proteins))

my_netsurf_txt %<>%
    dplyr::select(
        ., Proteins, AA, Position, Class, AlphaProb, BetaProb, CoilProb)

my_netsurf <- dplyr::left_join(
    x = my_netsurf_csv, y = my_netsurf_txt,
    by = c("id" = "Proteins", "seq" = "AA", "n" = "Position"))

data.table::fwrite(
    x = my_netsurf,
    file = sub(".txt", "_formatted.txt", my_netsurf_txt_f),
    append = FALSE, quote = FALSE, sep = "\t",
    row.names = FALSE, col.names = TRUE)


