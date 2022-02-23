


library(magrittr)
library(ggplot2)

my_netsurf_csv_f <- "H:/data/Synechocystis_6frame/2022-02-23_NetsurfP/Phosphorylated_proteins.csv"

my_netsurf_txt_f <- "H:/data/Synechocystis_6frame/2022-02-23_NetsurfP/Phosphorylated_proteins_netsurfp.txt"

my_netsurf_csv <- data.table::fread(
    input = my_netsurf_csv_f, sep = ",",
    header = TRUE, stringsAsFactors = FALSE)

my_netsurf_txt <- data.table::fread(
    input = my_netsurf_txt_f, sep = "\t",
    header = TRUE, stringsAsFactors = FALSE)

my_netsurf_txt %<>%
    dplyr::select(
        ., Proteins, AA, Position, Class, AlphaProb, BetaProb, CoilProb)

my_netsurf <- dplyr::left_join(
    x = my_netsurf_csv, y = my_netsurf_txt,
    by = c("id" = "Proteins", "seq" = "AA", "n" = "Position"))

data.table::fwrite(
    x = my_netsurf,
    file = "H:/data/Synechocystis_6frame/2022-02-23_NetsurfP/Phosphorylated_proteins_formatted.txt",
    append = FALSE, quote = FALSE, sep = "\t",
    row.names = FALSE, col.names = TRUE)


