



library(magrittr)

my_data <- data.table::fread(
    input = "C:/Users/kxmna01/Desktop/lastz_test/bovis_vs_bcg_for_plotting",
    sep = "\t", quote = "", header = TRUE, stringsAsFactors = FALSE)

plot(my_data, type = "l")

c("#033A61", "#035416", "#6B0705", "#945725")

my_data_dup <- data.table::fread(
    input = "H:/data/Synechocystis_6frame/Conservation/Synechocystis_sp_PCC_6803_vs_itself_for_plotting",
    sep = "\t", quote = "", header = TRUE, stringsAsFactors = FALSE)

plot(my_data_dup, type = "l")


