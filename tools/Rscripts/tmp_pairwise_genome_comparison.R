



library(magrittr)

my_data <- data.table::fread(
    input = "H:/data/Synechocystis_6frame/Conservation/Synechocystis_sp_PCC_6803_vs_Synechococcus_elongatus_PCC_7942_for_plotting",
    sep = "\t", quote = "", header = TRUE, stringsAsFactors = FALSE)

plot(my_data, type = "l")

c("#033A61", "#035416", "#6B0705", "#945725")

my_data_dup <- data.table::fread(
    input = "H:/data/Synechocystis_6frame/Conservation/Synechocystis_sp_PCC_6803_vs_itself_for_plotting",
    sep = "\t", quote = "", header = TRUE, stringsAsFactors = FALSE)

plot(my_data_dup, type = "l")




my_data_elongatus <- data.table::fread(
    input = "H:/data/Synechocystis_6frame/Conservation/Synechococcus_elongatus_PCC_6301_vs_Synechococcus_elongatus_PCC_7942_for_plotting",
    sep = "\t", quote = "", header = TRUE, stringsAsFactors = FALSE)

plot(my_data_elongatus, type = "l")


