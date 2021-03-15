



library(magrittr)

my_data <- data.table::fread(
    input = "C:/Users/kxmna01/Desktop/Conservation/IPPAS_B-1465_vs_ASM34078v1_for_plotting",
    sep = "\t", quote = "", header = TRUE, stringsAsFactors = FALSE)

plot(my_data, type = "l")



