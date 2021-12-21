


library(magrittr)

my_panther_f <- "D:/Local_databases/X-databases/Panther/PANTHER16.0_HMM_classifications"
panther_col <- c(
    "PANTHER ID",
    "Name",
    "Molecular function",
    "Biological process",
    "Cellular components",
    "Protein class",
    "Pathway"
)
my_panther <- data.table::fread(
    input = my_panther_f, sep = "\t", quote = "", header = FALSE,
    stringsAsFactors = FALSE, col.names = panther_col)

my_panther_format <- my_panther

data.table::fwrite(
    x = my_panther_format,
    file = "D:/Local_databases/X-databases/Panther/Panther_parsed.txt",
    append = FALSE, quote = FALSE, sep = "\t",
    row.names = FALSE, col.names = TRUE)


