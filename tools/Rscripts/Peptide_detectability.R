
library(magrittr)
library(ggplot2)

my_proteotypic_df <- data.table::fread(input = "H:/Reference_protease_cleavages_proteotypic.txt")
outfile <- "Reference_protease_cleavages"


data.table::fwrite(
    x = my_proteotypic_df %>% dplyr::filter(., Filter == "Pass") %>% as.data.frame(.) %>% .["Sequence"],
    file = paste0(outfile, "_sequence.txt"),
    append = FALSE, quote = FALSE, sep = "\t",
    row.names = FALSE, col.names = FALSE)


my_predict_f <- c(
    ref = "H:/Reference_protease_cleavages_sequence_Predictions.txt",
    hidden = "H:/Hidden_protease_cleavages_sequence_Predictions.txt")

my_predict <- lapply(my_predict_f, function(x) {
    data.table::fread(
        input = x, sep = "\t", quote = "",
        header = TRUE, stringsAsFactors = FALSE)
}) %>%
    plyr::ldply(., data.frame, .id = "Database")

summary(my_predict %>% dplyr::filter(., Database == "ref") %>% .[["Prob"]])
summary(my_predict %>% dplyr::filter(., Database == "hidden") %>% .[["Prob"]])

toplot <- my_predict %>%
    dplyr::group_by(., Database, Detectability) %>%
    dplyr::summarise(., Count = dplyr::n()) %>%
    dplyr::group_by(., Database) %>%
    dplyr::mutate(
        ., Percentage = Count * 100 / sum(Count))

my_pep <- data.table::fread(
    input = "H:/data/Synechocystis_6frame/PeptPosition/Peptides_coordinates.txt",
    sep = "\t", quote = "",
    header = TRUE, stringsAsFactors = FALSE)

toplot <- my_predict %>%
    dplyr::filter(., Peptide %in% unique(my_pep$id) & Database == "ref") %>%
    dplyr::group_by(., Database, Detectability) %>%
    dplyr::summarise(., Count = dplyr::n()) %>%
    dplyr::group_by(., Database) %>%
    dplyr::mutate(
        ., Percentage = Count * 100 / sum(Count),
        Database = "MS-identified") %>%
    dplyr::bind_rows(toplot, .)

toplot$Database <- factor(
    x = toplot$Database, levels = c("ref", "MS-identified", "hidden"),
    ordered = TRUE)

pdf("Hidden_detectability.pdf", 10, 10)
ggplot(
    toplot, aes(
        x = Detectability, y = Percentage,
        fill = Database, colour = Database, label = round(x = Percentage, 1))) +
    geom_bar(stat = "identity", position = "dodge", alpha = 0.7) +
    geom_text(position = position_dodge(width = 0.9), vjust = -0.3) +
    ggpubr::theme_pubr() +
    scale_fill_manual(values = c("#387eb8", "#d1d2d4", "#e21e25")) +
    scale_colour_manual(values = c("#387eb8", "#d1d2d4", "#e21e25"))
dev.off()


