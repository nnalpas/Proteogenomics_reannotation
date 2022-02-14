


rm(list = ls())

library(magrittr)
library(ggplot2)

my_plots <- list()

my_pg_f <- "H:/data/Synechocystis_6frame/2022-02-14_iBAQ/Scy004_resuscitation_Scy001_iBAQ_norm_t-test.txt"

my_data <- data.table::fread(
    input = my_pg_f, sep = "\t", quote = "", header = TRUE,
    stringsAsFactors = FALSE, colClasses = "character", na.strings = "NaN")

my_data_format <- my_data %>%
    dplyr::filter(., !grepl("chr1", Proteins)) %>%
    dplyr::select(., Proteins, tidyselect::starts_with("T-test Significant ")) %>%
    tidyr::pivot_longer(data = ., cols = -Proteins) %>%
    dplyr::filter(., value == "+") %>%
    dplyr::mutate(
        ., value = 1, name = sub("T-test Significant ", "", name)) %>%
    tidyr::pivot_wider(
        data = ., names_from = name, values_from = value)
my_data_format[is.na(my_data_format)] <- 0

my_sig_unchar <- my_data %>%
    dplyr::filter(
        ., !grepl("chr1", Proteins) &
            Characterization == "Uncharacterized") %>%
    dplyr::select(
        ., Proteins, Characterization,
        tidyselect::starts_with("T-test Difference "),
        tidyselect::starts_with("T-test Significant ")) %>%
    tidyr::pivot_longer(data = ., cols = c(-Proteins, -Characterization)) %>%
    dplyr::mutate(
        ., Parameter = sub("(Difference|Significant).+", "\\1", name),
        Condition = sub(".+(Difference|Significant) ", "", name)) %>%
    dplyr::select(., -name) %>%
    tidyr::pivot_wider(data = ., names_from = Parameter, values_from = value) %>%
    dplyr::filter(., `T-test Significant` == "+") %>%
    dplyr::mutate(
        ., Experiment = ifelse(
            grepl("Vegetative", Condition),
            "Resuscitation", "Carbon down-shift"),
        Vegetative = "0",
        HC = "0") %>%
    dplyr::select(., Proteins, Experiment, Condition, Vegetative, HC) %>%
    unique(.)

my_data_toplot <- my_data %>%
    dplyr::left_join(x = my_sig_unchar, y = .) %>%
    dplyr::select(
        ., Proteins, Experiment, Condition, Vegetative, HC,
        tidyselect::starts_with("T-test Difference ")) %>%
    tidyr::pivot_longer(data = ., cols = c(-Proteins, -Experiment, -Condition)) %>%
    dplyr::mutate(
        ., name = sub(".+(Difference|Significant) ", "", name),
        value = as.double(value)) %>%
    dplyr::mutate(
        ., Exp = ifelse(
            grepl("Vegetative", name),
            "Resuscitation", "Carbon down-shift")) %>%
    dplyr::filter(., Experiment == Exp) %>%
    dplyr::arrange(., Proteins, Experiment) %>%
    dplyr::mutate(., Significant = ifelse(Condition == name, "Yes", "No"))

my_data_toplot$name <- factor(
    x = my_data_toplot$name,
    levels = c(
        "HC", "LC3/HC", "LC24/HC", "Vegetative", "Chlorosis/Vegetative",
        "Stage1Early/Vegetative", "Stage1Mid/Vegetative",
        "Stage2/Vegetative"),
    ordered = TRUE)

my_cols <- c(Yes = "#387eb8", No = "black")

for (x in c("Resuscitation", "Carbon down-shift")) {
    my_plots[[x]] <- ggplot(
        my_data_toplot %>% dplyr::filter(., Experiment == x),
        aes(
            x = name, y = value, group = Proteins,
            shape = Proteins, colour = Significant)) +
        geom_point(size = 5) +
        geom_line(linetype = "dashed", colour = "black") +
        ggpubr::theme_pubr() +
        scale_colour_manual(values = my_cols)
}

pdf("Scy004_resuscitation_Scy001_iBAQ_norm_t-test_profile.pdf", 10, 10)
my_plots
dev.off()

write.table(
    x = my_data_format,
    file = sub(".txt", "_forOA.txt", my_pg_f),
    append = FALSE, quote = FALSE, sep = "\t",
    row.names = FALSE, col.names = TRUE)

save.image("Protein_iBAQ_t-test.RData")


