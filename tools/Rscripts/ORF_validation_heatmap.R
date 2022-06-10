


rm(list = ls())

library(ggplot2)
library(magrittr)

my_plots <- list()

my_data_f <- "H:/data/Synechocystis_6frame/2022-03-04_ORF_validation/Venn_ORF_validation_fmeasure.RDS"

feature_cols <- c(
    "Proteins", "Novel_peptide_count", "Novel_sum_MSMS_Count",
    "ORFNoveltyReason", "MainID", "Peptide 2+", "RBS valid",
    "TU valid", "Fmeasure valid", "PX valid", "Novel_min_PEP")

my_novel <- readRDS(my_data_f)

feature <- my_novel %>%
    dplyr::select(., tidyselect::all_of(feature_cols)) %>%
    dplyr::mutate(., Reason = dplyr::case_when(
        grepl("match elsewhere", ORFNoveltyReason) ~ "Conflicting annotation",
        TRUE ~ ORFNoveltyReason
    ) %>% sub(";.+", "", .)) %>%
    dplyr::arrange(., dplyr::desc(as.integer(Novel_sum_MSMS_Count)), Novel_min_PEP)

feature$Proteins <- factor(
    x = feature$Proteins, levels = feature$Proteins, ordered = TRUE)

quart_msms <- summary(as.integer(feature$Novel_sum_MSMS_Count))
feature_format_heat <- feature %>%
    dplyr::mutate(., MSMS_Count = dplyr::case_when(
        as.integer(Novel_sum_MSMS_Count) <= quart_msms[["1st Qu."]] ~ "25%",
        as.integer(Novel_sum_MSMS_Count) <= quart_msms[["Median"]] ~ "50%",
        as.integer(Novel_sum_MSMS_Count) <= quart_msms[["3rd Qu."]] ~ "75%",
        TRUE ~ "100%")) %>%
    dplyr::select(
        ., Proteins, Reason, MSMS_Count, `Peptide 2+`, `TU valid`,
        `Fmeasure valid`, `PX valid`) %>%
    tidyr::pivot_longer(data = ., cols = -Proteins)

feature_format_heat$name <- factor(
    x = feature_format_heat$name,
    levels = c("Reason", "MSMS_Count", "Fmeasure valid", "PX valid", "TU valid", "Peptide 2+"),
    ordered = TRUE)

feature_format_heat$Proteins <- factor(
    x = feature_format_heat$Proteins,
    levels = feature$Proteins,
    ordered = TRUE)

my_cols <- c(
    `Potential alternate start` = "#d1d2d4", `Potential alternate end` = "#246E39",
    `SAV` = "#fbaf3f", `Potentially novel` = "#e21e25",
    `Conflicting annotation` = "#387eb8",
    `TRUE` = "#404040", `FALSE` = "white",
    `25%` = "#3C5941", `50%` = "#B4BD93", `75%` = "#E3C28D", `100%` = "#C7522B")

my_plots[["heat_orf_properties"]] <- ggplot(
    feature_format_heat,
    aes(x = name, y = Proteins, fill = value)) +
    geom_tile(stat = "identity", colour = "white") +
    ggpubr::theme_pubr() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
    scale_fill_manual(values = my_cols)

pdf(file = "Novel_ORF_reason_heatmap.pdf", width = 10, height = 20)
my_plots
dev.off()


