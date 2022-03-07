


rm(list = ls())

library(ggplot2)
library(magrittr)

my_plots <- list()

my_novel_f <- "/mnt/storage/kxmna01/data/Synechocystis_6frame/2022-03-04_ORF_validation/Venn_ORF_validation.txt"
feature_cols <- c(
    "Proteins", "Novel_peptide_count", "Novel_sum_MSMS_Count",
    "ORFNoveltyReason", "MainID", "Peptide 2+", "RBS valid",
    "TU valid", "High quality", "PX valid")

my_novel <- data.table::fread(
    input = my_novel_f, sep = "\t", quote = "", header = TRUE,
    stringsAsFactors = FALSE, colClasses = "character")

feature <- my_novel %>%
    dplyr::select(., tidyselect::all_of(feature_cols)) %>%
    dplyr::mutate(., Reason = dplyr::case_when(
        ORFNoveltyReason %in% c("Potentially novel", "SAV", "Potential alternate start", "Potential alternate end") ~ ORFNoveltyReason,
        TRUE ~ "Conflicting annotation"
    )) %>%
    dplyr::mutate(., MainID = ifelse(MainID != "", MainID, Proteins)) %>%
    dplyr::arrange(., dplyr::desc(as.integer(Novel_sum_MSMS_Count)))

feature$MainID <- factor(
    x = feature$MainID, levels = feature$MainID, ordered = TRUE)

quart_msms <- summary(as.integer(feature$Novel_sum_MSMS_Count))
feature_format_heat <- feature %>%
    dplyr::mutate(., MSMS_Count = dplyr::case_when(
        as.integer(Novel_sum_MSMS_Count) <= quart_msms[["1st Qu."]] ~ "25%",
        as.integer(Novel_sum_MSMS_Count) <= quart_msms[["Median"]] ~ "50%",
        as.integer(Novel_sum_MSMS_Count) <= quart_msms[["3rd Qu."]] ~ "75%",
        TRUE ~ "100%")) %>%
    dplyr::select(
        ., MainID, Reason, MSMS_Count, `Peptide 2+`, `TU valid`,
        `High quality`, `PX valid`) %>%
    tidyr::pivot_longer(data = ., cols = -MainID)

feature_format_heat$name <- factor(
    x = feature_format_heat$name,
    levels = c("Reason", "MSMS_Count", "High quality", "PX valid", "TU valid", "Peptide 2+"),
    ordered = TRUE)

feature_format_heat$MainID <- factor(
    x = feature_format_heat$MainID,
    levels = feature$MainID,
    ordered = TRUE)

my_cols <- c(
    `Potential alternate start` = "#387eb8", `Potential alternate end` = "#fbaf3f",
    `SAV` = "#246E39", `Potentially novel` = "#e21e25",
    `Conflicting annotation` = "#d1d2d4",
    `TRUE` = "#404040", `FALSE` = "white",
    `25%` = "#3C5941", `50%` = "#B4BD93", `75%` = "#E3C28D", `100%` = "#C7522B")

my_plots[["heat_orf_properties"]] <- ggplot(
    feature_format_heat,
    aes(x = name, y = MainID, fill = value)) +
    geom_tile(stat = "identity", colour = "white") +
    ggpubr::theme_pubr() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
    scale_fill_manual(values = my_cols)

pdf(file = "Novel_ORF_reason_heatmap.pdf", width = 10, height = 20)
my_plots
dev.off()


