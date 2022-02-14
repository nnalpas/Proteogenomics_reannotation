


library(magrittr)
library(ggplot2)

my_plots <- list()
my_cols <- c("#387eb8", "#d1d2d4", "#e21e25", "#fbaf3f", "#404040")

my_data_f <- "H:/data/Synechocystis_6frame/2022-02-14_iBAQ/Scy004_resuscitation_Scy001_iBAQ_norm.txt"
my_annot_f <- "H:/data/Synechocystis_6frame/Custom_annotation/2021-12-21_Custom_Uniprot_Eggnog_annotations.txt"

my_data <- data.table::fread(
    input = my_data_f, sep = "\t", quote = "", header = TRUE,
    stringsAsFactors = FALSE, colClasses = "character", data.table = FALSE)

my_annot <- data.table::fread(
    input = my_annot_f, sep = "\t", quote = "", header = TRUE,
    stringsAsFactors = FALSE, colClasses = "character", data.table = FALSE)

my_data[, which(colnames(my_data) != "Proteins")] %<>% 
    lapply(., as.double)

my_data_format <- my_data %>%
    dplyr::select(., Proteins, tidyselect::starts_with("SCy001_L")) %>%
    tidyr::pivot_longer(
        data = ., cols = -Proteins,
        names_to = "Replicates", values_to = "iBAQ") %>%
    dplyr::group_by(., Replicates) %>%
    dplyr::mutate(
        ., Total = sum(iBAQ, na.rm = TRUE),
        iBAQ_perc = iBAQ*100/Total) %>%
    dplyr::ungroup(.)

my_ranks <- my_data_format %>%
    dplyr::group_by(., Proteins) %>%
    dplyr::summarise(
        ., iBAQ_perc_mean = mean(iBAQ_perc, na.rm = TRUE),
        iBAQ_perc_se = sd(iBAQ_perc, na.rm = TRUE)/sqrt(dplyr::n())) %>%
    dplyr::ungroup(.) %>%
    dplyr::filter(., !is.na(iBAQ_perc_mean)) %>%
    dplyr::arrange(., dplyr::desc(iBAQ_perc_mean)) %>%
    dplyr::mutate(
        ., Rank = 1:dplyr::n(),
        iBAQ_perc_cum = cumsum(iBAQ_perc_mean),
        iBAQ_perc_quantile = dplyr::case_when(
            iBAQ_perc_cum <= 2.5 ~ 2.5,
            iBAQ_perc_cum <= 25 ~ 25,
            iBAQ_perc_cum <= 50 ~ 50,
            iBAQ_perc_cum <= 75 ~ 75,
            iBAQ_perc_cum <= 97.5 ~ 97.5,
            TRUE ~ 100))

my_ranks_stats <- my_ranks %>%
    dplyr::group_by(., iBAQ_perc_quantile) %>%
    dplyr::summarise(
    ., iBAQ_perc_mean = max(iBAQ_perc_mean),
    Count = dplyr::n()) %>%
    dplyr::rowwise(.) %>%
    dplyr::mutate(
        ., Label = paste0(
            "Q", iBAQ_perc_quantile, ": ", Count,
            " proteins (max. = ", round(iBAQ_perc_mean, 2), ")")) %>%
    dplyr::ungroup(.)

my_plots[["iBAQ_perc_sort"]] <- ggplot(
    my_ranks, aes(x = Rank, y = iBAQ_perc_cum)) +
    geom_point(stat = "identity", shape = 1) +
    ggpubr::theme_pubr() +
    geom_hline(
        data = my_ranks_stats,
        mapping = aes(yintercept = iBAQ_perc_quantile), linetype = "dotted") +
    geom_text(
        data = my_ranks_stats,
        mapping = aes(
            x = 2000, y = iBAQ_perc_quantile, label = Label),
        stat = "identity", vjust = 1.5)

my_ranks_top <- my_ranks %>%
    dplyr::filter(., iBAQ_perc_quantile == 25) %>%
    dplyr::left_join(
        x = ., y = my_annot, by = c("Proteins" = "#query_name")) %>%
    dplyr::mutate(., fill = dplyr::case_when(
        grepl("Phycobilisome", Custom_classification) ~ "Phycobilisome",
        grepl("PSII", Custom_classification) ~ "PSII",
        grepl("PSI", Custom_classification) ~ "PSI",
        grepl("Calvin cycle", Custom_classification) ~ "Calvin cycle",
        TRUE ~ "Others"))

my_ranks_top$Proteins <- factor(
    x = my_ranks_top$Proteins, levels = my_ranks_top$Proteins, ordered = TRUE)

my_cols_fill <- my_cols %>%
    set_names(c("Phycobilisome", "PSI", "Calvin cycle", "PSII", "Others"))

my_plots[["iBAQ_perc_top25"]] <- ggplot(
    my_ranks_top,
    aes(x = Proteins, y = iBAQ_perc_mean, fill = fill)) +
    geom_bar(stat = "identity", position = "dodge") +
    ggpubr::theme_pubr() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
    scale_fill_manual(values = my_cols_fill)

pdf("Protein_iBAQ_quantiles.pdf", 10, 10)
my_plots
dev.off()


