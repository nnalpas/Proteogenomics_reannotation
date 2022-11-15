


rm(list = ls())

library(magrittr)
library(ggplot2)

my_big_cols <- c(
    "#387eb8", "#404040", "#e21e25", "#fbaf3f",
    "#d1d2d4", "#246E39", "#753B94", "#067A91")

set_include <- c(
    "Pel1J_accepted", "Pel1J_unique_accepted",
    "Pel4J_accepted", "Pel4J_unique_accepted",
    "Plk_accepted", "Plk_unique_accepted",
    "Pel_common_accepted", "Pel_Plk_common_accepted")

mandatory_path <- c()

mandatory_resource <- c()

my_plots <- list()

#my_oa_f <- "H:/data/Synechocystis_6frame/Phylostratigraphy/Phylostrata_proteins.txt"

#my_oa_f <- "H:/data/Synechocystis_6frame/2022-03-04_ORF_validation/2022-06-09_Phylostrata_for_OA.txt.OA.txt"

my_oa_f <- "C:/Users/nalpanic/SynologyDrive/Work/Abaumannii_trimeth/Analysis/Functional_enrichment/2022-11-15_Trimethylation (K)_multi_software_wide_for_OA.txt.OA.txt"

my_oa <- data.table::fread(
    input = my_oa_f, sep = "\t", quote = "",
    header = TRUE, stringsAsFactors = FALSE)

my_oa %<>%
    dplyr::rowwise(.) %>%
    dplyr::mutate(., overlap_perc = overlap / size * 100)

my_oa_filt <- my_oa %>%
    dplyr::filter(., set %in% set_include) %>%
    dplyr::filter(., !is.na(padj) & padj <= 0.2) %>%
    dplyr::group_by(., set) %>%
    dplyr::arrange(., dplyr::desc(score)) %>%
    dplyr::mutate(
        ., Rank = 1:dplyr::n(),
        Selection = dplyr::case_when(
            Rank %in% c(1:100) ~ "Top",
            TRUE ~ "")) %>%
    dplyr::ungroup(.) %>%
    dplyr::filter(., pathway %in% mandatory_path | Selection == "Top")

path_filt <- unique(my_oa_filt$pathway)

my_oa_final <- my_oa_filt %>%
    dplyr::filter(., pathway %in% path_filt & padj <= 0.5) %>%
    #dplyr::filter(., pathway %in% mandatory_path & padj <= 0.05) %>%
    dplyr::mutate(., padj_range = dplyr::case_when(
        padj <= 0.001 ~ "<= 0.001",
        padj <= 0.01 ~ "<= 0.01",
        padj <= 0.05 ~ "<= 0.05",
        padj <= 0.2 ~ "<= 0.2",
        TRUE ~ NA_character_
    ))

if (length(mandatory_resource) > 0) {
    my_oa_final %<>%
        dplyr::filter(
            ., resource %in% mandatory_resource)
}

my_oa_final$padj_range <- factor(
    x = my_oa_final$padj_range,
    levels = c("<= 0.001", "<= 0.01", "<= 0.05", "<= 0.2"),
    ordered= TRUE)
my_oa_final$resource <- factor(
    x = my_oa_final$resource,
    levels = sort(unique(my_oa_final$resource)),
    ordered= TRUE)
my_oa_final$pathway <- factor(
    x = my_oa_final$pathway,
    levels = sort(unique(my_oa_final$pathway)),
    ordered = TRUE)
my_oa_final$set <- factor(
    x = my_oa_final$set,
    levels = set_include,
    ordered = TRUE)

my_plots[["Heatmap_oa_sig"]] <- ggplot(
    my_oa_final,
    aes(x = set, y = pathway, fill = resource, alpha = padj_range)) +
    geom_tile(colour = "white") +
    ggpubr::theme_pubr() +
    theme(
        legend.position = "right",
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
    facet_grid(rows = vars(resource), scales = "free_y", space = "free_y") +
    scale_alpha_manual(
        values = c(
            `<= 0.001` = 1, `<= 0.01` = 0.8,
            `<= 0.05` = 0.6, `<= 0.05` = 0.4)) +
    scale_fill_manual(values = my_big_cols)

my_oa_final <- my_oa %>%
    dplyr::filter(., pathway %in% path_filt) %>%
    #dplyr::filter(., pathway %in% mandatory_path) %>%
    dplyr::mutate(., padj_range = dplyr::case_when(
        padj <= 0.001 ~ "<= 0.001",
        padj <= 0.01 ~ "<= 0.01",
        padj <= 0.05 ~ "<= 0.05",
        padj <= 0.2 ~ "<= 0.2",
        TRUE ~ NA_character_
    )) %>%
    #dplyr::filter(
    #    ., resource %in% c(
    #        "Miscellaneous", "best_og_Subcategory",
    #        "Characterization", "Keywords", "GOCC Term")) %>%
    dplyr::mutate(., Perc_bin = dplyr::case_when(
        overlap_perc == 0 ~ "0",
        overlap_perc <= 20 ~ "20",
        overlap_perc <= 40 ~ "40",
        overlap_perc <= 60 ~ "60",
        overlap_perc <= 80 ~ "80",
        TRUE ~ "100"
    ))

my_oa_final$padj_range <- factor(
    x = my_oa_final$padj_range,
    levels = c("<= 0.001", "<= 0.01", "<= 0.05", "<= 0.2"),
    ordered= TRUE)
my_oa_final$resource <- factor(
    x = my_oa_final$resource,
    levels = sort(unique(my_oa_final$resource)),
    ordered= TRUE)
my_oa_final$pathway <- factor(
    x = my_oa_final$pathway,
    levels = sort(unique(my_oa_final$pathway)),
    ordered = TRUE)
my_oa_final$set <- factor(
    x = my_oa_final$set,
    levels = set_include,
    ordered = TRUE)

my_fills <- c("#FFFFFF", "#FFE6E6", "#FFCCCC", "#FF9999", "#FF8080", "#FF4D4D") %>%
    set_names(seq(from = 0, to = 100, by = 20))

my_plots[["Heatmap_oa_perc"]] <- ggplot(
    my_oa_final,
    aes(x = set, y = pathway, fill = Perc_bin, colour = padj_range)) +
    geom_tile(height = 0.8, width = 0.8, size = 0.6) +
    ggpubr::theme_pubr() +
    theme(
        legend.position = "right",
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
    facet_grid(rows = vars(resource), scales = "free_y", space = "free_y") +
    scale_colour_manual(
        values = c(
            `<= 0.001` = "#545454", `<= 0.01` = "#545454",
            `<= 0.05` = "#545454", `<= 0.2` = "#545454")) +
     scale_fill_manual(values = my_fills)

my_date <- Sys.Date()
pdf(
    file = paste0(my_date, "_function_OA.pdf"),
    width = 14, height = 18)
my_plots
dev.off()


