


library(magrittr)
library(ggplot2)

my_plots <- list()

#my_oa_f <- "H:/data/Synechocystis_6frame/Phylostratigraphy/Phylostrata_proteins.txt"

my_oa_f <- "H:/data/Synechocystis_6frame/OA_GSEA/2021-12-23_Phylostrata_for_OA.txt.OA.txt"

my_oa <- data.table::fread(
    input = my_oa_f, sep = "\t", quote = "",
    header = TRUE, stringsAsFactors = FALSE)

my_oa %<>%
    dplyr::rowwise(.) %>%
    dplyr::mutate(., overlap_perc = overlap / size * 100)

mandatory_path <- my_oa[
    my_oa$resource == "Custom_annotation",][["pathway"]] %>%
    unique(.) %>%
    grep("pxd|ScyCode", ., value = TRUE, invert = TRUE)
mandatory_path <- my_oa[
        my_oa$resource == "best_og_Subcategory",][["pathway"]] %>%
    unique(.) %>%
    c(mandatory_path, .)

my_oa_filt <- my_oa %>%
    dplyr::filter(., !grepl("pxd|ScyCode", pathway)) %>%
    #dplyr::filter(., !resource %in% c(
    #    "CAZy", "PFAMs", "BiGG_Reaction",
    #    "EC level 2 name", "EC level 3 name",
    #    "KEGG_Reaction_Name", "KEGG_rclass_Name",
    #    "KEGG_Module_Name", "KEGG_brite_Name", "best_og_Category")) %>%
    dplyr::filter(., !is.na(padj) & padj <= 0.1) %>%
    dplyr::group_by(., set) %>%
    dplyr::arrange(., dplyr::desc(score)) %>%
    dplyr::mutate(
        ., Rank = 1:dplyr::n(),
        Selection = dplyr::case_when(
            Rank %in% c(1:20) ~ "Top",
            TRUE ~ "")) %>%
    dplyr::ungroup(.) %>%
    dplyr::filter(., pathway %in% mandatory_path | Selection == "Top")
    #dplyr::filter(., Selection == "Top")

path_filt <- unique(my_oa_filt$pathway)

my_oa_final <- my_oa_filt %>%
    dplyr::filter(., pathway %in% path_filt & padj <= 0.5) %>%
    dplyr::mutate(., padj_range = dplyr::case_when(
        padj <= 0.001 ~ "<= 0.001",
        padj <= 0.01 ~ "<= 0.01",
        padj <= 0.05 ~ "<= 0.05",
        padj <= 0.1 ~ "<= 0.1",
        padj <= 0.5 ~ "<= 0.5",
        TRUE ~ NA_character_
    ))

my_oa_final$padj_range <- factor(
    x = my_oa_final$padj_range,
    levels = c("<= 0.001", "<= 0.01", "<= 0.05", "<= 0.1", "<= 0.5"),
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
    levels = c(
        "ps 1 - cellular organisms",
        "ps 2 - Bacteria",
        "ps 3 - Terrabacteria group",
        "ps 4 - Cyanobacteria/Melainabacteria group",
        "ps 5 - Cyanobacteria",
        "ps 6 - Synechococcales",
        "ps 8 - Synechocystis",
        "ps 9 - unclassified Synechocystis",
        "ps 10 - Synechocystis sp. PCC 6803"),
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
            `<= 0.05` = 0.5, `<= 0.1` = 0.3,
            `<= 0.5` = 0.1))

my_oa_final <- my_oa %>%
    dplyr::filter(., pathway %in% path_filt) %>%
    dplyr::mutate(., padj_range = dplyr::case_when(
        padj <= 0.001 ~ "<= 0.001",
        padj <= 0.01 ~ "<= 0.01",
        padj <= 0.05 ~ "<= 0.05",
        TRUE ~ "> 0.05"
    ))

my_oa_final$padj_range <- factor(
    x = my_oa_final$padj_range,
    levels = c("<= 0.001", "<= 0.01", "<= 0.05", "> 0.05"),
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
    levels = c(
        "ps 1 - cellular organisms",
        "ps 2 - Bacteria",
        "ps 3 - Terrabacteria group",
        "ps 4 - Cyanobacteria/Melainabacteria group",
        "ps 5 - Cyanobacteria",
        "ps 6 - Synechococcales",
        "ps 8 - Synechocystis",
        "ps 9 - unclassified Synechocystis",
        "ps 10 - Synechocystis sp. PCC 6803"),
    ordered = TRUE)

my_plots[["Heatmap_oa_perc"]] <- ggplot(
    my_oa_final,
    aes(x = set, y = pathway, fill = overlap_perc, colour = padj_range)) +
    geom_tile() +
    ggpubr::theme_pubr() +
    theme(
        legend.position = "right",
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
    facet_grid(rows = vars(resource), scales = "free_y", space = "free_y") +
    scale_colour_manual(
        values = c(
            `<= 0.001` = "black", `<= 0.01` = "#545454",
            `<= 0.05` = "#A8A8A8", `> 0.05` = "white")) +
    scale_fill_gradient(low = "white", high = "darkred")

my_date <- Sys.Date()
pdf(
    file = paste0(my_date, "_phylostratigraphy_OA.pdf"),
    width = 14, height = 18)
my_plots
dev.off()


