


library(magrittr)
library(ggplot2)

my_plots <- list()

my_oa_f <- "H:/data/Synechocystis_6frame/OA_GSEA/Phylostrata_for_OA.txt.OA.txt"

my_oa <- data.table::fread(
    input = my_oa_f, sep = "\t", quote = "",
    header = TRUE, stringsAsFactors = FALSE)

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
    dplyr::filter(., Selection == "Top")

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

my_plots[["Heatmap_0.5"]] <- ggplot(my_oa_final, aes(x = set, y = pathway, fill = resource, alpha = padj_range)) +
    geom_tile(colour = "white") +
    ggpubr::theme_pubr() +
    theme(
        legend.position = "right",
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
    facet_grid(rows = vars(resource), scales = "free_y", space = "free_y") +
    scale_alpha_manual(
        values = c(`<= 0.001` = 1, `<= 0.01` = 0.8, `<= 0.05` = 0.5, `<= 0.1` = 0.3, `<= 0.5` = 0.1))

pdf(file = "Phylostratigraphy_OA.pdf", width = 9, height = 14)
my_plots
dev.off()


