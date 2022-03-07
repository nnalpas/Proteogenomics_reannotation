


library(magrittr)
library(ggplot2)

my_plots <- list()

#my_oa_f <- "H:/data/Synechocystis_6frame/Phylostratigraphy/Phylostrata_proteins.txt"

#my_oa_f <- "/mnt/storage/kxmna01/data/Synechocystis_6frame/OA_GSEA/2021-12-23_Phylostrata_for_OA.txt.OA.txt"

my_oa_f <- "/mnt/storage/kxmna01/data/Synechocystis_6frame/2022-03-02_Phylostrata_codon_characteristics/my_codon_stats.txt.GSEA.txt"

my_oa <- data.table::fread(
    input = my_oa_f, sep = "\t", quote = "",
    header = TRUE, stringsAsFactors = FALSE)

mandatory_path <- my_oa[
    my_oa$resource %in% c("Miscellaneous", "best_og_Subcategory"),][["pathway"]] %>%
    unique(.) %>%
    grep("pxd|ScyCode", ., value = TRUE, invert = TRUE)
mandatory_path <- my_oa[
    my_oa$resource == "Keywords",][["pathway"]] %>%
    unique(.) %>%
    grep("Signal|Plasmid|Carboxysome|Carbon|Ribosom|Photosyn|Thylak|Transposable", ., value = TRUE) %>%
    c(mandatory_path, .)

my_oa_filt <- my_oa %>%
    dplyr::filter(., !grepl("pxd|ScyCode", pathway)) %>%
    #dplyr::filter(., !resource %in% c(
    #    "CAZy", "PFAMs", "BiGG_Reaction",
    #    "EC level 2 name", "EC level 3 name",
    #    "KEGG_Reaction_Name", "KEGG_rclass_Name",
    #    "KEGG_Module_Name", "KEGG_brite_Name", "best_og_Category")) %>%
    dplyr::filter(., !is.na(padj) & padj <= 0.05) %>%
    dplyr::group_by(., set) %>%
    dplyr::arrange(., dplyr::desc(abs(NES))) %>%
    dplyr::mutate(
        ., Rank = 1:dplyr::n(),
        Selection = dplyr::case_when(
            Rank %in% c(1:100) ~ "Top",
            TRUE ~ "")) %>%
    dplyr::ungroup(.) %>%
    dplyr::filter(., pathway %in% mandatory_path | Selection == "Top")
#dplyr::filter(., Selection == "Top")

path_filt <- unique(my_oa_filt$pathway)

my_oa_final <- my_oa_filt %>%
    #dplyr::filter(., pathway %in% path_filt & padj <= 0.5) %>%
    dplyr::filter(., pathway %in% mandatory_path & padj <= 0.05) %>%
    dplyr::mutate(., padj_range = dplyr::case_when(
        padj <= 0.001 ~ "<= 0.001",
        padj <= 0.01 ~ "<= 0.01",
        padj <= 0.05 ~ "<= 0.05",
        TRUE ~ NA_character_
    )) %>%
    dplyr::filter(
        ., resource %in% c(
            "Miscellaneous", "best_og_Subcategory",
            "Characterization", "Keywords", "GOCC Term"))

my_oa_final$padj_range <- factor(
    x = my_oa_final$padj_range,
    levels = c("<= 0.001", "<= 0.01", "<= 0.05"),
    ordered= TRUE)
my_oa_final$resource <- factor(
    x = my_oa_final$resource,
    levels = sort(unique(my_oa_final$resource)),
    ordered= TRUE)
my_oa_final$pathway <- factor(
    x = my_oa_final$pathway,
    levels = sort(unique(my_oa_final$pathway)),
    ordered = TRUE)

my_palette <- grDevices::hcl.colors(n = 10, palette = "Fall")

my_plots[["Heatmap_oa_sig"]] <- ggplot(
    my_oa_final %>% dplyr::filter(., set %in% c("GC%", "GC3", "bias", "cai")),
    aes(x = set, y = pathway, fill = NES, size = padj_range)) +
    geom_point(colour = "white", shape = 21) +
    ggpubr::theme_pubr() +
    theme(
        legend.position = "right",
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
        panel.grid.major = element_line(colour = "grey", size = 0.5, linetype = "dotted")) +
    #facet_grid(rows = vars(resource), scales = "free_y", space = "free_y") +
    scale_size_manual(
        values = c(
            `<= 0.001` = 8, `<= 0.01` = 5,
            `<= 0.05` = 2)) +
    scale_fill_gradientn(
        colours = grDevices::hcl.colors(n = 9, palette = "Fall"),
        limits = c(-max(abs(my_oa_final$NES)), max(abs(my_oa_final$NES)))) +
    facet_grid(rows = vars(resource), scales = "free_y", space = "free_y")

my_date <- Sys.Date()
pdf(
    file = paste0(my_date, "_codon_bias_OA.pdf"),
    width = 10, height = 10)
my_plots
dev.off()


