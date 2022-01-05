


rm(list = ls())

library(magrittr)
library(ggplot2)

my_plots <- list()

#my_oa_f <- "H:/data/Synechocystis_6frame/Phylostratigraphy/Phylostrata_proteins.txt"

my_oa_f <- "H:/data/Synechocystis_6frame/2022-01-04_PCA/ibaq_values.txt.GSEA.txt"

my_oa <- data.table::fread(
    input = my_oa_f, sep = "\t", quote = "",
    header = TRUE, stringsAsFactors = FALSE)

mandatory_path <- my_oa %>%
    dplyr::filter(., `-log10_padj` > 2 & resource == "KEGG_Pathway_Name") %>%
    .[["pathway"]] %>%
    unique(.)

mandatory_path <- my_oa %>%
    dplyr::filter(
        ., resource %in% c(
            "Pathway", "GOBP Term", "GOCC Term", "GOMF Term",
            "Panther Pathway", "KEGG_Pathway_Name")) %>%
    dplyr::group_by(., set) %>%
    dplyr::arrange(., desc(`-log10_padj`)) %>%
    dplyr::slice(., 1:10) %>%
    .[["pathway"]] %>%
    unique(.) %>%
    grep("pxd|ScyCode", ., value = TRUE, invert = TRUE)
#, "pfam_description"

my_oa_final <- my_oa %>%
    dplyr::filter(
        ., pathway %in% mandatory_path &
            #resource == "KEGG_Pathway_Name" &
            padj <= 0.05) %>%
    dplyr::mutate(., padj_range = dplyr::case_when(
        padj <= 0.001 ~ "<= 0.001",
        padj <= 0.01 ~ "<= 0.01",
        padj <= 0.05 ~ "<= 0.05",
        TRUE ~ NA_character_
    )) %>%
    dplyr::group_by(., pathway) %>%
    dplyr::mutate(., median = median(NES)) %>%
    dplyr::ungroup(.) %>%
    dplyr::arrange(., dplyr::desc(median))

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
    levels = unique(my_oa_final$pathway),
    ordered = TRUE)
my_oa_final$set <- factor(
    x = my_oa_final$set,
    levels = unique(my_oa$set),
    ordered = TRUE)

my_palette <- colorRampPalette(c("#387eb8", "#d1d2d4", "#e21e25"))(n = 10)

my_plots[["Heatmap_oa_sig"]] <- ggplot(
    my_oa_final,
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
            `<= 0.001` = 12, `<= 0.01` = 8,
            `<= 0.05` = 4)) +
    scale_fill_gradientn(
        colours = grDevices::hcl.colors(n = 9, palette = "Fall"))

my_date <- Sys.Date()
pdf(
    file = paste0(my_date, "_phylostratigraphy_OA.pdf"),
    width = 10, height = 20)
my_plots
dev.off()


