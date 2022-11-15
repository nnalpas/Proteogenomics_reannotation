


rm(list = ls())

library(magrittr)
library(ggplot2)

#score <- "NES"
score <- "score"

my_plots <- list()

#my_oa_f <- "H:/data/Synechocystis_6frame/Phylostratigraphy/Phylostrata_proteins.txt"

#my_oa_f <- "H:/data/Synechocystis_6frame/2022-01-07_iBAQ/PG_iBAQ.txt.GSEA.txt"

#my_oa_f <- "H:/data/Synechocystis_6frame/2022-01-20_PCA/pca_var_loadings.txt.GSEA.txt"

#my_oa_f <- "H:/data/Synechocystis_6frame/Kopf_2014_TU/Transcriptional_units_expr.txt.GSEA_withZeroes.txt"

#my_oa_f <- "H:/data/Synechocystis_6frame/2022-01-27_Copy_numbers/Scy004_copy_numbers_norm.txt.GSEA.txt"

#my_oa_f <- "/mnt/storage/kxmna01/data/Synechocystis_6frame/2022-03-02_Phylostrata_codon_characteristics/my_codon_stats.txt.GSEA.txt"

my_oa_f <- "C:/Users/nalpanic/SynologyDrive/Work/Abaumannii_trimeth/Analysis/Functional_enrichment/2022-11-15_Trimethylation (K)_multi_software_wide_for_OA.txt.OA.txt"

my_oa <- data.table::fread(
    input = my_oa_f, sep = "\t", quote = "",
    header = TRUE, stringsAsFactors = FALSE)

mandatory_path <- c()

my_oa_final <- my_oa %>%
    dplyr::filter(
        ., pathway %in% mandatory_path |
            padj <= 0.05) %>%
    dplyr::mutate(., padj_range = dplyr::case_when(
        padj <= 0.001 ~ "<= 0.001",
        padj <= 0.01 ~ "<= 0.01",
        padj <= 0.05 ~ "<= 0.05",
        TRUE ~ NA_character_
    )) %>%
    dplyr::group_by(., pathway) %>%
    dplyr::mutate(., median = median(!!as.name(score))) %>%
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

#my_palette <- colorRampPalette(c("#387eb8", "#d1d2d4", "#e21e25"))(n = 10)
my_palette <- grDevices::hcl.colors(n = 10, palette = "Fall")

my_plots <- lapply(unique(my_oa_final$resource), function(x) {
    ggplot(
        my_oa_final %>% dplyr::filter(., resource == x),
        aes(x = set, y = pathway, fill = !!as.name(score), size = padj_range)) +
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
            limits = c(-max(abs(my_oa_final[[score]])), max(abs(my_oa_final[[score]])))) +
        ggtitle(x)
}) %>%
    set_names(unique(my_oa_final$resource))

my_date <- Sys.Date()
pdf(
    file = paste0(my_date, "_enrichment_plots.pdf"),
    width = 10, height = 10)
my_plots
dev.off()


