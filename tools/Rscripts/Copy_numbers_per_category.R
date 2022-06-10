


rm(list = ls())

library(magrittr)
library(ggplot2)

my_plots <- list()

my_cols <- c("#387eb8", "#404040", "#e21e25", "#fbaf3f", "#d1d2d4")

#my_data_f <- "H:/data/Synechocystis_6frame/2022-01-27_Copy_numbers/Scy004_copy_numbers_norm.txt"
my_data_f <- "H:/data/Synechocystis_6frame/2022-02-14_iBAQ/Scy004_resuscitation_Scy001_iBAQ_norm.txt"

my_data <- data.table::fread(
    input = my_data_f, sep = "\t", quote = "",
    header = TRUE, stringsAsFactors = FALSE)

my_annot_f <- "H:/data/Synechocystis_6frame/Custom_annotation/2022-06-09_Custom_Uniprot_Eggnog_annotations.txt"

my_annot <- data.table::fread(
    input = my_annot_f, sep = "\t", quote = "",
    header = TRUE, stringsAsFactors = FALSE)

category <- my_annot$best_og_Subcategory %>%
    grep(";", ., invert = TRUE, value = TRUE) %>%
    table(.) %>%
    data.frame(.) %>%
    set_colnames(c("Category", "Count"))

my_data %<>%
    dplyr::filter(., Proteins %in% my_annot$`#query_name`)



# Focus on high abundant keyword

my_target <- c("Phycobilisome", "Ribosomal protein", "Thylakoid", "Carboxysome")

my_annot_format <- my_annot %>%
    dplyr::mutate(., Simple_Subcategory = dplyr::case_when(
        grepl(";", best_og_Subcategory) ~ "Multiple categories",
        best_og_Subcategory == "" ~ "Function unknown",
        best_og_Subcategory %in% category[category$Count >= 50, ][["Category"]] ~ best_og_Subcategory,
        TRUE ~ "Other categories"
    ),
    Simple_keyword = dplyr::case_when(
        grepl(my_target[[1]], Keywords) ~ my_target[[1]],
        grepl(my_target[[2]], Keywords) ~ my_target[[2]],
        grepl(my_target[[3]], Keywords) ~ my_target[[3]],
        grepl(my_target[[4]], Keywords) ~ my_target[[4]],
        #grepl(my_target[[5]], Keywords) ~ my_target[[5]],
        TRUE ~ "Others"),
    Preferred_name = Preferred_name)

my_data_format <- my_data %>%
    dplyr::filter(., !grepl("^chr|pca|pcb|psy", Proteins)) %>%
    tidyr::pivot_longer(data = ., cols = tidyselect::starts_with("SCy")) %>%
    dplyr::group_by(., Proteins) %>%
    dplyr::summarise(
        ., Median = median(value, na.rm = TRUE),
        Mean = mean(value, na.rm = TRUE)) %>%
    dplyr::ungroup(.) %>%
    dplyr::filter(., Mean > 0) %>%
    dplyr::left_join(
        x = ., y = my_annot_format[, c(
            "#query_name", "Simple_Subcategory", "Simple_keyword", "Preferred_name")],
        by = c("Proteins" = "#query_name"))

my_data_format$Simple_Subcategory <- factor(
    x = my_data_format$Simple_Subcategory,
    levels = sort(unique(my_data_format$Simple_Subcategory), decreasing = TRUE),
    ordered = TRUE)

my_target_col <- my_cols[1:length(my_target)] %>%
    set_names(my_target)

pl1 <- ggplot(
    my_data_format, aes(x = Mean, y = Simple_Subcategory)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(
        data = my_data_format %>%
            dplyr::filter(., Simple_keyword != "Others"),
        mapping = aes(x = Mean, y = Simple_Subcategory, fill = Simple_keyword),
        position = position_jitter(height = 0.1),
        shape = 21, size = 2) +
    #ggrepel::geom_text_repel(
    geom_text(
        data = my_data_format %>%
            dplyr::arrange(., dplyr::desc(Mean)) %>%
            dplyr::filter(., Simple_keyword != "Others") %>%
            dplyr::group_by(., Simple_Subcategory) %>%
            dplyr::filter(., Proteins != "") %>%
            dplyr::slice(., 1),
        mapping = aes(
            x = Mean, y = Simple_Subcategory,
            colour = Simple_keyword, label = Proteins)) +
    geom_text(
        data = my_data_format %>%
            dplyr::arrange(., Mean) %>%
            dplyr::filter(., Simple_keyword != "Others") %>%
            dplyr::group_by(., Simple_Subcategory) %>%
            dplyr::filter(., Proteins != "") %>%
            dplyr::slice(., 1),
        mapping = aes(
            x = Mean, y = Simple_Subcategory,
            colour = Simple_keyword, label = Proteins)) +
    ggpubr::theme_pubr() +
    scale_x_log10() +
    scale_fill_manual(values = my_target_col) +
    scale_colour_manual(values = my_target_col) +
    theme(legend.position = "none")

pl2 <- ggplot(
    my_data_format,
    aes(x = Mean, fill = Simple_keyword, colour = Simple_keyword)) +
    #geom_density(alpha = 0.3) +
    geom_density(aes(y = ..scaled..), alpha = 0.3) +
    ggpubr::theme_pubr() +
    scale_x_log10() +
    scale_fill_manual(values = c(my_target_col, `Others` = "#246E39")) +
    scale_colour_manual(values = c(my_target_col, `Others` = "#246E39"))

my_plots[["Copy_number_high"]] <- cowplot::plot_grid(
    pl2, pl1, nrow = 2,
    align = "hv", axis = "tblr", rel_heights = c(1, 4))



# Focus on low abundant or other keyword

my_target_other <- c(
    "Transposable element", "Protein kinase",
    "Potential alternate start", "Potentially novel", "Conflicting annotation")

my_annot_format_other <- my_annot %>%
    dplyr::mutate(., Simple_Subcategory = dplyr::case_when(
        grepl(";", best_og_Subcategory) ~ "Multiple categories",
        best_og_Subcategory == "" ~ "Function unknown",
        best_og_Subcategory %in% category[category$Count >= 50, ][["Category"]] ~ best_og_Subcategory,
        TRUE ~ "Other categories"
    ),
    Simple_keyword = dplyr::case_when(
        grepl(my_target_other[[1]], Keywords, fixed = T) ~ my_target_other[[1]],
        grepl(my_target_other[[2]], Miscellaneous, fixed = T) ~ my_target_other[[2]],
        grepl(my_target_other[[3]], Miscellaneous, fixed = T) ~ my_target_other[[3]],
        grepl(my_target_other[[4]], Miscellaneous, fixed = T) ~ my_target_other[[4]],
        grepl(my_target_other[[5]], Miscellaneous, fixed = T) ~ my_target_other[[5]],
        TRUE ~ "Others"),
    Preferred_name = Preferred_name)

my_data_format_other <- my_data %>%
    tidyr::pivot_longer(data = ., cols = tidyselect::starts_with("SCy")) %>%
    dplyr::group_by(., Proteins) %>%
    dplyr::summarise(
        ., Median = median(value, na.rm = TRUE),
        Mean = mean(value, na.rm = TRUE)) %>%
    dplyr::ungroup(.) %>%
    dplyr::filter(., Mean > 0) %>%
    dplyr::left_join(
        x = ., y = my_annot_format_other[, c(
            "#query_name", "Simple_Subcategory", "Simple_keyword", "Preferred_name", "Miscellaneous")],
        by = c("Proteins" = "#query_name")) %>%
    dplyr::filter(., !grepl("Low quality", Miscellaneous))

my_data_format_other$Simple_Subcategory <- factor(
    x = my_data_format_other$Simple_Subcategory,
    levels = sort(unique(my_data_format_other$Simple_Subcategory), decreasing = TRUE),
    ordered = TRUE)

my_target_other_col <- my_cols[1:length(my_target_other)] %>%
    set_names(my_target_other)

pl1 <- ggplot(
    my_data_format_other, aes(x = Mean, y = Simple_Subcategory)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(
        data = my_data_format_other %>%
            dplyr::filter(., Simple_keyword != "Others"),
        mapping = aes(x = Mean, y = Simple_Subcategory, fill = Simple_keyword),
        position = position_jitter(height = 0.1),
        shape = 21, size = 2) +
    #ggrepel::geom_text_repel(
    geom_text(
        data = my_data_format_other %>%
            dplyr::arrange(., Mean) %>%
            dplyr::filter(., Simple_keyword != "Others") %>%
            dplyr::group_by(., Simple_Subcategory) %>%
            dplyr::filter(., Proteins != "") %>%
            dplyr::slice(., 1),
        mapping = aes(
            x = Mean, y = Simple_Subcategory,
            colour = Simple_keyword, label = Proteins)) +#,
        #max.overlaps = 100) +
    geom_text(
        data = my_data_format_other %>%
            dplyr::arrange(., dplyr::desc(Mean)) %>%
            dplyr::filter(., Simple_keyword != "Others") %>%
            dplyr::group_by(., Simple_Subcategory) %>%
            dplyr::filter(., Proteins != "") %>%
            dplyr::slice(., 1),
        mapping = aes(
            x = Mean, y = Simple_Subcategory,
            colour = Simple_keyword, label = Proteins)) +
    ggpubr::theme_pubr() +
    scale_x_log10() +
    scale_fill_manual(values = my_target_other_col) +
    scale_colour_manual(values = my_target_other_col) +
    theme(legend.position = "none")

pl2 <- ggplot(
    my_data_format_other,
    aes(x = Mean, fill = Simple_keyword, colour = Simple_keyword)) +
    geom_density(aes(y = ..scaled..), alpha = 0.3) +
    ggpubr::theme_pubr() +
    scale_x_log10() +
    scale_fill_manual(values = c(my_target_other_col, `Others` = "#246E39")) +
    scale_colour_manual(values = c(my_target_other_col, `Others` = "#246E39"))

my_plots[["Copy_number_low"]] <- cowplot::plot_grid(
    pl2, pl1, nrow = 2,
    align = "hv", axis = "tblr", rel_heights = c(1, 4))

pdf("Copy_numbers_per_category.pdf", 12, 12)
my_plots
dev.off()


