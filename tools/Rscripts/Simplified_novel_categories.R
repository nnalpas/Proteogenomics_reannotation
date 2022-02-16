


library(magrittr)
library(ggplot2)

my_cols <- c("#387eb8", "#d1d2d4", "#e21e25", "#fbaf3f", "#3B3B3B", "#834F96")

my_data_f <- "H:/data/Synechocystis_6frame/2021-12-29_ORF_validation/Venn_ORF_validation.txt"

my_data <- data.table::fread(
    input = my_data_f, sep = "\t", quote = "",
    header = TRUE, stringsAsFactors = FALSE)

my_data_format <- my_data %>%
    dplyr::mutate(., FinalCategory = dplyr::case_when(
        ORFNoveltyReason %in% c("Potentially novel", "SAV", "Potential alternate start", "Potential alternate start") ~ ORFNoveltyReason,
        TRUE ~ "Conflicting annotation"
    ))

my_data_format$FinalCategory <- factor(
    x = my_data_format$FinalCategory,
    levels = unique(my_data_format$FinalCategory),
    ordered = TRUE)

toplot <- my_data_format %>%
    plyr::ddply(
        .data = ., .variables = "FinalCategory",
        .fun = dplyr::summarise, Count = dplyr::n(), .drop = FALSE) %>%
    dplyr::mutate(., Type = "All novel")

toplot <- my_data_format %>%
    dplyr::filter(., `High quality` == TRUE) %>%
    plyr::ddply(
        .data = ., .variables = "FinalCategory",
        .fun = dplyr::summarise, Count = dplyr::n(), .drop = FALSE) %>%
    dplyr::mutate(., Type = "Quality filtered") %>%
    dplyr::bind_rows(toplot, .)

pl <- ggplot(
    toplot, aes(
        x = Count, y = FinalCategory, fill = FinalCategory,
        colour = FinalCategory, alpha = Type, label = Count)) +
    geom_bar(stat = "identity", position = "dodge") +
    geom_text(position = position_dodge(width = 0.9), hjust = -0.3) +
    ggpubr::theme_pubr() +
    scale_fill_manual(
        values = my_cols[1:length(unique(toplot$FinalCategory))]) +
    scale_colour_manual(
        values = my_cols[1:length(unique(toplot$FinalCategory))]) +
    scale_alpha_manual(values = c(0.5, 0.8))

pdf("Updated_category_novels.pdf", 10, 10)
pl
dev.off()

