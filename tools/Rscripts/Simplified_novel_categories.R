


rm(list = ls())

library(magrittr)
library(ggplot2)

my_cols <- c("#387eb8", "#d1d2d4", "#e21e25", "#fbaf3f", "#3B3B3B", "#834F96")

my_data_f <- "H:/data/Pathogens_6frame/2022-07-21_ORF_validation/Venn_ORF_validation_fmeasure.RDS"

my_data_format <- readRDS(my_data_f)

#my_data_format <- my_data %>%
#    dplyr::mutate(., FinalCategory = dplyr::case_when(
#        grepl("match elsewhere", ORFNoveltyReason) ~ "Conflicting annotation",
#        TRUE ~ ORFNoveltyReason
#    ) %>% sub(";.+", "", .))

my_data_format$Reason <- factor(
    x = my_data_format$Reason,
    levels = c("Conflicting annotation", "Potential alternate start", "Potentially novel", "SAV"),
    ordered = TRUE)

toplot <- my_data_format %>%
    plyr::ddply(
        .data = ., .variables = "Reason",
        .fun = dplyr::summarise, Count = dplyr::n(), .drop = FALSE) %>%
    dplyr::mutate(., Type = "All novel")

toplot <- my_data_format %>%
    dplyr::filter(., `Fmeasure valid` == TRUE) %>%
    plyr::ddply(
        .data = ., .variables = "Reason",
        .fun = dplyr::summarise, Count = dplyr::n(), .drop = FALSE) %>%
    dplyr::mutate(., Type = "Quality filtered") %>%
    dplyr::bind_rows(toplot, .)

pl <- ggplot(
    toplot, aes(
        x = Count, y = Reason, fill = Reason,
        colour = Reason, alpha = Type, label = Count)) +
    geom_bar(stat = "identity", position = "dodge") +
    geom_text(position = position_dodge(width = 0.9), hjust = -0.3) +
    ggpubr::theme_pubr() +
    scale_fill_manual(
        values = my_cols[1:length(unique(toplot$Reason))]) +
    scale_colour_manual(
        values = my_cols[1:length(unique(toplot$Reason))]) +
    scale_alpha_manual(values = c(0.5, 0.8))

pdf(paste(dirname(my_data_f), "Updated_category_novels.pdf", sep = "/"), 10, 10)
pl
dev.off()


