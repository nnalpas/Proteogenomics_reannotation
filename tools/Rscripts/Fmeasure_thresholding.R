


rm(list = ls())

library(magrittr)
library(ggplot2)

my_data <- data.table::fread(
    input = "H:/data/Synechocystis_6frame/2022-03-04_ORF_validation/Venn_ORF_validation.txt",
    sep = "\t", quote = "\"", header = TRUE,
    stringsAsFactors = FALSE, colClasses = "character")

my_data %<>%
    dplyr::mutate(., Manual_validation = dplyr::case_when(
        grepl("^Likely false positive", Comment) ~ "FP",
        grepl("^Fit", Comment) ~ "TP",
        TRUE ~ "Uncertain"
    ))

my_data$Novel_min_PEP %<>%
    as.numeric(.)

my_data$Novel_max_Score %<>%
    as.numeric(.)

my_data %<>%
    dplyr::arrange(., dplyr::desc(Novel_min_PEP))

my_plots <- list()

my_plots[["pep_hist"]] <- ggplot(
    my_data, aes(
        x = -log(x = Novel_min_PEP, base = 10),
        fill = Manual_validation, colour = Manual_validation)) +
    geom_histogram(stat = "bin", bins = 100) +
    scale_x_log10() +
    ggpubr::theme_pubr()

my_plots[["pep_dens"]] <- ggplot(
    my_data, aes(
        x = -log(x = Novel_min_PEP, base = 10),
        fill = Manual_validation, colour = Manual_validation)) +
    geom_density() +
    scale_x_log10() +
    ggpubr::theme_pubr()

total_tp <- sum(my_data$Manual_validation == "TP")
total_fp <- sum(my_data$Manual_validation == "FP")

my_data_format <- data.frame()
for (x in 1:nrow(my_data)) {
    
    my_data_format <- my_data %>%
        dplyr::slice(., x:nrow(my_data)) %>%
        dplyr::mutate(
            ., 
            TP = sum(Manual_validation == "TP"),
            FP = sum(Manual_validation == "FP"),
            FN = total_tp - TP,
            TN = total_fp - FP) %>%
        dplyr::slice(., 1) %>%
        dplyr::bind_rows(my_data_format, .)

}

my_data_format %<>%
    dplyr::rowwise(.) %>%
    dplyr::mutate(
        Precision = TP / (TP + FP),
        Accuracy = (TP + TN) / (TP + TN + FP + FN),
        Sensitivity = TP / (TP + FN),
        Specificity = TN / (TN + FP),
        Fmeasure = (2 * Precision * Sensitivity) / (Precision + Sensitivity)) %>%
    dplyr::ungroup(.)

pep_threshold <- my_data_format[
    my_data_format$Fmeasure == max(my_data_format$Fmeasure), ][["Novel_min_PEP"]]

my_data_format %<>%
    dplyr::mutate(., `Fmeasure valid` = ifelse(
        Novel_min_PEP <= pep_threshold, "TRUE", "FALSE"))

my_plots[["Precision"]] <- ggplot(
    my_data_format, aes(x = -log(x = Novel_min_PEP, base = 10), y = Precision)) +
    geom_line(size = 1) +
    scale_x_log10() +
    ggpubr::theme_pubr()

my_plots[["Accuracy"]] <- ggplot(
    my_data_format, aes(x = -log(x = Novel_min_PEP, base = 10), y = Accuracy)) +
    geom_line(size = 1) +
    scale_x_log10() +
    ggpubr::theme_pubr()

my_plots[["Sensitivity"]] <- ggplot(
    my_data_format, aes(x = -log(x = Novel_min_PEP, base = 10), y = Sensitivity)) +
    geom_line(size = 1) +
    scale_x_log10() +
    ggpubr::theme_pubr()

my_plots[["Specificity"]] <- ggplot(
    my_data_format, aes(x = -log(x = Novel_min_PEP, base = 10), y = Specificity)) +
    geom_line(size = 1) +
    scale_x_log10() +
    ggpubr::theme_pubr()

my_plots[["Fmeasure"]] <- ggplot(
    my_data_format, aes(x = -log(x = Novel_min_PEP, base = 10), y = Fmeasure)) +
    geom_line(size = 1) +
    scale_x_log10() +
    ggpubr::theme_pubr()

my_data_format$`Peptide 2+` %<>%
    as.logical(.)
my_data_format$`PX valid` %<>%
    as.logical(.)
my_data_format$`TU valid` %>%
    as.logical(.)

my_plots[["venn_noRBS"]] <- ggvenn::ggvenn(
    data = my_data_format,
    columns = c("Peptide 2+", "PX valid", "TU valid"),
    fill_color = c("#387eb8", "#404040", "#e21e25"))

pdf(file = "H:/data/Synechocystis_6frame/2022-03-04_ORF_validation/Fmeasure_threshold.pdf", 8, 8)
my_plots
dev.off()

saveRDS(
    object = my_data_format,
    file = "H:/data/Synechocystis_6frame/2022-03-04_ORF_validation/Venn_ORF_validation_fmeasure.RDS")

data.table::fwrite(
    x = my_data_format,
    file = "H:/data/Synechocystis_6frame/2022-03-04_ORF_validation/Venn_ORF_validation_fmeasure.txt",
    append = FALSE, quote = FALSE, sep = "\t",
    row.names = FALSE, col.names = TRUE)


