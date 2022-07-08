


rm(list = ls())

library(magrittr)
library(ggplot2)

my_data <- data.table::fread(
    input = "H:/data/Srim_6frame_3rd_analysis/Novel_res_validation/ORF_validation.txt",
    sep = "\t", quote = "\"", header = TRUE,
    stringsAsFactors = FALSE)

my_data %<>%
    dplyr::mutate(., Reason = dplyr::case_when(
        grepl("match elsewhere", ORFNoveltyReason) ~ "Conflicting annotation",
        TRUE ~ ORFNoveltyReason
    ) %>% sub(";.+", "", .)) %>%
    dplyr::arrange(., dplyr::desc(as.integer(Novel_sum_MSMS_Count)), Novel_min_PEP)

my_data %<>%
    dplyr::mutate(., Manual_validation = dplyr::case_when(
        Novel_sum_Intensity == 0 | Novel_sum_MSMS_Count <= 2 | `Start valid` == FALSE ~ "FP",
        Reason == "Conflicting annotation" ~ "TP",
        TRUE ~ "Uncertain"
    ))
    #dplyr::mutate(., Manual_validation = dplyr::case_when(
    #    grepl("^Likely false positive", Comment) ~ "FP",
    #    grepl("^Fit", Comment) ~ "TP",
    #    TRUE ~ "Uncertain"
    #))

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
    !is.na(my_data_format$Fmeasure) &
    my_data_format$Fmeasure == max(my_data_format$Fmeasure, na.rm = TRUE), ] %>%
    .[["Novel_min_PEP"]] %>%
    max(.)

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

table(my_data_format$`Fmeasure valid`)

my_data_format$`Peptide 2+` %<>%
    as.logical(.)
my_data_format$`PX valid` %<>%
    as.logical(.)
my_data_format$`Start valid` %>%
    as.logical(.)

my_plots[["venn"]] <- ggvenn::ggvenn(
    data = my_data_format,
    columns = c("Peptide 2+", "PX valid", "Start valid"),
    fill_color = c("#387eb8", "#404040", "#e21e25"))

my_date <- Sys.Date()
outdir <- paste0(
    "H:/data/Srim_6frame_3rd_analysis/",
    my_date,
    "_ORF_validation/")
dir.create(outdir)

pdf(file = paste0(outdir, "Fmeasure_threshold.pdf"), 8, 8)
my_plots
dev.off()

saveRDS(
    object = my_data_format,
    file = paste0(outdir, "Venn_ORF_validation_fmeasure.RDS"))

data.table::fwrite(
    x = my_data_format,
    file = paste0(outdir, "Venn_ORF_validation_fmeasure.txt"),
    append = FALSE, quote = FALSE, sep = "\t",
    row.names = FALSE, col.names = TRUE)


