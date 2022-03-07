


library(magrittr)
library(ggplot2)

my_plots <- list()
my_date <- Sys.Date()

my_f <- "/mnt/storage/kxmna01/data/Synechocystis_6frame/NoveltyExplain/ORF_novelty_reason_valid.RDS"
my_dir <- paste0("/mnt/storage/kxmna01/data/Synechocystis_6frame/", my_date, "_ORF_validation")

dir.create(my_dir)

my_data <- readRDS(
    file = my_f)

toplot <- my_data %>%
    dplyr::select(
        ., Proteins, Novel_peptide_count, RBS_fits_PeptideId,
        OperonID, OnlyIdBySite, PEPfilter, PX_novel_sequence) %>%
    dplyr::transmute(
        ., Proteins = Proteins,
        `Peptide 2+` = (!is.na(Novel_peptide_count) & Novel_peptide_count >= 2),
        `RBS valid` = tidyr::replace_na(data = RBS_fits_PeptideId, replace = FALSE),
        `TU valid` = (!is.na(OperonID) & OperonID != ""),
        `High quality` = (OnlyIdBySite & PEPfilter %in% c("class 1", "class 2")),
        `PX valid` = !is.na(PX_novel_sequence))

my_plots[["venn_noRBS"]] <- ggvenn::ggvenn(
    data = toplot,
    columns = c("Peptide 2+", "PX valid", "TU valid", "High quality"),
    fill_color = c("#387eb8", "#d1d2d4", "#e21e25", "#fbaf3f"))
my_plots[["venn_noPepCount"]] <- ggvenn::ggvenn(
    data = toplot,
    columns = c("RBS valid", "PX valid", "TU valid", "High quality"),
    fill_color = c("#387eb8", "#d1d2d4", "#e21e25", "#fbaf3f"))
my_plots[["venn_noPepCountorRNS"]] <- ggvenn::ggvenn(
    data = toplot,
    columns = c("PX valid", "TU valid", "High quality"),
    fill_color = c("#387eb8", "#d1d2d4", "#e21e25"))

my_data_final <- dplyr::left_join(x = my_data, y = toplot, by = "Proteins")

pdf(file = paste(my_dir, "Venn_ORF_validation.pdf", sep = "/"), width = 8, height = 8)
my_plots
dev.off()

data.table::fwrite(
    x = my_data_final,
    file = paste(my_dir, "Venn_ORF_validation.txt", sep = "/"),
    append = FALSE, quote = FALSE, sep = "\t",
    row.names = FALSE, col.names = TRUE)

save.image(paste(my_dir, "Venn_ORF_validation.RData", sep = "/"))


