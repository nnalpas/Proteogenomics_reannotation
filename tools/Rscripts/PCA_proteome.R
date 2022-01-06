


rm(list = ls())

library(magrittr)
library(ggplot2)

source("C:/Users/kxmna01/Documents/GitHub/Metaproteomics/tools/Rscripts/helper_functions.R")
source("C:/Users/kxmna01/Documents/GitHub/Metaproteomics/tools/Rscripts/helper_methods.R")

my_plots <- list()

#my_data_f <- "H:/data/Synechocystis_6frame/MQ_6frame_valid/combined/txt/proteinGroups.txt"
my_data_f <- "H:/data/Synechocystis_6frame/2022-01-05_Normalisation_pg/PG_normalised.txt"
ids_col <- "Protein IDs"

my_data <- data.table::fread(
    input = my_data_f, sep = "\t", quote = "", header = TRUE,
    stringsAsFactors = FALSE, colClasses = "character", na.strings = "NaN")

my_data_format <- my_data
my_data_format <- my_data_format[
    !grepl("CON__|REV__", my_data_format$`Protein IDs`), ]
rows <- my_data_format[[ids_col]]
my_data_format[[ids_col]] <- NULL
my_data_format <- as.data.frame(apply(my_data_format, 2, as.double))
my_data_format[!is.na(my_data_format) & my_data_format == 0] <- NA
rownames(my_data_format) <- rows

my_pheno <- data.frame(
    Samples = colnames(my_data_format),
    stringsAsFactors = FALSE
) %>%
    set_rownames(.[["Samples"]])

x <- as.matrix(my_data_format)
my_valid_entries <- data.frame(matrix(nrow = ncol(x), ncol = ncol(x))) %>%
    set_colnames(colnames(my_data_format)) %>%
    set_rownames(colnames(my_data_format))

for (i in seq_len(ncol(x))) {
    for (j in seq_len(ncol(x))) {
        x2 <- x[, i]
        y2 <- x[, j]
        ok <- complete.cases(x2, y2)
        my_valid_entries[i, j] <- sum(ok)
    }
}
my_valid_entries %<>%
    tibble::rownames_to_column(.data = ., var = "Rows") %>%
    tidyr::pivot_longer(data = ., cols = -Rows, names_to = "Cols")

my_plots[["Entries_valid"]] <- ggplot(
    my_valid_entries, aes(x = Cols, y = Rows, fill = value, label = value)) + 
    geom_tile() +
    geom_text() +
    ggpubr::theme_pubr() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
    scale_fill_gradientn(
        colours = grDevices::hcl.colors(n = 20, palette = "Fall"))

my_palette <- colorRampPalette(c("#387eb8", "#d1d2d4", "#e21e25"))(n = 30)
my_palette <- grDevices::hcl.colors(n = 20, palette = "Fall")

ms_cor <- cor(
    x = my_data_format, use = "pairwise.complete.obs", method = "spearman")
corrplot::corrplot(
    corr = ms_cor, type = "upper", method = "ellipse",
    addCoef.col = "darkgrey", tl.col = "black", diag = FALSE,
    main = paste(
        "\nSpearman corr. (",
        nrow(my_data_format), " entries)", sep = ""),
    tl.cex = 1.5,
    cl.cex = 1.5,
    number.cex = 1.2,
    col = my_palette)
my_plots[["Corr_Spearman"]] <- recordPlot()

my_samples_charact <- data.frame(
    name = colnames(my_data_format),
    shape = rep(
        x = c("solid", "dashed", "dotted"),
        length.out = length(colnames(my_data_format))),
    colour = rep(
        x = c("#387eb8", "#d1d2d4", "#e21e25", "#fbaf3f", "#3B3B3B", "#834F96"),
        each = 3, length.out = length(colnames(my_data_format))))

my_freq_total_entries <- apply(
    my_data_format, 2, function(x) { sum(!is.na(x)) }) %>%
    unlist(.) %>%
    data.table::data.table(name = names(.), Count = .) %>%
    dplyr::left_join(x = ., y = my_samples_charact) %>%
    dplyr::mutate(., Label = paste0(name, " (", shape, ")"))

my_plots[["Entries_total_freq"]] <- ggplot(
    my_freq_total_entries,
    aes(x = Label, y = Count, fill = colour, colour = colour, label = Count)) + 
    geom_bar(stat = "identity", position = "dodge") +
    geom_text(
        stat = "identity",
        position = position_dodge(width = 0.9), hjust = -0.2) +
    coord_flip() +
    ggpubr::theme_pubr() +
    scale_fill_identity() +
    scale_colour_identity()

my_freq_entries <- my_data_format %>%
    tibble::rownames_to_column(.data = ., var = "id") %>%
    tidyr::pivot_longer(data = ., cols = -id) %>%
    dplyr::filter(., !is.na(value)) %>%
    dplyr::left_join(x = ., y = my_samples_charact)

my_plots[["Entries_freq"]] <- ggplot(
    my_freq_entries,
    aes(x = value, group = name, colour = colour, linetype = shape)) + 
    geom_line(stat = "density", size = 1) +
    scale_x_log10() +
    ggpubr::theme_pubr() +
    scale_linetype_identity() +
    scale_color_identity()

my_data_filt <- my_data_format

my_data_filt <- my_data_filt[,
    apply(
        X = my_data_filt,
        MARGIN = 2,
        FUN = function(x) { sum(!is.na(x)) >= 1000})]

threshold <- 0.70
my_data_filt <- my_data_filt[
    apply(
        X = my_data_filt,
        MARGIN = 1,
        FUN = function(x) { sum(!is.na(x)) > floor(threshold*length(x)) }),]

#my_data_filt <- my_data_filt[
#    apply(
#        X = my_data_filt,
#        MARGIN = 1,
#        FUN = function(x) { !is.na(var(x, na.rm = TRUE)) })
#    ,]

my_freq_filt_entries <- my_data_filt %>%
    tibble::rownames_to_column(.data = ., var = "id") %>%
    tidyr::pivot_longer(data = ., cols = -id) %>%
    dplyr::filter(., !is.na(value)) %>%
    dplyr::left_join(x = ., y = my_samples_charact)

my_plots[["Entries_filt_freq"]] <- ggplot(
    my_freq_filt_entries,
    aes(x = value, group = name, colour = colour, linetype = shape)) + 
    geom_line(stat = "density", size = 1) +
    scale_x_log10() +
    ggpubr::theme_pubr() +
    scale_linetype_identity() +
    scale_color_identity()

my_freq_na_entries <- apply(
    my_data_filt, 2, function(x) { sum(is.na(x)) }) %>%
    unlist(.) %>%
    data.table::data.table(name = names(.), Count = .) %>%
    dplyr::left_join(x = ., y = my_samples_charact) %>%
    dplyr::mutate(., Label = paste0(name, " (", shape, ")"))

my_plots[["Entries_na_freq"]] <- ggplot(
    my_freq_na_entries,
    aes(x = Label, y = Count, fill = colour, colour = colour, label = Count)) + 
    geom_bar(stat = "identity", position = "dodge") +
    geom_text(
        stat = "identity",
        position = position_dodge(width = 0.9), hjust = -0.2) +
    coord_flip() +
    ggpubr::theme_pubr() +
    scale_fill_identity() +
    scale_colour_identity()

my_data_filt %<>%
    log10(.)# %>%
    #scale(x = ., center = TRUE, scale = FALSE)


my_freq_filt_entries_transform <- my_data_filt %>%
    tibble::rownames_to_column(.data = ., var = "id") %>%
    tidyr::pivot_longer(data = ., cols = -id) %>%
    dplyr::filter(., !is.na(value)) %>%
    dplyr::left_join(x = ., y = my_samples_charact)

my_plots[["Entries_filt_transform_freq"]] <- ggplot(
    my_freq_filt_entries_transform,
    aes(x = value, group = name, colour = colour, linetype = shape)) + 
    geom_line(stat = "density", size = 1) +
    scale_x_log10() +
    ggpubr::theme_pubr() +
    scale_linetype_identity() +
    scale_color_identity()

nb <- missMDA::estim_ncpPCA(X = t(my_data_filt), ncp.max = 5)
my_data_imput <- missMDA::imputePCA(X = t(my_data_filt), ncp = 3, scale = TRUE)

my_freq_imput_entries <- my_data_imput$completeObs %>%
    t(.) %>%
    data.table::data.table(.) %>%
    tibble::rownames_to_column(.data = ., var = "id") %>%
    tidyr::pivot_longer(data = ., cols = -id) %>%
    dplyr::filter(., !is.na(value)) %>%
    dplyr::left_join(x = ., y = my_samples_charact)

my_plots[["Entries_imput_freq"]] <- ggplot(
    my_freq_imput_entries,
    aes(x = value, group = name, colour = colour, linetype = shape)) + 
    geom_line(stat = "density", size = 1) +
    scale_x_log10() +
    ggpubr::theme_pubr() +
    scale_linetype_identity() +
    scale_color_identity()

res_pca <- FactoMineR::PCA(X = my_data_imput$completeObs, ncp = 5)
eigen <- factoextra::get_eig(res_pca)
factoextra::fviz_screeplot(res_pca, addlabels = TRUE)
var <- factoextra::get_pca_var(res_pca)
factoextra::fviz_pca_var(
    res_pca, col.var = "contrib", axes = c(1, 2),
    gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
    repel = TRUE, select.var = list(contrib = 100))
factoextra::fviz_contrib(res_pca, choice = "var", axes = 1, top = 100)
factoextra::fviz_contrib(res_pca, choice = "var", axes = 2, top = 100)
factoextra::fviz_contrib(res_pca, choice = "var", axes = 3, top = 100)
ind <- factoextra::get_pca_ind(res_pca)
factoextra::fviz_pca_ind(
    res_pca, col.ind = "cos2", axes = c(1, 2),
    gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
    repel = TRUE, max.overlaps = 5)
axis_lim <- max(abs(ind$coord))
my_plots[["PCA_biplot"]] <- factoextra::fviz_pca_biplot(
    res_pca, repel = FALSE, axes = c(1, 2),
    col.var = "cos2",
    gradient.cols = my_palette,
    select.var = list(cos2 = 30), max.overlaps = 5) +
    xlim(-axis_lim, axis_lim) +
    ylim(-axis_lim, axis_lim)

pdf(file = "QC_corr_pca.pdf", width = 8, height = 8)
my_plots
dev.off()

my_var_cor <- var$cor %>%
    data.frame(.) %>%
    tibble::rownames_to_column(.data = ., var = "Protein IDs") %>%
    tidyr::separate_rows(data = ., `Protein IDs`, sep = ";", convert = FALSE)
data.table::fwrite(
    x = my_var_cor, file = "pca_var_loadings.txt",
    append = FALSE, quote = FALSE, sep = "\t",
    row.names = FALSE, col.names = TRUE)

my_data_split <- my_data %>%
    tidyr::separate_rows(data = ., `Protein IDs`, sep = ";", convert = FALSE)
data.table::fwrite(
    x = my_data_split, file = "pca_values_imputed.txt",
    append = FALSE, quote = FALSE, sep = "\t",
    row.names = FALSE, col.names = TRUE)

save.image("QC_corr_pca.RData")


