


### Define working directory ---------------------------------------------

#
rm(list = ls())

#
library(magrittr)
library(dendextend)
library(ggplot2)

my_plots <- list()
my_cols <- c("#387eb8", "#d1d2d4", "#e21e25", "#fbaf3f", "#3B3B3B", "#834F96")

# 
my_pg_f <- "H:/data/Synechocystis_6frame/Kopf_2014_TU/Transcriptional_units_expr.txt"
my_annot_f <- "H:/data/Synechocystis_6frame/Custom_annotation/2021-12-21_Custom_Uniprot_Eggnog_annotations.txt"
#my_novel_f <- "H:/data/Synechocystis_6frame/2021-12-29_ORF_validation/Venn_ORF_validation.txt"
normalisation <- "zscore"
nmb_k <- 5
nmb_h <- 1
ids_col <- "id"
x_axis <- "Experiment"
y_axis <- "Normalised intensity"
feature_cols <- c(
    "#query_name", "Characterization", "chromosome", "aa_length",
    "TU ID", "Protein families", "GOBP Term", "GOCC Term", "GOMF Term",
    "EC number", "Interpro_NAME", "best_og_name", "best_og_Category", "EC",
    "GOBP Term.EggNOG", "GOCC Term.EggNOG", "GOMF Term.EggNOG")

my_annot <- data.table::fread(
    input = my_annot_f, sep = "\t", quote = "", header = TRUE,
    stringsAsFactors = FALSE, colClasses = "character")

feature <- my_annot %>%
    dplyr::filter(., Characterization == "Uncharacterized") %>%
    dplyr::select(., tidyselect::all_of(feature_cols)) %>%
    dplyr::rename(., Proteins = `#query_name`) %>%
    dplyr::mutate(., MainID = Proteins)

feature %<>%
    dplyr::mutate(., 
                  TU_identified = dplyr::case_when(
                      !is.na(`TU ID`) & `TU ID` != "" ~ "Annotated",
                      TRUE ~ "Not annotated"))

my_data <- data.table::fread(
    input = my_pg_f, sep = "\t", quote = "", header = TRUE,
    stringsAsFactors = FALSE, colClasses = "character", na.strings = "NaN")

#my_novel <- data.table::fread(
#    input = my_novel_f, sep = "\t", quote = "", header = TRUE,
#    stringsAsFactors = FALSE, colClasses = "character")



### Format the data ------------------------------------------------------

my_data_format <- my_data %>%
    dplyr::select(
        ., id, `15C`, `42C`, `-CO2`, Darkness, `-Fe`, `Hl`, `-N`, `-P`)
rows <- my_data_format[[ids_col]]
my_data_format[[ids_col]] <- NULL
my_data_format <- as.data.frame(apply(my_data_format, 2, as.double))
my_data_format[is.na(my_data_format)] <- 0
rownames(my_data_format) <- rows

my_freq_entries_all <- my_data_format %>%
    #t(.) %>%
    data.frame(.) %>%
    tibble::rownames_to_column(.data = ., var = "id") %>%
    tidyr::pivot_longer(data = ., cols = -id) %>%
    dplyr::filter(., !is.na(value) & value != 0)

my_plots[["Entries_freq_all"]] <- ggplot(
    my_freq_entries_all,
    aes(x = value, group = name, colour = name)) + 
    geom_line(stat = "density", size = 1) +
    scale_x_log10() +
    ggpubr::theme_pubr()# +
#scale_linetype_identity() +
#scale_color_identity()


#my_data_format <- my_data_format[
#    , apply(
#        X = my_data_format,
#        MARGIN = 2,
#        FUN = function(x) { sum(x > 0) >= 1500})]

patterns <- feature$Proteins %>%
    unique(.) %>%
    paste0(., collapse = "|") %>%
    paste0("(", ., ")(,|~)")
my_data_format <- my_data_format[grep(patterns, rownames(my_data_format)), ]

#my_data_format[my_data_format == 0] <- NA



### Tabulate dendogram data ----------------------------------------------

# Data normalisation
my_data_scaled <- log10(my_data_format+1)

#
if (normalisation == "zscore") {
    my_data_scaled %<>%
        t(.) %>%
        scale(., center = T, scale = T) %>%
        t(.)
}

my_freq_entries <- my_data_scaled %>%
    #t(.) %>%
    data.frame(.) %>%
    tibble::rownames_to_column(.data = ., var = "id") %>%
    tidyr::pivot_longer(data = ., cols = -id) %>%
    dplyr::filter(., !is.na(value) & value != 0)

my_plots[["Entries_freq"]] <- ggplot(
    my_freq_entries,
    aes(x = value, group = name, colour = name)) + 
    geom_line(stat = "density", size = 1) +
    #scale_x_log10() +
    ggpubr::theme_pubr()# +
#scale_linetype_identity() +
#scale_color_identity()

# 
pg_distance <- dist(x = my_data_scaled, method ="euclidean")
pg_hcluster <- hclust(d = pg_distance, method ="ward.D2")

#
pg_dendo <- as.dendrogram(pg_hcluster)

# Set the colors of branches
#cols_branches <- rainbow(n = nmb_k)
cols_branches <- c("#387eb8", "#d1d2d4", "#e21e25", "#fbaf3f", "#3B3B3B", "#834F96")[1:nmb_k]
pg_dendo <- color_branches(pg_dendo, k = nmb_k, col = cols_branches)
col_labels <- get_leaves_branches_col(pg_dendo)
col_labels <- col_labels[order(order.dendrogram(pg_dendo))]
leaves_label <- get_leaves_attr(pg_dendo, "label")

#
#my_palette <- colorRampPalette(c("blue", "purple", "red"))(n = 30)
#my_palette <- colorRampPalette(c("#2D9296", "#EDCA56", "#960835"))(n = 30)
#my_palette <- colorRampPalette(c("#387eb8", "#d1d2d4", "#e21e25"))(n = 30)

# 
h_cols <- data.frame(
    id = labels(pg_dendo),
    col = get_leaves_branches_col(pg_dendo),
    stringsAsFactors = FALSE)
h_cutree_k <- cutree(tree = pg_dendo, k = nmb_k) %>%
    as.data.frame(.) %>%
    set_colnames("k_cluster")
h_cutree_k$id <- row.names(h_cutree_k)
h_cutree_h <- cutree(tree = pg_dendo, h = nmb_h) %>%
    as.data.frame(.) %>%
    set_colnames("height")
h_cutree_h$id <- row.names(h_cutree_h)
h_cutree_k %<>%
    dplyr::left_join(x = ., y = h_cols) %>%
    dplyr::left_join(x = ., y = h_cutree_h)

#
pg_toplot <- my_data_scaled %>%
    as.data.frame(.) %>%
    tibble::rownames_to_column(.data = ., var = "id") %>%
    #dplyr::mutate(., id = as.integer(id)) %>%
    tidyr::gather(
        data = ., key = "Experiment", value = "scaled_value", -id,
        na.rm = FALSE, convert = TRUE) %>%
    #tidyr::separate(
    #    data = ., col = "key", into = c(x_axis, "Label"),
    #    sep = " ", remove = TRUE, convert = TRUE) %>%
    #dplyr::left_join(
    #    x = ., y = pg_format, by = c("id", x_axis, "Label")) %>%
    dplyr::left_join(
        x = ., y = h_cutree_k, by = "id") %>%
    #tidyr::separate_rows(data = ., id, sep = ";")# %>%
    dplyr::mutate(., Proteins = sub("~.+", "", id) %>% gsub(" |\"", "", .)) %>%
    tidyr::separate_rows(data = ., Proteins, sep = ",") %>%
    dplyr::left_join(
        x = .,
        y = my_annot,
        by = c("id" = "#query_name"))# %>%
    #dplyr::mutate(., size = 0.4)

#
pg_toplot %<>%
    dplyr::group_by(., k_cluster) %>%
    dplyr::mutate(
        ., count = dplyr::n_distinct(id)) %>%
    dplyr::ungroup(.) %>%
    dplyr::mutate(
        ., Cluster_Label = paste0("Cluster ", k_cluster, " (", count, " proteins)"))

# 
prot_label <- data.frame(Proteins = rownames(my_data_scaled)) %>%
    dplyr::mutate(., MainID = gsub("\"", "", Proteins))
    #dplyr::left_join(
    #    x = ., y = feature[, c("Proteins", "MainID")],
    #    by = "Proteins")# %>%
    #dplyr::mutate(., MainID = dplyr::case_when(
    #    is.na(MainID) ~ sub("(.*;|^)(s.+?)(;.*|$)", "\\2", Proteins),
    #    TRUE ~ MainID))
prot_label <- prot_label[["MainID"]] %>%
    set_names(prot_label[["Proteins"]])

#my_data_scaled[my_data_scaled == 0] <- NA



### Export graphs and data -----------------------------------------------

# 
pdf(
    file = paste0(
        "Clustering_k", nmb_k, "_", x_axis, "_",
        y_axis, "_", basename(my_pg_f), ".pdf"),
    width = 10, height = 10)

my_palette <- grDevices::hcl.colors(n = 20, palette = "Fall")

my_plots[["heatmap_scaled"]] <- gplots::heatmap.2(
    x = as.matrix(my_data_scaled),
    main = paste("Euclidean distance - Ward.D2 clustering"),
    trace = "none",
    margins = c(5, 7),
    col = my_palette,
    #breaks = col_breaks,
    dendrogram = "row",
    Rowv = pg_dendo,
    Colv = FALSE,
    labRow = prot_label,
    key.xlab = paste("Scaled", y_axis),
    cexRow = 0.4,
    cexCol = 1,
    na.rm = TRUE,
    RowSideColors = col_labels,
    colRow = col_labels)

my_plots

dev.off()

data.table::fwrite(
    x = function_annot_stats,
    file = "Uncharacterised_function_reannotation.txt",
    append = FALSE, quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)

#
to_export <- list(
    data = pg_toplot,
    dist = pg_distance,
    cluster = pg_hcluster,
    dendo = pg_dendo
)
saveRDS(
    object = to_export,
    file = paste0(
        "Clustering_k", nmb_k, "_", x_axis, "_",
        y_axis, "_", basename(my_pg_f)))

#
write.table(
    x = pg_toplot,
    file = paste0(
        "Clustering_k", nmb_k, "_", x_axis,
        "_", y_axis, "_", basename(my_pg_f)),
    append = FALSE, quote = FALSE, sep = "\t",
    row.names = FALSE, col.names = TRUE)

pg_toplot_oa <- pg_toplot %>%
    dplyr::select(., Proteins, k_cluster) %>%
    unique(.) %>%
    dplyr::mutate(., value = 1) %>%
    tidyr::pivot_wider(
        data = ., names_from = k_cluster, values_from = value)
pg_toplot_oa[is.na(pg_toplot_oa)] <- 0

#
write.table(
    x = pg_toplot_oa,
    file = paste0(
        "Clustering_k", nmb_k, "_", x_axis,
        "_", y_axis, "_for_OA_", basename(my_pg_f)),
    append = FALSE, quote = FALSE, sep = "\t",
    row.names = FALSE, col.names = TRUE)

save.image(sub(
    "\\.txt", ".RData",
    paste0(
        "Clustering_k", nmb_k, "_", x_axis,
        "_", y_axis, "_", basename(my_pg_f))))


