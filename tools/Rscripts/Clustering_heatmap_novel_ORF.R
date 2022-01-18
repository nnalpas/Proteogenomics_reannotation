



### Define working directory ---------------------------------------------

#
rm(list = ls())

#
library(dplyr)
library(tidyr)
library(plyr)
library(magrittr)
library(MSnbase)
library(data.table)
library(gplots)
library(dendextend)
library(ggplot2)
library(ggrepel)
library(ggpubr)
library(mgcv)

my_plots <- list()

# 
my_pg_f <- "H:/data/Synechocystis_6frame/MQ_6frame/combined/txt/proteinGroups.txt"
my_evid_f <- "H:/data/Synechocystis_6frame/Novel_res/Group_evidence.RDS"
my_novel_f <- "H:/data/Synechocystis_6frame/2021-12-29_ORF_validation/Venn_ORF_validation.txt"
normalisation <- "zscore"
nmb_k <- 5
nmb_h <- 1
x_axis <- "Experiment"
y_axis <- "iBAQ"
feature_cols <- c(
    "Proteins", "Novel_peptide_count", "Novel_sum_MSMS_Count",
    "ORFNoveltyReason", "MainID", "Peptide 2+", "RBS valid",
    "TU valid", "High quality", "PX valid")

my_pg <- data.table::fread(
    input = my_pg_f, sep = "\t", quote = "", header = TRUE,
    stringsAsFactors = FALSE, colClasses = "character")

my_evid <- readRDS(my_evid_f)

my_evid_format <- my_evid %>%
    dplyr::filter(
        ., Reverse != "+" & `Potential contaminant` != "+" &
            Database == "Novel") %>%
    dplyr::select(., `Protein group IDs`, Intensity, Experiment) %>%
    tidyr::separate_rows(data = ., `Protein group IDs`, sep = ";") %>%
    dplyr::group_by(., `Protein group IDs`, Experiment) %>%
    dplyr::summarise(., Intensity = sum(Intensity, na.rm = TRUE)) %>%
    dplyr::ungroup(.) %>%
    dplyr::mutate(., Intensity = ifelse(Intensity > 0, 1, 0))

my_novel <- data.table::fread(
    input = my_novel_f, sep = "\t", quote = "", header = TRUE,
    stringsAsFactors = FALSE, colClasses = "character")

feature <- my_novel %>%
    dplyr::select(., tidyselect::all_of(feature_cols)) %>%
    dplyr::mutate(., Reason = dplyr::case_when(
        grepl("start", ORFNoveltyReason) ~ "New start",
        grepl("termination", ORFNoveltyReason) ~ "New stop",
        grepl("novel", ORFNoveltyReason) ~ "Novel",
        TRUE ~ ORFNoveltyReason)) %>%
    dplyr::mutate(., MainID = ifelse(MainID != "", MainID, Proteins))

pg_abundance <- my_pg %>%
    dplyr::filter(., Reverse != "+" & `Potential contaminant` != "+") %>%
    dplyr::select(., id, `Protein IDs`, tidyselect::starts_with("iBAQ ")) %>%
    dplyr::select(., -tidyselect::matches(" (L|M|H)"), -`iBAQ peptides`) %>%
    tidyr::separate_rows(data = ., `Protein IDs`, sep = ";") %>%
    dplyr::filter(., `Protein IDs` %in% my_novel$Proteins) %>%
    tidyr::pivot_longer(
        data = ., cols = -tidyselect::all_of(c("id", "Protein IDs")),
        names_to = "Param", values_to = "iBAQ") %>%
    tidyr::separate(
        data = ., col = "Param",
        into = c("Quantification", "Experiment"),
        sep = " ") %>%
    dplyr::mutate(., iBAQ = as.integer(iBAQ)) %>%
    dplyr::left_join(
        x = ., y = my_evid_format,
        by = c("id" = "Protein group IDs", "Experiment")) %>%
    dplyr::mutate(
        ., iBAQ = tidyr::replace_na(data = iBAQ, replace = 0),
        Intensity = tidyr::replace_na(data = Intensity, replace = 0)) %>%
    dplyr::mutate(., iBAQ_norm = iBAQ * Intensity) %>%
    dplyr::select(., `Protein IDs`, Experiment, iBAQ_norm) %>%
    tidyr::pivot_wider(
        data = ., names_from = "Experiment", values_from = "iBAQ_norm") %>%
    as.data.frame(.)



### Format the data ------------------------------------------------------

rownames(pg_abundance) <- pg_abundance[["Protein IDs"]]
pg_abundance$`Protein IDs` <- NULL
pg_abundance %<>%
    dplyr::mutate_all(~as.double(.)) %>%
    base::as.matrix(.)
#pg_abundance[pg_abundance == "NaN"] <- 0
no_var <- apply(pg_abundance, 1, var)
pg_abundance <- pg_abundance[no_var != 0, ]



### Tabulate dendogram data ----------------------------------------------

# Data normalisation
#pg_scaled <- pg_abundance

#
#if (!is.null(normalisation)) {
#    pg_scaled %<>%
#        t(.) %>%
#        scale(., center = T, scale = T) %>%
#        t(.)
#}

pg_scaled <- log10(pg_abundance+1)

# 
pg_distance <- dist(x = pg_scaled, method ="euclidean")
pg_hcluster <- hclust(d = pg_distance, method ="ward.D2")

#
pg_dendo <- as.dendrogram(pg_hcluster)

# Set the colors of branches
cols_branches <- rainbow(n = nmb_k)
pg_dendo <- color_branches(pg_dendo, k = nmb_k, col = cols_branches)
col_labels <- get_leaves_branches_col(pg_dendo)
col_labels <- col_labels[order(order.dendrogram(pg_dendo))]
leaves_label <- get_leaves_attr(pg_dendo, "label")

#
#my_palette <- colorRampPalette(c("blue", "purple", "red"))(n = 30)
#my_palette <- colorRampPalette(c("#2D9296", "#EDCA56", "#960835"))(n = 30)
my_palette <- colorRampPalette(c("#387eb8", "#d1d2d4", "#e21e25"))(n = 30)

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
    dplyr::left_join(x = ., y = h_cutree_h) %>%
    dplyr::left_join(
        x = .,
        y = feature, by = c("id" = "Proteins"))

#
pg_toplot <- pg_scaled %>%
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
    dplyr::mutate(., size = 0.4)

#
pg_toplot %<>%
    dplyr::group_by(., k_cluster) %>%
    dplyr::mutate(
        ., count = dplyr::n_distinct(id)) %>%
    dplyr::ungroup(.) %>%
    dplyr::mutate(
        ., Cluster_Label = paste0("Cluster ", k_cluster, " (", count, " proteins)"))

# 
prot_label <- data.frame(Proteins = rownames(pg_scaled)) %>%
    dplyr::left_join(
        x = ., y = feature[, c("Proteins", "MainID")],
        by = "Proteins")
prot_label <- prot_label[["MainID"]] %>%
    set_names(prot_label[["Proteins"]])

pg_scaled[pg_scaled == 0] <- NA



# Further novel ORF description ------------------------------------------

my_identifi_per_exp <- apply(pg_scaled, 2, function(x) {
    sum(!is.na(x))
}) %>%
    data.frame(Experiment = names(.), Count = .)

my_identifi_per_exp$Experiment <- factor(
    x = my_identifi_per_exp$Experiment,
    levels = colnames(pg_scaled), ordered = TRUE)

feature_format <- leaves_label %>%
    data.frame(Proteins = .) %>%
    dplyr::left_join(x = ., y = feature)
feature_format$MainID <- factor(
    x = feature_format$MainID,
    levels = unique(feature_format$MainID), ordered = TRUE)

quart_msms <- summary(as.integer(feature_format$Novel_sum_MSMS_Count))
feature_format_heat <- feature_format %>%
    dplyr::mutate(., MSMS_Count = dplyr::case_when(
        as.integer(Novel_sum_MSMS_Count) <= quart_msms[["1st Qu."]] ~ "25%",
        as.integer(Novel_sum_MSMS_Count) <= quart_msms[["Median"]] ~ "50%",
        as.integer(Novel_sum_MSMS_Count) <= quart_msms[["3rd Qu."]] ~ "75%",
        TRUE ~ "100%")) %>%
    dplyr::select(
        ., MainID, Reason, MSMS_Count, `Peptide 2+`, `RBS valid`, `TU valid`,
        `High quality`, `PX valid`) %>%
    tidyr::pivot_longer(data = ., cols = -MainID)

feature_format_heat$name <- factor(
    x = feature_format_heat$name,
    levels = unique(feature_format_heat$name), ordered = TRUE)

my_cols <- c(
    `New start` = "#387eb8", `New stop` = "#fbaf3f", `SAV` = "#d1d2d4", `Novel` = "#e21e25",
    `TRUE` = "darkgrey", `FALSE` = "white",
    `25%` = "#7F98EB", `50%` = "#263BBF", `75%` = "#EB7775", `100%` = "#B8110E")

my_plots[["histogram_exp"]] <- ggplot(
    my_identifi_per_exp, aes(x = Experiment, y = Count)) +
    geom_bar(
        stat = "identity", position = "dodge", colour = "black", alpha = 0.6) +
    ggpubr::theme_pubr() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

my_plots[["histogram_orf"]] <- ggplot(
    feature_format,
    aes(x = MainID, y = as.integer(Novel_sum_MSMS_Count))) +
    geom_bar(
        stat = "identity", position = "dodge", colour = "black", alpha = 0.6) +
    ggpubr::theme_pubr() +
    coord_flip() +
    ylab("MS/MS count")

my_plots[["heat_orf_properties"]] <- ggplot(
    feature_format_heat,
    aes(x = name, y = MainID, fill = value)) +
    geom_tile(colour = "white") +
    ggpubr::theme_pubr() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
    scale_fill_manual(values = my_cols)



### Export graphs and data -----------------------------------------------

# 
pdf(
    file = paste0(
        "Clustering_k", nmb_k, "_", x_axis, "_",
        y_axis, "_", basename(my_pg_f), ".pdf"),
    width = 10, height = 10)

my_plots[["heatmap"]] <- heatmap.2(
    x = pg_scaled,
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

#
to_export <- list(
    data = pg_toplot,
    heatmap = my_heatmap,
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

save.image(sub(
    "\\.txt", ".RData",
    paste0(
        "Clustering_k", nmb_k, "_", x_axis,
        "_", y_axis, "_", basename(my_pg_f))))


