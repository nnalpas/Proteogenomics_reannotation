


library(magrittr)
library(ggplot2)

my_fasta_f <- "H:/data/Synechocystis_6frame/Phylostratigraphy/1148.faa"
my_evid_f <- "H:/data/Synechocystis_6frame/MQ_6frame_valid/combined/txt/evidence.txt"
my_phenodata_f <- "H:/data/Synechocystis_6frame/Phenodata/Scy004_resuscitation_gradseq.txt"
my_peptidome_f <- "H:/data/Synechocystis_6frame/Genome/Tryptic_digest_peptidome.txt"
my_abundance_cols <- c("Intensity L", "Intensity M", "Intensity H")

my_fasta <- seqinr::read.fasta(file = my_fasta_f, seqtype = "AA", as.string = TRUE)

my_evid <- data.table::fread(
    input = my_evid_f, sep = "\t", header = TRUE,
    stringsAsFactors = FALSE, colClasses = "character") %>%
    dplyr::filter(., Reverse != "+" & `Potential contaminant` != "+")

my_phenodata <- data.table::fread(
    input = my_phenodata_f, sep = "\t", header = TRUE,
    stringsAsFactors = FALSE, colClasses = "character")

my_peptidome <- data.table::fread(
    input = my_peptidome_f, sep = "\t", header = TRUE,
    stringsAsFactors = FALSE)

my_plots <- list()

my_evid_format <- my_evid %>%
    dplyr::filter(., `Raw file` %in% unique(my_phenodata$`Raw file`)) %>%
    dplyr::select(
        ., Sequence, Proteins, `Raw file`, Experiment, Charge, `m/z`,
        `Retention time`, PEP, `MS/MS count`, Score,
        tidyselect::all_of(my_abundance_cols),
        Reverse, `Potential contaminant`, id)

my_evid_format %<>%
    dplyr::left_join(x = ., y = my_phenodata, by = c("Raw file", "Experiment"))

my_evid_long <- my_evid_format %>%
    dplyr::select(
        ., id, Experiment, Replicate, Fraction, Proteins,
        tidyselect::all_of(my_abundance_cols)) %>%
    tidyr::pivot_longer(data = ., cols = tidyselect::all_of(my_abundance_cols)) %>%
    dplyr::mutate(., Label = sub("Intensity ", "", name)) %>%
    tidyr::unite(
        data = ., col = "Sample",
        Experiment, Replicate, Fraction, Label,
        sep = "_", remove = FALSE) %>%
    dplyr::mutate(., value = as.integer(value))

my_plots[["Raw_inten_box"]] <- ggplot(
    data = my_evid_long, aes(x = Sample, y = value, fill = Label)) +
    geom_boxplot() +
    ggpubr::theme_pubr() +
    theme(axis.text = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
    scale_y_log10()

my_plots[["Count_evid_box"]] <- ggplot(
    data = my_evid_long %>% dplyr::filter(., value > 0),
    aes(x = Sample, fill = Label)) +
    geom_histogram(stat = "count") +
    ggpubr::theme_pubr() +
    theme(axis.text = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
    scale_y_log10()

# Normalisation by dividing peptide intensity based on
# how many proteins it maps to
my_evid_norm <- my_evid_long %>%
    tidyr::separate_rows(data = ., Proteins, sep = ";") %>%
    dplyr::group_by(., id) %>%
    dplyr::mutate(., Protein_count = dplyr::n_distinct(Proteins)) %>%
    dplyr::ungroup(.) %>%
    dplyr:::mutate(., value_norm = value / Protein_count)

#my_norm <- my_evid_long %>%
#    dplyr::group_by(., id, Experiment, Replicate, Fraction) %>%
#    dplyr::mutate(
#        ., Mean = mean(value, na.rm = TRUE),
#        Stdev = sd(value, na.rm = TRUE),
#        Factor = Mean/value) %>%
#    dplyr::filter(., !is.na(Factor) & is.finite(Factor)) %>%
#    dplyr::group_by(., Experiment, Replicate, Label) %>%
#    dplyr::summarise(., Factor_median = median(Factor)) %>%
#    dplyr::ungroup(.)

#my_evid_norm <- my_evid_long %>%
#    dplyr::left_join(x = ., y = my_norm) %>%
#    dplyr::mutate(., value_norm = value * Factor_median)

my_plots[["Norm_inten_box"]] <- ggplot(
    data = my_evid_norm, aes(x = Sample, y = value_norm, fill = Label)) +
    geom_boxplot() +
    ggpubr::theme_pubr() +
    theme(axis.text = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
    scale_y_log10()

#my_plots[["Prot_evid_count_bar"]] <- ggplot(
#    data = my_evid_norm, aes(x = Sample, y = value_norm, fill = Label)) +
#    geom_boxplot() +
#    ggpubr::theme_pubr() +
#    theme(axis.text = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
#    scale_y_log10()

my_pg_norm <- my_evid_norm %>%
    dplyr::group_by(
        ., Experiment, Replicate, Proteins, name) %>%
    dplyr::summarise(
        ., value = sum(value, na.rm = TRUE),
        value_norm = sum(value_norm, na.rm = TRUE),
        count = dplyr::n()) %>%
    dplyr::ungroup(.) %>%
    dplyr::mutate(., Label = sub("Intensity ", "", name)) %>%
    tidyr::unite(
        data = ., col = "Sample",
        Experiment, Label, Replicate,
        sep = "_", remove = FALSE)

my_plots[["Norm_inten_pg_box"]] <- ggplot(
    data = my_pg_norm, aes(x = Sample, y = value_norm, fill = Label)) +
    geom_boxplot() +
    ggpubr::theme_pubr() +
    theme(axis.text = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
    scale_y_log10()

my_plots[["Raw_inten_pg_box"]] <- ggplot(
    data = my_pg_norm, aes(x = Sample, y = value, fill = Label)) +
    geom_boxplot() +
    ggpubr::theme_pubr() +
    theme(axis.text = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
    scale_y_log10()

my_plots[["Inten_pg_corr_box"]] <- ggplot(
    data = my_pg_norm, aes(x = value, y = value_norm, colour = count > 1)) +
    geom_point() +
    ggpubr::theme_pubr() +
    theme(axis.text = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
    scale_x_log10() +
    scale_y_log10() +
    facet_wrap(facets = vars(Sample))

#my_mw <- Peptides::mw(seq = my_fasta[[1]], monoisotopic = TRUE) %>%
#    data.frame(Proteins = names(my_fasta), MW = .)

my_ibaq <- my_pg_norm %>%
    #dplyr::mutate(., id = 1:dplyr::n()) %>%
    #tidyr::separate_rows(data = ., Proteins, sep = ";") %>%
    dplyr::left_join(x = ., y = my_peptidome) %>%
    #dplyr::group_by(., id) %>%
    #dplyr::mutate(
    #    ., Proteins = paste0(unique(Proteins), collapse = ";"),
    #    MW = max(MW, na.rm = TRUE)) %>%
    #dplyr::ungroup(.) %>%
    #unique(.) %>%
    dplyr::filter(., is.finite(Tryptic_count))

my_ibaq_custom <- my_ibaq %>%
    dplyr::mutate(
        ., iBAQ = value/Tryptic_count,
        iBAQ_norm = value_norm/Tryptic_count) %>%
    dplyr::group_by(., Sample, Experiment, Replicate, name) %>%
    dplyr::mutate(
        ., Total = sum(iBAQ, na.rm = TRUE),
        Total_norm = sum(iBAQ_norm, na.rm = TRUE)) %>%
    dplyr::ungroup(.)
my_ibaq_custom$value[
    !is.na(my_ibaq_custom$value) & my_ibaq_custom$value == 0] <- NA
my_ibaq_custom$value_norm[
    !is.na(my_ibaq_custom$value_norm) & my_ibaq_custom$value_norm == 0] <- NA
my_ibaq_custom$iBAQ[
    !is.na(my_ibaq_custom$iBAQ) & my_ibaq_custom$iBAQ == 0] <- NA
my_ibaq_custom$iBAQ_norm[
    !is.na(my_ibaq_custom$iBAQ_norm) & my_ibaq_custom$iBAQ_norm == 0] <- NA

my_plots[["Norm_ibaq_pg_box"]] <- ggplot(
    data = my_ibaq_custom, aes(x = Sample, y = iBAQ_norm, fill = Label)) +
    geom_boxplot() +
    ggpubr::theme_pubr() +
    theme(axis.text = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
    scale_y_log10()

my_plots[["Raw_ibaq_pg_box"]] <- ggplot(
    data = my_ibaq_custom, aes(x = Sample, y = iBAQ, fill = Label)) +
    geom_boxplot() +
    ggpubr::theme_pubr() +
    theme(axis.text = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
    scale_y_log10()

my_plots[["Norm_ibaq_pg_sum"]] <- ggplot(
    data = unique(my_ibaq_custom[, c("Sample", "Label", "Total_norm")]),
    aes(x = Sample, y = Total_norm, fill = Label)) +
    geom_bar(stat = "identity", position = "dodge") +
    ggpubr::theme_pubr() +
    theme(axis.text = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
    scale_y_log10()

my_plots[["Raw_ibaq_pg_sum"]] <- ggplot(
    data = unique(my_ibaq_custom[, c("Sample", "Label", "Total")]),
    aes(x = Sample, y = Total, fill = Label)) +
    geom_bar(stat = "identity", position = "dodge") +
    ggpubr::theme_pubr() +
    theme(axis.text = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
    scale_y_log10()

my_ibaq_custom_wide <- my_ibaq_custom %>%
    dplyr::select(., Sample, Proteins, Tryptic_count, value) %>%
    tidyr::pivot_wider(data = ., names_from = Sample, values_from = value)

data.table::fwrite(
    x = my_ibaq_custom_wide,
    file = "Scy004_resuscitation_Scy001_intensities_raw.txt",
    append = FALSE, quote = FALSE, sep = "\t",
    row.names = FALSE, col.names = TRUE)

my_ibaq_custom_wide <- my_ibaq_custom %>%
    dplyr::select(., Sample, Proteins, iBAQ_norm) %>%
    tidyr::pivot_wider(
        data = ., names_from = Sample, values_from = iBAQ_norm)

my_relative_ibaq_custom <- my_ibaq_custom_wide %>%
    tibble::column_to_rownames(.data = ., var = "Proteins") %>%
    as.matrix(.) %>%
    preprocessCore::normalize.quantiles(x = ., keep.names = TRUE) %>%
    as.data.frame(.) %>%
    tibble::rownames_to_column(.data = ., var = "Proteins")

data.table::fwrite(
    x = my_ibaq_custom_wide, file = "Scy004_resuscitation_Scy001_iBAQ_norm.txt",
    append = FALSE, quote = FALSE, sep = "\t",
    row.names = FALSE, col.names = TRUE)

data.table::fwrite(
    x = my_relative_ibaq_custom,
    file = "Scy004_resuscitation_Scy001_iBAQ_norm_quantile.txt",
    append = FALSE, quote = FALSE, sep = "\t",
    row.names = FALSE, col.names = TRUE)

my_ibaq_custom_wide <- my_ibaq_custom %>%
    dplyr::select(., Sample, Proteins, iBAQ) %>%
    tidyr::pivot_wider(data = ., names_from = Sample, values_from = iBAQ)

my_relative_ibaq_custom <- my_ibaq_custom_wide %>%
    tibble::column_to_rownames(.data = ., var = "Proteins") %>%
    as.matrix(.) %>%
    preprocessCore::normalize.quantiles(x = ., keep.names = TRUE) %>%
    as.data.frame(.) %>%
    tibble::rownames_to_column(.data = ., var = "Proteins")

data.table::fwrite(
    x = my_ibaq_custom_wide, file = "Scy004_resuscitation_Scy001_iBAQ_raw.txt",
    append = FALSE, quote = FALSE, sep = "\t",
    row.names = FALSE, col.names = TRUE)

data.table::fwrite(
    x = my_relative_ibaq_custom,
    file = "Scy004_resuscitation_Scy001_iBAQ_raw_quantile.txt",
    append = FALSE, quote = FALSE, sep = "\t",
    row.names = FALSE, col.names = TRUE)

data.table::fwrite(
    x = my_ibaq_custom,
    file = "Scy004_resuscitation_Scy001_intensities_long.txt",
    append = FALSE, quote = FALSE, sep = "\t",
    row.names = FALSE, col.names = TRUE)



