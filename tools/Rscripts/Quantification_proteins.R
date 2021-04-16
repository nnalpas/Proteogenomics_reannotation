


rm(list = ls())

library(magrittr)
library(MultiAssayExperiment)
library(ggplot2)

var_to_ignore <- c("Intensity", "Intensity L", "Intensity M", "Intensity H")

my_plots <- list()

my_file <- "H:/data/Synechocystis_6frame/SummarizedExp/MultiAssay.RDS"
my_f_files <- c(
    ref = "H:/data/Synechocystis_6frame/Genome/Synechocystis_sp_PCC_6803_cds_aa.fasta",
    micro_prot = "H:/data/Synechocystis_6frame/Genome/micro_proteins_Synechocystis_sp_PCC6803_20180419.fasta",
    never_id = "H:/data/Synechocystis_6frame/2020-12-02_Never_identified/Never_identified.fasta",
    novel_id = "H:/data/Synechocystis_6frame/Nuc_translation/Find0_Synechocystis_sp_PCC_6803_genome_FIXED_filter.fasta",
    wh_id = "H:/data/Synechocystis_6frame/Nuc_translation/Find0_Synechocystis_sp_PCC_6803_genome_FIXED_filter_Wolfgang_Hess_ORFs.fasta")

my_data <- readRDS(my_file)

my_fastas <- lapply(X = my_f_files, FUN = function(x) {
    seqinr::read.fasta(file = x, seqtype = "AA", as.string = T) %>%
        set_names(sub(".+\\|(.+)\\|.+", "\\1", names(.)))
}) %>%
    set_names(names(my_f_files))

non_sites_assays <- grep("Sites", names(experiments(my_data)), invert = TRUE)

abundance_df <- list()
#for (x in non_sites_assays) {
for (x in 1:length(experiments(my_data))) {
    
    toplot <- assays(experiments(my_data)[[x]])[[1]] %>%
        as.data.frame(.) %>%
        tibble::rownames_to_column(.data = ., var = "ID") %>%
        tidyr::pivot_longer(data = ., cols = -ID) %>%
        dplyr::filter(., !name %in% var_to_ignore) %>%
        dplyr::filter(., !is.na(value) & value != 0) %>%
        dplyr::group_by(., name) %>%
        dplyr::mutate(
            ., zscore = (value - mean(value, na.rm = TRUE))/sd(value, na.rm = TRUE)) %>%
        dplyr::ungroup(.)
    
    my_plots[[paste(names(experiments(my_data)[x]), "abund raw")]] <- ggplot(
        toplot,
        aes(x = name, y = value)) +
        geom_boxplot(position = "dodge") +
        ggpubr::theme_pubr() +
        #scale_y_continuous(trans = "log10") +
        theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
        ggtitle(paste0(names(experiments(my_data)[x]), ": raw abundance"))
    my_plots[[paste(names(experiments(my_data)[x]), "abund zscr")]] <- ggplot(
        toplot,
        aes(x = name, y = zscore)) +
        geom_boxplot(position = "dodge") +
        ggpubr::theme_pubr() +
        scale_y_continuous(trans = "log10") +
        theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
        ggtitle(paste0(names(experiments(my_data)[x]), ": zscored abundance"))
    
    toplot_summ <- toplot %>%
        dplyr::filter(., !name %in% c("iBAQ", "Intensity"))
    if (any(grepl("iBAQ", toplot$name))) {
        toplot_summ %<>%
            dplyr::filter(., grepl("iBAQ", name))
    }
    abundance_df[[names(experiments(my_data)[x])]] <- toplot_summ %>%
        dplyr::filter(., !is.na(value)) %>%
        tidyr::separate_rows(data = ., ID, sep = ";") %>%
        dplyr::group_by(., ID) %>%
        dplyr::summarise(
            ., Conditions = paste0(name, collapse = ";"),
            `Condition count` = dplyr::n_distinct(name),
            `Max abundance` = max(value),
            `Max condition` = paste0(name[value == max(value)], collapse = ";"))
    
}

abundance_df %<>%
    plyr::ldply(., data.table::data.table, .id = "Assay") %>%
    tidyr::separate(
        data = ., col = ID, into = c("protein id"),
        sep = "~", remove = FALSE) %>%
    tidyr::separate_rows(data = ., `protein id`, sep = ",|;")

protein_df <- lapply(my_fastas, names) %>%
    plyr::ldply(., data.table::data.table) %>%
    set_colnames(c("Fasta", "protein id")) %>%
    dplyr::left_join(x = ., y = abundance_df)

protein_df_spread <- protein_df %>%
    dplyr::filter(., !is.na(Assay)) %>%
    dplyr::mutate_all(~as.character(.)) %>%
    tidyr::pivot_longer(
        data = ., cols = c(
            `Conditions`, `Condition count`,
            `Max abundance`, `Max condition`)) %>%
    tidyr::unite(data = ., col = "key", Assay, name, sep = "__") %>%
    tidyr::pivot_wider(data = ., names_from = "key", values_from = "value")

# histogram best condition per fasta
toplot <- protein_df %>%
    dplyr::filter(., !is.na(Assay)) %>%
    tidyr::separate_rows(data = ., `Max condition`, sep = ";") %>%
    dplyr::group_by(., Fasta, Assay, `Max condition`) %>%
    dplyr::summarise(., Count = dplyr::n_distinct(`protein id`))

for (x in unique(toplot$Fasta)) {
    toplot_tmp <- toplot %>%
        dplyr::filter(., Fasta == x) %>%
        dplyr::arrange(., Assay, Count)
    toplot_tmp$`Max condition` <- factor(
        x = toplot_tmp$`Max condition`,
        levels = unique(toplot_tmp$`Max condition`),
        ordered = TRUE)
    my_plots[[paste(x, "best cond")]] <- ggplot(
        toplot_tmp,
        aes(x = `Max condition`, y = Count, fill = Assay, colour = Assay)) +
        geom_bar(stat = "identity", position = "dodge", alpha = 0.8) +
        ggpubr::theme_pubr() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
        ggtitle(paste0(x, ": all conditions"))
    my_plots[[paste(x, "top cond")]] <- ggplot(
        toplot_tmp %>%
            dplyr::group_by(., Assay) %>%
            dplyr::arrange(., dplyr::desc(Count)) %>%
            dplyr::slice(., 1:20),
        aes(x = `Max condition`, y = Count, fill = Assay, colour = Assay)) +
        geom_bar(stat = "identity", position = "dodge", alpha = 0.8) +
        ggpubr::theme_pubr() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
        ggtitle(paste0(x, ": top 20 conditions"))
}

pdf("Quantif_2_conditions.pdf", width = 18, height = 10)
my_plots
dev.off()

data.table::fwrite(
    x = protein_df, file = "Quantif_2_conditions.txt",
    append = FALSE, quote = FALSE, sep = "\t",
    row.names = FALSE, col.names = TRUE)

data.table::fwrite(
    x = protein_df_spread, file = "Quantif_2_conditions_spread.txt",
    append = FALSE, quote = FALSE, sep = "\t",
    row.names = FALSE, col.names = TRUE)


