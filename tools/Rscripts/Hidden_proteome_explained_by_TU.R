


library(magrittr)
library(ggplot2)

my_plots <- list()

my_transcri_f <- "H:/data/Synechocystis_6frame/Kopf_2014_TU/Transcriptional_units_expr.txt"

my_transcri <- data.table::fread(
    input = my_transcri_f, sep = "\t", quote = "", header = TRUE,
    stringsAsFactors = FALSE)

my_fasta_f <- c(
    "H:/data/Synechocystis_6frame/Genome/micro_proteins_Synechocystis_sp_PCC6803_20180419.fasta",
    "H:/data/Synechocystis_6frame/Genome/Synechocystis_sp_PCC_6803_cds_aa.fasta"
)

my_fasta <- lapply(my_fasta_f, function(x) {
    seqinr::read.fasta(file = x, seqtype = "AA", as.string = TRUE)
}) %>%
    unlist(., recursive = FALSE)

complete_proteome <- data.frame(Proteins = names(my_fasta), Proteome = TRUE)

my_pg_files <- c(
    SCyCode = "H:/data/Synechocystis_6frame/MQ_6frame_valid/combined/txt/proteinGroups.txt",
    pxd_nolabel = "H:/data/Synechocystis_6frame/ORF_validation/pxd_nolabel/combined/txt/proteinGroups.txt",
    pxd_005085 = "H:/data/Synechocystis_6frame/ORF_validation/pxd005085/combined/txt/proteinGroups.txt",
    pxd_014662_1 = "H:/data/Synechocystis_6frame/ORF_validation/pxd014662_1/combined/txt/proteinGroups.txt",
    pxd_014662_2 = "H:/data/Synechocystis_6frame/ORF_validation/pxd014662_2/combined/txt/proteinGroups.txt"
)

my_data <- lapply(my_pg_files, function(x) {
    data.table::fread(
        input = x, sep = "\t", quote = "", header = TRUE,
        stringsAsFactors = FALSE, colClasses = "character",
        na.strings = "NaN") %>%
        tidyr::separate_rows(data = ., `Protein IDs`, sep = ";") %>%
        dplyr::filter(., !grepl("^(REV__|CON__|chr|pca|pcb|psy)", `Protein IDs`))
        #.[["Protein IDs"]] %>%
        #grep("^(REV|CON|chr|pca|pcb|psy)__", ., invert = TRUE, value = TRUE) %>%
        #unique(.)
}) %>%
    plyr::ldply(., data.table::data.table, .id = "Processing") %>%
    dplyr::select(
        ., Processing, `Protein IDs`, Intensity, `Intensity L`,
        `Intensity M`, `Intensity H`) %>%
    dplyr::rowwise(.) %>%
    dplyr::mutate(., Intensity = dplyr::case_when(
        !is.na(Intensity) ~ as.double(Intensity),
        TRUE ~ sum(as.double(c(
            `Intensity L`, `Intensity M`, `Intensity H`)), na.rm = TRUE)
    )) %>%
    dplyr::select(., Processing, Proteins = `Protein IDs`, Intensity) %>%
    tidyr::separate(
        data = ., col = Processing, into = c("Type"),
        sep = "_", remove = FALSE, extra = "drop")

toplot <- my_data %>%
    dplyr::group_by(., Processing) %>%
    dplyr::summarise(., Count = dplyr::n_distinct(Proteins)) %>%
    dplyr::ungroup(.)

my_cols <- c("#387eb8", "#404040", "#e21e25", "#fbaf3f", "#d1d2d4") %>%
    set_names(toplot$Processing)

my_plots[["histogram_proc_identification"]] <- ggplot(
    data = toplot,
    mapping =  aes(x = Processing, y = Count, fill = Processing, colour = Processing)) +
    geom_bar(stat = "identity", position="dodge", alpha = 0.5) +
    ggpubr::theme_pubr() +
    scale_fill_manual(values = my_cols) +
    scale_colour_manual(values = my_cols)

toplot <- my_data %>%
    dplyr::select(., -Processing, -Intensity) %>%
    unique(.) %>%
    dplyr::mutate(., value = TRUE) %>%
    tidyr::pivot_wider(data = ., names_from = Type, values_from = value) %>%
    dplyr::left_join(x = complete_proteome, y = .) %>%
    dplyr::mutate_all(~tidyr::replace_na(data = ., replace = FALSE))

my_plots[["venn_identification_processing"]] <- ggvenn::ggvenn(
    data = toplot,
    columns = c("Proteome", "SCyCode", "pxd"),
    fill_color = c("#387eb8", "#d1d2d4", "#e21e25"))

never_identified <- list(
    never = toplot %>%
        dplyr::filter(Proteome == TRUE & SCyCode == FALSE & pxd == FALSE) %>%
        .[["Proteins"]],
    unique_SCyCode = toplot %>%
        dplyr::filter(Proteome == TRUE & SCyCode == TRUE & pxd == FALSE) %>%
        .[["Proteins"]],
    unique_PXD = toplot %>%
        dplyr::filter(Proteome == TRUE & SCyCode == FALSE & pxd == TRUE) %>%
        .[["Proteins"]]
)

for (x in names(never_identified)) {
    
    tmp <- my_fasta[never_identified[[x]]]
    tmp_names <- lapply(tmp, function(y) {attr(x = y, which = "Annot")}) %>%
        unlist(.) %>%
        sub("^>", "", .)
    seqinr::write.fasta(
        sequences = tmp, names = tmp_names,
        file.out = paste0("Never_identified_in_", x, ".fasta"),
        open = "w", as.string = TRUE)
    
}

toplot <- my_data %>%
    dplyr::filter(., Proteins %in% never_identified$unique_PXD) %>%
    dplyr::select(., -Type, -Intensity) %>%
    #unique(.) %>%
    dplyr::mutate(., value = TRUE) %>%
    tidyr::pivot_wider(data = ., names_from = Processing, values_from = value) %>%
    #dplyr::left_join(x = complete_proteome, y = .) %>%
    dplyr::mutate_all(~tidyr::replace_na(data = ., replace = FALSE))

my_plots[["venn_identification_PXD"]] <- ggvenn::ggvenn(
    data = toplot,
    columns = c("pxd_nolabel", "pxd_005085", "pxd_014662_1", "pxd_014662_2"),
    fill_color = c("#387eb8", "#d1d2d4", "#e21e25", "#fbaf3f"))

my_data_format <- my_data %>%
    dplyr::mutate(., Unique = dplyr::case_when(
        Proteins %in% c(
            never_identified$unique_SCyCode, never_identified$unique_PXD) ~ TRUE,
        TRUE ~ FALSE
    )) %>%
    dplyr::group_by(., Processing) %>%
    dplyr::arrange(., dplyr::desc(Intensity)) %>%
    dplyr::mutate(., Rank = 1:dplyr::n()) %>%
    dplyr::ungroup(.)

my_plots[["density_unique_identification"]] <- ggplot(data = my_data_format,
       mapping =  aes(x = Intensity, fill = Unique)) +
    geom_density(alpha = 0.5) +
    ggpubr::theme_pubr() +
    scale_x_log10() +
    facet_grid(rows = vars(Processing), scales = "free_y") +
    scale_fill_manual(values = c(`FALSE` = "#d1d2d4", `TRUE` = "#387eb8"))

my_plots[["histogram_unique_identification"]] <- ggplot(
    data = my_data_format,
    mapping =  aes(x = Intensity, fill = Unique)) +
    geom_histogram(position="identity", alpha = 0.5, colour = "black") +
    ggpubr::theme_pubr() +
    scale_x_log10() +
    facet_grid(rows = vars(Processing), scales = "free_y") +
    scale_fill_manual(values = c(`FALSE` = "#d1d2d4", `TRUE` = "#387eb8"))

my_data_oa <- my_data_format %>%
    dplyr::bind_rows(
        .,
        complete_proteome %>%
            dplyr::mutate(
                ., 
                Processing = "Proteome",
                Type = "Proteome",
                Unique = dplyr::case_when(
                    Proteins %in% never_identified$never ~ TRUE,
                    TRUE ~ FALSE
                ))) %>%
    dplyr::filter(., Unique == TRUE) %>%
    dplyr::select(., Processing, Type, Proteins) %>%
    dplyr::mutate(., value = 1)

my_data_oa_proc <- my_data_oa %>%
    dplyr::select(., -Type) %>%
    tidyr::pivot_wider(data = ., names_from = Processing, values_from = value) %>%
    dplyr::mutate_all(~tidyr::replace_na(data = ., replace = 0))

my_data_oa_type <- my_data_oa %>%
    dplyr::select(., -Processing) %>%
    unique(.) %>%
    tidyr::pivot_wider(data = ., names_from = Type, values_from = value) %>%
    dplyr::mutate_all(~tidyr::replace_na(data = ., replace = 0))

data.table::fwrite(
    x = my_data_oa_proc, file = "Hidden_proteome_for_oa_by_processing.txt",
    append = FALSE, quote = FALSE, sep = "\t",
    row.names = FALSE, col.names = TRUE)

data.table::fwrite(
    x = my_data_oa_type, file = "Hidden_proteome_for_oa_by_type.txt",
    append = FALSE, quote = FALSE, sep = "\t",
    row.names = FALSE, col.names = TRUE)

my_transcri_format <- my_transcri %>%
    dplyr::mutate(
        ., Hidden = grepl(
            paste0(never_identified$never, collapse = "|"), id)) %>%
    dplyr::filter(., grepl(".+~", id))

data.table::fwrite(
    x = my_transcri_format, file = "Hidden_proteome_found_by_TU.txt",
    append = FALSE, quote = FALSE, sep = "\t",
    row.names = FALSE, col.names = TRUE)

# sigB performed in Perseus
my_transcri_format_sig <- data.table::fread(
    input = "H:/data/Synechocystis_6frame/2022-02-02_Hidden_proteome/Hidden_proteome_found_by_TU_sigB.txt",
    sep = "\t", quote = "", header = TRUE,
    stringsAsFactors = FALSE) %>%
    dplyr::mutate(., Label = ifelse(
        `Max/Avg B significant` == "+" & Hidden == TRUE, id, NA))

my_plots[["density_unique_TU_identification"]] <- ggplot(
    data = my_transcri_format,
    mapping =  aes(x = Average, fill = Hidden)) +
    geom_density(alpha = 0.5) +
    ggpubr::theme_pubr() +
    scale_x_log10() +
    scale_fill_manual(values = c(`FALSE` = "#d1d2d4", `TRUE` = "#387eb8"))

my_plots[["histogram_unique_TU_identification"]] <- ggplot(
    data = my_transcri_format,
    mapping =  aes(x = Average, fill = Hidden)) +
    geom_histogram(position="identity", alpha = 0.5, colour = "black") +
    ggpubr::theme_pubr() +
    scale_x_log10() +
    scale_fill_manual(values = c(`FALSE` = "#d1d2d4", `TRUE` = "#387eb8"))


toplot <- data.frame(
    Proteins = never_identified$never) %>%
    dplyr::rowwise(.) %>%
    dplyr::mutate(
        ., TUid = dplyr::case_when(
            !any(grepl(Proteins, my_transcri_format$id)) ~ "No TU",
            any(grepl(paste0("^", Proteins, "~"), my_transcri_format$id)) ~ "Single gene in TU",
            TRUE ~ "Multiple genes in TU"
        ))

toplot$TUid <- factor(
    x = toplot$TUid, levels = c(
        "No TU", "Single gene in TU", "Multiple genes in TU"),
    ordered = TRUE)

my_cols <- c(
    `Single gene in TU` = "#d1d2d4", `Multiple genes in TU` = "#e21e25",
    `No TU` = "#387eb8")

my_plots[["histogram_proteins_per_TU"]] <- ggplot(
    data = toplot,
    mapping =  aes(x = TUid, fill = TUid, colour = TUid)) +
    geom_bar(stat = "count", position = "dodge", alpha = 0.9) +
    geom_text(stat = "count", aes(label = ..count..), vjust = -0.3) +
    ggpubr::theme_pubr() +
    scale_fill_manual(values = my_cols) +
    scale_colour_manual(values = my_cols)

my_plots[["volcano_TU_abundance"]] <- ggplot(
    data = my_transcri_format_sig,
    mapping =  aes(
        x = `Max/Avg`, y = Average,
        fill = Hidden, 
        size = `Max/Avg B significant`, shape = `Max/Avg B significant`,
        label = Label)) +
    geom_point() +
    ggrepel::geom_text_repel() +
    ggpubr::theme_pubr() +
    #scale_y_log10() +
    scale_fill_manual(values = c(`FALSE` = "#d1d2d4", `TRUE` = "#387eb8")) +
    scale_shape_manual(values = c(21, 22)) +
    scale_size_manual(values = c(2, 3.5))

pdf("Hidden_proteome_explained_by_TU.pdf", 10, 10)
my_plots
dev.off()

save.image("Hidden_proteome_explained_by_TU.RData")


