


library(magrittr)
library(ggplot2)

my_plots <- list()
my_cols <- c("#387eb8", "#d1d2d4", "#e21e25", "#fbaf3f", "#3B3B3B", "#834F96")

wkdir <- "/mnt/storage/kxmna01/data/Synechocystis_6frame/"



### Phosphosites ---------------------------------------------------------

my_data_f <- paste0(wkdir, "MQ_6frame_valid/combined/txt/summary/PhosphoSites_nonredundant.txt")

my_data <- data.table::fread(
    input = my_data_f, sep = "\t", quote = "", header = TRUE,
    stringsAsFactors = FALSE, colClasses = "character")

my_data_filt <- my_data %>%
    #dplyr::filter(., Reverse != "+" & `Potential contaminant` != "+")
    dplyr::filter(., Reverse != "+" & Contaminant != "+")

count_phospho_prot <- length(unique(unlist(strsplit(
    x = my_data_filt$Proteins, split = ";"))))

my_plots[["Localisation_distribution"]] <- ggplot(
    my_data_filt, aes(x = as.numeric(`Localization.prob`))) +
    geom_histogram() +
    ggpubr::theme_pubr()

phospho_stats <- function(x) {
    
    my_total_stats <- data.frame()
    
    my_total_stats <- x %>%
        dplyr::mutate(., Category = "Total") %>%
        dplyr::group_by(., Category) %>%
        dplyr::summarise(., Count = dplyr::n()) %>%
        dplyr::bind_rows(my_total_stats, .)
    
    my_total_stats <- x %>%
        dplyr::filter(
            ., !is.na(`Localization.prob`) & `Localization.prob` > 0.95) %>%
        dplyr::mutate(., Category = "Localised >0.95") %>%
        dplyr::group_by(., Category) %>%
        dplyr::summarise(., Count = dplyr::n()) %>%
        dplyr::bind_rows(my_total_stats, .)
    
    my_total_stats <- x %>%
        dplyr::mutate(., Category = `Amino.acid`) %>%
        dplyr::group_by(., Category) %>%
        dplyr::summarise(., Count = dplyr::n()) %>%
        dplyr::bind_rows(my_total_stats, .)
    
    my_total_stats <- x %>%
        dplyr::mutate(., Category = `Amino.acid`) %>%
        dplyr::group_by(., Category) %>%
        dplyr::summarise(., Count = dplyr::n()) %>%
        dplyr::bind_rows(my_total_stats, .)
    
    my_total_stats %<>%
        dplyr::ungroup(.) %>%
        dplyr::arrange(., dplyr::desc(Count))
    
    my_total_stats$Category <- factor(
        x = my_total_stats$Category,
        levels = unique(my_total_stats$Category),
        ordered = TRUE)
    
    pl <- ggplot(
        my_total_stats, aes(
            x = Category, y = Count, fill = Category,
            colour = Category, label = Count)) +
        geom_bar(stat = "identity", position = "dodge", alpha = 0.8) +
        geom_text(
            stat = "identity", position = position_dodge(width = 0.9),
            vjust = -0.5) +
        ggpubr::theme_pubr() +
        scale_fill_manual(values = my_cols) +
        scale_colour_manual(values = my_cols)
    
    return(list(data = my_total_stats, plot = pl))
    
}

all_stats <- phospho_stats(x = my_data_filt)

my_plots[["Phospho_stats_all"]] <- all_stats[["plot"]]

my_gr_phospho_f <- paste0(wkdir, "GRanges/Phospho (STY)Sites_grange.RDS")

my_gr_phospho <- readRDS(my_gr_phospho_f)

my_gr_ref_f <- paste0(wkdir, "GRanges/Ref_prot_grange.RDS")

my_gr_ref <- readRDS(my_gr_ref_f)

my_gr_phospho_notref <- IRanges::subsetByOverlaps(
    x = my_gr_phospho, ranges = my_gr_ref, invert = TRUE)

my_novel_stats <- data.frame(
    Names = names(my_gr_phospho_notref)) %>%
    tidyr::separate(
        data = ., col = Names,
        into = c("Proteins", "Positions.within.proteins", "Amino.acid"),
        sep = "~") %>%
    dplyr::inner_join(x = ., my_data_filt)

all_stats_novel <- phospho_stats(x = my_novel_stats)

my_plots[["Phospho_stats_novel"]] <- all_stats_novel[["plot"]]

my_reason_f <- paste0(wkdir, "2021-12-29_ORF_validation/Venn_ORF_validation.txt")

my_reason <- data.table::fread(
    input = my_reason_f, sep = "\t", quote = "", header = TRUE,
    stringsAsFactors = FALSE, colClasses = "character")

my_novel_hq <- my_reason %>%
    dplyr::filter(., `High quality` == TRUE) %>%
    .[["Proteins"]]

my_novel_stats_hq <- my_novel_stats %>%
    dplyr::filter(., Proteins %in% my_novel_hq)

all_stats_novel_hq <- phospho_stats(x = my_novel_stats_hq)

my_plots[["Phospho_stats_novel_HQ"]] <- all_stats_novel_hq[["plot"]]



### Phosphopeptides ------------------------------------------------------

my_evid_f <- paste0(wkdir, "MQ_6frame_valid/combined/txt/evidence.txt")

my_evid <- data.table::fread(
    input = my_evid_f, sep = "\t", quote = "", header = TRUE,
    stringsAsFactors = FALSE, colClasses = "character")

my_evid_filt <- my_evid %>%
    dplyr::filter(., grepl("Phospho", Modifications)) %>%
    dplyr::select(
        ., Sequence, `Modified sequence`,
        `Phospho (STY)`, `Phospho (STY) site IDs`,
        tidyselect::starts_with("Missed cleavages")) %>%
    unique(.)

my_evid_filt_long <- my_evid_filt %>%
    dplyr::mutate(., ID = 1:dplyr::n()) %>%
    tidyr::separate_rows(data = ., `Phospho (STY) site IDs`, sep = ";") %>%
    dplyr::filter(., `Phospho (STY) site IDs` %in% my_data$id) %>%
    dplyr::arrange(., `Phospho (STY) site IDs`)

phospho_evid_stats <- my_evid_filt_long %>%
    dplyr::group_by(., `Phospho (STY)`) %>%
    dplyr::summarise(., Count = dplyr::n_distinct(`Modified sequence`)) %>%
    dplyr::ungroup(.)

my_plots[["Phospho_multi_all"]] <- ggplot(
    phospho_evid_stats, aes(
        x = `Phospho (STY)`, y = Count, fill = `Phospho (STY)`,
        colour = `Phospho (STY)`, label = Count)) +
    geom_bar(stat = "identity", position = "dodge", alpha = 0.8) +
    geom_text(
        stat = "identity", position = position_dodge(width = 0.9),
        vjust = -0.5) +
    ggpubr::theme_pubr() +
    scale_fill_manual(values = my_cols) +
    scale_colour_manual(values = my_cols)

novel_evid_stats <- my_evid_filt_long %>%
    dplyr::filter(., `Phospho (STY) site IDs` %in% my_novel_stats$id) %>%
    dplyr::group_by(., `Phospho (STY)`) %>%
    dplyr::summarise(., Count = dplyr::n_distinct(`Modified sequence`)) %>%
    dplyr::ungroup(.)

my_plots[["Phospho_multi_novel"]] <- ggplot(
    novel_evid_stats, aes(
        x = `Phospho (STY)`, y = Count, fill = `Phospho (STY)`,
        colour = `Phospho (STY)`, label = Count)) +
    geom_bar(stat = "identity", position = "dodge", alpha = 0.8) +
    geom_text(
        stat = "identity", position = position_dodge(width = 0.9),
        vjust = -0.5) +
    ggpubr::theme_pubr() +
    scale_fill_manual(values = my_cols) +
    scale_colour_manual(values = my_cols)

novel_evid_stats_hq <- my_evid_filt_long %>%
    dplyr::filter(., `Phospho (STY) site IDs` %in% my_novel_stats_hq$id) %>%
    dplyr::group_by(., `Phospho (STY)`) %>%
    dplyr::summarise(., Count = dplyr::n_distinct(`Modified sequence`)) %>%
    dplyr::ungroup(.)

my_plots[["Phospho_multi_novel_hq"]] <- ggplot(
    novel_evid_stats_hq, aes(
        x = `Phospho (STY)`, y = Count, fill = `Phospho (STY)`,
        colour = `Phospho (STY)`, label = Count)) +
    geom_bar(stat = "identity", position = "dodge", alpha = 0.8) +
    geom_text(
        stat = "identity", position = position_dodge(width = 0.9),
        vjust = -0.5) +
    ggpubr::theme_pubr() +
    scale_fill_manual(values = my_cols) +
    scale_colour_manual(values = my_cols)



### Phosphorylated proteins ----------------------------------------------

my_pg_phospho <- my_data_filt %>%
    tidyr::separate_rows(data = ., Proteins, sep = ";") %>%
    dplyr::left_join(
        x = .,
        y = my_reason[, c("Proteins", "ORFNoveltyReason", "High quality")],
        by = "Proteins")

my_pg_phospho_count <- my_pg_phospho %>%
    dplyr::group_by(., Proteins) %>%
    dplyr::summarise(., Count = dplyr::n()) %>%
    dplyr::ungroup(.) %>%
    dplyr::arrange(., dplyr::desc(Count))

my_phospho_intens <- my_pg_phospho %>%
    dplyr::select(
        ., Proteins, `Positions.within.proteins`, `Amino.acid`,
        ORFNoveltyReason, `High quality`, `Localization.prob`,
        tidyselect::matches("Intensity.(Resuscitation|SCy004|SCy015)"))

my_phospho_identi_ref <- my_phospho_intens %>%
    dplyr::filter(., !grepl("^(c|p)", Proteins)) %>%
    tidyr::pivot_longer(
        data = .,
        cols = tidyselect::matches("Intensity.(Resuscitation|SCy004|SCy015)")) %>%
    dplyr::filter(., value > 0) %>%
    dplyr::mutate(., value = TRUE) %>%
    tidyr::pivot_wider(data = ., names_from = name, values_from = value) %>%
    #dplyr::left_join(x = complete_proteome, y = .) %>%
    dplyr::mutate(
        ., `Intensity.Resuscitation_chlorosis` = tidyr::replace_na(data = `Intensity.Resuscitation_chlorosis`, replace = FALSE),
        `Intensity.SCy004` = tidyr::replace_na(data = `Intensity.SCy004`, replace = FALSE),
        `Intensity.SCy015` = tidyr::replace_na(data = `Intensity.SCy015`, replace = FALSE))

my_plots[["venn_identification_ref"]] <- ggvenn::ggvenn(
    data = my_phospho_identi_ref,
    columns = c("Intensity.Resuscitation_chlorosis", "Intensity.SCy004", "Intensity.SCy015"),
    fill_color = c("#387eb8", "#d1d2d4", "#e21e25"))

my_phospho_intens_ref <- my_phospho_intens %>%
    dplyr::filter(., !grepl("^(c|p)", Proteins)) %>%
    tidyr::pivot_longer(
        data = .,
        cols = tidyselect::matches("Intensity.(Resuscitation|SCy004|SCy015)")) %>%
    dplyr::filter(., value > 0) %>%
    dplyr::group_by(., Proteins, `Positions.within.proteins`, `Amino.acid`) %>%
    dplyr::mutate(
        ., Resuscitation_chlorosis = dplyr::case_when(
            !any(name == "Intensity.Resuscitation_chlorosis") ~ "Not quanti.",
            all(name == "Intensity.Resuscitation_chlorosis") ~ "Unique quanti.",
            TRUE ~ "Shared id."
        ),
        SCy004 = dplyr::case_when(
            !any(name == "Intensity.SCy004") ~ "Not quanti.",
            all(name == "Intensity.SCy004") ~ "Unique quanti.",
            TRUE ~ "Shared id."
        ),
        SCy015 = dplyr::case_when(
            !any(name == "Intensity.SCy015") ~ "Not quanti.",
            all(name == "Intensity.SCy015") ~ "Unique quanti.",
            TRUE ~ "Shared id."
        ))

my_pg_f <- paste0(wkdir, "MQ_6frame_valid/combined/txt/proteinGroups.txt")

my_pg <- data.table::fread(
    input = my_pg_f, sep = "\t", quote = "", header = TRUE,
    stringsAsFactors = FALSE, colClasses = "character")

my_pg_format <- my_pg %>%
    dplyr::filter(., Reverse != "+" & `Potential contaminant` != "+") %>%
    dplyr::select(
        ., `Protein IDs`,
        tidyselect::matches("Intensity . (Resuscitation|SCy004|SCy015)")) %>%
    tidyr::separate_rows(data = ., `Protein IDs`, sep = ";")

my_pg_ibaq <- my_pg_format %>%
    tidyr::pivot_longer(data = ., cols = -`Protein IDs`) %>%
    tidyr::separate(
        data = ., col = name,
        into = c("Type", "Label", "Experiment"),
        sep = " ") %>%
    dplyr::group_by(., `Protein IDs`, Experiment) %>%
    dplyr::summarise(., Intensity = sum(as.double(value), na.rm = TRUE)) %>%
    dplyr::ungroup(.)

my_peptidome_f <- paste0(wkdir, "Genome/Tryptic_digest_peptidome.txt")

my_peptidome <- data.table::fread(
    input = my_peptidome_f, sep = "\t", quote = "", header = TRUE,
    stringsAsFactors = FALSE)

my_pg_ibaq %<>%
    dplyr::rename(., Proteins = `Protein IDs`) %>%
    dplyr::left_join(
        x = ., y = my_peptidome) %>%
    dplyr::mutate(., iBAQ = Intensity/Tryptic_count)

my_phospho_pg_intens_ref <- my_phospho_intens_ref %>%
    dplyr::select(., -name, -value) %>%
    unique(.) %>%
    tidyr::pivot_longer(
        data = ., cols = c(Resuscitation_chlorosis, SCy004, SCy015),
        names_to = "Experiment", values_to = "Group") %>%
    dplyr::left_join(x = ., y = my_pg_ibaq)

my_phospho_pg_intens_ref$Group <- factor(
    x = my_phospho_pg_intens_ref$Group,
    levels = unique(my_phospho_pg_intens_ref$Group))

my_plots[["Phospho_pg_density"]] <- ggplot(
    my_phospho_pg_intens_ref, aes(
        x = iBAQ, fill = Group, colour = Group)) +
    geom_density(alpha = 0.3) +
    #geom_density(aes(y = ..scaled..), alpha = 0.3) +
    ggpubr::theme_pubr() +
    facet_grid(rows = vars(Experiment)) +
    scale_x_log10() +
    scale_fill_manual(
        values = my_cols[1:length(unique(my_phospho_pg_intens_ref$Group))]) +
    scale_colour_manual(
        values = my_cols[1:length(unique(my_phospho_pg_intens_ref$Group))])

my_plots[["Phospho_pg_hist"]] <- ggplot(
    my_phospho_pg_intens_ref, aes(
        x = iBAQ, fill = Group, colour = Group)) +
    geom_histogram(position = "identity", alpha = 0.3) +
    ggpubr::theme_pubr() +
    facet_grid(rows = vars(Experiment)) +
    scale_x_log10() +
    scale_fill_manual(
        values = my_cols[1:length(unique(my_phospho_pg_intens_ref$Group))]) +
    scale_colour_manual(
        values = my_cols[1:length(unique(my_phospho_pg_intens_ref$Group))])

my_plots[["Phospho_pg_box"]] <- ggplot(
    my_phospho_pg_intens_ref, aes(
        x = Experiment, y = iBAQ, fill = Group, colour = Group)) +
    geom_boxplot(alpha = 0.6) +
    ggpubr::theme_pubr() +
    scale_y_log10() +
    scale_fill_manual(
        values = my_cols[1:length(unique(my_phospho_pg_intens_ref$Group))]) +
    scale_colour_manual(
        values = my_cols[1:length(unique(my_phospho_pg_intens_ref$Group))])



### Compiled data --------------------------------------------------------

my_motifs <- my_data_filt %>%
    dplyr::select(
        ., id, Proteins, `Amino.acid`, `Positions.within.proteins`,
        `Sequence.window`, `Localization.prob`) %>%
    dplyr::left_join(
        x = .,
        y = my_reason[, c("Proteins", "ORFNoveltyReason", "High quality")],
        by = "Proteins") %>%
    dplyr::mutate(
        ., Novel = ifelse(id %in% my_novel_stats$id, "+", ""),
        Novel_hq = ifelse(
            id %in% my_novel_stats$id & `High quality` == TRUE, "+", ""))

data.table::fwrite(
    x = my_motifs, file = "Phosphosites_categories.txt", append = FALSE,
    quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)



### Functional representation --------------------------------------------

my_annot_f <- paste0(wkdir, "Custom_annotation/2021-12-21_Custom_Uniprot_Eggnog_annotations.txt")

my_annot <- data.table::fread(
    input = my_annot_f, sep = "\t", quote = "", header = TRUE,
    stringsAsFactors = FALSE)

category <- unique(unlist(strsplit(my_annot$best_og_Subcategory, split = ";")))

my_annot_format <- my_annot %>%
    dplyr::mutate(., Simple_Subcategory = dplyr::case_when(
        grepl(";", best_og_Subcategory) ~ "Multiple categories",
        is.na(best_og_Subcategory) | best_og_Subcategory == "" ~ "Function unknown",
        best_og_Subcategory %in% category ~ best_og_Subcategory,
        TRUE ~ "Other categories"
    ))

my_motifs_annot <- my_motifs %>%
    tidyr::separate_rows(data = ., Proteins, sep = ";") %>%
    dplyr::left_join(
        x = .,
        y = my_annot_format[, c("#query_name", "Simple_Subcategory")],
        by = c("Proteins" = "#query_name")) %>%
    dplyr::group_by(., id) %>%
    dplyr::summarise_all(~paste0(unique(.), collapse = ";")) %>%
    dplyr::ungroup(.)

my_motifs_annot_count <- my_motifs_annot %>%
    dplyr::group_by(., Simple_Subcategory) %>%
    dplyr::summarise(
        ., Phospho_count = dplyr::n_distinct(id),
        Protein_count = dplyr::n_distinct(Proteins)) %>%
    dplyr::ungroup(.) %>%
    dplyr::mutate(., Ratio = Phospho_count/Protein_count) %>%
    dplyr::arrange(., Phospho_count)

my_motifs_annot_count$Simple_Subcategory <- factor(
    x = my_motifs_annot_count$Simple_Subcategory,
    levels = unique(my_motifs_annot_count$Simple_Subcategory), ordered = TRUE)

my_motifs_annot_count %<>%
    tidyr::pivot_longer(data = ., cols = c(Phospho_count, Protein_count))

my_plots[["Phospho_func_bar"]] <- ggplot(
    my_motifs_annot_count,
    aes(
        x = value, y = Simple_Subcategory,
        fill = name, colour = name, label = round(x = Ratio, digits = 2))) +
    geom_bar(
        stat = "identity",
        position = position_dodge(width = -0.3), alpha = 0.7) +
    geom_text(stat = "identity", position = "identity", colour = "black") +
    ggpubr::theme_pubr() +
    scale_fill_manual(
        values = my_cols[c(1, 5)]) +
    scale_colour_manual(
        values = my_cols[c(1, 5)])

my_data_top <- my_data %>%
    tidyr::separate_rows(data = ., Proteins, sep = ";") %>%
    dplyr::filter(., !grepl("^(c|p)", Proteins)) %>%
    dplyr::group_by(., Proteins) %>%
    dplyr::summarise(
        ., Count = dplyr::n_distinct(id),
        Count_loc = sum(Localization.prob > 0.95)) %>%
    dplyr::ungroup(.) %>%
    dplyr::arrange(., dplyr::desc(Count), dplyr::desc(Count_loc))

top_prot_phos <- my_data_top$Proteins[1:5]

my_data_top$Proteins <- factor(
    x = my_data_top$Proteins,
    levels = unique(my_data_top$Proteins),
    ordered = TRUE)

my_data_top %<>%
    tidyr::pivot_longer(data = ., cols = c(-Proteins))

my_plots[["Phospho_pg_top"]] <- ggplot(
    my_data_top %>% dplyr::filter(., Proteins %in% top_prot_phos),
    aes(
        x = Proteins, y = value,
        fill = name, colour = name, label = value)) +
    geom_bar(
        stat = "identity",
        position = "dodge", alpha = 0.7) +
    geom_text(
        stat = "identity",
        position = position_dodge(width = 0.9), vjust = -0.3) +
    ggpubr::theme_pubr() +
    scale_fill_manual(
        values = my_cols[c(1, 5)]) +
    scale_colour_manual(
        values = my_cols[c(1, 5)])

my_data_loc <- my_data %>%
    tidyr::separate_rows(
        data = ., Proteins, `Positions.within.proteins`, sep = ";") %>%
    dplyr::filter(., !grepl("^(c|p)", Proteins)) %>%
    dplyr::select(
        ., id, Proteins, `Positions.within.proteins`,
        `Amino.acid`, `Score.for.localization`, `Localization.prob`) %>%
    dplyr::mutate(
        ., `Score.for.localization` = as.integer(`Score.for.localization`)) %>%
    dplyr::arrange(., dplyr::desc(`Score.for.localization`)) %>%
    dplyr::mutate(
        ., Rank = 1:dplyr::n(),
        Localised = ifelse(
            `Localization.prob` > 0.95, "Localised", "Not localised"),
        Top_prot = ifelse(Proteins %in% top_prot_phos, "Top", "Others"),
        Label = ifelse(Proteins %in% top_prot_phos, Proteins, "Others"))

my_plots[["Phospho_localisation_rank"]] <- ggplot(
    my_data_loc, aes(
        x = Rank, y = `Score.for.localization`,
        colour = Localised, fill = Localised,
        shape = Label,
        size = Top_prot)) +
    geom_point() +
    ggpubr::theme_pubr() +
    scale_fill_manual(
        values = my_cols[c(1, 5)]) +
    scale_colour_manual(
        values = my_cols[c(1, 5)]) +
    scale_size_manual(values = c(2, 3)) +
    scale_shape_manual(values = c(
        Others = 1,
        sll1578 = 21,
        sll1031 = 22,
        sll0103 = 23,
        slr1841 = 24,
        slr0599 = 25
    ))



### Interesting phosphorylated target identification ---------------------

tmp <- my_annot_format %>%
    dplyr::filter(
        ., Simple_Subcategory == "Transcription")

my_target_phos <- my_data_top %>%
    tidyr::pivot_wider(data = ., names_from = name, values_from = value)

my_target_phos_exp <- my_phospho_intens_ref %>%
    dplyr::select(., -name, -value) %>%
    unique(.) %>%
    dplyr::group_by(., Proteins) %>%
    dplyr::summarise(
        ., Resuscitation_chlorosis_count = sum(
            Resuscitation_chlorosis %in% c("Unique quanti.", "Shared id.")),
        Resuscitation_chlorosis_uniq_count = sum(
            Resuscitation_chlorosis %in% c("Unique quanti.")),
        SCy004_count = sum(
            SCy004 %in% c("Unique quanti.", "Shared id.")),
        SCy004_uniq_count = sum(
            SCy004 %in% c("Unique quanti.")),
        SCy015_count = sum(
            SCy015 %in% c("Unique quanti.", "Shared id.")),
        SCy015_uniq_count = sum(
            SCy015 %in% c("Unique quanti.")))

my_target_phos %<>%
    dplyr::left_join(x = ., y = my_target_phos_exp) %>%
    dplyr::left_join(x = ., y = my_annot_format, by = c("Proteins" = "#query_name")) %>%
    dplyr::arrange(., dplyr::desc(Count))

my_uniq_exp_phos <- my_target_phos %>%
    dplyr::filter(
        ., Count >= 2 & (Count == Resuscitation_chlorosis_uniq_count |
            Count == SCy004_uniq_count | Count == SCy015_uniq_count))

data.table::fwrite(
    x = my_target_phos, file = "Phospho_per_conditions.txt",
    append = FALSE, quote = FALSE, sep = "\t",
    row.names = FALSE, col.names = TRUE)



### Export fasta files for phosphorylated proteins -----------------------

my_fasta_f <- paste0(wkdir, "Genome/Synechocystis_sp_PCC_6803_cds_aa.fasta")

my_fasta <- seqinr::read.fasta(
    file = my_fasta_f, seqtype = "AA",
    as.string = TRUE,
    whole.header = FALSE)

my_phospho_prot_list <- split(
    x = my_target_phos$Proteins,
    f = floor(1:length(my_target_phos$Proteins)/20))

for (x in names(my_phospho_prot_list)) {
    
    my_fasta_filt <- my_fasta[names(my_fasta) %in% my_phospho_prot_list[[x]]]
    headers <- lapply(my_fasta_filt, function(y) {
        attr(y, "Annot")}) %>%
        unlist(.) %>%
        sub("^>", "", .)
    seqinr::write.fasta(
        sequences = my_fasta_filt, names = headers,
        file.out = paste0("Phosphorylated_proteins_20_", x, ".fasta"),
        open = "w", nbchar = 60, as.string = TRUE)
    
}



### Location of phospho in protein structure -----------------------------

my_netsurfp_f <- paste0(wkdir, "2022-02-23_NetsurfP/Phosphorylated_proteins_formatted.txt")

my_netsurfp <- data.table::fread(
    input = my_netsurfp_f, sep = "\t",
    header = TRUE, stringsAsFactors = FALSE)

my_netsurfp_stats <- my_data_loc %>%
    dplyr::mutate(
        ., `Positions.within.proteins` = as.integer(Positions.within.proteins)) %>%
    dplyr::left_join(
        x = my_netsurfp, y = ., by = c(
            "id" = "Proteins", "seq" = "Amino.acid",
            "n" = "Positions.within.proteins"))

my_big_cols <- c(
    "#387eb8", "#404040", "#e21e25", "#fbaf3f", "#d1d2d4",
    "#246E39", "#753B94", "#009DB5")

for (x in c("Class", "q3", "q8")) {
    
    my_netsurfp_stats_toplot <- my_netsurfp_stats %>%
        dplyr::filter(., seq %in% c("S", "T", "Y")) %>%
        dplyr::group_by(., !!as.name(x)) %>%
        dplyr::summarise(
            ., Total = dplyr::n(),
            Phospho = sum(!is.na(Localised)),
            Phospho_loc = sum(!is.na(Localised) & Localised == "Localised")) %>%
        tidyr::pivot_longer(
            data = ., cols = -!!as.name(x), names_to = "Type", values_to = "Count")
    
    my_netsurfp_stats_toplot$Type <- factor(
        x = my_netsurfp_stats_toplot$Type,
        levels = unique(my_netsurfp_stats_toplot$Type),
        ordered = TRUE)
    
    my_plots[[paste0("netsufrp_dodge_", x)]] <- ggplot(
        my_netsurfp_stats_toplot, aes(
            x = !!as.name(x), y = Count, fill = Type,
            colour = Type, label = Count)) +
        geom_bar(stat = "identity", position = "dodge", alpha = 0.8) +
        geom_text(
            stat = "identity", position = position_dodge(width = 0.9),
            vjust = -0.5) +
        ggpubr::theme_pubr() +
        scale_fill_manual(values = my_big_cols) +
        scale_colour_manual(values = my_big_cols)
    
    my_plots[[paste0("netsufrp_stack_", x)]] <- ggplot(
        my_netsurfp_stats_toplot, aes(
            x = Type, y = Count, fill = !!as.name(x),
            colour = !!as.name(x), label = Count)) +
        geom_bar(stat = "identity", position = "fill", alpha = 0.8) +
        geom_text(
            stat = "identity", position = "fill",
            vjust = -0.5) +
        ggpubr::theme_pubr() +
        scale_fill_manual(values = my_big_cols) +
        scale_colour_manual(values = my_big_cols)
    
}

for (x in c("rsa", "AlphaProb", "BetaProb", "CoilProb", "disorder")) {
    
    my_netsurfp_stats_toplot <- my_netsurfp_stats %>%
        dplyr::filter(., seq %in% c("S", "T", "Y")) %>%
        dplyr::mutate(., Type = dplyr::case_when(
            is.na(Localised) ~ "Not phosphorylated",
            TRUE ~ Localised
        ))
    
    my_netsurfp_stats_toplot$Type <- factor(
        x = my_netsurfp_stats_toplot$Type,
        levels = unique(my_netsurfp_stats_toplot$Type),
        ordered = TRUE)
    
    my_plots[[paste0("netsufrp_violin_", x)]] <- ggplot(
        my_netsurfp_stats_toplot, aes(
            x = Type, y = !!as.name(x), fill = Type,
            colour = Type)) +
        geom_violin(alpha = 0.4) +
        geom_boxplot(width = 0.2, alpha = 0.7, size = 1.2) +
        ggpubr::theme_pubr() +
        scale_fill_manual(values = my_big_cols) +
        scale_colour_manual(values = my_big_cols)
    
}



### Phosphorylated protein coverage --------------------------------------

pl_list <- list()

prot_to_plot <- my_netsurfp %>%
    dplyr::mutate(., start = n, end = n, y = 0, Sequence = seq, Structure = q3)

phos_to_plot <- my_data %>%
    tidyr::separate_rows(
        data = ., Proteins, `Positions.within.proteins`, sep = ";") %>%
    dplyr::filter(., !grepl("^(c|p)", Proteins)) %>%
    dplyr::select(
        ., id, Proteins, `Positions.within.proteins`,
        `Amino.acid`, `Score.for.localization`, `Localization.prob`,
        tidyselect::matches("Intensity.(Resuscitation|SCy004|SCy015)")) %>%
    dplyr::mutate(
        ., `Score.for.localization` = as.integer(`Score.for.localization`)) %>%
    tidyr::pivot_longer(
        data = .,
        cols = tidyselect::matches("Intensity.(Resuscitation|SCy004|SCy015)"),
        names_to = "samples", values_to = "Intensity") %>%
    dplyr::group_by(
        ., id, Proteins, `Positions.within.proteins`, `Amino.acid`,
        `Score.for.localization`, `Localization.prob`) %>%
    dplyr::mutate(
        ., samples = ifelse(
            as.double(Intensity) > 0,
            sub("Intensity\\.", "",
                samples), NA)) %>%
    dplyr::summarise(., Quantified = paste0(na.omit(samples), collapse = ";")) %>%
    dplyr::ungroup(.) %>%
    dplyr::mutate(., Quantified = dplyr::case_when(
        is.na(Quantified) | Quantified == "" ~ "Not quantified",
        grepl(";", Quantified) ~ "Shared",
        TRUE ~ Quantified),
        Localised = ifelse(
            `Localization.prob` > 0.95, "Localised", "Not localised"))

for (x in unique(prot_to_plot$id)) {
    
    prot_to_plot_tmp <- prot_to_plot %>%
        dplyr::filter(., id == x)
    
    phos_to_plot_tmp <- phos_to_plot %>%
        dplyr::filter(., Proteins == x)
    
    pl_list[[paste0("Phospho_", x)]] <- ggplot() +
        geom_rect(
            data = prot_to_plot_tmp,
            mapping = aes(
                xmin = start-0.5,
                xmax = end+0.5,
                ymin = y,
                ymax = (y + -4),
                fill = Structure),
            colour = "black") +
        geom_segment(
            data = phos_to_plot_tmp,
            mapping = aes(
                x = as.integer(`Positions.within.proteins`),
                xend = as.integer(`Positions.within.proteins`),
                y = 0,
                yend = as.numeric(`Score.for.localization`),
                linetype = Localised),
            colour = "#3B3B3B") +
        geom_point(
            data = phos_to_plot_tmp,
            mapping = aes(
                x = as.integer(`Positions.within.proteins`),
                y = as.numeric(`Score.for.localization`),
                fill = Quantified),
            colour = "#3B3B3B",
            size = 7, shape = 21) +
        geom_text(
            data = phos_to_plot_tmp,
            mapping = aes(
                x = as.integer(`Positions.within.proteins`),
                y = as.numeric(`Score.for.localization`),
                label = `Amino.acid`),
            colour = "white",
            size = 5) +
        ggpubr::theme_pubr() +
        xlab("Amino acid position") +
        ylab("Localisation probability") +
        scale_fill_manual(
            values = c(
                `Resuscitation_chlorosis` = "#387eb8", `SCy004` = "#e21e25",
                `SCy015` = "#fbaf3f", `Not quantified` = "#d1d2d4",
                `Shared` = "#3B3B3B", `C` = "white", `E` = "#d1d2d4",
                `H` = "#3B3B3B")) +
        scale_linetype_manual(
            values = c(`Localised` = "solid", `Not localised` = "dotted")) +
        ggtitle(x)
    
}



### Export files ---------------------------------------------------------

pdf("Phosphorylation_statistics.pdf", 10, 10)
my_plots
dev.off()

pdf("Phosphorylation_coverage.pdf", 10, 10)
pl_list
dev.off()

save.image("Phosphorylation_statistics.RData")


