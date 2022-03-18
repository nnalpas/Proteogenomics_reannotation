


rm(list = ls())

library(magrittr)
library(ggplot2)

my_plots <- list()

my_ibaq_f <- "/mnt/storage/kxmna01/data/Synechocystis_6frame/2022-02-14_iBAQ/Scy004_resuscitation_Scy001_iBAQ_norm.txt"

my_ibaq <- data.table::fread(
    input = my_ibaq_f, sep = "\t", quote = "", header = TRUE,
    stringsAsFactors = FALSE, colClasses = "character", data.table = FALSE)

my_ibaq_filt <- my_ibaq %>%
    dplyr::select(., Proteins, tidyselect::starts_with("SCy001_L")) %>%
    tidyr::pivot_longer(data = ., cols = -Proteins) %>%
    dplyr::group_by(., Proteins) %>%
    dplyr::summarise(., Mean = mean(as.double(value), na.rm = TRUE)) %>%
    dplyr::ungroup(.) %>%
    dplyr::mutate(., Mean = ifelse(!is.na(Mean) & Mean > 0, Mean, NA))

my_path <- "/mnt/storage/kxmna01/data/Synechocystis_6frame/Codon_usage_other_species/iBAQ"

my_files <- list.files(
    path = my_path, pattern = "*iBAQ.txt", full.names = TRUE, recursive = TRUE) %>%
    set_names(sub("_protein_iBAQ.txt", "", basename(.)))

my_exp <- c(
    `Arabidopsis_thaliana` = "^Intensity",
    `Bacillus_subtilis` = "^Intensity",
    `Escherichia_coli` = "^Intensity",
    `Euglena_gracilis` = "^Intensity (1|2|3)",
    `Phaeodactylum_tricornutum` = "^Intensity (1|2|3|4|5|6)",
    `Saccharomyces_cerevisiae` = "^Intensity J11_(A|B|C|D)",
    `Synechococcus_elongatus` = "^Intensity"
)

my_other_ibaq <- lapply(my_files, function(x) {
    data.table::fread(
        input = x, sep = "\t", quote = "", header = TRUE,
        stringsAsFactors = FALSE, colClasses = "character", data.table = FALSE)
})

my_other_ibaq_filt <- lapply(names(my_other_ibaq), function(x) {
    my_other_ibaq[[x]] %>%
        dplyr::select(
            ., Proteins = `Protein IDs`, tidyselect::matches(my_exp[[x]])) %>%
        tidyr::pivot_longer(data = ., cols = -Proteins) %>%
        dplyr::group_by(., Proteins) %>%
        dplyr::summarise(., Mean = mean(as.double(value), na.rm = TRUE)) %>%
        dplyr::ungroup(.) %>%
        dplyr::mutate(., Mean = ifelse(!is.na(Mean) & Mean > 0, Mean, NA)) %>%
        tidyr::separate_rows(data = ., Proteins, sep = ";")
}) %>%
    set_names(names(my_other_ibaq)) %>%
    plyr::ldply(., data.table::data.table, .id = "Species")

my_eggnog_f <- "/mnt/storage/kxmna01/data/Synechocystis_6frame/Custom_annotation/2021-12-21_Custom_Uniprot_Eggnog_annotations.txt"

my_eggnog <- data.table::fread(
    input = my_eggnog_f, sep = "\t", quote = "", header = TRUE,
    stringsAsFactors = FALSE, colClasses = "character", data.table = FALSE) %>%
    dplyr::select(., Proteins = `#query_name`, KEGG_ko, best_og_Subcategory) %>%
    tidyr::separate_rows(data = ., KEGG_ko, sep = ";") %>%
    unique(.)

my_en_path <- "/mnt/storage/kxmna01/data/Synechocystis_6frame/Codon_usage_other_species/eggnog"

my_en_files <- list.files(
    path = my_en_path, pattern = "*.annotations", full.names = TRUE, recursive = TRUE) %>%
    set_names(sub("_protein.emapper.annotations", "", basename(.)))

my_other_eggnog <- lapply(my_en_files, function(x) {
    data.table::fread(
        input = x, sep = "\t", quote = "", header = TRUE,
        stringsAsFactors = FALSE, colClasses = "character", data.table = FALSE) %>%
        dplyr::select(., Proteins = `#query_name`, KEGG_ko, best_og_cat)
})

my_other_eggnog$Arabidopsis_thaliana$Proteins %<>% sub("\\..*", "", .)
my_other_eggnog$Bacillus_subtilis$Proteins %<>% sub("^.+\\|(.+)\\|.+", "\\1", .)
my_other_eggnog$Escherichia_coli$Proteins %<>% sub("^.+\\|(.+)\\|.+", "\\1", .)
my_other_eggnog$Synechococcus_elongatus$Proteins %<>% sub("^.+_prot_", "", .) %>% sub("(\\..+)_.+$", "\\1", .)

my_fasta <- seqinr::read.fasta(
    file = "/mnt/storage/kxmna01/data/Synechocystis_6frame/Codon_usage_other_species/Remmers-2018/Phaeodactylum_tricornutum_protein.fasta",
    seqtype = "AA", as.string = TRUE, whole.header = TRUE)
my_cross <- data.frame(Header = names(my_fasta)) %>%
    tidyr::separate(
        data = ., col = Header, into = c("ID", "Proteins"),
        sep = " ", extra = "merge") %>%
    dplyr::mutate(
        ., Proteins = sub(".*\\[locus_tag=(.+?)\\].*", "\\1", Proteins))
my_other_eggnog[["Phaeodactylum_tricornutum"]] %<>%
    dplyr::right_join(x = my_cross, y = ., by = c("ID" = "Proteins")) %>%
    dplyr::select(., -ID) %>%
    dplyr::filter(., !is.na(Proteins))

my_other_eggnog_format <- my_other_eggnog %>%
    plyr::ldply(., data.table::data.table, .id = "Species") %>%
    tidyr::separate_rows(data = ., KEGG_ko, sep = ",") %>%
    unique(.)

my_ibaq_filt_eg <- dplyr::left_join(
    x = my_ibaq_filt,
    y = my_eggnog %>% dplyr::select(., -best_og_Subcategory),
    by = "Proteins")

my_other_ibaq_filt_eg <- dplyr::left_join(
    x = my_other_ibaq_filt,
    y = my_other_eggnog_format %>% dplyr::select(., -best_og_cat),
    by = c("Species", "Proteins"))

my_eg_quanti <- my_ibaq_filt_eg %>%
    dplyr::mutate(., Species = "Synechocystis sp.") %>%
    dplyr::bind_rows(., my_other_ibaq_filt_eg) %>%
    dplyr::filter(., !is.na(Mean) & Mean > 0) %>%
    dplyr::mutate(., KEGG_ko = dplyr::case_when(
        is.na(KEGG_ko) ~ NA_character_,
        KEGG_ko %in% c("", "-") ~ NA_character_,
        TRUE ~ KEGG_ko
    ))

my_eg_quanti_format <- my_eg_quanti %>%
    dplyr::mutate(., KEGG_ko = ifelse(is.na(KEGG_ko), Proteins, KEGG_ko)) %>%
    dplyr::group_by(., KEGG_ko, Species) %>%
    dplyr::summarise(
        ., Count = dplyr::n_distinct(Proteins),
        Mean = mean(Mean, na.rm = TRUE)) %>%
    dplyr::ungroup(.)
my_eg_quanti_format$Log_mean <- log10(my_eg_quanti_format$Mean)

my_eg_quanti_corr <- my_eg_quanti_format %>%
    dplyr::select(., -Count, -Mean) %>%
    tidyr::pivot_wider(data = ., names_from = Species, values_from = Log_mean)

my_stats <- my_eg_quanti %>%
    dplyr::group_by(., Species) %>%
    dplyr::summarise(
        ., Quantified = dplyr::n_distinct(Proteins),
        Ortholog = dplyr::n_distinct(na.omit(KEGG_ko)))

my_eg_individual <- my_eg_quanti_corr %>%
    dplyr::filter(., !is.na(`Synechocystis sp.`)) %>%
    tidyr::pivot_longer(
        data = ., cols = -KEGG_ko,
        names_to = "Species", values_to = "Log_mean") %>%
    dplyr::mutate(., Photosynthetic = dplyr::case_when(
        Species %in% c(
            "Arabidopsis_thaliana", "Euglena_gracilis",
            "Phaeodactylum_tricornutum", "Synechococcus_elongatus",
            "Synechocystis sp.") ~ 1,
        TRUE ~ 0
    )) %>%
    dplyr::filter(., !is.na(Log_mean)) %>%
    dplyr::group_by(., KEGG_ko) %>%
    dplyr::summarise(
        ., All = TRUE,
        Ortholog = unique(grepl("^ko", KEGG_ko)),
        Count = dplyr::n_distinct(Species),
        Count_photo = sum(Photosynthetic),
        Species = paste0(unique(Species), collapse = ";")) %>%
    dplyr::ungroup(.) %>%
    dplyr::mutate(
        ., Shared = Count > 1,
        Core = (Count == 8),
        Photosynthesis = (Count_photo == 5))

my_plots[["venn"]] <- ggvenn::ggvenn(
    data = my_eg_individual,
    columns = c("Ortholog", "Shared", "Core", "Photosynthesis"),
    fill_color = c("#387eb8", "#d1d2d4", "#e21e25", "#fbaf3f"))

my_residuals <- data.frame()

for (i in unique(my_stats$Species)) {
    
    x <- my_eg_quanti_corr[["Synechocystis sp."]]
    y <- my_eg_quanti_corr[[i]]
    ok <- complete.cases(x, y)
    corr <- cor.test(
        x = x, y = y, exact = TRUE,
        alternative = "two.sided", method = "spearman")
    my_stats[my_stats$Species == i, "Shared"] <- sum(ok)
    my_stats[my_stats$Species == i, "Rho"] <- corr$estimate
    my_stats[my_stats$Species == i, "Pvalue"] <- corr$p.value
    
    toplot <- my_eg_quanti_corr %>%
        dplyr::select(., KEGG_ko, `Synechocystis sp.`, !!as.name(i)) %>%
        dplyr::filter(., !is.na(`Synechocystis sp.`) & !is.na(!!as.name(i)))
    
    toplot$Residuals <- 0
    toplot$Category <- "q10-q90"
    if (i != "Synechocystis sp.") {
        my_fit <- lm(
            formula = as.formula(paste(i, "~ `Synechocystis sp.`")),
            data = toplot, x = TRUE, y = TRUE)
        toplot$Residuals <- my_fit$residuals
        quants <- quantile(toplot$Residuals, c(0.1, 0.9))
        toplot %<>%
            dplyr::mutate(., Category = dplyr::case_when(
                Residuals < quants[["10%"]] ~ "q10",
                Residuals > quants[["90%"]] ~ "q90",
                TRUE ~ "q10-q90"
            ))
    }
    
    my_plots[[paste0("Scatter_", i)]] <- ggplot(toplot, aes(
        x = `Synechocystis sp.`, y = !!as.name(i), colour = Category)) +
        geom_point() +
        geom_smooth(method = "lm", formula = y~x, colour = "black") +
        ggpubr::theme_pubr() +
        scale_colour_manual(
            values = c(
                `q10` = "#387eb8", `q10-q90` = "#404040", `q90` = "#e21e25"))
    
    my_residuals <- toplot %>%
        dplyr::select(., KEGG_ko, Residuals, Category) %>%
        dplyr::mutate(., Species = i) %>%
        dplyr::bind_rows(my_residuals, .)
    
}

my_eg_quanti_final <- my_eg_quanti_format %>%
    dplyr::left_join(x = ., y = my_residuals, by = c("KEGG_ko", "Species"))

og <- data.table::fread(
    input = "/mnt/storage/kxmna01/Func_categories.txt",
    sep = "\t", quote = "", header = TRUE,
    stringsAsFactors = FALSE)

my_ogcat <- og %>%
    dplyr::rename(
        ., best_og_code = ID,
        best_og_Category = Category,
        best_og_Subcategory = Subcategory)

my_eg_annot <- my_other_eggnog_format %>%
    tidyr::separate_rows(data = ., best_og_cat, sep = "(?<=.)(?=.)") %>%
    dplyr::left_join(
        x = ., y = my_ogcat, by = c("best_og_cat" = "best_og_code")) %>%
    dplyr::select(., Proteins, KEGG_ko, best_og_Subcategory)

my_eg_annot <- my_eggnog %>%
    dplyr::select(., Proteins, KEGG_ko, best_og_Subcategory) %>%
    dplyr::bind_rows(my_eg_annot, .)

my_eg_annot_format <- my_eg_annot %>%
    tidyr::pivot_longer(
        data = ., cols = -best_og_Subcategory,
        names_to = "Type", values_to = "ID") %>%
    tidyr::separate_rows(data = ., best_og_Subcategory, sep = ";") %>%
    dplyr::mutate(., best_og_Subcategory = dplyr::case_when(
        is.na(best_og_Subcategory) | best_og_Subcategory %in% c("", "-") ~ "Function unknown",
        TRUE ~ best_og_Subcategory
    )) %>%
    dplyr::filter(., !is.na(ID) & !ID %in% c("", "-"))

my_eg_annot_final <- my_eg_annot_format %>%
    dplyr::group_by(., ID) %>%
    dplyr::summarise(
        ., best_og_Subcategory = paste0(
            unique(best_og_Subcategory), collapse = ";")) %>%
    dplyr::ungroup(.)

my_stats_cat_all <- data.frame()

for (p in unique(my_eg_annot_format$best_og_Subcategory)) {
    
    path_ids <- my_eg_annot_format %>%
        dplyr::filter(., best_og_Subcategory == p) %>%
        .[["ID"]] %>%
        unique(.)
    
    my_stats_cat <- my_eg_quanti %>%
        dplyr::filter(., Proteins %in% path_ids | KEGG_ko %in% path_ids) %>%
        dplyr::group_by(., Species) %>%
        dplyr::summarise(
            ., Quantified = dplyr::n_distinct(Proteins),
            Ortholog = dplyr::n_distinct(na.omit(KEGG_ko))) %>%
        dplyr::ungroup(.) %>%
        dplyr::mutate(., Category = p)
    
    for (i in unique(my_stats_cat$Species)) {
        
        x <- my_eg_quanti_corr[my_eg_quanti_corr$KEGG_ko %in% path_ids, ][["Synechocystis sp."]]
        y <- my_eg_quanti_corr[my_eg_quanti_corr$KEGG_ko %in% path_ids, ][[i]]
        ok <- complete.cases(x, y)
        if (sum(ok) > 2) {
            corr <- cor.test(
                x = x, y = y, exact = TRUE,
                alternative = "two.sided", method = "spearman")
        } else {
            corr <- list(estimate = NA, p.value = NA)
        }
        my_stats_cat[my_stats_cat$Species == i, "Shared"] <- sum(ok)
        my_stats_cat[my_stats_cat$Species == i, "Rho"] <- corr$estimate
        my_stats_cat[my_stats_cat$Species == i, "Pvalue"] <- corr$p.value
        
    }
    
    my_stats_cat_all %<>%
        dplyr::bind_rows(., my_stats_cat)
    
}

toplot <- my_stats_cat_all %>%
    dplyr::filter(., Species != "Synechocystis sp." & Pvalue < 0.05)
toplot$Category <- factor(
    x = toplot$Category,
    levels = sort(unique(toplot$Category), decreasing = TRUE),
    ordered = TRUE)

my_palette <- grDevices::hcl.colors(n = 11, palette = "Fall")

my_plots[["Heatmap"]] <- ggplot(
    toplot, aes(x = Species, y = Category, colour = Rho, size = Shared)) +
    geom_point() +
    ggpubr::theme_pubr() +
    scale_color_gradientn(colours = my_palette, limits = c(0, 1)) +
    scale_size(range = c(2, 20))

toplot <- my_stats %>%
    dplyr::select(., -Pvalue) %>%
    tidyr::pivot_longer(data = ., cols = -Species) %>%
    dplyr::filter(., Species != "Synechocystis sp.")
toplot$name <- factor(
    x = toplot$name, levels = c("Quantified", "Ortholog", "Shared", "Rho"),
    ordered = TRUE)

my_plots[["Comparison_stats"]] <- ggplot(
    toplot, aes(x = Species, y = value, fill = name, colour = name)) +
    geom_bar(stat = "identity", position = "dodge", alpha = 0.7) +
    ggpubr::theme_pubr() +
    facet_grid(rows = vars(name), scales = "free_y") +
    scale_colour_manual(
        values = c("#387eb8", "#404040", "#e21e25", "#fbaf3f")) +
    scale_fill_manual(
        values = c("#387eb8", "#404040", "#e21e25", "#fbaf3f"))

data.table::fwrite(
    x = my_eg_annot_final,
    file = "2022-03-18_Cross_species_Eggnog_annotations.txt",
    append = FALSE, quote = FALSE, sep = "\t",
    row.names = FALSE, col.names = TRUE)

data.table::fwrite(
    x = my_eg_individual,
    file = "2022-03-18_Synechocystis_ortholog_comparison.txt",
    append = FALSE, quote = FALSE, sep = "\t",
    row.names = FALSE, col.names = TRUE)

my_eg_individual_oa <- my_eg_individual %>%
    dplyr::select(., -Count, -Count_photo, -Species) %>%
    dplyr::rowwise(.) %>%
    dplyr::mutate(
        .,
        All_unique = dplyr::case_when(
            All & sum(c(All, Ortholog, Shared, Photosynthesis, Core)) == 1 ~ TRUE,
            TRUE ~ FALSE),
        Ortholog_unique = dplyr::case_when(
            Ortholog & sum(c(Ortholog, Shared, Photosynthesis, Core)) == 1 ~ TRUE,
            TRUE ~ FALSE),
        Shared_unique = dplyr::case_when(
            Shared & sum(c(Shared, Photosynthesis, Core)) == 1 ~ TRUE,
            TRUE ~ FALSE),
        Photosynthesis_unique = dplyr::case_when(
            Photosynthesis & sum(c(Photosynthesis, Core)) == 1 ~ TRUE,
            TRUE ~ FALSE),
        Core_unique = dplyr::case_when(
            Core & sum(c(Core)) == 1 ~ TRUE,
            TRUE ~ FALSE))

data.table::fwrite(
    x = my_eg_individual_oa,
    file = "2022-03-18_Synechocystis_ortholog_comparison_forOA.txt",
    append = FALSE, quote = FALSE, sep = "\t",
    row.names = FALSE, col.names = TRUE)

data.table::fwrite(
    x = my_eg_quanti_final,
    file = "2022-03-18_Synechocystis_comparison_all.txt",
    append = FALSE, quote = FALSE, sep = "\t",
    row.names = FALSE, col.names = TRUE)

save.image("2022-03-18_Synechocystis_comparison.RData")


