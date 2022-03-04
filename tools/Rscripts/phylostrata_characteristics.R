


rm(list = ls())

library(magrittr)
library(ggplot2)
library(ggplot2)

my_plots <- list()

my_nucl_f <- "/mnt/storage/kxmna01/data/Synechocystis_6frame/Genome/Synechocystis_sp_PCC_6803_gene_nucl.fasta"

my_nucl <- seqinr::read.fasta(
    file = my_nucl_f, seqtype = "DNA", as.string = TRUE, whole.header = FALSE)

my_annot_f <- "/mnt/storage/kxmna01/data/Synechocystis_6frame/Custom_annotation/2021-12-21_Custom_Uniprot_Eggnog_annotations.txt"

my_annot <- data.table::fread(
    input = my_annot_f, sep = "\t", quote = "",
    header = TRUE, stringsAsFactors = FALSE, colClasses = "character")

my_phylo_f <- "/mnt/storage/kxmna01/data/Synechocystis_6frame/Phylostratigraphy/Phylostrata_proteins.txt"

my_phylo <- data.table::fread(
    input = my_phylo_f, sep = "\t", quote = "",
    header = TRUE, stringsAsFactors = FALSE)

my_ibaq_f <- "/mnt/storage/kxmna01/data/Synechocystis_6frame/2022-02-14_iBAQ/Scy004_resuscitation_Scy001_iBAQ_raw.txt"

my_ibaq <- data.table::fread(
    input = my_ibaq_f, sep = "\t", quote = "",
    header = TRUE, stringsAsFactors = FALSE, colClasses = "character")

my_ibaq_format <- my_ibaq %>%
    tidyr::pivot_longer(data = ., cols = -Proteins) %>%
    dplyr::filter(., grepl("SCy001_L", name)) %>%
    dplyr::mutate(., Experiment = sub("^(.+?)_.+$", "\\1", name)) %>%
    dplyr::group_by(., Experiment, Proteins) %>%
    dplyr::summarise(., iBAQ = median(as.double(value), na.rm = TRUE))

my_phylo_format <- my_phylo %>%
    dplyr::left_join(x = ., y = my_ibaq_format, by = c("qseqid" = "Proteins")) %>%
    dplyr::filter(., !is.na(Experiment))

my_phylo_format$ps <- factor(
    x = my_phylo_format$ps, levels = c(1:10), ordered = TRUE)

my_gc <- lapply(my_nucl, function(x) {
    stringr::str_count(string = x, pattern = "g|c")/nchar(x)*100
}) %>%
    plyr::ldply(., data.table::data.table) %>%
    set_colnames(c("qseqid", "GC%"))

my_phylo_format_gc <- my_phylo_format %>%
    dplyr::full_join(x = ., y = my_gc)

my_nucl_set <- coRdon::readSet(file = my_nucl_f)
my_codon <- coRdon::codonTable(x = my_nucl_set)

my_codon_stats <- my_codon@counts %>%
    as.data.frame(.) %>%
    dplyr::mutate(
        ., qseqid = sub(" .+", "", my_codon@ID),
        length = my_codon@len) %>%
    tidyr::pivot_longer(data = ., cols = c(-qseqid, -length)) %>%
    dplyr::mutate(., `Position3` = dplyr::case_when(
        grepl("..(G|C)", name, ignore.case = TRUE) ~ "GC3",
        TRUE ~ "AT3")) %>%
    dplyr::group_by(., qseqid, `Position3`) %>%
    dplyr::summarise(
        ., Percentage = sum(value, na.rm = TRUE)/unique(length)*100) %>%
    dplyr::ungroup(.) %>%
    tidyr::pivot_wider(
        data = ., names_from = `Position3`, values_from = Percentage) %>%
    dplyr::full_join(x = my_phylo_format_gc, y = .)

my_codon_stats <- coRdon::MILC(cTobject = my_codon, id_or_name2 = "11") %>%
    as.data.frame(.) %>%
    set_colnames(c("milc")) %>%
    dplyr::mutate(., qseqid = sub(" .+", "", names(my_nucl_set))) %>%
    dplyr::full_join(x = my_codon_stats, y = .)

my_codon_stats <- coRdon::B(cTobject = my_codon, id_or_name2 = "11") %>%
    as.data.frame(.) %>%
    set_colnames(c("bias")) %>%
    dplyr::mutate(., qseqid = sub(" .+", "", names(my_nucl_set))) %>%
    dplyr::full_join(x = my_codon_stats, y = .)

my_codon_stats <- coRdon::ENC(cTobject = my_codon, id_or_name2 = "11") %>%
    as.data.frame(.) %>%
    set_colnames(c("enc")) %>%
    dplyr::mutate(., qseqid = sub(" .+", "", names(my_nucl_set))) %>%
    dplyr::full_join(x = my_codon_stats, y = .)

my_codon_stats <- coRdon::ENCprime(cTobject = my_codon, id_or_name2 = "11") %>%
    as.data.frame(.) %>%
    set_colnames(c("encprime")) %>%
    dplyr::mutate(., qseqid = sub(" .+", "", names(my_nucl_set))) %>%
    dplyr::full_join(x = my_codon_stats, y = .)

my_codon_stats <- coRdon::MCB(cTobject = my_codon, id_or_name2 = "11") %>%
    as.data.frame(.) %>%
    set_colnames(c("mcb")) %>%
    dplyr::mutate(., qseqid = sub(" .+", "", names(my_nucl_set))) %>%
    dplyr::full_join(x = my_codon_stats, y = .)

my_codon_stats <- coRdon::SCUO(cTobject = my_codon, id_or_name2 = "11") %>%
    as.data.frame(.) %>%
    set_colnames(c("scuo")) %>%
    dplyr::mutate(., qseqid = sub(" .+", "", names(my_nucl_set))) %>%
    dplyr::full_join(x = my_codon_stats, y = .)

my_nucl_char <- seqinr::read.fasta(
    file = my_nucl_f, seqtype = "DNA", as.string = FALSE, whole.header = FALSE)

my_high_exp <- my_ibaq_format %>%
    dplyr::filter(., !grepl("^(c|p)", Proteins)) %>%
    dplyr::mutate(., Total = sum(iBAQ, na.rm = TRUE), Perc = iBAQ/Total*100) %>%
    dplyr::arrange(., dplyr::desc(Perc)) %>%
    dplyr::mutate(., Cumsum = cumsum(Perc)) %>%
    dplyr::slice(., 1:70) %>%
    .[["Proteins"]]

my_nucl_char_high <- my_nucl_char[names(my_nucl_char) %in% my_high_exp]
my_w <- lapply(my_nucl_char_high, function(x) {
    seqinr::uco(
        seq = x, frame = 0, index = "eff")
}) %>%
    plyr::ldply(., data.table::data.table) %>%
    set_colnames(c("qseqid", "CODON", "Count")) %>%
    dplyr::full_join(x = ., y = seqinr::SEQINR.UTIL$CODON.AA) %>%
    dplyr::group_by(., AA, L, CODON) %>%
    dplyr::summarise(., Sum = sum(Count, na.rm = TRUE)) %>%
    dplyr::group_by(., AA, L) %>%
    dplyr::mutate(., w = Sum/max(Sum)) %>%
    dplyr::ungroup(.) %>%
    dplyr::arrange(., CODON)

my_cai <- lapply(my_nucl_char, function(x) {
    seqinr::cai(seq = x, numcode = 11, w = my_w$w)
}) %>%
    plyr::ldply(., data.table::data.table) %>%
    set_colnames(c("qseqid", "cai"))

my_uco <- lapply(my_nucl_char, function(x) {
    seqinr::uco(
        seq = x, frame = 0, as.data.frame = TRUE)
}) %>%
    plyr::ldply(., data.table::data.table, .id = "qseqid")

my_codon_stats %<>%
    dplyr::full_join(x = ., my_cai)

for (x in c(
    "iBAQ", "AT3", "GC3", "GC%", "milc", "bias",
    "enc", "encprime", "mcb", "scuo", "cai")) {
    
    my_plots[[paste0(x, "_violin")]] <- ggplot(
        my_codon_stats, aes(
            x = ps, y = !!as.name(x), group = ps, fill = ps, colour = ps)) +
        geom_violin(alpha = 0.3) +
        geom_boxplot(alpha = 0.7, width = 0.2) +
        scale_y_log10() +
        ggpubr::theme_pubr()
    
    my_plots[[paste0(x, "_scatter")]] <- ggplot(
        my_codon_stats, aes(x = as.integer(ps), y = !!as.name(x))) +
        geom_point() +
        scale_y_log10() +
        geom_smooth(method='lm', formula= y~x) +
        ggpubr::theme_pubr()
    
    my_plots[[paste0(x, "_vs_iBAQ_scatter")]] <- ggplot(
        my_codon_stats, aes(x = iBAQ, y = !!as.name(x))) +
        geom_point() +
        scale_x_log10() +
        scale_y_log10() +
        geom_smooth(method='lm', formula= y~x) +
        ggpubr::theme_pubr()
    
}

#library(seqinr)
#data(caitab)
#my_preferred_codon <- caitab %>%
#    data.frame(.) %>%
#    tibble::rownames_to_column(.data = ., var = "CODON") %>%
#    dplyr::full_join(
#        x = ., y = my_w %>% dplyr::select(., CODON, AA, L, syn = w))
my_other_caitab <- data.table::fread(
    input = "/mnt/storage/kxmna01/data/Synechocystis_6frame/Codon_usage_other_species/Codon_preference_other_species.txt",
    sep = "\t", quote = "", header = TRUE, stringsAsFactors = FALSE)
my_preferred_codon <- my_other_caitab %>%
    dplyr::full_join(
        x = my_w %>% dplyr::select(., CODON, AA, L, Synechocystis_sp = w), y = .)
#        x = my_preferred_codon, y = .)
#detach("package:seqinr", unload = TRUE)

my_preferred_codon_long <- my_preferred_codon %>%
    tidyr::pivot_longer(
        data = ., cols = c(-CODON, -AA, -L),
        names_to = "Species", values_to = "W")
my_preferred_codon_long$Species <- factor(
    x = my_preferred_codon_long$Species,
    levels = unique(my_preferred_codon_long$Species),
    ordered = TRUE)

my_palette <- grDevices::hcl.colors(n = 9, palette = "Fall")
my_big_cols <- c(
    "#387eb8", "#404040", "#e21e25", "#fbaf3f", "#d1d2d4", "#246E39", "#753B94")

my_plots[["w_across_species"]] <- ggplot(
    my_preferred_codon_long,
    aes(x = Species, y = CODON, fill = W)) +
    geom_tile(colour = "white") +
    facet_grid(rows = vars(AA), scales = "free_y", space = "free_y") +
    scale_fill_gradientn(colours = my_palette) +
    ggpubr::theme_pubr() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))

my_plots[["w_across_species_violin"]] <- ggplot(
    my_preferred_codon_long,
    aes(x = Species, y = W, fill = Species)) +
    geom_violin(alpha = 0.3) +
    geom_boxplot(alpha = 0.7, width = 0.1) +
    scale_fill_manual(values = my_big_cols) +
    ggpubr::theme_pubr() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))

for (x in unique(my_preferred_codon_long$Species)) {
    my_cor <- cor.test(
        x = my_preferred_codon$Synechocystis_sp,
        y = my_preferred_codon[[x]], method = "spearman")
    my_plots[[paste0(x, "_vs_syn_scatter")]] <- ggplot(
        my_preferred_codon, aes(x = Synechocystis_sp, y = !!as.name(x))) +
        geom_point() +
        geom_smooth(method='lm', formula= y~x) +
        ggpubr::theme_pubr() +
        ggtitle(paste0("Spearman rho = ", my_cor$estimate, " (p = ", my_cor$p.value, ")"))
}

pdf("Phylostrata_characteristics.pdf", 10, 10)
my_plots
dev.off()

data.table::fwrite(
    x = my_preferred_codon, file = "my_preferred_codon.txt",
    append = FALSE, quote = FALSE, sep = "\t",
    row.names = FALSE, col.names = TRUE)

data.table::fwrite(
    x = my_codon_stats, file = "my_codon_stats.txt",
    append = FALSE, quote = FALSE, sep = "\t",
    row.names = FALSE, col.names = TRUE)

data.table::fwrite(
    x = my_uco, file = "codon_usage_indices.txt",
    append = FALSE, quote = FALSE, sep = "\t",
    row.names = FALSE, col.names = TRUE)


