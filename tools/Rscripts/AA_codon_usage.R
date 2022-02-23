


library(magrittr)
library(ggplot2)
library(Biostrings)

my_plots <- list()
my_cols <- c(
    `Yes` = "#387eb8",
    `No` = "#3B3B3B")
my_big_cols <- c(
    "#387eb8", "#404040", "#e21e25", "#fbaf3f", "#d1d2d4", "#246E39", "#753B94")
my_palette <- grDevices::hcl.colors(n = 9, palette = "Fall")

my_nucl_f <- "H:/data/Synechocystis_6frame/Genome/Synechocystis_sp_PCC_6803_gene_nucl.fasta"

my_prot_f <- "H:/data/Synechocystis_6frame/Genome/Synechocystis_sp_PCC_6803_cds_aa.fasta"

my_nucl <- seqinr::read.fasta(
    file = my_nucl_f, seqtype = "DNA", as.string = FALSE)

my_prot <- seqinr::read.fasta(
    file = my_prot_f, seqtype = "AA", as.string = FALSE)

my_aa_freq <- unlist(my_prot) %>%
    table(.) %>%
    as.data.frame(.) %>%
    set_colnames(c("AA", "Freq")) %>%
    dplyr::arrange(., dplyr::desc(Freq)) %>%
    dplyr::mutate(
        ., `O-phosphorylable` = ifelse(AA %in% c("S", "T", "Y"), "Yes", "No"))

my_aa_freq$AA <- factor(
    x = my_aa_freq$AA, levels = my_aa_freq$AA, ordered = TRUE)

my_plots[["aa_freq"]] <- ggplot(
    my_aa_freq,
    aes(x = AA, y = Freq, fill = `O-phosphorylable`,
        colour = `O-phosphorylable`, label = Freq)) +
    geom_bar(stat = "identity", position = "dodge", alpha = 0.7) +
    geom_text(vjust = -0.3) +
    ggpubr::theme_pubr() +
    scale_fill_manual(values = my_cols) +
    scale_colour_manual(values = my_cols)

my_codons <- lapply(my_nucl, function(x) {
    codon <- c()
    res <- ""
    if (length(x) %% 3 != 0) {
        warning("Sequence not multiple of 3...")
    } else {
        for (y in 1:length(x)) {
            res %<>%
                paste0(., x[[y]])
            if ((y %% 3) == 0) {
                codon %<>%
                    c(., res)
                res <- ""
            }
        }
    }
    return(codon)
})

my_codon_freq <- unlist(my_codons) %>%
    data.frame(Name = names(.), Codon = .) %>%
    tidyr::separate(
        data = ., col = Name, into = c("Proteins", "Position"), sep = c(7))

my_trans <- data.frame(
    Codon = unique(my_codon_freq$Codon))
my_trans$AA_noInit <- translate(
    x = DNAStringSet(my_trans$Codon),
    genetic.code = getGeneticCode(id_or_name2 = "11"),
    no.init.codon = TRUE) %>%
    as.character(.)
my_trans$AA_Init <- translate(
    x = DNAStringSet(my_trans$Codon),
    genetic.code = getGeneticCode(id_or_name2 = "11"),
    no.init.codon = FALSE) %>%
    as.character(.)

my_codon_freq %<>%
    dplyr::left_join(x = ., y = my_trans, by = "Codon") %>%
    dplyr::mutate(., AA = ifelse(Position == 1, AA_Init, AA_noInit)) %>%
    dplyr::group_by(., Codon, AA) %>%
    dplyr::summarise(., Freq = dplyr::n()) %>%
    dplyr::ungroup(.) %>%
    dplyr::mutate(
        ., Total = sum(Freq),
        Perc = Freq*100/Total) %>%
    tidyr::separate(
        data = ., col = Codon,
        into = c("First position", "Second position", "Third position"),
        sep = c(1,2), remove = FALSE) %>%
    dplyr::arrange(., dplyr::desc(Freq)) %>%
    dplyr::group_by(., AA) %>%
    dplyr::mutate(., Rank = as.character(1:dplyr::n())) %>%
    dplyr::ungroup(.)

my_codon_freq$AA <- factor(
    x = my_codon_freq$AA,
    levels = c(as.character(my_aa_freq$AA), "*"),
    ordered = TRUE)

my_plots[["codon_freq"]] <- ggplot(
    my_codon_freq,
    aes(x = `Second position`, y = `Third position`, fill = Freq,
        label = round(Perc, 2))) +
    geom_tile(colour = "white") +
    geom_text() +
    facet_wrap(facets = vars(`First position`)) +
    ggpubr::theme_pubr() +
    scale_fill_gradientn(
        colours = my_palette)

my_plots[["aa_codon_freq"]] <- ggplot(
    my_codon_freq,
    aes(
        x = Codon, y = Freq, fill = Rank,
        colour = Rank, label = round(Perc, 1))) +
    geom_bar(stat = "identity", position = "dodge", alpha = 0.7) +
    geom_text(vjust = -0.3) +
    ggpubr::theme_pubr() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
    facet_grid(cols = vars(AA), scales = "free_x", space = "free_x") +
    scale_fill_manual(values = my_big_cols) +
    scale_colour_manual(values = my_big_cols)

my_plots[["aa_codon_freq_stack"]] <- ggplot(
    my_codon_freq,
    aes(
        x = AA, y = Freq, fill = Rank,
        colour = Rank, label = Codon)) +
    geom_bar(stat = "identity", position = "stack", alpha = 0.7) +
    geom_text(position = "stack", vjust = -0.3, colour = "black") +
    ggpubr::theme_pubr() +
    scale_fill_manual(values = my_big_cols) +
    scale_colour_manual(values = my_big_cols)

pdf("AA_codon_usage.pdf", 10, 10)
my_plots
dev.off()


