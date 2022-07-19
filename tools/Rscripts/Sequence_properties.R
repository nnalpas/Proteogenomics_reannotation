


rm(list = ls())

library(magrittr)
library(ggplot2)
library(data.table)

my_plots <- list()
my_cols <- c("#387eb8", "#3B3B3B", "#e21e25", "#fbaf3f", "#d1d2d4", "#834F96")

#my_files <- c(
#    ref = "H:/data/Synechocystis_6frame/Genome/Synechocystis_sp_PCC_6803_cds_aa.fasta",
#    micro_prot = "H:/data/Synechocystis_6frame/Genome/micro_proteins_Synechocystis_sp_PCC6803_20180419.fasta",
#    never_id = "H:/data/Synechocystis_6frame/2020-12-02_Never_identified/Never_identified.fasta",
#    novel_id = "H:/data/Synechocystis_6frame/Nuc_translation/Find0_Synechocystis_sp_PCC_6803_genome_FIXED_filter.fasta")

#my_files <- c(
#    ref = "H:/data/Synechocystis_6frame/Genome/Synechocystis_sp_PCC_6803_cds_aa.fasta",
#    ref_micro = "H:/data/Synechocystis_6frame/Genome/micro_proteins_Synechocystis_sp_PCC6803_20180419.fasta",
#    never_id = "H:/data/Synechocystis_6frame/2022-02-02_Hidden_proteome/Never_identified_in_never.fasta")

my_files <- c(
    ref = "H:/data/Srim_6frame_3rd_analysis/Genome/CP048261_CP048262_prot_sequence_FIXED3.fasta",
    novel = "H:/data/Srim_6frame_3rd_analysis/Novel_res/Find0_CP048261_CP048262_genome_sequence_FIXED_identified.fasta")

my_fastas <- lapply(X = my_files, FUN = function(x) {
    seqinr::read.fasta(file = x, seqtype = "AA", as.string = T) %>%
        set_names(sub(".+\\|(.+)\\|.+", "\\1", names(.)))
}) %>%
    set_names(names(my_files))

aaProp <- function(seq) {
    seq <- Peptides:::aaCheck(seq)
    seq <- lapply(seq, function(seq) {
        table(unlist(seq))/length(seq) %>%
            as.vector(.)
    })
    return(seq)
}

aaComp <- function(seq) {
    seq <- Peptides:::aaCheck(seq)
    seq <- lapply(seq, function(seq) {
        table(unlist(seq))
    })
    aacomp <- lapply(seq, function(seq) {
        AA <- matrix(ncol = 2, nrow = 9)
        rownames(AA) <- c("Tiny", "Small", "Aliphatic", 
                          "Aromatic", "NonPolar", "Polar", 
                          "Charged", "Basic", "Acidic")
        colnames(AA) <- c("Number", "Mole%")
        AA[1, 1] <- sum(seq[c("A", "C", "G", 
                              "S", "T")], na.rm = TRUE)
        AA[2, 1] <- sum(seq[c("A", "B", "C", 
                              "D", "G", "N", "P", "S", 
                              "T", "V")], na.rm = TRUE)
        AA[3, 1] <- sum(seq[c("A", "I", "L", 
                              "V")], na.rm = TRUE)
        AA[4, 1] <- sum(seq[c("F", "H", "W", 
                              "Y")], na.rm = TRUE)
        AA[5, 1] <- sum(seq[c("A", "C", "F", 
                              "G", "I", "L", "M", "P", 
                              "V", "W", "Y")], na.rm = TRUE)
        AA[6, 1] <- sum(seq[c("D", "E", "H", 
                              "K", "N", "Q", "R", "S", 
                              "T", "Z")], na.rm = TRUE)
        AA[7, 1] <- sum(seq[c("B", "D", "E", 
                              "H", "K", "R", "Z")], na.rm = TRUE)
        AA[8, 1] <- sum(seq[c("H", "K", "R")], 
                        na.rm = TRUE)
        AA[9, 1] <- sum(seq[c("B", "D", "E", 
                              "Z")], na.rm = TRUE)
        AA[, 2] <- (AA[, 1]/sum(seq) * 100)
        AA <- as.data.frame(round(AA, 3))
        return(setNames(object = AA[[2]], nm = rownames(AA)))
    })
    return(aacomp)
}

#aIndex <- function(seq) {
#    seq <- Peptides:::aaCheck(seq)
#    seq <- lapply(seq, function(seq) {
#        table(unlist(seq))/length(seq)
#    })
#    aliphatic <- unlist(lapply(seq, function(seq) {
#        c(aliphatic = sum(
#            c(
#                seq["A"],
#                (2.9 * seq["V"]),
#                3.9 * seq[c("L", "I")]),
#            na.rm = T) * 100)
#    }))
#    return(aliphatic)
#}

prot_stats <- lapply(my_fastas, function(x) {
    list(
        aaComp = aaComp(seq = x) %>% set_names(names(x)),
        aIndex = Peptides::aIndex(seq = x) %>% set_names(names(x)),
        aaProp = aaProp(seq = x) %>% set_names(names(x)),
        blosumIndices = Peptides::blosumIndices(seq = x) %>% set_names(names(x)),
        boman = Peptides::boman(seq = x) %>% set_names(names(x)),
        charge = Peptides::charge(seq = x) %>% set_names(names(x)),
        crucianiProperties = Peptides::crucianiProperties(seq = x) %>% set_names(names(x)),
        fasgaiVectors = Peptides::fasgaiVectors(seq = x) %>% set_names(names(x)),
        hydrophobicity = Peptides::hydrophobicity(seq = x, scale = "KyteDoolittle") %>% set_names(names(x)),
        instaIndex = Peptides::instaIndex(seq = x) %>% set_names(names(x)),
        kideraFactors = Peptides::kideraFactors(seq = x) %>% set_names(names(x)),
        lengthpep = Peptides::lengthpep(seq = x) %>% set_names(names(x)),
        mw = Peptides::mw(seq = x, monoisotopic = FALSE) %>% set_names(names(x)),
        pI = Peptides::pI(seq = x, pKscale = "EMBOSS") %>% set_names(names(x))
    )
})

prot_stats_df <- prot_stats %>%
    unlist(.) %>%
    data.frame(Names = names(.), Value = ., stringsAsFactors = FALSE) %>%
    tidyr::separate(
        data = ., col = "Names",
        into = c("Database", "Index", "Protein", "Param"),
        sep = "\\.", fill = "right") %>%
    dplyr::mutate(., Param = ifelse(is.na(Param), Index, Param))

prot_stats_df$Database <- factor(
    x = prot_stats_df$Database, levels = names(my_fastas), ordered = TRUE)

for (x in unique(prot_stats_df$Index)) {
    
    toplot <- prot_stats_df %>%
        dplyr::filter(., Index == x)
    
    my_plots[[paste0(x, "_violin")]] <- ggplot(
        data = toplot,
        mapping = aes(
           x = Database, y = Value, fill = Database, colour = Database)) +
        geom_violin(position = "dodge", alpha = 0.5, draw_quantiles = c(0.25, 0.5, 0.75)) +
        facet_wrap(facets = "Param") +
        ggtitle(x) +
        ggpubr::theme_pubr() +
        scale_fill_manual(values = my_cols[1:length(my_fastas)]) +
        scale_colour_manual(values = my_cols[1:length(my_fastas)])
    
    my_plots[[paste0(x, "_box")]] <- ggplot(
        data = toplot,
        mapping = aes(
            x = Database, y = Value, fill = Database, colour = Database)) +
        geom_boxplot(position = "dodge", alpha = 0.5) +
        facet_wrap(facets = "Param") +
        ggtitle(x) +
        ggpubr::theme_pubr() +
        scale_fill_manual(values = my_cols[1:length(my_fastas)]) +
        scale_colour_manual(values = my_cols[1:length(my_fastas)])
    
}

pdf("Sequence_AA_statistics.pdf", 10, 10)
my_plots
dev.off()

prot_stats_df_pivot <- prot_stats_df %>%
    tidyr::unite(data = ., col = "Names", Index, Param, sep = " ") %>%
    tidyr::pivot_wider(data = ., names_from = "Names", values_from = "Value")

data.table::fwrite(
    x = prot_stats_df_pivot, file = "Sequence_AA_statistics.txt",
    append = FALSE, quote = FALSE, sep = "\t",
    row.names = FALSE, col.names = TRUE)


