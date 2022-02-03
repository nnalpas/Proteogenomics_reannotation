


library(magrittr)
library(ggplot2)

my_plots <- list()

my_grange_f <- "H:/data/Synechocystis_6frame/GRanges/Ref_prot_grange.RDS"

my_grange <- readRDS(my_grange_f)

my_fasta_f <- "H:/data/Synechocystis_6frame/2022-02-02_Hidden_proteome/Never_identified_in_never.fasta"

my_fasta <- seqinr::read.fasta(
    file = my_fasta_f, seqtype = "AA", as.string = TRUE)

my_data <- as.data.frame(my_grange)

toplot <- my_data %>%
    dplyr::mutate(
        .,
        start_round = plyr::round_any(x = start, accuracy = 40000, f = floor),
        Type = ifelse(id %in% names(my_fasta), "Hidden", "Others")) %>%
    dplyr::group_by(., Type) %>%
    dplyr::mutate(
        ., Total = dplyr::n()) %>%
    dplyr::group_by(., seqnames, start_round, Type) %>%
    dplyr::summarise(
        ., Count = dplyr::n(), Percentage = Count*100/unique(Total),
        Loci_names = paste0(id, collapse = ", ")) %>%
    dplyr::ungroup(.) %>%
    dplyr::mutate(., Label = dplyr::case_when(
        Percentage > 2.5 ~ Loci_names,
        TRUE ~ NA_character_
    ))

my_plots[["Count"]] <- ggplot(
    toplot, aes(
        x = start_round, y = Count, fill = Type, colour = Type, label = Label)) +
    geom_point(shape = 21) +
    ggrepel::geom_text_repel() +
    ggpubr:::theme_pubr() +
    facet_grid(cols = vars(seqnames), scales = "free", space = "free_x") +
    scale_fill_manual(values = c("#387eb8", "#404040")) +
    scale_colour_manual(values = c("#387eb8", "#404040"))

my_plots[["Percentage"]] <- ggplot(
    toplot, aes(
        x = start_round, y = Percentage, fill = Type, colour = Type, label = Label)) +
    geom_point(shape = 21) +
    ggrepel::geom_text_repel() +
    ggpubr:::theme_pubr() +
    facet_grid(cols = vars(seqnames), scales = "free", space = "free_x") +
    scale_fill_manual(values = c("#387eb8", "#404040")) +
    scale_colour_manual(values = c("#387eb8", "#404040"))

pdf("Hidden_proteins_locations.pdf", 10, 10)
my_plots
dev.off()


