


library(magrittr)
library(ggplot2)
library(ggbio)
library(GenomicRanges)

my_plots <- list()

bin <- 10000

my_grange_f <- "/mnt/storage/kxmna01/data/Synechocystis_6frame/GRanges/Ref_prot_grange.RDS"

my_grange <- readRDS(my_grange_f)

my_genome_f <- "/mnt/storage/kxmna01/data/Synechocystis_6frame/GRanges/Genome_grange.RDS"

genome_grange <- readRDS(my_genome_f)

my_fasta_f <- "/mnt/storage/kxmna01/data/Synechocystis_6frame/2022-02-02_Hidden_proteome/Never_identified_in_never.fasta"

my_fasta <- seqinr::read.fasta(
    file = my_fasta_f, seqtype = "AA", as.string = TRUE)

my_data <- as.data.frame(my_grange)

toplot <- my_data %>%
    dplyr::mutate(
        .,
        start_round = plyr::round_any(x = start, accuracy = bin, f = floor)+1,
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
    ),
    end = start_round+(bin-1), strand = "*")
toplot$threshold <- quantile(x = toplot$Percentage, 0.975)

new_gr <- with(
    data = toplot,
    expr = GenomicRanges::GRanges(
        seqnames = seqnames,
        ranges = IRanges::IRanges(start_round, end),
        strand = strand))
seqlevels(new_gr) <- seqlevels(genome_grange) %>% as.character(.)
seqinfo(new_gr) <- seqinfo(genome_grange)
values(new_gr) <- toplot %>%
    dplyr::select(
        ., -seqnames, -start_round, -end, -strand)

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

my_plots[["Percentage_bar"]] <- ggplot(
    toplot, aes(
        x = start_round, y = Percentage, fill = Type, colour = Type, label = Label)) +
    geom_bar(stat = "identity", position = "dodge") +
    ggrepel::geom_text_repel() +
    ggpubr:::theme_pubr() +
    facet_grid(cols = vars(seqnames), scales = "free", space = "free_x") +
    scale_fill_manual(values = c("#387eb8", "#404040")) +
    scale_colour_manual(values = c("#387eb8", "#404040"))

my_plots[["circos"]] <- ggplot() +
    layout_circle(
        genome_grange, geom = "ideo", fill = "gray70",
        radius = 30, trackWidth = 2) +
    layout_circle(
        genome_grange, geom = "scale", size = 4,
        radius = 33, trackWidth = 2) +
    layout_circle(
        genome_grange, geom = "text", size = 7, aes(label = seqnames),
        angle = 0, radius = 37, trackWidth = 5) +
    layout_circle(
        new_gr,
        aes(y = Percentage, fill = Type, group = Type),
        geom = "bar",
        radius = 22, trackWidth = 8) +
    layout_circle(
        new_gr,
        aes(y = threshold),
        geom = "line",
        radius = 22, trackWidth = 8, linetype = "dashed", colour = "#e21e25") +
    scale_fill_manual(values = c("#387eb8", "#404040"))

pdf(paste0("Hidden_proteins_locations_bin", bin, ".pdf"), 10, 10)
my_plots
dev.off()

data.table::fwrite(
    x = toplot, file = paste0("Hidden_proteins_locations_bin", bin, ".txt"),
    quote = FALSE, append = FALSE, sep = "\t",
    row.names = FALSE, col.names = TRUE)


