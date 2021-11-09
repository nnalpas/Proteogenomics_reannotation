


library(magrittr)
library(ggplot2)

my_plots <- list()

my_f_tu <- "H:/data/Synechocystis_6frame/Kopf_2014_TU/Transcriptional_units_formatted_20201022.txt"

my_f_gr <- "H:/data/Synechocystis_6frame/GRanges/Ref_prot_grange.RDS"

my_f_fasta <- "H:/data/Synechocystis_6frame/Genome/Synechocystis_sp_PCC_6803_cds_aa.fasta"

my_f_assay <- "H:/data/Synechocystis_6frame/SummarizedExp/MultiAssay.RDS"

my_f_pg <- "H:/data/Synechocystis_6frame/MaxQuant/combined/txt/proteinGroups.txt"

my_tu <- data.table::fread(
    input = my_f_tu, sep = "\t", header = TRUE,
    stringsAsFactors = FALSE, colClasses = "character")

my_tu_map <- my_tu %>%
    dplyr::select(., id, Gene_name, Gene_count) %>%
    tidyr::separate_rows(data = ., Gene_name, sep = ", ", convert = FALSE) %>%
    dplyr::filter(., Gene_name != "")

my_fasta <- seqinr::read.fasta(
    file = my_f_fasta, seqtype = "AA", as.string = TRUE, whole.header = FALSE)

my_assay <- readRDS(my_f_assay)

my_pg <- data.table::fread(
    input = my_f_pg, sep = "\t", header = TRUE, quote = "",
    stringsAsFactors = FALSE, colClasses = "character")

my_gr <- readRDS(my_f_gr)

my_data <- data.frame(ranges(my_gr), strand = strand(my_gr)) %>%
    dplyr::full_join(x = ., y = my_tu_map, by = c("names" = "Gene_name"))

table(is.na(my_data$id))
table(is.na(my_data$start))

problem_ids <- unique(my_data[is.na(my_data$start), ][["names"]])

problem_ids[problem_ids %in% names(my_fasta)]

my_data_filt <- my_data %>%
    dplyr::filter(., !is.na(id) & !is.na(start)) %>%
    dplyr::group_by(., id) %>%
    dplyr::arrange(., start) %>%
    dplyr::mutate(
        ., `rank_+` = 1:dplyr::n(), `rank_-` = dplyr::n():1,
        Count_valid = Gene_count == dplyr::n(),
        Rank = ifelse(strand == "+", `rank_+`, `rank_-`))

table(is.na(as.double(my_pg$Intensity)))
table(is.na(as.double(my_pg$iBAQ)))

#my_abund <- assays(my_assay)[["proteinGroups"]] %>%
    #as.data.frame(.) %>%
    #tibble::rownames_to_column(.data = ., var = "id") %>%
    #tidyr::pivot_longer(
        #data = ., cols = -id,
        #names_to = "Conditions", values_to = "Intensity")
my_abund <- dplyr::select(my_pg, id = `Protein IDs`, dplyr::starts_with("iBAQ")) %>%
    dplyr::select(., -`iBAQ peptides`) %>%
    as.data.frame(.) %>%
    tidyr::pivot_longer(
        data = ., cols = -id,
        names_to = "Conditions", values_to = "Intensity") %>%
    dplyr::mutate(., Intensity = as.double(Intensity))

toplot <- my_abund %>%
    dplyr::mutate(., Quantified = !is.na(Intensity)) %>%
    plyr::ddply(
        .data = ., .variables = c("Conditions", "Quantified"),
        .fun = summarise, Count = dplyr::n(), .drop = FALSE)

ggplot(toplot, aes(x = Conditions, y = Count, fill = Quantified)) +
    geom_bar(stat = "identity", position = "dodge") +
    ggpubr::theme_pubr() +
    coord_flip()

my_abund_corr <- my_abund %>%
    tidyr::separate_rows(data = ., id, sep = ";") %>%
    dplyr::left_join(
        x = .,
        y = my_data_filt[, c("names", "Count_valid", "Rank", "Gene_count")],
        by = c("id" = "names")) %>%
    dplyr::filter(., !is.na(Rank))

for (x in unique(my_abund_corr$Conditions)) {
    
    toplot <- my_abund_corr %>%
        dplyr::filter(., Conditions == x)
    
    my_plots[[paste(x, "scatter all")]] <- ggplot(
        toplot, aes(x = log10(Intensity), y = Rank, fill = Count_valid, colour = Count_valid)) +
        geom_point() +
        ggpubr::theme_pubr() +
        xlab(x) +
        ggtitle("All TU")
    
    my_plots[[paste(x, "box all")]] <- ggplot(
        toplot, aes(x = Rank, y = log10(Intensity), group = Rank)) +
        geom_boxplot(position = "dodge") +
        ggpubr::theme_pubr() +
        xlab(x) +
        ggtitle("All TU")
    
}

my_abund_corr_filt <- my_abund_corr %>%
    dplyr::filter(., Gene_count > 1)

for (x in unique(my_abund_corr$Conditions)) {
    
    toplot <- my_abund_corr_filt %>%
        dplyr::filter(., Conditions == x)
    
    my_plots[[paste(x, "scatter gene > 1")]] <- ggplot(
        toplot, aes(x = log10(Intensity), y = Rank, fill = Count_valid, colour = Count_valid)) +
        geom_point() +
        ggpubr::theme_pubr() +
        xlab(x) +
        ggtitle("Tu with more than 1 gene")
    
    my_plots[[paste(x, "box gene > 1")]] <- ggplot(
        toplot, aes(x = Rank, y = log10(Intensity), group = Rank)) +
        geom_boxplot(position = "dodge") +
        ggpubr::theme_pubr() +
        xlab(x) +
        ggtitle("Tu with more than 1 gene")
    
}

pdf("Correlation_gene_order_versus_protein_abundance.pdf", 8, 8)
my_plots
dev.off()


