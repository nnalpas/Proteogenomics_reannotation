


library(rBLAST)
library(magrittr)
library(GenomicRanges)

db <- blast(
    db = "H:/data/Synechocystis_6frame/Phylostratigraphy/1148.faa",
    type = "blastp")

micro_prot <- readAAStringSet(
    filepath = "C:/Users/kxmna01/Desktop/sORFS-Ribo-seq.fasta",
    format = "fasta")

my_res <- predict(db, micro_prot)

my_best <- my_res %>%
    dplyr::group_by(., QueryID) %>%
    dplyr::mutate(., Count = dplyr::n_distinct(SubjectID)) %>%
    dplyr::arrange(., E, dplyr::desc(`Perc.Ident`)) %>%
    dplyr::slice(., 1) %>%
    dplyr::ungroup(.)

pep_gr <- readRDS("H:/data/Synechocystis_6frame/GRanges/Peptides_grange.RDS")

micro_prot_gr <- data.frame(id = names(micro_prot)) %>%
    tidyr::separate(
        data = ., col = id,
        into = c("chr_id", "range", "strand"),
        sep = ":", remove = FALSE) %>%
    tidyr::separate(
        data = ., col = range,
        into = c("start", "end"), sep = "-") %>%
    dplyr::mutate(
        ., chromosome = dplyr::case_when(
            chr_id == "NC_000911.1" ~ "chr1",
            TRUE ~ sub("\\.", "_", tolower(chr_id))),
        start = as.integer(start),
        end = as.integer(end))

# Create a GRanges object for all entries with coordinates
micro_grange <- with(
    data = micro_prot_gr,
    expr = GRanges(
        seqnames = chromosome,
        ranges = IRanges(start, end),
        strand = strand))

# Add seqinfo to the created GRanges object
seqlevels(micro_grange) <- seqlevels(pep_gr) %>% as.character(.)
seqinfo(micro_grange) <- seqinfo(pep_gr)

# Add values to the created GRanges object
values(micro_grange) <- micro_prot_gr %>%
    dplyr::select(
        ., -chromosome, -start, -end, -strand)

micro_prot_pep <- subsetByOverlaps(
    pep_gr, micro_grange, type = "within") %>%
    mcols(.) %>%
    as.data.frame(.) %>%
    dplyr::mutate(., Protein = Proteins) %>%
    tidyr::separate_rows(data = ., Protein, sep = ", ")

my_best_pep <- dplyr::left_join(
    x = data.frame(QueryID = names(micro_prot)), y = my_best) %>%
    dplyr::select(., QueryID, SubjectID, `Perc.Ident`, E) %>%
    dplyr::left_join(
        x = .,
        y = micro_prot_pep %>% dplyr::select(., Protein, Sequence, group, NoveltyReason),
        by = c("SubjectID" = "Protein")) %>%
    dplyr::group_by(., QueryID) %>%
    dplyr::mutate(., Peptide_count = dplyr::n_distinct(na.omit(Sequence))) %>%
    dplyr::summarise_all(~paste0(unique(.), collapse = ";")) %>%
    dplyr::ungroup(.) %>%
    dplyr::mutate_all(~sub("NA", "", .))

data.table::fwrite(
    x = my_best_pep, file = "sORFS-Ribo-seq_MSidentification.txt",
    append = FALSE, quote = FALSE, sep = "\t",
    row.names = FALSE, col.names = TRUE)


