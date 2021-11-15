



my_summyrised <- readRDS(
    file = "H:/data/Synechocystis_6frame/SummarizedExp/MultiAssay.RDS")

my_fasta_f <- list(
    "H:/data/Synechocystis_6frame/Genome/micro_proteins_Synechocystis_sp_PCC6803_20180419.fasta",
    "H:/data/Synechocystis_6frame/Genome/Synechocystis_sp_PCC_6803_cds_aa.fasta"
)

all_ids <- lapply(
    c(
        "proteinGroups"),# "pxd_nolabel_proteinGroups", "pxd005085_proteinGroups",
        #"pxd014662_1_proteinGroups", "pxd014662_2_proteinGroups"),
    function(x) {
        rownames(assays(my_summyrised)[[x]]) %>%
            base::strsplit(x = ., split = ";", fixed = TRUE) %>%
            unlist(.) %>%
            unique(.)
    }) %>%
    unlist(.) %>%
    unique(.)

my_fasta <- lapply(my_fasta_f, function(x) {
    seqinr::read.fasta(
        file = x, seqtype = "AA", as.string = TRUE, whole.header = FALSE)
}) %>%
    unlist(., recursive = FALSE)

my_fasta_filt <- my_fasta[!names(my_fasta) %in% all_ids]

my_headers <- lapply(my_fasta_filt, function(x) {
    attr(x = x, which = "Annot") %>%
        sub("^>", "", .)
}) %>%
    unlist(.)

seqinr::write.fasta(
    sequences = my_fasta_filt, names = my_headers,
    file.out = "2021-11-15_Never_identified_in_Scycode.fasta",
    open = "w", nbchar = 60, as.string = TRUE)




