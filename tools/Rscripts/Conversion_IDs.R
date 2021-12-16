


my_fasta_f <- "H:/data/Synechocystis_6frame/Genome/Synechocystis_sp_PCC_6803_cds_aa.fasta"

my_f <- "H:/data/Synechocystis_6frame/EggnogMapper/Synechocystis_UniProt_annotation.tab"

my_fasta <- seqinr::read.fasta(
    file = my_fasta_f, seqtype = "AA", as.string = TRUE)

my_custom_annot <- data.table::fread(
    input = my_f, sep = "\t", quote = "", header = TRUE,
    stringsAsFactors = FALSE, colClasses = "character")

head(my_custom_annot$Entry)
grep(",|;|-| ", my_custom_annot$Entry, value = TRUE)

head(my_custom_annot$`Gene names`)
grep(",|;|-| ", my_custom_annot$`Gene names`, value = TRUE)
table(is.na(my_custom_annot$`Gene names`) | my_custom_annot$`Gene names` == "" | my_custom_annot$`Gene names` == "-")
my_custom_annot$`Gene names`[is.na(my_custom_annot$`Gene names`) | my_custom_annot$`Gene names` == "" | my_custom_annot$`Gene names` == "-"]


head(my_custom_annot$`Gene names  (ordered locus )`)
grep(",|;|-| ", my_custom_annot$`Gene names  (ordered locus )`, value = TRUE)


gene_crossmap <- my_custom_annot %>%
    dplyr::select(., `Entry`, `Gene names`, `Gene names  (ordered locus )`, Sequence) %>%
    tidyr::pivot_longer(
        data = ., cols = c(`Gene names`, `Gene names  (ordered locus )`),
        values_to = "genes") %>%
    tidyr::separate_rows(data = ., genes, sep = "(; | )") %>%
    dplyr::filter(., genes != "") %>%
    dplyr::select(., -name) %>%
    unique(.)

my_db <- data.table::data.table(
    ID = names(my_fasta), Sequence = unlist(my_fasta))

annot_to_gene <- my_db %>%
    dplyr::inner_join(
        x = ., y = gene_crossmap, by = c("ID" = "genes", "Sequence"))

annot_to_gene <- my_db %>%
    dplyr::filter(., !ID %in% annot_to_gene$ID) %>%
    dplyr::inner_join(
        x = ., y = dplyr::select(gene_crossmap, -genes), by = c("Sequence")) %>%
    dplyr::bind_rows(annot_to_gene, .)

annot_to_gene <- my_db %>%
    dplyr::filter(., !ID %in% annot_to_gene$ID) %>%
    dplyr::inner_join(
        x = ., y = dplyr::select(gene_crossmap, -Sequence), by = c("ID" = "genes")) %>%
    dplyr::bind_rows(annot_to_gene, .)

annot_to_gene %<>%
    dplyr::group_by(., ID) %>%
    dplyr::mutate(
        ., UniProtID = paste0(na.omit(unique(Entry)), collapse = ";"),
        Duplication = dplyr::n_distinct(Entry)) %>%
    dplyr::group_by(., Entry) %>%
    dplyr::mutate(
        ., ReverseDuplication = dplyr::n_distinct(ID)) %>%
    dplyr::ungroup(.)

annot_to_gene <- my_custom_annot %>%
    dplyr::select(., Entry_kegg = Entry, `Cross-reference (KEGG)`) %>%
    tidyr::separate_rows(data = ., `Cross-reference (KEGG)`, sep = ";") %>%
    dplyr::mutate(
        ., `Cross-reference (KEGG)` = sub(
            "sy.:", "", `Cross-reference (KEGG)`)) %>%
    dplyr::filter(
        ., !is.na(`Cross-reference (KEGG)`) &
            `Cross-reference (KEGG)` != "") %>%
    dplyr::full_join(
        x = annot_to_gene, y = ., by = c("ID" = "Cross-reference (KEGG)")) %>%
    dplyr::mutate_all(~tidyr::replace_na(data = ., replace = ""))

dim(my_db[!my_db$ID %in% annot_to_gene$ID, ])
dim(annot_to_gene[annot_to_gene$Entry != annot_to_gene$Entry_kegg, ])

annot_to_gene %<>%
    dplyr::mutate(., Entry = dplyr::case_when(
        Entry_kegg == "" ~ Entry,
        Entry == "" ~ Entry_kegg,
        TRUE ~ Entry
    )) %>%
    dplyr::filter(., ID %in% names(my_fasta)) %>%
    dplyr::select(., ID, MainUniProtID = Entry, DupUniProtID = UniProtID) %>%
    unique(.)

missing_ids <- my_fasta[!names(my_fasta) %in% annot_to_gene$ID]

seqinr::write.fasta(
    sequences = missing_ids, names = names(missing_ids),
    file.out = "Syn_proteins_without_Uniprot.fasta",
    open = "w", nbchar = 60, as.string = TRUE)

my_final_annot <- annot_to_gene %>%
    dplyr::left_join(
        x = ., y = my_custom_annot, by = c("MainUniProtID" = "Entry")) %>%
    unique(.) %>%
    dplyr::select(., -Sequence)

data.table::fwrite(
    x = my_final_annot,
    file = sub(".tab", "_formatted_2021-12-14.txt", my_f),
    append = FALSE, quote = FALSE, sep = "\t",
    row.names = FALSE, col.names = TRUE)


