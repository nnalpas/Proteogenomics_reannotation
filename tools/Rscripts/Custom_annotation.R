


library(magrittr)

my_out <- "H:/data/Synechocystis_6frame/Custom_annotation"
my_fasta_f <- "H:/data/Synechocystis_6frame/Phylostratigraphy/1148.faa"
my_micro_fasta_f <- "H:/data/Synechocystis_6frame/Genome/micro_proteins_Synechocystis_sp_PCC6803_20180419.fasta"
my_mq_f <- "H:/data/Synechocystis_6frame/MQ_6frame_valid/combined/txt/"
my_mq_pride_f <- "H:/data/Synechocystis_6frame/ORF_validation"
my_novel_f <- "H:/data/Synechocystis_6frame/NoveltyExplain/ORF_novelty_reason.RDS"
my_position_f <- "H:/data/Synechocystis_6frame/ProtPosition/Ref_and_ORF_prot_coordinates.txt"
my_tu_f <- "H:/data/Synechocystis_6frame/Kopf_2014_TU/Transcriptional_units_per_gene.txt"
my_ps_f <- "H:/data/Synechocystis_6frame/EggnogMapper/EggNOG_annotation_2_perseus_PS_parsed.txt"

dir.create(my_out)



### Import all known and novel ORFs --------------------------------------

my_fasta <- seqinr::read.fasta(
    file = my_fasta_f, seqtype = "AA",
    as.string = TRUE, whole.header = FALSE)



### Novel ORF annotation -------------------------------------------------

my_novelty <- readRDS(my_novel_f)

my_novelty_final <- my_novelty %>%
    dplyr::select(
        ., Proteins, Database, PEPfilter,
        OnlyIdBySite, ORFNoveltyReason) %>%
    dplyr::mutate(
        ., Quality = dplyr::case_when(
            OnlyIdBySite & PEPfilter %in% c("class 1", "class 2") ~ "High quality",
            TRUE ~ "Low quality"),
        Explanation = sub("( \\(|;|\\||s\\/).*", "", ORFNoveltyReason)) %>%
    tidyr::pivot_longer(data = ., cols = c(Database, Quality, Explanation))



### Phosphoprotein identification ----------------------------------------

my_phospho_f <- "Phospho (STY)Sites.txt"

my_phospho <- data.table::fread(
    input = paste0(my_mq_f, my_phospho_f),
    sep = "\t", quote = "",
    header = TRUE, stringsAsFactors = FALSE,
    colClasses = "character")

my_phospho_final <- my_phospho %>%
    dplyr::filter(., Reverse != "+" & `Potential contaminant` != "+") %>%
    dplyr::select(., Proteins, `Localization prob`, Intensity) %>%
    tidyr::separate_rows(data = ., Proteins, sep = ";") %>%
    dplyr::mutate(
        ., Identified = "Identified phospho (STY)",
        #Localised = dplyr::case_when(
        #    `Localization prob` > 0.95 ~ "Localisation prob. > 0.95",
        #    `Localization prob` > 0.75 ~ "Localisation prob. > 0.75",
        #    TRUE ~ "Not localised"),
        Quantified = dplyr::case_when(
            !is.na(Intensity) & Intensity > 0 ~ "Quantified phospho (STY)",
            TRUE ~ "Not quantified phospho (STY)")) %>%
    tidyr::pivot_longer(
        data = ., cols = c(Identified, Quantified))#c(Identified, Localised, Quantified))



### Kinases annotation ---------------------------------------------------

my_headers <- lapply(my_fasta, function(x) {
    attr(x, "Annot")
}) %>%
    unlist(.)

my_kinases <- my_headers %>%
    grep(
        "kinase|sll1574 |sll1575 |slr1697 |slr0559 |sll0776 |slr1443 |slr1225 |slr0152 |sll0005 |sll1770 |slr0868 |slr1919 |sll0095 |slr0889 ",
        ., value = TRUE, ignore.case = TRUE) %>%
    sub("^>", "", .) %>%
    data.frame(Header = .) %>%
    tidyr::separate(
        data = ., col = Header,
        into = c("Proteins", "Genome", "Description", "Location"),
        sep = " \\| ", extra = "merge") %>%
    tidyr::separate_rows(data = ., Description, sep = ", ") %>%
    dplyr::mutate(
        ., Kinase = dplyr::case_when(
            grepl("kinase", Description, ignore.case = TRUE) ~ Description,
            TRUE ~ NA_character_)) %>%
    dplyr::group_by(., Proteins) %>%
    dplyr::summarise(
        ., Kinase = paste0(na.omit(unique(Kinase)), collapse = ";")) %>%
    dplyr::mutate(
        ., Kinase = dplyr::case_when(
            Kinase != "" ~ Kinase,
            TRUE ~ "serine/threonine kinase")) %>%
    dplyr::filter(., !grepl("inhibitor", Kinase))
my_kinases$Kinase %<>%
    sub("two-componen t sensor|two-component system sensory|two-component sen sor", "two-component sensor", .) %>%
    sub("histidin e", "histidine", .) %>%
    sub("protein kinase", "kinase", .)
my_kinases_final <- my_kinases %>%
    dplyr::mutate(
        ., Type = "Kinase",
        Target = dplyr::case_when(
            grepl("serine|threonine|tyrosine|histidine|aspartate|glutamate", Kinase) ~ "Protein kinase",
            TRUE ~ "Other kinase")) %>%
    tidyr::pivot_longer(data = ., cols = c(Type, Target, Kinase))



### Protein identification -----------------------------------------------

my_pg_f <- "proteinGroups.txt"

my_pg_files <- list.dirs(
    path = my_mq_pride_f, full.names = TRUE, recursive = FALSE) %>%
    set_names(basename(.)) %>%
    lapply(., function(x) paste(x, "/combined/txt/", my_pg_f, sep = "")) %>%
    unlist(.) %>%
    c(ScyCode = paste0(my_mq_f, my_pg_f))

my_pg <- lapply(my_pg_files, function(x) {
    data.table::fread(
        input = x, sep = "\t", quote = "",
        header = TRUE, stringsAsFactors = FALSE,
        colClasses = "character")
}) %>%
    plyr::ldply(., data.table::data.table, .id = "Processing")

all_prot <- data.frame(
    Proteins = rep(names(my_fasta), each = length(my_pg_files)),
    Processing = rep(names(my_pg_files), times = length(my_fasta)))

my_pg_final <- my_pg %>%
    dplyr::filter(., Reverse != "+" & `Potential contaminant` != "+") %>%
    dplyr::select(., Proteins = `Protein IDs`, Intensity, Processing) %>%
    tidyr::separate_rows(data = ., Proteins, sep = ";") %>%
    dplyr::left_join(x = all_prot, y = .) %>%
    dplyr::mutate(
        ., Identified = dplyr::case_when(
            !is.na(Intensity) ~ paste("Identified", Processing),
            TRUE ~ paste("Not identified", Processing)),
        Quantified = dplyr::case_when(
            !is.na(Intensity) & Intensity > 0 ~ paste("Quantified", Processing),
            TRUE ~ paste("Not quantified", Processing))) %>%
    dplyr::group_by(., Proteins) %>%
    dplyr::mutate(
        ., Identified_all = dplyr::case_when(
            any(!is.na(Intensity)) ~ "Identified",
            TRUE ~ "Not identified"),
        Quantified_all = dplyr::case_when(
            any(!is.na(Intensity) & Intensity > 0) ~ "Quantified",
            TRUE ~ "Not quantified")) %>%
    tidyr::pivot_longer(
        data = .,
        cols = c(Identified_all, Quantified_all, Identified, Quantified))



### Microproteins annotation ---------------------------------------------

my_micro <- seqinr::read.fasta(
    file = my_micro_fasta_f, seqtype = "AA",
    as.string = TRUE, whole.header = FALSE) %>%
    names(.) %>%
    data.table::data.table(Proteins = ., value = "Microproteins")



### Protein genomic position ---------------------------------------------

my_position <- data.table::fread(
    input = my_position_f,
    sep = "\t", quote = "",
    header = TRUE, stringsAsFactors = FALSE,
    colClasses = "character") %>%
    dplyr::select(., -ExactCoord, -Comment)



### Transcriptional units annotation -------------------------------------

my_tu <- data.table::fread(
    input = my_tu_f,
    sep = "\t", quote = "",
    header = TRUE, stringsAsFactors = FALSE,
    colClasses = "character") %>%
    tidyr::separate_rows(data = ., `TU ID`, `TU type`, sep = ";") %>%
    dplyr::group_by(., `Sense tags`) %>%
    dplyr::summarise_all(~paste0(unique(.), collapse = ";")) %>%
    dplyr::ungroup(.)



### SCyCode functional annotation ----------------------------------------

my_ps <- data.table::fread(
    input = my_ps_f,
    sep = "\t", quote = "",
    header = TRUE, stringsAsFactors = FALSE,
    colClasses = "character")



### Compile the custom annotation ----------------------------------------

my_custom_annot <- dplyr::bind_rows(
    my_pg_final, my_phospho_final, my_kinases_final,
    my_novelty_final, my_micro) %>%
    dplyr::group_by(., Proteins) %>%
    dplyr::summarise(
        ., Miscellaneous = paste0(
            na.omit(unique(value)), collapse = ";")) %>%
    dplyr::ungroup(.) %>%
    dplyr::rename(., `#query_name` = Proteins) %>%
    dplyr::left_join(x = ., y = my_position, by = c("#query_name" = "id")) %>%
    dplyr::left_join(x = ., y = my_tu, by = c("#query_name" = "Sense tags")) %>%
    dplyr::left_join(x = ., y = my_ps, by = "#query_name") %>%
    dplyr::mutate_all(~replace(x = ., list = is.na(.), values = ""))



### Export the custom annotation -----------------------------------------

data.table::fwrite(
    x = my_custom_annot,
    file = paste(my_out, "2021-12-21_annotation_session.txt", sep = "/"),
    append = FALSE, quote = FALSE, sep = "\t",
    row.names = FALSE, col.names = TRUE)

save.image(paste(my_out, "2021-12-21_annotation_session.RData", sep = "/"))



### END ------------------------------------------------------------------


