


library(magrittr)
library(GO.db)

opt <- list()
opt$annotation <- "H:/data/Synechocystis_6frame/EggnogMapper/Synechocystis_UniProt_annotation_formatted_2021-12-21.txt"
opt$ecdb <- "D:/Local_databases/X-databases/EC/enzclass.txt"
opt$interpro = "D:/Local_databases/X-databases/Interpro/interpro_parsed.txt"
opt$panther = "D:/Local_databases/X-databases/Panther/Panther_parsed.txt"
opt$pirsf <- "D:/Local_databases/X-databases/PIRSF/pirsfinfo_parsed.txt"
opt$prosite <- "D:/Local_databases/X-databases/Prosite/2021-12/prosite_parsed.txt"
opt$tigrfam <- "D:/Local_databases/X-databases/Tigrfam/hmm_PGAP.tsv"
opt$output = "H:/data/Synechocystis_6frame/EggnogMapper"



### Data import ----------------------------------------------------------

annot <- myfilelist <- opt$annotation %>%
    base::strsplit(x = ., split = ",", fixed = TRUE) %>%
    unlist(.) %>%
    lapply(X = ., FUN = function(x) {
        data.table::fread(
            input = x,
            sep = "\t",
            header = TRUE,
            stringsAsFactors = FALSE,
            colClasses = "character")
    }) %>%
    plyr::ldply(., dplyr::bind_rows, .id = NULL)

# Import the different parsed database
ec <- data.table::fread(
    input = opt$ecdb, sep = NULL, quote = "", header = FALSE,
    stringsAsFactors = FALSE, blank.lines.skip = TRUE)
interpro <- data.table::fread(
    input = opt$interpro, sep = "\t", quote = "", header = TRUE,
    stringsAsFactors = FALSE)
panther <- data.table::fread(
    input = opt$panther, sep = "\t", quote = "", header = TRUE,
    stringsAsFactors = FALSE)
pirsf <- data.table::fread(
    input = opt$pirsf, sep = "\t", quote = "", header = TRUE,
    stringsAsFactors = FALSE)
prosite <- data.table::fread(
    input = opt$prosite, sep = "\t", quote = "", header = TRUE,
    stringsAsFactors = FALSE)
tigrfam <- data.table::fread(
    input = opt$tigrfam, sep = "\t", quote = "", header = TRUE,
    stringsAsFactors = FALSE)



### Reorganising UniProt download annotation -----------------------------

unique_id <- length(unique(annot$ID))

delimiters <- apply(annot, 2, function(x) {
    unique(unlist(na.omit(stringr::str_extract_all(x, " |;|,|:|-|\\[|\\.|\\/|\\="))))
})

if (any(lengths(delimiters) > 1)) {
    warning("Possible multiple delimiters!")
}

col_to_split <- names(delimiters)[lengths(delimiters) > 0]
delimiter_check <- lapply(col_to_split, function(c) {
    lapply(delimiters[[c]], function(x) {
        head(grep(x, annot[[c]], value = TRUE, fixed = TRUE), n = 1)
    }) %>%
        set_names(delimiters[[c]]) %>%
        plyr::ldply(., data.table::data.table) %>%
        set_colnames(c("delimiter", "example_match"))
}) %>%
    set_names(col_to_split) %>%
    plyr::ldply(., data.table::data.table, .id = "column")

data.table::fwrite(
    x = delimiter_check, file = "delimiter_check.txt",
    append = FALSE, quote = FALSE, sep = "\t",
    row.names = FALSE, col.names = TRUE)

# Perform action as described in file
delimiter_actions <- data.table::fread(
    input = "H:/data/Synechocystis_6frame/EggnogMapper/UniProt_delimiter_check.txt",
    sep = "\t", quote = "", header = TRUE,
    stringsAsFactors = FALSE, strip.white = FALSE) %>%
    dplyr::filter(., action != "")
delimiter_actions %<>%
    split(x = ., f = delimiter_actions$column)

annot_format <- data.table::data.table(ID = unique(annot$ID))

for (x in names(delimiters)[-1]) {
    
    my_col_format <- annot %>%
        dplyr::select(., ID, !!as.name(x))
    
    if (x %in% names(delimiter_actions)) {
        for (a in 1:nrow(delimiter_actions[[x]])) {
            
            if (delimiter_actions[[x]][a, ][["action"]] == "separate_rows") {
                my_col_format %<>%
                    tidyr::separate_rows(
                        data = ., -ID,
                        sep = delimiter_actions[[x]][a, ][["delimiter"]])
            } else if (delimiter_actions[[x]][a, ][["action"]] == "separate") {
                my_col_format %<>%
                    tidyr::separate(
                        data = ., col = -ID, into = c("name", "value"),
                        sep = delimiter_actions[[x]][a, ][["delimiter"]],
                        extra = "merge") %>%
                    dplyr::mutate(., name = paste(x, name, sep = "-")) %>%
                    dplyr::mutate_all(~gsub(";", ",", .)) %>%
                    dplyr::group_by(., ID, name) %>%
                    dplyr::summarise(., value = paste0(value, collapse = ";")) %>%
                    dplyr::ungroup(.) %>%
                    tidyr::pivot_wider(
                        data = ., names_from = "name", values_from = "value")
            } else if (delimiter_actions[[x]][a, ][["action"]] == "gsub") {
                my_col_format[[x]] %<>%
                    gsub(
                        pattern = delimiter_actions[[x]][a, ][["delimiter"]],
                        replacement = delimiter_actions[[x]][a, ][["substitution"]],
                        x = .)
            } else {
                stop(paste0(
                    "Action: '", delimiter_actions[[x]][a, ][["action"]],
                    "' not allowed!"))
            }
            
            my_col_format %<>%
                dplyr::mutate_all(~sub("^ +", "", .)) %>%
                dplyr::mutate_all(~sub(" +$", "", .)) %>%
                dplyr::mutate_all(~gsub('\"', "", .)) %>%
                dplyr::mutate_all(~sub(";$", "", .)) %>%
                dplyr::mutate_all(~sub("\\.$", "", .))
            
            if (x %in% colnames(my_col_format)) {
                my_col_format %<>%
                    dplyr::filter(., !is.na(!!as.name(x)) & !!as.name(x) != "")
            }
            
        }
    }
    
    if (nrow(my_col_format) != length(unique(my_col_format$ID))) {
        if (any(grepl(";", my_col_format[[x]]))) {
            warning("There are semi-colon remaining!")
        }
        my_col_format %<>%
            dplyr::group_by(., ID) %>%
            dplyr::summarise_all(~paste0(unique(.), collapse = ";")) %>%
            dplyr::ungroup(.)
    }
    
    annot_format %<>%
        dplyr::left_join(x = ., y = my_col_format, by = "ID")
    
}

annot_format %<>%
    dplyr::rowwise(.) %>%
    dplyr::mutate(., Characterization = dplyr::case_when(
        grepl("uncharact|unknow", `Protein names`, ignore.case = T) ~ "Uncharacterized",
        grepl(paste(ID, "protein"), `Protein names`, ignore.case = T) ~ "Probably uncharacterized",
        TRUE ~ "Characterized"
    )) %>%
    dplyr::ungroup(.) %>%
    dplyr::mutate_all(~tidyr::replace_na(data = ., replace = ""))



### Formatting EC categories ---------------------------------------------

# 
my_ec <- ec %>%
    dplyr::filter(., !grepl("^#", `V1`)) %>%
    tidyr::separate(
        data = ., col = "V1", into = c("EC", "EC name"), sep = "\\.-  +") %>%
    dplyr::mutate(
        ., EC = gsub(" ", "", EC) %>% gsub("\\.-", "", .),
        `EC name` = gsub(";", ",", `EC name`)) %>%
    dplyr::mutate(
        ., Level = (1+stringr::str_count(string = EC, pattern = "\\.")),
        `EC split` = EC) %>%
    tidyr::separate(
        data = ., col = `EC split`,
        into = c("1", "2", "3"), sep = "\\.") %>%
    dplyr::mutate(., `EC1` = `1`) %>%
    tidyr::unite(
        data = ., col = "EC2", `1`, `2`, sep = ".",
        remove = FALSE, na.rm = TRUE) %>%
    tidyr::unite(
        data = ., col = "EC3", `1`, `2`, `3`, sep = ".",
        remove = FALSE, na.rm = TRUE) %>%
    dplyr::select(., -`1`, -`2`, -`3`)
my_ec_format <- my_ec %>%
    dplyr::filter(., Level == 1) %>%
    dplyr::select(., `EC1` = EC, `EC1 name` = `EC name`) %>%
    dplyr::left_join(
        x = my_ec,
        y = .)
my_ec_format <- my_ec %>%
    dplyr::filter(., Level == 2) %>%
    dplyr::select(., `EC2` = EC, `EC2 name` = `EC name`) %>%
    dplyr::left_join(
        x = my_ec_format,
        y = .)
my_ec_format <- my_ec %>%
    dplyr::filter(., Level == 3) %>%
    dplyr::select(., `EC3` = EC, `EC3 name` = `EC name`) %>%
    dplyr::left_join(
        x = my_ec_format,
        y = .)
my_ec_format %<>%
    tidyr::unite(
        data = ., col = "EC2 label", `EC1 name`, `EC2 name`,
        sep = " ", remove = FALSE, na.rm = TRUE) %>%
    tidyr::unite(
        data = ., col = "EC3 label", `EC1 name`, `EC2 name`, `EC3 name`,
        sep = " ", remove = FALSE, na.rm = TRUE) %>%
    dplyr::select(
        ., EC, `EC level 1` = `EC1`, `EC level 2` = `EC2`,
        `EC level 3` = `EC3`, `EC level 1 name` = `EC1 name`,
        `EC level 2 name` = `EC2 label`, `EC level 3 name` = `EC3 label`)

annot_to_ec <- annot_format %>%
    dplyr::select(., ID, `EC number`) %>%
    dplyr::filter(
        ., !is.na(`EC number`) & `EC number` != "-" & `EC number` != "") %>%
    tidyr::separate_rows(data = ., `EC number`, sep = ";") %>%
    dplyr::mutate(., `EC short` = sub("(.+)\\..+", "\\1", `EC number`))

annot_to_ec_final <- dplyr::left_join(
    x = annot_to_ec, y = my_ec_format, by = c("EC short" = "EC")) %>%
    dplyr::select(., -`EC short`) %>%
    dplyr::group_by(., ID) %>%
    dplyr::summarise_all(~paste0(na.omit(.), collapse = ";")) %>%
    dplyr::ungroup(.)



### Formatting GO categories ---------------------------------------------

# 
go_ancestors <- c(
    as.list(GOBPANCESTOR), as.list(GOCCANCESTOR), as.list(GOMFANCESTOR))

go_df <- as.list(GOTERM) %>%
    lapply(., function(x) {
        data.table::data.table(
            Ontology = Ontology(x), Term = Term(x), stringsAsFactors = FALSE)
    }) %>%
    plyr::ldply(., data.table::data.table, .id = "IDs") %>%
    dplyr::mutate(., Term = gsub(";", ",", Term))

# Identify GO ancestor
annot_to_go <- annot_format %>%
    dplyr::select(., ID, `Gene ontology IDs`) %>%
    dplyr::filter(
        ., !is.na(`Gene ontology IDs`) & `Gene ontology IDs` != "-" &
            `Gene ontology IDs` != "") %>%
    tidyr::separate_rows(data = ., `Gene ontology IDs`, sep = ",") %>%
    dplyr::rowwise(.) %>%
    dplyr::mutate(
        ., 
        Ancestor = tryCatch(
            expr = {paste0(c(`Gene ontology IDs`, go_ancestors[[`Gene ontology IDs`]]), collapse = ";")},
            error = function(cond) {return(`Gene ontology IDs`)}))

# For all GO (incl. ancestor) get ontology and terms
annot_to_go %<>%
    tidyr::separate_rows(data = ., Ancestor, sep = ";") %>%
    dplyr::left_join(x = ., y = go_df, by = c("Ancestor" = "IDs")) %>%
    dplyr::filter(., Ancestor != "all" & !is.na(Term))

annot_to_go_final <- annot_to_go %>%
    dplyr::select(., -`Gene ontology IDs`) %>%
    unique(.) %>%
    tidyr::pivot_longer(data = ., cols = c(Ancestor, Term)) %>%
    dplyr::mutate(
        ., Ontology = paste0("GO", Ontology),
        name = ifelse(name == "Ancestor", NA_character_, name)) %>%
    tidyr::unite(
        data = ., col = "key", Ontology, name, sep = " ", na.rm = TRUE) %>%
    dplyr::group_by(., ID, `key`) %>%
    dplyr::summarise_all(~paste0(na.omit(.), collapse = ";")) %>%
    dplyr::ungroup(.) %>%
    tidyr::pivot_wider(data = ., names_from = "key", values_from = "value")



### Formatting Interpro categories ---------------------------------------

my_interpro_format <- interpro %>%
    set_colnames(sub("ENTRY", "Interpro", colnames(.)))

annot_to_interpro <- annot_format %>%
    dplyr::select(., ID, InterPro = `Cross-reference (InterPro)`) %>%
    dplyr::filter(
        ., !is.na(InterPro) & InterPro != "-" & InterPro != "") %>%
    tidyr::separate_rows(data = ., InterPro, sep = ";")

annot_to_interpro_final <- dplyr::left_join(
    x = annot_to_interpro, y = my_interpro_format,
    by = c("InterPro" = "Interpro_AC")) %>%
    dplyr::group_by(., ID) %>%
    dplyr::summarise_all(~paste0(na.omit(.), collapse = ";")) %>%
    dplyr::ungroup(.)



### Formatting Panther categories ----------------------------------------

my_panther_format <- panther %>%
    set_colnames(sub("^", "Panther ", colnames(.)))

annot_to_panther <- annot_format %>%
    dplyr::select(., ID, Panther = `Cross-reference (PANTHER)`) %>%
    dplyr::filter(
        ., !is.na(Panther) & Panther != "-" & Panther != "") %>%
    tidyr::separate_rows(data = ., Panther, sep = ";")

annot_to_panther_final <- dplyr::left_join(
    x = annot_to_panther, y = my_panther_format,
    by = c("Panther" = "Panther PANTHER ID")) %>%
    dplyr::group_by(., ID) %>%
    dplyr::summarise_all(~paste0(na.omit(.), collapse = ";")) %>%
    dplyr::ungroup(.)



### Formatting PIRSF categories ------------------------------------------

my_pirsf_format <- pirsf %>%
    dplyr::select(., -curation_status, -parent) %>%
    set_colnames(sub("^", "PIRSF ", colnames(.)))

annot_to_pirsf <- annot_format %>%
    dplyr::select(., ID, PIRSF = `Cross-reference (PIRSF)`) %>%
    dplyr::filter(
        ., !is.na(PIRSF) & PIRSF != "-" & PIRSF != "") %>%
    tidyr::separate_rows(data = ., PIRSF, sep = ";")

annot_to_pirsf_final <- dplyr::left_join(
    x = annot_to_pirsf, y = my_pirsf_format,
    by = c("PIRSF" = "PIRSF PIRSF_number")) %>%
    dplyr::group_by(., ID) %>%
    dplyr::summarise_all(~paste0(na.omit(.), collapse = ";")) %>%
    dplyr::ungroup(.)



### Formatting Prosite categories ------------------------------------------

my_prosite_format <- prosite %>%
    dplyr::select(., -DT) %>%
    set_colnames(sub("^", "Prosite ", colnames(.)))

annot_to_prosite <- annot_format %>%
    dplyr::select(., ID, Prosite = `Cross-reference (PROSITE)`) %>%
    dplyr::filter(
        ., !is.na(Prosite) & Prosite != "-" & Prosite != "") %>%
    tidyr::separate_rows(data = ., Prosite, sep = ";")

annot_to_prosite_final <- dplyr::left_join(
    x = annot_to_prosite, y = my_prosite_format,
    by = c("Prosite" = "Prosite AC")) %>%
    dplyr::group_by(., ID) %>%
    dplyr::summarise_all(~paste0(na.omit(.), collapse = ";")) %>%
    dplyr::ungroup(.)



### Formatting TIGRFAM categories ----------------------------------------

my_tigrfam_format <- tigrfam %>%
    dplyr::select(
        ., source_identifier, label, family_type, product_name) %>%
    set_colnames(sub("^", "TIGRFAM ", colnames(.)))

annot_to_tigrfam <- annot_format %>%
    dplyr::select(., ID, TIGRFAM = `Cross-reference (TIGRFAMs)`) %>%
    dplyr::filter(
        ., !is.na(TIGRFAM) & TIGRFAM != "-" & TIGRFAM != "") %>%
    tidyr::separate_rows(data = ., TIGRFAM, sep = ";")

annot_to_tigrfam_final <- dplyr::left_join(
    x = annot_to_tigrfam, y = my_tigrfam_format,
    by = c("TIGRFAM" = "TIGRFAM source_identifier")) %>%
    dplyr::group_by(., ID) %>%
    dplyr::summarise_all(~paste0(na.omit(.), collapse = ";")) %>%
    dplyr::ungroup(.)



### Compile all annotations ----------------------------------------------

annot_final <- annot_format %>%
    dplyr::select(
        ., -`Gene ontology IDs`, -`EC number`,
        -`Cross-reference (InterPro)`, -`Cross-reference (PANTHER)`,
        -`Cross-reference (PIRSF)`, -`Cross-reference (PROSITE)`,
        -`Cross-reference (TIGRFAMs)`) %>%
    dplyr::left_join(x = ., y = annot_to_ec_final, by = "ID") %>%
    dplyr::left_join(x = ., y = annot_to_go_final, by = "ID") %>%
    dplyr::left_join(x = ., y = annot_to_interpro_final, by = "ID") %>%
    dplyr::left_join(x = ., y = annot_to_panther_final, by = "ID") %>%
    dplyr::left_join(x = ., y = annot_to_pirsf_final, by = "ID") %>%
    dplyr::left_join(x = ., y = annot_to_prosite_final, by = "ID") %>%
    dplyr::left_join(x = ., y = annot_to_tigrfam_final, by = "ID") %>%
    dplyr::mutate_all(~tidyr::replace_na(data = ., replace = ""))

if (nrow(annot_final) != unique_id) {
    stop("Final number of rows not equal to input!")
}



### Export the results data ----------------------------------------------

# Export annotation usable in Perseus
write.table(
    x = annot_final,
    file = paste0(opt$output, "/UniProt_annotation_2_perseus.txt"),
    quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)

# Save session
save.image(paste0(
    opt$output, "/UniProt_annotation_2_perseus.RData"))



### END ------------------------------------------------------------------

# Time scripts end
print(paste("Completed ", format(Sys.time(), "%Y-%m-%d"), sep = ""))


