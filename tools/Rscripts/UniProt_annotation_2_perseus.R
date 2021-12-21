


library(magrittr)
library(GO.db)

opt <- list()
opt$annotation <- "H:/data/Synechocystis_6frame/EggnogMapper/Synechocystis_UniProt_annotation_formatted_2021-12-21.txt"
opt$ecdb <- "D:/Local_databases/X-databases/EC/enzclass.txt"
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

# Import the enzyme commission database
ec <- data.table::fread(
    input = opt$ecdb, sep = NULL, quote = "", header = FALSE,
    stringsAsFactors = FALSE, blank.lines.skip = TRUE)



### Reorganising UniProt download annotation -----------------------------

unique_id <- length(unique(annot$ID))

delimiters <- apply(annot, 2, function(x) {
    unique(unlist(na.omit(str_extract_all(x, " |;|,|:|-|\\[|\\.|\\/|\\="))))
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
    input = "delimiter_check.txt", sep = "\t", quote = "",
    header = TRUE, stringsAsFactors = FALSE, strip.white = FALSE) %>%
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



### Compile all annotations ----------------------------------------------

annot_final <- annot_format %>%
    dplyr::select(., -`Gene ontology IDs`, -`EC number`) %>%
    dplyr::left_join(x = ., y = annot_to_ec_final, by = "ID") %>%
    dplyr::left_join(x = ., y = annot_to_go_final, by = "ID") %>%
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


