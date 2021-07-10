#!/usr/bin/env Rscript

# This script identify all novel peptides and ORFs from a MaxQuant search
# as input as well as the fasta files used for the search



### Environment set-up ---------------------------------------------------

# Start with clean environment
rm(list = ls())

# Define current time
date_str <- format(Sys.time(), "%Y-%m-%d")
print(paste("Start", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))

# Define the current user
user <- Sys.info()[["user"]]



### List of required packages -----------------------------------------------

# Source the custom user functions
if (interactive()) {
    source(
        file = paste(
            "C:/Users",
            user,
            "Documents/GitHub/Proteogenomics_reannotation/",
            "tools/Rscripts/helper.R",
            sep = "/"))
} else {
    source(
        file = paste(
            Sys.getenv("HOME"),
            "bin/helper.R",
            sep = "/"))
}

# Load the required packages (or install if not already in library)
# Note: the select function from dplyr and UniProt.ws is overlapping,
# therefore order of package loading is important
library(plyr)
library(dplyr)
library(seqinr)
library(bit64)
library(magrittr)
library(data.table)
library(splitstackshape)
library(stringr)
library(optparse)
library(ggplot2)
library(gtable)
library(ggpubr)
library(UpSetR)
library(plotly)
library(grid)
library(gridExtra)
library(purrr)
library(foreach)
library(doParallel)
library(rmarkdown)
library(cleaver)



### Parameters setting up ------------------------------------------------

# Define input parameters (interactively or from command line)
if (interactive()) {
    
    # Define the list of input parameters
    opt <- list(
        output = readline(
            prompt = "Define the output directory!"),
        maxquant = choose.dir(
            caption = "Choose the MaxQuant txt folder!"),
        reference = choose.files(
            caption = "Choose fasta file containing known proteins!",
            multi = FALSE),
        novel = choose.files(
            caption = "Choose fasta file containing ORF proteins!",
            multi = FALSE),
        threads = readline(prompt = "How many cores to use?") %>% as.integer())
    
} else {
    
    # Define the list of command line parameters
    option_list <- list(
        make_option(
            opt_str = c("-o", "--output"),
            type = "character", default = NULL,
            help = "Output directory", metavar = "character"),
        make_option(
            opt_str = c("-m", "--maxquant"),
            type = "character", default = NULL,
            help = "MaxQuant txt results folder", metavar = "character"),
        make_option(
            opt_str = c("-r", "--reference"),
            type = "character", default = NULL,
            help = "Reference protein fasta file", metavar = "character"),
        make_option(
            opt_str = c("-n", "--novel"),
            type = "character", default = NULL,
            help = "ORF fasta file", metavar = "character"),
        make_option(
            opt_str = c("-t", "--threads"),
            type = "integer", default = "",
            help = "Number of cores to use",
            metavar = "character"))
    
    # Parse the parameters provided on command line by user
    opt_parser <- OptionParser(option_list = option_list)
    opt <- parse_args(opt_parser)
    
}

# Check whether inputs parameters were provided
if (
    is.null(opt$output) | is.null(opt$maxquant) |
    is.null(opt$reference) | is.null(opt$novel)) {
    
    print_help(opt_parser)
    stop("All arguments must be supplied!")
    
}
if (
    identical(opt$threads, NULL) |
    identical(opt$threads, "") |
    identical(opt$threads, integer(0))) {
    
    warning("Default number of threads will be used!")
    opt$threads <- 1
    
}
registerDoParallel(cores = opt$threads)
print(paste("Number of threads registered:", getDoParWorkers()))

# Create output directory if not already existing
dir.create(opt$output)



### Data import ----------------------------------------------------------

# Import the maxquant evidence table
evid <- mq_read(
    path = opt$maxquant,
    name = "evidence.txt",
    integer64 = "double")

# Keep only evidences with at least 1 MS/MS
warning(paste(
    "There are ",
    table(c(evid$`MS/MS count` == 0, TRUE))[["TRUE"]]-1,
    "PSM with no MS/MS (matching), these will be ignored!"))
evid %<>%
    dplyr::filter(., `MS/MS count` >= 1)

# Import the maxquant peptide table
pep <- mq_read(
    path = opt$maxquant,
    name = "peptides.txt",
    integer64 = "double")

# Import the maxquant proteingroups table
pg <- mq_read(
    path = opt$maxquant,
    name = "proteinGroups.txt",
    integer64 = "double")

# Keep evidence IDs for all protein groups that are NOT only identified by site
pg_not_sites <- unique(as.integer(
    pg[pg$`Only identified by site` == "", ][["id"]]))

# Import the fasta files
fasta <- c(Known = opt$reference, Novel = opt$novel) %>%
    purrr::map(
        .x = ., .f = seqinr::read.fasta, seqtype = "AA", as.string = TRUE)
names(fasta$Known) %<>%
    uni_id_clean(.)

# Import all sites data
site_files <- list.files(
    path = opt$maxquant, pattern = "Sites.txt", full.names = TRUE) %>%
    set_names(sub("Sites\\.txt", "", basename(.)))
if (length(site_files) > 0) {
    sites <- lapply(X = site_files, function(x) {
        mq_read(
            path = dirname(x),
            name = basename(x) %>%
                sub("\\(", "\\\\(", .) %>%
                sub("\\)", "\\\\)", .),
            integer64 = "double")
    })
}



### Novel peptide identification -----------------------------------------

# Define most represented cleave sites and the enzymatic rule
digest_rule <- digestion_rule(
    x = pep[["Amino acid before"]])

# Define missed cleavages rule
mc <- c(min(
    pep[, grep("Missed cleavages", colnames(pep))]):(max(
        pep[, grep("Missed cleavages", colnames(pep))])+2))

# Digest all proteins into peptide to get their location and
# their database of origin
digest_db <- prot_digest_foreach(
    databases = fasta,
    custom = digest_rule, mc = mc)

# Get all identified protein IDs
prot_ids <- evid %>%
    strsplit(x = .[["Proteins"]], split = ";", fixed = TRUE) %>%
    unlist(.) %>%
    grep("CON__|REV__", ., invert = TRUE, value = TRUE) %>%
    unique(.)

# Check whether all protein IDs can be found
if (all(!prot_ids %in% digest_db$Proteins)) {
    stop("ID clean-up failed, check peptide and digested table!")
} else {
    warning(paste0(
        "Not all protein IDs are represented in digested table",
        " files provided, ",
        length(which(
            !prot_ids %in% digest_db$Proteins)),
        " IDs are missing!"))
}

# Make sure that the 'amino acid before' column contains a single character
# since a single character is removed via substring later on
if (any(nchar(pep$`Amino acid before`) > 1)) {
    stop(paste(
        "The peptide table column 'Amino acid before' contains more than",
        "one character, this is not allowed!"))
}
evid <- pep %>%
    dplyr::select(., Sequence, `Amino acid before`) %>%
    dplyr::left_join(x = evid, y = ., by = "Sequence")

# List all peptide sequence (including sequence with aa before)
# this is the MAIN data.frame that will be used to define each peptide type
pep_comp <- dplyr::bind_rows(
    data.table(
        table = "evid",
        evid[, c(
            "Sequence", "Amino acid before", "Proteins", "Leading proteins",
            "Reverse", "Potential contaminant", "Protein group IDs")],
        stringsAsFactors = FALSE),
    data.table(
        table = "pep",
        pep[!pep$Sequence %in% evid$Sequence, c(
            "Sequence", "Amino acid before", "Proteins", "Leading razor protein",
            "Reverse", "Potential contaminant", "Protein group IDs")],
        stringsAsFactors = FALSE) %>%
        dplyr::rename(., `Leading proteins` = `Leading razor protein`)
) %>%
    unique(.) %>%
    dplyr::mutate(
        ., DigestSeq = paste0(`Amino acid before`, Sequence))

# Digest all proteins into peptide to get their location and
# their database of origin
digest_datab <- unique_to_database(
    digest = digest_db,
    pep = pep_comp$Sequence)

# Retrieve the peptide location for those peptide that are missing
# the amino acid before (typically a methionine)
pep_comp %<>%
    dplyr::left_join(
        x = ., y = unique(digest_datab[, c("Sequence", "Dbuniqueness")]))

# Check if some peptides are missing
digest_datab_missing <- data.table::data.table()
if (any(is.na(pep_comp$Dbuniqueness))) {
    missing_seq <- pep_comp %>%
        dplyr::filter(., is.na(Dbuniqueness)) %>%
        .[["DigestSeq"]] %>%
        unique(.)
    digest_datab_missing <- unique_to_database(
        digest = digest_db,
        pep = missing_seq)
    pep_comp <- digest_datab_missing %>%
        dplyr::select(., Sequence, Dbuniqueness_miss = Dbuniqueness) %>%
        unique(.) %>%
        dplyr::left_join(
            x = pep_comp, y = ., by = c("DigestSeq" = "Sequence")) %>%
        tidyr::unite(
            data = ., col = "Dbuniqueness", Dbuniqueness, Dbuniqueness_miss,
            sep = ";", remove = TRUE, na.rm = TRUE) %>%
        dplyr::mutate(
            ., Dbuniqueness = ifelse(
                Dbuniqueness == "", NA_character_, Dbuniqueness))
}

# Check that the value in Dbuniqueness are valid
if (any(
    !na.omit(pep_comp$Dbuniqueness) %in% unique(digest_datab$Dbuniqueness))) {
    not_allowed <- pep_comp$Dbuniqueness %>%
        unique(.) %>%
        na.omit(.) %>%
        .[!. %in% unique(digest_datab$Dbuniqueness)] %>%
        paste0(., collapse = " - ")
    stop(paste0("The Dbuniqueness value is not allowed: ", not_allowed, "!"))
}

# Compile all peptide positions
pep_pos <- digest_datab_missing %>%
    dplyr::mutate(., Sequence = substring(Sequence, 2), start = start + 1) %>%
    dplyr::bind_rows(digest_datab, .) %>%
    dplyr::select(., -id) %>%
    dplyr::mutate(., id = Sequence)

# Compile peptide type info for each peptide sequence
pep_comp %<>%
    dplyr::ungroup(.) %>%
    dplyr::rename(., DatabID = Dbuniqueness)
if (length(pep_comp$Sequence) != length(unique(pep_comp$Sequence))) {
    stop("Duplicate sequence value in 'pep_comp'!")
}
pep_comp %<>%
    dplyr::mutate(
        ., OnlyIdBySite = ifelse(
            `Protein group IDs` %in% pg_not_sites, TRUE, FALSE),
        group = dplyr::case_when(
            grepl("Known", DatabID) ~ "Known",
            grepl("CON__", Proteins) ~ "Contaminant",
            grepl("Novel", DatabID) ~ "Novel",
            grepl("REV__", Proteins) |
                grepl("REV__", `Leading proteins`) ~ "Reverse",
            TRUE ~ NA_character_
        ),
        Database = dplyr::case_when(
            group %in% c("Known", "Contaminant") ~ "Target",
            group == "Reverse" ~ "Decoy",
            group == "Novel" ~ "Novel",
            TRUE ~ NA_character_
        ))

# Add the group and database info back into all main tables
evid_match <- pep_comp %>%
    dplyr::select(., Sequence, OnlyIdBySite, group, Database) %>%
    dplyr::left_join(x = evid, y = ., by = "Sequence")
pep_match <- pep_comp %>%
    dplyr::select(., Sequence, OnlyIdBySite, group, Database) %>%
    dplyr::left_join(x = pep, y = ., by = "Sequence")
pep_pos <- pep_comp %>%
    dplyr::select(., Sequence, OnlyIdBySite, group, Database) %>%
    dplyr::left_join(x = pep_pos, y = ., by = "Sequence")

# Print warning for unidentified sequence origin
print(paste(
    "There are", length(which(is.na(evid_match$group))),
    "NA values, these need to be checked!", sep = " "))

# Print the repartition of peptide per group
print(table(evid_match$group, useNA = "always"))

# Print warning for unidentified sequence origin
print(paste(
    "There are", length(which(is.na(pep_match$group))),
    "NA values, these need to be checked!", sep = " "))

# Print the repartition of peptide per group
print(table(pep_match$group, useNA = "always"))
pep_match %<>%
    dplyr::filter(., !is.na(group))



### Sites identification -------------------------------------------------

# In case there are sites
if (exists("sites")) {
    
    for (x in names(sites)) {
        
        site_id <- paste(x, "site IDs")
        ids_cross_map <- evid %>%
            dplyr::select(
                ., `Evidence IDs` = id, `Site IDs` = !!as.name(site_id)) %>%
            dplyr::filter(
                ., !is.na(`Site IDs`) & `Site IDs` != "") %>%
            tidyr::separate_rows(
                data = ., `Site IDs`, sep = ";", convert = TRUE)
        ids_cross_map <- sites[[x]] %>%
            dplyr::select(., `Site IDs` = id, `Evidence IDs`) %>%
            dplyr::filter(
                ., !is.na(`Evidence IDs`) & `Evidence IDs` != "") %>%
            tidyr::separate_rows(
                data = ., `Evidence IDs`, sep = ";", convert = TRUE) %>%
            dplyr::bind_rows(ids_cross_map, .) %>%
            unique(.)
        ids_cross_map_format <- evid_match %>%
            dplyr::select(
                ., `Evidence IDs` = id, Sequence, OnlyIdBySite,
                group, Database) %>%
            unique(.) %>%
            dplyr::left_join(x = ids_cross_map, y = ., by = "Evidence IDs") %>%
            dplyr::mutate(
                ., `Peptide details` = paste0(Sequence, " (", group, ")")) %>%
            dplyr::group_by(., `Site IDs`) %>%
            dplyr::summarise(
                ., `Peptide details` = paste0(
                    unique(`Peptide details`), collapse = ";"),
                group = dplyr::case_when(
                    any(group == "Known") ~ "Known",
                    any(group == "Contaminant") ~ "Contaminant",
                    any(group == "Reverse") ~ "Reverse",
                    any(group == "Novel") ~ "Novel",
                    TRUE ~ NA_character_
                ),
                Database = dplyr::case_when(
                    any(Database == "Target") ~ "Target",
                    any(Database == "Decoy") ~ "Decoy",
                    any(Database == "Novel") ~ "Novel",
                    TRUE ~ NA_character_
                ))
        
        site_pos <- ids_cross_map_format %>%
            dplyr::rename(., `id` = `Site IDs`) %>%
            dplyr::left_join(
                x = sites[[x]], y = ., by = c("id"))
        
        # Prepare sites location info
        site_pos %<>%
            dplyr::mutate(
                ., Proteins = ifelse(Proteins == "", Protein, Proteins),
                `Positions within proteins` = ifelse(
                    `Positions within proteins` == "",
                    Position, `Positions within proteins`)) %>%
            tidyr::separate_rows(
                data = ., Proteins, `Positions within proteins`, sep = ";") %>%
            dplyr::mutate(
                ., `Site ID` = paste(
                    Proteins, `Positions within proteins`,
                    `Amino acid`, sep = "~")) %>%
            dplyr::mutate(
                ., start = `Positions within proteins`,
                end = `Positions within proteins`) %>%
            dplyr::select(
                ., id = `Site ID`, Proteins, start, end, `Peptide details`,
                group, Database, Sequence = `Amino acid`,
                `Localization prob`, PEP, Score, `Sequence window`,
                dplyr::ends_with("Probabilities"),
                dplyr::starts_with("Intensity"))
        
        # Save the peptide location data
        saveRDS(
            object = site_pos,
            file = paste(
                opt$output, "/", x, "Sites_location.RDS", sep = ""))
        
    }
    
}



### Results visualisation ------------------------------------------------

# Define path to markdown file
if (interactive()) {
    rmd_file <- paste(
        "C:/Users",
        user,
        "Documents/GitHub/Proteogenomics_reannotation",
        "tools/Rscripts/GetNovelEntries.rmd",
        sep = "/")
} else {
    rmd_file <- paste(
        Sys.getenv("HOME"),
        "bin/GetNovelEntries.rmd",
        sep = "/")
}

# Generate the markdown report
report_markdown(
    rmd_file = rmd_file,
    params = list(
        fastas = fasta,
        pep_match = pep_match,
        evid_match = evid_match %>%
            dplyr::filter(., !is.na(group) & !is.na(Database))),
    format = "html_document",
    ext = "html")



### Export the results data ----------------------------------------------

# Save the group mapping data
saveRDS(
    object = evid_match,
    file = paste(
        opt$output, "/", "Group_evidence.RDS", sep = ""))

# Save the peptide location data
saveRDS(
    object = pep_pos,
    file = paste(
        opt$output, "/", "Peptides_location.RDS", sep = ""))

# Export complete evidence info for these novel evidence
write.table(
    x = evid_match,
    file = paste(
        opt$output, "/", "Group_evidence.txt", sep = ""),
    quote = FALSE,
    sep = "\t",
    row.names = FALSE,
    col.names = TRUE)

# Export complete evidence info for these novel evidence
my_id_novels <- evid_match %>%
    dplyr::filter(., group == "Novel") %>%
    cSplit(
        indt = ., splitCols = "Proteins", sep = ";",
        direction = "long", fixed = TRUE) %>%
    .[["Proteins"]] %>%
    unique(.) %>%
    as.character(.)
write.table(
    x = my_id_novels,
    file = paste(
        opt$output, "/", "Novel_ORF.txt", sep = ""),
    quote = FALSE,
    sep = "\t",
    row.names = FALSE,
    col.names = FALSE)

# Save subset fasta file for identified novel entries
my_novel_fasta <- seqinr::read.fasta(
    file = opt$novel, seqtype = "AA", as.string = TRUE)
my_novel_fasta_id <- my_novel_fasta[my_id_novels]
seqinr::write.fasta(
    sequences = my_novel_fasta_id,
    names = unlist(lapply(my_novel_fasta_id, function(x) {
        sub("^>", "", attr(x = x, which = "Annot"))
    })),
    file.out = paste(
        opt$output, "/",
        sub(".fasta", "_identified.fasta", basename(opt$novel)),
        sep = ""),
    open = "w", nbchar = 60, as.string = TRUE)

# Save session
save.image(paste(
    opt$output, "/", "GetNovelEntries.RData", sep = ""))



### END ------------------------------------------------------------------

# Time scripts end
print(paste("Completed ", format(Sys.time(), "%Y-%m-%d"), sep = ""))


