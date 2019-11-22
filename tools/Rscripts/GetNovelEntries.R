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
library(plotly)
library(grid)
library(gridExtra)
library(purrr)
library(foreach)
library(doParallel)
library(shiny)



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
    table(evid$`MS/MS count` == 0)[["TRUE"]],
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
pg_not_sites <- pg %>%
    dplyr::filter(., `Only identified by site` == "") %>%
    cSplit(
        indt = ., splitCols = "Evidence IDs",
        sep = ";", direction = "long") %>%
    .[["Evidence IDs"]] %>%
    as.integer(.) %>%
    unique(.)

# Import the fasta files
fasta <- c(Known = opt$reference, Novel = opt$novel) %>%
    purrr::map(
        .x = ., .f = seqinr::read.fasta, seqtype = "AA", as.string = TRUE)
names(fasta$Known) %<>%
    uni_id_clean(.)



### Novel peptide identification -----------------------------------------

# Define most represented cleave sites and the enzymatic rule
digest_rule <- digestion_rule(
    x = pep[["Amino acid before"]])

# Define missed cleavages rule
mc <- c(min(
    pep[["Missed cleavages"]]):max(
        pep[["Missed cleavages"]]))

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

# Define all peptide sequence identified and include peptide
# that lost the first methionine
evid <- pep %>%
    dplyr::select(., Sequence, `Amino acid before`) %>%
    dplyr::left_join(x = evid, y = ., by = "Sequence")

# List all peptide sequence (including sequence with aa before)
all_possible_pep <- unique(c(
    evid[["Sequence"]],
    paste0(
        evid[["Amino acid before"]],
        evid[["Sequence"]])))

# Digest all proteins into peptide to get their location and
# their database of origin
digest_datab <- unique_to_database(
    digest = digest_db,
    pep = all_possible_pep)

# Clean up the environment of most memory hungry variables
# should be uncommented in case of RAM run out
#rm(digest_db)

# Retrieve the peptide location for those peptide that are missing
# the amino acid before (typically a methionine)
missing_pep <- evid %>%
    dplyr::filter(., !Sequence %in% digest_datab$Sequence) %>%
    mq_rev_con_filt(.) %>%
    dplyr::mutate(
        ., DigestSeq = paste0(`Amino acid before`, Sequence)) %>%
    dplyr::select(., Sequence, DigestSeq) %>%
    dplyr::left_join(
        x = ., y = digest_datab,
        by = c("DigestSeq" = "Sequence")) %>%
    dplyr::filter(., !is.na(id)) %>%
    dplyr::mutate(., start = start + 1) %>%
    dplyr::select(., -DigestSeq)
pep_pos <- dplyr::bind_rows(digest_datab, missing_pep) %>%
    dplyr::select(., -id) %>%
    dplyr::filter(
        ., Sequence %in% evid[["Sequence"]])

# Compile peptide type info for each peptide sequence
pep_comp <- pep_pos %>%
    dplyr::group_by(., Sequence) %>%
    dplyr::summarise(
        ., DatabID = paste(sort(unique(Dbuniqueness)), collapse = "-"))

# New dataframe to hold info about fasta of origin for each sequence
evid_match <- evid %>%
    dplyr::left_join(x = ., y = pep_comp, by = "Sequence") %>%
    dplyr::ungroup(.) %>%
    dplyr::mutate(
        ., 
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

# Print warning for unidentified sequence origin
print(paste(
    "There are", length(which(is.na(evid_match$group))),
    "NA values, these need to be checked!", sep = " "))

# Print the repartition of peptide per group
print(table(evid_match$group, useNA = "always"))



### Results export and visualisation -------------------------------------

# Open a file for plot visualisation
pdf(
    file = paste0(opt$output, "/", "GetNovelEntries.pdf"),
    width = 10, height = 10)

# Define colour coding for consistency across plots
group_colours <- c(
    `Known` = "#2180FC", `Contaminant` = "#fc9a21",
    `Reverse` = "#ABA9A9", `Novel` = "#F23D3D")
database_colours <- c(
    `Target` = "#2180FC", `Decoy` = "#ABA9A9", `Novel` = "#F23D3D")

# Histogram of evidence counts
# Histogram of evidence counts
toplot <- evid_match %>% #params[["evid_match"]] %>%
    dplyr::group_by(., group) %>%
    dplyr::summarise(., evid_count = n())
toplot$group <- factor(
    x = toplot$group, levels = names(group_colours), ordered = TRUE)
pl <- ggplot(
    toplot,
    aes(
        x = group, y = evid_count, fill = group,
        colour = group, label = evid_count)) +
    geom_bar(stat = "identity", position = "dodge", alpha = 0.7) +
    geom_text(vjust = -0.5) +
    xlab("Group") +
    ylab("Count") +
    labs(fill = "Group") +
    theme_pubr() +
    scale_y_log10() +
    scale_fill_manual(values = group_colours) +
    scale_colour_manual(guide = FALSE, values = group_colours)




pl <- plots_hist(
    data = toplot,
    key = "group",
    value = "evid_count",
    group = "group",
    fill = "group",
    main = "PSM count",
    xlabel = "Groups",
    ylabel = "Count (log scale)",
    textsize = 25,
    label = "evid_count",
    transf = "log10",
    bw = TRUE)
pl[[1]]

# Histogram of evidence counts
toplot <- evid_match %>%
    dplyr::group_by(., Database) %>%
    dplyr::summarise(., evid_count = n())
pl <- plots_hist(
    data = toplot,
    key = "Database",
    value = "evid_count",
    group = "Database",
    fill = "Database",
    main = "PSM count",
    xlabel = "Databases",
    ylabel = "Count (log scale)",
    textsize = 25,
    label = "evid_count",
    transf = "log10",
    bw = TRUE)
pl[[1]]

# Boxplot of evidence resolution
toplot <- evid_match %>%
    dplyr::select(., group, Resolution)
toplot <- filter_na(my_data = toplot, my_col = "Resolution")
pl <- plots_box(
    data = toplot,
    key = "group",
    value = "Resolution",
    main = "PSM resolution",
    textsize = 25,
    fill = "grey",
    xlabel = "Databases",
    ylabel = "Resolution",
    outlier_simplify = TRUE)
plot(pl[[1]])

# Boxplot of evidence resolution
toplot <- evid_match %>%
    dplyr::select(., Database, Resolution)
toplot <- filter_na(my_data = toplot, my_col = "Resolution")
pl <- plots_box(
    data = toplot,
    key = "Database",
    value = "Resolution",
    main = "PSM resolution",
    textsize = 25,
    fill = "grey",
    xlabel = "Databases",
    ylabel = "Resolution",
    outlier_simplify = TRUE)
plot(pl[[1]])

# Boxplot of evidence mass error (ppm)
toplot <- evid_match %>%
    dplyr::select(., group, `Mass error [ppm]`) %>%
    set_colnames(c("group", "Masserror"))
toplot <- filter_na(my_data = toplot, my_col = "Masserror")
pl <- plots_box(
    data = toplot,
    key = "group",
    value = "Masserror",
    main = "PSM mass error (ppm)",
    textsize = 25,
    fill = "grey",
    xlabel = "Databases",
    ylabel = "Mass error (ppm)",
    outlier_simplify = TRUE)
plot(pl[[1]])

# Bowplot of evidence mass error (ppm)
toplot <- evid_match %>%
    dplyr::select(., Database, `Mass error [ppm]`) %>%
    set_colnames(c("Database", "Masserror"))
toplot <- filter_na(my_data = toplot, my_col = "Masserror")
pl <- plots_box(
    data = toplot,
    key = "Database",
    value = "Masserror",
    main = "PSM mass error (ppm)",
    textsize = 25,
    fill = "grey",
    xlabel = "Databases",
    ylabel = "Mass error (ppm)",
    outlier_simplify = TRUE)
plot(pl[[1]])

# Boxplot of evidence PEP
toplot <- evid_match %>%
    dplyr::select(., group, PEP)
toplot <- filter_na(my_data = toplot, my_col = "PEP")
pl <- plots_box(
    data = toplot,
    key = "group",
    value = "PEP",
    main = "PSM PEP",
    textsize = 25,
    zoom = c(-0.001, quantile(toplot$PEP, 0.95)[[1]]),
    fill = "grey",
    xlabel = "Databases",
    ylabel = "PEP",
    outlier_simplify = TRUE)
plot(pl[[1]])
plot(pl[[2]])

# Boxplot of evidence PEP
toplot <- evid_match %>%
    dplyr::select(., Database, PEP)
toplot <- filter_na(my_data = toplot, my_col = "PEP")
pl <- plots_box(
    data = toplot,
    key = "Database",
    value = "PEP",
    main = "PSM PEP",
    textsize = 25,
    zoom = c(-0.001, quantile(toplot$PEP, 0.95)[[1]]),
    fill = "grey",
    xlabel = "Databases",
    ylabel = "PEP",
    outlier_simplify = TRUE)
plot(pl[[1]])
plot(pl[[2]])

# Boxplot of evidence score
toplot <- evid_match %>%
    dplyr::select(., group, Score)
toplot <- filter_na(my_data = toplot, my_col = "Score")
pl <- plots_box(
    data = toplot,
    key = "group",
    value = "Score",
    main = "PSM score",
    textsize = 25,
    fill = "grey",
    xlabel = "Databases",
    ylabel = "Score",
    outlier_simplify = TRUE)
plot(pl[[1]])

# Bowplot of evidence score
toplot <- evid_match %>%
    dplyr::select(., Database, Score)
toplot <- filter_na(my_data = toplot, my_col = "Score")
pl <- plots_box(
    data = toplot,
    key = "Database",
    value = "Score",
    main = "PSM score",
    textsize = 25,
    fill = "grey",
    xlabel = "Databases",
    ylabel = "Score",
    outlier_simplify = TRUE)
plot(pl[[1]])

# Close the picture file
dev.off()



### Generate the report --------------------------------------------------

# Generate the markdown report
report_markdown(
    rmd_file = "GetNovelEntries.rmd",
    params = list(
        fastas = fasta,
        evid_match = evid_match),
    format = "html_document",
    ext = "html")

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
write.table(
    x = evid_match %>%
        dplyr::filter(., group == "Novel") %>%
        cSplit(
            indt = ., splitCols = "Proteins", sep = ";",
            direction = "long", fixed = TRUE) %>%
        .[["Proteins"]] %>%
        unique(.),
    file = paste(
        opt$output, "/", "Novel_ORF.txt", sep = ""),
    quote = FALSE,
    sep = "\t",
    row.names = FALSE,
    col.names = FALSE)



### END ------------------------------------------------------------------

# Time scripts end
print(paste("Completed ", format(Sys.time(), "%Y-%m-%d"), sep = ""))


