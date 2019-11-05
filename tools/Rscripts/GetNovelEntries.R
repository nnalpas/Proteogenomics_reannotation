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
            "R/6frames_proteogenomics/helper.R",
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
load_package("plyr")
load_package("dplyr")
load_package("seqinr")
load_package("bit64")
load_package("magrittr")
load_package("data.table")
load_package("splitstackshape")
load_package("stringr")
load_package("optparse")
load_package("ggplot2")
load_package("gtable")
load_package("grid")
load_package("gridExtra")
load_package("purrr")



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
            multi = FALSE))
    
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
            help = "ORF fasta file", metavar = "character"))
    
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

# Create output directory if not already existing
dir.create(opt$output)



### Data import ----------------------------------------------------------

# Import the maxquant evidence table
evid <- mq_read(
    path = opt$maxquant,
    name = "evidence.txt",
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

# Format association of peptide to protein
data <- evid %>%
    dplyr::select(., Sequence, Proteins) %>%
    cSplit(indt = ., splitCols = "Proteins", sep = ";", direction = "long") %>%
    unique(.) %>%
    base::as.data.frame(., stringsAsFactors = FALSE)

# Locate all the peptide within associated proteins
pep_loc <- list()
pep_loc[["Known"]] <- pept_locate(
    data = data, peptide = "Sequence",
    proteins = "Proteins", fasta = fasta$Known) %>%
    dplyr::filter(., !is.na(start))
pep_loc[["Novel"]] <- pept_locate(
    data = data, peptide = "Sequence",
    proteins = "Proteins", fasta = fasta$Novel) %>%
    dplyr::filter(., !is.na(start))
pep_loc[["Contaminant"]] <- data %>%
    dplyr::filter(., grepl("^CON__", Proteins)) %>%
    set_colnames(c("pep", "prot")) %>%
    dplyr::mutate(., start = NA_integer_, end = NA_integer_)
pep_loc[["Reverse"]] <- evid %>%
    dplyr::filter(., Reverse == "+") %>%
    dplyr::select(., Sequence, Proteins) %>%
    set_colnames(c("pep", "prot")) %>%
    dplyr::mutate(., start = NA_integer_, end = NA_integer_)

# New dataframe to hold info about fasta of origin for each sequence
evid_match <- plyr::ldply(
    .data = pep_loc, .fun = "data.frame", .id = "DatabID") %>%
    dplyr::group_by(., pep) %>%
    dplyr::summarise(
        ., 
        group = ifelse(
            test = any(DatabID == "Known"),
            yes = "Known",
            no = ifelse(
                test = any(DatabID == "Contaminant"),
                yes = "Contaminant",
                no = ifelse(
                    test = any(DatabID == "Novel"),
                    yes = "Novel",
                    no = ifelse(
                        test = any(DatabID == "Reverse"),
                        yes = "Reverse",
                        no = NA_character_)))),
        Database = ifelse(
            test = group %in% c("Known", "Contaminant"),
            yes = "Target",
            no = ifelse(
                test = group == "Reverse",
                yes = "Decoy",
                no = ifelse(
                    test = group == "Novel",
                    yes = "Novel",
                    no = NA_character_)))) %>%
    dplyr::left_join(x = evid, y = ., by = c("Sequence" = "pep")) %>%
	dplyr::mutate(
		., OnlyIdBySite = ifelse(id %in% pg_not_sites, TRUE, FALSE))

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

# Histogram of evidence counts
toplot <- evid_match %>%
    dplyr::group_by(., group) %>%
    dplyr::summarise(., evid_count = n())
pl <- plots_hist(
    data = toplot,
    key = "group",
    value = "evid_count",
    group = "group",
    fill = "group",
    main = "PSM count",
    xlabel = "Databases",
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
    dplyr::select(., group, `Mass Error [ppm]`) %>%
    set_colnames(c("group", "Masserror"))
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
    dplyr::select(., Database, `Mass Error [ppm]`) %>%
    set_colnames(c("Database", "Masserror"))
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

# Bowplot of evidence PEP
toplot <- evid_match %>%
    dplyr::select(., Database, PEP)
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

# Save the group mapping data
saveRDS(
    object = evid_match,
    file = paste(
        opt$output, "/", "Sequence_group_mapping.RDS", sep = ""))

# Save the peptide location data
saveRDS(
    object = pep_loc,
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

