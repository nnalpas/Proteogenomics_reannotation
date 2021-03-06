#!/usr/bin/env Rscript

# This script gets genomic coordinates using MaxQuant data for each entry, 
# such as peptide or modified site, that can be located within protein



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
library(plyr)
library(dplyr)
library(magrittr)
library(data.table)
library(splitstackshape)
library(stringr)
library(optparse)
library(seqinr)



### Parameters setting up ------------------------------------------------

# Define input parameters (interactively or from command line)
if (interactive()) {
    
    # Define the list of input parameters
    opt <- list(
        peptide = choose.files(
            caption = "Choose the peptide location file!",
            multi = FALSE),
        novel_reason = choose.files(
            caption = "Choose peptide novelty (.RDS) file!",
            multi = FALSE),
        known_prot = choose.files(
            caption = "Choose the known protein genomic coordinates file!",
            multi = FALSE),
        orf_prot = choose.files(
            caption = "Choose the ORF protein genomic coordinates file!",
            multi = FALSE),
        output = readline(
            prompt = "Define the output folder!"))
    
} else {
    
    # Define the list of command line parameters
    option_list <- list(
        make_option(
            opt_str = c("-p", "--peptide"),
            type = "character", default = NULL,
            help = "Peptide location file",
            metavar = "character"),
        make_option(
            opt_str = c("-r", "--novel_reason"),
            type = "character", default = NULL,
            help = "Peptide novelty info",
            metavar = "character"),
        make_option(
            opt_str = c("-k", "--known_prot"),
            type = "character", default = NULL,
            help = "Known protein genomic coordinates file",
            metavar = "character"),
        make_option(
            opt_str = c("-n", "--orf_prot"),
            type = "character", default = NULL,
            help = "ORF protein genomic coordinates file",
            metavar = "character"),
        make_option(
            opt_str = c("-o", "--output"),
            type = "character", default = "PeptPosition",
            help = "Output folder [default= %default]",
            metavar = "character"))
    
    # Parse the parameters provided on command line by user
    opt_parser <- OptionParser(option_list = option_list)
    opt <- parse_args(opt_parser)
    
}

# Check whether inputs parameter were provided
if (
    identical(opt$peptide, NULL) |
    identical(opt$peptide, "") |
    identical(opt$peptide, character(0))) {
    
    print_help(opt_parser)
    stop(paste(
        "Peptide location within proteins required!"))
    
}
if (
    identical(opt$novel_reason, NULL) |
    identical(opt$novel_reason, "") |
    identical(opt$novel_reason, character(0))) {
    
    print_help(opt_parser)
    stop("The input peptide novelty must be supplied!")
    
}
if (
    identical(opt$known_prot, NULL) |
    identical(opt$known_prot, "") |
    identical(opt$known_prot, character(0))) {
    
    print_help(opt_parser)
    stop(paste(
        "Known protein genomic coordinates file is required!"))
    
}
if (
    identical(opt$orf_prot, NULL) |
    identical(opt$orf_prot, "") |
    identical(opt$orf_prot, character(0))) {
    
    print_help(opt_parser)
    stop(paste(
        "ORF protein genomic coordinates file is required!"))
    
}

# Check whether output parameter was provided
if (identical(opt$output, NULL) |
    identical(opt$output, "") |
    identical(opt$output, character(0))) {
    
    opt$output <- "./PeptPosition"
    warning(paste0(
        "Output results to '",
        opt$output,
        "'!"))
    
}

# Create output directory if not already existing
dir.create(opt$output)



### Data import ----------------------------------------------------------

# Import the peptide location
pep_loc <- readRDS(file = opt$peptide)

# Import the peptide novelty explanation
evid_reason <- readRDS(file = opt$novel_reason)

# Import the known protein genomic coordinates
known_prot <- data.table::fread(
    input = opt$known_prot, sep = "\t", header = TRUE,
    stringsAsFactors = FALSE, quote = "")

# Import the ORF protein genomic coordinates
orf_prot <- data.table::fread(
    input = opt$orf_prot, sep = "\t", header = TRUE,
    stringsAsFactors = FALSE, quote = "")



### Compute genomic coordinates per peptides -----------------------------

# Check how many peptides were not located (NA values)
# and therefore will be filtered out
print(paste(
    "There are", nrow(pep_loc[is.na(pep_loc$start), ]), "peptides",
    "with no location value, these will be filtered out!"))

# Filter out all peptides that do not have a location within their associated
# proteins
pep_loc %<>%
    dplyr::filter(., !is.na(start))# %>%
    #dplyr::select(., -Database)

# Compile genomic coordinates of known and ORF proteins
coordinates <- dplyr::bind_rows(known_prot, orf_prot)

# Calculate nucleotide location within protein for each peptide
pep_loc %<>%
    dplyr::mutate(
        .,
        start_nucl = as.integer(start) * 3 - 2,
        end_nucl = as.integer(end) * 3) %>%
    dplyr::rename(
        .,
        startAA = start,
        endAA = end)

# Compute the genomic coordinates for each peptide
pep_coord <- coordinates %>%
    dplyr::rename(., start_prot = start, end_prot = end) %>%
    dplyr::left_join(
        x = pep_loc, y = .,
        by = c("Proteins" = "id")) %>%
    dplyr::mutate(
        .,
        start = ifelse(
            strand == "+",
            start_prot + start_nucl - 1,
            start_prot - start_nucl + 1),
        end = ifelse(
            strand == "+",
            start_prot + end_nucl - 1,
            start_prot - end_nucl + 1),
        nucl_length = nchar(Sequence) * 3,
        aa_length = nchar(Sequence))

# Keep only required columns and remove peptides without coordinates
pep_coord %<>%
    dplyr::select(
        ., -start_nucl, -end_nucl, -start_prot, -end_prot) %>%
    dplyr::filter(., !is.na(start)) %>%
    unique(.)

# Concatenate peptide with same id and coordinates (duplicate)
# we are now looking at genomic position so this is independent of protein name
pep_coord %<>%
    dplyr::group_by(., id, start, end) %>%
    dplyr::summarise_all(funs(toString(x = unique(.), width = NULL))) %>%
    dplyr::ungroup(.)

# Include for each peptide its novelty reason providing 'Sequence' exists
if ("Sequence" %in% colnames(pep_coord)) {
    pep_coord <- evid_reason %>%
        dplyr::select(., Sequence, NoveltyReason) %>%
        unique(.) %>%
        dplyr::left_join(
            x = pep_coord,
            y = .,
            by = "Sequence")
}



### Results export -------------------------------------------------------

# Export the peptide entries coordinates data
write.table(
    x = pep_coord,
    file = paste0(
        opt$output, "/",
        sub("_location.RDS$", "_coordinates.txt", basename(opt$peptide))),
    quote = FALSE, sep = "\t",
    row.names = FALSE, col.names = TRUE)

# Define end time
print(paste("Complete", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))


