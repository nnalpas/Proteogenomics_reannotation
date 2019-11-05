#!/usr/bin/env Rscript

# This script gets genomic coordinates for each peptide located within protein



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
load_package("plyr")
load_package("dplyr")
load_package("magrittr")
load_package("data.table")
load_package("splitstackshape")
load_package("stringr")
load_package("optparse")
load_package("seqinr")



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
        output = "peptide_location.txt")
    
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
            type = "character", default = "peptide_location.txt",
            help = "Output filename [default= %default]",
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
if (opt$output == "peptide_location.txt") {
    
    opt$output <- paste0("./", "peptide_location.txt")
    warning(paste0(
        "Output results to '",
        opt$output,
        "'!"))
    
}

# Create output directory if not already existing
dir.create(dirname(opt$output))



### Data import ----------------------------------------------------------

# Import the peptide location
pep_loc <- readRDS(file = opt$peptide) %>%
    plyr::ldply(.data = ., .fun = "data.frame", .id = "Database")

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
    dplyr::filter(., !is.na(start))

# Compile genomic coordinates of known and ORF proteins
coordinates <- dplyr::bind_rows(known_prot, orf_prot)

# Calculate nucleotide location within protein for each peptide
pep_loc %<>%
    dplyr::mutate(
        .,
        startAA = start,
        endAA = end) %>%
    dplyr::mutate(
        .,
        start_nucl = as.integer(start) * 3 - 2,
        end_nucl = as.integer(end) * 3)

# Compute the genomic coordinates for each peptide
pep_coord <- pep_loc %>%
    dplyr::left_join(
        x = ., y = coordinates,
        by = c("prot" = "id"), suffix = c("_pep", "_prot")) %>%
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
        nucl_length = nchar(pep) * 3,
        aa_length = nchar(pep))

# Keep only required columns and remove peptides without coordinates
pep_coord %<>%
    dplyr::select(
        ., pep, chromosome, strand, frame, start, end, nucl_length,
        aa_length, prot, startAA, endAA, ExactCoord, Comment, Database) %>%
    dplyr::filter(., !is.na(start)) %>%
    unique(.)

# Concatenate peptide with same sequence and coordinates (duplicate)
pep_coord %<>%
    dplyr::group_by(., pep, start, end) %>%
    dplyr::summarise_all(funs(toString(x = unique(.), width = NULL)))

# Rename the pep column to id (id is a mandatory column for GRange script)
colnames(pep_coord)[colnames(pep_coord) == "pep"] <- "id"

# Include for each peptide its novelty reason 
pep_coord %<>%
    dplyr::ungroup(.) %>%
    dplyr::left_join(
        x = .,
        y = evid_reason %>%
            dplyr::select(., Sequence, NoveltyReason) %>%
            unique(.),
        by = c("id" = "Sequence")) %>%
    dplyr::mutate(
        .,
        Database = ifelse(
            !is.na(NoveltyReason) & NoveltyReason == "Not novel",
            "Known", Database) %>%
            sub("^(Known).*", "\\1", .))



### Results export -------------------------------------------------------

# Export the peptide entries coordinates data
write.table(
    x = pep_coord, file = opt$output,
    quote = FALSE, sep = "\t",
    row.names = FALSE, col.names = TRUE)

# Define end time
print(paste("Complete", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))


