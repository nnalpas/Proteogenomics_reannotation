#!/usr/bin/env Rscript

# This script determines the ORF genomic coordinates based on fasta header



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
#source(
#    file = paste(
#        "C:/Users",
#        user,
#        "Documents/GitHub/Miscellaneous/R/General/General_function.R",
#        sep = "/"))
source(
    file = paste(
        Sys.getenv("HOME"),
        "bin/General_function.R",
        sep = "/"))

# Load the required packages (or install if not already in library)
load_package("plyr")
load_package("dplyr")
load_package("tidyr")
load_package("magrittr")
load_package("data.table")
load_package("splitstackshape")
load_package("stringr")
load_package("optparse")
load_package("seqinr")



### Parameters setting up ------------------------------------------------

# Define the list of command line parameters
option_list <- list(
    make_option(
        opt_str = c("-f", "--fasta"),
        type = "character", default = NULL,
        help = "Fasta file", metavar = "character"),
    make_option(
        opt_str = c("-o", "--output"),
        type = "character", default = "orf_coordinates.txt", 
        help = "Output file name [default= %default]", metavar = "character"))

# Parse the parameters provided on command line by user
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# Check whether inputs parameter was provided
if (is.null(opt$fasta)) {
    
    print_help(opt_parser)
    stop("The fasta input argument must be supplied!")
    
}

# Check whether output parameter was provided
if (opt$output == "orf_coordinates.txt") {
    
    warning("Output results to orf_coordinates.txt!")
    
}

# For manual parameters set-up
#opt <- list(
#    fasta = "C:/Users/kxmna01/Dropbox/Home_work_sync/Work/Colleagues shared work/Vaishnavi Ravikumar/Bacillus_subtilis_6frame/Databases/Find0_GCA_000009045.1_FIXED.fasta",
#    output = "C:/Users/kxmna01/Dropbox/Home_work_sync/Work/Colleagues shared work/Vaishnavi Ravikumar/Bacillus_subtilis_6frame/Databases/Find0_GCA_000009045.1_FIXED_orf_coordinates.txt")



### Data import ----------------------------------------------------------

# Import the fasta file
fasta <- opt$fasta %>%
    as.character(.) %>%
    seqinr::read.fasta(
        file = ., seqtype = "AA", as.string = TRUE)



### Get genomic coordinates ----------------------------------------------

# Get the start-end genomic position for each ORFs
orf_pos <- lapply(
    X = fasta, FUN = function(x) { attr(x, "Annot") }) %>%
    unlist(.) %>%
    sub("^.*? \\[(.*?)\\] .*$", "\\1", .) %>%
    base::data.frame(id = names(.), pos = ., stringsAsFactors = FALSE) %>%
    tidyr::separate(
        data = ., col = pos, into = c("start", "end"), sep = " - ") %>%
    dplyr::mutate(., start = as.numeric(start), end = as.numeric(end))

# Filter out the ORF in the fasta that are not holding codons (multiple of 3)
orf_pos %<>%
    dplyr::mutate(
        ., lengthCheck = keep_decimals(x = (abs(x = (start - end)) + 1) / 3)) %>%
    dplyr::filter(., lengthCheck == 0) %>%
    dplyr::select(., -lengthCheck) %>%
    base::as.data.frame(., stringsAsFactors = FALSE)

# Determine nucleotide and aa length
orf_pos %<>%
    dplyr::mutate(
        .,
        nucl_length = (abs(x = (end - start)) + 1),
        aa_length = nucl_length / 3)

# Determine strand and frame of ORF
orf_pos$strand <- ifelse(
    test = orf_pos$start < orf_pos$end, yes = 1, no = -1)
orf_pos[orf_pos$strand == 1, "frame"] <- ((as.numeric(
    orf_pos[orf_pos$strand == 1, "start"]) + 2) / 3) %>%
    keep_decimals(.) %>%
    round(x = ((. + 0.33) * 3))
orf_pos[orf_pos$strand == -1, "frame"] <- ((as.numeric(
    orf_pos[orf_pos$strand == -1, "end"]) + 2) / 3) %>%
    keep_decimals(.) %>%
    round(x = ((. + 0.33) * -3))

# Export the ORF coordinates data
write.table(
    x = orf_pos, file = opt$output,
    quote = FALSE, sep = "\t",
    row.names = FALSE, col.names = TRUE)

# Define end time
print(paste("Complete", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))


