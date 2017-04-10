#!/usr/bin/env Rscript

# This script compiles blast and reciprocal blast hits from blasting results
# to identify reciprocal best hits and matching position



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
        "/home-link",
        user,
        "bin/General_function.R",
        sep = "/"))

# Load the required packages (or install if not already in library)
load_package(plyr)
load_package(dplyr)
load_package(magrittr)
load_package(data.table)
load_package(splitstackshape)
load_package(stringr)
load_package(optparse)
load_package(UniProt.ws)
load_package(seqinr)



### Parameters setting up ------------------------------------------------

# Define the list of command line parameters
option_list <- list(
    make_option(
        opt_str = c("-f", "--fasta"),
        type = "character", default = NULL,
        help = "Fasta file", metavar = "character"),
    make_option(
        opt_str = c("-b", "--blast"),
        type = "character", default = NULL,
        help = "Blast file", metavar = "character"),
    make_option(
        opt_str = c("-g", "--genomic_position"),
        type = "character", default = NULL,
        help = "Genomic position file", metavar = "character"),
    make_option(
        opt_str = c("-o", "--output"),
        type = "character", default = "protein_location.txt", 
        help = "Output file name [default= %default]", metavar = "character"))

# Parse the parameters provided on command line by user
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# Check whether inputs parameter was provided
if (is.null(opt$fasta) | is.null(opt$blast) | is.null(opt$genomic_position)) {
    
    print_help(opt_parser)
    stop(paste(
        "The three input arguments must be supplied",
        "(fasta, blast and genomic position files)!"))
    
}

# Check whether output parameter was provided
if (opt$output == "protein_location.txt") {
    
    warning("Output results to protein_location.txt!")
    
}



### Data import ----------------------------------------------------------

# Import the fasta file
fasta <- opt$fasta %>%
    as.character(.) %>%
    seqinr::read.fasta(
        file = ., seqtype = "AA", as.string = TRUE) %>%
    set_names(uni_id_clean(names(.)))
    
# Get all protein IDs
uni_ids <- names(fasta) %>%
    as.character(.) %>%
    unique(.)

# Import the Blast results
blast_data <- read.table(
    file = opt$blast,
    header = FALSE,
    sep = "\t",
    quote = "",
    col.names = c(
        "qseqid", "sseqid", "pident", "nident", "mismatch", "length",
        "gapopen", "qstart", "qend", "sstart", "send", "evalue",
        "bitscore", "score"),
    as.is = TRUE) %>%
    dplyr::mutate(
        .,
        qseqid = uni_id_clean(qseqid),
        sseqid = uni_id_clean(sseqid))

# Import the genomic position data
blast_data <- read.table(
    file = opt$genomic_position, header = TRUE, sep = "\t", quote = "")



### Get target ID position based on blast --------------------------------

#




