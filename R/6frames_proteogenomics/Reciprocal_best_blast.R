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



### Parameters setting up ------------------------------------------------

# Define the list of command line parameters
option_list <- list(
    make_option(
        opt_str = c("-b", "--blast"),
        type = "character", default = NULL, 
        help = "Blast data file name", metavar = "character"),
    make_option(
        opt_str = c("-r", "--reciprocal_blast"),
        type = "character", default = NULL, 
        help = "Reciprocal blast data file name", metavar = "character"),
    make_option(
        opt_str = c("-o", "--output"),
        type = "character", default = "reciprocal_best_hits.txt", 
        help = "Output file name [default= %default]", metavar = "character"))

# Parse the parameters provided on command line by user
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# Check whether inputs parameter was provided
if (is.null(opt$blast) | is.null(opt$reciprocal_blast)){
    
    print_help(opt_parser)
    stop(paste(
        "The two input arguments must be supplied",
        "(blast and reciprocal blast files)!"))
    
}

# Check whether output parameter was provided
if (is.null(opt$output)){
    
    opt$output <- paste(dirname(opt$input), opt$output, sep = "/")
    warning("Output results to reciprocal_best_hits.txt!")
    
}



### Data import and best blast computation -------------------------------

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

# Get the best blast match for each query
best_blast_data <- best_blast(data = blast_data, key = "qseqid")

# Import the reciprocal Blast results
reciproc_data <- read.table(
    file = opt$reciprocal_blast,
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

# Get the best reciprocal blast match for each query
best_reciproc_data <- best_blast(data = reciproc_data, key = "qseqid")



### Reciprocal best blast hits identification ----------------------------

# Merge the best blast and best reciprocal data
blast_merge <- best_blast_data %>%
    dplyr::full_join(
        x = .,
        y = best_reciproc_data,
        by = c("qseqid" = "sseqid"),
        suffix = c("_blast", "_reciproc"))

# Determine which entry have a reciprocal best hits




