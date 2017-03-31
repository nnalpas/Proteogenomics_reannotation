#!/usr/bin/Rscript

# This script identifies the best blast hits from blasting results
# for target entries



### Environment set-up ---------------------------------------------------

# Start with clean environment
rm(list = ls())

# Define current time
date_str <- format(Sys.time(), "%Y-%m-%d")
print(paste("Start ", date.str, sep = ""))

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
        "/home",
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
        opt_str = c("-i", "--input"), type = "character", default = NULL, 
        help = "Blast data file name", metavar = "character"),
    make_option(
        opt_str = c("-o", "--output"), type = "character", default = "out.txt", 
        help = "Output file name [default= %default]", metavar = "character"))

# Parse the parameters provided on command line by user
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# Check whether input parameter was provided
if (is.null(opt$input)){
    
    print_help(opt_parser)
    stop("At least one argument must be supplied (input file)!")
    
}

# Check whether output parameter was provided
if (is.null(opt$output)){
    
    opt$output <- paste(dirname(opt$input), opt$output, sep = "/")
    warning("Output results to out.txt!")
    
}



### Data import and processing -------------------------------------------

# Import the Blast results
blast_data <- read.table(
    file = opt$input,
    header = FALSE,
    sep = "\t",
    quote = "",
    col.names = c(
        "qseqid", "sseqid", "pident", "nident", "mismatch", "length",
        "gapopen", "qstart", "qend", "sstart", "send", "evalue",
        "bitscore", "score"),
    as.is = TRUE)

# Get the best blast match for each query
best_blast_data <- best_blast(data = blast_data, key = "qseqid")

# Export the best hits protein IDs that needs to be reciprocally blasted
write.table(
    x = unique(best_blast_data$sseqid),
    file = opt$output,
    quote = FALSE,
    sep = "\t",
    row.names = FALSE,
    col.names = FALSE)


