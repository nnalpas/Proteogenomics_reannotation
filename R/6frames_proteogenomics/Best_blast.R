#!/usr/bin/env Rscript

# This script identifies the best blast hits from blasting results
# for target entries



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
        opt_str = c("-i", "--input"), type = "character", default = NULL, 
        help = "Blast data file name", metavar = "character"),
    make_option(
        opt_str = c("-f", "--filter"), type = "character", default = NULL, 
        help = "Specific filtering to apply", metavar = "character"),
    make_option(
        opt_str = c("-m", "--multi_match"), type = "character", default = NULL, 
        help = "Filter for multi hits entry", metavar = "character"),
    make_option(
        opt_str = c("-o", "--out_path"), type = "character", default = NULL, 
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
if (is.null(opt$out_path)){
    
    opt$output <- dirname(opt$input)
    warning(paste("Output results to path: ", opt$output, "!", sep = ""))
    
}

# If filter and multi_match parameters are undefined, define as null
if (is.null(opt$filter)) {
    
    opt$filter <- NULL
    
}
if (is.null(opt$multi_match)) {
    
    opt$multi_match <- NULL
    
}

# For manual parameters set-up
#opt <- list(
#    input = "C:/Users/kxmna01/Dropbox/Home_work_sync/Work/Colleagues shared work/Vaishnavi Ravikumar/Bacillus_subtilis_6frame/blastp_ORF_Nicolas_vs_Karsten_11012017",
#    filter = "pident == 100 & nident == length & qstart == sstart & qend == send",
#    multi_match = "uniquify",
#    output = "C:/Users/kxmna01/Dropbox/Home_work_sync/Work/Colleagues shared work/Vaishnavi Ravikumar/Bacillus_subtilis_6frame/Databases")



### Data import and processing -------------------------------------------

# Import the Blast results
blast_data <- blast_read(file = opt$input, blast_format = "6")

# Get the best blast match for each query
best_blast_data <- best_blast(
    data = blast_data, key = "qseqid",
    filter = opt$filter, multi_match = opt$multi_match)

# Export the best hits protein IDs that needs to be reciprocally blasted
write.table(
    x = unique(best_blast_data$sseqid),
    file = paste(opt$output, "/Reciprocal_id_", basename(opt$input), sep = ""),
    quote = FALSE,
    sep = "\t",
    row.names = FALSE,
    col.names = FALSE)

# Export the best hits results
write.table(
    x = best_blast_data,
    file = paste(opt$output, "/Best_blast_", basename(opt$input), sep = ""),
    quote = FALSE,
    sep = "\t",
    row.names = FALSE,
    col.names = FALSE)

# Export the blast ID map results
write.table(
    x = best_blast_data %>%
        dplyr::select(., qseqid, sseqid),
    file = paste(opt$output, "/Blast_cross-map_", basename(opt$input), sep = ""),
    quote = FALSE,
    sep = "\t",
    row.names = FALSE,
    col.names = FALSE)

# Define end time
print(paste("Complete", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))


