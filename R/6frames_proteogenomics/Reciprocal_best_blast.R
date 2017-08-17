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
        type = "character", default = NULL, 
        help = "Output directory", metavar = "character"))

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
    
    opt$output <- dirname(opt$input)
    warning(paste("Output results to path: ", opt$output, "!", sep = ""))
    
}



### Data import and best blast computation -------------------------------

# Import the Blast results
blast_data <- blast_read(file = opt$blast, blast_format = "6") %>%
    dplyr::mutate(
        .,
        qseqid = uni_id_clean(qseqid),
        sseqid = uni_id_clean(sseqid))

# Get the best blast match for each query
best_blast_data <- best_blast(data = blast_data, key = "qseqid")

# Import the reciprocal Blast results
reciproc_data <- blast_read(
    file = opt$reciprocal_blast, blast_format = "6") %>%
    dplyr::mutate(
        .,
        qseqid = uni_id_clean(qseqid),
        sseqid = uni_id_clean(sseqid))

# Get the best reciprocal blast match for each query
best_reciproc_data <- best_blast(data = reciproc_data, key = "qseqid")



### Reciprocal best blast hits identification ----------------------------

# Merge the best blast and best reciprocal data
blast_merge <- best_blast_data %>%
    dplyr::left_join(
        x = .,
        y = best_reciproc_data,
        by = c("qseqid" = "sseqid"),
        suffix = c("_blast", "_reciproc"))

# Determine which entry have a reciprocal best hits
blast_merge_confirmed <- blast_merge %>%
    dplyr::filter(
        .,
        is.na(qseqid_reciproc) |
            sseqid == qseqid_reciproc) %>%
    dplyr::mutate(., best_reciprocal_type = case_when(
        is.na(.$qseqid_reciproc) ~ "None",
            .$best_count_reciproc == 1 ~ "Single",
        TRUE ~ "Multiple"))

# Export the best hits results
write.table(
    x = best_blast_data,
    file = paste0(
        opt$output, "/Best_Reciproc_Blast_", basename(opt$input)),
    quote = FALSE,
    sep = "\t",
    row.names = FALSE,
    col.names = TRUE)

# Export the blast ID map results
write.table(
    x = best_blast_data %>%
        dplyr::select(., qseqid, sseqid),
    file = paste0(
        opt$output, "/Reciproc_Blast_cross-map_", basename(opt$input)),
    quote = FALSE,
    sep = "\t",
    row.names = FALSE,
    col.names = FALSE)

# Define end time
print(paste("Complete", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))


