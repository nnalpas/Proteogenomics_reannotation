#!/usr/bin/env Rscript

# This script determines the novelty reasons for each ORF based on
# genomic coordinates, blast results...



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
            "/home-link",
            user,
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
load_package("bit64")
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
        novel_reason = choose.files(
            caption = "Choose peptide novelty (.RDS) file!",
            multi = FALSE),
        ref_grange = choose.files(
            caption = "Choose reference entries grange (.RDS) file!",
            multi = FALSE),
        orf_grange = choose.files(
            caption = "Choose ORF entries grange (.RDS) file!",
            multi = FALSE),
        output = readline(
            prompt = "Define the output directory!"))
    
} else {
    
    # Define the list of command line parameters
    option_list <- list(
        make_option(
            opt_str = c("-i", "--novel_reason"),
            type = "character", default = NULL,
            help = "Peptide novelty info",
            metavar = "character"),
        make_option(
            opt_str = c("-r", "--ref_grange"),
            type = "character", default = NULL,
            help = "Reference entries grange file",
            metavar = "character"),
        make_option(
            opt_str = c("-n", "--orf_grange"),
            type = "character", default = NULL,
            help = "ORF entries grange file",
            metavar = "character"),
        make_option(
            opt_str = c("-o", "--output"),
            type = "character", default = "",
            help = "Output directory", metavar = "character"))
    
    # Parse the parameters provided on command line by user
    opt_parser <- OptionParser(option_list = option_list)
    opt <- parse_args(opt_parser)
    
}

# Check whether inputs parameter was provided
if (is.null(opt$novel_reason)) {
    
    print_help(opt_parser)
    stop("The input peptide novelty must be supplied!")
    
}
if (is.null(opt$ref_grange)) {
    
    print_help(opt_parser)
    stop("The input reference entries grange must be supplied!")
    
}
if (is.null(opt$orf_grange)) {
    
    print_help(opt_parser)
    stop("The input ORF entries grange must be supplied!")
    
}

# Check whether output parameter was provided
if (opt$output == "") {
    
    opt$output <- dirname(opt$novel_reason)
    warning(paste("Output results to ", opt$output, "!"))
    
}

# Create output directory if not already existing
dir.create(opt$output)



### Data import ----------------------------------------------------------

# Import the novelty reason file containing peptide novelty info
evid_reason <- opt$novel_reason %>%
    as.character(.) %>%
    readRDS(file = .)

# Import the reference genomic ranges file
ref_grange <- opt$ref_grange %>%
    as.character(.) %>%
    readRDS(file = .)

# Import the ORF genomic ranges file
orf_grange <- opt$orf_grange %>%
    as.character(.) %>%
    readRDS(file = .)



### END ------------------------------------------------------------------

# Time scripts end
print(paste("Completed ", format(Sys.time(), "%Y-%m-%d"), sep = ""))


