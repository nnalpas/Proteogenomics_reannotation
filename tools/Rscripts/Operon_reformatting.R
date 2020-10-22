#!/usr/bin/env Rscript

# This script reformats operon data to be suitable for genomic range creation
# This script is specific for subtiwiki provided operon files



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
library(bit64)
library(ggplot2)
library(gtable)
library(grid)
library(gridExtra)
library(purrr)



### Parameters setting up ------------------------------------------------

# Define input parameters (interactively or from command line)
if (interactive()) {
    
    # Define the list of input parameters
    opt <- list(
        operon = choose.files(
            caption = "Choose operon file!",
            multi = FALSE),
        output = readline(
            prompt = "Define the output file!"))
    
} else {
    
    # Define the list of command line parameters
    option_list <- list(
        make_option(
            opt_str = c("-i", "--operon"),
            type = "character", default = "",
            help = "Operon file (SubtiWiki)",
            metavar = "character"),
        make_option(
            opt_str = c("-o", "--output"),
            type = "character", default = "",
            help = "Output file", metavar = "character"))
    
    # Parse the parameters provided on command line by user
    opt_parser <- OptionParser(option_list = option_list)
    opt <- parse_args(opt_parser)
    
}

# Check whether inputs parameter was provided
if (
    identical(opt$operon, NULL) |
    identical(opt$operon, "") |
    identical(opt$operon, character(0))) {
    
    print_help(opt_parser)
    stop("The input operon must be supplied!")
    
}

# Check whether output parameter was provided
if (
    identical(opt$output, NULL) |
    identical(opt$output, "") |
    identical(opt$output, character(0))) {
    
    print_help(opt_parser)
    stop("The output filename must be supplied!")
    
}

# Create output directory if not already existing
dir.create(dirname(opt$output))



### Data import, reformatting and export ---------------------------------

# Import the input data
operon <- data.table::fread(
    input = opt$operon, sep = "\t", header = TRUE,
    stringsAsFactors = FALSE, data.table = FALSE)

# Format the data to be suitable for genomic ranges object creation
operon_format <- operon %>%
    dplyr::group_by(., id) %>%
    dplyr::summarise(
        .,
        start = min(c(Start, End)),
        end = max(c(Start, End)),
        strand = unique(Strand),
        length = (end - start)+1,
        chromosome = Chromosome) %>%
    set_colnames(c(
        "id", "start", "end",
        "strand", "length", "chromosome"))
operon_format <- operon %>%
    dplyr::group_by(., id) %>%
    summarise_if(., is.character, toString, width = NULL) %>%
    dplyr::left_join(x = operon_format, y = ., by = "id")

# Export the reformatted data
write.table(
    x = operon_format, file = opt$output, quote = FALSE,
    sep = "\t", row.names = FALSE, col.names = TRUE)



### END ------------------------------------------------------------------

# Time scripts end
print(paste("Completed ", format(Sys.time(), "%Y-%m-%d"), sep = ""))


