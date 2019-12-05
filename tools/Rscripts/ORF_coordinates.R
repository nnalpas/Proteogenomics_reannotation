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
library(tidyr)
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
        fasta = choose.files(
            caption = "Choose fasta file containing ORF proteins!",
            multi = FALSE),
        output = readline(
            prompt = "Define the output directory!"))
    
} else {
    
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
    
}

# Check whether inputs parameter was provided
if (is.null(opt$fasta)) {
    
    print_help(opt_parser)
    stop("The fasta input argument must be supplied!")
    
}

# Check whether output parameter was provided
if (opt$output == "orf_coordinates.txt") {
    
    warning("Output results to orf_coordinates.txt!")
    
}

# Create output directory if not already existing
dir.create(dirname(opt$output))



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
    test = orf_pos$start < orf_pos$end, yes = "+", no = "-")
orf_pos[orf_pos$strand == "+", "frame"] <- ((as.numeric(
    orf_pos[orf_pos$strand == "+", "start"]) + 2) / 3) %>%
    keep_decimals(.) %>%
    round(x = ((. + 0.33) * 3))
orf_pos[orf_pos$strand == "-", "frame"] <- ((as.numeric(
    orf_pos[orf_pos$strand == "-", "end"]) + 2) / 3) %>%
    keep_decimals(.) %>%
    round(x = ((. + 0.33) * -3))

# Define the chromosome of origin
orf_pos$chromosome <- sub("_.+", "", orf_pos$id)

# Export the ORF coordinates data
write.table(
    x = orf_pos, file = opt$output,
    quote = FALSE, sep = "\t",
    row.names = FALSE, col.names = TRUE)

# Define end time
print(paste("Complete", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))


