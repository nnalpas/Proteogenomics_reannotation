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
        chr_pattern = readline(
            prompt = paste(
                "Pattern to use to extract chromosome name per sequence?")),
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
            opt_str = c("-c", "--chr_pattern"), type = "character",
            default = NULL, help = "Pattern to extract chromosome name",
            metavar = "character"),
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

# Possible values:
# opt[["chr_pattern"]] <- "(.+?)_.+"
# opt[["chr_pattern"]] <- ".+?\\| (.+?) \\|.+"
if (is.null(opt$chr_pattern)) {
    warning("The default chromosome pattern argument is '(.+?)_.+'")
    opt[["chr_pattern"]] <- "(.+?)_.+"
} else if (opt$chr_pattern == "") {
    warning("The default chromosome pattern argument is '(.+?)_.+'")
    opt[["chr_pattern"]] <- "(.+?)_.+"
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

# Extract the sequence headers
fasta_annot <- lapply(
    X = fasta, FUN = function(x) { attr(x, "Annot") }) %>%
    unlist(.)

# Get the start-end genomic position and chromosome for each ORFs
orf_pos <- fasta_annot %>%
    base::data.frame(
        id = names(.),
        header = .,
        stringsAsFactors = FALSE) %>%
    dplyr::mutate(
        ., pos = sub("^.+ \\[(.*?)\\].*?$", "\\1", header),
        chromosome = sub(opt$chr_pat, "\\1", header)) %>%
    dplyr::select(., -header) %>%
    tidyr::separate(
        data = ., col = pos, into = c("start_tmp", "end_tmp"), sep = "(-|:)") %>%
    dplyr::mutate(
        ., start_tmp = as.numeric(sub(" *", "", start_tmp)),
        end_tmp = as.numeric(sub(" *", "", end_tmp)))

# Check if strand info is stored in header or in position
if (any(grepl("\\((+|-)\\)", fasta_annot))) {
    orf_pos$strand <- sub(".+\\((.)\\).*", "\\1", fasta_annot)
} else {
    orf_pos$strand <- ifelse(
        test = orf_pos$start_tmp < orf_pos$end_tmp, yes = "+", no = "-")
}

# Stop if the strand identification did not work
if (all(orf_pos$strand == "+") | all(orf_pos$strand == "-")) {
    stop(paste0(
        "Sequences have been identified solely on '+' or '-',",
        " there must be a mixed!"))
}

# Assign start and end positions based on strand info
orf_pos %<>%
    dplyr::rowwise(.) %>%
    dplyr::mutate(
        ., start = ifelse(
            strand == "+", min(start_tmp, end_tmp), max(start_tmp, end_tmp)),
        end = ifelse(
            strand == "+", max(start_tmp, end_tmp), min(start_tmp, end_tmp))) %>%
    dplyr::select(., -start_tmp, -end_tmp) %>%
    dplyr::ungroup(.)

# Determine if start/end match the sequence length and adapt based
# on 0- or 1-based coordinates and presence of stop codon
orf_pos %<>%
    dplyr::mutate(., nucl_length = abs(start - end))
orf_pos$nchar <- unlist(lapply(fasta, nchar))*3

diff_length <- unique(abs(orf_pos$nucl_length - orf_pos$nchar))
if (length(diff_length) == 1) {
    
    if (diff_length == 1) {
        orf_pos$nucl_length <- orf_pos$nucl_length+1
    } else if (diff_length == 2) {
        orf_pos$end <- ifelse(
            orf_pos$strand == "+", orf_pos$end-3, orf_pos$end+3)
        orf_pos$nucl_length <- (abs(orf_pos$start - orf_pos$end))+1
    } else if (diff_length == 0) {
        orf_pos$start <- ifelse(
            orf_pos$strand == "+", orf_pos$start+1, orf_pos$start-1)
        orf_pos$nucl_length <- (abs(orf_pos$start - orf_pos$end))+1
    } else if (diff_length == 3) {
        orf_pos$start <- ifelse(
            orf_pos$strand == "+", orf_pos$start+1, orf_pos$start-1)
        orf_pos$end <- ifelse(
            orf_pos$strand == "+", orf_pos$end-3, orf_pos$end+3)
        orf_pos$nucl_length <- (abs(orf_pos$start - orf_pos$end))+1
    } else {
        stop(paste0(
            "The difference between calculated length and sequence",
            "length is not allowed: ", diff_length))
    }
    
} else {
    stop(paste0(
        "There seem to me a mixture of coordinates based: ",
        paste0(unique(abs(orf_pos$nucl_length - orf_pos$nchar)), collapse = ", ")))
}

# Filter out the ORF in the fasta that are not holding codons (multiple of 3)
orf_pos %<>%
    dplyr::mutate(
        ., lengthCheck = keep_decimals(x = nucl_length / 3)) %>%
    dplyr::filter(., lengthCheck == 0) %>%
    dplyr::select(., -lengthCheck) %>%
    base::as.data.frame(., stringsAsFactors = FALSE)

# Determine aa length
orf_pos %<>%
    dplyr::mutate(., aa_length = nucl_length / 3)

# Determine frame of ORF
orf_pos[orf_pos$strand == "+", "frame"] <- ((as.numeric(
    orf_pos[orf_pos$strand == "+", "start"]) + 2) / 3) %>%
    keep_decimals(.) %>%
    round(x = ((. + 0.33) * 3))
orf_pos[orf_pos$strand == "-", "frame"] <- ((as.numeric(
    orf_pos[orf_pos$strand == "-", "end"]) + 2) / 3) %>%
    keep_decimals(.) %>%
    round(x = ((. + 0.33) * -3))

# Export the ORF coordinates data
write.table(
    x = orf_pos, file = opt$output,
    quote = FALSE, sep = "\t",
    row.names = FALSE, col.names = TRUE)

# Define end time
print(paste("Complete", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))


