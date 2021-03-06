#!/usr/bin/env Rscript

# This script generates a BSgenome package



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
#if (interactive()) {
#    source(
#        file = paste(
#            "C:/Users",
#            user,
#            "Documents/GitHub/Proteogenomics_reannotation/",
#            "R/6frames_proteogenomics/helper.R",
#            sep = "/"))
#} else {
#    source(
#        file = paste(
#            Sys.getenv("HOME"),
#            "bin/helper.R",
#            sep = "/"))
#}

# Load the required packages (or install if not already in library)
library("plyr")
library("dplyr")
library("magrittr")
library("data.table")
library("splitstackshape")
library("stringr")
library("optparse")
library("seqinr")
library("BiocGenerics")
library("S4Vectors")
library("IRanges")
library("GenomeInfoDb")
library("GenomicRanges")
library("Biostrings")
library("rtracklayer")
library("BSgenome")
#library("devtools")



### Parameters setting up ------------------------------------------------

# Define input parameters (interactively or from command line)
if (interactive()) {
    
    # Define the list of input parameters
    opt <- list(
        seed = choose.files(
            caption = "Choose the seed file!",
            multi = FALSE),
        fasta = choose.dir(
            caption = "Choose the directory containing the fasta files!"),
        output = choose.dir(
            caption = "Choose the output directory!"))
    
} else {
    
    # Define the list of command line parameters
    option_list <- list(
        make_option(
            opt_str = c("-s", "--seed"),
            type = "character", default = NULL,
            help = "BSgenome seed file",
            metavar = "character"),
        make_option(
            opt_str = c("-f", "--fasta"),
            type = "character", default = NULL,
            help = "Fasta sequences directory",
            metavar = "character"),
        make_option(
            opt_str = c("-o", "--output"),
            type = "character", default = NULL,
            help = "Output directory",
            metavar = "character"))
    
    # Parse the parameters provided on command line by user
    opt_parser <- OptionParser(option_list = option_list)
    opt <- parse_args(opt_parser)
    
}

# Check whether inputs parameter were provided
if (is.null(opt$seed) | opt$seed == "") {
    
    print_help(opt_parser)
    stop("BSgenome seed file required!")
    
}
if (is.null(opt$fasta) | opt$fasta == "") {
    
    print_help(opt_parser)
    stop("Fasta sequence directory required!")
    
}
if (is.null(opt$output) | opt$output == "") {
    
    print_help(opt_parser)
    stop("Output directory required!")
    
}



### Forge target package -------------------------------------------------

# Make a BSgenome data package using the seed file
forgeBSgenomeDataPkg(
    x = opt$seed, seqs_srcdir = opt$fasta, destdir = opt$output)

# Get the name of the package directory
pkg_name <- readLines(con = opt$seed, n = 1) %>%
    sub(pattern = "^Package: ", replacement = "", x = .)

# Build the source package
orig_dir <- getwd()
setwd(opt$output)
print("Building package...")
system(
    command = paste(
        "R CMD build --no-manual --no-build-vignettes ",
        shQuote(pkg_name)),
    wait = TRUE)
setwd(orig_dir)
#devtools::build(
#    pkg = paste0(opt$output, "/", pkg_name),
#    path = opt$output, binary = FALSE)

# Check the package
print("Checking package...")
system(
    command = paste(
        "R CMD check --no-manual --no-build-vignettes ",
        shQuote(list.files(
            path = opt$output, pattern = paste0("^", pkg_name, ".+tar.gz$"), full.names = TRUE))),
    wait = TRUE)

# Install the package
print("Installing package...")
system(
    command = paste(
        "R CMD INSTALL ",
        shQuote(list.files(
            path = opt$output, pattern = paste0("^", pkg_name, ".+tar.gz$"), full.names = TRUE))),
    wait = TRUE)

# Define end time
print(paste("Complete", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))


