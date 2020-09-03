#!/usr/bin/env Rscript

# This script retrieve column annotations based on specific keys using
# uniprot.ws package



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
library(UniProt.ws)



### Parameters setting up ------------------------------------------------

# Define input parameters (interactively or from command line)
if (interactive()) {
    
    # Define the list of input parameters
    opt <- list(
        fasta = choose.files(
            caption = "Choose Fasta file of the known protein (UniProt)!",
            multi = FALSE),
        taxon = readline(
            prompt = "Provide the taxon identifier!"),
        columns = readline(
            prompt = "Give comma-separated annotation names to retrieve!"),
        key = readline(
            prompt = "Provide the key name to use for cross-annotation!"),
        separator = readline(
            prompt = "Provide the separator character(s) to parse annotation!"),
        output = readline(
            prompt = "Define the output file name!"))
    
} else {
    
    # Define the list of command line parameters
    option_list <- list(
        make_option(
            opt_str = c("-f", "--fasta"),
            type = "character", default = NULL,
            help = "Known protein fasta file",
            metavar = "character"),
        make_option(
            opt_str = c("-t", "--taxon"),
            type = "integer", default = NULL,
            help = "The taxon identifier for the species of interest",
            metavar = "character"),
        make_option(
            opt_str = c("-c", "--columns"),
            type = "character", default = NULL,
            help = "The annotation columns (comma-separated) to retrieve!",
            metavar = "character"),
        make_option(
            opt_str = c("-k", "--key"),
            type = "character", default = NULL,
            help = "The key name to use for cross-annotation!",
            metavar = "character"),
        make_option(
            opt_str = c("-s", "--separator"),
            type = "character", default = NULL,
            help = "The separator characters to parse annotation!",
            metavar = "character"),
        make_option(
            opt_str = c("-o", "--output"),
            type = "character", default = "", 
            help = "Output file name [default= %default]",
            metavar = "character"))
    
    # Parse the parameters provided on command line by user
    opt_parser <- OptionParser(option_list = option_list)
    opt <- parse_args(opt_parser)
    
}

# Check whether inputs parameter were provided
if (
    is.null(opt$fasta) | is.null(opt$taxon) |
    is.null(opt$columns) | is.null(opt$key)) {
    
    print_help(opt_parser)
    
}

# Store the annotation column name into vector
opt$columns %<>%
    strsplit(x = ., split = ",", fixed = TRUE) %>%
    unlist(.)

# Check whether output parameter was provided
if (opt$output == "") {
    
    opt$output <- sub(
        pattern = "^(.*)(\\..*)?",
        replacement = "\\1_annot.txt",
        x = opt$fasta)
    warning(paste0(
        "Output results to '",
        opt$output,
        "'!"))
    
}

# Create output directory if not already existing
dir.create(dirname(opt$output))



### Data import ----------------------------------------------------------

# Import the fasta files
fasta <- seqinr::read.fasta(
    file = opt$fasta, seqtype = "AA", as.string = TRUE)



### Obtain protein annotation from UniProt -------------------------------

# Clean up the UniProtKB
names(fasta) <- uni_id_clean(names(fasta))

#
if (grepl("GN=", attr(fasta[[1]], "Annot"))) {
    
    data_final <- lapply(fasta, function(x) {
        attr(x, "Annot")
    }) %>%
        set_names(names(fasta)) %>%
        plyr::ldply(., "data.frame", .id = opt[["key"]]) %>%
        dplyr::mutate(
            ., 
            `ENTRY-NAME` = sub(".+\\|([^ ]+).+", "\\1", `data.frame`),
            `GENE-NAME` = ifelse(
                grepl("GN=", `data.frame`, fixed = TRUE),
                sub(".+ GN=([^ ]+).+", "\\1", `data.frame`),
                NA),
            `PROTEIN-NAME` = sub("[^ ]+ (.+) OS=.+", "\\1", `data.frame`),
            TaxonID = sub(".+ OX=([^ ]+).+", "\\1", `data.frame`),
            Taxon = sub(".+ OS=(.+) OX=.+", "\\1", `data.frame`)) %>%
        dplyr::select(., -`data.frame`)
    
} else {
    
    warning("Missing gene name. Not a uniprot fasta?")
    data_final <- lapply(fasta, function(x) {
        attr(x, "Annot")
    }) %>%
        set_names(names(fasta)) %>%
        plyr::ldply(., "data.frame", .id = opt[["key"]]) %>%
        dplyr::mutate(., `data.frame` = sub("^>", "", data.frame)) %>%
        tidyr::separate(
            data = ., col = "data.frame",
            into = opt$columns, sep = opt$separator)
        
}



### Results export -------------------------------------------------------

# Export the bsu ideogram for reuse at later stage
write.table(
    x = data_final,
    file = opt$output,
    quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)

# Define end time
print(paste("Complete", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))


