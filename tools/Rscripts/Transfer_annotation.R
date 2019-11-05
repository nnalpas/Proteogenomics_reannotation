#!/usr/bin/env Rscript

# This script transfers annotation from known proteins to putative ORFs



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
            Sys.getenv("HOME"),
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



### Parameters setting up ------------------------------------------------

# Define input parameters (interactively or from command line)
if (interactive()) {
    
    # Define the list of input parameters
    opt <- list(
        crossmap = choose.files(
            caption = "Choose the cross-mapping entries file!",
            multi = FALSE),
        annotation = choose.files(
            caption = "Choose the annotation file!",
            multi = FALSE),
        key = readline(
            prompt = "The key column name to join on!"),
        rename = readline(
            prompt = "The other column name for joining to be renamed!"),
        output = readline(
            prompt = "Define the output file!"))
    
} else {
    
    # Define the list of command line parameters
    option_list <- list(
        make_option(
            opt_str = c("-c", "--crossmap"),
            type = "character", default = NULL,
            help = "Cross-mapping entries file",
            metavar = "character"),
        make_option(
            opt_str = c("-a", "--annotation"),
            type = "character", default = NULL,
            help = "Annotation file",
            metavar = "character"),
        make_option(
            opt_str = c("-k", "--key"),
            type = "character", default = "",
            help = "Key column name for joining",
            metavar = "character"),
        make_option(
            opt_str = c("-r", "--rename"),
            type = "character", default = "",
            help = "The other key column name for joining to be renamed",
            metavar = "character"),
        make_option(
            opt_str = c("-o", "--output"),
            type = "character", default = "",
            help = "Output filename [default= %default]",
            metavar = "character"))
    
    # Parse the parameters provided on command line by user
    opt_parser <- OptionParser(option_list = option_list)
    opt <- parse_args(opt_parser)
    
}

# Check whether inputs parameter were provided
if (is.null(opt$crossmap)) {
    
    print_help(opt_parser)
    stop(paste(
        "Cross-mapping entries required!"))
    
}
if (is.null(opt$annotation)) {
    
    print_help(opt_parser)
    stop(paste(
        "Annotation file is required!"))
    
}
if (opt$key == "") {
    
    print_help(opt_parser)
    stop(paste(
        "Key column name is required!"))
    
}
if (opt$rename == "") {
    opt <- list(rename = NULL)
}

# Check whether output parameter was provided
if (opt$output == "") {
    
    opt$output <- sub(
        pattern = "^(.*)(\\..*)?",
        replacement = "\\1_annot.txt",
        x = opt$crossmap)
    warning(paste0(
        "Output results to '",
        opt$output,
        "'!"))
    
}

# Create output directory if not already existing
dir.create(dirname(opt$output))



### Data import ----------------------------------------------------------

# Import the cross-mapping data
crossmap <- data.table::fread(
    input = opt$crossmap, header = TRUE, sep = "\t",
    quote = "", stringsAsFactors = TRUE)

# Import the annotation data
annotations <- data.table::fread(
    input = opt$annotation, header = TRUE, sep = "\t",
    quote = "", stringsAsFactors = TRUE)



### Transfer of annotation -----------------------------------------------

# Rename specific column to key for the two datasets
if (any(colnames(crossmap) %in% opt$rename)) {
    colnames(crossmap)[colnames(crossmap) == opt$rename] <- opt$key
}
if (any(colnames(annotations) %in% opt$rename)) {
    colnames(annotations)[colnames(annotations) == opt$rename] <- opt$key
}

# Merge the cross-map data with the annotation based on main key
data <- dplyr::left_join(x = crossmap, y = annotations, by = opt$key) %>%
    dplyr::group_by(., qseqid) %>%
    dplyr::summarise_all(funs(paste(unique(.), collapse = ";")))



### Results export -------------------------------------------------------

# Export the annotated data
write.table(
    x = data, file = opt$output, quote = FALSE, sep = "\t",
    row.names = FALSE, col.names = TRUE)

# Define end time
print(paste("Complete", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))


