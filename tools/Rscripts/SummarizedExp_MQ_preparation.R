#!/usr/bin/env Rscript

# This script generates the required expression and phenodata files
# from MaxQuant txt folder (using protein groups and sites)



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
library(magrittr)
library(optparse)



### Parameters setting up ------------------------------------------------

# Define input parameters (interactively or from command line)
if (interactive()) {
    
    # Define the list of input parameters
    opt <- list(
        input = choose.files(
            caption = "Choose input MaxQuant txt file!",
            multi = FALSE),
        prefix = readline(
            prompt = "Define a prefix for the output filenames!"),
        output = readline(
            prompt = "Define the output folder!"))
    
} else {
    
    # Define the list of command line parameters
    option_list <- list(
        make_option(
            opt_str = c("-i", "--input"),
            type = "character", default = "",
            help = "Input MaxQuant txt file (containing abundance columns)",
            metavar = "character"),
        make_option(
            opt_str = c("-p", "--prefix"),
            type = "character", default = "", 
            help = "Define a prefix for the output filenames",
            metavar = "character"),
        make_option(
            opt_str = c("-o", "--output"),
            type = "character", default = "SummarizedExp", 
            help = "Output folder name [default= %default]",
            metavar = "character"))
    
    # Parse the parameters provided on command line by user
    opt_parser <- OptionParser(option_list = option_list)
    opt <- parse_args(opt_parser)
    
}

# Check whether input parameter was provided
if (
    identical(opt$input, "") |
    identical(opt$input, character(0))) {
    
    print_help(opt_parser)
    stop(paste(
        "MaxQuant txt file required!"))
    
}

# Check whether output parameter was provided
if (
    identical(opt$output, NULL) |
    identical(opt$output, "") |
    identical(opt$output, character(0))) {
    
    opt$output <- dirname(opt$input)
    warning(paste0("Output results to ", opt$output, "!"))
    
}

# Create output directory if not already existing
dir.create(opt$output)



### Import data ----------------------------------------------------------

# Vector holding the file(s) to import and the quantification
# columns of interest (can be more than one type of columns)
what_to_look_for <- list(
    proteinGroups = c("iBAQ", "Intensity"),
    Sites = c("Intensity")
)

# Find out the file type for the expression table and therefore which
# columns to retrieve
my_col <- lapply(names(what_to_look_for), function(i) {
    if (grepl(i, opt[["input"]])) {
        what_to_look_for[i]
    }
}) %>%
    unlist(., recursive = FALSE)

# Check that the file type can be accurately defined
if (length(names(my_col)) > 1) {
    stop(paste(
        "The file type cannot be determined, could be either:",
        paste0(names(my_col), collapse = " or ")))
} else {
    my_col %<>%
        unlist(., use.names = FALSE)
}

# Import the expression data
my_data <- data.table::fread(
    input = opt[["input"]], sep = "\t", quote = "", header = TRUE,
    stringsAsFactors = FALSE, colClasses = "character")



### Data formatting ------------------------------------------------------

# Loop through each column type to retrieve
my_pheno <- lapply(my_col, function(j) {
    required_file_summExp(
        my_data = my_data, my_col = j, input_path = x)
}) %>%
    plyr::ldply(., "data.table", .id = NULL)



### Results export -------------------------------------------------------

# Define output filename
out_name <- opt[["input"]] %>%
    basename(.) %>%
    c(opt[["prefix"]], .) %>%
    paste0(., collapse = "_") %>%
    paste0(opt[["output"]], "/", .) %>%
    sub("\\.txt", "", .)

# Write the phenodata to output folder
data.table::fwrite(
    x = my_pheno, file = paste0(out_name, "_pheno.txt"),
    append = FALSE, quote = FALSE, sep = "\t",
    row.names = FALSE, col.names = TRUE)

# Copy the expression data to output folder with matching name
file.copy(
    from = opt[["input"]],
    to = paste0(out_name, "_expr.txt"),
    overwrite = FALSE)

# Define end time
print(paste("Complete", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))


