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
if (interactive()) {
    source(
        file = paste(
            "C:/Users",
            user,
            "Documents/GitHub/Miscellaneous/R/6frames_proteogenomics/helper.R",
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



### Parameters setting up ------------------------------------------------

# Define input parameters (interactively or from command line)
if (interactive()) {
    
    # Define the list of input parameters
    opt <- list(
        input = choose.files(
            caption = "Choose input Blast results",
            multi = FALSE),
        filter = readline(
            prompt = paste(
                "What filter to use to determine best blast",
                "(do not provide value for default)?")),
        multi_match = readline(
            prompt = paste(
                "What to do for multihit entries",
                "(either: remove, keep, uniquify)?")),
        output = ".")
    
} else {
    
    # Define the list of command line parameters
    option_list <- list(
        make_option(
            opt_str = c("-i", "--input"), type = "character",
            default = NULL, help = "Blast data file name",
            metavar = "character"),
        make_option(
            opt_str = c("-f", "--filter"), type = "character",
            default = NULL, help = "Specific filtering to apply",
            metavar = "character"),
        make_option(
            opt_str = c("-m", "--multi_match"), type = "character",
            default = NULL, help = "Filter for multi hits entry",
            metavar = "character"),
        make_option(
            opt_str = c("-o", "--output"), type = "character",
            default = ".", help = "Output directory",
            metavar = "character"))
    
    # Parse the parameters provided on command line by user
    opt_parser <- OptionParser(option_list = option_list)
    opt <- parse_args(opt_parser)
    
}

# Check whether input parameter was provided
if (is.null(opt$input)){
    
    print_help(opt_parser)
    stop("At least one argument must be supplied (input file)!")
    
}

# Check whether output parameter was provided
if (is.null(opt$output)){
    
    opt$output <- dirname(opt$input)
    warning(paste("Output results to path: ", opt$output, "!", sep = ""))
    
}

# If filter and multi_match parameters are undefined, define as null
if (is.null(opt$filter)) {
    opt["filter"] <- list(NULL)
} else if (opt$filter == "") {
    opt["filter"] <- list(NULL)
} else {
    opt$filter <- enquote(opt$filter)
}
if (is.null(opt$multi_match)) {
    opt["multi_match"] <- list(NULL)
} else if (opt$multi_match == "") {
    opt["multi_match"] <- list(NULL)
}

# Create output directory if not already existing
dir.create(opt$output)



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
    file = paste0(opt$output, "/Reciprocal_id_", basename(opt$input)),
    quote = FALSE,
    sep = "\t",
    row.names = FALSE,
    col.names = FALSE)

# Export the best hits results
write.table(
    x = best_blast_data,
    file = paste0(opt$output, "/Best_blast_", basename(opt$input)),
    quote = FALSE,
    sep = "\t",
    row.names = FALSE,
    col.names = TRUE)

# Export the blast ID map results
write.table(
    x = best_blast_data %>%
        dplyr::select(., qseqid, sseqid),
    file = paste0(opt$output, "/Blast_cross-map_", basename(opt$input)),
    quote = FALSE,
    sep = "\t",
    row.names = FALSE,
    col.names = TRUE)

# Define end time
print(paste("Complete", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))


