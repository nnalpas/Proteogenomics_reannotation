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
            caption = "Choose input MaxQuant txt file or GRange RDS file!",
            multi = FALSE),
        prefix = readline(
            prompt = "Define a prefix for the output filenames!"),
        grange = readline(
            prompt = "The input is a GRange object (logical)!") %>%
            as.logical(.),
        output = readline(
            prompt = "Define the output folder!"))
    
} else {
    
    # Define the list of command line parameters
    option_list <- list(
        make_option(
            opt_str = c("-i", "--input"),
            type = "character", default = "",
            help = "Input MaxQuant txt file or GRange RDS file (containing abundance columns)",
            metavar = "character"),
        make_option(
            opt_str = c("-p", "--prefix"),
            type = "character", default = "", 
            help = "Define a prefix for the output filenames",
            metavar = "character"),
        make_option(
            opt_str = c("-g", "--grange"),
            type = "logical", default = FALSE, 
            help = "The input is a GRange object [default= %default]",
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
        "MaxQuant txt or GRange RDS file is required!"))
    
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

# Which column to use as ID
what_col_for_id <- list(
    proteinGroups = c("Protein IDs")
)

# Find out the file type for the expression table and therefore which
# columns to retrieve
my_file_type <- lapply(names(what_to_look_for), function(i) {
    if (grepl(i, opt[["input"]])) {
        i
    }
}) %>%
    unlist(.)

# Check that the file type can be accurately defined
if (length(my_file_type) != 1) {
    stop(paste(
        "The file type cannot be determined, could be either:",
        paste0(my_file_type, collapse = " or ")))
}

# Import the expression data
if (opt[["grange"]]) {
    my_data <- readRDS(opt[["input"]]) %>%
        values(.)
} else {
    my_data <- data.table::fread(
        input = opt[["input"]], sep = "\t", quote = "", header = TRUE,
        stringsAsFactors = FALSE, colClasses = "character")
}



### Data formatting ------------------------------------------------------

# What column should be retrieved based on file type
my_cols <- what_to_look_for[[my_file_type]] %>%
    unlist(., use.names = FALSE)

# Loop through each column type to retrieve
my_pheno <- lapply(my_cols, function(j) {
    required_file_summExp(
        my_data = my_data, my_col = j, input_path = x)
}) %>%
    plyr::ldply(., "data.table", .id = NULL)

# In case imported file is not a grange object, add a 'id' column
if (!opt[["grange"]]) {
    
    my_data_format <- my_data
    
    # Before renaming column name to id, make sure an 'id' column
    # is not already existing
    if ("id" %in% colnames(my_data_format)) {
        
        my_data_format %<>%
            dplyr::rename(., id_OLD = id)
        warning("A column called 'id' was renamed to 'id_OLD'!")
        
    }
    
    # Define the column used as id
    my_id <- what_col_for_id[[my_file_type]]
    
    # If such column exist then 
    if (my_id %in% colnames(my_data_format)) {
        my_data_format %<>%
            dplyr::mutate(., id = !!as.name(my_id))
    } else {
        stop(paste0(
            "The column: ", my_id,
            " cannot be used as 'id', as it is missing!"))
    }
    
}



### Results export -------------------------------------------------------

# Define output filename
out_name <- opt[["input"]] %>%
    basename(.) %>%
    c(opt[["prefix"]], .) %>%
    paste0(., collapse = "_") %>%
    paste0(opt[["output"]], "/", .) %>%
    sub("\\.(txt|RDS)", "", .)

# Write the phenodata to output folder
data.table::fwrite(
    x = my_pheno, file = paste0(out_name, "_pheno.txt"),
    append = FALSE, quote = FALSE, sep = "\t",
    row.names = FALSE, col.names = TRUE)

# Write the expression data to output folder with matching name
if (opt[["grange"]]) {
    file.copy(
        from = opt[["input"]],
        to = paste0(out_name, "_expr.RDS"),
        overwrite = FALSE)
} else {
    data.table::fwrite(
        x = my_data_format, file = paste0(out_name, "_expr.txt"),
        append = FALSE, quote = FALSE, sep = "\t",
        row.names = FALSE, col.names = TRUE)
}

# Define end time
print(paste("Complete", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))


