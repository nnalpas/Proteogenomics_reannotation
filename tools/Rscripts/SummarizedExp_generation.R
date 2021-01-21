

#!/usr/bin/env Rscript

# This script generates a SummarizedExperiment object from inputs



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

# Load the required packages
library(magrittr)
library(MultiAssayExperiment)
library(optparse)



### Parameters setting up ------------------------------------------------

# Define input parameters (interactively or from command line)
if (interactive()) {
    
    # Define the list of input parameters
    opt <- list(
        my_folder = readline(
            prompt = "Define the input folder (containing multiple expr and pheno files!"),
        my_expr = choose.files(
            caption = "Choose an expression value file!",
            multi = FALSE),
        my_pheno = choose.files(
            caption = "Choose a phenodata file!",
            multi = FALSE),
        prefix = readline(
            prompt = "Define a prefix name for the output multi assay file!"),
        output = readline(
            prompt = "Define the output folder!"))
    
} else {
    
    # Define the list of command line parameters
    option_list <- list(
        make_option(
            opt_str = c("-i", "--my_folder"),
            type = "character", default = "",
            help = "Input folder (containing multiple expr and pheno files)",
            metavar = "character"),
        make_option(
            opt_str = c("-e", "--my_expr"),
            type = "character", default = "",
            help = "Input expr file",
            metavar = "character"),
        make_option(
            opt_str = c("-x", "--my_pheno"),
            type = "character", default = "",
            help = "Input pheno file",
            metavar = "character"),
        make_option(
            opt_str = c("-p", "--prefix"),
            type = "character", default = "", 
            help = "Define a prefix for the multi assay file",
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
if ((identical(opt$my_folder, NULL) |
     identical(opt$my_folder, "") |
     identical(opt$my_folder, character(0))) &
    ((identical(opt$my_expr, NULL) |
      identical(opt$my_expr, "") |
      identical(opt$my_expr, character(0))) |
     (identical(opt$my_pheno, NULL) |
      identical(opt$my_pheno, "") |
      identical(opt$my_pheno, character(0))))) {
    
    print_help(opt_parser)
    stop(paste(
        "Either a input folder or both expr and pheno inputs are required!"))
    
}

# Check whether output parameter was provided
if (
    identical(opt$output, NULL) |
    identical(opt$output, "") |
    identical(opt$output, character(0))) {
    
    opt$output <- ifelse(
        is.null(opt$my_folder), dirname(opt$my_expr), opt$my_folder)
    warning(paste0("Output results to ", opt$output, "!"))
    
}

# Create output directory if not already existing
dir.create(opt$output)



### Import data ----------------------------------------------------------

# Identify all files to import
if (!is.null(opt[["my_folder"]])) {
    my_exprs <- list.files(
        path = opt[["my_folder"]], pattern = "*_expr\\.(txt|RDS)", all.files = TRUE,
        full.names = TRUE, recursive = TRUE) %>%
        data.frame(
            Name = basename(.) %>% sub("_expr\\.(txt|RDS)", "", .),
            Expr = .,
            stringsAsFactors = FALSE)
    my_phenos <- list.files(
        path = opt[["my_folder"]], pattern = "*pheno\\.txt", all.files = TRUE,
        full.names = TRUE, recursive = TRUE) %>%
        data.frame(
            Name = basename(.) %>% sub("_pheno.txt", "", .),
            Pheno = .,
            stringsAsFactors = FALSE)
    my_files <- dplyr::full_join(x = my_exprs, y = my_phenos)
} else {
    my_files <- data.frame(
        Name = opt[["my_expr"]] %>%
            basename(.) %>%
            sub("_expr\\.(txt|RDS)", "", .),
        Expr = opt[["my_expr"]],
        Pheno = opt[["my_pheno"]],
        stringsAsFactors = FALSE)
}

# Check that expr and pheno data are matching
if (any(is.na(my_files))) {
    stop("Cannot identify all corresponding 'expr' and 'pheno' data!")
}



### Data formatting ------------------------------------------------------

# Loop through all datasets and generate for each a summarized experiment
all_summexp <- apply(X = my_files, MARGIN = 1, FUN = function(x) {
    make_SummarizedExperiment(name = x[[1]], assay = x[[2]], coldata = x[[3]])
}) %>%
    set_names(my_files[["Name"]])

# Provisional code to generate a multi assay object
multi_assay <- MultiAssayExperiment(all_summexp)



### Exporting results ----------------------------------------------------

# Export all summarized experiment object for later re-use
for (x in names(all_summexp)) {
    saveRDS(
        object = all_summexp[[x]],
        file = paste0(opt$output, "/", x, "_SummExp.RDS"))
}

# Export all summarized experiment object for later re-use
saveRDS(
    object = multi_assay,
    file = paste0(opt$output, "/", opt$prefix, "MultiAssay.RDS"))


