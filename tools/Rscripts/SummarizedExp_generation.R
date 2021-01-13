

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

my_plots <- list()



### Parameters setting up ------------------------------------------------

#
opt <- list(
    my_folder = "H:/data/Synechocystis_6frame/Kopf_2014_TU",
    my_expr = NULL,
    my_pheno = NULL,
    output = "C:/Users/kxmna01/Desktop/2021-01-13_SummarizedExp"
)

# Check whether output parameter was provided
if (
    identical(opt$output, NULL) |
    identical(opt$output, "") |
    identical(opt$output, character(0))) {
    
    opt$output <- dirname(opt$novel_reason)
    warning(paste0("Output results to ", opt$output, "!"))
    
}

# Create output directory if not already existing
dir.create(opt$output)



### Import data ----------------------------------------------------------

# Identify all files to import
if (!is.null(opt[["my_folder"]])) {
    my_exprs <- list.files(
        path = opt[["my_folder"]], pattern = "*_expr.txt", all.files = TRUE,
        full.names = TRUE, recursive = TRUE) %>%
        data.frame(
            Name = basename(.) %>% sub("_expr.txt", "", .),
            Expr = .,
            stringsAsFactors = FALSE)
    my_phenos <- list.files(
        path = opt[["my_folder"]], pattern = "*pheno.txt", all.files = TRUE,
        full.names = TRUE, recursive = TRUE) %>%
        data.frame(
            Name = basename(.) %>% sub("_pheno.txt", "", .),
            Pheno = .,
            stringsAsFactors = FALSE)
    my_files <- dplyr::full_join(x = my_exprs, y = my_phenos)
} else {
    my_files <- data.frame(
        Name = opt[["my_expr"]] %>% basename(.) %>% sub("_expr.txt", "", .),
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
# tmp <- MultiAssayExperiment(all_summexp)

