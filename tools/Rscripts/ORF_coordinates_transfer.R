#!/usr/bin/env Rscript

# This script determines the ORF genomic coordinates of a set of entries
# based on blast results against known coordinate entries



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
#source(
#    file = paste(
#        "C:/Users",
#        user,
#        "Documents/GitHub/Miscellaneous/R/General/General_function.R",
#        sep = "/"))
source(
    file = paste(
        Sys.getenv("HOME"),
        "bin/General_function.R",
        sep = "/"))

# Load the required packages (or install if not already in library)
load_package("plyr")
load_package("dplyr")
load_package("tidyr")
load_package("magrittr")
load_package("data.table")
load_package("splitstackshape")
load_package("stringr")
load_package("optparse")
load_package("seqinr")



### Parameters setting up ------------------------------------------------

# Define the list of command line parameters
option_list <- list(
    make_option(
        opt_str = c("-i", "--idmap"),
        type = "character", default = NULL,
        help = "ID map file", metavar = "character"),
    make_option(
        opt_str = c("-c", "--coordinates"),
        type = "character", default = NULL,
        help = "Coordinates file", metavar = "character"),
    make_option(
        opt_str = c("-o", "--output"),
        type = "character", default = "orf_coordinates.txt",
        help = "Output file name [default= %default]", metavar = "character"))

# Parse the parameters provided on command line by user
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# Check whether inputs parameter was provided
if (is.null(opt$idmap) | is.null(opt$coordinates)) {
    
    print_help(opt_parser)
    stop("The idmap and coordinates arguments must be supplied!")
    
}

# Check whether output parameter was provided
if (opt$output == "orf_coordinates.txt") {
    
    warning("Output results to orf_coordinates.txt!")
    
}

# For manual parameters set-up
#opt <- list(
#    idmap = "C:/Users/kxmna01/Dropbox/Home_work_sync/Work/Colleagues shared work/Vaishnavi Ravikumar/Bacillus_subtilis_6frame/Databases/Blast_cross-map_blastp_ORF_Nicolas_vs_Karsten_11012017",
#    coordinates = "C:/Users/kxmna01/Dropbox/Home_work_sync/Work/Colleagues shared work/Vaishnavi Ravikumar/Bacillus_subtilis_6frame/Databases/Find0_GCA_000009045.1_FIXED_orf_coordinates.txt",
#    output = "C:/Users/kxmna01/Dropbox/Home_work_sync/Work/Colleagues shared work/Vaishnavi Ravikumar/Bacillus_subtilis_6frame/Databases/6frames_orf_coordinates.txt")



### Data import ----------------------------------------------------------

# Import the ID mapping file
idmap <- opt$idmap %>%
    as.character(.) %>%
    read.table(
        file = ., header = FALSE,
        sep = "\t", as.is = TRUE)

# Import the coordinates file
coordinates <- opt$coordinates %>%
    as.character(.) %>%
    read.table(
        file = ., header = TRUE,
        sep = "\t", as.is = TRUE)



### Get genomic coordinates ----------------------------------------------

# Detect which columns should be used for merging coordinates and IDmap
if (!any(idmap$V1 %in% coordinates$id) & !any(idmap$V2 %in% coordinates$id)) {
    stop("Impossible to differentiate which column to use for merging!")
} else if (
    !all(idmap$V1 %in% coordinates$id) & !all(idmap$V2 %in% coordinates$id)) {
    stop("Not all ID for merging can be found in coordinates variable!")
}
if (all(idmap$V1 %in% coordinates$id)) {
    print(paste(
        "Merging on column 1, starting with:",
        paste(head(idmap$V1), collapse = ";"),
        "!"))
    colnames(idmap) <- c("association_id", "id")
    orf_pos <- dplyr::left_join(
        x = idmap, y = coordinates, by = c("association_id" = "id"))
} else if (all(idmap$V2 %in% coordinates$id)) {
    print(paste(
        "Merging on column 2, starting with:",
        paste(head(idmap$V2), collapse = ";"),
        "!"))
    colnames(idmap) <- c("id", "association_id")
    orf_pos <- dplyr::left_join(
        x = idmap, y = coordinates, by = c("association_id" = "id"))
}

# Export the ORF coordinates data
write.table(
    x = orf_pos, file = opt$output,
    quote = FALSE, sep = "\t",
    row.names = FALSE, col.names = TRUE)

# Define end time
print(paste("Complete", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))


