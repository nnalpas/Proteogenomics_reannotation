#!/usr/bin/env Rscript

# This script compiles blast and reciprocal blast hits from blasting results
# to identify reciprocal best hits and matching position



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



### Parameters setting up ------------------------------------------------

# Define input parameters (interactively or from command line)
if (interactive()) {
    
    # Define the list of input parameters
    opt <- list(
        blast = choose.files(
            caption = "Choose input Blast results",
            multi = FALSE),
        reciprocal_blast = choose.files(
            caption = "Choose input reciprocal Blast results",
            multi = FALSE),
        output = NULL)
    
} else {
    
    # Define the list of command line parameters
    option_list <- list(
        make_option(
            opt_str = c("-b", "--blast"),
            type = "character", default = NULL, 
            help = "Blast data file name", metavar = "character"),
        make_option(
            opt_str = c("-r", "--reciprocal_blast"),
            type = "character", default = NULL, 
            help = "Reciprocal blast data file name", metavar = "character"),
        make_option(
            opt_str = c("-o", "--output"),
            type = "character", default = NULL, 
            help = "Output directory", metavar = "character"))
    
    # Parse the parameters provided on command line by user
    opt_parser <- OptionParser(option_list = option_list)
    opt <- parse_args(opt_parser)
    
}

# Check whether inputs parameter was provided
if (is.null(opt$blast) | is.null(opt$reciprocal_blast)){
    
    print_help(opt_parser)
    stop(paste(
        "The two input arguments must be supplied",
        "(blast and reciprocal blast files)!"))
    
}

# Check whether output parameter was provided
if (is.null(opt$output)){
    
    opt$output <- dirname(opt$reciprocal_blast)
    warning(paste("Output results to path: ", opt$output, "!", sep = ""))
    
}

# Create output directory if not already existing
dir.create(opt$output)



### Data import ----------------------------------------------------------

# Import the Blast results
blast_data <- blast_read(file = opt$blast, blast_format = "6") %>%
    dplyr::mutate(
        .,
        qseqid = uni_id_clean(qseqid),
        sseqid = uni_id_clean(sseqid))

# Import the reciprocal Blast results
reciproc_data <- blast_read(
    file = opt$reciprocal_blast, blast_format = "6") %>%
    dplyr::mutate(
        .,
        qseqid = uni_id_clean(qseqid),
        sseqid = uni_id_clean(sseqid))



### Reciprocal best blast hits identification ----------------------------

# Get the best blast match for each query
best_blast_data <- best_blast(data = blast_data, key = "qseqid")

# Get the best reciprocal blast match for each query
best_reciproc_data <- best_blast(data = reciproc_data, key = "qseqid")

# Merge the best blast and best reciprocal data
blast_merge <- best_blast_data %>%
    dplyr::left_join(
        x = .,
        y = best_reciproc_data,
        by = c("qseqid" = "sseqid", "sseqid" = "qseqid"),
        suffix = c("_blast", "_reciproc"))

# Determine which entry have a reciprocal best hits
# this requires filtering out all query ID without a reciprocal best hit
# and checking how many distinct subject ID are matched per query id
# as well as checking how many distinct subject sequence exist
# (the same sequence can have different IDs across different taxon)
blast_merge_confirmed <- blast_merge %>%
    dplyr::filter(., !is.na(qseq_reciproc)) %>%
    dplyr::group_by(., qseqid) %>%
    dplyr::mutate(
        .,
        final_count = n_distinct(sseqid),
        seq_count = n_distinct(qseq_reciproc)) %>%
    dplyr::ungroup(.) %>%
    dplyr::mutate(
        ., 
        best_reciprocal_type = dplyr::case_when(
            is.na(sseqid) ~ "None",
            final_count == 1 ~ "Single",
            TRUE ~ "Multiple"),
        comment = dplyr::case_when(
            seq_count == 1 ~ "Single",
            seq_count > 1 ~ "Multiple"))



### Results export -------------------------------------------------------

# Export the best hits results
write.table(
    x = blast_merge_confirmed,
    file = paste0(
        opt$output, "/Best_Reciproc_Blast_",
        basename(opt$reciprocal_blast)),
    quote = FALSE,
    sep = "\t",
    row.names = FALSE,
    col.names = TRUE)

# Export the blast ID map results
write.table(
    x = blast_merge_confirmed %>%
        dplyr::select(., qseqid, sseqid, best_reciprocal_type, comment),
    file = paste0(
        opt$output, "/Best_Reciproc_Blast_cross-map_",
        basename(opt$reciprocal_blast)),
    quote = FALSE,
    sep = "\t",
    row.names = FALSE,
    col.names = TRUE)

# Define end time
print(paste("Complete", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))


