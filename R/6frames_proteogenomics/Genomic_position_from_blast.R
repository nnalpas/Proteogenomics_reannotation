#!/usr/bin/env Rscript

# This script determines the reference protein genomic coordinates based on
# blast results against ORF of known coordinates



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
        "/home-link",
        user,
        "bin/General_function.R",
        sep = "/"))

# Load the required packages (or install if not already in library)
load_package(plyr)
load_package(dplyr)
load_package(magrittr)
load_package(data.table)
load_package(splitstackshape)
load_package(stringr)
load_package(optparse)
load_package(seqinr)



### Parameters setting up ------------------------------------------------

# Define the list of command line parameters
option_list <- list(
    make_option(
        opt_str = c("-f", "--fasta"),
        type = "character", default = NULL,
        help = "Fasta file", metavar = "character"),
    make_option(
        opt_str = c("-b", "--blast"),
        type = "character", default = NULL,
        help = "Blast file", metavar = "character"),
    make_option(
        opt_str = c("-g", "--genomic_position"),
        type = "character", default = NULL,
        help = "Genomic position file", metavar = "character"),
    make_option(
        opt_str = c("-o", "--output"),
        type = "character", default = "protein_location.txt", 
        help = "Output file name [default= %default]", metavar = "character"))

# Parse the parameters provided on command line by user
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# Check whether inputs parameter was provided
if (is.null(opt$fasta) | is.null(opt$blast) | is.null(opt$genomic_position)) {
    
    print_help(opt_parser)
    stop(paste(
        "The three input arguments must be supplied",
        "(fasta, blast and genomic position files)!"))
    
}

# Check whether output parameter was provided
if (opt$output == "protein_location.txt") {
    
    warning("Output results to protein_location.txt!")
    
}



### Data import ----------------------------------------------------------

# Import the fasta file
fasta <- opt$fasta %>%
    as.character(.) %>%
    seqinr::read.fasta(
        file = ., seqtype = "AA", as.string = TRUE) %>%
    set_names(uni_id_clean(names(.)))

# Get all protein IDs
uni_ids <- names(fasta) %>%
    as.character(.) %>%
    unique(.)

# Import the Blast results
blast_data <- opt$blast%>%
    as.character(.) %>%
    read.table(
        file = .,
        header = FALSE,
        sep = "\t",
        quote = "",
        col.names = c(
            "qseqid", "sseqid", "pident", "nident", "mismatch", "length",
            "gapopen", "qstart", "qend", "sstart", "send", "evalue",
            "bitscore", "score"),
        as.is = TRUE) %>%
    dplyr::mutate(
        .,
        qseqid = uni_id_clean(qseqid),
        sseqid = uni_id_clean(sseqid))

# Import the genomic position data
orf_pos <- opt$genomic_position %>%
    as.character(.) %>%
    read.table(
        file = ., header = TRUE, sep = "\t", quote = "")



### Get target ID position based on blast --------------------------------

# Determine the best blast for each queried protein
best_blast_data <- best_blast(data = blast_data, key = "qseqid") %>%
    dplyr::filter(., evalue < 0.0001)

# Give warning when a single best blast could not be found per query
if (any(best_blast_data$best_count != 1)) {
    
    warning(paste(
        "There are",
        length(unique(
            best_blast_data[best_blast_data$best_count != 1, "qseqid"])),
        "queries without a single best blast!", sep = " "))
    
}

# Add the orf position to the blast data
best_blast_pos <- best_blast_data %>%
    dplyr::left_join(x = ., y = orf_pos, by = c("sseqid" = "id"))

# Compute aa and nucleotide length of reference entries
ref_pos <- lapply(X = fasta, FUN = nchar) %>%
    plyr::ldply(., "data.frame") %>%
    set_colnames(c("id", "aa_length")) %>%
    dplyr::mutate(., nucl_length = aa_length * 3)

# Get the genomic start position for reference entries
best_blast_pos %<>%
    dplyr::mutate(., ref_start = case_when(
        .$strand == 1 ~ (.$start + (.$sstart * 3) - (.$qstart * 3)),
        .$strand == -1 ~ (.$start - (.$sstart * 3) + (.$qstart * 3))))

# Add the reference entries start position
ref_pos %<>%
    dplyr::left_join(
        x = .,
        y = best_blast_pos %>%
            dplyr::select(., qseqid, ref_start, strand, frame, best_count),
        by = c("id" = "qseqid"))

# Compute the reference entries end position
ref_pos %<>%
    dplyr::mutate(., ref_end = case_when(
        is.na(.$strand) ~ NA_integer_,
        .$strand == 1 ~ as.integer(.$ref_start + (.$aa_length * 3) - 1),
        .$strand == -1 ~ as.integer(.$ref_start - (.$aa_length * 3) + 1))) %>%
    dplyr::select(
        ., id, ref_start, ref_end, nucl_length, aa_length, strand, frame)
colnames(ref_pos) %<>%
    sub("^ref_", "", .)

# Export the reference entries coordinates data
write.table(
    x = ref_pos, file = opt$output,
    quote = FALSE, sep = "\t",
    row.names = FALSE, col.names = TRUE)

# Define end time
print(paste("Complete", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))


