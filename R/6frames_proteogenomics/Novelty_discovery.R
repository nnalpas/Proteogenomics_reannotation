#!/usr/bin/env Rscript

# This script determines the novelty reasons for each peptide and ORF based on
# genomic coordinates, blast results...



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
        opt_str = c("-e", "--evidence"),
        type = "character", default = NULL,
        help = "Evidence with peptide group info", metavar = "character"),
    make_option(
        opt_str = c("-p", "--peptide_location"),
        type = "character", default = NULL,
        help = "Peptide location file", metavar = "character"),
    make_option(
        opt_str = c("-r", "--reference_fasta"),
        type = "character", default = NULL,
        help = "Reference protein fasta file", metavar = "character"),
    make_option(
        opt_str = c("-n", "--novel_fasta"),
        type = "character", default = NULL,
        help = "ORF fasta file", metavar = "character"),
    make_option(
        opt_str = c("-rr", "--reciprocal_blast_ref"),
        type = "character", default = NULL,
        help = "Reciprocal blast against reference", metavar = "character"),
    make_option(
        opt_str = c("-ru", "--reciprocal_blast_uniprot"),
        type = "character", default = NULL,
        help = "Reciprocal blast against all uniprot", metavar = "character"),
    make_option(
        opt_str = c("-o", "--output"),
        type = "character", default = "", 
        help = "Output file name [default= %default]", metavar = "character"))

# Parse the parameters provided on command line by user
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# Check whether inputs parameter was provided
if (is.null(opt$fasta) | is.null(opt$blast) | is.null(opt$genomic_position)) {
    
    print_help(opt_parser)
    stop(paste(
        "The  input arguments must be supplied",
        "(fasta, blast and genomic position files)!"))
    
}

# Check whether output parameter was provided
if (opt$output == "") {
    
    warning("Output results to !")
    
}



### Data import ----------------------------------------------------------

# Import the evidence file containing peptide group info
evid <- opt$evidence %>%
    as.character(.) %>%
    readRDS(file = .)

# Import the peptide location within proteins
pep_loc <- opt$peptide_location %>%
    as.character(.) %>%
    readRDS(file = .)

# Import the fasta files
fasta <- c(Known = opt$reference_fasta, Novel = opt$novel_fasta) %>%
    purrr::map(
        .x = ., .f = seqinr::read.fasta, seqtype = "AA", as.string = TRUE)

# Import the reciprocal blast best hits between novel proteins and reference
# proteins
reciprocal_blast_ref <- opt$reciprocal_blast_ref %>%
    as.character(.) %>%
    read.table(file = ., header = TRUE, sep = "\t", quote = "")

# Import the reciprocal blast best hits between novel proteins and all
# uniprot proteins
reciprocal_blast_uniprot <- opt$reciprocal_blast_uniprot %>%
    as.character(.) %>%
    read.table(file = ., header = TRUE, sep = "\t", quote = "")



### Levenshtein distance -------------------------------------------------

# Keep only sequence for novel peptide
peptides <- evid %>%
    dplyr::filter(., group == "Novel") %>%
    .[["Sequence"]] %>%
    unique(.)

# Compute the levenshtein distance for all novel peptide and keep
# the minimum leven score result per peptide
leven_data <- adist(
    x = peptides,
    y = fasta$Known,
    partial = TRUE,
    ignore.case = TRUE) %>%
    set_rownames(peptides) %>%
    t(.) %>%
    base::data.frame(
        id = rownames(.),
        .,
        stringsAsFactors = FALSE) %>%
    tidyr::gather(data = ., key = "Sequence", value = "leven", -id) %>%
    dplyr::group_by(., Sequence) %>%
    dplyr::filter(., leven == min(leven)) %>%
    base::as.data.frame(., stringsAsFActors = FALSE)



