#!/usr/bin/env Rscript

# This script retrieve column annotations based on specific keys using
# uniprot.ws package



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
load_package("seqinr")
load_package("UniProt.ws")



### Parameters setting up ------------------------------------------------

# Define input parameters (interactively or from command line)
if (interactive()) {
    
    # Define the list of input parameters
    opt <- list(
        fasta = choose.files(
            caption = "Choose Fasta file of the known protein (UniProt)!",
            multi = FALSE),
        taxon = readline(
            prompt = "Provide the taxon identifier!"),
        columns = readline(
            prompt = "Give comma-separated annotation names to retrieve!"),
        key = readline(
            prompt = "Provide the key name to use for cross-annotation!"),
        output = readline(
            prompt = "Define the output file name!"))
    
} else {
    
    # Define the list of command line parameters
    option_list <- list(
        make_option(
            opt_str = c("-f", "--fasta"),
            type = "character", default = NULL,
            help = "Known protein fasta file",
            metavar = "character"),
        make_option(
            opt_str = c("-t", "--taxon"),
            type = "integer", default = NULL,
            help = "The taxon identifier for the species of interest",
            metavar = "character"),
        make_option(
            opt_str = c("-c", "--columns"),
            type = "character", default = NULL,
            help = "The annotation columns (comma-separated) to retrieve!",
            metavar = "character"),
        make_option(
            opt_str = c("-k", "--key"),
            type = "character", default = NULL,
            help = "The key name to use for cross-annotation!",
            metavar = "character"),
        make_option(
            opt_str = c("-o", "--output"),
            type = "character", default = "", 
            help = "Output file name [default= %default]",
            metavar = "character"))
    
    # Parse the parameters provided on command line by user
    opt_parser <- OptionParser(option_list = option_list)
    opt <- parse_args(opt_parser)
    
}

# Check whether inputs parameter were provided
if (
    is.null(opt$fasta) | is.null(opt$taxon) |
    is.null(opt$columns) | is.null(opt$key)) {
    
    print_help(opt_parser)
    
}

# Store the annotation column name into vector
opt$columns %<>%
    strsplit(x = ., split = ",", fixed = TRUE) %>%
    unlist(.)

# Check whether output parameter was provided
if (opt$output == "") {
    
    opt$output <- sub(
        pattern = "^(.*)(\\..*)?",
        replacement = "\\1_annot.txt",
        x = opt$fasta)
    warning(paste0(
        "Output results to '",
        opt$output,
        "'!"))
    
}

# Create output directory if not already existing
dir.create(dirname(opt$output))



### Data import ----------------------------------------------------------

# Import the fasta files
fasta <- seqinr::read.fasta(
    file = opt$fasta, seqtype = "AA", as.string = TRUE)



### Obtain protein annotation from UniProt -------------------------------

# Clean up the UniProtKB
names(fasta) <- uni_id_clean(names(fasta))

# Create UniProt.ws object for submitted taxon ID
up <- UniProt.ws(taxId = as.numeric(as.character(opt$taxon)))

# Get all the columns specified by user
data <- UniProt.ws::select(
    x = up, keys = names(fasta), columns = opt$columns, keytype = opt$key)

# Extract the gene name and locus from the GENES column if existing
if (any(colnames(data) %in% c("GENES"))) {
    data_final <- data %>%
        dplyr::mutate(
            .,
            `GENE-NAME` = sub(
                "^(.*?) .*$", "\\1", GENES, perl = TRUE),
            `ALT-GENE-NAME` = sub(
                "^.*? (.*)$", "\\1", GENES, perl = TRUE))
} else {
    data_final <- data
}



### Results export -------------------------------------------------------

# Export the bsu ideogram for reuse at later stage
write.table(
    x = data_final,
    file = opt$output,
    quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)

# Define end time
print(paste("Complete", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))


