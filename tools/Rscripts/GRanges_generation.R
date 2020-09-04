#!/usr/bin/env Rscript

# This script generates an ideogram for the genome of interest



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
library(plyr)
library(dplyr)
library(magrittr)
library(data.table)
library(splitstackshape)
library(stringr)
library(optparse)
library(seqinr)
library(GenomicRanges)



### Parameters setting up ------------------------------------------------

# Define input parameters (interactively or from command line)
if (interactive()) {
    
    # Define the list of input parameters
    opt <- list(
        coordinates = choose.files(
            caption = "Choose entries coordinate file!",
            multi = FALSE),
        genome = choose.files(
            caption = "Choose Fasta file of the genome!",
            multi = FALSE),
        bsgenome = readline(
            prompt = "Provide the BSgenome package name") %>%
            as.character(.),
        annotations = choose.files(
            caption = "Choose an annotations file!",
            multi = FALSE),
        output = readline(
            prompt = "Define the output filename (with .RDS extension)!"))
    
} else {
    
    # Define the list of command line parameters
    option_list <- list(
        make_option(
            opt_str = c("-c", "--coordinates"),
            type = "character", default = "",
            help = "Entries coordinate file (first column must be the key)",
            metavar = "character"),
        make_option(
            opt_str = c("-g", "--genome"),
            type = "character", default = "",
            help = "Genome fasta file",
            metavar = "character"),
        make_option(
            opt_str = c("-b", "--bsgenome"),
            type = "logical", default = "",
            help = "Provide the BSgenome package name",
            metavar = "character"),
        make_option(
            opt_str = c("-a", "--annotations"),
            type = "character", default = "",
            help = "Entries annotations file (first column must be the key)",
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

# Check whether input parameters were provided
if (
    identical(opt$coordinates, "") |
    identical(opt$coordinates, character(0))) {
    
    opt["coordinates"] <- list(NULL)
    warning(paste(
        "No coordinates provided, skip these steps!!"))
    
}
if (
    identical(opt$genome, "") |
    identical(opt$genome, character(0))) {
    
    print_help(opt_parser)
    stop(paste(
        "Genome fasta file required!"))
    
}
if (
    identical(opt$bsgenome, NULL) |
    identical(opt$bsgenome, "") |
    identical(opt$bsgenome, character(0))) {
    
    print_help(opt_parser)
    stop("The input BSgenome package must be supplied!")
    
}
warning(paste("Loading package", opt$bsgenome))
library(
    package = eval(opt$bsgenome),
    character.only = TRUE)

if (
    identical(opt$annotations, "") |
    identical(opt$annotations, character(0))) {
    
    opt["annotations"] <- list(NULL)
    warning(paste(
        "No annotations provided, skip these steps!!"))
    
}

# Check whether output parameter was provided
if (
    identical(opt$output, "") |
    identical(opt$output, character(0)) |
    length(grep("\\.RDS$", opt$output)) == 0) {
    
    print_help(opt_parser)
    stop(paste("Output filename is required (with .RDS extension)!"))
    
}

# Create output directory if not already existing
dir.create(dirname(opt$output))



### Data import ----------------------------------------------------------

# Import the fasta files (the genome is not used anymore, can be removed in future)
genome <- seqinr::read.fasta(
    file = opt$genome, seqtype = "DNA", as.string = TRUE)

# Check whether coordinates data are present
if (!is.null(opt$coordinates)) {
    coordinates <- data.table::fread(
        input = opt$coordinates, sep = "\t", header = TRUE,
        stringsAsFactors = FALSE, quote = "")
}

# Check whether annotations data are present
if (!is.null(opt$annotations)) {
    annotations <- data.table::fread(
        input = opt$annotations, sep = "\t", header = TRUE,
        stringsAsFactors = FALSE, quote = "")
}



### Ideogram generation --------------------------------------------------

# Get the BSgenome object
bsgeno <- eval(parse(text = opt$bsgenome))

# Chromosomes information based on genome fasta sequences
chromos <- seqinfo(bsgeno) %>%
    as.data.frame(.) %>%
    tibble::rownames_to_column(.data = ., var = "id") %>%
    dplyr::mutate(., chromosome = id, strand = "*", start = 1) %>%
    dplyr::select(
        ., id, chromosome, strand, start,
        end = seqlengths, isCircular, genome) %>%
    dplyr::arrange(., chromosome)

# Check whether coordinates are present
if (!exists("coordinates")) {
    
    # Data holding genomic ranges information for the genome only
    grange_data <- chromos %>%
        dplyr::select(., -isCircular, -genome)
    
} else {
    
    # Check that required columns are present in the coordinates data
    if (!all(c(
        "id", "strand", "chromosome", "start", "end") %in% colnames(
            coordinates))) {
        stop(paste0(
            'Columns: "id", "strand", "chromosome", "start", "end" are',
            ' missing from coordinates'))
    }
    
    # Data holding genomic ranges information for submitted entries
    grange_data <- coordinates
    
    # Check whether annotations are present
    if (exists("annotations")) {
        
        # Rename the first column name from annotations to match the one fro
        # the coordinates
        colnames(annotations)[1] <- colnames(grange_data)[1]
        
        # Check that the column renaming did not create duplicate
        if (any(duplicated(colnames(annotations)))) {
            stop("Duplicated column name in annotations after renaming!")
        }
        
        # Include the annotations to the coordinates
        grange_data %<>%
            dplyr::left_join(x = ., y = annotations)
        
    }
    
    # Format the start-end position so that start is always inferior to end
    grange_data %<>%
        dplyr::mutate(., start_tmp = start, end_tmp = end) %>%
        dplyr::rowwise(.) %>%
        dplyr::mutate(
            .,
            start = min(start_tmp, end_tmp),
            end = max(start_tmp, end_tmp))
    
    # Check whether entries on + strand had their start-end inverted
    if (
        grange_data %>%
        dplyr::filter(., strand == "+" & start_tmp != start) %>%
        nrow(.) %>%
        is_greater_than(., 0)) {
        
        warning("Positive strand entries had their start-end inverted!")
        
    }
    
    # Clean up the grange data
    grange_data %<>%
        dplyr::select(., -end_tmp, -start_tmp)
    
}

# Format columns to their right class (numeric)
grange_data %<>%
    dplyr::mutate(
        .,
        id = as.character(id),
        strand = as.character(strand),
        chromosome = as.character(chromosome) %>%
            sub(".*\\|(.+)\\|.*", "\\1", .),
        start = as.integer(start),
        end = as.integer(end)) %>%
    dplyr::arrange(., chromosome) %>%
    base::as.data.frame(., stringsAsFactors = FALSE)

# Create a GRanges object for all entries with coordinates
grange <- with(
    data = grange_data,
    expr = GRanges(
        seqnames = chromosome,
        ranges = IRanges(start, end),
        strand = strand))

# Add seqinfo to the created GRanges object
seqlevels(grange) <- seqlevels(bsgeno) %>% as.character(.)
seqinfo(grange) <- seqinfo(bsgeno)

# Add values to the created GRanges object
values(grange) <- grange_data %>%
    dplyr::select(
        ., -chromosome, -start, -end, -strand)

# Define names of the GRange entries
names(grange) <- grange_data$id



### Results export -------------------------------------------------------

# Export the grange for reuse at later stage
saveRDS(
    object = grange,
    file = opt$output)

# Define end time
print(paste("Complete", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))


