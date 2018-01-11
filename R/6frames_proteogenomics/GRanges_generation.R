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
            "R/6frames_proteogenomics/helper.R",
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
load_package("GenomicRanges")



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
        geno_name = readline(
            prompt = "Provide the genome name!"),
        circular = readline(
            prompt = "Is the chromosome circular (logical)?") %>%
            as.logical(.),
        annotation = choose.files(
            caption = "Choose an annotation file!",
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
            opt_str = c("-n", "--geno_name"),
            type = "character", default = "",
            help = "Genome name",
            metavar = "character"),
        make_option(
            opt_str = c("-t", "--circular"),
            type = "logical", default = "",
            help = "Genome type is circular?",
            metavar = "character"),
        make_option(
            opt_str = c("-a", "--annotation"),
            type = "character", default = "",
            help = "Entries annotation file (first column must be the key)",
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
    identical(opt$geno_name, "") |
    identical(opt$geno_name, character(0))) {
    
    print_help(opt_parser)
    stop(paste(
        "Genome name is required!"))
    
}
if (
    identical(opt$circular, "") |
    identical(opt$circular, logical(0))) {
    
    print_help(opt_parser)
    stop(paste(
        "Chromosome circular or linear (logical) required!"))
    
}
if (
    identical(opt$annotation, "") |
    identical(opt$annotation, character(0))) {
    
    opt["annotation"] <- list(NULL)
    warning(paste(
        "No annotation provided, skip these steps!!"))
    
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

# Import the fasta files
genome <- seqinr::read.fasta(
    file = opt$genome, seqtype = "DNA", as.string = TRUE)

# Check whether coordinates data are present
if (!is.null(opt$coordinates)) {
    coordinates <- data.table::fread(
        input = opt$coordinates, sep = "\t", header = TRUE,
        stringsAsFactors = FALSE, quote = "")
}

# Check whether annotation data are present
if (!is.null(opt$annotation)) {
    annotation <- data.table::fread(
        input = opt$annotation, sep = "\t", header = TRUE,
        stringsAsFactors = FALSE, quote = "")
}



### Ideogram generation --------------------------------------------------

# Get genome length
geno_size <- seqinr::getLength(genome)

# Chromosomes information based on genome fasta sequences
chromos <- base::data.frame(
    id = names(genome),
    chromosome = names(genome),
    strand = "*",
    start = 1,
    end = geno_size,
    type = opt$circular,
    genome = opt$geno_name,
    stringsAsFactors = FALSE)

# Check whether coordinates are present
if (!exists("coordinates")) {
    
    # Data holding genomic ranges information for the genome only
    grange_data <- chromos %>%
        dplyr::select(., -type, -genome)
    
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
    
    # Check whether annotation are present
    if (exists("annotation")) {
        
        # Rename the first column name from annotation to match the one fro
        # the coordinates
        colnames(annotation)[1] <- colnames(grange_data)[1]
        
        # Check that the column renaming did not create duplicate
        if (any(duplicated(colnames(annotation)))) {
            stop("Duplicated column name in annotation after renaming!")
        }
        
        # Include the annotation to the coordinates
        grange_data %<>%
            dplyr::left_join(x = ., y = annotation)
        
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
        chromosome = as.character(chromosome),
        start = as.integer(start),
        end = as.integer(end)) %>%
    base::as.data.frame(., stringsAsFactors = FALSE)

# Create a GRanges object for all entries with coordinates
grange <- with(
    data = grange_data,
    expr = GRanges(
        seqnames = chromosome,
        ranges = IRanges(start, end),
        strand = strand))

# Add seqinfo to the created GRanges object
seqinfo(grange) <- Seqinfo(
    seqnames = chromos$chromosome %>% as.character(.),
    seqlengths = chromos$end %>% as.integer(.),
    isCircular = chromos$type %>% as.logical(.),
    genome = chromos$genome %>% as.character(.))

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


