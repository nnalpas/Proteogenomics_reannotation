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
            type = "character", default = NULL,
            help = "Entries coordinate file",
            metavar = "character"),
        make_option(
            opt_str = c("-g", "--genome"),
            type = "character", default = NULL,
            help = "Genome fasta file",
            metavar = "character"),
        make_option(
            opt_str = c("-n", "--geno_name"),
            type = "character", default = NULL,
            help = "Genome name",
            metavar = "character"),
        make_option(
            opt_str = c("-a", "--annotation"),
            type = "character", default = NULL,
            help = "Entries annotation file",
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
if (is.null(opt$genome)) {
    
    print_help(opt_parser)
    stop(paste(
        "Genome fasta file required!"))
    
}
if (is.null(opt$geno_name)) {
    
    print_help(opt_parser)
    stop(paste(
        "Genome name is required!"))
    
}

# Check whether output parameter was provided
if (opt$output == "" | length(grep("\\.RDS$", opt$output)) == 0) {
    
    print_help(opt_parser)
    stop(paste("Output filename is required (with .RDS extension)!"))
    
}

# Create output directory if not already existing
dir.create(dirname(opt$output))



### Data import ----------------------------------------------------------

# Import the fasta files
genome <- seqinr::read.fasta(
    file = opt$genome, seqtype = "DNA", as.string = TRUE)



### Ideogram generation --------------------------------------------------

# Get genome length
geno_size <- seqinr::getLength(genome)

# Dataframe holding genome information for bacillus subtilis
data <- base::data.frame(
    Chromosome = names(genome),
    Strand = "*",
    Start = 1,
    End = geno_size,
    name = names(genome),
    length = geno_size,
    Type = TRUE,
    geno = opt$geno_name,
    stringsAsFactors = FALSE)

# Create an Ideogram (GRanges object)
grange <- with(
    data = data,
    expr = GRanges(
        seqnames = name,
        ranges = IRanges(Start, End),
        strand = Strand,
        Chromosome = Chromosome,
        seqinfo = Seqinfo(
            seqnames = name,
            seqlengths = length,
            isCircular = Type,
            genome = geno)))
names(grange) <- data$Chromosome



### Results export -------------------------------------------------------

# Export the bsu grange for reuse at later stage
saveRDS(
    object = grange,
    file = opt$output)

# Define end time
print(paste("Complete", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))


