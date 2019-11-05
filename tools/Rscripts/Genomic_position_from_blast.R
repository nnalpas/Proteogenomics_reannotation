#!/usr/bin/env Rscript

# This script determines the protein genomic coordinates based on
# blast results against genome



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
load_package("seqinr")
load_package("Biostrings")



### Parameters setting up ------------------------------------------------

# Define input parameters (interactively or from command line)
if (interactive()) {
    
    # Define the list of input parameters
    opt <- list()
    opt["fasta"] <- choose.files(
        caption = "Choose Fasta file containing entries to get coordinates",
        multi = FALSE) %>%
        list(.)
    opt["blast"] <- choose.files(
        caption = "Choose Blast results from which to get coordinates",
        multi = FALSE) %>%
        list(.)
    opt["genome"] <- choose.files(
        caption = "Choose Fasta file of the genome",
        multi = FALSE) %>%
        list(.)
    opt["output"] <- "protein_location.txt"
    
} else {
    
    # Define the list of command line parameters
    option_list <- list(
        make_option(
            opt_str = c("-f", "--fasta"),
            type = "character", default = NULL,
            help = "Fasta file",
            metavar = "character"),
        make_option(
            opt_str = c("-b", "--blast"),
            type = "character", default = NULL,
            help = "Blast file",
            metavar = "character"),
        make_option(
            opt_str = c("-g", "--genome"),
            type = "character", default = NULL,
            help = "Genome fasta file",
            metavar = "character"),
        make_option(
            opt_str = c("-o", "--output"),
            type = "character", default = "protein_location.txt", 
            help = "Output file name [default= %default]",
            metavar = "character"))
    
    # Parse the parameters provided on command line by user
    opt_parser <- OptionParser(option_list = option_list)
    opt <- parse_args(opt_parser)
    
}

# Check whether inputs parameter were provided
if (is.null(opt$fasta)) {
    
    print_help(opt_parser)
    stop(paste(
        "Fasta file required!"))
    
}
if (is.null(opt$blast)) {
    
    print_help(opt_parser)
    stop(paste(
        "Blast results required!"))
    
}
if (is.null(opt$genome)) {
    
    print_help(opt_parser)
    stop(paste(
        "Genome fasta file required!"))
    
}

# Check whether output parameter was provided
if (opt$output == "protein_location.txt") {
    
    opt$output <- paste0("./", date_str, "_protein_location.txt")
    warning(paste0(
        "Output results to '",
        opt$output,
        "'!"))
    
}

# Create output directory if not already existing
dir.create(dirname(opt$output))



### Data import ----------------------------------------------------------

# Import the fasta file
fasta <- opt$fasta %>%
    as.character(.) %>%
    seqinr::read.fasta(
        file = ., seqtype = "AA", as.string = TRUE) %>%
    set_names(uni_id_clean(names(.)))

# Import the fasta file
geno_size <- opt$genome %>%
    as.character(.) %>%
    seqinr::read.fasta(
        file = ., seqtype = "DNA", as.string = TRUE) %>%
    seqinr::getLength(.)

# Import the Blast results
blast_data <- blast_read(file = opt$blast, blast_format = "6") %>%
    dplyr::mutate(
        .,
        qseqid = uni_id_clean(qseqid),
        sseqid = uni_id_clean(sseqid))



### Compute strand, frame and genomic coordinate per entry ---------------

# Get the best blast and format sequence as stringset
blast_data_best <- best_blast(
    data = blast_data, key = "qseqid", multi_match = "uniquify") %>%
    dplyr::rowwise() %>%
    dplyr::mutate(
        .,
        qseq = list(Biostrings::AAString(qseq)),
        sseq = list(Biostrings::AAString(sseq)))

# Define genomic strand
blast_data_best %<>%
    dplyr::mutate(
        .,
        ExactCoord = ifelse(
            qstart == 1 & qend == qlen,
            "Exact",
            "Approximate"),
        strand = ifelse(
            test = sstart < send, yes = "+", no = "-"))

# Check that all entries in fasta file have a blast results
missing_entries <- names(fasta)[
    !names(fasta) %in% unique(blast_data_best$qseqid)]
warning(paste0(
    "There are ", length(missing_entries), " entries without blast!\n",
    "Please check the missing following entries and manually",
    " record coordinates if necessary:\n",
    paste(missing_entries, collapse = "\n")))

# Define translation frame
blast_data_best <- get_frame(
    data = blast_data_best,
    start = "sstart",
    strand = "strand",
    genome_size = geno_size)

# Adjust genomic coordinate if blast is not from start to end of query
blast_data_best %<>%
    dplyr::mutate(
        .,
        start = ifelse(
            strand == "+",
            (sstart - ((qstart - 1) * 3)),
            (sstart + ((qstart - 1) * 3))),
        end = ifelse(
            strand == "+",
            (send + ((qlen - qend) * 3)),
            (send - ((qlen - qend) * 3))),
        nucl_length_qlen = qlen * 3,
        nucl_length_str_end = abs(end - start) + 1,
        nucl_length = nucl_length_str_end,
        Comment = ifelse(
            nucl_length_qlen != nucl_length_str_end, "Warning", "OK")) %>%
    base::as.data.frame(., stringsAsFactors = FALSE)

# Keep only required columns (related to genomic coordinates)
data_final <- blast_data_best %>%
    dplyr::select(
        ., qseqid, sseqid, strand, frame, start, end,
        nucl_length, qlen, ExactCoord, Comment) %>%
    set_colnames(c(
        "id", "chromosome", "strand", "frame", "start", "end",
        "nucl_length", "aa_length", "ExactCoord", "Comment"))



### Results export -------------------------------------------------------

# Export the reference entries coordinates data
write.table(
    x = data_final, file = opt$output,
    quote = FALSE, sep = "\t",
    row.names = FALSE, col.names = TRUE)

# Define end time
print(paste("Complete", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))


