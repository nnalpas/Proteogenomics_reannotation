#!/usr/bin/env Rscript

# This script identify all novel peptides and ORFs from a MaxQuant search
# as input as well as the fasta files used for the search



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
# Note: the select function from dplyr and UniProt.ws is overlapping,
# therefore order of package loading is important
load_package(plyr)
load_package(dplyr)
load_package(seqinr)
load_package(bit64)
load_package(magrittr)
load_package(data.table)
load_package(splitstackshape)
load_package(stringr)
load_package(optparse)



### Parameters setting up ------------------------------------------------

# Define the list of command line parameters
option_list <- list(
    make_option(
        opt_str = c("-o", "--output"),
        type = "character", default = NULL,
        help = "Output directory", metavar = "character"),
    make_option(
        opt_str = c("-m", "--maxquant"),
        type = "character", default = NULL,
        help = "MaxQuant txt results folder", metavar = "character"),
    make_option(
        opt_str = c("-r", "--reference"),
        type = "character", default = NULL,
        help = "Reference protein fasta file", metavar = "character"),
    make_option(
        opt_str = c("-n", "--novel"),
        type = "character", default = NULL,
        help = "ORF fasta file", metavar = "character"))

# Parse the parameters provided on command line by user
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# Check whether inputs parameter was provided
if (
    is.null(opt$output) | is.null(opt$maxquant) |
    is.null(opt$reference) | is.null(opt$novel)) {
    
    print_help(opt_parser)
    stop("All arguments must be supplied!")
    
}

# For manual parameters set-up
#opt <- list(
#    output = "C:/Users/kxmna01/Dropbox/Home_work_sync/Work/Colleagues shared work/Vaishnavi Ravikumar/Bacillus_subtilis_6frame",
#    maxquant = "G:/data/Vaishnavi/combined - 6 frame translation/txt",
#    reference = "C:/Users/kxmna01/Dropbox/Home_work_sync/Work/Colleagues shared work/Vaishnavi Ravikumar/Bacillus_subtilis_6frame/Databases/uniprot-proteome_Bacillus_subtilis_168_UP000001570_20150318.fasta",
#    novel = "C:/Users/kxmna01/Dropbox/Home_work_sync/Work/Colleagues shared work/Vaishnavi Ravikumar/Bacillus_subtilis_6frame/Databases/Bsu_genome_assembly_GCA_000009045.1.out_FIXED_HEADER.fasta")



### Data import ----------------------------------------------------------

# Import the maxquant evidence table
evid <- mq_read(
    path = opt$maxquant,
    name = "evidence.txt",
    integer64 = "double")

# Import the fasta files
fasta <- c(Known = opt$reference, Novel = opt$novel) %>%
    purrr::map(
        .x = ., .f = seqinr::read.fasta, seqtype = "AA", as.string = TRUE)
names(fasta$Known) %<>%
    uni_id_clean(.)



### Novel peptide identification -----------------------------------------

# Format association of peptide to protein
data <- evid %>%
    cSplit(indt = ., splitCols = "Proteins", sep = ";", direction = "long") %>%
    dplyr::select(., Sequence, Proteins) %>%
    unique(.) %>%
    base::as.data.frame(., stringsAsFactors = FALSE)

# Format the fasta files
names(fasta_list$Known) %<>%
    sub(pattern = ".+\\|(.+)\\|.+", replacement = "\\1", x = .)

# Locate all the peptide within associated proteins
pep_loc <- list()
pep_loc[["Known"]] <- pept.locate(
    data = data, peptide = "Sequence",
    proteins = "Proteins", fasta = fasta_list$Known) %>%
    dplyr::filter(., !is.na(start))
pep_loc[["Novel"]] <- pept.locate(
    data = data, peptide = "Sequence",
    proteins = "Proteins", fasta = fasta_list$Novel) %>%
    dplyr::filter(., !is.na(start))
pep_loc[["Contaminant"]] <- data %>%
    dplyr::filter(., grepl("^CON__", Proteins)) %>%
    set_colnames(c("pep", "prot")) %>%
    dplyr::mutate(., start = NA_integer_, end = NA_integer_)
pep_loc[["Reverse"]] <- evid %>%
    dplyr::filter(., Reverse == "+") %>%
    dplyr::select(., Sequence, Proteins) %>%
    set_colnames(c("pep", "prot")) %>%
    dplyr::mutate(., start = NA_integer_, end = NA_integer_)

# New dataframe to hold info about fasta of origin for each sequence
evid_match <- evid %>%
    dplyr::mutate(
        .,
        group = ifelse(
            test = Sequence %in% pep_loc$Known$pep,
            yes = "Known",
            no = ifelse(
                test = Sequence %in% pep_loc$Contaminant$pep,
                yes = "Contaminant",
                no = ifelse(
                    test = Sequence %in% pep_loc$Novel$pep,
                    yes = "Novel",
                    no = ifelse(
                        test = Sequence %in% pep_loc$Reverse$pep,
                        yes = "Reverse",
                        no = NA_character_))))) %>%
    base::as.data.frame(., stringAsFactors = TRUE)

# Print warning for unidentified sequence origin
print(paste(
    "There are", length(which(is.na(evid_match$group))),
    "NA values, these need to be checked!", sep = " "))

# Print the repartition of peptide per group
print(table(evid_match$group, useNA = "always"))

# Save the group mapping data
saveRDS(
    object = evid_match,
    file = "Sequence_group_mapping.RDS")

# Save the peptide location data
saveRDS(
    object = pep_loc,
    file = "Peptides_location.RDS")

# Export complete evidence info for these novel evidence
write.table(
    x = evid_match[evid_match$group == "Novel", ],
    file = "Novel_evidence.txt",
    quote = FALSE,
    sep = "\t",
    row.names = FALSE,
    col.names = TRUE)

# Export complete evidence info for these novel evidence
write.table(
    x = evid_match %>%
        dplyr::filter(., group == "Novel") %>%
        cSplit(
            indt = ., splitCols = "Proteins", sep = ";",
            direction = "long", fixed = TRUE) %>%
        .[["Proteins"]] %>%
        unique(.),
    file = "Novel_ORF.txt",
    quote = FALSE,
    sep = "\t",
    row.names = FALSE,
    col.names = FALSE)



### END ------------------------------------------------------------------

# Time scripts end
print(paste("Completed ", format(Sys.time(), "%Y-%m-%d"), sep = ""))


