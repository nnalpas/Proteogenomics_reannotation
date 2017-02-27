#!/usr/bin/Rscript

# This script identify all novel peptides and ORFs from a MaxQuant search
# as input as well as the fasta files used for the search



### Parameters setting up ------------------------------------------------

# Start with clean environment
rm(list = ls())

# Define current time
date.str <- format(Sys.time(), "%Y-%m-%d")
print(paste("Start ", date.str, sep = ""))



### Define working directory ---------------------------------------------

# Import command line arguments
args <- commandArgs(trailingOnly = TRUE)

# Check for appropriate number of command line arguments
if (length(args) != 4) {
    stop("Usage: <WorkSpace> <MaxQuantSearch> <RefProteome> <6FrameProteome>!")
}

# Define the work space
work.space <- args[1]
setwd(work.space)
print(paste("Working directory: ", work.space, sep = ""))

# Select the maxquant txt folder
txt.dir <- args[2]

# List the fasta files that need to be imported
fasta.file <- c(args[3], args[4])
names(fasta.file) <- c("Known", "Novel")

# Define the current user
user <- Sys.info()[["user"]]



### List of required packages -----------------------------------------------

# Source the custom user functions
source(
    file = paste(
        "C:/Users",
        user,
        "Documents/GitHub/Miscellaneous/R/General/General_function.R",
        sep = "/"))
source(
    file = paste(
        "/Media/sf_C_DRIVE/Users",
        user,
        "Documents/GitHub/Miscellaneous/R/General/General_function.R",
        sep = "/"))

# Load the required packages (or install if not already in library)
# Note: the select function from dplyr and UniProt.ws is overlapping,
# therefore order of package loading is important
loadpackage(plyr)
loadpackage(dplyr)
loadpackage(seqinr)
loadpackage(bit64)
loadpackage(magrittr)
loadpackage(data.table)
loadpackage(splitstackshape)
loadpackage(stringr)



### Data import ----------------------------------------------------------

# Import the maxquant evidence table
evid <- maxquant.read(
    path = txt.dir,
    name = "evidence.txt",
    integer64 = "double")

# Import all fasta file data and store into list
fasta.list <- list()
for (x in 1:length(fasta.file)) {
    
    # Import the current fasta file
    tmp <- read.fasta(
        file = fasta.file[x], seqtype = "AA", as.string = TRUE)
    
    # Include imported fasta into the list
    fasta.list[names(fasta.file)[x]] <- list(tmp)
    
}

# Clean-up
rm(fasta.file)



### Novel peptide identification -----------------------------------------

# Format association of peptide to protein
data <- evid %>%
    cSplit(indt = ., splitCols = "Proteins", sep = ";", direction = "long") %>%
    dplyr::select(., Sequence, Proteins) %>%
    unique(.) %>%
    base::as.data.frame(., stringsAsFactors = FALSE)

# Format the fasta files
names(fasta.list$Known) %<>%
    sub(pattern = ".+\\|(.+)\\|.+", replacement = "\\1", x = .)

# Locate all the peptide within associated proteins
pep.loc <- list()
pep.loc[["Known"]] <- pept.locate(
    data = data, peptide = "Sequence",
    proteins = "Proteins", fasta = fasta.list$Known) %>%
    dplyr::filter(., !is.na(start))
pep.loc[["Novel"]] <- pept.locate(
    data = data, peptide = "Sequence",
    proteins = "Proteins", fasta = fasta.list$Novel) %>%
    dplyr::filter(., !is.na(start))
pep.loc[["Contaminant"]] <- data %>%
    dplyr::filter(., grepl("^CON__", Proteins)) %>%
    set_colnames(c("pep", "prot")) %>%
    dplyr::mutate(., start = NA_integer_, end = NA_integer_)
pep.loc[["Reverse"]] <- evid %>%
    dplyr::filter(., Reverse == "+") %>%
    dplyr::select(., Sequence, Proteins) %>%
    set_colnames(c("pep", "prot")) %>%
    dplyr::mutate(., start = NA_integer_, end = NA_integer_)

# New dataframe to hold info about fasta of origin for each sequence
evid.match <- evid %>%
    dplyr::mutate(
        .,
        group = ifelse(
            test = Sequence %in% pep.loc$Known$pep,
            yes = "Known",
            no = ifelse(
                test = Sequence %in% pep.loc$Contaminant$pep,
                yes = "Contaminant",
                no = ifelse(
                    test = Sequence %in% pep.loc$Novel$pep,
                    yes = "Novel",
                    no = ifelse(
                        test = Sequence %in% pep.loc$Reverse$pep,
                        yes = "Reverse",
                        no = NA_character_))))) %>%
    base::as.data.frame(., stringAsFactors = TRUE)

# Print warning for unidentified sequence origin
print(paste(
    "There are", length(which(is.na(evid.match$group))),
    "NA values, these need to be checked!", sep = " "))

# Print the repartition of peptide per group
print(table(evid.match$group, useNA = "always"))

# Save the group mapping data
saveRDS(
    object = evid.match,
    file = "Sequence_group_mapping.RDS")

# Save the peptide location data
saveRDS(
    object = pep.loc,
    file = "Peptides_location.RDS")

# Export complete evidence info for these novel evidence
write.table(
    x = evid.match[evid.match$group == "Novel", ],
    file = "Novel_evidence.txt",
    quote = FALSE,
    sep = "\t",
    row.names = FALSE,
    col.names = TRUE)

# Export complete evidence info for these novel evidence
write.table(
    x = evid.match %>%
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


