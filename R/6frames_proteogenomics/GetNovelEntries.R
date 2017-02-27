#!/usr/bin/Rscript

# This script identify all novel peptides and ORFs from a MaxQuant search
# as input as well as the fasta files used for the search



### Parameters setting up ------------------------------------------------

# Start with clean environment
rm(list = ls())



### Define working directory ---------------------------------------------

# Import command line arguments
args <- commandArgs(trailingOnly = TRUE)

# Check for appropriate number of command line arguments
if (length(args) != 3) {
    stop("Usage: <WorkSpace> <> <>!")
}

# Define the work space
work.space <- args[1]
setwd(work.space)

# Define the current user
user <- Sys.info()[["user"]]

# Define current time
date.str <- format(Sys.time(), "%Y-%m-%d")



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


