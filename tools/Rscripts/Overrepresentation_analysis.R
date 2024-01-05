#!/usr/bin/env Rscript

# This script performs overrepresentation analysis from a submitted fasta file



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
#if (interactive()) {
#    source(
#        file = paste(
#            "C:/Users",
#            user,
#            "Documents/GitHub/Proteogenomics_reannotation/",
#            "tools/Rscripts/helper.R",
#            sep = "/"))
#} else {
#    source(
#        file = paste(
#            Sys.getenv("HOME"),
#            "bin/helper.R",
#            sep = "/"))
#}

# Load the required packages
library(magrittr)
library(Easy)



### Parameters setting up ------------------------------------------------

opt <- list(
    annotation = "C:/Users/nalpanic/SynologyDrive/Work/Abaumannii_trimeth/Annotation/2023-05-16/Acinetobacter_baumannii_ATCC_17978_full_annotation_2023-05-16.txt",
    foreground = "C:/Users/nalpanic/SynologyDrive/Work/Colleagues shared work/Brandon_Robin/Abaumannii_mutants/Analysis/Condition_explanation/Succinyl_proteins_for_OA.txt",
    background = NULL,
    resource = "COG_function,GOBP Tree Term,GOMF Tree Term,GOCC Tree Term,KEGG Pathway Name,Virulence Database,EC level 1 name,EC level 2 name,EC level 3 name,InterPro Description,Other family Description,Subcellular Localization [b2g]",
    gene = "Locus Tag",
    idcol = "Accessions ABYAL",
    substitute = c("\\|.*", ""),
    pval = 1,
    padj = 1,
    minsize = 1,
    maxsize = Inf,
    threads = 1,
    output = "C:/Users/nalpanic/SynologyDrive/Work/Colleagues shared work/Brandon_Robin/Abaumannii_mutants/Analysis/Functional_annotation/")

# Check whether inputs parameters were provided
if (
    is.null(opt$output) | is.null(opt$annotation) |
    is.null(opt$foreground)) {
    
    print_help(opt_parser)
    stop("All arguments must be supplied!")
    
}
if (
    identical(opt$threads, NULL) |
    identical(opt$threads, "") |
    identical(opt$threads, integer(0))) {
    
    warning("Default number of threads will be used!")
    opt$threads <- 1
    
}
doParallel::registerDoParallel(cores = opt$threads)
print(paste("Number of threads registered:", foreach::getDoParWorkers()))

# Create output directory if not already existing
dir.create(opt$output)



### Data import ----------------------------------------------------------

my_annotation <- data.table::fread(
    input = opt$annotation, sep = "\t", quote = "", header = TRUE,
    stringsAsFactors = FALSE, colClasses = "character")

if (!is.na(opt$idcol) & !is.null(opt$idcol) & opt$idcol != "") {
    my_ranking <- data.table::fread(
        input = opt$foreground, sep = "\t", quote = "",
        header = TRUE, stringsAsFactors = FALSE)
    check_cols <- colnames(my_ranking)[!colnames(my_ranking) %in% opt$idcol]
    for (a in check_cols) {
        if (
            any(!my_ranking[[a]] %in% c(0, 1, NA_integer_)) &
            any(!my_ranking[[a]] %in% c(FALSE, TRUE, NA)) &
            any(!my_ranking[[a]] %in% c("FALSE", "TRUE", NA)) &
            any(!my_ranking[[a]] %in% c("No", "Yes", NA)) &
            any(!my_ranking[[a]] %in% c("", "+", NA))) {
            my_ranking[[a]] <- NULL
        }
    }
    my_ranking_long <- my_ranking %>%
        tidyr::pivot_longer(data = ., cols = -as.name(opt$idcol)) %>%
        dplyr::filter(., !is.na(value) & (value == 1 | isTRUE(value) | value == "TRUE" | value == "Yes" | value == "+")) %>%
        dplyr::mutate(., name = as.factor(name)) %>%
        unique(.)
    my_foreground <- split(
        x = my_ranking_long[[opt$idcol]], f = my_ranking_long$name)
} else {
    my_foreground <- seqinr::read.fasta(
        file = opt$foreground, seqtype = "AA", as.string = TRUE) %>%
        names(.) %>%
        list(.) %>%
        set_names(basename(opt$foreground))
}

if (!is.null(opt$background)) {
    my_background <- seqinr::read.fasta(
        file = opt$background, seqtype = "AA", as.string = TRUE)
} else {
    my_background <- unique(my_ranking[[opt$idcol]]) %>%
        set_names(unique(my_ranking[[opt$idcol]]))
}

if (!is.null(opt$substitute)) {
    my_annotation[[opt$gene]] %<>%
        sub(opt$substitute[1], opt$substitute[2], .)
    my_foreground %<>%
        lapply(X = ., function(x) sub(opt$substitute[1], opt$substitute[2], x))
    names(my_background) %<>%
        sub(opt$substitute[1], opt$substitute[2], .)
}

my_resource <- opt$resource %>%
    base::strsplit(x = ., split = ",", fixed = TRUE) %>%
    unlist(.) %>%
    unique(.)



### Perform overrepresentation analysis ----------------------------------

# Loop through selected resources and perform separate
# overrepresentation analysis
my_oa_combined <- lapply(X = my_resource, FUN = function(x) {
    #for (x in my_resource) {

    pathways <- fgsea_pathways(
        annotation = my_annotation, gene = opt$gene, resource = x)
    
    my_results <- lapply(my_foreground, function(stats) {
        fora_scored(
            pathways = pathways, genes = stats,
            universe = names(my_background), minSize = opt$minsize,
            maxSize = opt$maxsize, pval = opt$pval, padj = opt$padj)
    }) %>%
        plyr::ldply(., data.table::data.table, .id = "set")
    
}) %>%
    set_names(my_resource) %>%
    plyr::ldply(., dplyr::bind_rows, .id = "resource")



### Export the results data ----------------------------------------------

# Export annotation usable in Perseus
write.table(
    x = my_oa_combined,
    file = paste0(
        opt$output, "/", date_str, "_", basename(opt$foreground), ".OA.txt"),
    quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)

# Save session
save.image(paste0(
    opt$output, "/", date_str, "_", basename(opt$foreground), ".OA.RData"))



### END ------------------------------------------------------------------

# Time scripts end
print(paste("Completed ", format(Sys.time(), "%Y-%m-%d"), sep = ""))


