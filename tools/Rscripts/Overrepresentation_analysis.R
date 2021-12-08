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

# Load the required packages
library(magrittr)



### Parameters setting up ------------------------------------------------

opt <- list(
    annotation = "H:/data/Synechocystis_6frame/EggnogMapper/EggNOG_annotation_2_perseus.txt",
    foreground = "H:/data/Synechocystis_6frame/Phylostratigraphy/Phylostrata_for_OA.txt",
    background = "H:/data/Synechocystis_6frame/Phylostratigraphy/1148.faa",
    resource = "best_og_Category,best_og_Subcategory,GOBP Term,GOCC Term,GOMF Term,EC level 1 name,EC level 2 name,EC level 3 name,KEGG_Pathway_Name,KEGG_Module_Name,KEGG_Reaction_Name,KEGG_rclass_Name,KEGG_brite_Name,PFAMs,CAZy,BiGG_Reaction,Custom_annotation",
    gene = "#query_name",
    idcol = "qseqid",
    pval = 1,
    padj = 1,
    minsize = 1,
    maxsize = Inf,
    threads = 1,
    output = "H:/data/Synechocystis_6frame/OA_GSEA")

# Check whether inputs parameters were provided
if (
    is.null(opt$output) | is.null(opt$annotation) |
    is.null(opt$foreground) | is.null(opt$background)) {
    
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
        header = TRUE, stringsAsFactors = FALSE, colClasses = "character") %>%
        tidyr::pivot_longer(data = ., cols = -as.name(opt$idcol)) %>%
        dplyr::filter(., !is.na(value) & value == 1) %>%
        dplyr::mutate(., name = as.factor(name))
    my_foreground <- split(
        x = my_ranking[[opt$idcol]], f = my_ranking$name)
} else {
    my_foreground <- seqinr::read.fasta(
        file = opt$foreground, seqtype = "AA", as.string = TRUE) %>%
        names(.) %>%
        list(.) %>%
        set_names(basename(opt$foreground))
}

my_background <- seqinr::read.fasta(
    file = opt$background, seqtype = "AA", as.string = TRUE)

my_resource <- opt$resource %>%
    base::strsplit(x = ., split = ",", fixed = TRUE) %>%
    unlist(.) %>%
    unique(.)



### Perform overrepresentation analysis ----------------------------------

# Loop through selected resources and perform separate
# overrepresentation analysis
my_oa_combined <- lapply(X = my_resource, FUN = function(x) {

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
    file = paste0(opt$output, "/", basename(opt$foreground), ".OA.txt"),
    quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)

# Save session
save.image(paste0(opt$output, "/", basename(opt$foreground), ".OA.RData"))



### END ------------------------------------------------------------------

# Time scripts end
print(paste("Completed ", format(Sys.time(), "%Y-%m-%d"), sep = ""))


