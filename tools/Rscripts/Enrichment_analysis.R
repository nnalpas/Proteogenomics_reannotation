#!/usr/bin/env Rscript

# This script performs gene set enrichment analysis from a submitted file



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
    annotation = "H:/data/0989LaAn/eggnog/EggNOG_annotation_2_perseus.txt",
    ranking = "H:/data/0989LaAn/Ranalysis/2021-04_Perseus_analysis/Mix1_and_2_iBAQ.txt",
    resource = "Preferred_name,Custom_Pathway,best_og_Category,best_og_Subcategory,GOBP Term,GOCC Term,GOMF Term,EC level 1 name,EC level 2 name,EC level 3 name,KEGG_Pathway_Name,KEGG_Module_Name,KEGG_Reaction_Name,KEGG_rclass_Name,KEGG_brite_Name,PFAMs,CAZy,BiGG_Reaction",
    gene = "#query_name",
    leading = as.logical(TRUE),
    nozero = as.logical(TRUE),
    idcol = "Protein IDs",
    pval = 1,
    padj = 1,
    minsize = 1,
    maxsize = Inf,
    threads = 1,
    output = "H:/data/0989LaAn/Ranalysis/2021-11-18_OA_GSEA")

# Check whether inputs parameters were provided
if (
    is.null(opt$output) | is.null(opt$annotation) |
    is.null(opt$ranking)) {
    
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

my_ranking <- data.table::fread(
    input = opt$ranking, sep = "\t", quote = "",
    header = TRUE, stringsAsFactors = FALSE, colClasses = "character")

my_resource <- opt$resource %>%
    base::strsplit(x = ., split = ",", fixed = TRUE) %>%
    unlist(.) %>%
    unique(.)



### Perform gene-set enrichment analysis ---------------------------------

# Check format of the ranking data
if (any(colnames(my_ranking) != opt$idcol)) {
    
    my_cols <- colnames(my_ranking)[colnames(my_ranking) != opt$idcol]
    if (opt$leading) {
        my_ranking[[opt$idcol]] %<>%
            sub(";.+", "", .)
    } else {
        my_ranking %<>%
            tidyr::separate_rows(data = ., opt$idcol, sep = ";")
    }
    
    stats_list <- lapply(my_cols, function (x) {
        stats <- my_ranking %>%
            dplyr::select(., !!as.name(opt$idcol), !!as.name(x))
        stats[[x]] <- as.double(stats[[x]])
        stats %<>%
            dplyr::filter(., !is.na(!!as.name(x)))
        if (opt$nozero) {
            stats %<>%
                dplyr::filter(., !!as.name(x) != 0)
        }
        stats[[x]] %>%
            set_names(stats[[opt$idcol]])
    }) %>%
        set_names(my_cols)
    
}

# Loop through selected resources and perform separate
# gene-set enrichment analysis
my_gsea_combined <- lapply(X = my_resource, FUN = function(x) {
    
    pathways <- fgsea_pathways(
        annotation = my_annotation, gene = opt$gene, resource = x)
    
    my_results <- lapply(stats_list, function(stats) {
        fgsea_scored(
            pathways = pathways, stats = stats,
            minSize = opt$minsize, maxSize = opt$maxsize,
            pval = opt$pval, padj = opt$padj)
    }) %>%
        plyr::ldply(., data.table::data.table, .id = "set")
    
}) %>%
    set_names(my_resource) %>%
    plyr::ldply(., dplyr::bind_rows, .id = "resource")



### Export the results data ----------------------------------------------

# Export annotation usable in Perseus
write.table(
    x = my_gsea_combined,
    file = paste0(opt$output, "/", basename(opt$ranking), ".GSEA.txt"),
    quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)

# Save session
save.image(paste0(opt$output, "/", basename(opt$ranking), ".GSEA.RData"))



### END ------------------------------------------------------------------

# Time scripts end
print(paste("Completed ", format(Sys.time(), "%Y-%m-%d"), sep = ""))

