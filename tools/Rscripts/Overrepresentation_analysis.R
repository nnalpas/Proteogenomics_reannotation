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
    annotation = "H:/data/Synechocystis_6frame/Custom_annotation/2021-12-21_Custom_Uniprot_Eggnog_annotations.txt",
    foreground = "H:/data/Synechocystis_6frame/2022-02-14_iBAQ/Scy004_resuscitation_Scy001_iBAQ_norm_t-test_forOA.txt",
    background = "H:/data/Synechocystis_6frame/Phylostratigraphy/1148.faa",
    resource = "Miscellaneous,TU ID,Custom_classification,Active site-note,Binding site-note,Catalytic activity-Reaction,ChEBI,Pathway,Site-note,Keywords,Protein existence,Status,Developmental stage,Induction,Tissue specificity,Subcellular location [CC],Intramembrane,Topological domain-note,Transmembrane-note,Post-translational modification,Modified residue-note,Propeptide-id,Signal peptide-note,Transit peptide,Protein families,Domain [FT]-note,Motif-note,Characterization,EC level 1 name,EC level 2 name,EC level 3 name,GOBP Term,GOCC Term,GOMF Term,Interpro_NAME,Panther Name,Panther Protein class,Panther Pathway,PIRSF name,Prosite,Prosite DE,TIGRFAM label,TIGRFAM product_name,best_og_Category,best_og_Subcategory,GOBP Term.EggNOG,GOCC Term.EggNOG,GOMF Term.EggNOG,EC level 1 name.EggNOG,EC level 2 name.EggNOG,EC level 3 name.EggNOG,KEGG_Pathway_Name,KEGG_Module_Name,KEGG_Reaction_Name,KEGG_rclass_Name,KEGG_brite_Name,pfam_id,pfam_description,pfam_clan_id,pfam_clan_description,CAZy,BiGG_Reaction",
    gene = "#query_name",
    idcol = "Proteins",
    pval = 1,
    padj = 1,
    minsize = 1,
    maxsize = Inf,
    threads = 1,
    output = "H:/data/Synechocystis_6frame/2022-02-14_iBAQ")

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

my_date <- Sys.Date()

# Export annotation usable in Perseus
write.table(
    x = my_oa_combined,
    file = paste0(
        opt$output, "/", my_date, "_", basename(opt$foreground), ".OA.txt"),
    quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)

# Save session
save.image(paste0(
    opt$output, "/", my_date, "_", basename(opt$foreground), ".OA.RData"))



### END ------------------------------------------------------------------

# Time scripts end
print(paste("Completed ", format(Sys.time(), "%Y-%m-%d"), sep = ""))


