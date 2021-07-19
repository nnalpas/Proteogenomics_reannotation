#!/usr/bin/env Rscript

# This script incorporates the results of interpro functional annotation



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

# Load the required packages (or install if not already in library)
library(magrittr)
library(optparse)
library(foreach)
library(doParallel)
library(ggplot2)



### Parameters setting up ------------------------------------------------

# Define input parameters (interactively or from command line)
if (interactive()) {
    
    # Define the list of input parameters
    opt <- list(
        novel_reason = choose.files(
            caption = "Choose orf novelty (.RDS) file!",
            multi = FALSE),
        interpro = choose.dir(caption = ""),
        ref_grange = choose.files(
            caption = "Choose reference entries grange (.RDS) file!",
            multi = FALSE),
        orf_grange = choose.files(
            caption = "Choose ORF entries grange (.RDS) file!",
            multi = FALSE),
        pep_grange = choose.files(
            caption = "Choose peptide entries grange (.RDS) file!",
            multi = FALSE),
        genome_grange = choose.files(
            caption = "Choose genome entry grange (.RDS) file!",
            multi = FALSE),
        bsgenome = readline(prompt = paste0(
            "Provide name of the BSgenome package")) %>%
            as.character(),
        pep_class = readline(prompt = paste0(
            "Which PEP class to keep",
            " (separated by space; e.g. 'class 1')")) %>%
            as.character(),
        threads = readline(prompt = "How many cores to use?") %>% as.integer(),
        output = readline(
            prompt = "Define the output directory!"))
    
} else {
    
    # Define the list of command line parameters
    option_list <- list(
        make_option(
            opt_str = c("-i", "--novel_reason"),
            type = "character", default = NULL,
            help = "Peptide novelty info",
            metavar = "character"),
        make_option(
            opt_str = c("-p", "--interpro"),
            type = "character", default = NULL,
            help = "Interpro annotation folder",
            metavar = "character"),
        make_option(
            opt_str = c("-r", "--ref_grange"),
            type = "character", default = NULL,
            help = "Reference entries grange file",
            metavar = "character"),
        make_option(
            opt_str = c("-n", "--orf_grange"),
            type = "character", default = NULL,
            help = "ORF entries grange file",
            metavar = "character"),
        make_option(
            opt_str = c("-d", "--pep_grange"),
            type = "character", default = NULL,
            help = "Peptide entries grange file",
            metavar = "character"),
        make_option(
            opt_str = c("-g", "--genome_grange"),
            type = "character", default = NULL,
            help = "Genome entry grange file",
            metavar = "character"),
        make_option(
            opt_str = c("-b", "--bsgenome"),
            type = "character", default = NULL,
            help = "The BSgenome package name for genomic visualisation",
            metavar = "character"),
        make_option(
            opt_str = c("-c", "--pep_class"),
            type = "character", default = NULL,
            help = "Which PEP class to keep (separated by space; e.g. 'class 1')",
            metavar = "character"),
        make_option(
            opt_str = c("-t", "--threads"),
            type = "integer", default = NULL,
            help = "Number of cores to use",
            metavar = "character"),
        make_option(
            opt_str = c("-o", "--output"),
            type = "character", default = NULL,
            help = "Output directory", metavar = "character"))
    
    # Parse the parameters provided on command line by user
    opt_parser <- OptionParser(option_list = option_list)
    opt <- parse_args(opt_parser)
    
}

# Check whether inputs parameter was provided
if (
    identical(opt$novel_reason, NULL) |
    identical(opt$novel_reason, "") |
    identical(opt$novel_reason, character(0))) {
    
    print_help(opt_parser)
    stop("The input peptide novelty must be supplied!")
    
}
if (
    identical(opt$interpro, NULL) |
    identical(opt$interpro, "") |
    identical(opt$interpro, character(0))) {
    
    print_help(opt_parser)
    warning("No interpro annotation folder provided by user!")
    
}
if (
    identical(opt$ref_grange, NULL) |
    identical(opt$ref_grange, "") |
    identical(opt$ref_grange, character(0))) {
    
    print_help(opt_parser)
    stop("The input reference entries grange must be supplied!")
    
}
if (
    identical(opt$orf_grange, NULL) |
    identical(opt$orf_grange, "") |
    identical(opt$orf_grange, character(0))) {
    
    print_help(opt_parser)
    stop("The input ORF entries grange must be supplied!")
    
}
if (
    identical(opt$pep_grange, NULL) |
    identical(opt$pep_grange, "") |
    identical(opt$pep_grange, character(0))) {
    
    print_help(opt_parser)
    stop("The input peptide entries grange must be supplied!")
    
}
if (
    identical(opt$genome_grange, NULL) |
    identical(opt$genome_grange, "") |
    identical(opt$genome_grange, character(0))) {
    
    print_help(opt_parser)
    stop("The input genome entry grange must be supplied!")
    
}
if (
    identical(opt$pep_class, NULL) |
    identical(opt$pep_class, "") |
    identical(opt$pep_class, character(0))) {
    
    opt["pep_class"] <- list("class 1")
    warning("No PEP class provided by user, default to 'class1'!")
    
}
if (
    identical(opt$threads, NULL) |
    identical(opt$threads, "") |
    identical(opt$threads, integer(0))) {
    
    warning("Default number of threads will be used!")
    
}
registerDoParallel(cores = opt$threads)
print(paste("Number of threads registered:", getDoParWorkers()))

# Load the BSgenome
if (
    identical(opt$bsgenome, NULL) |
    identical(opt$bsgenome, "") |
    identical(opt$bsgenome, character(0))) {
    
    print_help(opt_parser)
    stop("The input BSgenome package must be supplied!")
    
}
library(
    package = eval(opt$bsgenome),
    character.only = TRUE)

# Check whether output parameter was provided
if (
    identical(opt$output, NULL) |
    identical(opt$output, "") |
    identical(opt$output, character(0))) {
    
    opt$output <- dirname(opt$novel_reason)
    warning(paste0("Output results to ", opt$output, "!"))
    
}

# Create output directory if not already existing
dir.create(opt$output)



### Data import ----------------------------------------------------------

# Import the novelty reason file containing peptide novelty info
orf_reason_final <- opt$novel_reason %>%
    as.character(.) %>%
    readRDS(file = .)

# Import all interpro annotation files
my_interpro <- opt$interpro %>%
    list.files(path = ., pattern = "*.tsv$", full.names = TRUE) %>%
    lapply(., function(x) {
        data.table::fread(
            input = x, sep = "\t", quote = "", header = TRUE,
            stringsAsFactors = FALSE, colClasses = "character", fill = TRUE)
    }) %>%
    plyr::ldply(.data = ., .fun = data.table::data.table, .id = NULL)

# Import the genome genomic ranges file
genome_grange <- opt$genome_grange %>%
    as.character(.) %>%
    readRDS(file = .)

# Import the reference genomic ranges file
ref_grange <- opt$ref_grange %>%
    as.character(.) %>%
    readRDS(file = .)

# Import the ORF genomic ranges file
orf_grange <- opt$orf_grange %>%
    as.character(.) %>%
    readRDS(file = .)

# Import the peptides genomic ranges file
pep_grange <- opt$pep_grange %>%
    as.character(.) %>%
    readRDS(file = .)



### Interpro function annotation analysis --------------------------------

# List to hold all plots
my_plots <- list()

# Replace empty entries with "-"
my_interpro[my_interpro == ""] <- "-"

# 
my_interpro[my_interpro$`Protein accession` %in% names(ref_grange), "Database"] <- "Target"
my_interpro[my_interpro$`Protein accession` %in% names(orf_grange), "Database"] <- "Novel"

for (x in c("Overall", unique(my_interpro$Analysis))) {
    
    target_analysis <- x
    if (x == "Overall") {
        target_analysis <- unique(my_interpro$Analysis)
    }
    
    toplot <- my_interpro %>%
        dplyr::filter(., Analysis %in% target_analysis) %>%
        dplyr::group_by(., `Protein accession`, Database) %>%
        dplyr::summarise(., Nmb_annot = dplyr::n()) %>%
        dplyr::ungroup(.) %>%
        dplyr::mutate(
            ., Nmb_annot = ifelse(Nmb_annot >= 20, "20+", Nmb_annot)) %>%
        plyr::ddply(
            .data = ., .variables = c("Database", "Nmb_annot"),
            .fun = dplyr::summarise,
            Freq = dplyr::n_distinct(`Protein accession`),
            .drop = FALSE)
    toplot$Nmb_annot <- factor(
        x = toplot$Nmb_annot, levels = c(1:19, "20+"), ordered = TRUE)
    
    my_plots[[paste("Interpro distrib", x)]] <- ggplot(
        toplot, aes(
            x = Nmb_annot, y = Freq, fill = Database, colour = Database)) +
        geom_bar(
            stat = "identity", position = "dodge", alpha = 0.3) +
        ggpubr::theme_pubr() +
        theme(legend.position = "right") +
        facet_grid(rows = vars(Database), scales = "free_y") +
        ggtitle(paste(x, "distribution of interpro annotation"))
    
}

my_interpro_novel <- orf_reason_final %>%
    dplyr::select(., Proteins, OnlyIdBySite, PEPfilter, ORFNoveltyReason) %>%
    dplyr::left_join(
        x = ., y = my_interpro, by = c("Proteins" = "Protein accession"))

for (x in c("Overall", unique(my_interpro$Analysis))) {
    
    target_analysis <- x
    if (x == "Overall") {
        target_analysis <- unique(my_interpro$Analysis)
    }
    
    toplot <- my_interpro_novel %>%
        dplyr::filter(., !is.na(Analysis) & Analysis %in% target_analysis) %>%
        dplyr::group_by(., ORFNoveltyReason) %>%
        dplyr::summarise(., Count = dplyr::n_distinct(Proteins)) %>%
        dplyr::ungroup(.) %>%
        dplyr::mutate(., Type = "All annotated")
    
    toplot <- my_interpro_novel %>%
        dplyr::filter(., OnlyIdBySite & PEPfilter %in% opt[["pep_class"]]) %>%
        dplyr::filter(., !is.na(Analysis) & Analysis %in% target_analysis) %>%
        dplyr::group_by(., ORFNoveltyReason) %>%
        dplyr::summarise(., Count = dplyr::n_distinct(Proteins)) %>%
        dplyr::ungroup(.) %>%
        dplyr::mutate(., Type = "Quality filtered annotated") %>%
        dplyr::bind_rows(toplot, .)
    
    my_plots[[paste("Annotated novel", x)]] <- ggplot(
        toplot, aes(
            x = ORFNoveltyReason, y = Count, fill = Type, colour = Type)) +
        geom_bar(
            stat = "identity", position = "dodge") +
        ggpubr::theme_pubr() +
        theme(legend.position = "right") +
        coord_flip() +
        ggtitle(paste(x, "novel ORF with interpro annotation"))
    
}






