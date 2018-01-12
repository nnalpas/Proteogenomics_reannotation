#!/usr/bin/env Rscript

# This script determines the novelty reasons for each ORF based on
# genomic coordinates, blast results...



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
load_package("bit64")
load_package("ggplot2")
load_package("gtable")
load_package("grid")
load_package("gridExtra")
load_package("purrr")
load_package("foreach")
load_package("doParallel")



### Parameters setting up ------------------------------------------------

# Define input parameters (interactively or from command line)
if (interactive()) {
    
    # Define the list of input parameters
    opt <- list(
        novel_reason = choose.files(
            caption = "Choose peptide novelty (.RDS) file!",
            multi = FALSE),
        ref_grange = choose.files(
            caption = "Choose reference entries grange (.RDS) file!",
            multi = FALSE),
        orf_grange = choose.files(
            caption = "Choose ORF entries grange (.RDS) file!",
            multi = FALSE),
        operon_grange = choose.files(
            caption = "Choose operon entries grange (.RDS) file!",
            multi = FALSE),
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
            opt_str = c("-p", "--operon_grange"),
            type = "character", default = NULL,
            help = "Operon entries grange file",
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
    identical(opt$operon_grange, NULL) |
    identical(opt$operon_grange, "") |
    identical(opt$operon_grange, character(0))) {
    
    opt["operon_grange"] <- list(NULL)
    warning("Input operon not supplied, operon processes ignored!")
    
}
if (
    identical(opt$threads, NULL) |
    identical(opt$threads, "") |
    identical(opt$threads, integer(0))) {
    
    warning("Default number of threads will be used!")
    
}
registerDoParallel(cores = opt$threads)
print(paste("Number of threads registered:", getDoParWorkers()))

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
evid_reason <- opt$novel_reason %>%
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

# Import the operon genomic ranges file if defined
if (!is.null(opt$operon)) {
    operon_grange <- opt$operon %>%
        as.character(.) %>%
        readRDS(file = .)
}



### ORF neighbour analysis -----------------------------------------------

# Identify the putative novel ORFs
targets <- evid_reason %>%
    dplyr::filter(., Database == "Novel") %>%
    .[["Proteins"]] %>%
    strsplit(x = ., split = ";", fixed = TRUE) %>%
    unlist(.) %>%
    unique(.)

# Remove target IDs not present in genomic range
targets_filt <- targets[targets %in% names(orf_grange)]

# Define the subset of novel ORF that were found expressed
sub_orf_grange <- orf_grange[targets_filt]

# Identify all reference entries that are neighbour of the novel ORFs
neighbours_list <- gr_near_dist(
    query = sub_orf_grange,
    subject = ref_grange,
    near_type = c("nearest","precede", "follow"),
    select = "all")

# Reformat the neighbour list to dataframe
neighbours_cat <- neighbours_list %>%
    plyr::ldply(.data = ., .fun = "data.frame", .id = "Neighbour")

# Filter the neighbour for in-frame (distance is multiple of 3) and
# close proximity (inferior to 100 bp)
neighbours_inframe <- neighbours_cat %>%
    dplyr::filter(., (dist %% 3) == 0 & dist <= 100)

# Interprete the results of the neighbouring entries analysis
neighbours_analysis <- neighbours_inframe %>%
    dplyr::mutate(., interpr = dplyr::case_when(
        queryStrand == "+" & dist <= 3 &
            Neighbour == "precede" ~ "Five neighbour",
        queryStrand == "+" & dist <= 3 &
            Neighbour == "follow" ~ "Three neighbour",
        queryStrand == "-" & dist <= 3 &
            Neighbour == "precede" ~ "Three neighbour",
        queryStrand == "-" & dist <= 3 &
            Neighbour == "follow" ~ "Five neighbour",
        queryStrand == "+" & dist > 3 &
            Neighbour == "precede" ~ "Nearby Five neighbour",
        queryStrand == "+" & dist > 3 &
            Neighbour == "follow" ~ "Nearby Three neighbour",
        queryStrand == "-" & dist > 3 &
            Neighbour == "precede" ~ "Nearby Three neighbour",
        queryStrand == "-" & dist > 3 &
            Neighbour == "follow" ~ "Nearby Five neighbour",
        TRUE ~ "Unexplained"))



### Data visualisation and export ----------------------------------------

# Loop through each neighbour type
for (x in names(neighbours_list)) {
    
    # Generate frequency of neighbouring entries distance (between ORF and ref)
    toplot <- neighbours_list[[x]] %>%
        #dplyr::mutate(
        #    .,
        #    dist = factor(
        #        x = dist,
        #        levels = seq(0, max(dist), by = 1),
        #        labels = seq(0, max(dist), by = 1),
        #        ordered = TRUE)) %>%
        plyr::ddply(
            .data = .,
            .variables = .(queryStrand, dist),
            .fun = function(piece) {
                summarise(piece, freq = length(queryID)) },
            .drop = FALSE,
            .parallel = TRUE)
    
    # Visualise frequency as histogram
    pl <- plots_hist(
        data = toplot,
        key = "dist",
        value = "freq",
        group = "dist",
        fill = "queryStrand",
        posit = "stack",
        main = paste("Frequency of distance between", x, "entries"),
        xlabel = "Distance (bp)",
        ylabel = "Count",
        legend = "bottom",
        bw = TRUE,
        xdir = "vertical",
        auto_scale = 10)
    print(pl[[1]])
    pl <- plots_hist(
        data = toplot %>% dplyr::filter(., as.numeric(dist) < 101),
        key = "dist",
        value = "freq",
        group = "dist",
        fill = "queryStrand",
        posit = "stack",
        main = paste("Frequency of  0 - 100 bp between", x, "entries"),
        xlabel = "Distance (bp)",
        ylabel = "Count",
        legend = "bottom",
        bw = TRUE,
        xdir = "vertical")
    print(pl[[1]])
    
}



### END ------------------------------------------------------------------

# Close the cluster
stopImplicitCluster()

# Time scripts end
print(paste("Completed ", format(Sys.time(), "%Y-%m-%d"), sep = ""))


