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
library(plyr)
library(dplyr)
library(magrittr)
library(data.table)
library(splitstackshape)
library(stringr)
library(optparse)
library(seqinr)
library(bit64)
library(ggplot2)
library(ggbio)
library(cowplot)
library(gtable)
library(grid)
library(gridExtra)
library(purrr)
library(foreach)
library(doParallel)
#library(motifRG)
#library(ggseqlogo)



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
        pep_grange = choose.files(
            caption = "Choose peptide entries grange (.RDS) file!",
            multi = FALSE),
        sanger_grange = choose.files(
            caption = "Choose Sanger seq entries grange (.RDS) file!",
            multi = FALSE),
        genome_grange = choose.files(
            caption = "Choose genome entry grange (.RDS) file!",
            multi = FALSE),
        pep_class = readline(prompt = paste0(
            "Which PEP class to keep",
            " (separated by comma; e.g. 'class 1')")) %>%
            as.character(),
        bsgenome = readline(prompt = paste0(
            "Provide name of the BSgenome package")) %>%
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
            opt_str = c("-d", "--pep_grange"),
            type = "character", default = NULL,
            help = "Peptide entries grange file",
            metavar = "character"),
        make_option(
            opt_str = c("-s", "--sanger_grange"),
            type = "character", default = NULL,
            help = "Sanger seq entries grange file",
            metavar = "character"),
        make_option(
            opt_str = c("-g", "--genome_grange"),
            type = "character", default = NULL,
            help = "Genome entry grange file",
            metavar = "character"),
        make_option(
            opt_str = c("-c", "--pep_class"),
            type = "character", default = NULL,
            help = "Which PEP class to keep (separated by comma; e.g. 'class 1')",
            metavar = "character"),
        make_option(
            opt_str = c("-b", "--bsgenome"),
            type = "character", default = NULL,
            help = "The BSgenome package name for genomic visualisation",
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
    identical(opt$pep_grange, NULL) |
    identical(opt$pep_grange, "") |
    identical(opt$pep_grange, character(0))) {
    
    print_help(opt_parser)
    stop("The input peptide entries grange must be supplied!")
    
}
if (
    identical(opt$sanger_grange, NULL) |
    identical(opt$sanger_grange, "") |
    identical(opt$sanger_grange, character(0))) {
    
    opt["sanger_grange"] <- list(NULL)
    warning("No Sanger seq entries grange provided by user!")
    
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
evid_reason <- opt$novel_reason %>%
    as.character(.) %>%
    readRDS(file = .)

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

# Import the Sanger seq genomic ranges file
if (
    identical(opt$sanger_grange, NULL) ||
    identical(opt$sanger_grange, "") ||
    identical(opt$sanger_grange, character(0)) ||
    !file.exists(opt$sanger_grange)) {
    sanger_grange <- GRanges()
    seqinfo(sanger_grange) <- seqinfo(genome_grange)
} else {
    sanger_grange <- opt$sanger_grange %>%
        as.character(.) %>%
        readRDS(file = .)
}

# Import the operon genomic ranges file if defined
if (!is.null(opt$operon) && file.exists(opt$operon)) {
    operon_grange <- opt$operon %>%
        as.character(.) %>%
        readRDS(file = .)
}

# Format the PEP classes if defined
if (!is.null(opt$pep_class)) {
    opt$pep_class %<>%
        as.character(.) %>%
        strsplit(x = ., split = ",", fixed = TRUE) %>%
        unlist(.) %>%
        unique(.)
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
        dist <= 3 & Neighbour == "precede" ~ "Three neighbour",
        dist <= 3 & Neighbour == "follow" ~ "Five neighbour",
        dist > 3 & Neighbour == "precede" ~ "Nearby Three neighbour",
        dist > 3 & Neighbour == "follow" ~ "Nearby Five neighbour",
        TRUE ~ "No neighbour"))



### Operon overlap analysis ----------------------------------------------

# Check whether operon data are present
if (exists("operon_grange")) {
    
    # Identify the operon that overlap with putative novel ORF
    overlap_operon <- findOverlaps(
        query = sub_orf_grange, subject = operon_grange,
        type = "any", select = "all")
    
    # Include the entries ID and strand in the dataframe
    overlap_operon %<>%
        as.data.frame(.) %>%
        dplyr::mutate(
            .,
            queryID = names(sub_orf_grange)[queryHits],
            queryStrand = as.character(strand(sub_orf_grange)[queryHits]),
            subjectID = names(operon_grange)[subjectHits],
            subjectStrand = as.character(strand(operon_grange)[subjectHits]))
    
}



### Possible start analysis ----------------------------------------------

# Get the BSgenome object
my_bsgeno <- eval(parse(text = opt$bsgenome))

# Define number of nucleotide after first nucleotide (to define start codon)
nucleotide_after <- 2

# Create new GRange object to study nucleotide frequence of start codon
# on positive strand
plus_grange <- subset(ref_grange, strand == "+")
names(plus_grange) %<>%
    sub("^", "start_", .)
end(plus_grange) <- (start(plus_grange) + nucleotide_after)

# Create new GRange object to study nucleotide frequence of start codon
# on negative strand
minus_grange <- subset(ref_grange, strand == "-")
names(minus_grange) %<>%
    sub("^", "start_", .)
start(minus_grange) <- (end(minus_grange) - nucleotide_after)

# Compile positive and negative strand entries
start_codon <- c(plus_grange, minus_grange)

# Get the nucleotide sequence associated with the start codon
start_codon_seq <- getSeq(x = my_bsgeno, names = start_codon)

# Check frequency of start codon
start_codon_freq <- trinucleotideFrequency(
    x = start_codon_seq, step = 1, as.prob = FALSE, as.array = FALSE,
    fast.moving.side = "right", with.labels = TRUE,
    simplify.as = "collapsed")

# Reformat and calculate percentage of count per start codon
start_codon_freq %<>%
    as.data.frame(.) %>%
    set_colnames("Count") %>%
    dplyr::mutate(
        .,
        Codon = rownames(.),
        Percentage = Count / sum(Count))

# Keep only the start codon with frequency higher than 1%
used_start <- start_codon_freq %>%
    dplyr::filter(., Percentage > 0.01) %>%
    .[["Codon"]]

# Locate start codon across genome
find_start_codon <- function(x) {
    res <- Biostrings::vmatchPattern(
        pattern = DNAString(x), subject = my_bsgeno)
    mcols(res) <- data.frame(Start_codon = rep(x = x, times = length(res)))
    res
}

start_match <- foreach::foreach(
    x = used_start, .combine = c, .packages = c("Biostrings")) %dopar%
    find_start_codon(x = x)

start_match_fr <- calculate_frame(start_match)
pep_grange_fr <- calculate_frame(pep_grange)
sub_orf_grange_fr <- calculate_frame(sub_orf_grange)
#ref_grange_fr <- calculate_frame(ref_grange)

all_starts <- data.frame()
for (x in c(-3, -2, -1, 1, 2, 3)) {
    start_match_fr_sub <- subset(start_match_fr, Frame == x)
    pep_grange_fr_sub <- subset(pep_grange_fr, Frame == x) %>%
        subset(., group == "Novel")
    sub_orf_grange_fr_sub <- subset(sub_orf_grange_fr, Frame == x)
    #ref_grange_fr_sub <- subset(ref_grange_fr, Frame == x)
    
    for (y in sub_orf_grange_fr_sub$id) {
        #curr_orf <- subset(ref_grange_fr_sub, id == y)
        curr_orf <- subset(sub_orf_grange_fr_sub, id == y)
        curr_pep <- subsetByOverlaps(
            pep_grange_fr_sub,
            curr_orf,
            type = "within")
        if (x < 0) {
            curr_express <- GRanges(
                seqnames = seqnames(curr_orf),
                ranges = IRanges(max(end(curr_pep))-3, end(curr_orf)),
                strand = "-")
            curr_start <- subsetByOverlaps(
                start_match_fr_sub,
                curr_express,
                type = "within") %>%
                as.data.frame(.) %>%
                dplyr::arrange(., start)
        } else {
            curr_express <- GRanges(
                seqnames = seqnames(curr_orf),
                ranges = IRanges(start(curr_orf), min(start(curr_pep))+3),
                strand = "+")
            curr_start <- subsetByOverlaps(
                start_match_fr_sub,
                curr_express,
                type = "within") %>%
                as.data.frame(.) %>%
                dplyr::arrange(., dplyr::desc(start))
        }
        
        all_starts <- curr_start %>%
            dplyr::mutate(
                ., id = y,
                Starts = paste0(Start_codon, ":", start, "-", end)) %>%
            dplyr::group_by(., id) %>%
            dplyr::summarise(
                ., Starts = paste0(Starts, collapse = ";")) %>%
            dplyr::ungroup(.) %>%
            dplyr::bind_rows(all_starts, .)
        
    }
    
}



### Novelty reason determination -----------------------------------------

# Get all ORF to explain
orf_reason <- evid_reason %>%
    dplyr::filter(., Database %in% c("Novel", "Target")) %>%
    dplyr::select(
        ., Sequence, Proteins, PEP, `MS/MS count`, Score, Intensity,
        group, Database, OnlyIdBySite, NoveltyReason, PEPfilter,
        blast_ref, blast_best) %>%
    cSplit(
        indt = ., splitCols = "Proteins", sep = ";",
        direction = "long", fixed = TRUE) %>%
    dplyr::filter(., Proteins %in% targets_filt)

# Compile info for each novel ORFs based on all sequences (target and novel)
#orf_reason_cat <- orf_reason %>%
#    dplyr::group_by(., Proteins) %>%
#    dplyr::mutate(
#        .,
#        all_peptide_count = n_distinct(Sequence),
#        all_Sequence = paste(unique(Sequence), collapse = ";"),
#        all_min_PEP = min(PEP, na.rm = TRUE),
#        all_sum_MSMS_Count = sum(`MS/MS count`, na.rm = TRUE),
#        all_max_Score = max(Score, na.rm = TRUE),
#        all_sum_Intensity = sum(Intensity, na.rm = TRUE))

# Compile peptide info for each novel ORFs
orf_reason_cat <- orf_reason %>%
    dplyr::group_by(., Proteins, Database) %>%
    dplyr::summarise(
        .,
        peptide_count = n_distinct(Sequence),
        Sequence = paste(unique(Sequence), collapse = ";"),
        min_PEP = min(PEP, na.rm = TRUE),
        sum_MSMS_Count = sum(`MS/MS count`, na.rm = TRUE),
        max_Score = max(Score, na.rm = TRUE),
        sum_Intensity = sum(Intensity, na.rm = TRUE)) %>%
    tidyr::gather(
        data = ., key = "Param", value = "Value",
        -Proteins, -Database, convert = TRUE) %>%
    tidyr::unite(
        data = ., col = "Key", Database, Param, sep = "_", remove = TRUE) %>%
    tidyr::spread(data = ., key = "Key", value = "Value", convert = TRUE) %>%
    ungroup()

# Compile ORF database and quality metrics
orf_reason_cat <- orf_reason %>%
    dplyr::filter(., Database == "Novel") %>%
    dplyr::group_by(Proteins) %>%
    dplyr::summarise(
        group = paste(unique(group), collapse = ";"),
        Database = paste(unique(Database), collapse = ";"),
        OnlyIdBySite = ifelse(all(OnlyIdBySite == FALSE), FALSE, TRUE),
        PEPfilter = dplyr::case_when(
            any(PEPfilter == "class 1") ~ "class 1",
            any(PEPfilter == "class 2") ~ "class 2",
            any(PEPfilter == "class 3") ~ "class 3",
            TRUE ~ NA_character_)) %>%
    dplyr::left_join(x = orf_reason_cat, y = ., by = "Proteins") %>%
    dplyr::arrange(., desc(Novel_peptide_count))

# Simplify the peptide novelty by prioritising the best reason
orf_reason_clean <- orf_reason %>%
    dplyr::filter(., Database == "Novel") %>%
    dplyr::select(., Proteins, NoveltyReason) %>%
    tidyr::separate(
        data = ., col = "NoveltyReason", into = c("ref_reason", "best_reason"),
        sep = " \\(", remove = TRUE, convert = FALSE) %>%
    dplyr::group_by(Proteins, ref_reason) %>%
    dplyr::summarise(
        ., reason = dplyr::case_when(
            any(is.na(best_reason)) ~ NA_character_,
            any(best_reason == "partial match elsewhere)") ~ "partial match elsewhere)",
            any(best_reason == "SAV match elsewhere)") ~ "SAV match elsewhere)",
            any(best_reason == "exact match elsewhere)") ~ "exact match elsewhere)")) %>%
    dplyr::ungroup() %>%
    tidyr::unite(data = ., col = "NoveltyReason", ref_reason, reason, sep = " (", remove = TRUE, na.rm = TRUE) %>%
    dplyr::group_by(Proteins) %>%
    dplyr::summarise(., ORFNoveltyReason = paste(
        sort(as.character(unique(NoveltyReason))), collapse = ";")) %>%
    dplyr::left_join(x = orf_reason_cat, y = ., by = "Proteins")

# Clean-up the ORF novelty reason by removing 'Not novel' if another reason
#orf_reason_clean <- orf_reason_cat %>%
#    dplyr::mutate(., ORFNoveltyReason = sub(
#        "(;Not novel|Not novel;)", "", PepNoveltyReason))

# Clean-up the ORF novelty reason by removing 'Potential alternate start'
# if 'Alternate start (known other species)' reason
#orf_reason_clean %<>%
#    dplyr::mutate(., ORFNoveltyReason = sub(
#        "(Alternate start \\(known other species\\));Potential alternate start",
#        "\\1",
#        ORFNoveltyReason))

# Clean-up the ORF novelty reason by removing 'SAV'
# if 'SAVs/InDels' reason
#orf_reason_clean %<>%
#    dplyr::mutate(., ORFNoveltyReason = sub(
#        "SAV;(SAVs/InDels)",
#        "\\1",
#        ORFNoveltyReason))

# Include the neighbour reference entry analysis
orf_reason_neighb <- neighbours_analysis %>%
    dplyr::filter(., interpr %in% c("Five neighbour", "Three neighbour")) %>%
    dplyr::select(., queryID, subjectID, interpr) %>%
    dplyr::group_by(., queryID, interpr) %>%
    dplyr::summarise(., subjectID = paste(subjectID, collapse = ";")) %>%
    dplyr::mutate(., interpr = factor(x = interpr, levels = c("Five neighbour", "Three neighbour"))) %>%
    tidyr::complete(., interpr = interpr) %>%
    tidyr::pivot_wider(
        data = ., names_from = interpr, values_from = subjectID) %>%
    set_colnames(sub("neighbour", "neighbour (+/-3bp)", colnames(.))) %>%
    dplyr::left_join(
        x = orf_reason_clean, y = ., by = c("Proteins" = "queryID"))

if (!"Five neighbour (+/-3bp)" %in% colnames(orf_reason_neighb)) {
    orf_reason_neighb[["Five neighbour (+/-3bp)"]] <- NA_character_
}

if (!"Three neighbour (+/-3bp)" %in% colnames(orf_reason_neighb)) {
    orf_reason_neighb[["Three neighbour (+/-3bp)"]] <- NA_character_
}

# Include the reference protein ID which match novel ORF if any
orf_reason_neighb <- as.data.frame(orf_grange, stringsAsFactors = FALSE) %>%
    dplyr::filter(., !is.na(MainID)) %>%
    dplyr::select(., id, MainID) %>%
    dplyr::left_join(
        x = orf_reason_neighb, y = ., by = c("Proteins" = "id"))

# Include the operon overlap analysis if possible
if (exists("overlap_operon")) {
    
    orf_reason_opr <- overlap_operon %>%
        dplyr::select(., queryID, subjectID) %>%
        set_colnames(c("id", "OperonID")) %>%
        dplyr::group_by(., id) %>%
        dplyr::summarise(
            ., OperonID = paste(unique(OperonID), collapse = ";")) %>%
        dplyr::left_join(
            x = orf_reason_neighb, y = ., by = c("Proteins" = "id"))
    
} else {
    
    orf_reason_opr <- orf_reason_neighb %>%
        dplyr::mutate(., OperonID = NA)
    
}

orf_reason_start <- orf_reason_opr %>%
    dplyr::left_join(
        x = ., y = all_starts,
        by = c("Proteins" = "id"))

# Compute final ORF novelty reason based on five/three prime neighbour
orf_reason_final <- orf_reason_start %>% 
    dplyr::mutate(., ORFNoveltyReason = dplyr::case_when(
        grepl("Potential alternate start|Potentially novel", ORFNoveltyReason) &
            !is.na(`Five neighbour (+/-3bp)`) ~ paste(
                ORFNoveltyReason, "Erroneous termination", sep = ";"),
        TRUE ~ ORFNoveltyReason)) %>%
    dplyr::mutate(., ORFNoveltyReason = dplyr::case_when(
        grepl("Potential alternate end|Potentially novel", ORFNoveltyReason) &
            !is.na(`Three neighbour (+/-3bp)`) ~ paste(
                ORFNoveltyReason, "Erroneous start", sep = ";"),
        TRUE ~ ORFNoveltyReason))

# Include the nucleotide and amino acid length
orf_reason_final <- data.frame(
    Proteins = names(sub_orf_grange),
    nucl_length = sub_orf_grange@elementMetadata@listData$nucl_length,
    aa_length = sub_orf_grange@elementMetadata@listData$aa_length) %>%
    dplyr::left_join(x = orf_reason_final, y = .)

# Include the ORF genomic coordinates and strand
orf_reason_final <- orf_grange %>%
    as.data.frame(., stringsAsFactors = FALSE) %>%
    dplyr::select(., seqnames, start, end, strand, id) %>%
    set_colnames(c("Chromosome", "start", "end", "strand", "Proteins")) %>%
    dplyr::left_join(x = orf_reason_final, y = ., by = "Proteins")

# Include the best blast data
blast_info <- orf_reason %>%
    dplyr::filter(., Database == "Novel") %>%
    dplyr::select(., Proteins, blast_ref, blast_best) %>%
    tidyr::pivot_longer(
        data = ., names_to = "key", values_to = "value", cols = -Proteins) %>%
    tidyr::separate_rows(
        data = ., value, sep = ";", convert = FALSE) %>%
    tidyr::separate(
        data = ., col = "value",
        into = c("id", "evalue", "score", "pident", "description", "taxon"),
        sep = "\\|\\|", remove = TRUE, convert = TRUE) %>%
    dplyr::group_by(Proteins, key) %>%
    dplyr::summarise_all(~paste(unique(.), collapse = ";")) %>%
    dplyr::ungroup(.)
orf_reason_final <- blast_info %>%
    dplyr::filter(., key == "blast_best") %>%
    dplyr::select(
        ., Proteins, best_blast_id = id,
        best_blast_description = description, best_blast_taxon = taxon) %>%
    dplyr::left_join(x = orf_reason_final, y = ., by = "Proteins")



### Coverage location evidences ------------------------------------------

# Loop through all novel ORFs
for (x in orf_reason_final$Proteins) {
    pep_gr_filt <- subset(pep_grange, grepl(paste0(x, "(,|$)"), Proteins))
    orf_reason_final[orf_reason_final$Proteins == x, "Coverage_start"] <- min(
        c(start(pep_gr_filt), end(pep_gr_filt)))
    orf_reason_final[orf_reason_final$Proteins == x, "Coverage_end"] <- max(
        c(start(pep_gr_filt), end(pep_gr_filt)))
}



### Focus on high quality novel ORF --------------------------------------

# Keep high quality novel ORF by filtering out the entries only identified
# by sites and entries above PEP threshold
orf_reason_highqual <- orf_reason_final %>%
    dplyr::filter(., OnlyIdBySite & PEPfilter %in% opt[["pep_class"]])



### Data description and coverage ----------------------------------------

# Open a file for plot visualisation
pdf(
    file = paste0(opt$output, "/", "ORF_novelty_reason.pdf"),
    width = 10, height = 10)

# Histogram of evidence counts
toplot <- evid_reason %>%
    dplyr::group_by(., Database) %>%
    dplyr::summarise(., evid_count = n())
pl_datab_count <- plots_hist(
    data = toplot,
    key = "Database",
    value = "evid_count",
    group = "Database",
    fill = "Database",
    main = "PSM count (re-annotated)",
    xlabel = "Databases",
    ylabel = "Count (log scale)",
    textsize = 25,
    label = "evid_count",
    transf = "log10",
    bw = TRUE)
plot(pl_datab_count[[1]])

# Boxplot of evidence PEP
toplot <- evid_reason %>%
    dplyr::select(., Database, PEP)
pl_datab_pep <- plots_box(
    data = toplot,
    key = "Database",
    value = "PEP",
    main = "PSM PEP (re-annotated)",
    textsize = 25,
    fill = "grey",
    xlabel = "Databases",
    ylabel = "PEP",
    outlier_simplify = TRUE)
plot(pl_datab_pep[[1]])

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
            .parallel = FALSE)
    
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

# Plot the neighbour results
toplot <- neighbours_analysis %>%
    dplyr::group_by(., interpr) %>%
    dplyr::summarise(., count = n_distinct(queryID)) %>%
    dplyr::mutate(., Type = "All novel")
toplot <- neighbours_analysis %>%
    dplyr::filter(., queryID %in% unique(orf_reason_highqual$Proteins)) %>%
    dplyr::group_by(., interpr) %>%
    dplyr::summarise(., count = n_distinct(queryID)) %>%
    dplyr::mutate(., Type = "Quality filtered") %>%
    dplyr::bind_rows(toplot, .)
plots_hist(
    data = toplot,
    key = "interpr",
    value = "count",
    group = "Type",
    fill = "Type",
    main = "Reference entry neighbouring ORF",
    xlabel = "Neighbour type",
    ylabel = "Count of ORF with neighbour",
    textsize = 15,
    label = "count",
    bw = TRUE,
    legend = "bottom")

# Plot the operon overlap results
toplot <- orf_reason_final %>%
    dplyr::mutate(., within_operon = !is.na(OperonID))
toplot$within_operon <- factor(
    x = toplot$within_operon, levels = c(TRUE, FALSE), ordered = TRUE)
toplot %<>%
    plyr::ddply(
        .data = ., .variables = .(within_operon), .fun = summarise,
        count = n_distinct(Proteins), .drop = FALSE) %>%
    #dplyr::group_by(., within_operon) %>%
    #dplyr::summarise(., count = n_distinct(Proteins)) %>%
    dplyr::mutate(., Type = "All novel")
toplot_tmp <- orf_reason_final %>%
    dplyr::filter(., Proteins %in% unique(orf_reason_highqual$Proteins)) %>%
    dplyr::mutate(., within_operon = !is.na(OperonID))
toplot_tmp$within_operon <- factor(
    x = toplot_tmp$within_operon, levels = c(TRUE, FALSE), ordered = TRUE)
toplot <- toplot_tmp %>%
    plyr::ddply(
        .data = ., .variables = .(within_operon), .fun = summarise,
        count = n_distinct(Proteins), .drop = FALSE) %>%
    #dplyr::group_by(., within_operon) %>%
    #dplyr::summarise(., count = n_distinct(Proteins)) %>%
    dplyr::mutate(., Type = "Quality filtered") %>%
    dplyr::bind_rows(toplot, .)
plots_hist(
    data = toplot,
    key = "within_operon",
    value = "count",
    group = "Type",
    fill = "Type",
    main = "ORF within known operon",
    xlabel = "Within operon",
    ylabel = "Count of ORF",
    textsize = 15,
    label = "count",
    bw = TRUE,
    legend = "bottom")

# Plot the possible start results
start_stats <- orf_reason_final %>%
    dplyr::mutate(., start_type = dplyr::case_when(
        is.na(Starts) ~ "No start",
        grepl(";", Starts) ~ "Multiple start",
        TRUE ~ "single start"))
start_stats$start_type <- factor(
    x = start_stats$start_type,
    levels = c("No start", "Multiple start", "single start"),
    ordered = TRUE)
toplot <- start_stats %>%
    plyr::ddply(
        .data = ., .variables = .(start_type), .fun = summarise,
        count = n_distinct(Proteins), .drop = FALSE) %>%
    dplyr::mutate(., Type = "All novel")
toplot <- start_stats %>%
    dplyr::filter(., Proteins %in% unique(orf_reason_highqual$Proteins)) %>%
    plyr::ddply(
        .data = ., .variables = .(start_type), .fun = summarise,
        count = n_distinct(Proteins), .drop = FALSE) %>%
    dplyr::mutate(., Type = "Quality filtered") %>%
    dplyr::bind_rows(toplot, .)
plots_hist(
    data = toplot,
    key = "start_type",
    value = "count",
    group = "Type",
    fill = "Type",
    main = "ORF with start",
    xlabel = "Start presence",
    ylabel = "Count of ORF",
    textsize = 15,
    label = "count",
    bw = TRUE,
    legend = "bottom")

# Visualise ORF amino acid length repartition
toplot <- data.frame(
    Proteins = names(ref_grange),
    aa_length = ref_grange@elementMetadata@listData$aa_length) %>%
    dplyr::mutate(., Type = "Known entries")
toplot <- orf_reason_final %>%
    dplyr::select(., Proteins, aa_length) %>%
    dplyr::mutate(., Type = "All novel") %>%
    dplyr::bind_rows(toplot, .)
toplot <- orf_reason_final %>%
    dplyr::filter(., Proteins %in% unique(orf_reason_highqual$Proteins)) %>%
    dplyr::select(., Proteins, aa_length) %>%
    dplyr::mutate(., Type = "Quality filtered") %>%
    dplyr::bind_rows(toplot, .)
plots_box(
    data = toplot,
    key = "Type",
    value = "aa_length",
    main = "ORF length",
    xlabel = "Type",
    ylabel = "Amino acid length",
    colour = "black",
    fill = "grey",
    textsize = 15)

# Visualise ORF per novelty type
toplot <- orf_reason_final %>%
    dplyr::select(., Proteins, ORFNoveltyReason) %>%
    dplyr::group_by(., ORFNoveltyReason) %>%
    dplyr::summarise(., count = n_distinct(Proteins)) %>%
    dplyr::mutate(., Type = "All novel")
toplot <- orf_reason_final %>%
    dplyr::filter(., Proteins %in% unique(orf_reason_highqual$Proteins)) %>%
    dplyr::select(., Proteins, ORFNoveltyReason) %>%
    dplyr::group_by(., ORFNoveltyReason) %>%
    dplyr::summarise(., count = n_distinct(Proteins)) %>%
    dplyr::mutate(., Type = "Quality filtered") %>%
    dplyr::bind_rows(toplot, .) %>%
    dplyr::mutate(., ORFNoveltyReason = stringr::str_trunc(
        string = ORFNoveltyReason, width = 60, side = "right"))
pl_orf_reason <- plots_hist(
    data = toplot,
    key = "ORFNoveltyReason",
    value = "count",
    group = "Type",
    fill = "Type",
    main = "ORF novelty explanation",
    xlabel = "Novelty reason type",
    ylabel = "Count of ORF",
    textsize = 15,
    label = "count",
    bw = TRUE,
    legend = "bottom",
    xdir = "vertical")
plot(pl_orf_reason[[1]] + coord_flip())

# Visualise ORF per novelty type and per operon status
toplot <- orf_reason_final %>%
    dplyr::select(., Proteins, ORFNoveltyReason, OperonID) %>%
    dplyr::mutate(., within_operon = !is.na(OperonID)) %>%
    plyr::ddply(
        .data = ., .variables = .(ORFNoveltyReason, within_operon),
        .fun = summarise, count = n_distinct(Proteins), .drop = FALSE) %>%
    dplyr::mutate(., Type = "All novel")
toplot <- orf_reason_final %>%
    dplyr::filter(., Proteins %in% unique(orf_reason_highqual$Proteins)) %>%
    dplyr::select(., Proteins, ORFNoveltyReason, OperonID) %>%
    dplyr::mutate(., within_operon = !is.na(OperonID)) %>%
    plyr::ddply(
        .data = ., .variables = .(ORFNoveltyReason, within_operon),
        .fun = summarise, count = n_distinct(Proteins), .drop = FALSE) %>%
    dplyr::mutate(., Type = "Quality filtered") %>%
    dplyr::bind_rows(toplot, .) %>%
    dplyr::mutate(., ORFNoveltyReason = stringr::str_trunc(
        string = ORFNoveltyReason, width = 60, side = "right"))
pl_all_orf_reason_operon <- plots_hist(
    data = toplot %>% dplyr::filter(., Type == "All novel"),
    key = "ORFNoveltyReason",
    value = "count",
    group = "ORFNoveltyReason",
    fill = "within_operon",
    posit = "stack",
    main = "All ORF novelty with operon",
    xlabel = "Novelty reason type",
    ylabel = "Count of ORF",
    textsize = 15,
    bw = TRUE,
    legend = "bottom",
    xdir = "vertical")
plot(pl_all_orf_reason_operon[[1]] + coord_flip())
pl_hq_orf_reason_operon <- plots_hist(
    data = toplot %>% dplyr::filter(., Type == "Quality filtered"),
    key = "ORFNoveltyReason",
    value = "count",
    group = "ORFNoveltyReason",
    fill = "within_operon",
    posit = "stack",
    main = "High quality ORF novelty with operon",
    xlabel = "Novelty reason type",
    ylabel = "Count of ORF",
    textsize = 15,
    bw = TRUE,
    legend = "bottom",
    xdir = "vertical")
plot(pl_hq_orf_reason_operon[[1]] + coord_flip())

# Visualise ORF per novelty type and start results
toplot <- start_stats %>%
    plyr::ddply(
        .data = ., .variables = .(ORFNoveltyReason, start_type),
        .fun = summarise, count = n_distinct(Proteins), .drop = FALSE) %>%
    dplyr::mutate(., Type = "All novel")
toplot <- start_stats %>%
    dplyr::filter(., Proteins %in% unique(orf_reason_highqual$Proteins)) %>%
    plyr::ddply(
        .data = ., .variables = .(ORFNoveltyReason, start_type),
        .fun = summarise, count = n_distinct(Proteins), .drop = FALSE) %>%
    dplyr::mutate(., Type = "Quality filtered") %>%
    dplyr::bind_rows(toplot, .) %>%
    dplyr::mutate(., ORFNoveltyReason = stringr::str_trunc(
        string = ORFNoveltyReason, width = 60, side = "right"))
pl_all_orf_reason_start <- plots_hist(
    data = toplot %>% dplyr::filter(., Type == "All novel"),
    key = "ORFNoveltyReason",
    value = "count",
    group = "ORFNoveltyReason",
    fill = "start_type",
    posit = "stack",
    main = "All ORF novelty with start",
    xlabel = "Novelty reason type",
    ylabel = "Count of ORF",
    textsize = 15,
    bw = TRUE,
    legend = "bottom",
    xdir = "vertical")
plot(pl_all_orf_reason_start[[1]] + coord_flip())
pl_hq_orf_reason_start <- plots_hist(
    data = toplot %>% dplyr::filter(., Type == "Quality filtered"),
    key = "ORFNoveltyReason",
    value = "count",
    group = "ORFNoveltyReason",
    fill = "start_type",
    posit = "stack",
    main = "High quality ORF novelty with start",
    xlabel = "Novelty reason type",
    ylabel = "Count of ORF",
    textsize = 15,
    bw = TRUE,
    legend = "bottom",
    xdir = "vertical")
plot(pl_hq_orf_reason_start[[1]] + coord_flip())

# Frequency of ORF per number of novel peptide per novelty type
toplot <- orf_reason_final %>%
    dplyr::select(., Proteins, ORFNoveltyReason, Novel_peptide_count) %>%
    plyr::ddply(
        .data = .,
        .variables = .(ORFNoveltyReason, Novel_peptide_count),
        .fun = summarise,
        count = dplyr::n_distinct(Proteins),
        .drop = FALSE) %>%
    dplyr::mutate(., Type = "All novel")
toplot <- orf_reason_final %>%
    dplyr::filter(., Proteins %in% unique(orf_reason_highqual$Proteins)) %>%
    dplyr::select(., Proteins, ORFNoveltyReason, Novel_peptide_count) %>%
    plyr::ddply(
        .data = .,
        .variables = .(ORFNoveltyReason, Novel_peptide_count),
        .fun = summarise,
        count = dplyr::n_distinct(Proteins),
        .drop = FALSE) %>%
    dplyr::mutate(., Type = "Quality filtered") %>%
    dplyr::bind_rows(toplot, .)
pl_all_orf_freq <- plots_hist(
    data = toplot %>%
        dplyr::filter(., Type == "All novel"),
    key = "Novel_peptide_count",
    value = "count",
    group = "ORFNoveltyReason",
    fill = "ORFNoveltyReason", posit = "dodge",
    main = "All ORF frequency per peptide count",
    xlabel = "Number of novel peptides",
    ylabel = "Count of ORF",
    textsize = 15,
    label = "count",
    legend = "none",
    xdir = "horizontal")
plot(pl_all_orf_freq[[1]] + facet_wrap(facets = "ORFNoveltyReason"))
pl_qc_orf_freq <- plots_hist(
    data = toplot %>%
        dplyr::filter(., Type == "Quality filtered"),
    key = "Novel_peptide_count",
    value = "count",
    group = "ORFNoveltyReason",
    fill = "ORFNoveltyReason", posit = "dodge",
    main = "High quality ORF frequency per peptide count",
    xlabel = "Number of novel peptides",
    ylabel = "Count of ORF",
    textsize = 15,
    label = "count",
    legend = "none",
    xdir = "horizontal")
plot(pl_qc_orf_freq[[1]] + facet_wrap(facets = "ORFNoveltyReason"))

# Close the device
dev.off()

# Export neighbour results (as txt and RDS files)
saveRDS(
    object = neighbours_analysis,
    file = paste(
        opt$output, "/", "Neighbours_analysis.RDS", sep = ""))
write.table(
    x = neighbours_analysis,
    file = paste(
        opt$output, "/", "Neighbours_analysis.txt", sep = ""),
    quote = FALSE,
    sep = "\t",
    row.names = FALSE,
    col.names = TRUE)

# Export operon overlap results (as txt and RDS files)
if (exists("overlap_operon")) {
    saveRDS(
        object = overlap_operon,
        file = paste(
            opt$output, "/", "Overlap_operon.RDS", sep = ""))
    write.table(
        x = overlap_operon,
        file = paste(
            opt$output, "/", "Overlap_operon.txt", sep = ""),
        quote = FALSE,
        sep = "\t",
        row.names = FALSE,
        col.names = TRUE)
}

# Export all ORF novelty reason results (as txt and RDS files)
saveRDS(
    object = orf_reason_final,
    file = paste(
        opt$output, "/", "ORF_novelty_reason.RDS", sep = ""))
write.table(
    x = orf_reason_final,
    file = paste(
        opt$output, "/", "ORF_novelty_reason.txt", sep = ""),
    quote = FALSE,
    sep = "\t",
    row.names = FALSE,
    col.names = TRUE)

# Export all ORF novelty reason results (as txt and RDS files)
saveRDS(
    object = orf_reason_highqual,
    file = paste(
        opt$output, "/", "HighQual_ORF_novelty_reason.RDS", sep = ""))
write.table(
    x = orf_reason_highqual,
    file = paste(
        opt$output, "/", "HighQual_ORF_novelty_reason.txt", sep = ""),
    quote = FALSE,
    sep = "\t",
    row.names = FALSE,
    col.names = TRUE)



### Genomic visualisation and ORF coverage -------------------------------

# Open a file for plot visualisation
pdf(
    file = paste0(opt$output, "/", "ORF_genomic_visualisation.pdf"),
    width = 10, height = 10)

# Define reference entries that are expressed
expr_known <- evid_reason %>%
    dplyr::filter(., group == "Known") %>%
    strsplit(x = .[["Proteins"]], split = ";", fixed = TRUE) %>%
    unlist(.) %>%
    unique(.) %>%
    data.frame(id = ., Expressed = TRUE) %>%
    dplyr::left_join(
        x = data.frame(id = names(ref_grange)), y = ., by = "id") %>%
    dplyr::mutate(
        ., Expressed = ifelse(is.na(Expressed), FALSE, Expressed)) %>%
    set_rownames(.[["id"]]) %>%
    dplyr::select(., -id)

# Include the entries epression pattern into the GRange metadata
ref_grange_expr <- ref_grange
values(ref_grange_expr) <- cbind(
    values(ref_grange_expr), expr_known)

# Compute and visualise the genomic coverage based on nucleotide coverage
pl_rectvenn <- plots_rectvenn(
    ideo = genome_grange, ref = ref_grange_expr, pep = pep_grange)
pl_rectvenn

# Define novel entries that are expressed
expr_novel <- evid_reason %>%
    dplyr::filter(., group == "Novel") %>%
    strsplit(x = .[["Proteins"]], split = ";", fixed = TRUE) %>%
    unlist(.) %>%
    unique(.) %>%
    data.frame(id = ., Expressed = TRUE) %>%
    dplyr::left_join(
        x = data.frame(id = names(orf_grange)), y = ., by = "id") %>%
    dplyr::mutate(
        ., Expressed = ifelse(is.na(Expressed), FALSE, Expressed)) %>%
    set_rownames(.[["id"]]) %>%
    dplyr::select(., -id)

# Include the entries expression pattern into the GRange metadata
orf_grange_expr <- orf_grange
values(orf_grange_expr) <- cbind(
    values(orf_grange_expr), expr_novel)

# Use ggbio extension to plot ORF location on genome as a circos graph
colou <- c("#4682B4", "#BD5E5E", "#437A3C", "#F57C36", "#D58DEB", "#B2B83F")
pl_circos <- ggplot() +
    ggtitle(label = organism(my_bsgeno)) +
    layout_circle(
        genome_grange, geom = "ideo", fill = "gray70",
        radius = 30, trackWidth = 2) +
    layout_circle(
        genome_grange, geom = "scale", size = 4,
        radius = 33, trackWidth = 2) +
    layout_circle(
        genome_grange, geom = "text", size = 7, aes(label = seqnames),
        angle = 0, radius = 37, trackWidth = 5) +
    layout_circle(
        subset(x = ref_grange_expr, strand == "+" & Expressed),
        geom = "rect", color = colou[1],
        radius = 26, trackWidth = 4) +
    layout_circle(
        subset(x = ref_grange_expr, strand == "-" & Expressed),
        geom = "rect", color = colou[2],
        radius = 22, trackWidth = 4) +
    layout_circle(
        subset(x = orf_grange_expr, strand == "+" & Expressed),
        geom = "rect", color = colou[3],
        radius = 18, trackWidth = 4) +
    layout_circle(
        subset(x = orf_grange_expr, strand == "-" & Expressed),
        geom = "rect", color = colou[4],
        radius = 14, trackWidth = 4) +
    annotate(
        geom = "text", x = -12, y = 6, hjust = 0,
        label = "1. Annotated ORF (+ strand)", colour = colou[1]) +
    annotate(
        geom = "text", x = -12, y = 2, hjust = 0,
        label = "2. Annotated ORF (- strand)", colour = colou[2]) +
    annotate(
        geom = "text", x = -12, y = -2, hjust = 0,
        label = "3. Putative novel ORF (+ strand)", colour = colou[3]) +
    annotate(
        geom = "text", x = -12, y = -6, hjust = 0,
        label = "4. Putative novel ORF (- strand)", colour = colou[4])
plot(pl_circos)

# Get all peptide associated nucleotide position
coverage_pep <- gr_nucleotide_pos(
    grange = pep_grange, filter = 'grepl("Known|Target", Database)')

# Convert to dataframe
coverage_nucl <- coverage_pep %>%
    ldply(., "data.frame") %>%
    set_colnames(c("Sequence", "Nucl_pos")) %>%
    dplyr::filter(., !is.na(Nucl_pos))

# Compute msms count for each peptide sequence
coverage_nucl_count <- evid_reason %>%
    dplyr::select(., Sequence, `MS/MS count`) %>%
    dplyr::group_by(., Sequence) %>%
    dplyr::summarise(., Count = sum(`MS/MS count`)) %>%
    dplyr::left_join(coverage_nucl, ., by = "Sequence")

# Format dataframe count column to factor
coverage_nucl_count$Count <- factor(
    x = coverage_nucl_count$Count,
    levels = seq(1, max(coverage_nucl_count$Count, na.rm = TRUE)),
    labels = seq(1, max(coverage_nucl_count$Count, na.rm = TRUE)),
    ordered = TRUE)

# Calculate overall nucleotide coverage frequencies
toplot <- coverage_nucl_count %>%
    plyr::ddply(
        .data = ., .variables = c("Count"),
        .fun = summarise, Freq = length(Sequence),
        .drop = FALSE)

# Calculate the quantiles of count frequencies
quantiles_toplot <- quantile(
    x = as.integer(as.character(coverage_nucl_count$Count)),
    probs = c(0, 0.25, 0.5, 0.75, 1), na.rm = TRUE) %>%
    as.data.frame(.) %>%
    set_colnames("Values") %>%
    dplyr::mutate(., Quartiles = row.names(.))
quantiles_toplot <- data.frame(
    Quartiles = "Mean",
    Values = base::mean(as.integer(as.character(coverage_nucl_count$Count)), na.rm = TRUE)) %>%
    dplyr::bind_rows(quantiles_toplot, .) %>%
    dplyr::select(., Quartiles, Values) %>%
    dplyr::mutate(., Values = round(x = Values, digits = 1))

# Generate the histogram frequency of MS/MS count per nucleotide
pl_coverage <- plots_hist(
    data = toplot %>%
        dplyr::filter(., Count %in% c(1:160)),
    key = "Count",
    value = "Freq",
    group = "Count",
    fill = "grey",
    main = "Coverage per nucleotide",
    xlabel = "MS/MS counts",
    ylabel = "Nucleotide counts",
    textsize = 25)

# Add quantiles table to the graph
pl_coverage <- pl_coverage[[1]] +
    annotation_custom(
        grob = tableGrob(
            d = quantiles_toplot, theme = ttheme_minimal(), rows = NULL),
        xmin = 120, xmax = 150, ymin = 1E5, ymax = 3E5) +
    scale_x_discrete(
        breaks = c(1, seq(20, 160, by = 20)),
        labels = c(1, seq(20, 160, by = 20)))
plot(pl_coverage)
dev.off()

# Define a high quality target list
all_targets <- orf_reason_final %>%
    #dplyr::filter(., Novel_peptide_count >= 2) %>%
    .[["Proteins"]]

# Filter the peptide GRange for duplicate ID
pep_grange_unique <- subset(
    pep_grange,
    !(id %in% unique(names(pep_grange)[duplicated(names(pep_grange))])))

# Loop through all high quality candidate ORFs
#warning("Generating genomic visualisation only for high quality ORFs!")
pl_genome_list <- list()
for (i in all_targets) {
    
    # Generate all genomic representation for reference and novel entries
    # as well as peptide and sanger sequences
    genomic_vis_data <- plots_orf_genomic(
        x = i,
        bsgeno = my_bsgeno,
        ref_gr = ref_grange_expr,
        orf_gr = orf_grange_expr,
        pep_gr = pep_grange_unique,
        sanger_gr = sanger_grange,
        #ref_label = "GENE.NAME")
        ref_label = "protein")
        #ref_label = "PROTEIN.ID")
    
    # Plot the known entries and 6-frames ORFs tracks
    my_plots <- genomic_vis_data[["plots"]]
    pl_genome <- tracks(
        `Known ORFs` = my_plots$`Known ORFs`,
        `Frame 1` = my_plots$`Frame 1`,
        `Frame 2` = my_plots$`Frame 2`,
        `Frame 3` = my_plots$`Frame 3`,
        `Frame -1` = my_plots$`Frame -1`,
        `Frame -2` = my_plots$`Frame -2`,
        `Frame -3` = my_plots$`Frame -3`,
        heights = c(0.5, rep(0.1, times = 6)),
        title = paste0(
            i,
            ": ",
            as.character(orf_reason_final[
                orf_reason_final$Proteins == i, "ORFNoveltyReason"])),
        #xlab = "Genomic position",
        label.bg.fill = "white",
        label.text.angle = 45) +
        theme_bw() +
        theme(
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            axis.text.y =  element_blank(),
            axis.ticks.y = element_blank(),
            legend.position = "none")
    xlim(pl_genome) <- c(
        genomic_vis_data$region_coordinates[["start"]],
        genomic_vis_data$region_coordinates[["end"]])
    #print(pl_genome)
    
    # Plot the peptides and Sanger sequence tracks
    pl_genome_test <- tracks(`Genome` = my_plots$Genome)
    my_test <- try(print(pl_genome_test), silent = TRUE)
    
    if (attr(my_test,"class") == "try-error") {
        warning("Nucleotide sequence not drawn due to incompatible seqname")
        pl_genome_zoom <- tracks(
            `Peptide` = my_plots$Peptide,
            `Sanger` = my_plots$`Sequenced PCR`,
            heights = c(2, 0.5),
            title = paste0(
                i,
                ": zoomed in"),
            xlab = "Genomic position",
            label.bg.fill = "grey",
            label.text.angle = 45,
            track.bg.color = "grey") +
            theme_bw() +
            theme(
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                axis.text.y =  element_blank(),
                axis.title.y =  element_blank(),
                axis.ticks.y = element_blank(),
                legend.position = "right")
    } else {
        pl_genome_zoom <- tracks(
            `Peptide` = my_plots$Peptide,
            `Sanger` = my_plots$`Sequenced PCR`,
            `Genome` = my_plots$Genome,
            heights = c(2, 0.5, 0.2),
            title = paste0(
                i,
                ": zoomed in"),
            xlab = "Genomic position",
            label.bg.fill = "grey",
            label.text.angle = 45,
            track.bg.color = "grey") +
            theme_bw() +
            theme(
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                axis.text.y =  element_blank(),
                axis.title.y =  element_blank(),
                axis.ticks.y = element_blank(),
                legend.position = "right")
    }
    
    xlim(pl_genome_zoom) <- c(
        genomic_vis_data$region_zoom[["start"]],
        genomic_vis_data$region_zoom[["end"]])

    # Use cowplot to combine the two genome tracks
    g <- as(pl_genome, "grob")
    gz <- as(pl_genome_zoom, "grob")
    pl_genome_merge <- cowplot::plot_grid(
        g, gz, align = "h", axis = "lr", nrow = 2)
    #print(pl_genome_merge)
    
    # Include the plots for current candidate into list
    pl_genome_list[i] <- list(pl_genome_merge)
    
}

# Open a file for plot visualisation
pdf(
    file = paste0(opt$output, "/", "ORF_visualisation.pdf"),
    width = 10, height = 10)
pl_genome_list
dev.off()

# Export all plots as RDS file
all_figures <- list(
    datab_count = pl_datab_count[[1]],
    datab_pep = pl_datab_pep[[1]],
    reason = pl_orf_reason[[1]],
    orf_freq = pl_all_orf_freq[[1]],
    orf_qc_freq = pl_qc_orf_freq[[1]],
    venn = pl_rectvenn,
    circos = pl_circos,
    coverage = pl_coverage,
    genome_vis = pl_genome_list)
saveRDS(
    object = all_figures,
    file = paste(
        opt$output, "/", "Novelty_ORF_figures.RDS", sep = ""))



### END ------------------------------------------------------------------

# Close the cluster
stopImplicitCluster()

# Time scripts end
print(paste("Completed ", format(Sys.time(), "%Y-%m-%d"), sep = ""))


