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
library(gtable)
library(grid)
library(gridExtra)
library(purrr)
library(foreach)
library(doParallel)
library(motifRG)
library(ggseqlogo)
library(BSgenome.Srimosus.Strathclyde.062019)



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
        pep_pos = choose.files(
            caption = "Choose peptide position within protein (.RDS) file!",
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
        add_rbs = readline(prompt = paste0(
            "Provide additional RBS sequence",
            " (separated by space)?")) %>%
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
            opt_str = c("-e", "--pep_pos"),
            type = "character", default = NULL,
            help = "Peptide position file",
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
            opt_str = c("-a", "--add_rbs"),
            type = "character", default = NULL,
            help = "Additional RBS sequence (separated by space)",
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
    identical(opt$pep_pos, NULL) |
    identical(opt$pep_pos, "") |
    identical(opt$pep_pos, character(0))) {
    
    print_help(opt_parser)
    stop("The input peptide position within proteins must be supplied!")
    
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
    
    print_help(opt_parser)
    stop("The input Sanger seq entries grange must be supplied!")
    
}
if (
    identical(opt$genome_grange, NULL) |
    identical(opt$genome_grange, "") |
    identical(opt$genome_grange, character(0))) {
    
    print_help(opt_parser)
    stop("The input genome entry grange must be supplied!")
    
}
if (
    identical(opt$add_rbs, NULL) |
    identical(opt$add_rbs, "") |
    identical(opt$add_rbs, character(0))) {
    
    opt["add_rbs"] <- list(NULL)
    warning("No RBS sequence provided by user!")
    
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

# Import the peptide position file
pep_loc_data <- opt$pep_pos %>%
    as.character(.) %>%
    readRDS(file = .)
 
# Import the peptides genomic ranges file
pep_grange <- opt$pep_grange %>%
    as.character(.) %>%
    readRDS(file = .)

# Import the Sanger seq genomic ranges file
sanger_grange <- opt$sanger_grange %>%
    as.character(.) %>%
    readRDS(file = .)

# Import the genome genomic ranges file
genome_grange <- opt$genome_grange %>%
    as.character(.) %>%
    readRDS(file = .)

# Import the operon genomic ranges file if defined
if (!is.null(opt$operon)) {
    operon_grange <- opt$operon %>%
        as.character(.) %>%
        readRDS(file = .)
}

# Format the RBS sequences if defined
if (!is.null(opt$add_rbs)) {
    user_rbs <- opt$add_rbs %>%
        as.character(.) %>%
        strsplit(x = ., split = " ", fixed = TRUE) %>%
        unlist(.) %>%
        toupper(.) %>%
        unique(.)
} else {
    user_rbs <- character(0)
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



### RBS motif analysis ---------------------------------------------------

# Get the BSgenome object
bsgeno <- BSgenome.Srimosus.Strathclyde.062019

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
start_codon_seq <- getSeq(x = bsgeno, names = start_codon)

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

# Compute the MS/MS count for all known protein and ORF within range of
# the target and for each raw file
msms_count <- evid_reason %>%
    dplyr::filter(., group == "Known") %>%
    dplyr::group_by(., `Leading Proteins`, `Raw file`) %>%
    dplyr::summarise(., Count = sum(`MS/MS Count`)) %>%
    set_colnames(c("Proteins", "Raw", "Count")) %>%
    dplyr::ungroup(.) %>%
    dplyr::group_by(., Proteins) %>%
    dplyr::summarise(., Count = median(Count)) %>%
    dplyr::arrange(., Count) %>%
    cSplit(
        indt = ., splitCols = "Proteins", sep = ";",
        direction = "long", fixed = TRUE) %>%
    dplyr::filter(., !grepl("^seq_", Proteins))

# Get the highest and lowest expressed protein
low_express <- msms_count %>%
    dplyr::filter(., Count <= median(Count)) %>%
    .[["Proteins"]] %>%
    as.character(.)
high_express <- msms_count %>%
    dplyr::filter(., Count > median(Count)) %>%
    .[["Proteins"]] %>%
    as.character(.)

# Create list of protein stratified by abundance
ref_expr_list <- list(
    `Low abundant` = low_express, `High abundant` = high_express)

# Define number of nucleotide before first nucleotide (to define RBS region)
nucleotide_before <- 30

# Loop through the stratified list of proteins
rbs_motifs <- list()
rbs_motifs_score <- list()
rbs_plot <- list()
for (i in 1:length(ref_expr_list)) {
    
    # Create foreground GRange object (+/- X bp of start)
    # to study RBS nucleotide frequence for entries on + and - strand
    pos_grange <- subset(
        ref_grange, strand == "+" & id %in% ref_expr_list[[i]])
    ranges(pos_grange) <- IRanges(
        start = (start(pos_grange) - nucleotide_before),
        end = (start(pos_grange) + nucleotide_after))
    neg_grange <- subset(
        ref_grange, strand == "-" & id %in% ref_expr_list[[i]])
    ranges(neg_grange) <- IRanges(
        start = (end(neg_grange) - nucleotide_after),
        end = (end(neg_grange) + nucleotide_before))
    
    # Concatenate the entries from positive and negative strand
    fg_grange <- c(pos_grange, neg_grange)
    names(fg_grange) <- sub("^", "fg_", fg_grange$id)
    
    # Get the nucleotide sequence associated with the grange object
    fg_seq <- getSeq(x = bsgeno, names = fg_grange)
    fg_seq_list <- as.vector(fg_seq)
    
    # Create background GRange object (+/- X bp of end)
    # to study RBS nucleotide frequence for entries on + and - strand
    pos_grange <- subset(
        ref_grange, strand == "+" & id %in% ref_expr_list[[i]])
    ranges(pos_grange) <- IRanges(
        start = (end(pos_grange) - nucleotide_before),
        end = (end(pos_grange) + nucleotide_after))
    neg_grange <- subset(
        ref_grange, strand == "-" & id %in% ref_expr_list[[i]])
    ranges(neg_grange) <- IRanges(
        start = (start(neg_grange) - nucleotide_after),
        end = (start(neg_grange) + nucleotide_before))
    
    # Concatenate the entries from positive and negative strand
    bg_grange <- c(pos_grange, neg_grange)
    names(bg_grange) <- sub("^", "bg_", bg_grange$id)
    
    # Get the nucleotide sequence associated with the grange object
    bg_seq <- getSeq(x = bsgeno, names = bg_grange)
    bg_seq_list <- as.vector(bg_seq)
    
    # Loop through nucleotide length (for RBS motif)
    tmp_rbs_motifs <- list()
    tmp_rbs_motifs_score <- c()
    for (x in c(6:12)) {
        
        # Find the motif between foreground and background
        tmp <- findMotifFgBg(
            fg.seq = fg_seq, bg.seq = bg_seq,
            start.width = x, both.strand = FALSE, flank = 1,
            max.width = 16, enriched.only = TRUE)
        
        # Compile all motifs data if existing
        if (length(tmp$motifs) != 0) {
            tmp_rbs_motifs[paste0("length_", x)] <- list(tmp)
            for (y in names(tmp$motifs)) {
                tmp_rbs_motifs_score <- c(
                    tmp_rbs_motifs_score, tmp$motifs[[y]]@score)
            }
        }
        
    }
    
    # Create dataframe from the motif list
    tmp_rbs_motifs_score %<>%
        ldply(., "data.frame") %>%
        set_colnames(c("Motifs", "Score"))
    
    # Draw consensus logo sequence of RBS for expressed protein
    rbs_plot[[names(ref_expr_list)[i]]] <- ggplot() +
        geom_logo(
            data = fg_seq_list,
            method = "bits",
            seq_type = "dna") +
        theme_logo() +
        ggtitle(paste(names(ref_expr_list)[i], "proteins")) +
        scale_x_continuous(
            name = "Nucleotide position prior CDS",
            breaks = c(seq(
                from = 1,
                to = (1 + nucleotide_after + nucleotide_before),
                by = 1)),
            labels = c(seq(
                from = -nucleotide_before,
                to = nucleotide_after,
                by = 1))) +
        annotation_custom(
            grob = tableGrob(
                d = tmp_rbs_motifs_score %>%
                    dplyr::arrange(., desc(Score)) %>%
                    dplyr::slice(., 1:10),
                theme = ttheme_minimal(base_size = 8), rows = NULL),
            xmin = 1, xmax = 10, ymin = 0.5, ymax = 2)
    
    # Compile all RBS data into list
    rbs_motifs[names(ref_expr_list)[i]] <- list(tmp_rbs_motifs)
    rbs_motifs_score[names(ref_expr_list)[i]] <- list(tmp_rbs_motifs_score)
    
}

# Concatenate all RBS motifs sequence in single vector
all_rbs_motifs <- c(
    rbs_motifs_score$`Low abundant`$Motifs,
    rbs_motifs_score$`High abundant`$Motifs,
    user_rbs) %>%
    unique(.)

# Order RBS motifs by decreasing length
all_rbs_motifs <- all_rbs_motifs[
    order(nchar(all_rbs_motifs), all_rbs_motifs, decreasing = T)]

# Create the RBS motifs regular expression
pattern_rbs_motifs <- paste(all_rbs_motifs, collapse = "|") %>%
    paste0("(", ., ")") %>%
    paste0(., ".{1,14}(", paste(used_start, collapse = "|"), ")")

# Loop through all target ORFs
rbs_results <- list() 
for (x in names(sub_orf_grange)) {
    
    # Get only current ORF genomic position
    tmp_grange <- subset(sub_orf_grange, id == x)
    
    # Check if multiple entry with same name
    if (length(tmp_grange) == 1) {
        
        # Include X nucleotides prior to start codon in grange
        if (as.character(strand(tmp_grange)) == "+") {
            ranges(tmp_grange) <- IRanges(
                start = (start(tmp_grange) - nucleotide_before),
                end = (end(tmp_grange)))
        } else if (as.character(strand(tmp_grange)) == "-") {
            ranges(tmp_grange) <- IRanges(
                start = (start(tmp_grange)),
                end = (end(tmp_grange) + nucleotide_before))
        } else {
            stop("Strand unknown!")
        }
        names(tmp_grange) <- x
        
        # Get sequence for updated grange
        tmp_seq <- getSeq(x = bsgeno, names = tmp_grange)
        
        # Attempt to locate pattern within the sequence
        tmp <- tmp_seq %>%
            str_locate_all(
                string = .,
                pattern = pattern_rbs_motifs)
        
        # Check whether locations were obtained
        if (
            length(tmp[[1]][, "start"]) == 0 |
            length(tmp[[1]][, "end"]) == 0) {
            
            # Add explicitely NA if there is no match
            tmp <- base::data.frame(
                id = names(tmp_seq),
                start = NA_integer_,
                end = NA_integer_,
                seq = NA_character_,
                stringsAsFactors = FALSE)
            
        } else {
            
            # Format as dataframe the current peptide positions
            tmp %<>%
                set_names(x) %>%
                ldply(., "data.frame", .id = "id") %>%
                dplyr::rowwise(.) %>%
                dplyr::mutate(
                    .,
                    seq = substring(
                        text = tmp_seq,
                        first = start,
                        last = end))
            
        }
        
        # Compile all RBS search result in list
        rbs_results[x] <- list(tmp)
        
    } else {
        
        warning("Multiple entry with same name, these will be ignored!")
        
    }
    
}

# Transform list into dataframe
rbs_results %<>%
    plyr::ldply(., "data.frame", .id = "id")

# For each novel ORF get the minimum peptide start position and
# maximum peptide end position
pep_loc_data_cat <- pep_loc_data[["Novel"]] %>%
    dplyr::group_by(., prot) %>%
    dplyr::summarise(
        .,
        minStartPep = min(start, na.rm = T),
        maxEndPep = max(end, na.rm = T))

# For each novel ORF get the minimum start and minimum end RBS position
rbs_results_cat <- rbs_results %>%
    dplyr::filter(., !is.na(start)) %>%
    dplyr::mutate(
        .,
        start = (start - nucleotide_before - 1),
        end = (end - nucleotide_before - 3)) %>%
    dplyr::group_by(., id) %>%
    dplyr::summarise(
        .,
        minStartRBS = min(start, na.rm = T),
        minEndRBS = min(end, na.rm = T))

# Combine the peptide and RBS data for each novel ORF and
# determine whether putative RBS fits the identified peptides position
rbs_explain <- dplyr::left_join(
    x = rbs_results_cat, y = pep_loc_data_cat, by = c("id" = "prot")) %>%
    dplyr::mutate(
        ., RBS_fits_PeptideId = ifelse(minEndRBS <= minStartPep, TRUE, FALSE))



### Novelty reason determination -----------------------------------------

# Get all ORF to explain
orf_reason <- evid_reason %>%
    dplyr::select(
        ., Sequence, Proteins, PEP, `MS/MS Count`, Score, Intensity,
        group, Database, OnlyIdBySite, NoveltyReason, PEPfilter) %>%
    cSplit(
        indt = ., splitCols = "Proteins", sep = ";",
        direction = "long", fixed = TRUE) %>%
    dplyr::filter(., Proteins %in% targets)

# Compile info for each novel ORFs based on all sequences (target and novel)
orf_reason_cat <- orf_reason %>%
    dplyr::group_by(., Proteins) %>%
    dplyr::mutate(
        .,
        all_peptide_count = n_distinct(Sequence),
        all_Sequence = paste(unique(Sequence), collapse = ";"),
        all_min_PEP = min(PEP, na.rm = TRUE),
        all_sum_MSMS_Count = sum(`MS/MS Count`, na.rm = TRUE),
        all_max_Score = max(Score, na.rm = TRUE),
        all_sum_Intensity = sum(Intensity, na.rm = TRUE))

# Compile info for each novel ORFs based on novel only sequences
orf_reason_cat %<>%
    dplyr::filter(., Database == "Novel") %>%
    dplyr::group_by(., Proteins) %>%
    dplyr::summarise(
        .,
        novel_peptide_count = n_distinct(Sequence),
        novel_Sequence = paste(unique(Sequence), collapse = ";"),
        novel_min_PEP = min(PEP, na.rm = TRUE),
        novel_sum_MSMS_Count = sum(`MS/MS Count`, na.rm = TRUE),
        novel_max_Score = max(Score, na.rm = TRUE),
        novel_sum_Intensity = sum(Intensity, na.rm = TRUE),
        all_peptide_count = unique(all_peptide_count),
        all_Sequence = unique(all_Sequence),
        all_min_PEP = unique(all_min_PEP),
        all_sum_MSMS_Count = unique(all_sum_MSMS_Count),
        all_max_Score = unique(all_max_Score),
        all_sum_Intensity = unique(all_sum_Intensity),
        group = paste(unique(group), collapse = ";"),
        Database = paste(unique(Database), collapse = ";"),
        PepNoveltyReason = paste(
            sort(as.character(unique(NoveltyReason))), collapse = ";"),
        OnlyIdBySite = ifelse(all(OnlyIdBySite == FALSE), FALSE, TRUE),
        PEPfilter = ifelse(all(PEPfilter == FALSE), FALSE, TRUE)) %>%
    dplyr::arrange(., desc(novel_peptide_count))

# Clean-up the ORF novelty reason by removing 'Not novel' if another reason
orf_reason_clean <- orf_reason_cat %>%
    dplyr::mutate(., ORFNoveltyReason = sub(
        "(;Not novel|Not novel;)", "", PepNoveltyReason))

# Clean-up the ORF novelty reason by removing 'Potential alternate start'
# if 'Alternate start (known other species)' reason
orf_reason_clean %<>%
    dplyr::mutate(., ORFNoveltyReason = sub(
        "(Alternate start \\(known other species\\));Potential alternate start",
        "\\1",
        ORFNoveltyReason))

# Clean-up the ORF novelty reason by removing 'SAV'
# if 'SAVs/InDels' reason
orf_reason_clean %<>%
    dplyr::mutate(., ORFNoveltyReason = sub(
        "SAV;(SAVs/InDels)",
        "\\1",
        ORFNoveltyReason))

# Include the neighbour reference entry analysis
orf_reason_neighb <- neighbours_analysis %>%
    dplyr::filter(., interpr %in% c("Five neighbour", "Three neighbour")) %>%
    dplyr::select(., queryID, subjectID, interpr) %>%
    tidyr::spread(
        data = ., key = interpr, value = subjectID, convert = FALSE) %>%
    set_colnames(sub("neighbour", "neighbour (+/-3bp)", colnames(.))) %>%
    dplyr::left_join(
        x = orf_reason_clean, y = ., by = c("Proteins" = "queryID"))

# Include the reference protein ID which match novel ORF if any
orf_reason_neighb <- as.data.frame(orf_grange, stringsAsFactors = FALSE) %>%
    dplyr::filter(., !is.na(UNIPROTKB)) %>%
    dplyr::select(., id, UNIPROTKB) %>%
    dplyr::left_join(
        x = orf_reason_neighb, y = ., by = c("Proteins" = "id"))

# Include the operon overlap analysis
orf_reason_opr <- overlap_operon %>%
    dplyr::select(., queryID, subjectID) %>%
    set_colnames(c("id", "OperonID")) %>%
    dplyr::group_by(., id) %>%
    dplyr::summarise(., OperonID = paste(unique(OperonID), collapse = ";")) %>%
    dplyr::left_join(
        x = orf_reason_neighb, y = ., by = c("Proteins" = "id"))

# Include the RBS motif presence analysis
orf_reason_rbs <- rbs_results %>%
    dplyr::filter(., !is.na(start)) %>%
    dplyr::mutate(
        .,
        start = (start - nucleotide_before - 1),
        end = (end - nucleotide_before - 1)) %>%
    tidyr::unite(data = ., col = rbs_position, start, end, sep = "/") %>%
    dplyr::group_by(., id) %>%
    dplyr::summarise(
        ., rbs_position = paste(unique(rbs_position), collapse = ";")) %>%
    dplyr::select(., id, rbs_position) %>%
    dplyr::left_join(
        x = orf_reason_opr, y = ., by = c("Proteins" = "id")) %>%
    dplyr::left_join(
        x = ., y = rbs_explain %>% dplyr::select(., id, RBS_fits_PeptideId),
        by = c("Proteins" = "id"))

# Compute final ORF novelty reason based on five prime neighbour and RBS motif
orf_reason_final <- orf_reason_rbs %>% 
    dplyr::mutate(., ORFNoveltyReason = dplyr::case_when(
        ORFNoveltyReason == "Potential alternate start" &
            !is.na(`Five neighbour (+/-3bp)`) ~ paste(
                ORFNoveltyReason, "Erroneous termination", sep = ";"),
        ORFNoveltyReason == "Potentially novel" &
            !is.na(`Five neighbour (+/-3bp)`) ~ "Erroneous termination",
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



### Focus on high quality novel ORF --------------------------------------

# Keep high quality novel ORF by filtering out the entries only identified
# by sites and entries above PEP threshold
orf_reason_highqual <- orf_reason_final %>%
    dplyr::filter(., OnlyIdBySite & PEPfilter)



### Data visualisation and export ----------------------------------------

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
    dplyr::mutate(., within_operon = !is.na(OperonID)) %>%
    dplyr::group_by(., within_operon) %>%
    dplyr::summarise(., count = n_distinct(Proteins)) %>%
    dplyr::mutate(., Type = "All novel")
toplot <- orf_reason_final %>%
    dplyr::filter(., Proteins %in% unique(orf_reason_highqual$Proteins)) %>%
    dplyr::mutate(., within_operon = !is.na(OperonID)) %>%
    dplyr::group_by(., within_operon) %>%
    dplyr::summarise(., count = n_distinct(Proteins)) %>%
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

# Display the RBS seqlogo results
marrangeGrob(grobs = rbs_plot, ncol = 1, nrow = 1, top = NULL)

# Plot the RBS motif results
toplot <- orf_reason_final %>%
    dplyr::mutate(., rbs_type = dplyr::case_when(
        is.na(rbs_position) ~ "No RBS",
        RBS_fits_PeptideId ~ "Good fit RBS",
        TRUE ~ "Putative RBS")) %>%
    dplyr::group_by(., rbs_type) %>%
    dplyr::summarise(., count = n_distinct(Proteins)) %>%
    dplyr::mutate(., Type = "All novel")
toplot <- orf_reason_final %>%
    dplyr::filter(., Proteins %in% unique(orf_reason_highqual$Proteins)) %>%
    dplyr::mutate(., rbs_type = dplyr::case_when(
        is.na(rbs_position) ~ "No RBS",
        RBS_fits_PeptideId ~ "Good fit RBS",
        TRUE ~ "Putative RBS")) %>%
    dplyr::group_by(., rbs_type) %>%
    dplyr::summarise(., count = n_distinct(Proteins)) %>%
    dplyr::mutate(., Type = "Quality filtered") %>%
    dplyr::bind_rows(toplot, .)
plots_hist(
    data = toplot,
    key = "rbs_type",
    value = "count",
    group = "Type",
    fill = "Type",
    main = "ORF with RBS motifs",
    xlabel = "RBS motif presence",
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
    dplyr::bind_rows(toplot, .)
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
    legend = "right",
    xdir = "vertical")
plot(pl_orf_reason[[1]])

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
    dplyr::bind_rows(toplot, .)
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
    legend = "right",
    xdir = "vertical")
plot(pl_all_orf_reason_operon[[1]])
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
    legend = "right",
    xdir = "vertical")
plot(pl_hq_orf_reason_operon[[1]])

# Visualise ORF per novelty type and RBS motif results
toplot <- orf_reason_final %>%
    dplyr::mutate(., rbs_type = dplyr::case_when(
        is.na(rbs_position) ~ "No RBS",
        RBS_fits_PeptideId ~ "Good fit RBS",
        TRUE ~ "Putative RBS")) %>%
    plyr::ddply(
        .data = ., .variables = .(ORFNoveltyReason, rbs_type),
        .fun = summarise, count = n_distinct(Proteins), .drop = FALSE) %>%
    dplyr::mutate(., Type = "All novel")
toplot <- orf_reason_final %>%
    dplyr::filter(., Proteins %in% unique(orf_reason_highqual$Proteins)) %>%
    dplyr::mutate(., rbs_type = dplyr::case_when(
        is.na(rbs_position) ~ "No RBS",
        RBS_fits_PeptideId ~ "Good fit RBS",
        TRUE ~ "Putative RBS")) %>%
    plyr::ddply(
        .data = ., .variables = .(ORFNoveltyReason, rbs_type),
        .fun = summarise, count = n_distinct(Proteins), .drop = FALSE) %>%
    dplyr::mutate(., Type = "Quality filtered") %>%
    dplyr::bind_rows(toplot, .)
pl_all_orf_reason_rbs <- plots_hist(
    data = toplot %>% dplyr::filter(., Type == "All novel"),
    key = "ORFNoveltyReason",
    value = "count",
    group = "ORFNoveltyReason",
    fill = "rbs_type",
    posit = "stack",
    main = "All ORF novelty with RBS",
    xlabel = "Novelty reason type",
    ylabel = "Count of ORF",
    textsize = 15,
    bw = TRUE,
    legend = "right",
    xdir = "vertical")
plot(pl_all_orf_reason_rbs[[1]])
pl_hq_orf_reason_rbs <- plots_hist(
    data = toplot %>% dplyr::filter(., Type == "Quality filtered"),
    key = "ORFNoveltyReason",
    value = "count",
    group = "ORFNoveltyReason",
    fill = "rbs_type",
    posit = "stack",
    main = "High quality ORF novelty with operon",
    xlabel = "Novelty reason type",
    ylabel = "Count of ORF",
    textsize = 15,
    bw = TRUE,
    legend = "right",
    xdir = "vertical")
plot(pl_hq_orf_reason_rbs[[1]])

# Frequency of ORF per number of novel peptide per novelty type
toplot <- orf_reason_final %>%
    dplyr::select(., Proteins, ORFNoveltyReason, novel_peptide_count) %>%
    plyr::ddply(
        .data = .,
        .variables = .(ORFNoveltyReason, novel_peptide_count),
        .fun = summarise,
        count = dplyr::n_distinct(Proteins),
        .drop = FALSE) %>%
    dplyr::mutate(., Type = "All novel")
toplot <- orf_reason_final %>%
    dplyr::filter(., Proteins %in% unique(orf_reason_highqual$Proteins)) %>%
    dplyr::select(., Proteins, ORFNoveltyReason, novel_peptide_count) %>%
    plyr::ddply(
        .data = .,
        .variables = .(ORFNoveltyReason, novel_peptide_count),
        .fun = summarise,
        count = dplyr::n_distinct(Proteins),
        .drop = FALSE) %>%
    dplyr::mutate(., Type = "Quality filtered") %>%
    dplyr::bind_rows(toplot, .)
pl_all_orf_freq <- plots_hist(
    data = toplot %>%
        dplyr::filter(., Type == "All novel"),
    key = "novel_peptide_count",
    value = "count",
    group = "ORFNoveltyReason",
    fill = "ORFNoveltyReason", posit = "dodge",
    main = "All ORF frequency per peptide count",
    xlabel = "Number of novel peptides",
    ylabel = "Count of ORF",
    textsize = 15,
    label = "count",
    legend = "bottom",
    xdir = "horizontal")
plot(pl_all_orf_freq[[1]])
pl_qc_orf_freq <- plots_hist(
    data = toplot %>%
        dplyr::filter(., Type == "Quality filtered"),
    key = "novel_peptide_count",
    value = "count",
    group = "ORFNoveltyReason",
    fill = "ORFNoveltyReason", posit = "dodge",
    main = "High quality ORF frequency per peptide count",
    xlabel = "Number of novel peptides",
    ylabel = "Count of ORF",
    textsize = 15,
    label = "count",
    legend = "bottom",
    xdir = "horizontal")
plot(pl_qc_orf_freq[[1]])

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

# Export RBS motifs results (as txt and RDS files)
saveRDS(
    object = rbs_results,
    file = paste(
        opt$output, "/", "RBS_motifs_results.RDS", sep = ""))
write.table(
    x = rbs_results,
    file = paste(
        opt$output, "/", "RBS_motifs_results.txt", sep = ""),
    quote = FALSE,
    sep = "\t",
    row.names = FALSE,
    col.names = TRUE)

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



### Genomic visualisation and plots export -------------------------------

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

# Close the device
dev.off()

# Open a file for plot visualisation
pdf(
    file = paste0(opt$output, "/", "ORF_genomic_visualisation.pdf"),
    width = 10, height = 10)

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

# Include the entries epression pattern into the GRange metadata
orf_grange_expr <- orf_grange
values(orf_grange_expr) <- cbind(
    values(orf_grange_expr), expr_novel)

# Use ggbio extension to plot ORF location on genome as a circos graph
colou <- c("#4682B4", "#BD5E5E", "#437A3C", "#F57C36", "#D58DEB", "#B2B83F")
pl_circos <- ggplot() +
    ggtitle(label = organism(bsgeno)) +
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
    grange = pep_grange, filter = 'grepl("Known", Database)')

# Convert to dataframe
coverage_nucl <- coverage_pep %>%
    ldply(., "data.frame") %>%
    set_colnames(c("Sequence", "Nucl_pos")) %>%
    dplyr::filter(., !is.na(Nucl_pos))

# Compute msms count for each peptide sequence
coverage_nucl_count <- evid_reason %>%
    dplyr::select(., Sequence, `MS/MS Count`) %>%
    dplyr::group_by(., Sequence) %>%
    dplyr::summarise(., Count = sum(`MS/MS Count`)) %>%
    dplyr::left_join(coverage_nucl, ., by = "Sequence")

# Format dataframe count column to factor
coverage_nucl_count$Count <- factor(
    x = coverage_nucl_count$Count,
    levels = seq(1, max(coverage_nucl_count$Count)),
    labels = seq(1, max(coverage_nucl_count$Count)),
    ordered = TRUE)

# Calculate overall nucleotide coverage frequencies
toplot <- coverage_nucl_count %>%
    plyr::ddply(
        .data = ., .variables = c("Count"),
        .fun = summarise, Freq = n(),
        .drop = FALSE)

# Calculate the quantiles of count frequencies
quantiles_toplot <- quantile(
    x = as.integer(as.character(coverage_nucl_count$Count)),
    probs = c(0, 0.25, 0.5, 0.75, 1)) %>%
    as.data.frame(.) %>%
    set_colnames("Values") %>%
    dplyr::mutate(., Quartiles = row.names(.))
quantiles_toplot <- data.frame(
    Quartiles = "Mean",
    Values = mean(as.integer(as.character(coverage_nucl_count$Count)))) %>%
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
        xmin = 130, xmax = 150, ymin = 2E5, ymax = 4E5) +
    scale_x_discrete(
        breaks = c(1, seq(20, 160, by = 20)),
        labels = c(1, seq(20, 160, by = 20)))
plot(pl_coverage)

# Define a high quality target list
high_qual_targets <- orf_reason_highqual %>%
    dplyr::filter(., novel_peptide_count >= 2) %>%
    .[["Proteins"]]

# Filter the peptide GRange for duplicate ID
pep_grange_unique <- subset(
    pep_grange,
    !(id %in% unique(names(pep_grange)[duplicated(names(pep_grange))])))

# Loop through all high quality candidate ORFs
warning("Generating genomic visualisation only for high quality ORFs!")
pl_genome_list <- list()
for (i in high_qual_targets) {
    
    # Generate all genomic representation for reference and novel entries
    # as well as peptide and sanger sequences
    genomic_vis_data <- plots_orf_genomic(
        x = i,
        ref_gr = ref_grange_expr,
        orf_gr = orf_grange_expr,
        pep_gr = pep_grange_unique,
        sanger_gr = sanger_grange,
        ref_label = "GENE.NAME")
    
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
            as.character(orf_reason_highqual[
                orf_reason_highqual$Proteins == i, "ORFNoveltyReason"])),
        xlab = "Genomic position",
        label.bg.fill = "white") +
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
    print(pl_genome)
    
    # Plot the peptides and Sanger sequence tracks
    pl_genome_zoom <- tracks(
        `Peptide` = my_plots$Peptide,
        `Sequenced PCR` = my_plots$`Sequenced PCR`,
        `Genome` = my_plots$Genome,
        heights = c(2, 0.5, 0.2),
        title = paste0(
            i,
            ": ",
            as.character(orf_reason_highqual[
                orf_reason_highqual$Proteins == i, "ORFNoveltyReason"])),
        xlab = "Genomic position",
        label.bg.fill = "white") +
        theme_bw() +
        theme(
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            axis.text.y =  element_blank(),
            axis.title.y =  element_blank(),
            axis.ticks.y = element_blank(),
            legend.position = "right")
    xlim(pl_genome_zoom) <- c(
        genomic_vis_data$region_zoom[["start"]],
        genomic_vis_data$region_zoom[["end"]])
    print(pl_genome_zoom)
    
    # Include the plots for current candidate into list
    pl_genome_list[i] <- list(list(
        "pl_genome" = pl_genome, "pl_genome_zoom" = pl_genome_zoom))
    
}

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

# Close the pdf file
dev.off()

# Time scripts end
print(paste("Completed ", format(Sys.time(), "%Y-%m-%d"), sep = ""))


