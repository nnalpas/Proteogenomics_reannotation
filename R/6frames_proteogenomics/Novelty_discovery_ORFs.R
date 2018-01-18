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
load_package("motifRG")
load_package("ggseqlogo")
load_package("BSgenome.Bsubtilis.EMBL.AL0091263")



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

# Get the Bsu genome object
bsu <- BSgenome.Bsubtilis.EMBL.AL0091263

# Create new GRange object to study nucleotide frequence of start codon
# on positive strand
plus_grange <- subset(ref_grange, strand == "+")
names(plus_grange) %<>%
    sub("^", "start_", .)
end(plus_grange) <- (start(plus_grange) + 2)

# Create new GRange object to study nucleotide frequence of start codon
# on negative strand
minus_grange <- subset(ref_grange, strand == "-")
names(minus_grange) %<>%
    sub("^", "start_", .)
start(minus_grange) <- (end(minus_grange) - 2)

# Compile positive and negative strand entries
start_codon <- c(plus_grange, minus_grange)

# Get the nucleotide sequence associated with the start codon
start_codon_seq <- getSeq(x = bsu, names = start_codon)

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

# Loop through the stratified list of proteins
rbs_motifs <- list()
rbs_motifs_score <- list()
for (i in 1:length(ref_expr_list)) {
    
    # Create foreground GRange object (+/- X bp of start)
    # to study RBS nucleotide frequence for entries on + and - strand
    pos_grange <- subset(
        ref_grange, strand == "+" & id %in% ref_expr_list[[i]])
    ranges(pos_grange) <- IRanges(
        start = (start(pos_grange) - 30),
        end = (start(pos_grange) + 2))
    neg_grange <- subset(
        ref_grange, strand == "-" & id %in% ref_expr_list[[i]])
    ranges(neg_grange) <- IRanges(
        start = (end(neg_grange) - 2),
        end = (end(neg_grange) + 30))
    
    # Concatenate the entries from positive and negative strand
    fg_grange <- c(pos_grange, neg_grange)
    names(fg_grange) <- sub("^", "fg_", fg_grange$id)
    
    # Get the nucleotide sequence associated with the grange object
    fg_seq <- getSeq(x = bsu, names = fg_grange)
    fg_seq_list <- as.vector(fg_seq)
    
    # Create background GRange object (+/- X bp of end)
    # to study RBS nucleotide frequence for entries on + and - strand
    pos_grange <- subset(
        ref_grange, strand == "+" & id %in% ref_expr_list[[i]])
    ranges(pos_grange) <- IRanges(
        start = (end(pos_grange) - 30),
        end = (end(pos_grange) + 2))
    neg_grange <- subset(
        ref_grange, strand == "-" & id %in% ref_expr_list[[i]])
    ranges(neg_grange) <- IRanges(
        start = (start(neg_grange) - 2),
        end = (start(neg_grange) + 30))
    
    # Concatenate the entries from positive and negative strand
    bg_grange <- c(pos_grange, neg_grange)
    names(bg_grange) <- sub("^", "bg_", bg_grange$id)
    
    # Get the nucleotide sequence associated with the grange object
    bg_seq <- getSeq(x = bsu, names = bg_grange)
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
    ggplot() +
        geom_logo(
            data = fg_seq_list,
            method = "bits",
            seq_type = "dna") +
        theme_logo() +
        ggtitle(paste(names(ref_expr_list)[i], "proteins")) +
        scale_x_continuous(
            name = "Nucleotide position prior CDS",
            breaks = c(seq(from = 1, to = 33, by = 1)),
            labels = c(seq(from = -30, to = 2, by = 1))) +
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


# User provided additional RBS motifs sequence
user_rbs <- c(
    "aaaggaggtgt", "agaggtggtgt", "atattaagaggaggag", "agagaacaaggagggg")


# Concatenate all RBS motifs sequence in single vector
all_rbs_motifs <- c(
    rbs_motifs_score$`Low abundant`$Motifs,
    rbs_motifs_score$`High abundant`$Motifs) %>%
    c(., toupper(user_rbs)) %>%
    unique(.)

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
        
        # Include 30 nucleotides prior to start codon in grange
        if (as.character(strand(tmp_grange)) == "+") {
            ranges(tmp_grange) <- IRanges(
                start = (start(tmp_grange) - 30),
                end = (end(tmp_grange)))
        } else if (as.character(strand(tmp_grange)) == "-") {
            ranges(tmp_grange) <- IRanges(
                start = (start(tmp_grange)),
                end = (end(tmp_grange) + 30))
        } else {
            stop("Strand unknown!")
        }
        names(tmp_grange) <- x
        
        # Get sequence for updated grange
        tmp_seq <- getSeq(x = bsu, names = tmp_grange)
        
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


