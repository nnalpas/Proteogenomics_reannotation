#!/usr/bin/env Rscript

# This script allows plotting of all ORF based on genomic coordinates
# as well as strand and frame of translation



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
#source(
#    file = paste(
#        "C:/Users",
#        user,
#        "Documents/GitHub/Miscellaneous/R/General/General_function.R",
#        sep = "/"))
source(
    file = paste(
        "/home-link",
        user,
        "bin/General_function.R",
        sep = "/"))

# Load the required packages (or install if not already in library)
load_package(plyr)
load_package(dplyr)
load_package(tidyr)
load_package(seqinr)
load_package(UniProt.ws)
load_package(magrittr)
load_package(WriteXLS)
load_package(data.table)
load_package(splitstackshape)
load_package(VennDiagram)
load_package(ggplot2)
load_package(grid)
load_package(gridExtra)
load_package(RColorBrewer)
load_package(stringr)
load_package(Biostrings)
load_package(RecordLinkage)
load_package(VariantAnnotation)
load_package(cgdsr)
load_package(bit64)
load_package(cleaver)
load_package(plotly)
load_package(GenomicRanges)
load_package(biovizBase)
load_package(ggbio)
load_package(ggradar)



### Parameters setting up ------------------------------------------------

# Define the list of command line parameters
option_list <- list(
    make_option(
        opt_str = c("-i", "--idmap"),
        type = "character", default = NULL,
        help = "ID map file", metavar = "character"),
    make_option(
        opt_str = c("-c", "--coordinates"),
        type = "character", default = NULL,
        help = "Coordinates file", metavar = "character"),
    make_option(
        opt_str = c("-o", "--out_path"),
        type = "character", default = NULL, 
        help = "Output directory", metavar = "character"))

# Parse the parameters provided on command line by user
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# Check whether inputs parameter was provided
if (is.null(opt$coordinates)) {
    
    print_help(opt_parser)
    stop("The coordinates argument must be supplied!")
    
}

# Check whether output parameter was provided
if (is.null(opt$out_path)){
    
    opt$output <- dirname(opt$input)
    warning(paste("Output results to path: ", opt$output, "!", sep = ""))
    
}

# For manual parameters set-up
#opt <- list(
#    coordinates = "C:/Users/kxmna01/Dropbox/Home_work_sync/Work/Colleagues shared work/Vaishnavi Ravikumar/Bacillus_subtilis_6frame/Databases/6frames_orf_coordinates.txt",
#    out_path = "C:/Users/kxmna01/Dropbox/Home_work_sync/Work/Colleagues shared work/Vaishnavi Ravikumar/Bacillus_subtilis_6frame/Circos")



### Data import ----------------------------------------------------------

# Import the orf coordinates
orf_coord <- read.table(
    file = opt$coordinates, header = TRUE,
    sep = "\t", quote = "", as.is = TRUE)



### Ideogram generation --------------------------------------------------

# Dataframe holding genome information for bacillus subtilis
tmp <- base::data.frame(
    Chromosome = 1, Strand = "*", Start = 1, End = 4215606,
    name = "chr1", length = 4215606, Type = TRUE, geno = "AL009126.3",
    stringsAsFactors = FALSE)

# Create an Ideogram (GRanges object) for B. subtilis
bsu_ideo <- with(
    data = tmp,
    expr = GRanges(
        seqnames = name,
        ranges = IRanges(Start, End),
        strand = Strand,
        Chromosome = Chromosome,
        seqinfo = Seqinfo(
            seqnames = name,
            seqlengths = length,
            isCircular = Type,
            genome = geno)))

# Export the bsu ideogram for reuse at later stage
saveRDS(
    object = bsu_ideo,
    file = paste0(opt$out_path, "/bsu_AL009126.3_ideo.RDS"))



### ORF visualisation ----------------------------------------------------

# Format dataframe as genomic position and other info for each reference ORF
tmp <- orf_coord_filt %>%
    dplyr::filter(., !is.na(UniProtID)) %>%
    base::as.data.frame(., stringsAsFactors = FALSE)
tmp %<>%
    dplyr::mutate(
        .,
        Start = ifelse(test = strand == 1, yes = start, no = end),
        End = ifelse(test = strand == 1, yes = end, no = start),
        Strand = ifelse(test = strand == 1, yes = "+", no = "-"),
        UniProtKBID = uni_id_clean(UniProtID),
        Chromosome = 1,
        chr.name = "chr1",
        length = 4215606,
        Type = TRUE,
        geno = "AL009126.3") %>%
    dplyr::select(., -start, -end, -strand, -UniProtID) %>%
    dplyr::group_by(., id, Start, End, Strand, frame, Chromosome, chr.name) %>%
    summarise_each(
        funs(toString(x = unique(.), width = NULL))) %>%
    base::as.data.frame(., stringsAsFactors = FALSE)

# Create a GRanges object for all reference and identified novel ORF 
ref_grange <- with(
    data = tmp,
    expr = GRanges(
        seqnames = chr.name,
        ranges = IRanges(Start, End),
        strand = Strand))

# Add values and seqinfo to the created GRanges object
values(ref_grange) <- tmp %>%
    dplyr::select(., id, UniProtKBID, association_id, frame, Chromosome)
seqinfo(ref_grange) <- Seqinfo(
    seqnames = tmp$chr.name %>% unique(.) %>% as.character(.),
    seqlengths = tmp$length %>% unique(.) %>% as.integer(.),
    isCircular = tmp$Type %>% unique(.) %>% as.logical(.),
    genome = tmp$geno %>% unique(.) %>% as.character(.))

# Format dataframe as genomic position and other info for each novel ORF
tmp <- orf_coord_filt %>%
    dplyr::filter(., is.na(UniProtID)) %>%
    base::as.data.frame(., stringsAsFactors = FALSE)
tmp %<>%
    dplyr::mutate(
        .,
        Start = ifelse(test = strand == 1, yes = start, no = end),
        End = ifelse(test = strand == 1, yes = end, no = start),
        Strand = ifelse(test = strand == 1, yes = "+", no = "-"),
        UniProtKBID = uni_id_clean(UniProtID),
        Chromosome = 1,
        chr.name = "chr1",
        length = 4215606,
        Type = TRUE,
        geno = "AL009126.3") %>%
    dplyr::select(., -start, -end, -strand, -UniProtID) %>%
    dplyr::group_by(., id, Start, End, Strand, frame, Chromosome, chr.name) %>%
    summarise_each(
        funs(toString(x = unique(.), width = NULL))) %>%
    base::as.data.frame(., stringsAsFactors = FALSE)

# Create a GRanges object for all reference and identified novel ORF 
orf_grange <- with(
    data = tmp,
    expr = GRanges(
        seqnames = chr.name,
        ranges = IRanges(Start, End),
        strand = Strand))

# Add values and seqinfo to the created GRanges object
values(orf_grange) <- tmp %>%
    dplyr::select(., id, UniProtKBID, association_id, frame, Chromosome)
seqinfo(orf_grange) <- Seqinfo(
    seqnames = tmp$chr.name %>% unique(.) %>% as.character(.),
    seqlengths = tmp$length %>% unique(.) %>% as.integer(.),
    isCircular = tmp$Type %>% unique(.) %>% as.logical(.),
    genome = tmp$geno %>% unique(.) %>% as.character(.))

# Use ggbio extension to plot ORF location on genome as a circos graph
colou <- c("#4682B4", "#BD5E5E", "#437A3C", "#F57C36", "#D58DEB", "#B2B83F")
pl <- ggplot() +
    ggtitle(label = "Bacillus subtilis ORFs") +
    layout_circle(
        bsu_ideo, geom = "ideo", fill = "gray70",
        radius = 30, trackWidth = 3) +
    layout_circle(
        bsu_ideo, geom = "scale", size = 2, radius = 33, trackWidth = 2) +
    layout_circle(
        bsu_ideo, geom = "text", aes(label = seqnames),
        angle = 0, radius = 36, trackWidth = 5) +
    layout_circle(
        subset(x = ref_grange, strand == "+"), geom = "rect", color = colou[1],
        radius = 26, trackWidth = 3) +
    layout_circle(
        subset(x = ref_grange, strand == "-"), geom = "rect", color = colou[2],
        radius = 23, trackWidth = 3) +
    layout_circle(
        subset(x = orf_grange, strand == "+"), geom = "rect", color = colou[3],
        radius = 20, trackWidth = 3) +
    layout_circle(
        subset(x = orf_grange, strand == "-"), geom = "rect", color = colou[4],
        radius = 17, trackWidth = 3)
pl



