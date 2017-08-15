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



# temporary



# Import the blast results of Nicolas' ORF versus the reference proteome
blast_NN_vs_ref <- read.table(
    file = "C:/Users/kxmna01/Dropbox/Home_work_sync/Work/Colleagues shared work/Vaishnavi Ravikumar/Bacillus_subtilis_6frame/blastp_ORF_Nicolas_vs_RefProt_19012017",
    header = FALSE,
    sep = "\t",
    quote = "",
    col.names = c(
        "qseqid", "sseqid", "pident", "nident", "mismatch", "length",
        "gapopen", "qstart", "qend", "sstart", "send", "evalue",
        "bitscore", "score"),
    as.is = TRUE)

# Find the best blast hit for each query id
blast_NN_vs_ref_best <- best_blast(
    data = blast_NN_vs_ref, key = "qseqid", multi_match = "remove")

# Keep the mapping of ORF to uniprotID
orf_to_uniprot <- blast_NN_vs_ref_best %>%
    dplyr::select(., qseqid, sseqid) %>%
    set_colnames(c("qseqid", "UniProtID")) %>%
    base::as.data.frame(., stringsAsFactors = FALSE)
orf_to_uniprot$qseqid %<>%
    gsub("^lcl\\|", "", .)
orf_to_uniprot %<>%
    set_colnames(c("id", "UniProtID"))


orf_coord %<>%
    dplyr::left_join(., orf_to_uniprot, by = c("association_id" = "id"))


evid_match <- readRDS("C:/Users/kxmna01/Dropbox/Home_work_sync/Work/Colleagues shared work/Vaishnavi Ravikumar/Bacillus_subtilis_6frame/15082017/2017-08-15_Sequence_group_mapping.RDS")
evid_expr <- evid_match %>%
    strsplit(x = .[["Proteins"]], split = ";", fixed = TRUE) %>%
    unlist(.) %>%
    unique(.)



#orf_coord_filt <- orf_coord %>%
#    dplyr::filter(., id %in% evid_expr | UniProtID %in% evid_expr)






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
tmp <- orf_coord %>%
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
tmp %<>%
    dplyr::filter(., !duplicated(UniProtKBID))

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
tmp <- orf_coord %>%
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
tmp %<>%
    dplyr::filter(., !duplicated(id))

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
ref_grange_expr <- subset(x = ref_grange, UniProtKBID %in% evid_expr)
orf_grange_expr <- subset(x = orf_grange, id %in% evid_expr)
colou <- c("#4682B4", "#BD5E5E", "#437A3C", "#F57C36", "#D58DEB", "#B2B83F")
pl <- ggplot() +
    ggtitle(label = "Bacillus subtilis ORFs") +
    layout_circle(
        bsu_ideo, geom = "ideo", fill = "gray70",
        radius = 30, trackWidth = 2) +
    layout_circle(
        bsu_ideo, geom = "scale", size = 4,
        radius = 33, trackWidth = 2) +
    layout_circle(
        bsu_ideo, geom = "text", size = 7, aes(label = seqnames),
        angle = 0, radius = 37, trackWidth = 5) +
    layout_circle(
        subset(x = ref_grange_expr, strand == "+"),
        geom = "rect", color = colou[1],
        radius = 26, trackWidth = 4) +
    layout_circle(
        subset(x = ref_grange_expr, strand == "-"),
        geom = "rect", color = colou[2],
        radius = 22, trackWidth = 4) +
    layout_circle(
        subset(x = orf_grange_expr, strand == "+"),
        geom = "rect", color = colou[3],
        radius = 18, trackWidth = 4) +
    layout_circle(
        subset(x = orf_grange_expr, strand == "-"),
        geom = "rect", color = colou[4],
        radius = 14, trackWidth = 4)


pdf(file = paste0(opt$out_path, "/ORFs_circos.pdf"), width = 10, height = 10)
plot(pl)
dev.off()


### Genomic coverage -----------------------------------------------------



# Get all chromosome nucleotide position
chrom_nuc <- bsu_ideo@ranges[[1]]

# Get all protein-coding associated nucleotide position
coding_nuc <- lapply(X = ref_grange@ranges, FUN = function(x) {
    x
})
names(coding_nuc) <- ref_grange@elementMetadata@listData$UniProtKBID

# Get all expressed protein associated nucleotide position
exprs_nuc <- lapply(X = ref_grange_expr@ranges, FUN = function(x) {
    x
})
names(exprs_nuc) <- ref_grange_expr@elementMetadata@listData$UniProtKBID
exprs_nuc_stranded <- list()
orf_coord_filt <- orf_coord[!is.na(orf_coord$UniProtID), ] %>%
    dplyr::filter(., !duplicated(UniProtID))
orf_coord_filt$UniProtID %<>%
    uni_id_clean(.)
for (x in names(exprs_nuc)) {
    if (orf_coord_filt[orf_coord_filt$UniProtID == x, "strand"] == -1) {
        exprs_nuc_stranded[[x]] <- rev(exprs_nuc[[x]])
    } else {
        exprs_nuc_stranded[[x]] <- exprs_nuc[[x]]
    }
}

# Get all peptide associated nucleotide position
pep_loc <- readRDS(file = "C:/Users/kxmna01/Dropbox/Home_work_sync/Work/Colleagues shared work/Vaishnavi Ravikumar/Bacillus_subtilis_6frame/15082017/2017-08-15_Peptides_location.RDS") %>%
    ldply(., "data.frame", .id = "Database")
pep_loc_filt <- pep_loc %>%
    dplyr::filter(., Database == "Known")

tmp <- pep_loc_filt[
    !is.na(pep_loc_filt$start) & !duplicated(pep_loc_filt$pep), ]
cover_nuc <- apply(
    X = tmp,
    MARGIN = 1,
    FUN = function(x) {
        
        val <- seq(from = as.integer(x["start"]) * 3 - 2, to = as.integer(x["end"]) * 3)
        tmp <- exprs_nuc_stranded[[x["prot"]]][val]
        tmp
        
    }
)
names(cover_nuc) <- tmp$pep

# Format data into a standard plotting dataframe
toplot <- data.frame(
    Param = paste(
        "Chromosome", round(length(chrom_nuc) / 1000000, 1),
        "Mb", sep = " "),
    Position = as.integer(unique(chrom_nuc)),
    stringsAsFactors = FALSE)
toplot <- data.frame(
    Param = paste(
        "Protein-coding",
        round(length(unique(unlist(coding_nuc))) / 1000000, 1),
        "Mb", sep = " "),
    Position = as.integer(unique(unlist(coding_nuc))),
    stringsAsFactors = FALSE) %>%
    base::rbind(toplot, .)
toplot <- data.frame(
    Param = paste(
        "Expressed protein",
        round(length(unique(unlist(exprs_nuc))) / 1000000, 1),
        "Mb", sep = " "),
    Position = as.integer(unique(unlist(exprs_nuc))),
    stringsAsFactors = FALSE) %>%
    base::rbind(toplot, .)
toplot <- data.frame(
    Param = paste(
        "Detected peptide",
        round(length(unique(unlist(cover_nuc))) / 1000000, 1),
        "Mb", sep = " "),
    Position = as.integer(unique(unlist(cover_nuc))),
    stringsAsFactors = FALSE) %>%
    base::rbind(toplot, .)

pdf(file = paste0(opt$out_path, "/Genomic_coverage.pdf"), width = 10, height = 10)

# Plot a square venn diagram of chromosome coverage
plot.new()
rect(
    xleft = 0,
    ybottom = 0,
    xright = 1,
    ytop = 1,
    border = "red",
    lwd = 2)
text(
    x = 0.015,
    y = 1.02,
    labels = paste(
        "Chromosome", round(length(chrom_nuc) / 1000000, 1),
        "Mb", sep = " "),
    col = "red", cex = 1.0, adj = 0)
rect(
    xleft = 0.01,
    ybottom = 0.01,
    xright = sqrt(
        length(unique(unlist(coding_nuc))) / length(chrom_nuc)) + 0.01,
    ytop = sqrt(
        length(unique(unlist(coding_nuc))) / length(chrom_nuc)) + 0.01,
    border = "green",
    lwd = 2)
text(
    x = 0.025,
    y = sqrt(
        length(unique(unlist(coding_nuc))) / length(chrom_nuc)) + 0.032,
    labels = paste(
        "Protein-coding",
        round(length(unique(unlist(coding_nuc))) / 1000000, 1),
        "Mb", sep = " "),
    col = "green", cex = 1.0, adj = 0)
rect(
    xleft = 0.02,
    ybottom = 0.02,
    xright = sqrt(
        length(unique(unlist(exprs_nuc))) / length(chrom_nuc)) + 0.02,
    ytop = sqrt(
        length(unique(unlist(exprs_nuc))) / length(chrom_nuc)) + 0.02,
    border = "blue",
    lwd = 2)
text(
    x = 0.035,
    y = sqrt(
        length(unique(unlist(exprs_nuc))) / length(chrom_nuc)) + 0.04,
    labels = paste(
        "Expressed protein",
        round(length(unique(unlist(exprs_nuc))) / 1000000, 1),
        "Mb", sep = " "),
    col = "blue", cex = 1.0, adj = 0)
rect(
    xleft = 0.03,
    ybottom = 0.03,
    xright = sqrt(
        length(unique(unlist(cover_nuc))) / length(chrom_nuc)) + 0.03,
    ytop = sqrt(
        length(unique(unlist(cover_nuc))) / length(chrom_nuc)) + 0.03,
    border = "gold",
    lwd = 2)
text(
    x = 0.045,
    y = sqrt(
        length(unique(unlist(cover_nuc))) / length(chrom_nuc)) + 0.05,
    labels = paste(
        "Detected peptide",
        round(length(unique(unlist(cover_nuc))) / 1000000, 1),
        "Mb", sep = " "),
    col = "gold", cex = 1.0, adj = 0)

dev.off()



### Coverage per nucleotide ----------------------------------------------

# 
peptide_to_nucl <- cover_nuc %>%
    ldply(., "data.frame") %>%
    set_colnames(c("Sequence", "Nucl_pos")) %>%
    dplyr::filter(., !is.na(Nucl_pos))

# Merge with evid and the msms count



