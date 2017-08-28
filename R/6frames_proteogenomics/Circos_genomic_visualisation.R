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
    dplyr::left_join(., orf_to_uniprot, by = c("association_id" = "id")) %>%
    dplyr::mutate(., UniProtID = uni_id_clean(UniProtID))


evid_match <- readRDS("C:/Users/kxmna01/Dropbox/Home_work_sync/Work/Colleagues shared work/Vaishnavi Ravikumar/Bacillus_subtilis_6frame/15082017/2017-08-15_Sequence_group_mapping.RDS")
evid_expr_known <- evid_match %>%
    dplyr::filter(., group == "Known") %>%
    strsplit(x = .[["Proteins"]], split = ";", fixed = TRUE) %>%
    unlist(.) %>%
    unique(.)
evid_expr_novel <- evid_match %>%
    dplyr::filter(., group == "Novel") %>%
    strsplit(x = .[["Proteins"]], split = ";", fixed = TRUE) %>%
    unlist(.) %>%
    unique(.)

operon <- read.table(
    file = "C:/Users/kxmna01/Dropbox/Home_work_sync/Work/Colleagues shared work/Vaishnavi Ravikumar/Bacillus_subtilis_6frame/Databases/Bsu_operon_01062017.opr",
    header = TRUE,
    sep = "\t",
    quote = "",
    as.is = TRUE)


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
    dplyr::filter(., !is.na(id)) %>%
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
ref_grange_expr <- subset(x = ref_grange, UniProtKBID %in% evid_expr_known)
orf_grange_expr <- subset(x = orf_grange, id %in% evid_expr_novel)
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




# Format dataframe as genomic position and other info for each novel ORF
tmp <- operon %>%
    dplyr::group_by(., OperonID) %>%
    dplyr::summarise(
        .,
        GeneCount = n(),
        Start = min(Start),
        End = max(End),
        Strand = unique(Strand),
        Length = max(End) - min(Start)) %>%
    dplyr::filter(., !is.na(OperonID)) %>%
    base::as.data.frame(., stringsAsFactors = FALSE)
tmp %<>%
    dplyr::mutate(
        .,
        Chromosome = 1,
        chr.name = "chr1",
        length = 4215606,
        Type = TRUE,
        geno = "AL009126.3") %>%
    base::as.data.frame(., stringsAsFactors = FALSE)
tmp %<>%
    dplyr::filter(., !duplicated(OperonID))

# Create a GRanges object for all reference and identified novel ORF 
operon_grange <- with(
    data = tmp,
    expr = GRanges(
        seqnames = chr.name,
        ranges = IRanges(Start, End),
        strand = Strand))

# Add values and seqinfo to the created GRanges object
values(operon_grange) <- tmp %>%
    dplyr::select(., OperonID, GeneCount, Chromosome)
seqinfo(operon_grange) <- Seqinfo(
    seqnames = tmp$chr.name %>% unique(.) %>% as.character(.),
    seqlengths = tmp$length %>% unique(.) %>% as.integer(.),
    isCircular = tmp$Type %>% unique(.) %>% as.logical(.),
    genome = tmp$geno %>% unique(.) %>% as.character(.))


colou <- c("#4682B4", "#BD5E5E", "#437A3C", "#F57C36", "#D58DEB", "#B2B83F")
pl <- ggplot() +
    ggtitle(label = "Bacillus subtilis ORFs") +
    layout_circle(
        bsu_ideo, geom = "ideo", fill = "#e6e6e6",
        radius = 30.5, trackWidth = 2) +
    layout_circle(
        bsu_ideo, geom = "scale", size = 4,
        radius = 33, trackWidth = 2) +
    layout_circle(
        bsu_ideo, geom = "text", size = 7, aes(label = seqnames),
        angle = 0, radius = 38, trackWidth = 5) +
    layout_circle(
        subset(x = operon_grange, strand == "+" & GeneCount > 2),
        geom = "rect", color = colou[5],
        radius = 30.5, trackWidth = 1.9) +
    layout_circle(
        subset(x = operon_grange, strand == "-" & GeneCount > 2),
        geom = "rect", color = colou[6],
        radius = 30.5, trackWidth = 1.9) +
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


pdf(file = paste0(opt$out_path, "/ORFs_circos_with_operon.pdf"), width = 10, height = 10)
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
peptide_to_nucl_count <- evid_match %>%
    dplyr::select(., Sequence, `MS/MS Count`) %>%
    dplyr::group_by(., Sequence) %>%
    dplyr::summarise(., Count = sum(`MS/MS Count`)) %>%
    dplyr::left_join(peptide_to_nucl, ., by = "Sequence")

# 
peptide_to_nucl_count$Count <- factor(
    x = peptide_to_nucl_count$Count,
    levels = seq(1, max(peptide_to_nucl_count$Count)),
    labels = seq(1, max(peptide_to_nucl_count$Count)),
    ordered = TRUE)
toplot <- peptide_to_nucl_count %>%
    plyr::ddply(
        .data = ., .variables = c("Count"),
        .fun = summarise, Freq = n(),
        .drop = FALSE)

quantiles_toplot <- quantile(
    x = as.integer(as.character(peptide_to_nucl_count$Count)),
    probs = c(0, 0.25, 0.5, 0.75, 1)) %>%
    as.data.frame(.) %>%
    set_colnames("Quartiles")


pl <- plots_hist(
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

pl <- pl[[1]] +
    annotation_custom(
        grob = tableGrob(d = quantiles_toplot, theme = ttheme_minimal()),
        xmin = 130, xmax = 150, ymin = 2E5, ymax = 4E5)


pdf(file = paste0(opt$out_path, "/Coverage_per_nucleotide.pdf"), width = 10, height = 10)
plot(pl)
dev.off()


peptide_count <- peptide_to_nucl_count %>%
    dplyr::select(., -Nucl_pos) %>%
    unique(.)
peptide_count$Count %<>%
    as.character(.) %>%
    as.numeric(.)

peptide_count <- evid_match %>%
    dplyr::select(
        ., Sequence, Length, Modifications, `Modified sequence`,
        starts_with("Missed cleavages"), Proteins, `Leading Proteins`,
        `Leading Razor Protein`, `Gene Names`, `Protein Names`,
        Type, `Labeling State`, `Raw file`, Experiment, `MS/MS m/z`,
        Charge, `m/z`, Mass, Resolution, PEP, `MS/MS Count`, Score,
        starts_with("Intensity"), Reverse, `Potential contaminant`, id,
        group, Database) %>%
    dplyr::left_join(., peptide_count, by = "Sequence")


write.table(
    x = peptide_count,
    file = paste0(opt$out_path, "/Coverage_per_peptide.txt"),
    quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)



### Visualise specific genomic region


start_pos <- start(orf_grange_expr[orf_grange_expr@elementMetadata@listData$id == "seq_51322"]) - 2000
end_pos <- end(orf_grange_expr[orf_grange_expr@elementMetadata@listData$id == "seq_51322"]) + 2000


tmp <- operon_grange[
    start(operon_grange) %in% c(start_pos:end_pos) |
        end(operon_grange) %in% c(start_pos:end_pos)] %>%
    as.data.frame(.)
tmp$value <- rep(x = 1, times = nrow(tmp))

pl1 <- autoplot(
    operon_grange[
        start(operon_grange) %in% c(start_pos:end_pos) |
            end(operon_grange) %in% c(start_pos:end_pos)],
    mapping = aes(fill = strand),
    geom = "arrowrect", layout = "linear", colour = "black") +
    geom_text(
        data = tmp,
        mapping = aes(
            x = start + ((end - start) / 2), y = value, label = OperonID),
        nudge_y = 0.45, check_overlap = TRUE)
    scale_fill_manual(values = c(`+` = colou[5], `-` = colou[6]))
pl1


tmp <- ref_grange_expr[
    start(ref_grange_expr) %in% c(start_pos:end_pos) |
        end(ref_grange_expr) %in% c(start_pos:end_pos)] %>%
    as.data.frame(.)
tmp$value <- rep(x = 1, times = nrow(tmp))

pl2 <- autoplot(
    ref_grange_expr[
        start(ref_grange_expr) %in% c(start_pos:end_pos) |
            end(ref_grange_expr) %in% c(start_pos:end_pos)],
    mapping = aes(fill = strand),
    geom = "arrowrect", layout = "linear", colour = "black") +
    geom_text(
        data = tmp,
        mapping = aes(
            x = start + ((end - start) / 2), y = value, label = UniProtKBID),
        nudge_y = 0.45, check_overlap = TRUE) +
    scale_fill_manual(values = c(`+` = colou[1], `-` = colou[2]))
pl2



tmp <- orf_grange[
    start(orf_grange) %in% c(start_pos:end_pos) |
        end(orf_grange) %in% c(start_pos:end_pos)] %>%
    as.data.frame(.)
tmp$value <- rep(x = 1, times = nrow(tmp))

pl3 <- autoplot(
    orf_grange[
        (start(orf_grange) %in% c(start_pos:end_pos) |
             end(orf_grange) %in% c(start_pos:end_pos)) &
            orf_grange$frame == 1],
    mapping = aes(fill = strand),
    geom = "arrowrect", layout = "linear", colour = "black") +
    scale_fill_manual(values = c(`+` = colou[3], `-` = colou[4]))
pl3
pl4 <- autoplot(
    orf_grange[
        (start(orf_grange) %in% c(start_pos:end_pos) |
            end(orf_grange) %in% c(start_pos:end_pos)) &
            orf_grange$frame == 2],
    mapping = aes(fill = strand),
    geom = "arrowrect", layout = "linear", colour = "black") +
    scale_fill_manual(values = c(`+` = colou[3], `-` = colou[4]))
pl4
pl5 <- autoplot(
    orf_grange[
        (start(orf_grange) %in% c(start_pos:end_pos) |
             end(orf_grange) %in% c(start_pos:end_pos)) &
            orf_grange$frame == 3],
    mapping = aes(fill = strand),
    geom = "arrowrect", layout = "linear", colour = "black") +
    scale_fill_manual(values = c(`+` = colou[3], `-` = colou[4]))
pl5
pl6 <- autoplot(
    orf_grange[
        (start(orf_grange) %in% c(start_pos:end_pos) |
             end(orf_grange) %in% c(start_pos:end_pos)) &
            orf_grange$frame == -1],
    mapping = aes(fill = strand),
    geom = "arrowrect", layout = "linear", colour = "black") +
    scale_fill_manual(values = c(`+` = colou[3], `-` = colou[4]))
pl6
pl7 <- autoplot(
    orf_grange[
        (start(orf_grange) %in% c(start_pos:end_pos) |
             end(orf_grange) %in% c(start_pos:end_pos)) &
            orf_grange$frame == -2],
    mapping = aes(fill = strand),
    geom = "arrowrect", layout = "linear", colour = "black") +
    scale_fill_manual(values = c(`+` = colou[3], `-` = colou[4]))
pl7
pl8 <- autoplot(
    orf_grange[
        (start(orf_grange) %in% c(start_pos:end_pos) |
             end(orf_grange) %in% c(start_pos:end_pos)) &
            orf_grange$frame == -3],
    mapping = aes(fill = strand),
    geom = "arrowrect", layout = "linear", colour = "black") +
    scale_fill_manual(values = c(`+` = colou[3], `-` = colou[4]))
pl8


pl <- tracks(
    `Operon` = pl1, `Known ORFs` = pl2,
    `Frame 1` = pl3, `Frame 2` = pl4, `Frame 3` = pl5,
    `Frame -1` = pl6, `Frame -2` = pl7, `Frame -3` = pl8,
    heights = c(2, 2, rep(0.5, times = 6)),
    xlab = "Genomic position",
    label.bg.fill = "white") +
    theme_bw() +
    theme(
        axis.text.y =  element_blank(),
        axis.ticks.y = element_blank(),
        legend.position = "none")
xlim(pl) <- c(start_pos:end_pos)
pl







# Get all expressed protein associated nucleotide position
exprs_nuc_orf <- lapply(X = ranges(orf_grange_expr), FUN = function(x) {
    x
})
names(exprs_nuc_orf) <- orf_grange_expr$id
exprs_nuc_stranded_orf <- list()
orf_coord_filt <- orf_coord[!is.na(orf_coord$id), ] %>%
    dplyr::filter(., !duplicated(id))
for (x in names(exprs_nuc_orf)) {
    if (orf_coord_filt[orf_coord_filt$id == x, "strand"] == -1) {
        exprs_nuc_stranded_orf[[x]] <- rev(exprs_nuc_orf[[x]])
    } else {
        exprs_nuc_stranded_orf[[x]] <- exprs_nuc_orf[[x]]
    }
}


all_exprs_nuc <- c(exprs_nuc_stranded, exprs_nuc_stranded_orf)
tmp <- pep_loc_filt[
    !is.na(pep_loc_filt$start) & !duplicated(pep_loc_filt$pep), ]
cover_nuc <- apply(
    X = tmp,
    MARGIN = 1,
    FUN = function(x) {
        
        val <- c(as.integer(x["start"]) * 3 - 2, as.integer(x["end"]) * 3)
        tmp <- c(
            x[["pep"]],
            all_exprs_nuc[[x["prot"]]][val[1]],
            all_exprs_nuc[[x["prot"]]][val[2]],
            x[["prot"]],
            x[["Database"]])
        if (length(tmp) < 4) {
            warning(paste0("Wrong length: ", paste(tmp, collapse = ";")))
            tmp <- c(
                x[["pep"]],
                NA,
                NA,
                x[["prot"]],
                x[["Database"]])
        }
        names(tmp) <- c("pep", "Start", "End", "Prot", "Database")
        tmp
        
    }
)
cover_nuc <- t(cover_nuc) %>%
    as.data.frame(., stringAsFactors = FALSE) %>%
    dplyr::mutate(
        .,
        Start = as.integer(as.character(Start)),
        End = as.integer(as.character(End)))
cover_nuc_final <- cover_nuc %>%
    dplyr::filter(., Database == "Known") %>%
    dplyr::left_join(
        x = .,
        y = orf_coord[!is.na(orf_coord$UniProtID), ],
        by = c("Prot" = "UniProtID")) %>%
    dplyr::mutate(., UniProtID = Prot)
cover_nuc_final <- cover_nuc %>%
    dplyr::filter(., Database == "Novel") %>%
    dplyr::left_join(
        x = .,
        y = orf_coord[!is.na(orf_coord$id), ],
        by = c("Prot" = "id")) %>%
    dplyr::mutate(., id = Prot) %>%
    dplyr::bind_rows(cover_nuc_final, .)

# Format dataframe as genomic position and other info for each peptide
tmp <- cover_nuc_final %>%
    dplyr::filter(., !is.na(pep) & !is.na(Start) & !is.na(End)) %>%
    base::as.data.frame(., stringsAsFactors = FALSE)
tmp %<>%
    dplyr::select(., pep, Start, End, Prot, Database, strand, frame) %>%
    dplyr::mutate(
        .,
        start = ifelse(test = strand == 1, yes = Start, no = End),
        end = ifelse(test = strand == 1, yes = End, no = Start),
        Strand = ifelse(test = strand == 1, yes = "+", no = "-"),
        ProteinID = Prot,
        Chromosome = 1,
        chr.name = "chr1",
        length = 4215606,
        Type = TRUE,
        geno = "AL009126.3") %>%
    dplyr::select(., -Start, -End, -Prot, -strand) %>%
    base::as.data.frame(., stringsAsFactors = FALSE)
tmp %<>%
    dplyr::filter(., !duplicated(pep))

# Create a GRanges object for all reference and identified novel ORF 
pep_grange <- with(
    data = tmp,
    expr = GRanges(
        seqnames = chr.name,
        ranges = IRanges(start, end),
        strand = Strand))

# Add values and seqinfo to the created GRanges object
values(pep_grange) <- tmp %>%
    dplyr::select(., pep, Database, frame, Chromosome, ProteinID)
seqinfo(pep_grange) <- Seqinfo(
    seqnames = tmp$chr.name %>% unique(.) %>% as.character(.),
    seqlengths = tmp$length %>% unique(.) %>% as.integer(.),
    isCircular = tmp$Type %>% unique(.) %>% as.logical(.),
    genome = tmp$geno %>% unique(.) %>% as.character(.))
