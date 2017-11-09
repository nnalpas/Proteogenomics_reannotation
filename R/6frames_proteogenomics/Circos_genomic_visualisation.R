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
load_package("plyr")
load_package("dplyr")
load_package("tidyr")
load_package("seqinr")
load_package("UniProt.ws")
load_package("magrittr")
load_package("WriteXLS")
load_package("data.table")
load_package("splitstackshape")
load_package("VennDiagram")
load_package("ggplot2")
load_package("grid")
load_package("gridExtra")
load_package("RColorBrewer")
load_package("stringr")
load_package("Biostrings")
load_package("RecordLinkage")
load_package("VariantAnnotation")
load_package("cgdsr")
load_package("bit64")
load_package("cleaver")
load_package("plotly")
load_package("GenomicRanges")
load_package("biovizBase")
load_package("ggbio")
load_package("ggradar")
load_package("BSgenome.Bsubtilis.EMBL.AL0091263")
load_package("GenomicFeatures")



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
#    out_path = "C:/Users/kxmna01/Dropbox/Home_work_sync/Work/Colleagues shared work/Vaishnavi Ravikumar/Bacillus_subtilis_6frame/Genome_plot")



### Data import ----------------------------------------------------------

# Import the orf coordinates
orf_coord <- read.table(
    file = opt$coordinates, header = TRUE,
    sep = "\t", quote = "", as.is = TRUE)

# Genome size
geno_size <- 4215606



# temporary



# Import the experimental design (with conditions)
exp_design <- read.table(
    file = "C:/Users/kxmna01/Dropbox/Home_work_sync/Work/Colleagues shared work/Vaishnavi Ravikumar/Bacillus_subtilis_6frame/Bsu_Conditions_06022017.txt", header = TRUE,
    sep = "\t", quote = "", as.is = TRUE)

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
    data = blast_NN_vs_ref, key = "qseqid", multi_match = "uniquify")

# Keep the mapping of ORF to uniprotID
orf_to_uniprot <- blast_NN_vs_ref_best %>%
    dplyr::select(., qseqid, sseqid) %>%
    set_colnames(c("association_id", "UniProtID")) %>%
    dplyr::mutate(
        .,
        association_id = gsub("^lcl\\|", "", association_id),
        UniProtID = uni_id_clean(UniProtID)) %>%
    base::as.data.frame(., stringsAsFactors = FALSE)

# Include this ID cross-mapping to the ORF coordinates
orf_coord %<>%
    dplyr::left_join(., orf_to_uniprot)


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

ref_details <- read.table(
    file = "C:/Users/kxmna01/Dropbox/Home_work_sync/Work/Colleagues shared work/Vaishnavi Ravikumar/Bacillus_subtilis_6frame/Databases/uniprot-taxonomy-Bsu.tab",
    header = TRUE,
    sep = "\t",
    quote = "",
    as.is = TRUE) %>%
    dplyr::mutate(
        .,
        Gene_name = sub("([^ ]*).*", "\\1", `Gene.names`),
        Locus = sub(".* (.*)$", "\\1", `Gene.names`))

orf_coord %<>%
    dplyr::left_join(., ref_details, by = c("UniProtID" = "Entry"))

# Import peptide location within proteins
pep_loc <- readRDS(file = "C:/Users/kxmna01/Dropbox/Home_work_sync/Work/Colleagues shared work/Vaishnavi Ravikumar/Bacillus_subtilis_6frame/15082017/2017-08-15_Peptides_location.RDS") %>%
    ldply(., "data.frame", .id = "Database")

# Import the blast result of PCR sanger sequence against Bsu genome
sanger_vs_genome <- blast_read(
    file = "C:/Users/kxmna01/Dropbox/Home_work_sync/Work/Colleagues shared work/Vaishnavi Ravikumar/PCR_validation/Sanger_seq/Sanger_vs_Genome",
    blast_format = "6")

# Identify best blast for each sanger sequence
best_sanger_vs_genome <- best_blast(data = sanger_vs_genome, key = "qseqid")

# Import the blast results of Nicolas' ORF versus the reference proteome
tblastn_ref_vs_genome <- blast_read(
    file = "C:/Users/kxmna01/Dropbox/Home_work_sync/Work/Colleagues shared work/Vaishnavi Ravikumar/Bacillus_subtilis_6frame/tblastn_Refprot_vs_Genome_08112017",
    blast_format = "6")

# Get the best blast and format sequence as stringset
tblastn_ref_vs_genome_best <- best_blast(
    data = tblastn_ref_vs_genome, key = "qseqid", multi_match = "uniquify") %>%
    dplyr::rowwise() %>%
    dplyr::mutate(
        .,
        qseq = list(AAString(qseq)),
        sseq = list(AAString(sseq))) %<>%
    dplyr::left_join(., ref_details, by = c("qseqid" = "Entry"))

# Define genomic strand
tblastn_ref_vs_genome_best %<>%
    dplyr::mutate(
        .,
        ExactCoord = ifelse(
            qstart == 1 & qend == qlen,
            "Exact",
            "Approximate"),
        strand = ifelse(
            test = sstart < send, yes = 1, no = -1))

# Define translation frame
tblastn_ref_vs_genome_best <- get_frame(
    data = tblastn_ref_vs_genome_best,
    start = "sstart",
    strand = "strand",
    genome_size = geno_size)

# Adjust genomic coordinate if blast is not from start to end of query
tblastn_ref_vs_genome_best %<>%
    dplyr::mutate(
        .,
        start = ifelse(
            strand == "1",
            (sstart - ((qstart - 1) * 3)),
            (sstart + ((qstart - 1) * 3))),
        end = ifelse(
            strand == "1",
            (send + ((qlen - qend) * 3)),
            (send - ((qlen - qend) * 3))),
        nucl_length_qlen = qlen * 3,
        nucl_length_str_end = abs(end - start) + 1,
        nucl_length = nucl_length_str_end,
        Comment = ifelse(
            nucl_length_qlen != nucl_length_str_end, "Warning", "OK")) %>%
    base::as.data.frame(., stringsAsFactors = FALSE)

# Include info on the ORF ID 
tblastn_ref_vs_genome_best <- orf_coord %>%
    dplyr::select(., "association_id", "id", "UniProtID") %>%
    unique(.) %>%
    dplyr::left_join(
        x = tblastn_ref_vs_genome_best, y = ., by = c("qseqid" = "UniProtID")) %>%
    dplyr::group_by(., qseqid) %>%
    dplyr::summarise_all(funs(toString(x = unique(.), width = NULL)))



### Ideogram generation --------------------------------------------------

# Dataframe holding genome information for bacillus subtilis
tmp <- base::data.frame(
    Chromosome = 1, Strand = "*", Start = 1, End = geno_size,
    name = "chr1", length = geno_size, Type = TRUE, geno = "AL0091263",
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



### Create GRange of ORFs and other entries ------------------------------

# Format dataframe as genomic position and other info for each reference ORF
tmp <- tblastn_ref_vs_genome_best %>%
    dplyr::select(
        ., qseqid, association_id, id, start, end, strand, frame,
        `Entry.name`, Status, `Protein.names`, `Gene.names`, Organism,
        Length, `Subcellular.location..CC.`, Gene_name, Locus,
        ExactCoord, qlen, nucl_length, Comment) %>%
    plyr::rename(., c("qseqid" = "UniProtID", "qlen" = "aa_length")) %>%
    base::as.data.frame(., stringsAsFactors = FALSE)
tmp %<>%
    dplyr::mutate(
        .,
        Start = as.integer(ifelse(test = strand == 1, yes = start, no = end)),
        End = as.integer(ifelse(test = strand == 1, yes = end, no = start)),
        Strand = ifelse(test = strand == 1, yes = "+", no = "-"),
        UniProtKBID = UniProtID,
        Expressed = ifelse(UniProtID %in% evid_expr_known, TRUE, FALSE),
        Chromosome = 1,
        chr.name = "chr1",
        length = geno_size,
        Type = TRUE,
        geno = "AL0091263") %>%
    dplyr::select(
        .,
        UniProtKBID, id, association_id, Start, End, Strand, frame, Expressed,
        Chromosome, `chr.name`, length, Type, geno, Gene_name, Locus,
        ExactCoord, Comment) %>%
    dplyr::group_by(
        ., UniProtKBID, Start, End, Strand, frame, Chromosome, chr.name) %>%
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
    dplyr::select(
        ., id, UniProtKBID, association_id, frame, Expressed,
        Chromosome, Gene_name, Locus, ExactCoord, Comment)
seqinfo(ref_grange) <- Seqinfo(
    seqnames = tmp$chr.name %>% unique(.) %>% as.character(.),
    seqlengths = tmp$length %>% unique(.) %>% as.integer(.),
    isCircular = tmp$Type %>% unique(.) %>% as.logical(.),
    genome = tmp$geno %>% unique(.) %>% as.character(.))
names(ref_grange) <- ref_grange$UniProtKBID

# Export the reference GRanges for reuse at later stage
saveRDS(
    object = ref_grange,
    file = paste0(opt$out_path, "/bsu_AL009126.3_ref_grange.RDS"))

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
        UniProtKBID = UniProtID,
        Expressed = ifelse(UniProtID %in% evid_expr_known | id %in% evid_expr_novel, TRUE, FALSE),
        Chromosome = 1,
        chr.name = "chr1",
        length = geno_size,
        Type = TRUE,
        geno = "AL0091263") %>%
    dplyr::select(
        .,
        UniProtKBID, id, association_id, Start, End, Strand, frame, Expressed,
        Chromosome, `chr.name`, length, Type, geno, Gene_name, Locus) %>%
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
    dplyr::select(
        ., id, UniProtKBID, association_id, frame, Expressed,
        Chromosome, Gene_name, Locus)
seqinfo(orf_grange) <- Seqinfo(
    seqnames = tmp$chr.name %>% unique(.) %>% as.character(.),
    seqlengths = tmp$length %>% unique(.) %>% as.integer(.),
    isCircular = tmp$Type %>% unique(.) %>% as.logical(.),
    genome = tmp$geno %>% unique(.) %>% as.character(.))
names(orf_grange) <- orf_grange$id

# Export the reference GRanges for reuse at later stage
saveRDS(
    object = orf_grange,
    file = paste0(opt$out_path, "/bsu_AL009126.3_orf_grange.RDS"))

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
    dplyr::filter(., !is.na(OperonID))%>%
    dplyr::mutate(
        .,
        Chromosome = 1,
        chr.name = "chr1",
        length = geno_size,
        Type = TRUE,
        geno = "AL0091263") %>%
    base::as.data.frame(., stringsAsFactors = FALSE)
tmp %<>%
    dplyr::filter(., !duplicated(OperonID))

# Create a GRanges object for all operons
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
names(operon_grange) <- operon_grange$OperonID

# Export the reference GRanges for reuse at later stage
saveRDS(
    object = operon_grange,
    file = paste0(opt$out_path, "/bsu_AL009126.3_operon_grange.RDS"))

# Get all known protein associated nucleotide position
tmp <- lapply(
    X = ranges(ref_grange),
    FUN = function(x) {
        x
    }) %>%
    tibble::tibble(UniProtKBID = names(.), nucl_pos = .)
tmp <- lapply(
    X = tmp$nucl_pos,
    FUN = function(x) {
        rev(x)
    }) %>%
    tibble::tibble(UniProtKBID = tmp$UniProtKBID, rev_nucl_pos = .) %>%
    dplyr::left_join(tmp, .)
tmp <- ref_grange %>%
    tibble::as.tibble(.) %>%
    dplyr::left_join(., tmp)
nuc_stranded_ref <- tmp %>%
    dplyr::mutate(
        .,
        mainID = UniProtKBID,
        nucl_pos_stranded = case_when(
            as.character(.[["strand"]]) == "+" ~ .[["nucl_pos"]],
            as.character(.[["strand"]]) == "-" ~ .[["rev_nucl_pos"]]))

# Get all ORF associated nucleotide position
tmp <- lapply(
    X = ranges(orf_grange),
    FUN = function(x) {
        x
    }) %>%
    tibble::tibble(id = names(.), nucl_pos = .)
tmp <- lapply(
    X = tmp$nucl_pos,
    FUN = function(x) {
        rev(x)
    }) %>%
    tibble::tibble(id = tmp$id, rev_nucl_pos = .) %>%
    dplyr::left_join(tmp, .)
tmp <- orf_grange %>%
    tibble::as.tibble(.) %>%
    dplyr::left_join(., tmp)
nuc_stranded_orf <- tmp %>%
    dplyr::mutate(
        .,
        mainID = id,
        nucl_pos_stranded = case_when(
            as.character(.[["strand"]]) == "+" ~ .[["nucl_pos"]],
            as.character(.[["strand"]]) == "-" ~ .[["rev_nucl_pos"]]))

# Compile all nucleotide position of all known genes and ORFs
nuc_stranded_all <- dplyr::bind_rows(nuc_stranded_ref, nuc_stranded_orf)

# Associate each peptide to its genes/orfs
nuc_stranded_pep <- pep_loc %>%
    set_colnames(c("Database", "pep", "mainID", "pep_start", "pep_end")) %>%
    dplyr::left_join(nuc_stranded_all, .) %>%
    dplyr::filter(., !is.na(pep) & !is.na(pep_start) & !is.na(pep_end))

# For each peptide transform aa position to nucleotide position
nuc_stranded_pep %<>%
    dplyr::mutate(
        .,
        pep_nucl_start = as.integer(pep_start) * 3 - 2,
        pep_nucl_end = as.integer(pep_end) * 3) %>%
    rowwise(.) %>%
    dplyr::mutate(
        .,
        pep_nucl_pos_start = nucl_pos_stranded[pep_nucl_start],
        pep_nucl_pos_end = nucl_pos_stranded[pep_nucl_end])

# Get vector of all peptide coming from known proteins
known_pep <- nuc_stranded_pep %>%
    dplyr::filter(., Database == "Known") %>%
    .[["pep"]] %>%
    unique(.)

# Keep peptide entries as coming exclusively from Known database and
# the remaining as Novel
tmp <- nuc_stranded_pep %>%
    dplyr::filter(
        .,
        (pep %in% known_pep & Database == "Known") |
            (!(pep %in% known_pep) & Database == "Novel"))

# Format dataframe as genomic position and other info for each peptide
tmp %<>%
    dplyr::filter(
        ., !is.na(pep) &
            !is.na(pep_nucl_pos_start) &
            !is.na(pep_nucl_pos_end)) %>%
    base::as.data.frame(., stringsAsFactors = FALSE)
tmp %<>%
    dplyr::select(
        ., pep, pep_nucl_pos_start, pep_nucl_pos_end,
        mainID, id, UniProtKBID, Database, strand, frame) %>%
    dplyr::mutate(
        .,
        start = ifelse(
            test = strand == "+",
            yes = pep_nucl_pos_start,
            no = pep_nucl_pos_end),
        end = ifelse(
            test = strand == "+",
            yes = pep_nucl_pos_end,
            no = pep_nucl_pos_start),
        Chromosome = 1,
        chr.name = "chr1",
        length = geno_size,
        Type = TRUE,
        geno = "AL0091263") %>%
    dplyr::select(., -pep_nucl_pos_start, -pep_nucl_pos_end) %>%
    base::as.data.frame(., stringsAsFactors = FALSE)
tmp %<>%
    dplyr::filter(., !duplicated(pep))

# Create a GRanges object for all reference and identified novel ORF 
pep_grange <- with(
    data = tmp,
    expr = GRanges(
        seqnames = chr.name,
        ranges = IRanges(start, end),
        strand = strand))

# Add values and seqinfo to the created GRanges object
values(pep_grange) <- tmp %>%
    dplyr::select(., pep, mainID, id, UniProtKBID, Database, frame, Chromosome)
seqinfo(pep_grange) <- Seqinfo(
    seqnames = tmp$chr.name %>% unique(.) %>% as.character(.),
    seqlengths = tmp$length %>% unique(.) %>% as.integer(.),
    isCircular = tmp$Type %>% unique(.) %>% as.logical(.),
    genome = tmp$geno %>% unique(.) %>% as.character(.))
names(pep_grange) <- pep_grange$pep

# Export the reference GRanges for reuse at later stage
saveRDS(
    object = pep_grange,
    file = paste0(opt$out_path, "/bsu_AL009126.3_pep_grange.RDS"))

# Format dataframe as genomic position and other info for each sanger sequence
tmp <- best_sanger_vs_genome %>%
    dplyr::select(., qseqid, sstart, send) %>%
    set_colnames(c("sangerid", "sstart", "send")) %>%
    dplyr::mutate(
        .,
        Start = ifelse(sstart < send, sstart, send),
        End = ifelse(sstart < send, send, sstart),
        Chromosome = 1,
        chr.name = "chr1",
        length = geno_size,
        Type = TRUE,
        geno = "AL0091263") %>%
    dplyr::select(., -sstart, -send) %>%
    base::as.data.frame(., stringsAsFactors = FALSE)

# Create a GRanges object for all reference and identified novel ORF 
sanger_grange <- with(
    data = tmp,
    expr = GRanges(
        seqnames = chr.name,
        ranges = IRanges(Start, End)))

# Add values and seqinfo to the created GRanges object
values(sanger_grange) <- tmp %>%
    dplyr::select(., sangerid, Chromosome)
seqinfo(sanger_grange) <- Seqinfo(
    seqnames = tmp$chr.name %>% unique(.) %>% as.character(.),
    seqlengths = tmp$length %>% unique(.) %>% as.integer(.),
    isCircular = tmp$Type %>% unique(.) %>% as.logical(.),
    genome = tmp$geno %>% unique(.) %>% as.character(.))
names(sanger_grange) <- sanger_grange$sangerid

# Export the reference GRanges for reuse at later stage
saveRDS(
    object = sanger_grange,
    file = paste0(opt$out_path, "/bsu_AL009126.3_sanger_grange.RDS"))



### Circos visualisation -------------------------------------------------

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

# Open pdf output and render the plot
pdf(file = paste0(opt$out_path, "/ORFs_circos.pdf"), width = 10, height = 10)
plot(pl)

# Use ggbio extension to plot ORF and operon location on genome
# as a circos graph
pl <- ggplot() +
    ggtitle(label = "Bacillus subtilis ORFs with operons") +
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

# Render the plot
plot(pl)
dev.off()



### Genomic coverage -----------------------------------------------------

# Compute and visualise the genomic coverage based on nucleotide coverage
pl <- plots_rectvenn(
    ideo = bsu_ideo, ref = ref_grange, pep = pep_grange)



### Coverage per nucleotide ----------------------------------------------

# Get all peptide associated nucleotide position
cover_nuc <- gr_nucleotide_pos(
    grange = pep_grange, filter = 'Database == "Known"')

# Convert to dataframe
peptide_to_nucl <- cover_nuc %>%
    ldply(., "data.frame") %>%
    set_colnames(c("Sequence", "Nucl_pos")) %>%
    dplyr::filter(., !is.na(Nucl_pos))

# Compute msms count for each peptide sequence
peptide_to_nucl_count <- evid_match %>%
    dplyr::select(., Sequence, `MS/MS Count`) %>%
    dplyr::group_by(., Sequence) %>%
    dplyr::summarise(., Count = sum(`MS/MS Count`)) %>%
    dplyr::left_join(peptide_to_nucl, ., by = "Sequence")

# Format dataframe count column to factor
peptide_to_nucl_count$Count <- factor(
    x = peptide_to_nucl_count$Count,
    levels = seq(1, max(peptide_to_nucl_count$Count)),
    labels = seq(1, max(peptide_to_nucl_count$Count)),
    ordered = TRUE)

# Calculate overall nucleotide coverage frequencies
toplot <- peptide_to_nucl_count %>%
    plyr::ddply(
        .data = ., .variables = c("Count"),
        .fun = summarise, Freq = n(),
        .drop = FALSE)

# Calculate the quantiles of count frequencies
quantiles_toplot <- quantile(
    x = as.integer(as.character(peptide_to_nucl_count$Count)),
    probs = c(0, 0.25, 0.5, 0.75, 1)) %>%
    as.data.frame(.) %>%
    set_colnames("Values") %>%
    dplyr::mutate(., Quartiles = row.names(.))
quantiles_toplot <- data.frame(
    Quartiles = "Mean",
    Values = mean(as.integer(as.character(peptide_to_nucl_count$Count)))) %>%
    dplyr::bind_rows(quantiles_toplot, .) %>%
    dplyr::select(., Quartiles, Values) %>%
    dplyr::mutate(., Values = round(x = Values, digits = 1))

# Generate the histogram frequency of MS/MS count per nucleotide
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

# Add quantiles table to the graph
pl <- pl[[1]] +
    annotation_custom(
        grob = tableGrob(
            d = quantiles_toplot, theme = ttheme_minimal(), rows = NULL),
        xmin = 130, xmax = 150, ymin = 2E5, ymax = 4E5)

# Output the plot to pdf
pdf(
    file = paste0(opt$out_path, "/", date_str, "_Coverage_per_nucleotide.pdf"),
    width = 10, height = 10)
plot(pl)
dev.off()

# Get only MS/MS count per peptide
peptide_count <- peptide_to_nucl_count %>%
    dplyr::select(., -Nucl_pos) %>%
    unique(.)
peptide_count$Count %<>%
    as.character(.) %>%
    as.numeric(.)

# Add the MS/MS count information to the evidence so that
# they can be sorted by abundance
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

# Export the evidence table together with MS/MS count values
write.table(
    x = peptide_count,
    file = paste0(opt$out_path, "/Coverage_per_peptide.txt"),
    quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)



### Visualise specific genomic region ------------------------------------

# Get the Bsu genome object
bsu <- BSgenome.Bsubtilis.EMBL.AL0091263

# Define the novel ORF of interest
Target_id <- "seq_223100"

# Define genomic region of interest
start_pos <- start(orf_grange_expr[Target_id]) - 1000
end_pos <- end(orf_grange_expr[Target_id]) + 1000

# The plot of operon for that genomic region
#tmp <- operon_grange[
#    (start(operon_grange) %in% c(start_pos:end_pos) |
#        end(operon_grange) %in% c(start_pos:end_pos)) &
#        operon_grange$GeneCount > 2] %>%
#    as.data.frame(.)
#tmp$value <- rep(x = 1, times = nrow(tmp))

#pl1 <- autoplot(
#    operon_grange[
#        (start(operon_grange) %in% c(start_pos:end_pos) |
#            end(operon_grange) %in% c(start_pos:end_pos)) &
#            operon_grange$GeneCount > 2],
#    mapping = aes(fill = strand),
#    geom = "arrowrect", layout = "linear", colour = "black") +
#    geom_text(
#        data = tmp,
#        mapping = aes(
#            x = start + ((end - start) / 2), y = value, label = OperonID),
#        nudge_y = 0.45, check_overlap = TRUE) + 
#    scale_fill_manual(values = c(`+` = colou[5], `-` = colou[6]))
#pl1

# The plot of known genes for that genomic region
tmp <- ref_grange[
    start(ref_grange) %in% c(start_pos:end_pos) |
        end(ref_grange) %in% c(start_pos:end_pos)] %>%
    as.data.frame(.)
tmp$value <- rep(x = 1, times = nrow(tmp))

pl2 <- autoplot(
    ref_grange[
        start(ref_grange) %in% c(start_pos:end_pos) |
            end(ref_grange) %in% c(start_pos:end_pos)],
    mapping = aes(fill = strand, alpha = Expressed),
    geom = "arrowrect", layout = "linear", colour = "black") +
    geom_text(
        data = tmp,
        mapping = aes(
            x = start + ((end - start) / 2), y = value, label = Gene_name),
        nudge_y = 0.45, check_overlap = TRUE) +
    scale_fill_manual(values = c(`+` = colou[1], `-` = colou[2])) +
    scale_alpha_manual(values = c("TRUE" = 1, "FALSE" = 0.5))
pl2

# The plot of novel ORFs for that genomic region
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
    mapping = aes(fill = strand, alpha = Expressed),
    geom = "arrowrect", layout = "linear", colour = "black") +
    scale_fill_manual(values = c(`+` = colou[3], `-` = colou[4])) +
    scale_alpha_manual(values = c("TRUE" = 1, "FALSE" = 0.5))
pl3
pl4 <- autoplot(
    orf_grange[
        (start(orf_grange) %in% c(start_pos:end_pos) |
            end(orf_grange) %in% c(start_pos:end_pos)) &
            orf_grange$frame == 2],
    mapping = aes(fill = strand, alpha = Expressed),
    geom = "arrowrect", layout = "linear", colour = "black") +
    scale_fill_manual(values = c(`+` = colou[3], `-` = colou[4])) +
    scale_alpha_manual(values = c("TRUE" = 1, "FALSE" = 0.5))
pl4
pl5 <- autoplot(
    orf_grange[
        (start(orf_grange) %in% c(start_pos:end_pos) |
             end(orf_grange) %in% c(start_pos:end_pos)) &
            orf_grange$frame == 3],
    mapping = aes(fill = strand, alpha = Expressed),
    geom = "arrowrect", layout = "linear", colour = "black") +
    scale_fill_manual(values = c(`+` = colou[3], `-` = colou[4])) +
    scale_alpha_manual(values = c("TRUE" = 1, "FALSE" = 0.5))
pl5
pl6 <- autoplot(
    orf_grange[
        (start(orf_grange) %in% c(start_pos:end_pos) |
             end(orf_grange) %in% c(start_pos:end_pos)) &
            orf_grange$frame == -1],
    mapping = aes(fill = strand, alpha = Expressed),
    geom = "arrowrect", layout = "linear", colour = "black") +
    scale_fill_manual(values = c(`+` = colou[3], `-` = colou[4])) +
    scale_alpha_manual(values = c("TRUE" = 1, "FALSE" = 0.5))
pl6
pl7 <- autoplot(
    orf_grange[
        (start(orf_grange) %in% c(start_pos:end_pos) |
             end(orf_grange) %in% c(start_pos:end_pos)) &
            orf_grange$frame == -2],
    mapping = aes(fill = strand, alpha = Expressed),
    geom = "arrowrect", layout = "linear", colour = "black") +
    scale_fill_manual(values = c(`+` = colou[3], `-` = colou[4])) +
    scale_alpha_manual(values = c("TRUE" = 1, "FALSE" = 0.5))
pl7
pl8 <- autoplot(
    orf_grange[
        (start(orf_grange) %in% c(start_pos:end_pos) |
             end(orf_grange) %in% c(start_pos:end_pos)) &
            orf_grange$frame == -3],
    mapping = aes(fill = strand, alpha = Expressed),
    geom = "arrowrect", layout = "linear", colour = "black") +
    scale_fill_manual(values = c(`+` = colou[3], `-` = colou[4])) +
    scale_alpha_manual(values = c("TRUE" = 1, "FALSE" = 0.5))
pl8



# Define genomic region of interest
start_pep <- start(orf_grange_expr[Target_id])
end_pep <- end(orf_grange_expr[Target_id])




tmp <- pep_grange[
    start(pep_grange) %in% c(start_pep:end_pep) |
        end(pep_grange) %in% c(start_pep:end_pep)] %>%
    as.data.frame(.)
tmp$value <- rep(x = 1, times = nrow(tmp))
tmp %<>%
    dplyr::arrange(., frame)
tmp$pep <- factor(
    x = tmp$pep,
    levels = as.character(tmp$pep),
    labels = as.character(tmp$pep),
    ordered = TRUE)



pl9 <- ggplot(data = tmp,
    mapping = aes(
        xmin = start, xmax = end,
        ymin = as.integer(pep), ymax = as.integer(pep)+0.5,
        fill = factor(frame), label = pep)) +
    geom_rect(colour = "black") +
    geom_text(mapping = aes(
        x = start + ((end - start) / 2), y = as.integer(pep)+0.275),
        colour = "white", size = 3, check_overlap = TRUE)
pl9



tmp <- sanger_grange[
    start(sanger_grange) %in% c(start_pep:end_pep) |
        end(sanger_grange) %in% c(start_pep:end_pep)] %>%
    as.data.frame(.)

pl10 <- ggplot(data = tmp,
              mapping = aes(
                  xmin = start, xmax = end,
                  ymin = as.integer(factor(sangerid)),
                  ymax = as.integer(factor(sangerid))+0.5,
                  colour = strand)) +
    geom_rect(fill = "grey")
pl10




wh <- orf_grange_expr[Target_id] %>%
    range(.)
pl11 <- autoplot(bsu, which = wh, geom = "rect")
pl11



# Plot the tracks
pl <- tracks(
    #`Operon` = pl1,
    `Known ORFs` = pl2,
    `Frame 1` = pl3, `Frame 2` = pl4, `Frame 3` = pl5,
    `Frame -1` = pl6, `Frame -2` = pl7, `Frame -3` = pl8,
    heights = c(
        #0.2,
        0.5, rep(0.1, times = 6)),
    xlab = "Genomic position",
    label.bg.fill = "white") +
    theme_bw() +
    theme(
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.y =  element_blank(),
        axis.ticks.y = element_blank(),
        legend.position = "none")
xlim(pl) <- c(start_pos:end_pos)
pl





pl_bis <- tracks(
    `Peptide` = pl9, `Sequenced PCR` = pl10, `Genome` = pl11,
    heights = c(2, 0.5, 0.2),
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
xlim(pl_bis) <- c(start_pep, end_pep)
pl_bis




pdf(
    file = paste0(opt$out_path, "/", Target_id, "_genomic_region.pdf"),
    width = 20, height = 10)
pl
pl_bis
dev.off()



### Temporary code -------------------------------------------------------

to_import <- list.files(
    opt$out_path, pattern = "^bsu.*\\.RDS$", full.names = T)

to_import %<>%
    set_names(sub("bsu_AL009126.3_(.*)\\.RDS", "\\1", basename(.)))

data <- purrr::map(.x = to_import, .f = readRDS)

list2env(x = data, envir = .GlobalEnv)
rm(data)

evid_match <- readRDS(
    "C:/Users/kxmna01/Dropbox/Home_work_sync/Work/Colleagues shared work/Vaishnavi Ravikumar/Bacillus_subtilis_6frame/15082017/2017-08-15_Sequence_group_mapping.RDS")

# Get the Bsu genome object
bsu <- BSgenome.Bsubtilis.EMBL.AL0091263

# Define the novel ORF of interest
Target_ids <- c(
    "seq_49263", "seq_51322", "seq_134853", "seq_145510",
    "seq_154909", "seq_163507", "seq_223100")

# Loop through all targets
for (Target_id in Target_ids) {
    
    # Define genomic region of interest
    start_pos <- start(orf_grange[Target_id]) - 1000
    end_pos <- end(orf_grange[Target_id]) + 1000
    
    # Select all reference within range of the target
    target_ref <- ref_grange[
        (start(ref_grange) %in% c(start_pos:end_pos) |
             end(ref_grange) %in% c(start_pos:end_pos)) &
            ref_grange$Expressed == TRUE]
    
    # Select all ORF within range of the target
    target_orf <- orf_grange[
        (start(orf_grange) %in% c(start_pos:end_pos) |
             end(orf_grange) %in% c(start_pos:end_pos)) &
            orf_grange$Expressed == TRUE & orf_grange$UniProtKBID == "NA"]$id
    
    # Get raw files associated to candidate novel ORF
    target_raw <- evid_match %>%
        dplyr::filter(
            .,
            grepl(
                paste0("(", paste(target_orf, collapse = "|"), ")"),
                Proteins)) %>%
        .[["Raw file"]] %>%
        unique(.)
    
    # Determine in which condition was the candidate novel ORF observed
    target_raw_condition <- exp_design %>%
        dplyr::filter(., Name %in% target_raw) %>%
        dplyr::select(., Name, Conditions)
    
    # Create the pdf to contain plots
    pdf(
        file = paste0(
            opt$out_path, "/", date_str, "_", Target_id, "_spectral_count.pdf"),
        width = 10, height = 10)
    
    # Generate the table of condition per raw file of interest
    pl <- ggplot() +
        theme_bw() +
        annotation_custom(
            grob = tableGrob(
                d = target_raw_condition, theme = ttheme_default(), rows = NULL))
    plot(pl)
    
    # Compute the MS/MS count for all known protein and ORF within range of
    # the target and for each raw file
    data <- evid_match %>%
        dplyr::filter(
            .,
            grepl(
                paste0(
                    "(",
                    paste(c(target_ref$UniProtKBID, target_orf), collapse = "|"),
                    ")"),
                Proteins)) %>%
        dplyr::group_by(., `Leading Proteins`, `Raw file`) %>%
        dplyr::summarise(., Count = sum(`MS/MS Count`)) %>%
        dplyr::mutate(
            ., Type = ifelse(`Raw file` %in% target_raw, "Novel", "Known")) %>%
        set_colnames(c("Proteins", "Raw", "Count", "Type")) %>%
        dplyr::mutate(
            .,
            Name = ifelse(
                Proteins %in% target_ref$UniProtKBID,
                target_ref[target_ref$UniProtKBID %in% Proteins]$Gene_name,
                Proteins))
    
    # Create the boxplot of MS/MS count per raw files and ORFs, also differentiate
    # raw files according to whether they belong to the candidate novel ORF or not
    pl <- plots_box(
        data = data, key = "Name", value = "Count", fill = "Type",
        main = paste0(Target_id, " region spectral count"),
        textsize = 20, bw = TRUE,
        label = c("N_label", "M_label", "Md_label"))
    plot(pl[[1]])
    
    # Close pdf output
    dev.off()
    
}



### Nucleotide sequence analysis -----------------------------------------

# Create new GRange object to study nucleotide frequence of start codon
start_grange <- subset(ref_grange, strand == "+")
names(start_grange) %<>%
    sub("^", "start_", .)
end(start_grange) <- (start(start_grange) + 2)

# Get the nucleotide sequence associated with the start codon
tmp <- getSeq(x = bsu, names = start_grange)

# After blasting of reference protein against genome, the genomic coordinates
# were corrected and the weird results previously observed here are no longer
# seen, indeed B. subtilis start codon should be as follows:
# ATG, TTG and GTG start codons are used in 78%, 13% and 9% of CDSs
nucleotideFrequencyAt(
    x = tmp, at = 1, as.prob = FALSE, as.array = TRUE,
    fast.moving.side = "right", with.labels = TRUE)
nucleotideFrequencyAt(
    x = tmp, at = 2, as.prob = FALSE, as.array = TRUE,
    fast.moving.side = "right", with.labels = TRUE)
nucleotideFrequencyAt(
    x = tmp, at = 3, as.prob = FALSE, as.array = TRUE,
    fast.moving.side = "right", with.labels = TRUE)
trinucleotideFrequency(
    x = tmp, step = 1, as.prob = FALSE, as.array = FALSE,
    fast.moving.side = "right", with.labels = TRUE,
    simplify.as = "collapsed")

# Compute the MS/MS count for all known protein and ORF within range of
# the target and for each raw file
data <- evid_match %>%
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

# Get the most and lowest expressed protein
low_express <- data %>%
    dplyr::slice(., 1:1525)
high_express <- data %>%
    dplyr::slice(., c((nrow(.)-1524):nrow(.)))

# Create GRange object to study RBS nucleotide frequence (+/- 50bp of start)
rbs_low_grange <- subset(
    ref_grange, strand == "+" & UniProtKBID %in% low_express$Proteins)
ranges(rbs_low_grange) <- IRanges(
    start = (start(rbs_low_grange) - 25),
    end = (start(rbs_low_grange) + 2))
names(rbs_low_grange) <- sub("^", "rbs_", rbs_low_grange$UniProtKBID)

# Get the nucleotide sequence associated with the start codon
rbs_low_seq <- getSeq(x = bsu, names = rbs_low_grange)
rbs_low_seq_list <- as.vector(rbs_low_seq)

# Draw consensus logo sequence of RBS for most highly expressed protein
ggplot() +
    geom_logo(
        data = rbs_low_seq_list,
        method = "bits",
        seq_type = "dna") +
    theme_logo()

# Export sequences to fasta file
Biostrings::writeXStringSet(
    x = rbs_low_seq,
    filepath = paste0(opt$out_path, "/", date_str, "_low_expr_RBS_seq.fasta"),
    append = FALSE, compress = FALSE,
    compression_level = NA, format = "fasta")

# Create GRange object to study RBS nucleotide frequence (+/- 50bp of start)
rbs_high_grange <- subset(
    ref_grange, strand == "+" & UniProtKBID %in% high_express$Proteins)
ranges(rbs_high_grange) <- IRanges(
    start = (start(rbs_high_grange) - 25),
    end = (start(rbs_high_grange) + 2))
names(rbs_high_grange) <- sub("^", "rbs_", rbs_high_grange$UniProtKBID)

# Get the nucleotide sequence associated with the start codon
rbs_high_seq <- getSeq(x = bsu, names = rbs_high_grange)
rbs_high_seq_list <- as.vector(rbs_high_seq)

# Draw consensus logo sequence of RBS for most highly expressed protein
ggplot() +
    geom_logo(
        data = rbs_high_seq_list,
        method = "bits",
        seq_type = "dna") +
    theme_logo()

# Export sequences to fasta file
Biostrings::writeXStringSet(
    x = rbs_high_seq,
    filepath = paste0(opt$out_path, "/", date_str, "_high_expr_RBS_seq.fasta"),
    append = FALSE, compress = FALSE,
    compression_level = NA, format = "fasta")




# Try to define the RBS motif sequence
#rbs_low_motifs <- oligonucleotideFrequency(
#    x = rbs_low_seq, width = 8, step = 1, as.prob = FALSE,
#    as.array = FALSE, fast.moving.side = "right", with.labels = TRUE,
#    simplify.as = "collapsed")


load_package("motifRG")

bg_low_grange <- subset(
    ref_grange, strand == "+" & UniProtKBID %in% low_express$Proteins)
ranges(bg_low_grange) <- IRanges(
    start = (end(bg_low_grange) - 25),
    end = (end(bg_low_grange) + 2))
names(bg_low_grange) <- sub("^", "bg_", bg_low_grange$UniProtKBID)

bg_low_seq <- getSeq(x = bsu, names = bg_low_grange)


rbs_low_motifs <- findMotifFgBg(
    fg.seq = rbs_low_seq, bg.seq = bg_low_seq,
    start.width = 6, both.strand = FALSE, flank = 1,
    enriched.only = TRUE, max.width = 16)



bg_high_grange <- subset(
    ref_grange, strand == "+" & UniProtKBID %in% high_express$Proteins)
ranges(bg_high_grange) <- IRanges(
    start = (end(bg_high_grange) - 25),
    end = (end(bg_high_grange) + 2))
names(bg_high_grange) <- sub("^", "bg_", bg_high_grange$UniProtKBID)

bg_high_seq <- getSeq(x = bsu, names = bg_high_grange)


rbs_high_motifs <- findMotifFgBg(
    fg.seq = rbs_high_seq, bg.seq = bg_high_seq,
    start.width = 6, both.strand = FALSE, flank = 1,
    enriched.only = TRUE, max.width = 16)


