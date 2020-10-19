#!/usr/bin/env Rscript

# This script determines the novelty reasons for each peptide based on
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
library(gtable)
library(ggpubr)
library(grid)
library(gridExtra)
library(purrr)
library(foreach)
library(doParallel)
library(taxize)



### Parameters setting up ------------------------------------------------

# Define input parameters (interactively or from command line)
if (interactive()) {
    
    # Define the list of input parameters
    opt <- list(
        evidence = choose.files(
            caption = "Choose evidence (.RDS) file containing database info!",
            multi = FALSE),
        reference_fasta = choose.files(
            caption = "Choose Fasta file containing known proteins!",
            multi = FALSE),
        reciprocal_blast_ref = choose.files(
            caption = paste(
                "Choose reciprocal best blast file of",
                "ORF versus Reference protein!"),
            multi = FALSE),
        reciprocal_blast_uniprot = choose.dir(
            caption = paste(
                "Choose reciprocal best blast folder!"),
            multi = FALSE),
        peptide_location = choose.files(
            caption = "Choose peptide position (.RDS) file!",
            multi = FALSE),
        threads = readline(prompt = "How many cores to use?") %>% as.integer(),
        output = readline(
            prompt = "Define the output directory!"))
    
} else {
    
    # Define the list of command line parameters
    option_list <- list(
        make_option(
            opt_str = c("-e", "--evidence"),
            type = "character", default = NULL,
            help = "Evidence with database group info",
            metavar = "character"),
        make_option(
            opt_str = c("-f", "--reference_fasta"),
            type = "character", default = NULL,
            help = "Reference protein fasta file",
            metavar = "character"),
        make_option(
            opt_str = c("-r", "--reciprocal_blast_ref"),
            type = "character", default = NULL,
            help = "Reciprocal best blast against reference",
            metavar = "character"),
        make_option(
            opt_str = c("-u", "--reciprocal_blast_uniprot"),
            type = "character", default = NULL,
            help = "Reciprocal best blast folder",
            metavar = "character"),
        make_option(
            opt_str = c("-p", "--peptide_location"),
            type = "character", default = NULL,
            help = "Peptide location file",
            metavar = "character"),
        make_option(
            opt_str = c("-t", "--threads"),
            type = "integer", default = "",
            help = "Number of cores to use",
            metavar = "character"),
        make_option(
            opt_str = c("-o", "--output"),
            type = "character", default = "",
            help = "Output directory", metavar = "character"))
    
    # Parse the parameters provided on command line by user
    opt_parser <- OptionParser(option_list = option_list)
    opt <- parse_args(opt_parser)
    
}

# Check whether inputs parameter was provided
if (
    identical(opt$evidence, NULL) |
    identical(opt$evidence, "") |
    identical(opt$evidence, character(0))) {
    
    print_help(opt_parser)
    stop("The input evidence must be supplied!")
    
}
if (
    identical(opt$reference_fasta, NULL) |
    identical(opt$reference_fasta, "") |
    identical(opt$reference_fasta, character(0))) {
    
    print_help(opt_parser)
    stop("The input reference fasta must be supplied!")
    
}
if (
    identical(opt$reciprocal_blast_ref, NULL) |
    identical(opt$reciprocal_blast_ref, "") |
    identical(opt$reciprocal_blast_ref, character(0))) {
    
    print_help(opt_parser)
    stop("The input reciprocal blast against reference must be supplied!")
    
}
if (
    identical(opt$reciprocal_blast_uniprot, NULL) |
    identical(opt$reciprocal_blast_uniprot, "") |
    identical(opt$reciprocal_blast_uniprot, character(0))) {
    
    print_help(opt_parser)
    stop("The input reciprocal blast folder must be supplied!")
    
}
if (
    identical(opt$peptide_location, NULL) |
    identical(opt$peptide_location, "") |
    identical(opt$peptide_location, character(0))) {
    
    print_help(opt_parser)
    stop("The input peptide position must be supplied!")
    
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
    
    opt$output <- dirname(opt$evidence)
    warning(paste("Output results to ", opt$output, "!"))
    
}

# Create output directory if not already existing
dir.create(opt$output)



### Data import ----------------------------------------------------------

# Import the evidence file containing peptide group info
evid <- opt$evidence %>%
    as.character(.) %>%
    readRDS(file = .)

# Import the fasta files
fasta <- c(Known = opt$reference_fasta) %>%
    purrr::map(
        .x = ., .f = seqinr::read.fasta, seqtype = "AA", as.string = TRUE)
names(fasta$Known) %<>%
    uni_id_clean(.)

# Import the reciprocal blast best hits between ORFs and
# reference proteins
reciprocal_blast_ref <- opt$reciprocal_blast_ref %>%
    as.character(.) %>%
    read.table(
        file = ., header = TRUE, sep = "\t", quote = "",
        as.is = TRUE, comment.char = "") %>%
    dplyr::mutate(
        ., staxid_blast = as.integer(staxid_blast),
        staxid_reciproc = as.integer(staxid_reciproc))

# Import the reciprocal blast best hits between novel ORFs and
# all uniprot proteins
reciprocal_blast_files <- opt$reciprocal_blast_uniprot %>%
    as.character(.) %>%
    list.files(
        path = .,
        pattern = "Best_Reciproc_Blast_.+_recip_annot",
        full.names = TRUE) %>%
    set_names(sub(".+vs_(.+)_recip_annot", "\\1", .)) %>%
    grep("cross", ., invert = TRUE, value = TRUE)
reciprocal_blast_all <- lapply(X = reciprocal_blast_files, FUN = function(x) {
        read.table(
            file = x, header = TRUE, sep = "\t", quote = "",
            as.is = TRUE, comment.char = "") %>%
        dplyr::mutate(
            ., staxid_blast = as.integer(staxid_blast),
            staxid_reciproc = as.integer(staxid_reciproc))
    })

# Import the peptide location within proteins
pep_loc <- opt$peptide_location %>%
    as.character(.) %>%
    readRDS(file = .)

# Keep only sequence for novel peptide
novel_pep <- evid %>%
    dplyr::filter(., group == "Novel") %>%
    .[["Sequence"]] %>%
    unique(.)



### Novelty reason determination -----------------------------------------

# Compile all reciprocal best hits
reciprocal_blast_all %<>%
    plyr::ldply(., "data.frame", .id = "DB", stringsAsFactors = F)

# Identify missing taxon name
all_tax_ids <- unique(
    c(reciprocal_blast_all$TaxonID, reciprocal_blast_ref$TaxonID)) %>%
    split(x = ., f = ceiling(seq_along(.) / 3))

# Query NCBI for taxon name, without API key only 3 queries per seconds
# are allowed, therefor include a pause of 2 seconds between query
my_taxon <- lapply(X = 1:length(all_tax_ids), FUN = function(x) {
    Sys.sleep(time = 2)
    taxize::id2name(all_tax_ids[[x]], db = "ncbi")}) %>%
    unlist(., recursive = FALSE)
my_taxon %<>%
    plyr::ldply(., "data.frame", .id = NULL) %>%
    dplyr::mutate(., TaxonID = as.integer(id)) %>%
    dplyr::select(., TaxonID, Taxon = name)

# Merge the taxon name with the reciprocal blast results
reciprocal_blast_ref$Taxon <- NULL
reciprocal_blast_ref %<>%
    dplyr::left_join(x = ., y = my_taxon)
reciprocal_blast_all$Taxon <- NULL
reciprocal_blast_all %<>%
    dplyr::left_join(x = ., y = my_taxon)

# Determine the novelty reasons for each peptide
novelty_reasons <- purrr::map(
    .x = novel_pep,
    .f = novel_pep_classify,
    coordinate = pep_loc,
    blast_ref = reciprocal_blast_ref,
    blast_all = reciprocal_blast_all) %>%
    plyr::ldply(., "data.frame", .id = NULL, stringsAsFactors = FALSE) %>%
    unique(.)

# Get the leading proteins
evid_leading <- evid %>%
    dplyr::filter(., Sequence %in% novel_pep) %>%
    dplyr::select(., Sequence, `Leading proteins`) %>%
    unique(.)

# Use the leading protein to decide between multi ORFs mapping
novelty_reasons_format <- novelty_reasons %>%
    dplyr::left_join(x = ., y = evid_leading, by = "Sequence") %>%
    dplyr::filter(., warning == "" | Proteins == `Leading proteins`)

# Include the novelty reason column to the evidence table
evid_reason <- novelty_reasons_format %>%
    dplyr::select(., -Proteins, -`Leading proteins`) %>%
    dplyr::left_join(x = evid, y = ., by = "Sequence")

# Update the evidence table based on novelty explanation, any "Not novel"
# entry must have its database field set back to "Known"
warning(paste(
    "Exactly", length(which(evid_reason$NoveltyReason == "Not novel")),
    "evidence entries were reclassified as 'Not novel'!"))
warning(paste(
    "Exactly", length(which(is.na(evid_reason$NoveltyReason) & evid_reason$group == "Novel")),
    "evidence entries were reclassified as 'Reverse'!"))
evid_reason %<>%
    dplyr::mutate(
        .,
        group = dplyr::case_when(
            is.na(group) ~ "Known",
            (!is.na(NoveltyReason) & NoveltyReason == "Not novel") ~ "Known",
            (is.na(NoveltyReason) & group == "Novel") ~ "Reverse",
            TRUE ~ group),
        Database = dplyr::case_when(
            is.na(Database) ~ "Target",
            (!is.na(NoveltyReason) & NoveltyReason == "Not novel") ~ "Target",
            (is.na(NoveltyReason) & Database == "Novel") ~ "Decoy",
            TRUE ~ Database))



### Focus on high quality novel peptide ----------------------------------

# Define value for PEP filterig of the novel evidences
pep_class1 <- median(x = evid[evid$Database == "Target", "PEP"], na.rm = TRUE)
pep_class2 <- median(x = evid[evid$Database == "Novel", "PEP"], na.rm = TRUE)

# Add a column for the soft PEP filtering (based on median known evidence PEP)
evid_reason %<>%
    dplyr::mutate(
        ., PEPfilter = dplyr::case_when(
            PEP <= pep_class1 ~ "class 1",
            PEP > pep_class1 & PEP <= pep_class2 ~ "class 2",
            TRUE ~ "class 3"))



### Data visualisation and export ----------------------------------------

# Save the evidence reasons data
saveRDS(
    object = evid_reason,
    file = paste(
        opt$output, "/", "Sequence_novelty_reason.RDS", sep = ""))

# Export complete evidence info for these novel evidence
write.table(
    x = evid_reason,
    file = paste(
        opt$output, "/", "Sequence_novelty_reason.txt", sep = ""),
    quote = FALSE,
    sep = "\t",
    row.names = FALSE,
    col.names = TRUE)

# Open a file for plot visualisation
pdf(
    file = paste0(opt$output, "/", "Peptide_novelty_reason.pdf"),
    width = 12, height = 9)

# Define size of text for all plots
textsize <- 20

# Visualise the PEP for evidences between databases as density
evid_reason$Database <- factor(
    x = evid_reason$Database,
    levels = c("Target", "Decoy", "Novel"),
    ordered = TRUE)

# Define colour code for database
database_colours <- c(
    `Target` = "#2180FC", `Decoy` = "#ABA9A9", `Novel` = "#F23D3D")

# 
toplot <- evid_reason %>%
    dplyr::filter(., PEP < 1)
toplot$Database <- factor(
    x = toplot$Database,
    levels = c("Decoy", "Novel", "Target"),
    ordered = TRUE)

pl <- ggplot(
    data = toplot,
    mapping = aes(x = PEP, fill = Database, colour = Database)) +
    geom_density(aes(y = ..scaled..), alpha = 0.8, size = 1) +
    geom_vline(
        xintercept = pep_class1,
        colour = "#4f6990", size = 1, linetype = "dashed") +
    annotation_custom(grob = grobTree(
        textGrob(
            paste("Class 1 <", pep_class1), x = pep_class1, y = 0.97,
            hjust = -0.3, gp = gpar(col = "#4f6990")))) +
    geom_vline(
        xintercept = pep_class2,
        colour = "#903333", size = 1, linetype = "dashed") +
    annotation_custom(grob = grobTree(
        textGrob(
            paste("Class 2 <", pep_class2), x = pep_class2, y = 0.92,
            hjust = -0.4, gp = gpar(col = "#903333")))) +
    coord_cartesian(
        xlim = c(0, quantile(evid$PEP, probs = 0.95)),
        ylim = c(0, 1.1)) +
    ylab(label = "Scaled density") +
    ggtitle(label = "PEP density (with PEP soft threshold)") +
    theme_bw() +
    theme(
        legend.position = "bottom",
        title = element_text(
            #face = "bold",
            size = (textsize * 1.2)),
        text = element_text(size = textsize),
        plot.title = element_text(
            #face = "bold",
            size = (textsize * 1.5))) +
    scale_fill_manual(values = database_colours) +
    scale_colour_manual(guide = FALSE, values = database_colours)
pl

# Split and gather the evidence data so that it can be analysed separately
# per feature (removing redundancy)
evid_reason_long <- evid_reason %>%
    dplyr::filter(., Database == "Novel") %>%
    cSplit(
        indt = ., splitCols = "Proteins", sep = ";",
        direction = "long", fixed = TRUE) %>%
    dplyr::select(
        ., evid = id, pep = Sequence, orf = Proteins,
        PEPfilter, OnlyIdBySite, NoveltyReason) %>%
    tidyr::gather(
        data = ., key = "feature", value = "id",
        evid, pep, orf, convert = TRUE) %>%
    dplyr::mutate(
        ., PEPfilter = ifelse(!OnlyIdBySite, "class 4", PEPfilter)) %>%
    unique(.)

# Compute novel feature count (peptide and ORF) in total
data <- evid_reason_long %>%
    dplyr::group_by(., feature) %>%
    dplyr::summarise(
        .,
        PEPfilter = "All novel",
        Count = n_distinct(id)) %>%
    base::as.data.frame(., stringsAsFactors = FALSE)

# Compute novel feature count (peptide and ORF) accross PEP classes
data <- evid_reason_long %>%
    dplyr::group_by(., feature, id) %>%
    dplyr::summarise(., PEPfilter = dplyr::case_when(
        any(PEPfilter == "class 1") ~ "class 1",
        any(PEPfilter == "class 2") ~ "class 2",
        any(PEPfilter == "class 3") ~ "class 3",
        any(PEPfilter == "class 4") ~ "class 4",
        TRUE ~ NA_character_)) %>%
    dplyr::ungroup(.) %>%
    dplyr::group_by(., feature, PEPfilter) %>%
    dplyr::summarise(
        ., Count = n_distinct(id)) %>%
    base::as.data.frame(., stringsAsFactors = FALSE) %>%
    dplyr::bind_rows(data, .)

# Feature as factor to keep order in plot
data$feature <- factor(
    x = data$feature,
    levels = c("evid", "pep", "orf"),
    labels = c("Evidence", "Peptide", "ORF"),
    ordered = TRUE)

# Plot the novel feature count
pl <- ggplot(
    data, aes(
        x = feature, y = Count, fill = PEPfilter,
        colour = PEPfilter, label = Count)) +
    geom_bar(stat = "identity", position = "dodge", alpha = 0.7, size = 1) +
    geom_text(
        position = position_dodge(width = 0.9), vjust = -0.3, hjust = 0.5) +
    ggtitle("Novel peptide/ORF count") +
    xlab(NULL) +
    theme_pubr() +
    theme(
        legend.position = "bottom",
        title = element_text(
            #face = "bold",
            size = (textsize * 1.2)),
        text = element_text(size = textsize),
        plot.title = element_text(
            #face = "bold",
            size = (textsize * 1.5)))
pl

# Compute novel feature count (peptide and evidence) in total
data <- evid_reason_long %>%
    dplyr::filter(., feature != "orf") %>%
    dplyr::group_by(., feature, NoveltyReason) %>%
    dplyr::summarise(
        .,
        PEPfilter = "All novel",
        Count = n_distinct(id)) %>%
    base::as.data.frame(., stringsAsFactors = FALSE)

# Compute novel feature count (peptide and evidence) accross PEP classes
data <- evid_reason_long %>%
    dplyr::filter(., feature != "orf") %>%
    dplyr::group_by(., feature, id, NoveltyReason) %>%
    dplyr::summarise(., PEPfilter = dplyr::case_when(
        any(PEPfilter == "class 1") ~ "class 1",
        any(PEPfilter == "class 2") ~ "class 2",
        any(PEPfilter == "class 3") ~ "class 3",
        any(PEPfilter == "class 4") ~ "class 4",
        TRUE ~ NA_character_)) %>%
    dplyr::ungroup(.) %>%
    dplyr::group_by(., feature, PEPfilter, NoveltyReason) %>%
    dplyr::summarise(
        ., Count = n_distinct(id)) %>%
    base::as.data.frame(., stringsAsFactors = FALSE) %>%
    dplyr::bind_rows(data, .)

# Feature as factor to keep order in plot
data$feature <- factor(
    x = data$feature,
    levels = c("evid", "pep", "orf"),
    labels = c("Evidence", "Peptide", "ORF"),
    ordered = TRUE)

# Plot the evidence count per novelty reasons
pl <- ggplot(
    data %>% dplyr::filter(., feature == "Evidence"),
    aes(
        x = NoveltyReason, y = Count, fill = PEPfilter,
        colour = PEPfilter, label = Count)) +
    geom_bar(stat = "identity", position = "dodge", alpha = 0.7, size = 1) +
    geom_text(
        position = position_dodge(width = 0.9), vjust = 0.5, hjust = -0.5) +
    ggtitle("Evidence novelty reasons") +
    xlab("Novelty reason") +
    theme_pubr() +
    theme(
        legend.position = "bottom",
        title = element_text(
            #face = "bold",
            size = (textsize * 1.2)),
        text = element_text(size = textsize),
        plot.title = element_text(
            #face = "bold",
            size = (textsize * 1.5))) +
    coord_flip() +
    facet_grid(~PEPfilter)
pl

# Plot the peptide count per novelty reasons
pl <- ggplot(
    data %>% dplyr::filter(., feature == "Peptide"),
    aes(
        x = NoveltyReason, y = Count, fill = PEPfilter,
        colour = PEPfilter, label = Count)) +
    geom_bar(stat = "identity", position = "dodge", alpha = 0.7, size = 1) +
    geom_text(
        position = position_dodge(width = 0.9), vjust = 0.5, hjust = -0.5) +
    ggtitle("Peptide novelty reasons") +
    xlab("Novelty reason") +
    theme_pubr() +
    theme(
        legend.position = "bottom",
        title = element_text(
            #face = "bold",
            size = (textsize * 1.2)),
        text = element_text(size = textsize),
        plot.title = element_text(
            #face = "bold",
            size = (textsize * 1.5))) +
    coord_flip() +
    facet_grid(~PEPfilter)
pl

# Close the picture file
dev.off()



### END ------------------------------------------------------------------

# Close the cluster
stopImplicitCluster()

# Time scripts end
print(paste("Completed ", format(Sys.time(), "%Y-%m-%d"), sep = ""))


