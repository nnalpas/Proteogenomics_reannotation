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
library(grid)
library(gridExtra)
library(purrr)
library(foreach)
library(doParallel)



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
        reciprocal_blast_uniprot = choose.files(
            caption = paste(
                "Choose reciprocal best blast file of",
                "ORF versus UniProt!"),
            multi = FALSE),
        reciprocal_blast_ncbi = choose.files(
            caption = paste(
                "Choose reciprocal best blast file of",
                "ORF versus NCBI!"),
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
            help = "Reciprocal best blast against all uniprot",
            metavar = "character"),
        make_option(
            opt_str = c("-n", "--reciprocal_blast_ncbi"),
            type = "character", default = NULL,
            help = "Reciprocal best blast against all uniprot",
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
    stop("The input reciprocal blast against UniProt must be supplied!")
    
}
if (
    identical(opt$reciprocal_blast_ncbi, NULL) |
    identical(opt$reciprocal_blast_ncbi, "") |
    identical(opt$reciprocal_blast_ncbi, character(0))) {
    
    print_help(opt_parser)
    stop("The input reciprocal blast against NCBI must be supplied!")
    
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
reciprocal_blast_uniprot <- opt$reciprocal_blast_uniprot %>%
    as.character(.) %>%
    read.table(
        file = ., header = TRUE, sep = "\t", quote = "",
        as.is = TRUE, comment.char = "") %>%
    dplyr::mutate(
        ., staxid_blast = as.integer(staxid_blast),
        staxid_reciproc = as.integer(staxid_reciproc))

# Import the reciprocal blast best hits between novel ORFs and
# all uniprot proteins
reciprocal_blast_ncbi <- opt$reciprocal_blast_ncbi %>%
    as.character(.) %>%
    read.table(
        file = ., header = TRUE, sep = "\t", quote = "",
        as.is = TRUE, comment.char = "") %>%
    dplyr::mutate(
        ., staxid_blast = as.integer(staxid_blast),
        staxid_reciproc = as.integer(staxid_reciproc))

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
reciprocal_blast_all <- dplyr::bind_rows(
    data.frame(DB = "Reference", reciprocal_blast_ref, stringsAsFactors = F),
    data.frame(DB = "UniProt", reciprocal_blast_uniprot, stringsAsFactors = F),
    data.frame(DB = "NCBI", reciprocal_blast_ncbi, stringsAsFactors = F))

# Determine the novelty reasons for each peptide
novelty_reasons <- purrr::map(
    .x = novel_pep,
    .f = novel_pep_classify,
    coordinate = pep_loc,
    blast_ref = reciprocal_blast_ref,
    blast_all = reciprocal_blast_all) %>%
    plyr::ldply(., "data.frame", .id = NULL, stringsAsFactors = FALSE)

# Include the novelty reason column to the evidence table
evid_reason <- evid %>%
    dplyr::left_join(x = ., y = novelty_reasons, by = "Sequence")

# Update the evidence table based on novelty explanation, any "Not novel"
# entry must have its database field set back to "Known"
warning(paste(
    "Exactly", length(which(evid_reason$NoveltyReason == "Not novel")),
    "evidence entries were reclassified as 'Not novel'!"))
evid_reason %<>%
    dplyr::mutate(
        .,
        group = ifelse(
            !is.na(NoveltyReason) & NoveltyReason == "Not novel",
            "Known", group),
        Database = ifelse(
            !is.na(NoveltyReason) & NoveltyReason == "Not novel",
            "Target", Database))



### Focus on high quality novel peptide ----------------------------------

# Define value for PEP filterig of the novel evidences
pep_class1 <- median(evid[evid$Database == "Target", "PEP"])
pep_class2 <- median(evid[evid$Database == "Novel", "PEP"])

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

# Compute novel feature count (peptide and ORF) in total
data <- evid_reason %>%
    dplyr::filter(., Database == "Novel") %>%
    cSplit(
        indt = ., splitCols = "Proteins", sep = ";",
        direction = "long", fixed = TRUE) %>%
    dplyr::select(., id, Sequence, Proteins) %>%
    unique(.) %>%
    dplyr::summarise(
        .,
        PEPfilter = "All novel",
        Number_evidence = n_distinct(id),
        Number_peptide = n_distinct(Sequence),
        Number_ORF = n_distinct(Proteins)) %>%
    tidyr::gather(
        data = ., key = "feature", value = "Count", -PEPfilter) %>%
    base::as.data.frame(., stringsAsFactors = FALSE)

# Compute novel feature count (peptide and ORF) that pass PEP filter
data <- evid_reason %>%
    dplyr::filter(., Database == "Novel" & OnlyIdBySite) %>%
    cSplit(
        indt = ., splitCols = "Proteins", sep = ";",
        direction = "long", fixed = TRUE) %>%
    dplyr::select(., id, Sequence, Proteins, PEPfilter) %>%
    unique(.) %>%
    dplyr::group_by(., PEPfilter) %>%
    dplyr::summarise(
        .,
        Number_evidence = n_distinct(id),
        Number_peptide = n_distinct(Sequence),
        Number_ORF = n_distinct(Proteins)) %>%
    tidyr::gather(
        data = ., key = "feature", value = "Count", -PEPfilter) %>%
    base::as.data.frame(., stringsAsFactors = FALSE) %>%
    dplyr::bind_rows(data, .)

# Compute novel feature count (peptide and ORF) for only identified by site
data <- evid_reason %>%
    dplyr::filter(., Database == "Novel"  & !OnlyIdBySite) %>%
    cSplit(
        indt = ., splitCols = "Proteins", sep = ";",
        direction = "long", fixed = TRUE) %>%
    dplyr::select(., id, Sequence, Proteins) %>%
    unique(.) %>%
    dplyr::summarise(
        .,
        PEPfilter = "class4 (OnlyIdBySite)",
        Number_evidence = n_distinct(id),
        Number_peptide = n_distinct(Sequence),
        Number_ORF = n_distinct(Proteins)) %>%
    tidyr::gather(
        data = ., key = "feature", value = "Count", -PEPfilter) %>%
    base::as.data.frame(., stringsAsFactors = FALSE) %>%
    dplyr::bind_rows(data, .)

# Feature as factor to keep order in plot
data$feature <- factor(
    x = data$feature,
    levels = c("Number_evidence", "Number_peptide", "Number_ORF"),
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

# Compute the peptide count per novelty reasons
data <- evid_reason %>%
    dplyr::filter(., Database == "Novel") %>%
    dplyr::select(., Sequence, Proteins, NoveltyReason) %>%
    dplyr::group_by(., NoveltyReason) %>%
    dplyr::summarise(
        .,
        PEPfilter = "All novel",
        Number_peptide = n_distinct(Sequence),
        Number_ORF = n_distinct(Proteins)) %>%
    tidyr::gather(
        data = ., key = "feature", value = "Count",
        -NoveltyReason, -PEPfilter) %>%
    base::as.data.frame(., stringsAsFactors = FALSE)

# Compute the peptide count per novelty reasons that pass PEP filter and
# pass only ID by sites
data <- evid_reason %>%
    dplyr::filter(., Database == "Novel" & OnlyIdBySite) %>%
    dplyr::select(., Sequence, Proteins, NoveltyReason, PEPfilter) %>%
    dplyr::group_by(., NoveltyReason, PEPfilter) %>%
    dplyr::summarise(
        .,
        Number_peptide = n_distinct(Sequence),
        Number_ORF = n_distinct(Proteins)) %>%
    tidyr::gather(
        data = ., key = "feature", value = "Count", -NoveltyReason, -PEPfilter) %>%
    base::as.data.frame(., stringsAsFactors = FALSE) %>%
    dplyr::bind_rows(data, .)

# Compute the peptide count per novelty reasons for only ID by sites
data <- evid_reason %>%
    dplyr::filter(., Database == "Novel" & !OnlyIdBySite) %>%
    dplyr::select(., Sequence, Proteins, NoveltyReason) %>%
    dplyr::group_by(., NoveltyReason) %>%
    dplyr::summarise(
        .,
        PEPfilter = "class4 (OnlyIdBySite)",
        Number_peptide = n_distinct(Sequence),
        Number_ORF = n_distinct(Proteins)) %>%
    tidyr::gather(
        data = ., key = "feature", value = "Count",
        -NoveltyReason, -PEPfilter) %>%
    base::as.data.frame(., stringsAsFactors = FALSE) %>%
    dplyr::bind_rows(data, .)

# Plot the peptide count per novelty reasons
pl <- ggplot(
    data %>% dplyr::filter(., feature == "Number_peptide"),
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

# Plot the ORF count per novelty reasons
pl <- ggplot(
    data %>% dplyr::filter(., feature == "Number_ORF"),
    aes(
        x = NoveltyReason, y = Count, fill = PEPfilter,
        colour = PEPfilter, label = Count)) +
    geom_bar(stat = "identity", position = "dodge", alpha = 0.7, size = 1) +
    geom_text(
        position = position_dodge(width = 0.9), vjust = 0.5, hjust = -0.5) +
    ggtitle("ORF novelty reasons") +
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


