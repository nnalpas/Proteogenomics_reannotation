#!/usr/bin/env Rscript

# This script determines the novelty reasons for each peptide and ORF based on
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
        output = readline(
            prompt = "Define the output directory!"))
    
    
    
    
        peptide_location = ,
        
        novel_fasta = choose.files(
            caption = "Choose Fasta file containing ORF proteins!",
            multi = FALSE),
        blast_ref = ,
        reciprocal_blast_uniprot = ,
        reciprocal_blast_ncbi = ,
        output = )
    
} else {
    
    # Define the list of command line parameters
    option_list <- list(
        make_option(
            opt_str = c("-e", "--evidence"),
            type = "character", default = NULL,
            help = "Evidence with database group info",
            metavar = "character"),
        make_option(
            opt_str = c("-r", "--reference_fasta"),
            type = "character", default = NULL,
            help = "Reference protein fasta file",
            metavar = "character"),
        make_option(
            opt_str = c("-o", "--output"),
            type = "character", default = "",
            help = "Output directory", metavar = "character"))
    
    
    
    
    
    
        make_option(
            opt_str = c("-p", "--peptide_location"),
            type = "character", default = NULL,
            help = "Peptide location file",
            metavar = "character"),
        
        make_option(
            opt_str = c("-n", "--novel_fasta"),
            type = "character", default = NULL,
            help = "ORF fasta file",
            metavar = "character"),
        make_option(
            opt_str = c("-b", "--blast_ref"),
            type = "character", default = NULL,
            help = "Blast against reference",
            metavar = "character"),
        make_option(
            opt_str = c("-ru", "--reciprocal_blast_uniprot"),
            type = "character", default = NULL,
            help = "Reciprocal blast against all uniprot",
            metavar = "character"),
        make_option(
            opt_str = c("-rn", "--reciprocal_blast_ncbi"),
            type = "character", default = NULL,
            help = "Reciprocal blast against all uniprot",
            metavar = "character"),
        make_option(
            opt_str = c("-o", "--output"),
            type = "character", default = "", 
            help = "Output file name [default= %default]",
            metavar = "character"))
    
    # Parse the parameters provided on command line by user
    opt_parser <- OptionParser(option_list = option_list)
    opt <- parse_args(opt_parser)
    
}

# Check whether inputs parameter was provided
if (is.null(opt$evidence)) {
    
    print_help(opt_parser)
    stop("The input evidence must be supplied!")
    
}
if (is.null(opt$reference_fasta)) {
    
    print_help(opt_parser)
    stop("The input reference fasta must be supplied!")
    
}

# Check whether output parameter was provided
if (opt$output == "") {
    
    opt$output <- dirname(opt$evidence)
    warning(paste("Output results to ", opt$output, "!"))
    
}






if (is.null(opt$fasta) | is.null(opt$blast) | is.null(opt$genomic_position)) {
    
    print_help(opt_parser)
    stop(paste(
        "The  input arguments must be supplied",
        "(fasta, blast and genomic position files)!"))
    
}



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







# Import the peptide location within proteins
pep_loc <- opt$peptide_location %>%
    as.character(.) %>%
    readRDS(file = .)

# Import the reciprocal blast best hits between novel proteins and reference
# proteins
reciprocal_blast_ref <- opt$reciprocal_blast_ref %>%
    as.character(.) %>%
    read.table(file = ., header = TRUE, sep = "\t", quote = "")

# Import the reciprocal blast best hits between novel proteins and all
# uniprot proteins
reciprocal_blast_uniprot <- opt$reciprocal_blast_uniprot %>%
    as.character(.) %>%
    read.table(file = ., header = TRUE, sep = "\t", quote = "")



ggplot(data = evid, mapping = aes(x = PEP, fill = Database)) +
    geom_density(aes(y = ..scaled..), alpha = 0.5) +
    geom_vline(
        xintercept = median(evid[evid$Database == "Target", "PEP"]),
        colour = "red", size = 1, linetype = "dashed") +
    coord_cartesian(xlim = c(0, quantile(evid$PEP, probs = 0.95))) +
    ylab(label = "Scaled density") +
    theme_bw()



### Levenshtein distance -------------------------------------------------

# Keep only sequence for novel peptide
novel_pep <- evid %>%
    dplyr::filter(., group == "Novel") %>%
    .[["Sequence"]] %>%
    unique(.)

# Compute the levenshtein distance for all novel peptide
leven_data <- adist(
    x = novel_pep,
    y = fasta$Known,
    partial = TRUE,
    ignore.case = TRUE)

# Reformat the data (transpose dataframe and gather)
leven_data %<>%
    set_rownames(novel_pep) %>%
    t(.) %>%
    base::data.frame(
        id = rownames(.),
        .,
        stringsAsFactors = FALSE) %>%
    tidyr::gather(data = ., key = "Sequence", value = "leven", -id)

# Keep the minimum levenshtein score result per peptide
leven_data %<>%
    dplyr::group_by(., Sequence) %>%
    dplyr::filter(., leven == min(leven)) %>%
    base::as.data.frame(., stringsAsFActors = FALSE)



### Reciprocal blast analysis --------------------------------------------

# Filter out hits that have e-value above 0.0001
reciprocal_blast_ref_filt <- reciprocal_blast_ref %>%
    dplyr::filter(., evalue_blast < 0.0001 & evalue_reciproc < 0.0001)
reciprocal_blast_uniprot_filt <- reciprocal_blast_uniprot %>%
    dplyr::filter(., evalue_blast < 0.0001 & evalue_reciproc < 0.0001)



### Novelty reason determination -----------------------------------------

# Get the list of novel peptides
novel_pep <- evid %>%
    dplyr::filter(., group == "Novel") %>%
    .[["Sequence"]] %>%
    unique(.)

# Determine the novelty reasons for each peptide
novelty_reasons <- purrr::map(
    .x = novel_pep,
    .f = novel_pep_classify,
    coordinate = pep_loc$Novel,
    levenshtein = leven_data,
    blast_ref = reciprocal_blast_ref_filt,
    blast_all = reciprocal_blast_uniprot_filt) %>%
    set_names(novel_pep) %>%
    plyr::ldply(., "data.frame") %>%
    set_colnames(c("id", "reason"))










# Compile a condensed dataframe of novel peptide info
data <- evid_match %>%
    dplyr::filter(., group == "Novel") %>%
    dplyr::select(
        .,
        Sequence, Length, Proteins, Mass, `Mass Error [ppm]`,
        `Mass Error [Da]`, PEP, Score, Intensity) %>%
    dplyr::group_by(., Sequence) %>%
    dplyr::summarise(
        .,
        Length = toString(x = unique(Length), width = NULL),
        Proteins = toString(x = unique(Proteins), width = NULL),
        Mass = toString(x = unique(Mass), width = NULL),
        minMassErrorppm = min(`Mass Error [ppm]`, na.rm = TRUE),
        minMassErrorDa = min(`Mass Error [Da]`, na.rm = TRUE),
        minPEP = min(PEP, na.rm = TRUE),
        maxScore = max(Score, na.rm = TRUE),
        maxIntensity = max(Intensity, na.rm = TRUE)) %>%
    base::as.data.frame(., stringsAsFactors = FALSE)

# Add the novel peptide info from evidence to position info
pep.pos.final <- dplyr::left_join(
    x = pep.pos, y = data, by = "Sequence")







