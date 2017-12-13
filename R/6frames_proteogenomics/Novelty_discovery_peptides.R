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
            opt_str = c("-r", "--reference_fasta"),
            type = "character", default = NULL,
            help = "Reference protein fasta file",
            metavar = "character"),
        make_option(
            opt_str = c("-rr", "--reciprocal_blast_ref"),
            type = "character", default = NULL,
            help = "Reciprocal best blast against reference",
            metavar = "character"),
        make_option(
            opt_str = c("-ru", "--reciprocal_blast_uniprot"),
            type = "character", default = NULL,
            help = "Reciprocal best blast against all uniprot",
            metavar = "character"),
        make_option(
            opt_str = c("-rn", "--reciprocal_blast_ncbi"),
            type = "character", default = NULL,
            help = "Reciprocal best blast against all uniprot",
            metavar = "character"),
        make_option(
            opt_str = c("-p", "--peptide_location"),
            type = "character", default = NULL,
            help = "Peptide location file",
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
if (is.null(opt$evidence)) {
    
    print_help(opt_parser)
    stop("The input evidence must be supplied!")
    
}
if (is.null(opt$reference_fasta)) {
    
    print_help(opt_parser)
    stop("The input reference fasta must be supplied!")
    
}
if (is.null(opt$reciprocal_blast_ref)) {
    
    print_help(opt_parser)
    stop("The input reciprocal blast against reference must be supplied!")
    
}
if (is.null(opt$reciprocal_blast_uniprot)) {
    
    print_help(opt_parser)
    stop("The input reciprocal blast against UniProt must be supplied!")
    
}
if (is.null(opt$reciprocal_blast_ncbi)) {
    
    print_help(opt_parser)
    stop("The input reciprocal blast against NCBI must be supplied!")
    
}
if (is.null(opt$peptide_location)) {
    
    print_help(opt_parser)
    stop("The input peptide position must be supplied!")
    
}

# Check whether output parameter was provided
if (opt$output == "") {
    
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
        as.is = TRUE, comment.char = "")

# Import the reciprocal blast best hits between novel ORFs and
# all uniprot proteins
reciprocal_blast_uniprot <- opt$reciprocal_blast_uniprot %>%
    as.character(.) %>%
    read.table(
        file = ., header = TRUE, sep = "\t", quote = "",
        as.is = TRUE, comment.char = "")

# Import the reciprocal blast best hits between novel ORFs and
# all uniprot proteins
reciprocal_blast_ncbi <- opt$reciprocal_blast_ncbi %>%
    as.character(.) %>%
    read.table(
        file = ., header = TRUE, sep = "\t", quote = "",
        as.is = TRUE, comment.char = "")

# Import the peptide location within proteins
pep_loc <- opt$peptide_location %>%
    as.character(.) %>%
    readRDS(file = .)



### Levenshtein distance -------------------------------------------------

# Keep only sequence for novel peptide
novel_pep <- evid %>%
    dplyr::filter(., group == "Novel") %>%
    .[["Sequence"]] %>%
    unique(.)

# Compute the levenshtein distance for all novel peptide
leven_dist <- adist(
    x = novel_pep,
    y = fasta$Known,
    partial = TRUE,
    ignore.case = TRUE)

# Reformat the data (transpose dataframe and gather)
leven_dist_format <- leven_dist %>%
    set_rownames(novel_pep) %>%
    t(.) %>%
    base::data.frame(
        id = rownames(.),
        .,
        stringsAsFactors = FALSE) %>%
    tidyr::gather(data = ., key = "Sequence", value = "leven", -id)

# Keep the minimum levenshtein score result per peptide
#leven_dist_format %<>%
#    dplyr::group_by(., Sequence) %>%
#    dplyr::filter(., leven == min(leven)) %>%
#    base::as.data.frame(., stringsAsFActors = FALSE)



### Novelty reason determination -----------------------------------------

# Compile all reciprocal best hits against UniProt and NCBI
reciprocal_blast_all <- dplyr::bind_rows(
    reciprocal_blast_uniprot, reciprocal_blast_ncbi)

# Determine the novelty reasons for each peptide
novelty_reasons <- purrr::map(
    .x = novel_pep,
    .f = novel_pep_classify,
    coordinate = pep_loc$Novel,
    levenshtein = leven_dist_format,
    blast_ref = reciprocal_blast_ref,
    blast_all = reciprocal_blast_all) %>%
    set_names(novel_pep) %>%
    plyr::ldply(., "data.frame") %>%
    set_colnames(c("Sequence", "NoveltyReason"))

# Include the novelty reason column to the evidence table
evid_reason <- evid %>%
    dplyr::left_join(x = ., y = novelty_reasons, by = "Sequence")



### Focus on high quality novel peptide ----------------------------------

# Define value for PEP filterig of the novel evidences
pep_threshold <- median(evid[evid$Database == "Target", "PEP"])

# Add a column for the soft PEP filtering (based on median known evidence PEP)
evid_reason %<>%
    dplyr::mutate(., PEPfilter = ifelse(PEP <= pep_threshold, TRUE, FALSE))



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




ggplot(data = evid, mapping = aes(x = PEP, fill = Database)) +
    geom_density(aes(y = ..scaled..), alpha = 0.5) +
    geom_vline(
        xintercept = pep_threshold,
        colour = "red", size = 1, linetype = "dashed") +
    coord_cartesian(xlim = c(0, quantile(evid$PEP, probs = 0.95))) +
    ylab(label = "Scaled density") +
    theme_bw()


# Source the custom user functions
if (interactive()) {
    mark_report <- paste(
        "C:/Users",
        user,
        "Documents/GitHub/Proteogenomics_reannotation/",
        "R/6frames_proteogenomics/Bsu_proteogenomics_report.rmd",
        sep = "/")
} else {
    mark_report <- paste(
        "/home-link",
        user,
        "bin/Bsu_proteogenomics_report.rmd",
        sep = "/")
}


report_markdown(rmd_file = mark_report)


# Define the report markdown file
report.file <- paste(
    "C:/Users",
    user, 
    "Documents/GitHub/Miscellaneous/R/6frames_proteogenomics",
    "Bsu_proteogenomics_report.rmd", sep = "/")

# Define temporary location for report generation
tempReport <- file.path(tempdir(), basename(report.file))
file.copy(
    from = report.file,
    to = tempReport,
    overwrite = TRUE)

# Define the required variables as markdown parameters
param <- list(
    evidences = evid_match,
    pept.posit = pep.pos,
    levenshtein = leven.data,
    pheno = exp.design,
    bsu.ideo = bsu.ideo,
    ref.grange = ref.grange,
    novel.grange = orf.grange,
    pep_loc = pep_located,
    manual.annot = manual.annot) 

# Define output file name
out_file <- paste(
    work_space,
    "/Bsu_proteogenomics_",
    format(Sys.time(), '%Y%m%d_%H-%M'),
    ".html",
    sep = "")

# Render the markdown report
rmarkdown::render(
    input = tempReport,
    output_format = "ioslides_presentation",
    output_file = out_file,
    params = param,
    envir = new.env(parent = globalenv()))


