#!/usr/bin/env Rscript

# This script determines the novel ORF that are validated through
# independant reprocessing in related organism



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
library(magrittr)
library(ggplot2)
library(msa)



### Parameters setting up ------------------------------------------------

# Define input parameters (interactively or from command line)
if (interactive()) {
    
    # Define the list of input parameters
    opt <- list(
        novel_reason = choose.files(
            caption = "Choose ORF novelty (.RDS) file!",
            multi = FALSE),
        pep_grange = choose.files(
            caption = "Choose peptide entries grange (.RDS) file!",
            multi = FALSE),
        reciprocal_blast_ref = choose.files(
            caption = paste(
                "Choose reciprocal best blast file of",
                "ORF versus Reference protein!"),
            multi = FALSE),
        reciprocal_blast_cross = choose.files(
            caption = paste(
                "Choose reciprocal best blast file of",
                "ORF versus related organism protein!"),
            multi = FALSE),
        evid_folder = choose.dir(
            caption = "Select folder with validation MaxQuant processing!"),
        fasta_ref = choose.files(
            caption = "Choose fasta file containing reference proteins!",
            multi = FALSE),
        fasta_novel = choose.files(
            caption = "Choose fasta file containing novel ORF!",
            multi = FALSE),
        fasta_cross = choose.files(
            caption = "Choose fasta file containing cross ORF!",
            multi = FALSE),
        pep_class = readline(prompt = paste0(
            "Which PEP class to keep",
            " (separated by comma; e.g. 'class 1')")) %>%
            as.character(),
        output = readline(
            prompt = "Define the output directory!"))
    
}

# Check whether inputs parameter was provided
if (
    identical(opt$novel_reason, NULL) |
    identical(opt$novel_reason, "") |
    identical(opt$novel_reason, character(0))) {
    
    print_help(opt_parser)
    stop("The input ORF novelty must be supplied!")
    
}
if (
    identical(opt$pep_grange, NULL) |
    identical(opt$pep_grange, "") |
    identical(opt$pep_grange, character(0))) {
    
    print_help(opt_parser)
    stop("The input peptide entries grange must be supplied!")
    
}
if (
    identical(opt$reciprocal_blast_ref, NULL) |
    identical(opt$reciprocal_blast_ref, "") |
    identical(opt$reciprocal_blast_ref, character(0))) {
    
    print_help(opt_parser)
    stop("The input reciprocal blast against reference must be supplied!")
    
}
if (
    identical(opt$reciprocal_blast_cross, NULL) |
    identical(opt$reciprocal_blast_cross, "") |
    identical(opt$reciprocal_blast_cross, character(0))) {
    
    print_help(opt_parser)
    stop("The input reciprocal blast against related organism must be supplied!")
    
}
if (
    identical(opt$evid_folder, NULL) |
    identical(opt$evid_folder, "") |
    identical(opt$evid_folder, character(0))) {
    
    print_help(opt_parser)
    stop("The input MaxQuant validation folder must be supplied!")
    
}
if (
    identical(opt$fasta_ref, NULL) |
    identical(opt$fasta_ref, "") |
    identical(opt$fasta_ref, character(0))) {
    
    print_help(opt_parser)
    stop("The input reference proteins fasta file must be supplied!")
    
}
if (
    identical(opt$fasta_novel, NULL) |
    identical(opt$fasta_novel, "") |
    identical(opt$fasta_novel, character(0))) {
    
    print_help(opt_parser)
    stop("The input novel ORF fasta file must be supplied!")
    
}
if (
    identical(opt$fasta_cross, NULL) |
    identical(opt$fasta_cross, "") |
    identical(opt$fasta_cross, character(0))) {
    
    print_help(opt_parser)
    stop("The input cross ORF fasta file must be supplied!")
    
}
if (
    identical(opt$pep_class, NULL) |
    identical(opt$pep_class, "") |
    identical(opt$pep_class, character(0))) {
    
    opt["pep_class"] <- list("class 1")
    warning("No PEP class provided by user, default to 'class1'!")
    
}

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

# Import the novelty reason file containing ORF novelty info
orf_reason <- opt$novel_reason %>%
    as.character(.) %>%
    readRDS(file = .)

# Import the peptides genomic ranges file
pep_grange <- opt$pep_grange %>%
    as.character(.) %>%
    readRDS(file = .)

# Import the reciprocal blast best hits between ORFs and
# reference proteins
reciprocal_blast_ref <- opt$reciprocal_blast_ref %>%
    strsplit(x = ., split = ",") %>%
    unlist(.) %>%
    as.character(.) %>%
    purrr::map(
        .x = ., .f = read.table, header = TRUE, sep = "\t", quote = "",
        as.is = TRUE, comment.char = "") %>%
    plyr::ldply(., "data.frame") %>%
    dplyr::mutate(
        ., staxid_blast = as.integer(staxid_blast),
        staxid_reciproc = as.integer(staxid_reciproc))

# Import the reciprocal blast best hits between ORFs and
# related organism proteins
reciprocal_blast_cross <- opt$reciprocal_blast_cross %>%
    strsplit(x = ., split = ",") %>%
    unlist(.) %>%
    as.character(.) %>%
    purrr::map(
        .x = ., .f = read.table, header = TRUE, sep = "\t", quote = "",
        as.is = TRUE, comment.char = "") %>%
    plyr::ldply(., "data.frame") %>%
    dplyr::mutate(
        ., staxid_blast = as.integer(staxid_blast),
        staxid_reciproc = as.integer(staxid_reciproc))

# Import the evidence.txt for all possible MaxQuant validation processing
evidence_files <- opt$evid_folder %>%
    as.character(.) %>%
    list.files(
        path = .,
        pattern = "PXD",
        recursive = FALSE,
        full.names = TRUE) %>%
    set_names(basename(.))
evidence_all <- lapply(X = evidence_files, FUN = function(x) {
    data.table::fread(
        input = paste0(x, "/combined/txt/evidence.txt"),
        sep = "\t", quote = "", header = TRUE,
        stringsAsFactors = FALSE, colClasses = "character")
}) %>%
    plyr::ldply(., dplyr::bind_rows, .id = "Processing")

# Import the fasta files
fasta <- c(
    Known = opt$fasta_ref, Novel = opt$fasta_novel,
    Cross = opt$fasta_cross) %>%
    purrr::map(
        .x = ., .f = seqinr::read.fasta, seqtype = "AA", as.string = TRUE)

# Format the PEP classes if defined
if (!is.null(opt$pep_class)) {
    opt$pep_class %<>%
        as.character(.) %>%
        strsplit(x = ., split = ",", fixed = TRUE) %>%
        unlist(.) %>%
        unique(.)
}



### Overlap novel peptides -----------------------------------------------

my_plots <- list()

toplot <- pep_grange %>%
    data.frame(.) %>%
    dplyr::select(., Sequence, group) %>%
    dplyr::mutate(., value = TRUE) %>%
    unique(.) %>%
    tidyr::pivot_wider(data = ., names_from = group, values_from = value)

toplot <- evidence_all %>%
    dplyr::select(., Sequence) %>%
    dplyr::mutate(., proteomeXchange = TRUE) %>%
    unique(.) %>%
    dplyr::full_join(x = toplot, y = .) %>%
    dplyr::mutate_all(~tidyr::replace_na(data = ., replace = FALSE))

my_plots[["venn_validation"]] <- ggvenn::ggvenn(
    data = toplot,
    columns = c("Known", "Novel", "proteomeXchange"),
    fill_color = c("#387eb8", "#d1d2d4", "#e21e25"))



### Localise peptides matching refined ORFs ------------------------------

pep_novel_uniq <- evidence_all %>%
    dplyr::select(., Processing, Sequence, Proteins) %>%
    tidyr::separate_rows(data = ., Proteins, sep = ";") %>%
    unique(.) %>%
    dplyr::filter(., Proteins %in% names(fasta$Cross))

my_pep_loc <- lapply(X = unique(pep_novel_uniq$Proteins), FUN = function(x) {
    curr_seqs <- unique(pep_novel_uniq[
        pep_novel_uniq$Proteins == x, ][["Sequence"]])
    my_pep_loc <- stringr::str_locate_all(
        string = fasta$Cross[[x]],
        pattern = curr_seqs) %>%
        set_names(curr_seqs) %>%
        plyr::ldply(., data.frame, .id = "Sequence")
}) %>%
    set_names(unique(pep_novel_uniq$Proteins)) %>%
    plyr::ldply(., data.frame, .id = "Proteins") %>%
    dplyr::left_join(x = pep_novel_uniq, y = .) %>%
    dplyr::mutate(., id = 1:nrow(.))



### Identify refined protein region --------------------------------------

reciprocal_blast_all <- reciprocal_blast_ref %>%
    dplyr::filter(., qseqid %in% orf_reason$Proteins) %>%
    dplyr::select(
        ., orfid = qseqid, refid = sseqid, orfstart = qstart_blast,
        orfend = qend_blast, orflen = qlen_blast,
        refstart = sstart_blast, refend = send_blast,
        reflen = slen_blast)

reciprocal_blast_all <- reciprocal_blast_cross %>%
    dplyr::select(
        ., orfid = qseqid, crossid = sseqid,
        orfcrossstart = qstart_blast, orfcrossend = qend_blast,
        crossstart = sstart_blast, crossend = send_blast,
        crosslen = slen_blast) %>%
    dplyr::full_join(x = reciprocal_blast_all, y = ., by = "orfid") %>%
    dplyr::left_join(
        x = ., y = orf_reason[, c("Proteins", "ORFNoveltyReason")],
        by = c("orfid" = "Proteins"))

reciprocal_blast_all %<>%
    dplyr::mutate(., refinestart = dplyr::case_when(
        grepl("alternate end", ORFNoveltyReason) & !is.na(crossid) & orfcrossend < orfend ~ NA_integer_,
        grepl("alternate end", ORFNoveltyReason) & !is.na(crossid) ~ as.integer((crossend - (orfcrossend - orfend))),
        TRUE ~ crossstart),
    refineend = dplyr::case_when(
        grepl("alternate start", ORFNoveltyReason) & !is.na(crossid) & orfcrossstart > orfstart ~ NA_integer_,
        grepl("alternate start", ORFNoveltyReason) & !is.na(crossid) ~ as.integer(crossstart + (orfstart - orfcrossstart)),
        TRUE ~ crossend
    )) %>%
    dplyr::filter(., !is.na(refinestart) & !is.na(refineend))



### Blast view of cross-validation ---------------------------------------

setwd(opt$output)

orf_reason$PX_datasets <- NA_character_
orf_reason$PX_novel_sequence <- NA_character_
orf_reason$PX_novel_peptide_count <- NA_integer_
my_pep_loc$Group <- "Known"
for (x in reciprocal_blast_all$crossid) {
    
    my_orf <- reciprocal_blast_all[
        reciprocal_blast_all$crossid == x, ][["orfid"]]
    my_ref <- reciprocal_blast_all[
        reciprocal_blast_all$crossid == x, ][["refid"]]
    my_start <- reciprocal_blast_all[
        reciprocal_blast_all$crossid == x, ][["refinestart"]]
    my_end <- reciprocal_blast_all[
        reciprocal_blast_all$crossid == x, ][["refineend"]]
    valid_pep <- my_pep_loc[
        my_pep_loc$Proteins == x &
            my_pep_loc$start < my_end &
            my_pep_loc$end > my_start, ]
    my_pep_loc[my_pep_loc$id %in% valid_pep$id, "Group"] <- "Novel"
    
    if (nrow(valid_pep) > 0) {
        
        orf_reason[
            orf_reason$Proteins == my_orf, "PX_datasets"] <- paste0(
                unique(valid_pep$Processing), collapse = ";")
        orf_reason[
            orf_reason$Proteins == my_orf, "PX_novel_sequence"] <- paste0(
            unique(valid_pep$Sequence), collapse = ";")
        orf_reason[
            orf_reason$Proteins == my_orf, "PX_novel_peptide_count"] <- length(
                unique(valid_pep[["Sequence"]]))
        
        seqs_string <- Biostrings::AAStringSet(x = c(
            fasta$Novel[[my_orf]],
            ifelse(is.na(my_ref), "FAKE", fasta$Known[[my_ref]]),
            fasta$Cross[[x]],
            as.character(unique(valid_pep[["Sequence"]]))
        ))
        
        my_alignment <- msa(inputSeqs = seqs_string, order = "input")
        
        msaPrettyPrint(
            x = my_alignment, output = "pdf", showNames = "none",
            file = paste0(opt$output, "/", my_orf, "_alignment.pdf"),
            showLogo = "none", askForOverwrite = FALSE, verbose = TRUE)
        
    }
    
}

orf_reason$`PX valid` <- ifelse(
    is.na(orf_reason$PX_novel_sequence), FALSE, TRUE)

orf_reason$`High quality` <- ifelse(
    orf_reason$PEPfilter %in% opt$pep_class & orf_reason$OnlyIdBySite == TRUE,
    TRUE, FALSE)

orf_reason$`Peptide 2+` <- ifelse(
    orf_reason$Novel_peptide_count > 1, TRUE, FALSE)

orf_reason$`Start valid` <- ifelse(
    !is.na(orf_reason$Starts), TRUE, FALSE)

data.table::fwrite(
    x = orf_reason,
    file = paste0(opt$output, "/ORF_validation.txt"),
    append = FALSE, quote = FALSE, sep = "\t",
    row.names = FALSE, col.names = TRUE)

saveRDS(object = orf_reason, file = paste0(opt$output, "/ORF_validation.RDS"))

data.table::fwrite(
    x = my_pep_loc,
    file = paste0(opt$output, "/Peptide_validation.txt"),
    append = FALSE, quote = FALSE, sep = "\t",
    row.names = FALSE, col.names = TRUE)

data.table::fwrite(
    x = reciprocal_blast_all,
    file = paste0(opt$output, "/Reciprocal_blast_validation.txt"),
    append = FALSE, quote = FALSE, sep = "\t",
    row.names = FALSE, col.names = TRUE)



### END ------------------------------------------------------------------

# Time scripts end
print(paste("Completed ", format(Sys.time(), "%Y-%m-%d"), sep = ""))


