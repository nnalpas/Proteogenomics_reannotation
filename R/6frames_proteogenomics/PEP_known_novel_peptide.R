


### Parameters setting up ------------------------------------------------

# Start with clean environment
rm(list = ls())



### Define working directory ---------------------------------------------

# Define the work space
work_space <- choose.dir()
setwd(work_space)

# Define the current user
user <- Sys.info()[["user"]]

# Define current time
date_str <- format(Sys.time(), "%Y-%m-%d")

# Variable for additional filter
filter <- FALSE



### List of required packages -----------------------------------------------

# Source the custom user functions
source(
    file = paste(
        "C:/Users/",
        user,
        "/Documents/GitHub/Miscellaneous/R/General/General_function.R",
        sep = ""))

# Load the required packages (or install if not already in library)
# Note: the select function from dplyr and UniProt.ws is overlapping,
# therefore order of package loading is important
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



### Data import ----------------------------------------------------------

# Select the maxquant txt folder
txt_dir <- choose.dir(caption = "Select the MaxQuant txt folder?")

# Import the maxquant evidence table
evid <- mq_read(
    path = txt_dir,
    name = "evidence.txt",
    integer64 = "double")

# Import the maxquant proteingroups table
pg <- mq_read(
    path = txt_dir,
    name = "proteinGroups.txt",
    integer64 = "double")
pg_not_sites <- pg %>%
    dplyr::filter(., `Only identified by site` == "") %>%
    cSplit(indt = ., splitCols = "Evidence IDs", sep = ";", direction = "long")

# List the fasta files that need to be imported
fasta_file <- c()
fasta_file["Known"] <- choose.files(
    caption = "Select KNOWN proteome Fasta file",
    multi = FALSE,
    filters = ".fasta")
fasta_file["Novel"] <- choose.files(
    caption = "Select NOVEL proteome Fasta file",
    multi = FALSE,
    filters = ".fasta")
fasta_file["Contaminant"] <- choose.files(
    caption = "Select CONTAMINANT proteome Fasta file",
    multi = FALSE,
    filters = ".fasta")

# Import the current fasta file
fasta_list <- purrr::map(
    .x = fasta_file, .f = read.fasta, seqtype = "AA", as.string = TRUE)

# Import the experimental design (with conditions)
exp.design <- read.table(
    file = "Bsu_Conditions_06022017.txt", header = TRUE,
    sep = "\t", quote = "", as.is = TRUE)

# Clean-up
rm(fasta_file)



### Novel peptide identification -----------------------------------------

# Additional filter
if (filter) {
    evid %<>%
        dplyr::filter(., id %in% unique(pg_not_sites$`Evidence IDs`))
}

# Format association of peptide to protein
data <- evid %>%
    cSplit(indt = ., splitCols = "Proteins", sep = ";", direction = "long") %>%
    dplyr::select(., Sequence, Proteins) %>%
    unique(.) %>%
    base::as.data.frame(., stringsAsFactors = FALSE)

# Format the fasta files
names(fasta_list$Contaminant) %<>%
    sub(pattern = "^", replacement = "CON__", x = .)
names(fasta_list$Known) %<>%
    sub(pattern = ".+\\|(.+)\\|.+", replacement = "\\1", x = .)
fastas <- c(fasta_list$Contaminant, fasta_list$Novel, fasta_list$Known)

# Locate all the peptide within associated proteins
pep_loc <- pept_locate(
    data = data,
    peptide = "Sequence",
    proteins = "Proteins",
    fasta = fastas)

# Find the source of the peptide from which fasta file they are originating
evid_match <- pep_loc %>%
    dplyr::group_by(., pep) %>%
    dplyr::summarise(
        .,
        group = ifelse(
            test = any(prot %in% names(fasta_list$Known)),
            yes = "Known",
            no = ifelse(
                test = any(prot %in% names(fasta_list$Contaminant)),
                yes = "Contaminant",
                no = ifelse(
                    test = any(prot %in% names(fasta_list$Novel)),
                    yes = "Novel",
                    no = NA_character_)))) %>%
    dplyr::left_join(x = evid, y = ., by = c("Sequence" = "pep"))

# Define the reverse hits as group
evid_match[evid_match$Reverse == "+", "group"] <- "Reverse"

# Print warning for unidentified sequence origin
print(paste(
    "There are", length(which(is.na(evid_match$group))),
    "NA values, these need to be checked!", sep = " "))

# Clean-up
rm(evid)

# Save the group mapping data
saveRDS(object = evid_match, file = "Sequence_group_mapping.RDS")

# Export complete evidence info for these novel evidence
write.table(
    x = evid_match[evid_match$group == "Novel", ],
    file = paste("Novel_evidence_", date_str, ".txt", sep = ""),
    quote = FALSE,
    sep = "\t",
    row.names = FALSE,
    col.names = TRUE)



### Focus on high quality novel peptide ----------------------------------

# Find the median PEP for known peptides
known.med.pep <- evid_match %>%
    dplyr::filter(., group == "Known") %>%
    dplyr::summarise(., median(PEP)) %>%
    as.numeric(.)

# Keep novel peptides that have lower PEP than known peptide median PEP
candidates <- evid_match %>%
    dplyr::filter(., PEP <= known.med.pep, group == "Novel") %>%
    base::as.data.frame(., stringsAsFactors = FALSE)

# Get position of the novel peptide within ORF
tmp <- evid_match %>%
    dplyr::filter(., group == "Novel") %>%
    dplyr::select(., Sequence) %>%
    unique(.) %>%
    dplyr::group_by(., Sequence) %>%
    dplyr::summarise(
        .,
        Proteins = grep(Sequence, fasta_list$Novel) %>%
            names(fasta_list$Novel)[.] %>%
            paste(., collapse = ";")) %>%
    cSplit(
        indt = ., splitCols = "Proteins", sep = ";", direction = "long") %>%
    base::as.data.frame(., stringsAsFactors = FALSE)
pep.pos <- apply(X = tmp, MARGIN = 1, FUN = function(x) {
    val <- str_locate_all(
        string = fasta_list$Novel[as.character(x[["Proteins"]])] %>% as.character(.),
        pattern = x[["Sequence"]] %>% as.character(.)) %>%
        unlist(.) %>%
        c(x[["Sequence"]], x[["Proteins"]], .)
}) %>%
    unlist(.) %>%
    t(.) %>%
    set_colnames(c("Sequence", "ORF", "start", "end")) %>%
    base::as.data.frame(., stringsAsFactors = FALSE)
pep.pos$start <- as.numeric(pep.pos$start)
pep.pos$end <- as.numeric(pep.pos$end)

# Add the reason for peptide novelty (or indicate if filtered out)
pep.pos <- pep.pos %>%
    dplyr::mutate(
        .,
        ReasonNovel = ifelse(
            test = Sequence %in% candidates$Sequence,
            yes = "",
            no = "PEPfilter")) %>%
    base::as.data.frame(., stringsAsFactors = FALSE)

# Clean-up
rm(tmp)



### Levenshtein distance -------------------------------------------------

# Keep only sequence for novel peptide
tmp <- evid_match %>%
    dplyr::filter(., group == "Novel") %>%
    dplyr::select(., Sequence) %>%
    unique(.)

# Compute the levenshtein distance for all novel peptide and keep
# the minimum leven score result per peptide
leven.data <- adist(
    x = tmp$Sequence,
    y = fasta_list$Known,
    partial = TRUE,
    ignore.case = TRUE) %>%
    set_rownames(tmp$Sequence) %>%
    t(.) %>%
    base::data.frame(
        id = rownames(.),
        .,
        stringsAsFactors = FALSE) %>%
    tidyr::gather(data = ., key = "Sequence", value = "leven", -id) %>%
    dplyr::group_by(., Sequence) %>%
    dplyr::filter(., leven == min(leven)) %>%
    base::as.data.frame(., stringsAsFActors = FALSE)

# Clean-up
rm(tmp)



### Novel ORF blast analysis ---------------------------------------------

# Read blast result from top candidates matched ORF against
# the reference proteome used in this study
blast.bsu.vs.ref <- read.table(
    file = "./blastp_BsuRef_novel_candidates_21122016",
    header = FALSE, sep = "\t", quote = "",
    col.names = c(
        "qseqid", "sseqid", "pident", "nident", "mismatch", "length",
        "gapopen", "qstart", "qend", "sstart", "send", "evalue",
        "bitscore", "score"),
    as.is = TRUE)

# Get the best blastp match for each query
blast.bsu.vs.ref.best <- best_blast(data = blast.bsu.vs.ref, key = "qseqid")

# Filter out the best match that have e-value above 0.0001
blast.bsu.vs.ref.best <- blast.bsu.vs.ref.best %>%
    dplyr::filter(., evalue < 0.0001)
blast.bsu.vs.ref.best$qseqid <- gsub(
    pattern = "^lcl\\|", replacement = "",
    x = blast.bsu.vs.ref.best$qseqid)

# Read blast result from top candidates matched ORF against
# all bacteria proteomes from uniprot
blast.bsu.vs.allbact <- read.table(
    file = "./blastp_allprot_novel_candidates_21122016",
    header = FALSE, sep = "\t", quote = "",
    col.names = c(
        "qseqid", "sseqid", "pident", "nident", "mismatch", "length",
        "gapopen", "qstart", "qend", "sstart", "send", "evalue",
        "bitscore", "score"),
    as.is = TRUE)

# Get the best blastp match for each query
blast.bsu.vs.allbact.best <- best_blast(
    data = blast.bsu.vs.allbact, key = "qseqid")

# Filter out the best match that have e-value above 0.0001
blast.bsu.vs.allbact.best <- blast.bsu.vs.allbact.best %>%
    dplyr::filter(., evalue < 0.0001)
blast.bsu.vs.allbact.best$qseqid <- gsub(
    pattern = "^lcl\\|", replacement = "",
    x = blast.bsu.vs.allbact.best$qseqid)

# Clean-up
rm(blast.bsu.vs.ref, blast.bsu.vs.allbact)



### Novelty reason determination -----------------------------------------

# Define the reason for novel peptide, second add the peptide that
# are novel due to SAV
for (x in 1:nrow(pep.pos)) {
    
    # Process only peptide with no reasons for novelty
    if (pep.pos$ReasonNovel[x] == "") {
        
        seq.pep <- pep.pos$Sequence[x] %>% as.character(.)
        start.pep <- pep.pos$start[x]
        end.pep <- pep.pos$end[x]
        orf.pep <- pep.pos$ORF[x]
        
        tmp.blast.bsu <- blast.bsu.vs.ref.best[
            blast.bsu.vs.ref.best$qseqid == orf.pep, ]
        tmp.blast.allbact <- blast.bsu.vs.allbact.best[
            blast.bsu.vs.allbact.best$qseqid == orf.pep, ]
        tmp.leven.data <- leven.data[
            leven.data$Sequence == seq.pep & leven.data$leven == 1, ]
        
        reason.pep <- ""
        
        # Check whether there are Bsu blast info for current entry
        if (nrow(tmp.blast.bsu) > 0) {
            for (y in 1:nrow(tmp.blast.bsu)) {
                
                # Check if novel peptide is SAV
                if (
                    start.pep >= tmp.blast.bsu$qstart &
                    start.pep <= tmp.blast.bsu$qend &
                    end.pep >= tmp.blast.bsu$qstart &
                    end.pep <= tmp.blast.bsu$qend) {
                    
                    reason.pep <- "Potential SAV"
                    
                    if (
                        nrow(tmp.leven.data) > 0 &
                        any(tmp.leven.data$id %in% tmp.blast.bsu$sseqid)) {
                        
                        reason.pep <- "SAV"
                        
                    }
                    
                    # Check if novel peptide is novel start site
                } else if (
                    (start.pep < tmp.blast.bsu$qstart &
                     end.pep >= tmp.blast.bsu$qstart) |
                    (start.pep < tmp.blast.bsu$qstart &
                     end.pep < tmp.blast.bsu$qstart)) {
                    
                    reason.pep <- "New start site"
                    
                    # Check if novel peptide is novel stop site
                } else if (
                    (start.pep <= tmp.blast.bsu$qend &
                     end.pep > tmp.blast.bsu$qend) |
                    (start.pep > tmp.blast.bsu$qend &
                     end.pep > tmp.blast.bsu$qend)) {
                    
                    reason.pep <- "New stop site"
                    
                }
                
            }
            
            # Check whether there are allbact blast info for current entry
        } else if (nrow(tmp.blast.allbact) > 0) {
            
            reason.pep <- "Known other bacteria"
            
            # If no blast info consider peptide as novel
        } else {
            
            reason.pep <- "Potentially novel"
            
        }
        
        pep.pos$ReasonNovel[x] <- reason.pep
        
    }
    
}

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

# Export condensed dataframe of novel peptide
write.table(
    x = pep.pos.final,
    file = paste("Novel_peptide_", date_str, ".txt", sep = ""),
    quote = FALSE,
    sep = "\t",
    row.names = FALSE,
    col.names = TRUE)

# Clean-up
rm(data)



### Biological explanation of novel ORF ----------------------------------

# Count the number of unique peptide matching novel ORF
filt <- pep.pos[pep.pos$ReasonNovel != "PEPfilter", "ORF"]
data <- pep.pos %>%
    dplyr::filter(., ORF %in% filt) %>%
    dplyr::group_by(., ORF) %>%
    dplyr::summarise(
        .,
        Sequences = toString(x = unique(Sequence), width = NULL),
        ReasonNovel = toString(x = unique(ReasonNovel), width = NULL),
        Unique_peptide_count = n_distinct(Sequence)) %>%
    base::as.data.frame(., stringsAsFactors = FALSE) %>%
    .[order(.$Unique_peptide_count, decreasing = TRUE), ]
tmp <- blast.bsu.vs.allbact.best
colnames(tmp)[colnames(tmp) == "qseqid"] <- "ORF"

# Look into the novel ORF that match known other bacterial entries or are
# completely uncharacterised
orf.candidates <- data[grep(
    pattern = paste(
        c(
            "Known other bacteria", "Potentially novel",
            "New start site", "SAV"),
        collapse = "|"),
    x = data$ReasonNovel), ] %>%
    dplyr::filter(., Unique_peptide_count > 1) %>%
    dplyr::left_join(x = ., y = tmp, by = "ORF") %>%
    .[order(.$Unique_peptide_count, decreasing = TRUE), ] %>%
    base::as.data.frame(., stringsAsFactors = FALSE)

# Export the table of novel ORF that needs validation
write.table(
    x = tmp,
    file = paste("Novel_ORF_toInterprete_", date_str, ".txt", sep = ""),
    quote = FALSE,
    sep = "\t",
    row.names = FALSE,
    col.names = TRUE)

# Clean-up
rm(data, tmp)



### ORF genomic position identification ----------------------------------

# List the fasta files that need to be imported
fasta_file <- choose.files(
    caption = "Select Fasta files",
    multi = TRUE,
    filters = ".fasta") %>%
    set_names(c("NicolasORF"))

# Import all fasta file data and store into list
for (x in 1:length(fasta_file)) {
    
    # Import the current fasta file
    tmp <- read.fasta(
        file = fasta_file[x], seqtype = "AA", as.string = TRUE)
    
    # Include imported fasta into the list
    fasta_list[names(fasta_file)[x]] <- list(tmp)
    
}

# Get the start-end genomic position for each ORFs
orf.pos <- lapply(
    X = fasta_list$NicolasORF, FUN = function(x) { attr(x, "Annot") }) %>%
    unlist(.) %>%
    sub("^.*? \\[(.*?)\\] .*$", "\\1", .) %>%
    base::data.frame(id = names(.), pos = ., stringsAsFactors = FALSE) %>%
    tidyr::separate(
        data = ., col = pos, into = c("start", "end"), sep = " - ")

# Format to numeric the ORFs position
orf.pos$start <- as.numeric(orf.pos$start)
orf.pos$end <- as.numeric(orf.pos$end)

# Filter out the ORF in the fasta that are not holding codons (multiple of 3)
orf.pos %<>%
    dplyr::mutate(
        ., lengthCheck = keep_decimals(x = (abs(x = (start - end)) + 1) / 3)) %>%
    dplyr::filter(., lengthCheck == 0) %>%
    dplyr::select(., -lengthCheck) %>%
    base::as.data.frame(., stringsAsFactors = FALSE)

# Determine strand and frame of ORF
orf.pos$strand <- ifelse(
    test = orf.pos$start < orf.pos$end, yes = 1, no = -1)
orf.pos[orf.pos$strand == 1, "frame"] <- ((as.numeric(
    orf.pos[orf.pos$strand == 1, "start"]) + 2) / 3) %>%
    keep_decimals(.) %>%
    round(x = ((. + 0.33) * 3))
orf.pos[orf.pos$strand == -1, "frame"] <- ((as.numeric(
    orf.pos[orf.pos$strand == -1, "end"]) + 2) / 3) %>%
    keep_decimals(.) %>%
    round(x = ((. + 0.33) * -3))

# Import the blast results of NN ORF versus KK ORF, NN ORF have positional
# information which can be used for KK ORF
blast.NN.vs.KK <- read.table(
    file = "blastp_ORF_Nicolas_vs_Karsten_11012017",
    header = FALSE,
    sep = "\t",
    quote = "",
    col.names = c(
        "qseqid", "sseqid", "pident", "nident", "mismatch", "length",
        "gapopen", "qstart", "qend", "sstart", "send", "evalue",
        "bitscore", "score"),
    as.is = TRUE)

# Get the best blastn match for each query
# (pident == 100%, nident == length, qstart == sstart, qend == send)
blast.NN.vs.KK.best <- blast.NN.vs.KK %>%
    dplyr::group_by(., qseqid) %>%
    unique(.) %>%
    dplyr::filter(
        .,
        pident == 100 & nident == length &
            qstart == sstart & qend == send) %>%
    dplyr::mutate(., count_entry = n()) %>%
    base::as.data.frame(., stringsAsFactors = FALSE)

# Add the entries that are missing from previous step to find best match
tmp <- blast.NN.vs.KK[!blast.NN.vs.KK$qseqid %in% blast.NN.vs.KK.best$qseqid, ]
tmp$count_entry <- 0
blast.NN.vs.KK.best <- rbind(blast.NN.vs.KK.best, tmp)

# Add start-end positions of ORFs to blast results
blast.NN.vs.KK.best %<>%
    dplyr::mutate(., id = sub("^lcl\\|", "", qseqid)) %>%
    dplyr::full_join(x = ., y = orf.pos, by = "id")

# Clean-up
rm(blast.NN.vs.KK, fasta_file, tmp)



### ORF matching to known protein ----------------------------------------

# Import the blast results of Nicolas' ORF versus the reference proteome
blast.NN.vs.ref <- read.table(
    file = "blastp_ORF_Nicolas_vs_RefProt_19012017", header = FALSE,
    sep = "\t",
    quote = "",
    col.names = c(
        "qseqid", "sseqid", "pident", "nident", "mismatch", "length",
        "gapopen", "qstart", "qend", "sstart", "send", "evalue",
        "bitscore", "score"),
    as.is = TRUE)

# Find the best blast hit for each query id
blast.NN.vs.ref.best <- best_blast(data = blast.NN.vs.ref, key = "qseqid")

# Keep the mapping of ORF to uniprotID
orf.uniprot <- blast.NN.vs.ref.best %>%
    dplyr::select(., qseqid, sseqid) %>%
    set_colnames(c("qseqid", "UniProtID")) %>%
    base::as.data.frame(., stringsAsFactors = FALSE)
orf.uniprot$qseqid %<>%
    gsub("^lcl\\|", "", .)
orf.uniprot %<>%
    set_colnames(c("id", "UniProtID"))

# Clean-up
rm(blast.NN.vs.ref)



### ORF neighbours identification ----------------------------------------

# Loop through all frames
neighbour.ORF <- data.frame()
for (x in unique(orf.pos$frame)) {
    
    # Check the strand and use appropriate code to find orf neighbours
    if (x %in% c(1, 2 , 3)) {
        tmp <- orf.pos %>%
            dplyr::filter(., frame == x) %>%
            dplyr::arrange(., start) %>%
            dplyr::select(., id) %>%
            .[["id"]] %>%
            base::data.frame(
                id = .,
                FiveNeighbID = c(.[length(.)], .[-length(.)]),
                ThreeNeighbID = c(.[-1], .[1]),
                stringsAsFactors = FALSE)
    } else {
        tmp <- orf.pos %>%
            dplyr::filter(., frame == x) %>%
            .[order(.[["start"]], decreasing = TRUE), ] %>%
            dplyr::select(., id) %>%
            .[["id"]] %>%
            base::data.frame(
                id = .,
                FiveNeighbID = c(.[length(.)], .[-length(.)]),
                ThreeNeighbID = c(.[-1], .[1]),
                stringsAsFactors = FALSE)
    }
    
    neighbour.ORF <- rbind(neighbour.ORF, tmp)
    
}

# Add positions and uniprotID for orf and its neighbours
tmp <- orf.pos[, c("id", "start", "end")] %>%
    dplyr::left_join(x = ., y = orf.uniprot, by = "id")
orf.neighb <- neighbour.ORF %>%
    dplyr::left_join(x = ., y = tmp, by = "id")
tmp %<>%
    set_colnames(c(
        "FiveNeighbID", "FiveNeighbstart",
        "FiveNeighbend", "FiveNeighbUniprot"))
orf.neighb %<>%
    dplyr::left_join(x = ., y = tmp, by = "FiveNeighbID")
tmp %<>%
    set_colnames(c(
        "ThreeNeighbID", "ThreeNeighbstart",
        "ThreeNeighbend", "ThreeNeighbUniprot"))
orf.neighb %<>%
    dplyr::left_join(x = ., y = tmp, by = "ThreeNeighbID") %>%
    dplyr::select(
        ., id, UniProtID, start, end, FiveNeighbID, FiveNeighbUniprot,
        FiveNeighbstart, FiveNeighbend, ThreeNeighbID, ThreeNeighbUniprot,
        ThreeNeighbstart, ThreeNeighbend)
tmp <- blast.NN.vs.KK.best %>%
    dplyr::select(., id, sseqid) %>%
    set_colnames(c("id", "ORF"))
orf.neighb %<>%
    dplyr::left_join(x = ., y = tmp, by = "id")

# Clean-up
rm(tmp)



### Novel ORF explanation by neighbours ----------------------------------

# Compile neighbouring info with the candidate ORF
orf.candidates.final <- orf.candidates %>%
    dplyr::select(
        ., ORF, Sequences, ReasonNovel, Unique_peptide_count) %>%
    unique(.) %>%
    dplyr::left_join(x = ., y = orf.neighb, by = "ORF") %>%
    dplyr::arrange(., start) %>%
    base::as.data.frame(., stringsAsFactors = FALSE)



### Check novel ORF per conditions ---------------------------------------

# Homogeneise the condition labels
exp.design$Conditions %<>%
    sub("Log/Stationary", "Log-Stationary Transition", .) %>%
    sub("Log/Transition", "Log", .) %>%
    sub(" Transition", "", .)

# Expand dataframe to get individual peptide sequence
tmp <- orf.candidates.final %>%
    cSplit(
        indt = ., splitCols = "Sequences", sep = ", ", direction = "long") %>%
    base::as.data.frame(., stringsAsFactors = FALSE)

# Obtain for each peptide all associated condition it was found in and
# compile per ORF all associated conditions
tmp2 <- evid_match %>%
    dplyr::select(., Sequence, `Raw file`) %>%
    unique(.) %>%
    dplyr::filter(., Sequence %in% tmp$Sequences) %>%
    dplyr::left_join(x = ., y = exp.design, by = c("Raw file" = "Name")) %>%
    dplyr::select(., Sequence, Conditions) %>%
    dplyr::left_join(
        x = ., y = tmp[, c("ORF", "Sequences")],
        by = c("Sequence" = "Sequences")) %>%
    dplyr::group_by(., ORF) %>%
    dplyr::summarise_each(
        funs(toString(x = unique(.), width = NULL)))

# Add the ORF-condition map to the top candidates
tmp %<>%
    dplyr::left_join(x = ., y = tmp2[, c("ORF", "Conditions")], by = "ORF") %>%
    dplyr::group_by(., ORF) %>%
    dplyr::summarise_each(
        funs(toString(x = unique(.), width = NULL))) %>%
    base::as.data.frame(., stringsAsFactors = FALSE)
orf.candidates.final <- tmp

# Export the data for easier examination
write.table(
    x = orf.candidates.final,
    file = paste("ORF_candidate_neighbour_", date_str, ".txt", sep = ""),
    quote = FALSE,
    sep = "\t",
    row.names = FALSE,
    col.names = TRUE)

# Clean-up
rm(tmp, tmp2)



### ORF visualisation ----------------------------------------------------

# Dataframe holding genome information for bacillus subtilis
tmp <- base::data.frame(
    Chromosome = 1, Strand = "*", Start = 1, End = 4215606,
    name = "chr1", length = 4215606, Type = TRUE, geno = "AL009126.3",
    stringsAsFactors = FALSE)

# Create an Ideogram (GRanges object) for B. subtilis
bsu.ideo <- with(
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
saveRDS(object = bsu.ideo, file = "./Databases/bsu_AL009126.3_ideo.RDS")

# Format dataframe as genomic position and other info for each reference ORF
tmp <- orf.neighb %>%
    dplyr::filter(., !is.na(UniProtID)) %>%
    dplyr::select(., id, UniProtID, start, end, ORF) %>%
    dplyr::left_join(
        x = .,
        y = orf.pos %>% dplyr::select(., -start, -end),
        by = "id") %>%
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
    dplyr::select(., -start, -end, -strand) %>%
    dplyr::group_by(., id, Start, End, Strand, frame, Chromosome, chr.name) %>%
    summarise_each(
        funs(toString(x = unique(.), width = NULL))) %>%
    base::as.data.frame(., stringsAsFactors = FALSE)

# Create a GRanges object for all reference and identified novel ORF 
ref.grange <- with(
    data = tmp,
    expr = GRanges(
        seqnames = chr.name,
        ranges = IRanges(Start, End),
        strand = Strand))

# Add values and seqinfo to the created GRanges object
values(ref.grange) <- tmp %>%
    dplyr::select(., id, UniProtKBID, ORF, frame, Chromosome)
seqinfo(ref.grange) <- Seqinfo(
    seqnames = tmp$chr.name %>% unique(.) %>% as.character(.),
    seqlengths = tmp$length %>% unique(.) %>% as.integer(.),
    isCircular = tmp$Type %>% unique(.) %>% as.logical(.),
    genome = tmp$geno %>% unique(.) %>% as.character(.))

# Format dataframe as genomic position and other info for each novel ORF
filt <- evid_match %>%
    dplyr::filter(., group == "Novel") %>%
    cSplit(indt = ., splitCols = "Proteins", sep = ";", direction = "long") %>%
    .[["Proteins"]] %>%
    as.character(.) %>%
    unique(.)
tmp <- orf.neighb %>%
    dplyr::filter(., (!is.na(ORF) & ORF %in% filt)) %>%
    dplyr::select(., id, UniProtID, start, end, ORF) %>%
    dplyr::left_join(
        x = .,
        y = orf.pos %>% dplyr::select(., -start, -end),
        by = "id") %>%
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
    dplyr::select(., -start, -end, -strand) %>%
    dplyr::group_by(., id, Start, End, Strand, frame, Chromosome, chr.name) %>%
    summarise_each(
        funs(toString(x = unique(.), width = NULL))) %>%
    base::as.data.frame(., stringsAsFactors = FALSE)

# Create a GRanges object for all reference and identified novel ORF 
orf.grange <- with(
    data = tmp,
    expr = GRanges(
        seqnames = chr.name,
        ranges = IRanges(Start, End),
        strand = Strand))

# Add values and seqinfo to the created GRanges object
values(orf.grange) <- tmp %>%
    dplyr::select(., id, UniProtKBID, ORF, frame, Chromosome)
seqinfo(orf.grange) <- Seqinfo(
    seqnames = tmp$chr.name %>% unique(.) %>% as.character(.),
    seqlengths = tmp$length %>% unique(.) %>% as.integer(.),
    isCircular = tmp$Type %>% unique(.) %>% as.logical(.),
    genome = tmp$geno %>% unique(.) %>% as.character(.))



### Get peptide position within proteins ---------------------------------

# Compile the required fasta files (novel and known)
fastas <- c(fasta_list$Novel, fasta_list$Known)

# Loop through all enzymatic cleavage rules
pep_list <- base::data.frame()
for (cleav in c("K|R", "K", "R")) {
    
    # Get all peptide sequence and position for all proteins
    pep_list <- prot_digest(fasta = fastas, custom = cleav) %>%
        base::rbind(pep_list, .)
    
}

# Clean up the ID to standard uniprot ID
pep_list$Proteins %<>%
    sub(
        pattern = ".+\\|(.+)\\|.+",
        replacement = "\\1",
        x = .)

# Get association list of peptide and protein (remove reverse and contaminant)
pep.prot <- evid_match %>%
    mq_rev_con_filt(.) %>%
    dplyr::select(., Sequence, Proteins) %>%
    cSplit(indt = ., splitCols = "Proteins", sep = ";", direction = "long") %>%
    unique(.) %>%
    base::as.data.frame(., stringsAsFactors = FALSE)

# Get peptide location within each associated proteins for identified peptides
pep_located <- pep.prot %>%
    dplyr::left_join(x = ., y = pep_list, by = c("Sequence", "Proteins")) %>%
    base::as.data.frame(., stringsAsFactors = FALSE)

# Select the peptide that could not be located and add Methionine
# at the start of the peptide sequence
tmp <- pep_located %>%
    dplyr::filter(., is.na(start)) %>%
    dplyr::select(., Sequence, Proteins) %>%
    set_colnames(c("OrigSequence", "Proteins")) %>%
    dplyr::mutate(., Sequence = sub("^", "M", OrigSequence))

# Keep only located peptide from this variable
pep_located %<>%
    dplyr::filter(., !is.na(start))

# Try to locate again the missing peptide (the Methionine addition should help)
# and add the new results to main variable
pep_located <- tmp %>%
    dplyr::left_join(x = ., y = pep_list, by = c("Sequence", "Proteins")) %>%
    dplyr::mutate(., Sequence = OrigSequence, start = (start + 1)) %>%
    dplyr::select(., -OrigSequence) %>%
    base::as.data.frame(., stringsAsFactors = FALSE) %>%
    base::rbind(pep_located, ., stringsAsFactors = FALSE)

# Remove any remaining contaminant entry (which also could not be located)
pep_located <- pep_located[
    grep("^CON__", pep_located$Proteins, invert = TRUE), ] %>%
    unique(.)



### Table of novel ORF manual annotation ---------------------------------

# Import the table of manul annotation
manual.annot <- read.table(
    file = "ORF_candidate_neighbour_2017-02-13.txt",
    header = TRUE, sep = "\t", quote = "", as.is = TRUE)
colnames(manual.annot)[
    grep("Unique_peptide", colnames(manual.annot))] <- "PeptideCount"



### Generate the report --------------------------------------------------

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

# Clean-up
rm(param)


