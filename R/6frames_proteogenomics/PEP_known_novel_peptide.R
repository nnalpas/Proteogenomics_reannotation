


### Parameters setting up ------------------------------------------------

# Start with clean environment
rm(list = ls())

# Define user input parameters
markdown <- TRUE



### Define working directory ---------------------------------------------

# Define the work space
work.space <- choose.dir()
setwd(work.space)

# Define the current user
user <- Sys.info()[["user"]]

# Define current time
date.str <- format(Sys.time(), "%Y-%m-%d")



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
loadpackage(plyr)
loadpackage(dplyr)
loadpackage(tidyr)
loadpackage(seqinr)
loadpackage(UniProt.ws)
loadpackage(magrittr)
loadpackage(WriteXLS)
loadpackage(data.table)
loadpackage(splitstackshape)
loadpackage(VennDiagram)
loadpackage(ggplot2)
loadpackage(grid)
loadpackage(gridExtra)
loadpackage(RColorBrewer)
loadpackage(stringr)
loadpackage(Biostrings)
loadpackage(RecordLinkage)
loadpackage(VariantAnnotation)
loadpackage(cgdsr)
loadpackage(bit64)
loadpackage(cleaver)
loadpackage(plotly)
loadpackage(GenomicRanges)
loadpackage(biovizBase)
loadpackage(ggbio)



### Data import ----------------------------------------------------------

# Import the maxquant evidence table
evid <- maxquant.read(
    path = ".",
    name = "evidence.txt",
    integer64 = "double")

# List the fasta files that need to be imported
fasta.file <- choose.files(
    caption = "Select Fasta files",
    multi = TRUE,
    filters = ".fasta") %>%
    set_names(c("Novel", "Contaminant", "Known"))

# Import all fasta file data and store into list
fasta.list <- list()
for (x in 1:length(fasta.file)) {
    
    # Import the current fasta file
    tmp <- read.fasta(
        file = fasta.file[x], seqtype = "AA", as.string = TRUE)
    
    # Include imported fasta into the list
    fasta.list[names(fasta.file)[x]] <- list(tmp)
    
}

# Clean-up
rm(fasta.file)



### Novel peptide identification -----------------------------------------

# Digest all protein and store peptide into list
pep.list <- lapply(X = fasta.list, FUN = function(x) {
    
    cliv <- lapply(
        X = c("K|R", "K", "R"),
        FUN = function(y) {
            tmp <- cleave(
                x = x %>% unlist(.),
                custom = y,
                missedCleavages = c(0:2)) %>%
                unlist(.) %>%
                unique(.)
        }) %>%
        unlist(.) %>%
        unique(.)
    
    cliv
})

# New dataframe to hold info about fasta of origin for each sequence
evid.match <- evid %>%
    dplyr::mutate(
        .,
        group = ifelse(
            test = Sequence %in% pep.list$Known,
            yes = "Known",
            no = ifelse(
                test = Sequence %in% pep.list$Contaminant,
                yes = "Contaminant",
                no = ifelse(
                    test = Sequence %in% pep.list$Novel,
                    yes = "Novel",
                    no = NA_character_)))) %>%
    base::as.data.frame(., stringAsFactors = TRUE)

# Reperform previous step to try find the sequence that miss first methionine
data <- evid.match %>%
    dplyr::filter(., is.na(group)) %>%
    dplyr::mutate(
        .,
        group = ifelse(
            test = paste(
                "M", Sequence, sep = "") %in% pep.list$Known,
            yes = "Known",
            no = ifelse(
                test = paste(
                    "M", Sequence, sep = "") %in% pep.list$Contaminant,
                yes = "Contaminant",
                no = ifelse(
                    test = paste(
                        "M", Sequence, sep = "") %in% pep.list$Novel,
                    yes = "Novel",
                    no = NA_character_)))) %>%
    base::as.data.frame(., stringAsFactors = TRUE)

# Combine the two peptide group type information
evid.match <- rbind(evid.match[!is.na(evid.match$group), ], data)

# Define the reverse hits as group
evid.match[evid.match$Reverse == "+", "group"] <- "Reverse"

# Print warning for unidentified sequence origin
print(paste(
    "There are", length(which(is.na(evid.match$group))),
    "NA values, these need to be checked!", sep = " "))

# Clean-up
rm(pep.list, evid)

# Save the group mapping data
saveRDS(object = evid.match, file = "Sequence_group_mapping.RDS")

# Export complete evidence info for these novel evidence
write.table(
    x = evid.match[evid.match$group == "Novel", ],
    file = paste("Novel_evidence_", date.str, ".txt", sep = ""),
    quote = FALSE,
    sep = "\t",
    row.names = FALSE,
    col.names = TRUE)



### Focus on high quality novel peptide ----------------------------------

# Find the median PEP for known peptides
known.med.pep <- evid.match %>%
    dplyr::filter(., group == "Known") %>%
    dplyr::summarise(., median(PEP)) %>%
    as.numeric(.)

# Keep novel peptides that have lower PEP than known peptide median PEP
candidates <- evid.match %>%
    dplyr::filter(., PEP <= known.med.pep, group == "Novel") %>%
    base::as.data.frame(., stringsAsFactors = FALSE)

# Get position of the novel peptide within ORF
tmp <- evid.match %>%
    dplyr::filter(., group == "Novel") %>%
    dplyr::select(., Sequence) %>%
    unique(.) %>%
    dplyr::group_by(., Sequence) %>%
    dplyr::summarise(
        .,
        Proteins = grep(Sequence, fasta.list$Novel) %>%
            names(fasta.list$Novel)[.] %>%
            paste(., collapse = ";")) %>%
    cSplit(
        indt = ., splitCols = "Proteins", sep = ";", direction = "long") %>%
    base::as.data.frame(., stringsAsFactors = FALSE)
pep.pos <- apply(X = tmp, MARGIN = 1, FUN = function(x) {
    val <- str_locate_all(
        string = fasta.list$Novel[as.character(x[["Proteins"]])] %>% as.character(.),
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
tmp <- evid.match %>%
    dplyr::filter(., group == "Novel") %>%
    dplyr::select(., Sequence) %>%
    unique(.)

# Compute the levenshtein distance for all novel peptide and keep
# the minimum leven score result per peptide
leven.data <- adist(
    x = tmp$Sequence,
    y = fasta.list$Known,
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
blast.bsu.vs.ref.best <- best.blast(data = blast.bsu.vs.ref, key = "qseqid")

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
blast.bsu.vs.allbact.best <- best.blast(
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
data <- evid.match %>%
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
    file = paste("Novel_peptide_", date.str, ".txt", sep = ""),
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
    file = paste("Novel_ORF_toInterprete_", date.str, ".txt", sep = ""),
    quote = FALSE,
    sep = "\t",
    row.names = FALSE,
    col.names = TRUE)

# Clean-up
rm(data, tmp)



### ORF genomic position identification ----------------------------------

# List the fasta files that need to be imported
fasta.file <- choose.files(
    caption = "Select Fasta files",
    multi = TRUE,
    filters = ".fasta") %>%
    set_names(c("NicolasORF"))

# Import all fasta file data and store into list
for (x in 1:length(fasta.file)) {
    
    # Import the current fasta file
    tmp <- read.fasta(
        file = fasta.file[x], seqtype = "AA", as.string = TRUE)
    
    # Include imported fasta into the list
    fasta.list[names(fasta.file)[x]] <- list(tmp)
    
}

# Get the start-end genomic position for each ORFs
orf.pos <- lapply(
    X = fasta.list$NicolasORF, FUN = function(x) { attr(x, "Annot") }) %>%
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
        ., lengthCheck = revtrunc(x = (abs(x = (start - end)) + 1) / 3)) %>%
    dplyr::filter(., lengthCheck == 0) %>%
    dplyr::select(., -lengthCheck) %>%
    base::as.data.frame(., stringsAsFactors = FALSE)

# Determine strand and frame of ORF
orf.pos$strand <- ifelse(
    test = orf.pos$start < orf.pos$end, yes = 1, no = -1)
orf.pos[orf.pos$strand == 1, "frame"] <- ((as.numeric(
    orf.pos[orf.pos$strand == 1, "start"]) + 2) / 3) %>%
    revtrunc(.) %>%
    round(x = ((. + 0.33) * 3))
orf.pos[orf.pos$strand == -1, "frame"] <- ((as.numeric(
    orf.pos[orf.pos$strand == -1, "end"]) + 2) / 3) %>%
    revtrunc(.) %>%
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
rm(blast.NN.vs.KK, fasta.file, tmp)



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
blast.NN.vs.ref.best <- best.blast(data = blast.NN.vs.ref, key = "qseqid")

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

# Export the data for easier examination
write.table(
    x = orf.candidates.final,
    file = paste("ORF_candidate_neighbour_", date.str, ".txt", sep = ""),
    quote = FALSE,
    sep = "\t",
    row.names = FALSE,
    col.names = TRUE)



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
        UniProtKBID = uniprotID.clean(UniProtID),
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
    seqnames = tmp$chr.name %>% unique(.),
    seqlengths = tmp$length %>% unique(.),
    isCircular = tmp$Type %>% unique(.),
    genome = tmp$geno %>% unique(.))

# Use ggbio extension to plot ORF location on genome as a circos graph
colou <- c("#4682B4", "#BD5E5E", "#437A3C", "#F57C36", "#D58DEB", "#B2B83F")
p <- ggplot() +
    ggtitle(label = "Reference ORF") +
    layout_circle(
        bsu.ideo, geom = "ideo", fill = "gray70",
        radius = 30, trackWidth = 3) +
    layout_circle(
        bsu.ideo, geom = "scale", size = 2, radius = 33, trackWidth = 2) +
    layout_circle(
        bsu.ideo, geom = "text", aes(label = seqnames),
        angle = 0, radius = 36, trackWidth = 5) +
    layout_circle(
        subset(x = ref.grange, frame == 1), geom = "rect", color = colou[1],
        radius = 26, trackWidth = 3) +
    layout_circle(
        subset(x = ref.grange, frame == 2), geom = "rect", color = colou[2],
        radius = 23, trackWidth = 3) +
    layout_circle(
        subset(x = ref.grange, frame == 3), geom = "rect", color = colou[3],
        radius = 20, trackWidth = 3) +
    layout_circle(
        subset(x = ref.grange, frame == -1), geom = "rect", color = colou[4],
        radius = 17, trackWidth = 3) +
    layout_circle(
        subset(x = ref.grange, frame == -2), geom = "rect", color = colou[5],
        radius = 14, trackWidth = 3) +
    layout_circle(
        subset(x = ref.grange, frame == -3), geom = "rect", color = colou[6],
        radius = 11, trackWidth = 3)
p

# Format dataframe as genomic position and other info for each novel ORF
filt <- evid.match %>%
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
        UniProtKBID = uniprotID.clean(UniProtID),
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
    seqnames = tmp$chr.name %>% unique(.),
    seqlengths = tmp$length %>% unique(.),
    isCircular = tmp$Type %>% unique(.),
    genome = tmp$geno %>% unique(.))

# Use ggbio extension to plot ORF location on genome as a circos graph
colou <- c("#4682B4", "#BD5E5E", "#437A3C", "#F57C36", "#D58DEB", "#B2B83F")
p <- ggplot() +
    ggtitle(label = "Novel ORF") +
    layout_circle(
        bsu.ideo, geom = "ideo", fill = "gray70",
        radius = 30, trackWidth = 3) +
    layout_circle(
        bsu.ideo, geom = "scale", size = 2, radius = 33, trackWidth = 2) +
    layout_circle(
        bsu.ideo, geom = "text", aes(label = seqnames),
        angle = 0, radius = 36, trackWidth = 5) +
    layout_circle(
        subset(x = orf.grange, frame == 1), geom = "rect", color = colou[1],
        radius = 26, trackWidth = 3) +
    layout_circle(
        subset(x = orf.grange, frame == 2), geom = "rect", color = colou[2],
        radius = 23, trackWidth = 3) +
    layout_circle(
        subset(x = orf.grange, frame == 3), geom = "rect", color = colou[3],
        radius = 20, trackWidth = 3) +
    layout_circle(
        subset(x = orf.grange, frame == -1), geom = "rect", color = colou[4],
        radius = 17, trackWidth = 3) +
    layout_circle(
        subset(x = orf.grange, frame == -2), geom = "rect", color = colou[5],
        radius = 14, trackWidth = 3) +
    layout_circle(
        subset(x = orf.grange, frame == -3), geom = "rect", color = colou[6],
        radius = 11, trackWidth = 3)
p






### Generate the report --------------------------------------------------

# Check whether markdown output is requested
if (markdown) {
    
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
        evidences = evid.match,
        pept.posit = pep.pos,
        levenshtein = leven.data)
    
    # Define output file name
    out.file <- paste(
        work.space,
        "/Bsu_proteogenomics_",
        format(Sys.time(), '%Y%m%d_%H-%M'),
        ".html",
        sep = "")
    
    # Render the markdown report
    rmarkdown::render(
        input = tempReport,
        output_format = "ioslides_presentation",
        output_file = out.file,
        params = param,
        envir = new.env(parent = globalenv()))
    
    # Clean-up
    rm(param)
    
}


