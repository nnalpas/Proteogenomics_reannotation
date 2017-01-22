


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

# Save the group mapping data
saveRDS(object = evid.match, file = "Sequence_group_mapping.RDS")



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






#
known.med.pep <- evid.match %>%
    dplyr::filter(., group == "Known") %>%
    dplyr::select(., Sequence, PEP) %>%
    dplyr::summarise(., median(PEP)) %>%
    as.numeric(.)

# Keep novel peptides with median PEP (based on known)
candidates <- evid.match %>%
    dplyr::filter(., PEP <= known.med.pep, group == "Novel") %>%
    base::as.data.frame(., stringsAsFactors = FALSE)

# Find the ORF where novel peptide map to
candidate.orf <- data.frame()
for (x in unique(evid.match[evid.match$group == "Novel", "Sequence"])) {
    
    orfs <- grep(pattern = x, x = fasta.6frame) %>%
        names(fasta.6frame)[.]
    candidate.orf %<>%
        rbind(
            ., data.frame(Sequence = x, ORF = orfs, stringsAsFactors = FALSE))
    
}

# Get position of the novel peptide within ORF
candidate.pos <- data.frame()
for (x in 1:nrow(candidate.orf)) {
    
    tmp <- str_locate_all(
        string = fasta.6frame[candidate.orf$ORF[x]],
        pattern = candidate.orf$Sequence[x]) %>%
        set_names(names(fasta.6frame[candidate.orf$ORF[x]])) %>%
        ldply(., data.frame) %>%
        set_colnames(c("ORF", "start", "end"))
    candidate.pos %<>%
        rbind(., data.frame(
            Sequence = candidate.orf$Sequence[x],
            tmp, stringsAsFactors = TRUE))
}
rm(candidate.orf)

# Export the list of ORF identified by a novel peptide
write.table(
    x = candidate.pos[
        candidate.pos$Sequence %in% unique(candidates$Sequence),
        "ORF"],
    file = paste("Novel_pep_orfs_", date.str, ".txt", sep = ""),
    quote = FALSE,
    sep = "\t",
    row.names = FALSE,
    col.names = FALSE)

# 
data <- candidate.pos %>%
    dplyr::summarise(
        .,
        type = "All novel",
        Number_peptide = n_distinct(Sequence),
        Number_ORF = n_distinct(ORF)) %>%
    tidyr::gather(
        data = ., key = "feature", value = "Count", -type) %>%
    base::as.data.frame(., stringsAsFactors = FALSE)

# 
data <- candidate.pos %>%
    dplyr::filter(., ReasonNovel != "PEPfilter") %>%
    dplyr::summarise(
        .,
        type = "PEP filtered",
        Number_peptide = n_distinct(Sequence),
        Number_ORF = n_distinct(ORF)) %>%
    tidyr::gather(
        data = ., key = "feature", value = "Count", -type) %>%
    base::as.data.frame(., stringsAsFactors = FALSE) %>%
    rbind(data, .)

# 
histPlots(
    data = data,
    key = "feature",
    value = "Count",
    group = "type",
    fill = "type",
    main = "Novel peptide summary")

# Read blast result from top candidates matched ORF against
# the reference proteome used in this study
blast.bsu <- read.table(
    file = "../blastp_BsuRef_novel_candidates_21122016",
    header = FALSE, sep = "\t", quote = "",
    col.names = c(
        "qseqid", "sseqid", "pident", "nident", "mismatch", "length",
        "gapopen", "qstart", "qend", "sstart", "send", "evalue",
        "bitscore", "score"),
    as.is = TRUE)

# Get the best blastp match for each query
blast.bsu.final <- data.frame()
for (id in unique(blast.bsu$qseqid)) {
    
    tmp <- blast.bsu[blast.bsu$qseqid == id, ] %>%
        unique(.)
    min.eval <- min(tmp$evalue)
    tmp <- tmp[tmp$evalue == min.eval, ]
        
    if (nrow(tmp) == 1) {
        blast.bsu.final <- rbind(blast.bsu.final, tmp)
    } else {
        max.score <- max(tmp$score)
        tmp <- tmp[tmp$score == max.score, ]
        if (nrow(tmp) == 1) {
            blast.bsu.final <- rbind(blast.bsu.final, tmp)
        } else {
            max.pident <- max(tmp$pident)
            tmp <- tmp[tmp$pident == max.pident, ]
            if (nrow(tmp) == 1) {
                blast.bsu.final <- rbind(blast.bsu.final, tmp)
            } else {
                print(paste("Cannot determine best match for: ", id, sep = ""))
            }
        }
    }
    
}

# Filter out the best match that have e-value above 0.0001
blast.bsu.final <- blast.bsu.final[blast.bsu.final$evalue < 0.0001, ]
blast.bsu.final$qseqid <- gsub(
    pattern = "^lcl\\|", replacement = "", x = blast.bsu.final$qseqid)

# Read blast result from top candidates matched ORF against
# all bacteria proteomes from uniprot
blast.allbact <- read.table(
    file = "../blastp_allprot_novel_candidates_21122016",
    header = FALSE, sep = "\t", quote = "",
    col.names = c(
        "qseqid", "sseqid", "pident", "nident", "mismatch", "length",
        "gapopen", "qstart", "qend", "sstart", "send", "evalue",
        "bitscore", "score"),
    as.is = TRUE)

# Get the best blastp match for each query
blast.allbact.final <- data.frame()
for (id in unique(blast.allbact$qseqid)) {
    
    tmp <- blast.allbact[blast.allbact$qseqid == id, ] %>%
        unique(.)
    min.eval <- min(tmp$evalue)
    tmp <- tmp[tmp$evalue == min.eval, ]
    
    if (nrow(tmp) == 1) {
        blast.allbact.final <- rbind(blast.allbact.final, tmp)
    } else {
        max.score <- max(tmp$score)
        tmp <- tmp[tmp$score == max.score, ]
        if (nrow(tmp) == 1) {
            blast.allbact.final <- rbind(blast.allbact.final, tmp)
        } else {
            max.pident <- max(tmp$pident)
            tmp <- tmp[tmp$pident == max.pident, ]
            blast.allbact.final <- rbind(blast.allbact.final, tmp)
            if (nrow(tmp) != 1) {
                print(paste("Cannot determine best match for: ", id, sep = ""))
            }
        }
    }
    
}

# Filter out the best match that have e-value above 0.0001
blast.allbact.final <- blast.allbact.final[blast.allbact.final$evalue < 0.0001, ]
blast.allbact.final$qseqid <- gsub(
    pattern = "^lcl\\|", replacement = "", x = blast.allbact.final$qseqid)

# Define the reason for novel peptide, first add the peptide that
# were filtered due to high PEP
candidate.pos$ReasonNovel <- NA
candidate.pos[
    !(candidate.pos$Sequence %in% unique(candidates$Sequence)),
    "ReasonNovel"] <- "PEPfilter"

# Define the reason for novel peptide, second add the peptide that
# are novel due to SAV
for (x in 1:nrow(candidate.pos)) {
    
    # Process only peptide with no reasons for novelty
    if (is.na(candidate.pos$ReasonNovel[x])) {
        
        seq.pep <- candidate.pos$Sequence[x] %>% as.character(.)
        start.pep <- candidate.pos$start[x]
        end.pep <- candidate.pos$end[x]
        orf.pep <- candidate.pos$ORF[x]
        
        tmp.blast.bsu <- blast.bsu.final[
            blast.bsu.final$qseqid == orf.pep, ]
        tmp.blast.allbact <- blast.allbact.final[
            blast.allbact.final$qseqid == orf.pep, ]
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
        
        candidate.pos$ReasonNovel[x] <- reason.pep
        
    }
    
}

# 
data <- candidate.pos %>%
    dplyr::filter(., ReasonNovel != "PEPfilter") %>%
    dplyr::group_by(., ReasonNovel) %>%
    dplyr::summarise(
        .,
        Number_peptide = n_distinct(Sequence),
        Number_ORF = n_distinct(ORF)) %>%
    tidyr::gather(
        data = ., key = "feature", value = "Count", -ReasonNovel) %>%
    base::as.data.frame(., stringsAsFactors = FALSE)

#
histPlots(
    data = data,
    key = "feature",
    value = "Count",
    group = "ReasonNovel",
    fill = "ReasonNovel",
    main = "Peptide novelty reasons")

# Investigate the potentially fully uncharacterised novel peptide
novelty <- unique(candidate.pos$ReasonNovel)[-1]
data <- data.frame()
for (key in novelty) {
    
    filt <- candidate.pos[candidate.pos$ReasonNovel == key, "ORF"]
    data <- candidate.pos %>%
        dplyr::filter(., ORF %in% filt) %>%
        dplyr::group_by(., ORF) %>%
        dplyr::summarise(., Unique_peptide_count = n_distinct(Sequence)) %>%
        dplyr::group_by(., Unique_peptide_count) %>%
        dplyr::summarise(., ORF_count = n_distinct(ORF)) %>%
        base::data.frame(., Novelty = key, stringsAsFActors = FALSE) %>%
        rbind(data, .)
    
}

#
histPlots(
    data = data,
    key = "Unique_peptide_count",
    value = "ORF_count",
    group = "Novelty",
    fill = "Novelty",
    main = "Count of novel ORF per count of unique novel peptide")

# Close the plot output
dev.off()

# Export complete evidence info for thee novel evidence
write.table(
    x = evid.match[evid.match$group == "Novel", ],
    file = paste("Novel_evidence_", date.str, ".txt", sep = ""),
    quote = FALSE,
    sep = "\t",
    row.names = FALSE,
    col.names = TRUE)

# Compile a condensed dataframe of novel peptide info and export it
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
candidate.pos.final <- merge(
    x = candidate.pos, y = data, by = "Sequence", all.x = TRUE)
write.table(
    x = candidate.pos.final,
    file = paste("Novel_peptide_", date.str, ".txt", sep = ""),
    quote = FALSE,
    sep = "\t",
    row.names = FALSE,
    col.names = TRUE)



### Biological explanation of novel ORF ----------------------------------

# 
filt <- candidate.pos[candidate.pos$ReasonNovel != "PEPfilter", "ORF"]
data <- candidate.pos %>%
    dplyr::filter(., ORF %in% filt) %>%
    dplyr::group_by(., ORF) %>%
    dplyr::summarise(
        .,
        Sequences = toString(x = unique(Sequence), width = NULL),
        ReasonNovel = toString(x = unique(ReasonNovel), width = NULL),
        Unique_peptide_count = n_distinct(Sequence)) %>%
    base::as.data.frame(., stringsAsFactors = FALSE) %>%
    .[order(.$Unique_peptide_count, decreasing = TRUE), ]

# Look first into the novel ORF that match known other bacterial entries
tmp <- data[grep(pattern = "Known other bacteria", x = data$ReasonNovel), ] %>%
    dplyr::filter(., Unique_peptide_count > 1) %>%
    merge(
        x = ., y = blast.allbact.final,
        by.x = "ORF", by.y = "qseqid", all.x = TRUE)

# 
tmp <- data[grep(pattern = "Novel", x = data$ReasonNovel), ] %>%
    dplyr::filter(., Unique_peptide_count > 1) %>%
    merge(
        x = ., y = blast.allbact.final,
        by.x = "ORF", by.y = "qseqid", all.x = TRUE) %>%
    rbind(tmp, .) %>%
    .[order(.$Unique_peptide_count, decreasing = TRUE), ] %>%
    base::as.data.frame(., stringsAsFactors = FALSE)

# 
write.table(
    x = tmp,
    file = paste("Novel_ORF_toInterprete_", date.str, ".txt", sep = ""),
    quote = FALSE,
    sep = "\t",
    row.names = FALSE,
    col.names = TRUE)

# to re-implement as wrapper function
#ggplot(
#    data = data.filt,
#    mapping = aes(x = PEP, fill = factor(group), colour = factor(group))) +
#    geom_density(alpha = 0.1) +
#    xlab(label = "Evidence PEP") +
#    ylab(label = "Density") +
#    ggtitle(label = "Density of evidence PEP") +
#    xlim(0, 0.05)



### ORF genomic position identification ----------------------------------

# Import the blast results
blast <- read.table(
    file = "blastp_ORF_Nicolas_vs_Karsten_11012017.txt",
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
blast.NN.KK.final <- blast %>%
    dplyr::group_by(., qseqid) %>%
    unique(.) %>%
    dplyr::filter(
        .,
        pident == 100 & nident == length &
            qstart == sstart & qend == send) %>%
    dplyr::mutate(., count_entry = n()) %>%
    base::as.data.frame(., stringsAsFactors = FALSE)

# Add the entries that are missing from previous step to find best match
tmp <- blast[!blast$qseqid %in% blast.NN.KK.final$qseqid, ]
tmp$count_entry <- 0
blast.NN.KK.final <- rbind(blast.NN.KK.final, tmp)

# Check how many match have the query entries
tmp <- blast.NN.KK.final %>%
    dplyr::select(., qseqid, count_entry) %>%
    dplyr::group_by(., qseqid) %>%
    unique(.) %>%
    dplyr::ungroup(.) %>%
    dplyr::group_by(., count_entry) %>%
    dplyr::summarise(., table = n()) %>%
    base::as.data.frame(., stringsAsFactors = FALSE)

# Import the fasta files
fasta.nn <- read.fasta(
    file = "G:/data/Vaishnavi/Databases/Find0_GCA_000009045.1.fasta",
    seqtype = "AA",
    as.string = TRUE)
fasta.kk <- read.fasta(
    file = paste(
        "G:/data/Vaishnavi/Databases/",
        "Bsu_genome_assembly_GCA_000009045.1.out_FIXED_HEADER.fasta",
        sep = ""),
    seqtype = "AA",
    as.string = TRUE)

# Get the start-end genomic position for each ORFs
positions <- lapply(X = fasta.nn, FUN = function(x) { attr(x, "Annot") }) %>%
    unlist(.) %>%
    sub("^.*? \\[(.*?)\\] .*$", "\\1", .) %>%
    base::data.frame(id = names(.), pos = ., stringsAsFactors = FALSE) %>%
    tidyr::separate(
        data = ., col = pos, into = c("start", "end"), sep = " - ") %>%
    base::as.data.frame(., stringsAsFactors = FALSE)

# Determine strand and frame of ORF
positions$strand <- ifelse(
    test = positions$start < positions$end, yes = 1, no = -1)
positions[positions$strand == 1, "frame"] <- ((as.numeric(
    positions[positions$strand == 1, "start"]) + 2) / 3) %>%
    revtrunc(.) %>%
    round(x = ((. + 0.33) * 3))
positions[positions$strand == -1, "frame"] <- ((as.numeric(
    positions[positions$strand == -1, "end"]) + 2) / 3) %>%
    revtrunc(.) %>%
    round(x = ((. + 0.33) * -3))

# Add start-end positions of ORFs to blast results
blast.NN.KK.final %<>%
    dplyr::mutate(., id = sub("^lcl\\|", "", qseqid)) %>%
    merge(
    x = ., y = positions,
    by = "id", all = TRUE)

# Import the blast results of Nicolas' ORF versus the reference proteome
blast.vs.ref <- read.table(
    file = "blastp_ORF_Nicolas_vs_RefProt_19012017.txt", header = FALSE,
    sep = "\t",
    quote = "",
    col.names = c(
        "qseqid", "sseqid", "pident", "nident", "mismatch", "length",
        "gapopen", "qstart", "qend", "sstart", "send", "evalue",
        "bitscore", "score"),
    as.is = TRUE)

# Find the best blast hit for each query id
blast.vs.ref.best <- best.blast(data = blast.vs.ref, key = "qseqid") %>%
    dplyr::select(., qseqid, sseqid) %>%
    set_colnames(c("qseqid", "UniProtID")) %>%
    base::as.data.frame(., stringsAsFactors = FALSE)

# Merge the Uniprot ID info to the blast map of ORF to UniprotID
blast.vs.ref.final <- merge(
    x = blast.NN.KK.final,
    y = blast.vs.ref.best,
    by = "qseqid",
    all = TRUE)

# Import the novel peptide table
novel.pep <- read.table(
    file = "Novel_peptide_2016-12-21.txt",
    header = TRUE, sep = "\t", quote = "")

# Check how many novel ORF are missing from the blast results
summary(unique(novel.pep$ORF) %in% unique(blast.NN.KK.final$sseqid))

# Merge the novel peptide with the ORF position and blast info
novel.pep.pos <- merge(
    x = blast.vs.ref.final, y = novel.pep,
    by.x = "sseqid", by.y = "ORF", all.y = TRUE) %>%
    dplyr::mutate(., id = sub("^lcl\\|", "", qseqid)) %>%
    merge(x = ., y = positions, by = "id", all.x = TRUE) %>%
    base::as.data.frame(., stringsAsFactors = FALSE)

# Keep only ORF with positions information and only these columns
data <- novel.pep.pos %>%
    dplyr::filter(., !is.na(`start.y`)) %>%
    dplyr::select(., sseqid, `start.y`, `end.y`, strand, frame) %>%
    unique(.) %>%
    base::as.data.frame(., stringsAsFactors = FALSE)

# Format to numeric the ORFs position
positions$start <- as.numeric(positions$start)
positions$end <- as.numeric(positions$end)
data$start.y <- as.numeric(data$start.y)
data$end.y <- as.numeric(data$end.y)
blast.vs.ref.best$qseqid %<>%
    sub("^lcl\\|", "", .)

# Loop through selected ORF
neighbour.ORF <- data.frame()
for (x in 1:nrow(data)) {
    
    # Filter the full ORF entries for neighbours
    tmp <- positions %>%
        dplyr::filter(
            .,
            strand == data[x, "strand"] & frame == data[x, "frame"])
    
    # Depending on strand detect 5' and 3' neighbour ORF
    if (data[x, "strand"] == 1) {
        
        FiveNeighb <- tmp %>%
            dplyr::filter(
                ., end < data[x, "start.y"]) %>%
            dplyr::filter(., end == max(end)) %>%
            base::as.data.frame(., stringsAsFActors = FALSE)
        ThreeNeighb <- tmp %>%
            dplyr::filter(., start > data[x, "end.y"]) %>%
            dplyr::filter(., start == min(start)) %>%
            base::as.data.frame(., stringsAsFActors = FALSE)
        
    } else {
        
        FiveNeighb <- tmp %>%
            dplyr::filter(
                ., end > data[x, "start.y"]) %>%
            dplyr::filter(., end == min(end)) %>%
            base::as.data.frame(., stringsAsFActors = FALSE)
        ThreeNeighb <- tmp %>%
            dplyr::filter(., start < data[x, "end.y"]) %>%
            dplyr::filter(., start == max(start)) %>%
            base::as.data.frame(., stringsAsFActors = FALSE)
        
    }
    
    # Check number of entry detected as neighbours
    if (nrow(FiveNeighb) == 1 & nrow(ThreeNeighb) == 1) {
        
        # Get UniProt ID for 5' and 3' neighbours if any
        FiveNeighbUniID <- ifelse(
            test = any(blast.vs.ref.best$qseqid == FiveNeighb[["id"]]),
            yes = blast.vs.ref.best[
                blast.vs.ref.best$qseqid == FiveNeighb[["id"]], "UniProtID"],
            no = "")
        ThreeNeighbUniID <- ifelse(
            test = any(blast.vs.ref.best$qseqid == ThreeNeighb[["id"]]),
            yes = blast.vs.ref.best[
                blast.vs.ref.best$qseqid == ThreeNeighb[["id"]], "UniProtID"],
            no = "")
        
        # Add the 5' and 3' ORF neighbours to the current novel ORF
        neighbour.ORF <- data.frame(
            sseqid = data[x, "sseqid"],
            FiveNeighbID = FiveNeighb[["id"]],
            FiveNeighbUniID = FiveNeighbUniID,
            FiveNeighbStart = FiveNeighb[["start"]],
            FiveNeighbEnd = FiveNeighb[["end"]],
            ThreeNeighbID = ThreeNeighb[["id"]],
            ThreeNeighbUniID = ThreeNeighbUniID,
            ThreeNeighbStart = ThreeNeighb[["start"]],
            ThreeNeighbEnd = ThreeNeighb[["end"]]) %>%
            rbind(neighbour.ORF, .)
        
    }
    
}

# Continue with adding neighbour info to novel ORF and calculating how many
# peptides confirm the ORF presence
tmp <- novel.pep.pos %>%
    merge(x = ., y = neighbour.ORF, by = "sseqid", all.x = TRUE) %>%
    dplyr::group_by(., sseqid, id) %>%
    dplyr::mutate(., ORFpepCount = n()) %>%
    base::as.data.frame(., stringsAsFactors = FALSE)



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
    
}


