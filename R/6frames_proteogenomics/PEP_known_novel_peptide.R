
rm(list = ls())

setwd("F:/data/Vaishnavi/combined - 6 frame translation/txt")

user <- Sys.info()[["user"]]

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

evid <- maxquant.read(
    path = ".",
    name = "evidence.txt",
    integer64 = "double")

# 
evid.filt <- evid
evid.filt$group <- "Known"
evid.filt[grep(pattern = "^seq_\\d+(;seq_\\d+)*$", x = evid.filt$Proteins), "group"] <- "Novel"
evid.filt[evid.filt$Reverse == "+", "group"] <- "Reverse"

# 
pdf(file = "Potential_novel_feature.pdf", width = 12, height = 10)

# Calculate quantile values of number of protein per groups field
data.filt <- evid.filt %>%
    dplyr::select(., group, PEP) %>%
    base::as.data.frame(., stringsAsFactors = FALSE)

# Create the boxplot using ggplot
textsize <- 15
ggplot(data = data.filt, mapping = aes(x = factor(group), y = PEP)) +
    geom_boxplot() +
    xlab(label = "Evidence group") +
    ylab(label = "PEP") +
    ggtitle(label = "PEP per evidence category") +
    theme(
        legend.position = "bottom",
        title = element_text(
            face = "bold",
            size = (textsize * 1.25)),
        text = element_text(size = textsize),
        plot.title = element_text(
            face = "bold",
            size = (textsize * 1.5)),
        axis.text.x = element_text(
            angle = 0, vjust = 1, hjust = 0))

# 
ggplot(
    data = data.filt,
    mapping = aes(x = PEP, fill = factor(group), colour = factor(group))) +
    geom_density(alpha = 0.1) +
    xlab(label = "Evidence PEP") +
    ylab(label = "Density") +
    ggtitle(label = "Density of evidence PEP") +
    xlim(0, 0.05)

# 
data <- table(evid.filt$group) %>%
    data.frame(., type = "evidence", stringsAsFactors = FALSE) %>%
    set_names(c("group", "Freq", "type"))

# 
data <- evid.filt %>%
    dplyr::select(., Sequence, group) %>%
    dplyr::group_by(., group) %>%
    dplyr::summarise(., Freq = n_distinct(Sequence)) %>%
    data.frame(., type = "peptide", stringsAsFactors = FALSE) %>%
    rbind(data, .)

# 
data <- evid.filt %>%
    dplyr::select(., `Protein group IDs`, group) %>%
    cSplit(
        indt = ., splitCols = "Protein group IDs",
        sep = ";", direction = "long") %>%
    dplyr::group_by(., group) %>%
    dplyr::summarise(., Freq = n_distinct(`Protein group IDs`)) %>%
    data.frame(., type = "protein group", stringsAsFactors = FALSE) %>%
    rbind(data, .)

tmp <- evid.filt %>%
    dplyr::select(., `Protein group IDs`, group) %>%
    cSplit(
        indt = ., splitCols = "Protein group IDs",
        sep = ";", direction = "long") %>%
    unique(.) %>%
    dplyr::count(`Protein group IDs`) %>%
    base::as.data.frame(., stringsAsFactors = FALSE)

# 
textsize <- 15
ggplot(
    data = data,
    mapping = aes(
        x = type, y = Freq, group = group, fill = group),
    environment = .GlobalEnv) +
    geom_bar(
        stat = "identity", position = "dodge",
        colour = "black") +
    geom_text(
        mapping = aes(label = Freq),
        position = position_dodge(width = 0.9),
        vjust = -0.25, size = (textsize * 0.25)) +
    ggtitle("Features identification") +
    xlab(label = "Feature type") +
    ylab(label = "Feature count (log10)") +
    theme(
        legend.position = "bottom",
        plot.background = element_rect(fill = "#f2f2f2"),
        panel.background = element_rect(fill = "white"),
        panel.grid.major = element_line(colour = "grey40"),
        title = element_text(
            face = "bold",
            size = (textsize * 1.25)),
        text = element_text(size = textsize),
        plot.title = element_text(
            face = "bold",
            size = (textsize * 1.5)),
        axis.text.x = element_text(
            angle = -45, vjust = 1, hjust = 0)) +
    scale_y_log10() +
    scale_fill_grey(start = 0.2, end = 0.8)

# Calculate quantile values of number of protein per groups field
data <- evid.filt %>%
    dplyr::select(., group, `Mass Error [ppm]`) %>%
    base::as.data.frame(., stringsAsFactors = FALSE)

# Create the boxplot using ggplot
textsize <- 15
ggplot(data = data, mapping = aes(x = factor(group), y = `Mass Error [ppm]`)) +
    geom_boxplot() +
    xlab(label = "Evidence group") +
    ylab(label = "Mass Error (ppm)") +
    ggtitle(label = "Mass Error per evidence category") +
    theme(
        legend.position = "bottom",
        title = element_text(
            face = "bold",
            size = (textsize * 1.25)),
        text = element_text(size = textsize),
        plot.title = element_text(
            face = "bold",
            size = (textsize * 1.5)),
        axis.text.x = element_text(
            angle = 0, vjust = 1, hjust = 0))

# 
ggplot(
    data = data,
    mapping = aes(
        x = `Mass Error [ppm]`,
        fill = factor(group),
        colour = factor(group))) +
    xlab(label = "Evidence mass error (ppm)") +
    ylab(label = "Density") +
    ggtitle(label = "Density of evidence mass error") +
    geom_density(alpha = 0.1)

# Calculate quantile values of number of protein per groups field
data <- evid.filt %>%
    dplyr::select(., group, Sequence, `Protein group IDs`) %>%
    cSplit(
        indt = ., splitCols = "Protein group IDs",
        sep = ";", direction = "long") %>%
    dplyr::group_by(., group, `Protein group IDs`) %>%
    dplyr::summarise(., n(), n_distinct(Sequence)) %>%
    tidyr::gather(
        data = .,
        key = "type",
        value = "count",
        `n()`,
        `n_distinct(Sequence)`) %>%
    base::as.data.frame(., stringsAsFactors = FALSE)
data[data$type == "n()", "type"] <- "evidence"
data[data$type == "n_distinct(Sequence)", "type"] <- "peptide"

# Create the boxplot using ggplot
textsize <- 15
ggplot(
    data = data,
    mapping = aes(
        x = factor(type),
        y = count,
        fill = factor(group))) +
    geom_boxplot() +
    xlab(label = "Feature type") +
    ylab(label = "Number of feature per protein group") +
    ggtitle(label = "Feature group per protein group (zoom -1 to 5000)") +
    theme(
        legend.position = "bottom",
        title = element_text(
            face = "bold",
            size = (textsize * 1.25)),
        text = element_text(size = textsize),
        plot.title = element_text(
            face = "bold",
            size = (textsize * 1.5)),
        axis.text.x = element_text(
            angle = 0, vjust = 1, hjust = 0)) +
    ylim(-1, 5000)

ggplot(
    data = data,
    mapping = aes(
        x = factor(type),
        y = count,
        fill = factor(group))) +
    geom_boxplot() +
    xlab(label = "Feature type") +
    ylab(label = "Number of feature per protein group") +
    ggtitle(label = "Feature group per protein group (zoom -1 to 50)") +
    theme(
        legend.position = "bottom",
        title = element_text(
            face = "bold",
            size = (textsize * 1.25)),
        text = element_text(size = textsize),
        plot.title = element_text(
            face = "bold",
            size = (textsize * 1.5)),
        axis.text.x = element_text(
            angle = 0, vjust = 1, hjust = 0)) +
    ylim(-1, 50)

# 
ggplot(
    data = data[data$type == "evidence", ],
    mapping = aes(
        x = `count`,
        fill = factor(group),
        colour = factor(group))) +
    geom_density(alpha = 0.1) +
    xlab(label = "Evidence count per protein group") +
    ylab(label = "Density") +
    ggtitle(label = "Density of evidence per protein group") +
    xlim(-1, 100)

# 
ggplot(
    data = data[data$type == "peptide", ],
    mapping = aes(
        x = `count`,
        fill = factor(group),
        colour = factor(group))) +
    geom_density(alpha = 0.1) +
    xlab(label = "Peptide count per protein group") +
    ylab(label = "Density") +
    ggtitle(label = "Density of peptide per protein group") +
    xlim(-1, 50)

#
tmp <- evid.filt %>%
    cSplit(
        indt = ., splitCols = "Protein group IDs",
        sep = ";", direction = "long") %>%
    base::merge(
        x = .,
        y = data %>%
            tidyr::spread(data = ., key = type, value = count) %>%
            dplyr::filter(., group == "Novel") %>%
            dplyr::select(., -group),
        by = "Protein group IDs",
        all.x = TRUE) %>%
    base::as.data.frame(., stringsAsFactors = FALSE)

#
evid.filt %>%
    dplyr::select(., group, PEP, peptide) %>%
    dplyr::group_by(., group) %>%
    dplyr::summarise(., median(PEP), quantile(PEP, 0.25)) %>%
    base::as.data.frame(., stringsAsFactors = FALSE)

tmp %>%
    dplyr::select(., group, `Protein group IDs`, peptide) %>%
    dplyr::filter(., group == "Novel") %>%
    unique(.) %>%
    dplyr::summarise(., median(peptide), quantile(peptide, 0.75), quantile(peptide, 0.9))

#
known.med.pep <- evid.filt %>%
    dplyr::select(., Sequence, PEP) %>%
    unique(.) %>%
    dplyr::summarise(., median(PEP)) %>%
    as.numeric(.)
prot <- tmp %>%
    dplyr::filter(., group == "Novel" & PEP <= known.med.pep & peptide >= 3) %>%
    dplyr::select(., `Protein group IDs`) %>%
    unique(.)

# 
pg <- maxquant.read(path = ".", name = "proteinGroups.txt")

# 
pg.filt <- pg %>%
    dplyr::filter(., id %in% prot$`Protein group IDs`) %>%
    base::as.data.frame(., stringsAsFactors = FALSE)

# 
pg.filt[, c(1:10,133:134)]
tmp[tmp$`Protein group IDs` == 3636 & tmp$group == "Novel", "Sequence"] %>%
    unique(.)

# 
pg.filt$novelReason <- NA
pg.filt[1, "novelReason"] <- "SAV"
pg.filt[2, "novelReason"] <- "SAV"
pg.filt[3, "novelReason"] <- "SAV"
pg.filt[4, "novelReason"] <- "Uncharacterised protein (L8ADC1) in B. subtilis BEST7613"
pg.filt[5, "novelReason"] <- "Novel start site for C0SPC1"
pg.filt[6, "novelReason"] <- "Novel start site for P50730"
pg.filt[7, "novelReason"] <- "Uncharacterized protein (A0A125UG60) B. sp. LM 4-2"
pg.filt[8, "novelReason"] <- "Novel start site for O31451"
pg.filt[9, "novelReason"] <- "Novel start site for P23449"
pg.filt[10, "novelReason"] <- "Novel start site for O31771"
pg.filt[11, "novelReason"] <- "SAV"
pg.filt[12, "novelReason"] <- "SAV"

#
grid.newpage()
grid.table(
    pg.filt[, c(
        "Protein IDs", "Peptide counts (all)",
        "Q-value", "novelReason")])

#
dev.off()

#
write.table(
    x = pg.filt,
    file = "Potential_novel.txt",
    quote = FALSE,
    sep = "\t",
    row.names = FALSE,
    col.names = TRUE)

out <- evid.filt %>%
    dplyr::filter(., PEP <= known.med.pep, group == "Novel") %>%
    base::as.data.frame(., stringsAsFactors = FALSE)
length(unique(out$Sequence))

out.format <- out %>%
    

write.table(
    x = out,
    file = "Peptide_novel_PEP_knownMedian.txt",
    quote = FALSE,
    sep = "\t",
    row.names = FALSE,
    col.names = TRUE)



### Novel peptide identification -----------------------------------------

# 
fasta.file <- c(
    "F:/data/Vaishnavi/Databases/Bsu_genome_assembly_GCA_000009045.1.out_FIXED_HEADER.fasta",
    "F:/data/Vaishnavi/Databases/uniprot-proteome_Bacillus_subtilis_168_UP000001570_20150318.fasta",
    "G:/MaxQuant/MaxQuant_1.5.1.0/bin/conf/contaminants.fasta"
)

# 
fasta.6frame <- read.fasta(file = fasta.file[1], seqtype = "AA", as.string = TRUE)
fasta.ref <- read.fasta(file = fasta.file[2], seqtype = "AA", as.string = TRUE)
fasta.cont <- read.fasta(file = fasta.file[3], seqtype = "AA", as.string = TRUE)


#
data <- base::data.frame(
    Sequence = unique(evid$Sequence),
    group = NA_character_,
    stringsAsFactors = FALSE)

# 
na.val <- 0
for (x in 1:nrow(data)) {
    if (length(grep(pattern = data$Sequence[x], x = fasta.ref)) > 0) {
        data$group[x] <- "Known"
    } else if (length(grep(pattern = data$Sequence[x], x = fasta.cont)) > 0) {
        data$group[x] <- "Contaminant"
    } else if (length(grep(pattern = data$Sequence[x], x = fasta.6frame)) > 0) {
        data$group[x] <- "Novel"
    } else {
        na.val <- na.val + 1
    }
}

# 
print(paste(
    "There are", na.val, " NA values, these need to be checked!", sep = " "))


#data <- rbind(data, data.orig[!is.na(data.orig$group), ])

#
saveRDS(object = data, file = "Sequence_group_mapping.RDS")

# 
evid.match <- base::merge(x = evid, y = data, by = "Sequence", all = TRUE)

# 
evid.match[evid.match$Reverse == "+", "group"] <- "Reverse"

#
saveRDS(object = evid.match, file = "evid_match.RDS")

tmp <- evid.match[evid.match$group == "Novel", "Sequence"] %>% unique(.)

# 
leven.data <- base::data.frame()
for (x in 1:length(tmp)) {
    
    # 
    dist.name <- adist(
        x = tmp[x],
        y = fasta.ref,
        partial = TRUE,
        ignore.case = TRUE) %>%
        t(.) %>%
        base::data.frame(
            id = rownames(.),
            leven = .,
            Sequence = tmp[x],
            stringsAsFactors = FALSE)
    
    # 
    leven.data <- rbind(
        leven.data,
        dist.name[dist.name$leven == min(dist.name$leven), ])
    
}


