

# Start with clean environment
rm(list = ls())

# Define parameters for script
wkdir <- "F:/data/Vaishnavi/combined - SILAC processing/txt/"
id.col <- "Protein IDs"
selec <- c("Protein IDs", "Q-value", "Intensity")
path.list <- c(
    "F:/data/Vaishnavi/combined - SILAC processing/txt/",
    "P:/Vaishnavi/220raw_05122016/combined_old/txt")
names <- c("Full processing", "Partial processing")
RevCon <- TRUE
ids <- "uniprot"

# Define workspace location
setwd(wkdir)

# Get current date
date.str <- Sys.time() %>%
    format(., "%Y-%m-%d")

# Define the user
user <- Sys.info()[["user"]]

# Source the mandatory custom functions
source(
    file = paste(
        "C:/Users/",
        user,
        "/Documents/GitHub/Miscellaneous/R/General/General_function.R",
        sep = ""))

# Load the required packages (or install if not already in library)
# Note: several functions are conflicting between packages, make sure to use
# the package name when using a function
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



### Look into proteins ---------------------------------------------------

# Create pdf for reporting
pdf(
    file = paste("Comparison_full_vs_220_", date.str, ".pdf", sep = ""),
    width = 10, height = 10)

# Create a dataframe containing all samples required data columns
names(path.list) <- names
selec %<>%
    lapply(X = ., FUN = as.symbol)
data <- data.frame()
for (x in names(path.list)) {
    
    # Read in and select required data
    tmp <- maxquant.read(
            path = path.list[[x]],
            name = "proteinGroups.txt")
    
    # Check whether reverse and contaminant should be filtered out
    if (RevCon) {
        tmp %<>%
            rev.con.filt(data = .)
    }
    
    # Keep only required columns
    tmp %<>%
        dplyr::select_(., .dots = selec) %>%
        base::cbind(., file = x, stringsAsFactors = FALSE)
    data <- rbind(data, tmp)
    rm(tmp)
    
}

# Check that ids are from uniprot and if yes clean ids
if (ids == "uniprot") {
    
    data[[id.col]] <- uniprotID.clean(x = data[[id.col]])
    
}

# Venn of protein group that match
vennPlots(
    data = data,
    key = "file",
    value = id.col,
    main = "Comparison of protein IDs")

# Correlation plot of Q-values
corrPlots(
    data = data %>% dplyr::select(., -Intensity),
    key = "file",
    value = "Q-value",
    main = "Q-value correlation")

# Boxplot of Q-values
boxPlots(
    data = data,
    key = "file",
    value = "Q-value",
    main = "Q-value boxplot",
    zoom = c(-0.0001, 0.0001))

# Correlation plot of Intensity
corrPlots(
    data = data %>%
        dplyr::select(., -`Q-value`),
    key = "file",
    value = "Intensity",
    main = "Intensity correlation")

# Boxplot of Intensity
boxPlots(
    data = data,
    key = "file",
    value = "Intensity",
    main = "Intensity boxplot",
    zoom = c(1E8, 1E10),
    transf = "log10")

# Get table of number of entries pper processing
nmb.table <- data %>%
    dplyr::group_by(., file) %>%
    dplyr::summarise(., n_distinct(`Protein IDs`)) %>%
    base::data.frame(., type = "protein", stringsAsFactors = FALSE) %>%
    set_colnames(c("Processing", "Number features", "Feature type"))



### Look into peptides ---------------------------------------------------

# Define new parameters
id.col <- "Sequence"
selec <- c("Sequence", "PEP", "Score", "Intensity")
ids <- ""

# Create a dataframe containing all samples required data columns
names(path.list) <- names
selec %<>%
    lapply(X = ., FUN = as.symbol)
data <- data.frame()
for (x in names(path.list)) {
    
    # Read in and select required data
    tmp <- maxquant.read(
        path = path.list[[x]],
        name = "peptides.txt",
        integer64 = "double")
    
    # Check whether reverse and contaminant should be filtered out
    if (RevCon) {
        tmp %<>%
            rev.con.filt(data = .)
    }
    
    # Keep only required columns
    tmp %<>%
        dplyr::select_(., .dots = selec) %>%
        base::cbind(., file = x, stringsAsFactors = FALSE)
    data <- rbind(data, tmp)
    rm(tmp)
    
}

# Check that ids are from uniprot and if yes clean ids
if (ids == "uniprot") {
    
    data[[id.col]] <- uniprotID.clean(x = data[[id.col]])
    
}

# Get table of number of entries pper processing
nmb.table <- data %>%
    dplyr::group_by(., file) %>%
    dplyr::summarise(., n_distinct(Sequence)) %>%
    base::data.frame(., type = "peptide", stringsAsFactors = FALSE) %>%
    set_colnames(c("Processing", "Number features", "Feature type")) %>%
    rbind(nmb.table, .)

# Histogram of number of feature
histPlots(
    data = nmb.table,
    key = "Feature type",
    value = "Number features",
    group = "Processing",
    fill = "Processing",
    main = "Number of feature per processing")

# Venn of peptide sequence that match
vennPlots(
    data = data,
    key = "file",
    value = id.col,
    main = "Comparison of peptide sequence")

# Correlation plot of PEP
corrPlots(
    data = data %>% dplyr::select(., -Intensity, -Score),
    key = "file",
    value = "PEP",
    main = "PEP correlation",
    zoom = c(-0.01, 0.25))

# Boxplot of PEP
boxPlots(
    data = data,
    key = "file",
    value = "PEP",
    main = "PEP boxplot",
    zoom = c(-0.0001, 0.0001))

# Correlation plot of Score
corrPlots(
    data = data %>% dplyr::select(., -Intensity, -PEP),
    key = "file",
    value = "Score",
    main = "Score correlation")

# Boxplot of Score
boxPlots(
    data = data,
    key = "file",
    value = "Score",
    main = "Score boxplot")

# Correlation plot of Intensity
corrPlots(
    data = data %>% dplyr::select(., -Score, -PEP),
    key = "file",
    value = "Intensity",
    main = "Intensity correlation",
    transf = "log10")

# Boxplot of Intensity
boxPlots(
    data = data,
    key = "file",
    value = "Intensity",
    main = "Intensity boxplot",
    zoom = c(5E5, 5E6),
    transf = "log10")

# Close the report
dev.off()


