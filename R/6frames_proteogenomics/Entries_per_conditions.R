


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
source(
    file = paste(
        "C:/Users",
        user,
        "Documents/GitHub/Miscellaneous/R/General/General_function.R",
        sep = "/"))

# Load the required packages (or install if not already in library)
load_package("plyr")
load_package("dplyr")
load_package("tidyr")
load_package("seqinr")
load_package("UniProt.ws")
load_package("magrittr")
load_package("WriteXLS")
load_package("data.table")
load_package("splitstackshape")
load_package("VennDiagram")
load_package("ggplot2")
load_package("grid")
load_package("gridExtra")
load_package("RColorBrewer")
load_package("stringr")
load_package("Biostrings")
load_package("RecordLinkage")
load_package("VariantAnnotation")
load_package("cgdsr")
load_package("bit64")
load_package("cleaver")
load_package("plotly")
load_package("GenomicRanges")
load_package("biovizBase")
load_package("ggbio")
load_package("ggradar")
load_package("BSgenome.Bsubtilis.EMBL.AL0091263")
load_package("GenomicFeatures")


# Import the experimental design (with conditions)
exp_design <- read.table(
    file = "E:/Processing/Nicolas/Vaishnavi/experimentalDesignTemplate - Conditions.txt",
    header = TRUE, sep = "\t", quote = "", as.is = TRUE)

evid <- mq_read(
    path = "E:/processing/Nicolas/Vaishnavi/",
    name = "evidence.txt")

prot <- mq_read(
    path = "E:/processing/Nicolas/Vaishnavi/",
    name = "proteinGroups.txt")

acet <- mq_read(
    path = "E:/processing/Nicolas/Vaishnavi/",
    name = "Acetyl \\(K\\)Sites_nonredundant.txt")

phospho <- mq_read(
    path = "E:/processing/Nicolas/Vaishnavi/",
    name = "Phospho \\(STY\\)Sites_nonredundant.txt")


all_prot <- strsplit(x = prot$`Evidence IDs`, split = ";", fixed = TRUE) %>%
    unlist(.) %>%
    unique(.) %>%
    as.integer(.)

all_acet <- strsplit(x = acet$`Evidence.IDs`, split = ";", fixed = TRUE) %>%
    unlist(.) %>%
    unique(.) %>%
    as.integer(.)

all_phospho <- strsplit(x = phospho$`Evidence.IDs`, split = ";", fixed = TRUE) %>%
    unlist(.) %>%
    unique(.) %>%
    as.integer(.)



all_filters <- list(
    "All" = all_prot,
    "Acetyl" = all_acet,
    "Phospho" =  all_phospho)

for (filters in c(1:length(all_filters))) {
    
    data_format <- evid %>%
        dplyr::filter(., id %in% all_filters[[filters]]) %>%
        dplyr::select(
            ., Proteins, `Raw file`, id,
            `Protein group IDs`) %>%
        unique(.) %>%
        cSplit(
            indt = ., splitCols = "Proteins", sep = ";",
            direction = "long", fixed = TRUE) %>%
        unique(.) %>%
        cSplit(
            indt = ., splitCols = "Protein group IDs", sep = ";",
            direction = "long", fixed = TRUE) %>%
        unique(.) %>%
        dplyr::left_join(
            x = ., y = exp_design[, c("Name", "Conditions")],
            by = c("Raw file" = "Name")) %>%
        dplyr::select(., -`Raw file`, -id) %>%
        unique(.)
    
    data_format_cat <- data_format %>%
        dplyr::group_by(., Proteins) %>%
        dplyr::summarise_all(funs(toString(x = unique(.), width = NULL)))
    
    write.table(
        x = data_format,
        file = paste(
            names(all_filters)[filters],
            "Proteins_per_conditions.txt", sep = "_"),
        quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
    
    write.table(
        x = data_format_cat,
        file = paste(
            names(all_filters)[filters],
            "Proteins_per_conditions_compiled.txt", sep = "_"),
        quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
    
}

for (filters in c(2:length(all_filters))) {
    
    data_format <- evid %>%
        dplyr::filter(., id %in% all_filters[[filters]]) %>%
        dplyr::select(
            ., `Modified sequence`, Proteins, `Raw file`, id,
            `Protein group IDs`,
            matches(paste0(names(all_filters)[filters], ".*IDs"))) %>%
        set_colnames(c(
            "Modified sequence", "Proteins", "Raw file",
            "id", "Protein group IDs", "site IDs")) %>%
        unique(.) %>%
        cSplit(
            indt = ., splitCols = "site IDs", sep = ";",
            direction = "long", fixed = TRUE) %>%
        unique(.) %>%
        dplyr::left_join(
            x = ., y = exp_design[, c("Name", "Conditions")],
            by = c("Raw file" = "Name")) %>%
        dplyr::select(., -`Raw file`, -id) %>%
        unique(.) %>%
        set_colnames(c(
            "Modified sequence", "Proteins", "Protein group IDs",
            paste(names(all_filters)[filters], "site IDs"), "Conditions"))
    
    data_format_cat <- data_format %>%
        dplyr::group_by(., `Modified sequence`) %>%
        dplyr::summarise_all(funs(toString(x = unique(.), width = NULL)))
        
    
    write.table(
        x = data_format,
        file = paste(
            names(all_filters)[filters],
            "Sites_per_conditions.txt", sep = "_"),
        quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
    
    write.table(
        x = data_format_cat,
        file = paste(
            names(all_filters)[filters],
            "Sites_per_conditions_compiled.txt", sep = "_"),
        quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
    
}



