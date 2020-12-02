


my_folder <- "H:/data/Synechocystis_6frame/MaxQuant/combined/txt/"
my_fasta_file <- "H:/data/Synechocystis_6frame/Genome/Synechocystis_sp_PCC_6803_cds_aa.fasta"
my_novel_file <- "H:/data/Synechocystis_6frame/NoveltyExplain/ORF_novelty_reason.RDS"
my_path_file <- "D:/Local_databases/X-databases/KEGG/2020-12-01_syn_kegg_pathways.txt"
my_tu_file <- "H:/data/Synechocystis_6frame/Kopf_2014_TU/Transcriptional_units_formatted_20201022.txt"
pep_class <- c("class 1", "class 2")



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

# Load the required packages
library(magrittr)
library(ggplot2)
library(dendextend)

my_plots <- list()



### Import data ----------------------------------------------------------

my_pg <- data.table::fread(
    input = paste0(my_folder, "proteinGroups.txt"),
    sep = "\t", quote = "", header = TRUE,
    stringsAsFactors = FALSE, colClasses = "character")

my_fasta <- seqinr::read.fasta(
    file = my_fasta_file, seqtype = "AA",
    as.string = TRUE, whole.header = TRUE)

my_novel <- readRDS(my_novel_file)

my_path <- data.table::fread(
    input = my_path_file,
    sep = "\t", quote = "", header = TRUE,
    stringsAsFactors = FALSE, colClasses = "character")

my_tu <- data.table::fread(
    input = my_tu_file,
    sep = "\t", header = TRUE,
    stringsAsFactors = FALSE, colClasses = "character")



### Data formatting ------------------------------------------------------

my_pg_filt <- my_pg %>%
    dplyr::filter(., Reverse != "+", `Potential contaminant` != "+")

my_path_format <- my_path %>%
    dplyr::mutate(
        ., KeggID = sub("syn:", "", KeggID),
        KeggID = sub("syn:", "", KeggID),
        UniProtID = sub("up:", "", UniProtID),
        PathwayID = sub("path:", "", PathwayID),
        PathwayName = sub(" - Synechocystis sp. PCC 6803", "", PathwayName)) %>%
    dplyr::group_by(., KeggID) %>%
    dplyr::summarise_all(~paste0(unique(na.omit(.)), collapse = "|")) %>%
    dplyr::ungroup(.) %>%
    dplyr::mutate(., PathwayName = ifelse(PathwayName == "", NA, PathwayName))

id_prot <- my_pg_filt %>%
    strsplit(x = .[["Protein IDs"]], split = ";", fixed = TRUE) %>%
    unlist(.) %>%
    unique(.)

my_identification <- data.table::data.table(
    Header = names(my_fasta)) %>%
    tidyr::separate(
        data = ., col = "Header",
        into = c("Protein IDs", "Genome", "Header_remain"),
        sep = " \\| ", extra = "merge") %>%
    tidyr::separate(
        data = ., col = "Header_remain",
        into = c("Description", "Location"),
        sep = " \\| ", fill = "left") %>%
    dplyr::mutate(
        ., Location = sub("join\\{(.+)\\}", "\\1", Location) %>%
            sub("\\:.+\\:", ":", .) %>%
            gsub("\\[|\\)", "", .) %>%
            sub("](", ":", ., fixed = TRUE)) %>%
    tidyr::separate(
        data = ., col = "Location",
        into = c("Start", "End", "Strand"),
        sep = "\\:", convert = TRUE) %>%
    dplyr::mutate(
        ., `Identification` = ifelse(
            !`Protein IDs` %in% id_prot, "Never identified", "Identified"),
        Location = round(x = Start, digits = -4)) %>%
    tidyr::unite(
        data = ., col = "Genome_loc", Genome, Location,
        sep = "_", remove = FALSE)

my_identification$Genome_loc <- factor(
    x = my_identification$Genome_loc,
    levels = unique(my_identification$Genome_loc),
    ordered = TRUE)

my_tu_format <- my_tu %>%
    tidyr::separate_rows(data = ., Gene_name, sep = ", ") %>%
    dplyr::select(
        ., TU_id = id, `Protein IDs` = Gene_name,
        TU_description = Gene_description, TU_type) %>%
    dplyr::group_by(., `Protein IDs`) %>%
    dplyr::summarise_all(~paste0(unique(na.omit(.)), collapse = "|")) %>%
    dplyr::filter(., !is.na(`Protein IDs`) & `Protein IDs` != "")



### Over-representation analysis -----------------------------------------

foreground <- my_identification[
    my_identification$Identification == "Never identified", ][["Protein IDs"]]

background <- my_identification[["Protein IDs"]]

my_oa_path <- clusterProfiler::enrichKEGG(
    gene = foreground, organism = "syn", keyType = "kegg",
    pAdjustMethod = "BH", universe = background,
    pvalueCutoff = 0.5, qvalueCutoff = 0.5)

my_oa_path_df <- my_oa_path@result %>%
    as.data.frame(.)

my_oa_mod <- clusterProfiler::enrichMKEGG(
    gene = foreground, organism = "syn", keyType = "kegg",
    pAdjustMethod = "BH", universe = background,
    minGSSize = 3,
    pvalueCutoff = 0.5, qvalueCutoff = 0.5)

my_oa_mod_df <- my_oa_mod@result %>%
    as.data.frame(.)



### Additional characterisation through histogram ------------------------

my_identification_annot <- my_identification %>%
    dplyr::left_join(
        x = ., y = my_path_format, by = c("Protein IDs" = "KeggID")) %>%
    dplyr::left_join(
        x = ., y = my_tu_format)

for (x in c("Description", "Genome", "Genome_loc", "TU_id")) {
    my_plots[[paste0(x, "_count")]] <- describe_count_plot(
        my_data = my_identification_annot, main_id = "Protein IDs",
        x = x, fill = "Identification", separator = "\\|",
        gtitle = x, flip = TRUE,
        ylabel = "Protein count")
}

for (x in c("PathwayName", "TU_type")) {
    my_plots[[paste0(x, "_count")]] <- describe_count_plot(
        my_data = my_identification_annot, main_id = "Protein IDs",
        x = x, fill = "Identification", separator = "\\|",
        gtitle = x, flip = TRUE, threshold = 20,
        ylabel = "Protein count")
}




