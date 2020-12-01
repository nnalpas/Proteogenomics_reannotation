

my_folder <- "H:/data/Synechocystis_6frame/MaxQuant/combined/txt/"
my_fasta_file <- "H:/data/Synechocystis_6frame/Genome/Synechocystis_sp_PCC_6803_cds_aa.fasta"
my_novel_file <- "H:/data/Synechocystis_6frame/NoveltyExplain/ORF_novelty_reason.RDS"
my_path_file <- "D:/Local_databases/X-databases/KEGG/2020-12-01_syn_kegg_pathways.txt"
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



### 

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
    dplyr::mutate(., `Identification` = ifelse(
        !`Protein IDs` %in% id_prot, "Never identified", "Identified"))

my_pg_format <- my_pg_filt %>%
    dplyr::select(
        ., `Protein IDs`, `Majority protein IDs`,
        dplyr::starts_with("iBAQ")) %>%
    tidyr::pivot_longer(
        data = ., cols = dplyr::starts_with("iBAQ"),
        names_to = "Experiment", values_to = "iBAQ") %>%
    dplyr::mutate(
        ., Label = sub(";.+", "", `Majority protein IDs`),
        Experiment = sub("^iBAQ ", "", Experiment),
        iBAQ = as.double(iBAQ)) %>%
    dplyr::filter(., !grepl("^peptides", Experiment)) %>%
    dplyr::group_by(., Experiment) %>%
    dplyr::mutate(
        ., iBAQ_sum = sum(iBAQ, na.rm = T),
        iBAQ_perc = iBAQ / iBAQ_sum * 100,
        iBAQ_log10 = log10(iBAQ + 1)) %>%
    dplyr::ungroup(.)

my_plots[["iBAQ_raw"]] <- ggplot(
    my_pg_format, aes(x = Experiment, y = iBAQ_log10)) +
    geom_boxplot() +
    ggpubr::theme_pubr() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

my_plots[["iBAQ_perc"]] <- ggplot(
    my_pg_format, aes(x = Experiment, y = iBAQ_perc)) +
    geom_boxplot() +
    ggpubr::theme_pubr() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
    scale_y_log10()

my_plots[["iBAQ_heat"]] <- ggplot(
    my_pg_format, aes(x = Experiment, y = Label, fill = iBAQ_log10)) +
    geom_tile() +
    ggpubr::theme_pubr() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
    scale_fill_gradient(low = "#17115C", high = "#C20808", na.value = "grey")

my_plots[["iBAQ_perc_heat"]] <- ggplot(
    my_pg_format, aes(x = Experiment, y = Label, fill = iBAQ_perc)) +
    geom_tile() +
    ggpubr::theme_pubr() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
    scale_fill_gradient(low = "#17115C", high = "#C20808", na.value = "grey")



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

my_identification_annot <- my_identification %>%
    dplyr::left_join(
        x = ., y = my_path_format, by = c("Protein IDs" = "KeggID"))

for (x in c("Description")) {
    my_plots[[paste0(x, "_count")]] <- describe_count_plot(
        my_data = my_identification_annot, main_id = "Protein IDs",
        x = x, fill = "Identification", separator = "\\|",
        gtitle = x, flip = TRUE,
        ylabel = "Protein count")
}

for (x in c("PathwayName")) {
    my_plots[[paste0(x, "_count")]] <- describe_count_plot(
        my_data = my_identification_annot, main_id = "Protein IDs",
        x = x, fill = "Identification", separator = "\\|",
        gtitle = x, flip = TRUE, threshold = 20,
        ylabel = "Protein count")
}



