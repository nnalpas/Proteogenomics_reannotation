

my_folder <- "H:/data/Synechocystis_6frame/MaxQuant/combined/txt/"
my_fasta_file <- "H:/data/Synechocystis_6frame/Genome/Synechocystis_sp_PCC_6803_cds_aa.fasta"
my_novel_file <- "H:/data/Synechocystis_6frame/NoveltyExplain/ORF_novelty_reason.RDS"
pep_class



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

my_plots <- list()



### Import data ----------------------------------------------------------

my_pg <- data.table::fread(
    input = paste0(my_folder, "proteinGroups.txt"),
    sep = "\t", quote = "", header = TRUE,
    stringsAsFactors = FALSE, colClasses = "character")

my_fasta <- seqinr::read.fasta(
    file = my_fasta_file, seqtype = "AA", as.string = TRUE)

my_novel <- readRDS(my_novel_file)



### 

my_pg_filt <- my_pg %>%
    dplyr::filter(., Reverse != "+", `Potential contaminant` != "+")

id_prot <- my_pg_filt %>%
    strsplit(x = .[["Protein IDs"]], split = ";", fixed = TRUE) %>%
    unlist(.) %>%
    unique(.)

my_identification <- data.table::data.table(
    `Protein IDs` = names(my_fasta)) %>%
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


