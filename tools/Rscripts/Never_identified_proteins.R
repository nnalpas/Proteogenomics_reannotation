#!/usr/bin/env Rscript

# This script determines the proteins that have never been identified
# and attempts to characterise them to provide a list of new experimental
# conditions to try...



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



### Parameters setting up ------------------------------------------------

#
opt <- list(
    my_folder = "H:/data/Synechocystis_6frame/MaxQuant/combined/txt/",
    my_fasta_file = "H:/data/Synechocystis_6frame/Genome/Synechocystis_sp_PCC_6803_cds_aa.fasta",
    my_novel_file = "H:/data/Synechocystis_6frame/NoveltyExplain/ORF_novelty_reason.RDS",
    my_path_file = "D:/Local_databases/X-databases/KEGG/2020-12-01_syn_kegg_pathways.txt",
    my_tu_file = "H:/data/Synechocystis_6frame/Kopf_2014_TU/Transcriptional_units_formatted_20201022.txt",
    #pep_class = c("class 1", "class 2"),
    max_oa = 10,
    output = "C:/Users/kxmna01/Desktop/2020-12-02_Never_identified"
)

# Define input parameters (interactively or from command line)
#if (interactive()) {
#    
#    # Define the list of input parameters
#    opt <- list(
#        novel_reason = choose.files(
#            caption = "Choose peptide novelty (.RDS) file!",
#            multi = FALSE),
#        ref_grange = choose.files(
#            caption = "Choose reference entries grange (.RDS) file!",
#            multi = FALSE),
#        orf_grange = choose.files(
#            caption = "Choose ORF entries grange (.RDS) file!",
#            multi = FALSE),
#        operon_grange = choose.files(
#            caption = "Choose operon entries grange (.RDS) file!",
#            multi = FALSE),
#        pep_pos = choose.files(
#            caption = "Choose peptide position within protein (.RDS) file!",
#            multi = FALSE),
#        pep_grange = choose.files(
#            caption = "Choose peptide entries grange (.RDS) file!",
#            multi = FALSE),
#        sanger_grange = choose.files(
#            caption = "Choose Sanger seq entries grange (.RDS) file!",
#            multi = FALSE),
#        genome_grange = choose.files(
#            caption = "Choose genome entry grange (.RDS) file!",
#            multi = FALSE),
#        add_rbs = readline(prompt = paste0(
#            "Provide additional RBS sequence",
#            " (separated by space)?")) %>%
#            as.character(),
#        pep_class = readline(prompt = paste0(
#            "Which PEP class to keep",
#            " (separated by space; e.g. 'class 1')")) %>%
#            as.character(),
#        bsgenome = readline(prompt = paste0(
#            "Provide name of the BSgenome package")) %>%
#            as.character(),
#        threads = readline(prompt = "How many cores to use?") %>% as.integer(),
#        output = readline(
#            prompt = "Define the output directory!"))
#    
#} else {
#    
#    # Define the list of command line parameters
#    option_list <- list(
#        make_option(
#            opt_str = c("-i", "--novel_reason"),
#            type = "character", default = NULL,
#            help = "Peptide novelty info",
#            metavar = "character"),
#        make_option(
#            opt_str = c("-r", "--ref_grange"),
#            type = "character", default = NULL,
#            help = "Reference entries grange file",
#            metavar = "character"),
#        make_option(
#            opt_str = c("-n", "--orf_grange"),
#            type = "character", default = NULL,
#            help = "ORF entries grange file",
#            metavar = "character"),
#        make_option(
#            opt_str = c("-p", "--operon_grange"),
#            type = "character", default = NULL,
#            help = "Operon entries grange file",
#            metavar = "character"),
#        make_option(
#            opt_str = c("-e", "--pep_pos"),
#            type = "character", default = NULL,
#            help = "Peptide position file",
#            metavar = "character"),
#        make_option(
#            opt_str = c("-d", "--pep_grange"),
#            type = "character", default = NULL,
#            help = "Peptide entries grange file",
#            metavar = "character"),
#        make_option(
#            opt_str = c("-s", "--sanger_grange"),
#            type = "character", default = NULL,
#            help = "Sanger seq entries grange file",
#            metavar = "character"),
#        make_option(
#            opt_str = c("-g", "--genome_grange"),
#            type = "character", default = NULL,
#            help = "Genome entry grange file",
#            metavar = "character"),
#        make_option(
#            opt_str = c("-a", "--add_rbs"),
#            type = "character", default = NULL,
#            help = "Additional RBS sequence (separated by space)",
#            metavar = "character"),
#        make_option(
#            opt_str = c("-c", "--pep_class"),
#            type = "character", default = NULL,
#            help = "Which PEP class to keep (separated by space; e.g. 'class 1')",
#            metavar = "character"),
#        make_option(
#            opt_str = c("-b", "--bsgenome"),
#            type = "character", default = NULL,
#            help = "The BSgenome package name for genomic visualisation",
#            metavar = "character"),
#        make_option(
#            opt_str = c("-t", "--threads"),
#            type = "integer", default = NULL,
#            help = "Number of cores to use",
#            metavar = "character"),
#        make_option(
#            opt_str = c("-o", "--output"),
#            type = "character", default = NULL,
#            help = "Output directory", metavar = "character"))
#    
#    # Parse the parameters provided on command line by user
#    opt_parser <- OptionParser(option_list = option_list)
#    opt <- parse_args(opt_parser)
#    
#}

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



### Import data ----------------------------------------------------------

my_pg <- data.table::fread(
    input = paste0(opt$my_folder, "proteinGroups.txt"),
    sep = "\t", quote = "", header = TRUE,
    stringsAsFactors = FALSE, colClasses = "character")

my_fasta <- seqinr::read.fasta(
    file = opt$my_fasta_file, seqtype = "AA",
    as.string = TRUE, whole.header = TRUE)

my_novel <- readRDS(opt$my_novel_file)

my_path <- data.table::fread(
    input = opt$my_path_file,
    sep = "\t", quote = "", header = TRUE,
    stringsAsFactors = FALSE, colClasses = "character")

my_tu <- data.table::fread(
    input = opt$my_tu_file,
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

never_id_prot <- my_identification[
    my_identification$Identification == "Never identified", ][["Protein IDs"]]
my_fasta_filt <- my_fasta[which(
    sub(" .+", "", names(my_fasta)) %in% never_id_prot)]



### Over-representation analysis -----------------------------------------

foreground <- my_identification[
    my_identification$Identification == "Never identified", ][["Protein IDs"]]

background <- my_identification[["Protein IDs"]]

my_oa_path <- clusterProfiler::enrichKEGG(
    gene = foreground, organism = "syn", keyType = "kegg",
    pAdjustMethod = "BH", universe = background,
    pvalueCutoff = 0.5, qvalueCutoff = 0.5)

my_oa_path_df <- my_oa_path@result %>%
    as.data.frame(.) %>%
    dplyr::rowwise(.) %>%
    dplyr::mutate(
        ., GeneRatioNum = eval(parse(text = GeneRatio)),
        padjust_log = -log10(p.adjust))

if (nrow(my_oa_path_df) > 0) {
    my_plots[["OA_path"]] <- ggplot(
        data = my_oa_path_df %>%
            dplyr::arrange(., dplyr::desc(padjust_log)) %>%
            dplyr::slice(., 1:opt$max_oa),
        mapping = aes(
            x = GeneRatioNum, y = padjust_log,
            fill = Count, label = Description)) +
        geom_point(shape = 21, size = 5) +
        ggrepel::geom_text_repel() +
        ggpubr::theme_pubr() +
        geom_hline(
            yintercept = -log10(0.05), colour = "gold",
            linetype = "dashed", size = 1) +
        scale_fill_gradient(low = "#3C3C82", high = "#A31717") +
        xlab("Gene ratio") +
        ylab("-log10 adj. p-value") +
        ggtitle(paste0("KEGG pathways (top ", opt$max_oa, ")"))
}

my_oa_mod <- clusterProfiler::enrichMKEGG(
    gene = foreground, organism = "syn", keyType = "kegg",
    pAdjustMethod = "BH", universe = background,
    minGSSize = 3,
    pvalueCutoff = 0.5, qvalueCutoff = 0.5)

my_oa_mod_df <- my_oa_mod@result %>%
    as.data.frame(.) %>%
    dplyr::rowwise(.) %>%
    dplyr::mutate(
        ., GeneRatioNum = eval(parse(text = GeneRatio)),
        padjust_log = -log10(p.adjust))

if (nrow(my_oa_mod_df) > 0) {
    my_plots[["OA_module"]] <- ggplot(
        data = my_oa_mod_df %>%
            dplyr::arrange(., dplyr::desc(padjust_log)) %>%
            dplyr::slice(., 1:opt$max_oa),
        mapping = aes(
            x = GeneRatioNum, y = padjust_log,
            fill = Count, label = Description)) +
        geom_point(shape = 21, size = 5) +
        ggrepel::geom_text_repel() +
        ggpubr::theme_pubr() +
        geom_hline(
            yintercept = -log10(0.05), colour = "gold",
            linetype = "dashed", size = 1) +
        scale_fill_gradient(low = "#3C3C82", high = "#A31717") +
        xlab("Gene ratio") +
        ylab("-log10 adj. p-value") +
        ggtitle(paste0("KEGG modules (top ", opt$max_oa, ")"))
}



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



### Exporting results ----------------------------------------------------

# Open a file for plot visualisation
pdf(
    file = paste0(opt$output, "/", "Never_identified.pdf"),
    width = 12, height = 12)
my_plots
dev.off()

# Export all tables
write.table(
    x = my_oa_path_df,
    file = paste0(opt$output, "/", "Never_identified_OA_pathways.txt"),
    append = FALSE, quote = FALSE, sep = "\t",
    row.names = FALSE, col.names = TRUE)
write.table(
    x = my_oa_mod_df,
    file = paste0(opt$output, "/", "Never_identified_OA_modules.txt"),
    append = FALSE, quote = FALSE, sep = "\t",
    row.names = FALSE, col.names = TRUE)
write.table(
    x = my_identification_annot,
    file = paste0(opt$output, "/", "Never_identified.txt"),
    append = FALSE, quote = FALSE, sep = "\t",
    row.names = FALSE, col.names = TRUE)

# Export fasta file of never identified proteins
seqinr::write.fasta(
    sequences = my_fasta_filt, names = names(my_fasta_filt),
    file.out = paste0(opt$output, "/", "Never_identified.fasta"),
    open = "w", nbchar = 60, as.string = TRUE)

# Save the session
save.image(paste0(opt$output, "/", "Never_identified.RData"))



### END ------------------------------------------------------------------

# Time scripts end
print(paste("Completed ", format(Sys.time(), "%Y-%m-%d"), sep = ""))


