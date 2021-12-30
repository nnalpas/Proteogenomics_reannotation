


library(magrittr)
library(ggplot2)

my_plots <- list()
my_date <- Sys.Date()

my_dir <- paste0("H:/data/Synechocystis_6frame/", my_date, "_Genomic_coverage")
dir.create(my_dir)

# Define the current user
user <- Sys.info()[["user"]]

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

#my_f <- "H:/data/Synechocystis_6frame/NoveltyExplain/ORF_novelty_reason_valid.RDS"
my_ggr_f <- "H:/data/Synechocystis_6frame/GRanges/Genome_grange.RDS"
my_rgr_f <- "H:/data/Synechocystis_6frame/GRanges/Ref_prot_grange.RDS"
my_ogr_f <- "H:/data/Synechocystis_6frame/GRanges/Orf_prot_grange.RDS"
my_pgr_f <- "H:/data/Synechocystis_6frame/GRanges/Peptides_grange.RDS"
my_evid_f <- "H:/data/Synechocystis_6frame/Novel_res/Group_evidence.RDS" # this parameter can be removed if MS/MS count is included in the peptide grange

my_ggr <- readRDS(my_ggr_f)
my_rgr <- readRDS(my_rgr_f)
my_ogr <- readRDS(my_ogr_f)
my_pgr <- readRDS(my_pgr_f)
my_evid <- readRDS(my_evid_f)

# Compute and visualise the genomic coverage based on nucleotide coverage
my_plots[["rect_venn"]] <- plots_rectvenn(
    ideo = my_ggr, ref = my_rgr, pep = my_pgr,
    colour = c("#387eb8", "#d1d2d4", "#e21e25", "#fbaf3f", "#404040"))

# Generate the histogram frequency of MS/MS count per nucleotide
my_plots[["msms_count"]] <- plot_nuc_coverage(pep = my_pgr, count = my_evid)

my_plots[["circos"]] <- plot_circos(
    ideo = my_ggr, ref = my_rgr, orf = my_ogr, pep = my_pgr)

pdf(paste(my_dir, "Ref_Novel_coverage.pdf", sep = "/"), 8, 8)
my_plots
dev.off()


