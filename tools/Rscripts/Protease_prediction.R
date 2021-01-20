


library(magrittr)
library(ggplot2)
library(data.table)

my_files <- c(
    ref = "L:/Software/AndromedaDatabases/Pseudomonas_aeruginosa/UP000002438_208964_complete_2019-12-11.fasta")

my_fastas <- lapply(X = my_files, FUN = function(x) {
    seqinr::read.fasta(file = x, seqtype = "AA", as.string = T) %>%
        set_names(sub(".+\\|(.+)\\|.+", "\\1", names(.)))
}) %>%
    set_names(names(my_files))

my_digest <- lapply(X = names(cleaver:::rules), FUN = function(enz) {
    prot_digest(
        fasta = my_fastas$ref, enzyme = enz,
        custom = NULL, miss_cleav = 0)
}) %>%
    set_names(names(cleaver:::rules))

my_proteotypic <- lapply(X = my_digest, FUN = function(x) {
    is_proteotypic(data = x)
})

my_plots <- list()

my_plots[["Overall"]] <- lapply(
    X = names(my_proteotypic), FUN = function(x) {
        plot_proteotypic(
            data = my_proteotypic[[x]], pep_seq = "Sequence",
            prot = "Proteins", stack = "Filter", title = x)
})

my_proteotypic_df <- plyr::ldply(
    my_proteotypic, "data.frame", .id = "Protease")

my_plots[["Coverage"]] <- protease_coverage(
    data = my_proteotypic_df, protease = "Protease")

my_plots[["Coverage_per_prot"]] <- protease_coverage(
    data = my_proteotypic_df, protease = "Protease", protein = "Proteins")

my_plots[["detectable_prot"]] <- plot_covered_protein(
    data = my_proteotypic_df, protease = "Protease", protein = "Proteins")

my_plots[["protease_combination"]] <- plot_protease_combination(
    my_proteotypic_df)

outfile <- basename(my_files) %>%
    sub("\\.fasta$", "", .)

save.image(file = paste0(outfile, ".RData"))

pdf(file = paste0(outfile, "_plots.pdf"), width = 10, height = 10)
my_plots
dev.off()

data.table::fwrite(
    x = my_proteotypic_df, file = paste0(outfile, "_proteotypic.txt"),
    append = FALSE, quote = FALSE, sep = "\t",
    row.names = FALSE, col.names = TRUE)


