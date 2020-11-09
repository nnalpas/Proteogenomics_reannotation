


library(magrittr)
library(ggplot2)
library(data.table)

my_files <- c(
    ref = "H:/data/Synechocystis_6frame/Genome/Synechocystis_sp_PCC_6803_cds_aa.fasta")

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
    X = my_proteotypic, FUN = function(x) {
        plot_proteotypic(
            data = x, pep_seq = "Sequence",
            prot = "Proteins", stack = "Filter")
})

my_proteotypic_df <- plyr::ldply(
    my_proteotypic, "data.frame", .id = "Protease")

my_coverage <- protease_coverage(
    data = my_proteotypic_df, protease = "Protease")

my_coverage_per_prot <- protease_coverage(
    data = my_proteotypic_df, protease = "Protease", protein = "Proteins")


