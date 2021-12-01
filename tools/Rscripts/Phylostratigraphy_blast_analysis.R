


library(magrittr)
library(phylostratr)
library(ggplot2)

focal_id <- "1148"
my_phylo_blast_folder <- "H:/data/Synechocystis_6frame/Phylostratigraphy/"

strata_uni <- uniprot_strata(taxid = focal_id, from = 1)

lineages <- lapply(
    1:nleafs(strata_uni@tree), lineage, type = "index", x = strata_uni@tree)

weights <- uniprot_weight_by_ref(clade = 2)
weights <- weights[strata_uni@tree$tip.label]
weights[is.na(weights)] <- 1

chosen <- FUN(lineages, 100, weights, strata_uni@tree)

my_blast_f <- list.files(
    path = my_phylo_blast_folder, pattern = "*_annot", full.names = TRUE) %>%
    set_names(sub("_annot", "", basename(.)))

chosen <- names(my_blast_f)

to_include <- chosen[!chosen %in% strata_uni@tree$tip.label]

strat_fill <- add_taxa(x = strata_uni, taxa = to_include)

if (any(!chosen %in% strat_fill@tree$tip.label)) {
    warning("Still missing taxid!")
}

strata_species <- subtree(
    x = strat_fill, id = chosen, type = "name")

# This part is important to remove all descendents of the focal id
strata_species %<>%
    prune(
        x = .,
        id = descendents(x = strata_app, id = "1148", type = "name"),
        type = "name")
strata_species %<>%
    add_taxa(x = ., taxa = focal_id)

# This code is equivalent to 'strata_best <- strata_besthits(strata_bla)'
strata_species@data$besthit <- lapply(my_blast_f[
    names(my_blast_f) %in% strata_species@tree$tip.label], function(x) {
    data.table::fread(
        input = x, sep = "\t", quote = "", header = TRUE,
        stringsAsFactors = FALSE,
        select = c(
            qseqid= "character", sseqid = "character",
            qstart = "integer", qend = "integer",
            sstart = "integer", send = "integer",
            evalue = "double", score = "double", staxid = "character")) %>%
        dplyr::group_by(., qseqid, staxid) %>%
        dplyr::filter(., score == max(score)) %>%
        dplyr::ungroup(.) %>%
        dplyr::distinct(., qseqid, staxid, .keep_all = TRUE)
})

#tmp <- strata_species
#strata_fold(tmp, function(s) {
#    do.call(what = rbind, s@data$besthit)})

strata_best <- merge_besthits(strata_species)

strata_final <- stratify(strata_best)

# Plot results
plot_heatmaps(
    hits = strata_best, filename = "Heatmap_phylostratr.pdf",
    tree = strata_species@tree, to_name = TRUE, focal_id = focal_id)
make_obo_pdf(
    d = strata_best, file = "Obo_phylostratr.pdf")




save.image("session.RData")


