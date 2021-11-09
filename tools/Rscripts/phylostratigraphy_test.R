


library(magrittr)
library(phylostratr)

focal_id <- "1148"

strata_uni <- uniprot_strata(taxid = focal_id, from = 1)

strata_app <- strata_apply(
    strata = strata_uni, f = diverse_subtree,
    n = 20, weights = uniprot_weight_by_ref())
#strata_app %<>%
#    use_recommended_prokaryotes(.)


strat_fill <- uniprot_fill_strata(strata_app)
strat_fill <- add_taxa(x = strat_fill, taxa = focal_id)
strat_fill@data$faa[[focal_id]] <- "uniprot-seqs/1148.faa"


strata_bla <- strata_blast(strata = strat_fill)

strata_best <- merge_besthits(best_hits(strata_bla))

strata_final <- stratify(strata_best)


