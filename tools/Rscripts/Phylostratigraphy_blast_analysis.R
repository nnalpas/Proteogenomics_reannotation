


library(magrittr)
library(phylostratr)
library(ggplot2)
library(ggtree)

focal_id <- "1148"
my_phylo_blast_folder <- "H:/data/Synechocystis_6frame/Phylostratigraphy/"

my_plots <- list()

strata_uni <- uniprot_strata(taxid = focal_id, from = 1)

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

strata_best <- merge_besthits(strata_species)

strata_final <- stratify(strata_best)

# Organise data for overrepresentation analysis
strata_final_oa <- strata_final %>%
    dplyr::mutate(., value = 1, Label = paste0("ps ", ps , " - ", mrca_name)) %>%
    tidyr::pivot_wider(data = ., names_from = Label, values_from = value) %>%
    dplyr::select(., -ps, -mrca, -mrca_name) %>%
    dplyr::mutate_all(~tidyr::replace_na(., replace = 0))

# Plot results
plot_heatmaps(
    hits = strata_best, filename = "Heatmap_phylostratr.pdf",
    tree = strata_species@tree, to_name = TRUE, focal_id = focal_id)

toplot <- strata_final %>%
    dplyr::group_by(., ps, mrca_name) %>%
    dplyr::summarise(., Count = dplyr::n_distinct(qseqid)) %>%
    dplyr::ungroup(.) %>%
    dplyr::mutate(., Label = paste0("ps ", ps , " - ", mrca_name)) %>%
    dplyr::arrange(., dplyr::desc(ps))
toplot$Label <- factor(
    x = toplot$Label, levels = unique(toplot$Label), ordered = TRUE)

my_plots[["ps_count"]] <- ggplot(toplot, aes(x = Label, y = Count, label = Count)) +
    geom_bar(stat = "identity", position = "dodge") +
    geom_text(
        stat = "identity", position = position_dodge(width = 0.9),
        hjust = -0.3) +
    ggpubr::theme_pubr() +
    coord_flip()

toplot_tree <- strata_species@tree
toplot_tree$tip.label <- partial_id_to_name(toplot_tree$tip.label)
toplot_tree$node.label <- partial_id_to_name(toplot_tree$node.label)
toplot_tree$node.label[30:length(toplot_tree$node.label)] <- ""

my_ps_list <- strata_best %>%
    dplyr::select(., staxid, mrca, ps) %>%
    unique(.) %>%
    dplyr::rowwise(.) %>%
    dplyr::mutate(
        ., TaxonName = partial_id_to_name(staxid),
        mrca_name = partial_id_to_name(mrca),
        Label = paste0("ps ", ps , " - ", mrca_name))

toplot <- my_ps_list %>%
    dplyr::group_by(., Label) %>%
    dplyr::summarise(., Count = dplyr::n_distinct(staxid)) %>%
    dplyr::ungroup(.)
toplot$Label <- factor(
    x = toplot$Label, levels = unique(toplot$Label), ordered = TRUE)

my_plots[["ps_org_count"]] <- ggplot(
    toplot, aes(x = Label, y = Count, label = Count)) +
    geom_bar(stat = "identity", position = "dodge") +
    geom_text(
        stat = "identity", position = position_dodge(width = 0.9),
        hjust = -0.3) +
    ggpubr::theme_pubr() +
    coord_flip()

my_group <- split(
    my_ps_list$TaxonName,
    my_ps_list$Label)

toplot_tree_grouped <- groupOTU(toplot_tree, my_group)

my_plots[["phylogeny_tree"]] <- ggtree(
    toplot_tree_grouped, aes(color = group)) +
    geom_text2(aes(label = label), hjust = .3) +
    geom_tiplab(size = 1)

my_plots[["phylogeny_circular"]] <- ggtree::ggtree(
    toplot_tree_grouped, aes(color = group), layout = 'circular') +
    ggtree::geom_tiplab(size = 1, aes(angle = angle))

pdf("Phylostrata_plots.pdf", 10, 10)
my_plots
dev.off()

data.table::fwrite(
    x = strata_final, file = "Phylostrata_proteins.txt",
    append = FALSE, quote = FALSE, sep = "\t",
    row.names = FALSE, col.names = TRUE)

data.table::fwrite(
    x = strata_final_oa, file = "Phylostrata_for_OA.txt",
    append = FALSE, quote = FALSE, sep = "\t",
    row.names = FALSE, col.names = TRUE)

save.image("Session_phylostratr.RData")


