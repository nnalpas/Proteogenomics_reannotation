


library(magrittr)
library(phylostratr)
library(ggplot2)

focal_id <- "1148"

strata_uni <- uniprot_strata(taxid = focal_id, from = 1)

strata_app <- strata_apply(
    strata = strata_uni, f = diverse_subtree,
    n = 200, weights = uniprot_weight_by_ref(clade = 2))
#strata_app %<>%
#    use_recommended_prokaryotes(.)


strat_fill <- uniprot_fill_strata(strata_app)
#strat_fill <- add_taxa(x = strat_fill, taxa = focal_id)
#strat_fill@data$faa[[focal_id]] <- "uniprot-seqs/1148.faa"


strata_bla <- strata_blast(strata = strat_fill)

strata_bla_best <- strata_besthits(strata_bla)
strata_best_merge <- merge_besthits(strata_bla_best)

strata_final <- stratify(strata_best_merge)






### Code to retrieve all proteome from NCBI ------------------------------


target_taxon <- strata_app@tree$tip.label


tree_statistics <- lapply(X = 1:8, function(x) {
    my_res <- treeio::tree_subset(
        tree = strata_app@tree, node = "1148", levels_back = x)
    data.frame(
        Count = length(my_res$tip.label),
        Tips = paste0(my_res$tip.label, collapse = ";"))
}) %>%
    plyr::ldply(., "data.frame") %>%
    set_colnames(c("Count", "Tips")) %>%
    dplyr::mutate(
        ., NodeLabel = rev(unique(strata_final$mrca_name)[1:8]),
        ps = rev(unique(strata_final$ps)[1:8])) %>%
    #dplyr::rowwise(.) %>%
    dplyr::mutate(
        ., CountDiverse = (
            Count - dplyr::lag(Count, default = 1))) %>%
    dplyr::ungroup(.)
tree_statistics


tree_statistics_split <- tree_statistics %>%
    tidyr::separate_rows(data = ., Tips, sep = ";", convert = TRUE) %>%
    dplyr::group_by(., Tips) %>%
    dplyr::filter(., ps == max(ps)) %>%
    dplyr::ungroup(.)


my_needed_taxa <- tree_statistics_split %>%
    dplyr::select(., Tips, NodeLabel) %>%
    dplyr::mutate(., Category = "Phylostratr")



tomislav_taxon_f <- "H:/data/Synechocystis_6frame/Futo-2020_msaa217_supplementary_data/Supplementary_file_S8.txt"
tomislav_taxon <- data.table::fread(input = tomislav_taxon_f, sep = "\t", quote = "", header = TRUE, stringsAsFactors = FALSE)
tomislav_taxon_filt <- tomislav_taxon %>%
    dplyr::filter(., `Phylostratum (ps)` %in% c(1:2))


table(tomislav_taxon_filt$Tax_ID %in% target_taxon)
View(tomislav_taxon_filt[!tomislav_taxon_filt$Tax_ID %in% target_taxon, ])


tomislav_taxon_final <- tomislav_taxon_filt %>%
    dplyr::filter(., !Tax_ID %in% target_taxon)


my_needed_taxa <- tomislav_taxon_filt %>%
    dplyr::select(., Tips = Tax_ID, NodeLabel = `Phylostratum name`) %>%
    dplyr::mutate(., Category = "Tomislav") %>%
    dplyr::bind_rows(my_needed_taxa, .)



require(XML)
data <- xmlParse("H:/data/Synechocystis_6frame/Phylostratigraphy/taxonomy_result_Cyanobacteria-Melainabacteria.xml")

xml_data <- xmlToList(data)

my_taxa <- lapply(1:length(xml_data), function(x) {
    my_res <- data.frame()
    if (all(c("TaxId", "ScientificName", "Lineage", "Rank") %in% names(xml_data[[x]]))) {
        my_res <- data.frame(
            TaxID = as.integer(xml_data[[x]][["TaxId"]]),
            ScientificName = xml_data[[x]][["ScientificName"]],
            Lineage = xml_data[[x]][["Lineage"]],
            Rank = xml_data[[x]][["Rank"]]
        )
    }
    my_res
}) %>%
    plyr::ldply(., dplyr::bind_rows, .id = NULL) %>%
    dplyr::mutate(., NodeLabel = dplyr::case_when(
        grepl("unclassified Synechocystis(;|$)", Lineage) ~ "unclassified Synechocystis",
        grepl("Merismopediaceae(;|$)", Lineage) ~ "Merismopediaceae",
        grepl("Synechococcales(;|$)", Lineage) ~ "Synechococcales",
        grepl("Cyanobacteria(;|$)", Lineage) ~ "Cyanobacteria",
        grepl("Cyanobacteria\\/Melainabacteria group(;|$)", Lineage) ~ "Cyanobacteria/Melainabacteria group",
        #grepl("Terrabacteria group(;|$)", Lineage) ~ "Terrabacteria group",
        #grepl("Bacteria(;|$)", Lineage) ~ "Bacteria",
        #grepl("cellular organisms(;|$)", Lineage) ~ "cellular organisms",
        TRUE ~ NA_character_
    )) %>%
    dplyr::filter(., !is.na(NodeLabel))


my_needed_taxa <- my_taxa %>%
    dplyr::select(., Tips = TaxID, NodeLabel) %>%
    dplyr::mutate(., Category = "NCBI_taxonomy") %>%
    dplyr::filter(
        ., !Tips %in% my_needed_taxa$Tips &
            NodeLabel %in% tree_statistics[
                tree_statistics$CountDiverse < 200, ][["NodeLabel"]]) %>%
    dplyr::bind_rows(my_needed_taxa, .)

my_needed_taxa <- my_taxa %>%
    dplyr::filter(., grepl(
        "PCC 6803|Gloeobacteria", Lineage, ignore.case = T)) %>%
    dplyr::select(., Tips = TaxID, NodeLabel) %>%
    dplyr::filter(., !Tips %in% my_needed_taxa$Tips) %>%
    dplyr::mutate(., Category = "NCBI_taxonomy") %>%
    dplyr::bind_rows(my_needed_taxa, .)



my_ncbi_refseq_f <- "H:/data/Synechocystis_6frame/Phylostratigraphy/assembly_summary_refseq.txt"
my_ncbi_genbank_f <- "H:/data/Synechocystis_6frame/Phylostratigraphy/assembly_summary_genbank.txt"
my_ncbi_f <- c(
    refseq = my_ncbi_refseq_f,
    genbank = my_ncbi_genbank_f
)

my_ncbi_list <- my_ncbi_f %>%
    lapply(., function(x) {
        data.table::fread(
            input = x, sep = "\t", quote = "",
            header = TRUE, stringsAsFactors = FALSE)
}) %>%
    plyr::ldply(., dplyr::bind_rows, .id = "Repository")


my_ncbi_list_filt <- my_ncbi_list %>%
    dplyr::inner_join(
        x = ., y = my_needed_taxa, by = c("species_taxid" = "Tips"))
my_ncbi_list_filt <- my_ncbi_list %>%
    dplyr::inner_join(
        x = ., y = my_needed_taxa, by = c("taxid" = "Tips")) %>%
    dplyr::bind_rows(my_ncbi_list_filt, my_ncbi_list_filt_tmp) %>%
    unique(.)

taxa_with_genome_count <- length(unique(my_ncbi_list_filt$species_taxid))

my_ncbi_list_filt$Repository <- factor(
    x = my_ncbi_list_filt$Repository,
    levels = c("refseq", "genbank"),
    ordered = TRUE)
my_ncbi_list_filt$refseq_category <- factor(
    x = my_ncbi_list_filt$refseq_category,
    levels = c("reference genome", "representative genome", "na"),
    ordered = TRUE)
my_ncbi_list_filt$assembly_level <- factor(
    x = my_ncbi_list_filt$assembly_level,
    levels = c("Complete Genome", "Chromosome", "Scaffold", "Contig"),
    ordered = TRUE)
my_ncbi_list_filt$genome_rep <- factor(
    x = my_ncbi_list_filt$genome_rep,
    levels = c("Full", "Partial"),
    ordered = TRUE)
my_ncbi_list_filt$Category <- factor(
    x = my_ncbi_list_filt$Category,
    levels = c("Phylostratr", "Tomislav", "NCBI_taxonomy"),
    ordered = TRUE)
my_ncbi_list_filt$seq_rel_date <- as.Date(x = my_ncbi_list_filt$seq_rel_date)

inclusion_reasons <- c(
    "", "partial", "missing strain identifier",
    "derived from single cell", "derived from environmental source",
    "derived from environmental source; derived from metagenome",
    "derived from metagenome", "missing tRNA genes",
    "derived from metagenome; not used as type",
    "unverified source organism; derived from metagenome",
    "derived from metagenome; derived from single cell",
    "derived from environmental source; partial",
    "unverified source organism; partial"
)

my_ncbi_list_final <- my_ncbi_list_filt %>%
    dplyr::group_by(., species_taxid) %>%
    dplyr::filter(., Repository == min(Repository)) %>%
    dplyr::filter(., refseq_category == min(refseq_category)) %>%
    dplyr::filter(., assembly_level == min(assembly_level)) %>%
    dplyr::filter(., genome_rep == min(genome_rep)) %>%
    dplyr::filter(., Category == min(Category)) %>%
    dplyr::filter(., excluded_from_refseq %in% inclusion_reasons) %>%
    dplyr::group_by(
        ., taxid, `species_taxid`, `organism_name`,
        `infraspecific_name`, isolate) %>%
    dplyr::filter(., seq_rel_date == max(seq_rel_date)) %>%
    dplyr::mutate(., DuplicationIsolate = dplyr::n()) %>%
    dplyr::group_by(., species_taxid) %>%
    dplyr::mutate(., DuplicationSpecies = dplyr::n()) %>%
    dplyr::ungroup(.) %>%
    dplyr::arrange(., species_taxid)

taxa_with_one_genome_count <- length(unique(my_ncbi_list_final$species_taxid))

if (taxa_with_genome_count != taxa_with_one_genome_count) {
    warning("Lost some taxa during genome filtering!")
}

table(my_ncbi_list_final$NodeLabel)

data.table::fwrite(
    x = my_ncbi_list_final,
    file = "H:/data/Synechocystis_6frame/Phylostratigraphy/Selected_representative_genome.txt",
    append = FALSE, quote = FALSE, sep = "\t",
    row.names = FALSE, col.names = TRUE)

my_genome_dir <- "T:/User/Nicolas/Phylostratigraphy"

my_genome_files <- list.files(
    path = my_genome_dir, pattern = "*.faa", full.names = TRUE) %>%
    data.table::data.table(path = .) %>%
    dplyr::mutate(
        ., assembly_accession = sub("^(.+?_.+?)_.+", "\\1", basename(path))) %>%
    dplyr::left_join(x = ., y = my_ncbi_list_final) %>%
    dplyr::group_by(., taxid) %>%
    dplyr::mutate(., Duplication = dplyr::n())

#data.table::fwrite(
#    x = my_genome_files,
#    file = "H:/data/Synechocystis_6frame/Phylostratigraphy/Downloaded_representative_genome.txt",
#    append = FALSE, quote = FALSE, sep = "\t",
#    row.names = FALSE, col.names = TRUE)

# uncompress and concatenate all proteins per taxon id (check for unique IDs)
out_dir <- "T:/User/Nicolas/Phylostratigraphy/Formatted"
dir.create(path = out_dir)
options(warn = 0)
for (t in unique(my_genome_files$taxid)) {
    
    t_fasta_f <- c(my_genome_files[my_genome_files$taxid == t, ][["path"]])
    
    if (length(t_fasta_f) == 1) {
        file.copy(from = t_fasta_f, to = paste0(out_dir, "/", t, ".faa"))
    } else {
        print(paste("Concatenating:", t, "..."))
        my_seqs <- lapply(t_fasta_f, function(f) {
            seqinr::read.fasta(
                file = f, seqtype = "AA",
                as.string = TRUE, whole.header = TRUE)
        }) %>%
            unlist(., recursive = FALSE)
        my_seqs_filt <- my_seqs[!duplicated(my_seqs)]
        seqinr::write.fasta(
            sequences = my_seqs_filt, names = names(my_seqs_filt),
            file.out = paste0(out_dir, "/", t, ".faa"),
            open = "w", nbchar = 80, as.string = TRUE)
        
    }
    
}

toplot <- my_genome_files %>%
    dplyr::select(., taxid, NodeLabel) %>%
    unique(.) %>%
    dplyr::group_by(., NodeLabel) %>%
    dplyr::summarise(., Count = dplyr::n())
toplot$NodeLabel <- factor(
    x = toplot$NodeLabel,
    levels = levels(tree_statistics$NodeLabel),
    ordered = TRUE)

ggplot(toplot, aes(x = NodeLabel, y = Count, label = Count)) +
    geom_bar(stat = "identity", position = "dodge") +
    geom_text(position = position_dodge(width = 0.9), hjust = -0.1) +
    ggpubr::theme_pubr() +
    coord_flip()

# create the serial blast phenodata
my_serial_blast <- my_genome_files %>%
    dplyr::ungroup(.) %>%
    dplyr::select(., taxid) %>%
    unique(.) %>%
    dplyr::mutate(
        ., Task = "blastp",
        Query = paste0("${PBS_O_HOME}/work/Synechocystis_6frame/Phylostratigraphy/", focal_id, ".faa"),
        Subject = paste0("${PBS_O_HOME}/work/Synechocystis_6frame/Phylostratigraphy/", taxid, ".faa"),
        Title = taxid, Selection = "all", threshold = 10, maxhits = 250,
        Param1 = "", Param2 = "") %>%
    dplyr::select(., -taxid)

data.table::fwrite(
    x = my_serial_blast,
    file = "T:/User/Nicolas/Phylostratigraphy/Phylostratigraphy_blasts.txt",
    append = FALSE, quote = FALSE, sep = "\t",
    row.names = FALSE, col.names = FALSE)

save.image("H:/data/Synechocystis_6frame/Phylostratigraphy/phylo_2021-11-25.RData")


