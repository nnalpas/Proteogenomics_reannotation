

library(magrittr)

my_ncbi_refseq_f <- "H:/data/Synechocystis_6frame/Phylostratigraphy/assembly_summary_refseq.txt"
my_ncbi_genbank_f <- "H:/data/Synechocystis_6frame/Phylostratigraphy/assembly_summary_genbank.txt"
my_ncbi_f <- c(
    refseq = my_ncbi_refseq_f,
    genbank = my_ncbi_genbank_f
)

my_needed_taxa <- data.table::fread(
    input = "C:/Users/kxmna01/Desktop/Pae_taxonomy_result.txt",
    sep = "\t", quote = "", header = FALSE, stringsAsFactors = FALSE,
    col.names = "taxid")

my_ncbi_list <- my_ncbi_f %>%
    lapply(., function(x) {
        data.table::fread(
            input = x, sep = "\t", quote = "",
            header = TRUE, stringsAsFactors = FALSE)
    }) %>%
    plyr::ldply(., dplyr::bind_rows, .id = "Repository")

my_ncbi_list_filt <- my_ncbi_list %>%
    dplyr::filter(., species_taxid %in% my_needed_taxa$taxid | taxid %in% my_needed_taxa$taxid)

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
    dplyr::group_by(., taxid) %>%
    dplyr::filter(., Repository == min(Repository)) %>%
    dplyr::filter(., refseq_category == min(refseq_category)) %>%
    dplyr::filter(., assembly_level == min(assembly_level)) %>%
    dplyr::filter(., genome_rep == min(genome_rep)) %>%
    dplyr::filter(., excluded_from_refseq %in% inclusion_reasons) %>%
    dplyr::group_by(
        ., taxid, `species_taxid`, `organism_name`,
        `infraspecific_name`, isolate) %>%
    dplyr::filter(., seq_rel_date == max(seq_rel_date)) %>%
    dplyr::mutate(., DuplicationIsolate = dplyr::n()) %>%
    dplyr::ungroup(.) %>%
    dplyr::arrange(., taxid)

nrow(my_ncbi_list_final)

my_ncbi_list_final %<>%
    dplyr::filter(., assembly_level %in% c("Complete Genome", "Chromosome"))

nrow(my_ncbi_list_final)

data.table::fwrite(
    x = my_ncbi_list_final,
    file = "H:/data/Pathogens_6frame/Genome/Selected_representative_genome.txt",
    append = FALSE, quote = FALSE, sep = "\t",
    row.names = FALSE, col.names = TRUE)






