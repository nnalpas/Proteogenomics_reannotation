

my_clan_f <- "D:/Local_databases/X-databases/PFAM/clan.txt/clan.txt"
clan_col <- c(
    "clan_acc",
    "clan_id",
    "previous_id",
    "clan_description",
    "clan_author",
    "deposited_by",
    "clan_comment",
    "updated",
    "created",
    "version",
    "number_structures",
    "number_archs",
    "number_species",
    "number_sequences",
    "competed",
    "uniprot_competed"
)
my_clan <- data.table::fread(
    input = my_clan_f, sep = "\t", quote = "", header = FALSE,
    stringsAsFactors = FALSE, col.names = clan_col)

my_pfam_clan_f <- "D:/Local_databases/X-databases/PFAM/Pfam-A.clans.tsv/Pfam-A.clans.tsv"
pfam_clan_col <- c(
    "pfam_acc",
    "clan_acc",
    "clan_id",
    "pfam_id",
    "description"
)
my_pfam_clan <- data.table::fread(
    input = my_pfam_clan_f, sep = "\t", quote = "", header = FALSE,
    stringsAsFactors = FALSE, col.names = pfam_clan_col)

my_pfam_f <- "D:/Local_databases/X-databases/PFAM/pfamA.txt/pfamA.txt"
pfam_col <- c(
    "pfam_acc",
    "pfam_id",
    "previous_id",
    "description",
    "deposited_by",
    "seed_source",
    "type",
    "comment",
    "sequence_GA",
    "domain_GA",
    "sequence_TC",
    "domain_TC",
    "sequence_NC",
    "domain_NC",
    "buildMethod",
    "model_length",
    "searchMethod",
    "msv_lambda",
    "msv_mu",
    "viterbi_lambda",
    "viterbi_mu",
    "forward_lambda",
    "forward_tau",
    "num_seed",
    "num_full",
    "updated",
    "created",
    "version",
    "number_archs",
    "number_species",
    "number_structures",
    "average_length",
    "percentage_id",
    "average_coverage",
    "change_status",
    "seed_consensus",
    "full_consensus",
    "number_shuffled_hits",
    "number_uniprot",
    "rp_seed",
    "number_rp15",
    "number_rp35",
    "number_rp55",
    "number_rp75"
)
my_pfam <- data.table::fread(
    input = my_pfam_f, sep = "\t", quote = "", header = FALSE,
    stringsAsFactors = FALSE, col.names = pfam_col)

my_pfam_final <- dplyr::full_join(
    x = my_pfam, y = my_pfam_clan,
    by = c("pfam_acc", "pfam_id", "description")) %>%
    dplyr::full_join(
        x = ., y = my_clan, by = c("clan_acc", "clan_id"),
        suffix = c("_pfam", "_clan"))

data.table::fwrite(
    x = my_pfam_final,
    file = "D:/Local_databases/X-databases/PFAM/pfam_and_clan_merge.txt",
    append = FALSE, quote = FALSE, sep = "\t",
    row.names = FALSE, col.names = TRUE)


