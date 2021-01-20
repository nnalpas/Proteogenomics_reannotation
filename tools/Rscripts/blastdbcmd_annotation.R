


### Working environment set-up -------------------------------------------

#
rm(list = ls())

#
library(dplyr)
library(tidyr)
library(plyr)
library(magrittr)
library(splitstackshape)
library(data.table)



### Import the data ------------------------------------------------------

#
my_uni <- fread(
    input = "C:\\Users\\kxmna01\\Dropbox\\Home_work_sync\\Work\\Srim_6frame\\ReciprocalBlast\\uniprot_annot.txt",
    sep = "\t", quote = "",
    header = FALSE, stringsAsFactors = FALSE, data.table = FALSE)

my_ref <- fread(
    input = "C:\\Users\\kxmna01\\Dropbox\\Home_work_sync\\Work\\Srim_6frame\\ReciprocalBlast\\uniref_annot.txt",
    sep = "\t", quote = "",
    header = FALSE, stringsAsFactors = FALSE, data.table = FALSE)

my_ncbi <- fread(
    input = "C:\\Users\\kxmna01\\Dropbox\\Home_work_sync\\Work\\Srim_6frame\\ReciprocalBlast\\refseq_annot.txt",
    sep = "\t", quote = "",
    header = FALSE, stringsAsFactors = FALSE, data.table = FALSE)

my_strath <- fread(
    input = "C:\\Users\\kxmna01\\Dropbox\\Home_work_sync\\Work\\Srim_6frame\\ReciprocalBlast\\Strathclyde_annot.txt",
    sep = ";", quote = "",
    header = FALSE, stringsAsFactors = FALSE, data.table = FALSE)

my_best_uni <- fread(
    input = "C:\\Users\\kxmna01\\Dropbox\\Home_work_sync\\Work\\Srim_6frame\\Blast\\Best_blast_ORFprot_vs_Uniprot",
    sep = "\t", quote = "",
    header = TRUE, stringsAsFactors = FALSE, data.table = FALSE)

my_best_ref <- fread(
    input = "C:\\Users\\kxmna01\\Dropbox\\Home_work_sync\\Work\\Srim_6frame\\Blast\\Best_blast_ORFprot_vs_Refprot",
    sep = "\t", quote = "",
    header = TRUE, stringsAsFactors = FALSE, data.table = FALSE)

my_best_ncbi <- fread(
    input = "C:\\Users\\kxmna01\\Dropbox\\Home_work_sync\\Work\\Srim_6frame\\Blast\\Best_blast_ORFprot_vs_NCBIprot",
    sep = "\t", quote = "",
    header = TRUE, stringsAsFactors = FALSE, data.table = FALSE)

my_best_strath <- fread(
    input = "C:\\Users\\kxmna01\\Dropbox\\Home_work_sync\\Work\\Srim_6frame\\Blast\\Best_blast_ORFprot_vs_Strathprot",
    sep = "\t", quote = "",
    header = TRUE, stringsAsFactors = FALSE, data.table = FALSE)

my_best_reci_uni <- fread(
    input = "C:\\Users\\kxmna01\\Dropbox\\Home_work_sync\\Work\\Srim_6frame\\ReciprocalBlast\\Best_Reciproc_Blast_Uniprot_vs_ORFprot",
    sep = "\t", quote = "",
    header = TRUE, stringsAsFactors = FALSE, data.table = FALSE)

my_best_reci_ref <- fread(
    input = "C:\\Users\\kxmna01\\Dropbox\\Home_work_sync\\Work\\Srim_6frame\\ReciprocalBlast\\Best_Reciproc_Blast_Refprot_vs_ORFprot",
    sep = "\t", quote = "",
    header = TRUE, stringsAsFactors = FALSE, data.table = FALSE)

my_best_reci_ncbi <- fread(
    input = "C:\\Users\\kxmna01\\Dropbox\\Home_work_sync\\Work\\Srim_6frame\\ReciprocalBlast\\Best_Reciproc_Blast_NCBIprot_vs_ORFprot",
    sep = "\t", quote = "",
    header = TRUE, stringsAsFactors = FALSE, data.table = FALSE)

my_best_reci_strath <- fread(
    input = "C:\\Users\\kxmna01\\Dropbox\\Home_work_sync\\Work\\Srim_6frame\\ReciprocalBlast\\Best_Reciproc_Blast_Strathprot_vs_ORFprot",
    sep = "\t", quote = "",
    header = TRUE, stringsAsFactors = FALSE, data.table = FALSE)



### Format the functional annotation -------------------------------------

#
my_uni_format <- my_uni %>%
    dplyr::mutate(., `V1` = sub(" OS=", ";", `V1`, fixed = TRUE)) %>%
    dplyr::mutate(., `V1` = sub(" OX=([^ ]+).+", ";\\1", `V1`)) %>%
    tidyr::separate(
        data = ., col = "V1",
        into = c("ID", "Description", "Taxon", "TaxonID"),
        sep = ";")

my_ref_format <- my_ref %>%
    dplyr::mutate(., `V1` = sub(" OS=", ";", `V1`, fixed = TRUE)) %>%
    dplyr::mutate(., `V1` = sub(" OX=([^ ]+).+", ";\\1", `V1`)) %>%
    tidyr::separate(
        data = ., col = "V1",
        into = c("ID", "Description", "Taxon", "TaxonID"),
        sep = ";")

my_ncbi_format <- my_ncbi %>%
    tidyr::separate(
        data = ., col = "V1",
        into = c("ID", "DescriptionTMP", "TaxonID"),
        sep = ";") %>%
    dplyr::mutate(
        ., DescriptionTMP = sub(
            "MULTISPECIES: ", "", DescriptionTMP, fixed = TRUE)) %>%
    dplyr::mutate(
        ., DescriptionTMP = sub("]", "", DescriptionTMP, fixed = TRUE)) %>%
    dplyr::mutate(
        ., DescriptionTMP = sub(" [", ";", DescriptionTMP, fixed = TRUE)) %>%
    tidyr::separate(
        data = ., col = "DescriptionTMP",
        into = c("Description", "Taxon"),
        sep = ";")

my_strath_format <- my_strath %>%
    set_colnames(c("ID", "header", "TaxonID")) %>%
    tidyr::separate(
        data = ., col = "header",
        into = c("Nothin", "GeneName", "Type", "DescriptionTMP", "Location"),
        sep = "\\|") %>%
    dplyr::mutate(
        ., Description = gsub(
            "(^ | $)", "", DescriptionTMP),
        Taxon = "Streptomyces rimosus subsp. rimosus ATCC 10970") %>%
    dplyr::select(., ID, Description, Taxon, TaxonID)



### Annotation for best hits ---------------------------------------------

#
my_best_uni %<>%
    dplyr::mutate(., sseqid = sub(".+\\|(.+?)\\|.+", "\\1", sseqid)) %>%
    dplyr::left_join(x = ., y = my_uni_format, by = c("sseqid" = "ID"))

my_best_ref %<>%
    dplyr::mutate(., sseqid = sub(".+\\|(.+?)\\|.+", "\\1", sseqid)) %>%
    dplyr::left_join(x = ., y = my_ref_format, by = c("sseqid" = "ID"))

my_best_ncbi %<>%
    dplyr::mutate(., sseqid = sub("\\|$", "", sseqid) %>% sub(".+\\|", "", .)) %>%
    dplyr::left_join(x = ., y = my_ncbi_format, by = c("sseqid" = "ID"))

my_best_strath %<>%
    dplyr::left_join(x = ., y = my_strath_format, by = c("sseqid" = "ID"))

#
data.table::fwrite(
    x = my_best_uni,
    file = "C:\\Users\\kxmna01\\Dropbox\\Home_work_sync\\Work\\Srim_6frame\\Blast\\Best_blast_ORFprot_vs_Uniprot_annot",
    append = FALSE, quote = FALSE, sep = "\t",
    row.names = FALSE, col.names = TRUE)

data.table::fwrite(
    x = my_best_ref,
    file = "C:\\Users\\kxmna01\\Dropbox\\Home_work_sync\\Work\\Srim_6frame\\Blast\\Best_blast_ORFprot_vs_Refprot_annot",
    append = FALSE, quote = FALSE, sep = "\t",
    row.names = FALSE, col.names = TRUE)

data.table::fwrite(
    x = my_best_ncbi,
    file = "C:\\Users\\kxmna01\\Dropbox\\Home_work_sync\\Work\\Srim_6frame\\Blast\\Best_blast_ORFprot_vs_NCBIprot_annot",
    append = FALSE, quote = FALSE, sep = "\t",
    row.names = FALSE, col.names = TRUE)

data.table::fwrite(
    x = my_best_strath,
    file = "C:\\Users\\kxmna01\\Dropbox\\Home_work_sync\\Work\\Srim_6frame\\Blast\\Best_blast_ORFprot_vs_Strathprot_annot",
    append = FALSE, quote = FALSE, sep = "\t",
    row.names = FALSE, col.names = TRUE)



### Annotation for reciprocal best hits ----------------------------------

#
my_best_reci_uni %<>%
    dplyr::left_join(x = ., y = my_uni_format, by = c("sseqid" = "ID"))

my_best_reci_ref %<>%
    dplyr::left_join(x = ., y = my_ref_format, by = c("sseqid" = "ID"))

my_best_reci_ncbi %<>%
    dplyr::left_join(x = ., y = my_ncbi_format, by = c("sseqid" = "ID"))

my_best_reci_strath %<>%
    dplyr::left_join(x = ., y = my_strath_format, by = c("sseqid" = "ID"))

#
data.table::fwrite(
    x = my_best_reci_uni,
    file = "C:\\Users\\kxmna01\\Dropbox\\Home_work_sync\\Work\\Srim_6frame\\ReciprocalBlast\\Best_Reciproc_Blast_Uniprot_vs_ORFprot_annot",
    append = FALSE, quote = FALSE, sep = "\t",
    row.names = FALSE, col.names = TRUE)

data.table::fwrite(
    x = my_best_reci_ref,
    file = "C:\\Users\\kxmna01\\Dropbox\\Home_work_sync\\Work\\Srim_6frame\\ReciprocalBlast\\Best_Reciproc_Blast_Refprot_vs_ORFprot_annot",
    append = FALSE, quote = FALSE, sep = "\t",
    row.names = FALSE, col.names = TRUE)

data.table::fwrite(
    x = my_best_reci_ncbi,
    file = "C:\\Users\\kxmna01\\Dropbox\\Home_work_sync\\Work\\Srim_6frame\\ReciprocalBlast\\Best_Reciproc_Blast_NCBIprot_vs_ORFprot_annot",
    append = FALSE, quote = FALSE, sep = "\t",
    row.names = FALSE, col.names = TRUE)

data.table::fwrite(
    x = my_best_reci_strath,
    file = "C:\\Users\\kxmna01\\Dropbox\\Home_work_sync\\Work\\Srim_6frame\\ReciprocalBlast\\Best_Reciproc_Blast_Strathprot_vs_ORFprot_annot",
    append = FALSE, quote = FALSE, sep = "\t",
    row.names = FALSE, col.names = TRUE)


