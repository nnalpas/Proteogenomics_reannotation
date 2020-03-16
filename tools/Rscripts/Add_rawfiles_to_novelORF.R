



library(magrittr)
library(dplyr)
library(data.table)

my_orfs <- fread(input = "C:/Users/kxmna01/Dropbox/Home_work_sync/Work/Srim_6frame/NoveltyExplain/ORF_novelty_reason.txt", sep = "\t", quote = "", header = TRUE, stringsAsFactors = FALSE, data.table = FALSE)

my_evid <- fread(input = "C:/Users/kxmna01/Dropbox/Home_work_sync/Work/Srim_6frame/Novel_res/Group_evidence.txt", sep = "\t", quote = "", header = TRUE, colClasses = "character", stringsAsFactors = FALSE, data.table = FALSE)

my_orfs_format <- my_orfs %>%
    dplyr::select(., Proteins, Novel_Sequence) %>%
    tidyr::separate_rows(data = ., Novel_Sequence, sep = ";") %>%
    dplyr::left_join(
        x = .,
        y = my_evid[, c(
            "Sequence", "Raw file", "Experiment",
            grep("Intensity ", colnames(my_evid), value = TRUE))],
        by = c("Novel_Sequence" = "Sequence")) %>%
    tidyr::pivot_longer(
        data = ., cols = grep("Intensity ", colnames(.), value = TRUE),
        names_to = "key", values_to = "Intensity") %>%
    dplyr::filter(., Intensity > 0) %>%
    dplyr::mutate(
        ., key = sub("Intensity ", "", key)) %>%
    dplyr::group_by(., Proteins, `Raw file`, Experiment) %>%
    dplyr::summarise(., key = paste(unique(key), collapse = ", ")) %>%
    dplyr::group_by(., Proteins) %>%
    dplyr::mutate(.,
        `Raw file` = paste0(`Raw file`, " (", key, ")")) %>%
    dplyr::select(., -key) %>%
    dplyr::group_by(., Proteins) %>%
    dplyr::summarise_all(~paste(unique(.), collapse = "; "))

my_orfs_final <- dplyr::left_join(
    x = my_orfs, y = my_orfs_format, by = "Proteins")

fwrite(x = my_orfs_final, file = "C:/Users/kxmna01/Dropbox/Home_work_sync/Work/Srim_6frame/NoveltyExplain/ORF_novelty_reason_RAWfiles.txt", append = FALSE, quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)


