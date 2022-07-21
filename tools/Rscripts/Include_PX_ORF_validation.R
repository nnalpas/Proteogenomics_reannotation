


rm(list = ls())

library(magrittr)

opt <- list()

opt$novel <- "H:/data/Pathogens_6frame/NoveltyExplain/ORF_novelty_reason.RDS"
opt$novel_valid <- "H:/data/Pathogens_6frame/Novel_res_validation"

orf_reason_final <- readRDS(opt$novel)

my_validation_f <- list.files(
    path = opt$novel_valid, pattern = "Peptides_grange.RDS",
    full.names = TRUE, recursive = TRUE) %>%
    set_names(sub(paste0(opt$novel_valid, "/"), "", dirname(.), fixed = TRUE))

if (length(my_validation_f) > 0) {
    
    pep_grange <- lapply(my_validation_f, readRDS) %>%
        as(., "GRangesList") %>%
        unlist(.)
    
    # Loop through all novel ORFs
    for (x in orf_reason_final$Proteins) {
        pep_gr_filt <- subset(pep_grange, grepl(paste0(x, "(,|$)"), Proteins)) %>%
            subset(., group == "Novel")
        if (length(pep_gr_filt) > 0) {
            orf_reason_final[orf_reason_final$Proteins == x, "PX_datasets"] <- names(pep_gr_filt) %>%
                sub("\\..+", "", .) %>%
                unique(.) %>%
                paste0(., collapse = ";")
            orf_reason_final[orf_reason_final$Proteins == x, "PX_novel_sequence"] <- names(pep_gr_filt) %>%
                sub(".+\\.", "", .) %>%
                unique(.) %>%
                paste0(., collapse = ";")
            orf_reason_final[orf_reason_final$Proteins == x, "PX_novel_peptide_count"] <- names(pep_gr_filt) %>%
                sub(".+\\.", "", .) %>%
                dplyr::n_distinct(.)
            orf_reason_final[orf_reason_final$Proteins == x, "PX_coverage_start"] <- min(
                c(start(pep_gr_filt), end(pep_gr_filt)))
            orf_reason_final[orf_reason_final$Proteins == x, "PX_Coverage_end"] <- max(
                c(start(pep_gr_filt), end(pep_gr_filt)))
        }
    }
    
} else {
    orf_reason_final$PX_datasets <- ""
    orf_reason_final$PX_novel_sequence <- ""
    orf_reason_final$PX_novel_peptide_count <- ""
    orf_reason_final$PX_coverage_start <- ""
    orf_reason_final$PX_Coverage_end <- ""
}

saveRDS(object = orf_reason_final, file = sub(".RDS", "_valid.RDS", opt$novel))

data.table::fwrite(
    x = orf_reason_final, file = sub(".RDS", "_valid.txt", opt$novel),
    append = FALSE, quote = FALSE, sep = "\t",
    row.names = FALSE, col.names = TRUE)


