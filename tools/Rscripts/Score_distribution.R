


library(magrittr)
library(ggplot2)

my_big_cols <- c(
    "#387eb8", "#404040", "#e21e25", "#fbaf3f", "#d1d2d4", "#246E39", "#753B94")

my_hq <- c(
    "chr1_105061", "chr1_121253", "chr1_30235", "chr1_66945", "chr1_282265", "pcb2_4_87", "chr1_143957", "chr1_302301", "chr1_11522", "pca2_4_128", "chr1_32540", "chr1_253927", "chr1_71049", "chr1_303974", "chr1_30", "chr1_92970", "chr1_89481", "chr1_98563", "chr1_92958", "chr1_21768", "chr1_158171", "chr1_8011", "chr1_205934", "chr1_16997", "chr1_214319", "chr1_195865", "chr1_93138", "chr1_16254", "chr1_9090", "chr1_236975", "chr1_16981", "chr1_211568", "chr1_232675", "chr1_168797", "chr1_114270", "chr1_157119", "chr1_15966", "pca2_4_103", "chr1_314381", "psysm_10543", "chr1_181836", "chr1_36421", "pcb2_4_104", "chr1_155520", "chr1_154008", "chr1_144958", "chr1_116115", "chr1_198514", "chr1_15210", "psysm_7517", "chr1_106542", "psysm_9217", "chr1_221462", "psysa_6401", "chr1_66981", "chr1_112570", "chr1_17047", "chr1_309173", "chr1_88209", "chr1_235525", "chr1_86903", "chr1_169827", "chr1_234166", "chr1_146410"
)

my_evid_group_f <- "T:/User/Nicolas/Synechocystis_6frame/Novel_res/Group_evidence.RDS"

my_evid_group <- readRDS(my_evid_group_f)

dim(my_evid_group)

my_evid_format <- my_evid_group %>%
    dplyr::select(., Proteins, group, Score)

my_evid_format %<>%
    dplyr::filter(
        ., grepl(paste0(my_hq, collapse = "|"), my_evid_format$Proteins)) %>%
    dplyr::mutate(., group = "High quality") %>%
    dplyr::bind_rows(my_evid_format, .)

my_evid_format %<>%
    dplyr::filter(., group %in% c("Known", "Novel", "High quality"))

my_evid_format$group <- factor(
    x = my_evid_format$group, levels = c("Known", "Novel", "High quality"), ordered = T)

pl <- ggplot(my_evid_format, aes(x = group, y = Score, fill = group, colour = group)) +
    geom_boxplot(alpha = 0.5) +
    ggpubr::theme_pubr() +
    scale_fill_manual(values = my_big_cols) +
    scale_colour_manual(values = my_big_cols) +
    xlab("Type") +
    ylab("Andromeda score")

pdf("Andromeda_score_distribution.pdf", 8, 8)
pl
dev.off()
    

