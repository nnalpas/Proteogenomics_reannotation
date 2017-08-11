


# Start with clean environment
rm(list = ls())



### Define working directory ---------------------------------------------

# Define the work space
work_space <- choose.dir()
setwd(work_space)

# Define the current user
user <- Sys.info()[["user"]]

# Define current time
date_str <- format(Sys.time(), "%Y-%m-%d")

# Variable for additional filter
filter <- TRUE



### List of required packages -----------------------------------------------

# Source the custom user functions
source(
    file = paste(
        "C:/Users/",
        user,
        "/Documents/GitHub/Miscellaneous/R/General/General_function.R",
        sep = ""))

# Load the required packages (or install if not already in library)
# Note: the select function from dplyr and UniProt.ws is overlapping,
# therefore order of package loading is important
load_package(plyr)
load_package(dplyr)
load_package(tidyr)
load_package(seqinr)
load_package(UniProt.ws)
load_package(magrittr)
load_package(WriteXLS)
load_package(data.table)
load_package(splitstackshape)
load_package(VennDiagram)
load_package(ggplot2)
load_package(grid)
load_package(gridExtra)
load_package(RColorBrewer)
load_package(stringr)
load_package(Biostrings)
load_package(RecordLinkage)
load_package(VariantAnnotation)
load_package(cgdsr)
load_package(bit64)
load_package(cleaver)
load_package(plotly)
load_package(GenomicRanges)
load_package(biovizBase)
load_package(ggbio)
load_package(ggradar)



### Data import ----------------------------------------------------------

# Select the maxquant txt folder
txt_dir <- choose.dir(caption = "Select the MaxQuant txt folder?")

# Import the maxquant evidence table
evid <- readRDS(file = "./2017-08-10_Sequence_group_mapping.RDS")

# Import the maxquant proteingroups table
pg <- mq_read(
    path = txt_dir,
    name = "proteinGroups.txt",
    integer64 = "double")
pg_not_sites <- pg %>%
    dplyr::filter(., `Only identified by site` == "") %>%
    cSplit(indt = ., splitCols = "Evidence IDs", sep = ";", direction = "long")


msms <- mq_read(
    path = txt_dir,
    name = "msms.txt",
    integer64 = "double")


evid_localfdr <- evid %>%
    dplyr::select(., id, PEP, Reverse) %>%
    dplyr::arrange(., PEP) %>%
    dplyr::mutate(., Reverse_val = ifelse(Reverse == "+", 1, 0)) %>%
    dplyr::mutate(., Reverse_cumsum = cumsum(Reverse_val)) %>%
    dplyr::mutate(., rank = base::rank(x = PEP, ties.method = "first")) %>%
    dplyr::mutate(., Local_FDR = (2 * Reverse_cumsum / rank * 100))
evid_localfdr_filt <- evid_localfdr %>%
    dplyr::filter(., Local_FDR <= (max(Local_FDR) / 3.5))




evid_count <- evid %>%
    dplyr::group_by(., group) %>%
    dplyr::summarise(., count = n_distinct(id)) %>%
    dplyr::mutate(., Type = "Without filtering")

evid_count <- evid %>%
    dplyr::filter(., id %in% unique(pg_not_sites$`Evidence IDs`)) %>%
    dplyr::group_by(., group) %>%
    dplyr::summarise(., count = n_distinct(id)) %>%
    dplyr::mutate(., Type = "No id by sites") %>%
    dplyr::bind_rows(evid_count, .)

evid_count <- evid %>%
    dplyr::filter(., id %in% unique(pg_not_sites$`Evidence IDs`)) %>%
    dplyr::filter(., id %in% unique(evid_localfdr_filt$id)) %>%
    dplyr::group_by(., group) %>%
    dplyr::summarise(., count = n_distinct(id)) %>%
    dplyr::mutate(., Type = "Adj. FDR") %>%
    dplyr::bind_rows(evid_count, .)




evid_pep <- evid %>%
    dplyr::select(., group, PEP) %>%
    dplyr::mutate(., Type = "Without filtering")

evid_pep <- evid %>%
    dplyr::filter(., id %in% unique(pg_not_sites$`Evidence IDs`)) %>%
    dplyr::select(., group, PEP) %>%
    dplyr::mutate(., Type = "No id by sites") %>%
    dplyr::bind_rows(evid_pep, .)

evid_pep <- evid %>%
    dplyr::filter(., id %in% unique(pg_not_sites$`Evidence IDs`)) %>%
    dplyr::filter(., id %in% unique(evid_localfdr_filt$id)) %>%
    dplyr::select(., group, PEP) %>%
    dplyr::mutate(., Type = "Adj. FDR") %>%
    dplyr::bind_rows(evid_pep, .)



pl <- plots_hist(
    data = evid_count,
    key = "group",
    value = "count",
    group = "Type",
    fill = "Type",
    main = "Number of evidence",
    transf = "log10",
    label = "count",
    posit = "dodge",
    legend = "bottom")
pl




pl <- plots_box(
    data = evid_pep,
    key = "Type",
    value = "PEP",
    fill = "group",
    main = "PEP distribution",
    zoom = c(-0.001, 0.04))
pl[[1]]
pl[[2]]










pep_table <- mq_read(
    path = "G:/data/Vaishnavi/combined - 6 frame translation/txt",
    name = "peptides.txt")
tmp <- evid %>%
    dplyr::filter(., `MS/MS Scan Number` > 0) %>%
    dplyr::left_join(., pep_table[, c("Sequence", "Start position")]) %>%
    dplyr::mutate(
        ., end_position = `Start position` + Length - 1) %>%
    dplyr::filter(., !is.na(`Start position`))

tmp2 <- apply(X = tmp, MARGIN = 1, FUN = function(x){
    seq(x[["Start position"]], x[["end_position"]], 1) %>%
        paste(., collapse = ";")
})
evid_withpos <- tmp2 %>%
    data.frame(
        Positions = .) %>%
    dplyr::bind_cols(tmp, .)

evid_coverage <- evid_withpos %>%
    dplyr::select(., Sequence, `Leading Proteins`, `MS/MS Count`, Positions) %>%
    dplyr::group_by(., Sequence, `Leading Proteins`, Positions) %>%
    dplyr::summarise(., MSMS_Count = sum(`MS/MS Count`)) %>%
    cSplit(indt = ., splitCols = "Positions", sep = ";", direction = "long")

aa_cover <- evid_coverage %>%
    dplyr::group_by(., MSMS_Count) %>%
    dplyr::summarise(., Freq = n())


pl <- plots_hist(
    data = aa_cover,
    key = "MSMS_Count",
    value = "Freq",
    group = "MSMS_Count",
    fill = "Freq",
    main = "MS/MS count per amino acid",
    transf = "log10",
    label = "Freq",
    legend = "bottom")

ggplot(
    data = data,
    mapping = aes_string(
        x = key,
        y = value)) +
    geom_bar(
        stat = "identity",
        position = "dodge",
        colour = "black") +
    ggtitle(label = main) +
    xlab(label = key) +
    ylab(label = value) +
    theme_bw() +
    theme(
        legend.position = leg,
        title = element_text(
            face = "bold",
            size = (textsize * 1.4)),
        text = element_text(size = textsize),
        plot.title = element_text(
            face = "bold",
            size = (textsize * 1.7)),
        axis.text.x = element_text(
            angle = 0, vjust = 0.5, hjust = 0.5)) +
    #scale_fill_grey(name = fill, start = 0, end = 0.95, na.value = "red") +
    scale_fill_discrete(name = fill) +
    scale_x_discrete(labels = cat.names)
pl




