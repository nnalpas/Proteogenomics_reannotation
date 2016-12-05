


setwd("F:/data/Vaishnavi/combined - 6 frame translation/txt")

user <- Sys.info()[["user"]]

source(
    file = paste(
        "C:/Users/",
        user,
        "/Documents/GitHub/Miscellaneous/R/General/General_function.R",
        sep = ""))

# Load the required packages (or install if not already in library)
# Note: the select function from dplyr and UniProt.ws is overlapping,
# therefore order of package loading is important
loadpackage(plyr)
loadpackage(dplyr)
loadpackage(tidyr)
loadpackage(seqinr)
loadpackage(UniProt.ws)
loadpackage(magrittr)
loadpackage(WriteXLS)
loadpackage(data.table)
loadpackage(splitstackshape)
loadpackage(VennDiagram)
loadpackage(ggplot2)
loadpackage(grid)
loadpackage(gridExtra)
loadpackage(RColorBrewer)
loadpackage(stringr)
loadpackage(Biostrings)
loadpackage(RecordLinkage)
loadpackage(VariantAnnotation)
loadpackage(cgdsr)
loadpackage(bit64)

evid <- maxquant.read(path = ".", name = "evidence.txt")


evid.filt <- evid
evid.filt$group <- "Known"
evid.filt[grep(pattern = "^seq_\\d+(;seq_\\d+)*$", x = evid.filt$Proteins), "group"] <- "Novel"


# Calculate quantile values of number of protein per groups field
data.filt <- evid.filt %>%
    dplyr::select(., group, PEP) %>%
    base::as.data.frame(., stringsAsFactors = FALSE)

# Create the boxplot using ggplot
textsize <- 15
ggplot(data = data.filt, mapping = aes(x = factor(group), y = PEP)) +
    geom_boxplot() +
    xlab(label = "Peptide group") +
    ylab(label = "PEP") +
    ggtitle(label = "PEP per peptide category") +
    theme(
        legend.position = "bottom",
        title = element_text(
            face = "bold",
            size = (textsize * 1.25)),
        text = element_text(size = textsize),
        plot.title = element_text(
            face = "bold",
            size = (textsize * 1.5)),
        axis.text.x = element_text(
            angle = 0, vjust = 1, hjust = 0))




plot <- ggplot(
    data = data.filt,
    mapping = aes(
        x = msrun, ymin = `0%`, lower = `25%`, middle = `50%`,
        upper = `75%`, ymax = `100%`,
        fill = species, alpha = level)) +
    #scale_alpha_manual(values = c(1, 0.5)) +
    geom_boxplot(stat = "identity") +
    xlab(label = "MS run") +
    ylab(label = ytitle) +
    ggtitle(label = type) +
    theme(
        legend.position = "bottom",
        title = element_text(
            face = "bold",
            size = (textsize * 1.25)),
        text = element_text(size = textsize),
        plot.title = element_text(
            face = "bold",
            size = (textsize * 1.5)),
        axis.text.x = element_text(
            angle = -45, vjust = 1, hjust = 0)) +
    guides(fill = guide_legend(ncol = 5)) +
    scale_x_discrete(
        labels = gsub(
            pattern = "^(.*?\\..*?)\\.(.*)$",
            replacement = "\\1.\n\\2",
            x = unique(data.filt$msrun))) +
    coord_cartesian(ylim = ylimit)


