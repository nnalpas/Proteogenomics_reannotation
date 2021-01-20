
library(gplots)
library(ggplot2)


#my_filt <- read.table(
#    file = paste0(dirname(opt$input), "/Novel_targets.txt"),
#    header = FALSE, sep = "\t", quote = "", as.is = TRUE) %>%
#    unlist(., use.names = FALSE)

my_filt <- read.table(
    file = paste0(dirname(opt$input), "/Novel_targets2.txt"),
    header = FALSE, sep = "\t", quote = "", as.is = TRUE) %>%
    unlist(., use.names = FALSE)

#
blast_data_filter <- blast_data %>%
    dplyr::filter(., qseqid %in% my_filt)

#
blast_data_format <- blast_data_filter %>%
    dplyr::group_by(., qseqid, staxid) %>%
    dplyr::summarise(
        ., pident = max(pident), nident = max(nident), evalue = min(evalue),
        Sig_label = ifelse(evalue <= 0.0001, "*", ""))


blast_data_format$qseqid <- factor(
    x = blast_data_format$qseqid,
    levels = unique(blast_data_format$qseqid),
    ordered = TRUE)



blast_data_format$staxid <- factor(
    x = blast_data_format$staxid,
    levels = unique(blast_data_format$staxid),
    ordered = TRUE)

# 
heatmap_pl <- ggplot(
    data = blast_data_format,
    #aes(x = qseqid, y = sseqid)) +
    aes(x = qseqid, y = staxid)) +
    geom_tile(aes(fill = pident)) +
    geom_text(aes(label = Sig_label), size = 0.5) +
    scale_fill_gradient2(low = "#232376", high = "#F43D00") +
    theme_bw() +
    theme(
        title = NULL,
        legend.position = "left",
        #axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        #axis.ticks.x = element_blank(),
        #axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        #axis.ticks.y = element_blank(),
        axis.text.x = element_text(
            angle = 90, vjust = 0.5, hjust = 1)) +
    xlab("Novel ORFs") +
    ylab("Blast results")



### Data taxon association -----------------------------------------------

library(RSQLite)
options(sqldf.driver = "SQLite")

# Define path to taxonomy SQL database
tax_datab <- "L:/Data/SQL_Taxonomy/Taxon.sqlite"


# Open/create database
tax <- dbConnect(SQLite(), dbname = tax_datab)
taxon_level <- "species"

# Get the level ID for specific level name
lev <- as.numeric(dbGetQuery(
    conn = tax,
    statement = paste(
        "SELECT LevelID FROM TaxonLevel WHERE Levelname = '",
        taxon_level, "'", sep = "")))

# Define variables to hold taxonomy annotation data
final_tax <- data.frame(matrix(ncol = 6))
colnames(final_tax) <- c(
    "TaxonID", "Taxonname", "Uniquename", "LevelID", "ParentID",
    "StartingTaxon")
final_tax <- final_tax[!is.na(final_tax$TaxonID),]

# 
val <- data.frame(TaxonID = unique(blast_data$staxid))
val$TargetTaxonID <- NA

# Loop through all taxon ID associated to all protein ID
for (targ in unique(val$TaxonID)) {
    
    # Get the complete ancestor tree for a taxon ID
    hierarc <- dbGetQuery(
        conn = tax,
        statement = paste(
            "WITH RECURSIVE
            find_parent(n) AS (
            VALUES(", targ, ")
            UNION
            SELECT ParentID FROM Taxonomy, find_parent
            WHERE Taxonomy.TaxonID = find_parent.n LIMIT 50)
            SELECT * FROM Taxonomy WHERE Taxonomy.TaxonID IN find_parent",
            sep = ""))
    
    if (nrow(hierarc) > 0) {
        
        # Check whether the target level is present in taxonomy tree
        if (any(hierarc$LevelID == lev)) {
            
            # Get the Taxon ID for target level ID
            final_targ <- hierarc[hierarc$LevelID == lev, "TaxonID"]
            
        } else {
            
            # Reorganise the taxonomy tree from branch to root
            tag <- targ
            reform_hierarc <- hierarc[hierarc$TaxonID == tag, ]
            while (tag != 1) {
                tag <- reform_hierarc[nrow(reform_hierarc), "ParentID"]
                tmp_hierar <- hierarc[hierarc$TaxonID == tag, ]
                reform_hierarc <- rbind(reform_hierarc, tmp_hierar)
            }
            
            # Define the closest level ID to the missing target level
            tmp_lev <- reform_hierarc[reform_hierarc$LevelID < lev, "LevelID"] %>%
                max(.)
            
            # Get the Taxon ID for newly defined level ID
            final_targ <- reform_hierarc[
                which(reform_hierarc$LevelID == tmp_lev), "TaxonID"]
            
        }
        
        # Associate target taxon ID to protein ID
        val[val$TaxonID == targ, "TargetTaxonID"] <- final_targ
        
        # Keep track of association between initial taxon ID and target taxon ID
        final_tax <- rbind(
            final_tax,
            data.frame(
                hierarc[hierarc$TaxonID == final_targ, ],
                StartingTaxon = targ))
        
    }
    
}
final_tax <- unique(final_tax)


#
blast_data_format$staxid <- as.integer(blast_data_format$staxid)
blast_data_species <- final_tax %>%
    dplyr::select(., Taxonname, StartingTaxon) %>%
    dplyr::left_join(
        x = blast_data_format, y = ., by = c("staxid" = "StartingTaxon")) %>%
    dplyr::mutate(
        ., staxidNEW = ifelse(
            is.na(Taxonname), as.character(staxid), Taxonname))


blast_data_species$qseqid <- factor(
    x = blast_data_species$qseqid,
    levels = unique(blast_data_species$qseqid),
    ordered = TRUE)


blast_data_species$staxidNEW <- factor(
    x = blast_data_species$staxidNEW,
    levels = unique(blast_data_species$staxidNEW),
    ordered = TRUE)

# 
heatmap_pl <- ggplot(
    data = blast_data_species,
    #aes(x = qseqid, y = sseqid)) +
    #aes(x = qseqid, y = staxid)) +
    aes(x = qseqid, y = staxidNEW)) +
    geom_tile(aes(fill = pident)) +
    geom_text(aes(label = Sig_label), size = 0.5) +
    scale_fill_gradient2(low = "#232376", high = "#F43D00") +
    theme_bw() +
    theme(
        title = NULL,
        legend.position = "left",
        #axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        #axis.ticks.x = element_blank(),
        #axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        #axis.ticks.y = element_blank(),
        axis.text.x = element_text(
            angle = 90, vjust = 0.5, hjust = 1)) +
    xlab("Novel ORFs") +
    ylab("Blast results")


