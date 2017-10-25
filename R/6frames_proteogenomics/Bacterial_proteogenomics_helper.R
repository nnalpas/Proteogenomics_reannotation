


### Define bacterial proteogenomics helper functions ---------------------

# Several functions are used across several script for the bacterial
# proteogenomics pipeline



### Genomic coverage -----------------------------------------------------

# Function to draw rectangular venn diagram of the genomic coverage
plots_rectvenn <- function(
    ideo = NULL,
    ref = NULL,
    pep = NULL) {
    
    # Check whether the correct variables have been submitted by users
    if (is.null(ideo)) {
        stop("Error: No ideogram GRange object was provided!")
    }
    if (is.null(ref)) {
        stop("Error: No reference GRange object was provided!")
    }
    if (is.null(pep)) {
        stop("Error: No peptide GRange object was provided!")
    }
    
    # Get all chromosome nucleotide position
    chrom_nuc <- gr_nucleotide_pos(
        grange = ideo)
    
    # Get all protein-coding associated nucleotide position
    coding_nuc <- gr_nucleotide_pos(
        grange = ref)
    
    # Get all expressed protein associated nucleotide position
    exprs_nuc <- gr_nucleotide_pos(
        grange = ref, filter = "Expressed == TRUE")
    
    # Get all peptide associated nucleotide position
    cover_nuc <- gr_nucleotide_pos(
        grange = pep, filter = 'Database == "Known"')
    
    # Plot a square venn diagram of chromosome coverage
    plot.new()
    rect(
        xleft = 0,
        ybottom = 0,
        xright = 1,
        ytop = 1,
        border = "red",
        lwd = 2)
    text(
        x = 0.015,
        y = 1.02,
        labels = paste(
            "Chromosome", round(length(unique(unlist(chrom_nuc))) / 1000000, 1),
            "Mb", sep = " "),
        col = "red", cex = 1.0, adj = 0)
    rect(
        xleft = 0.01,
        ybottom = 0.01,
        xright = sqrt(
            length(unique(unlist(coding_nuc))) /
                length(unique(unlist(chrom_nuc)))) + 0.01,
        ytop = sqrt(
            length(unique(unlist(coding_nuc))) /
                length(unique(unlist(chrom_nuc)))) + 0.01,
        border = "green",
        lwd = 2)
    text(
        x = 0.025,
        y = sqrt(
            length(unique(unlist(coding_nuc))) /
                length(unique(unlist(chrom_nuc)))) + 0.032,
        labels = paste(
            "Protein-coding",
            round(length(unique(unlist(coding_nuc))) / 1000000, 1),
            "Mb", sep = " "),
        col = "green", cex = 1.0, adj = 0)
    rect(
        xleft = 0.02,
        ybottom = 0.02,
        xright = sqrt(
            length(unique(unlist(exprs_nuc))) /
                length(unique(unlist(chrom_nuc)))) + 0.02,
        ytop = sqrt(
            length(unique(unlist(exprs_nuc))) /
                length(unique(unlist(chrom_nuc)))) + 0.02,
        border = "blue",
        lwd = 2)
    text(
        x = 0.035,
        y = sqrt(
            length(unique(unlist(exprs_nuc))) /
                length(unique(unlist(chrom_nuc)))) + 0.04,
        labels = paste(
            "Expressed protein",
            round(length(unique(unlist(exprs_nuc))) / 1000000, 1),
            "Mb", sep = " "),
        col = "blue", cex = 1.0, adj = 0)
    rect(
        xleft = 0.03,
        ybottom = 0.03,
        xright = sqrt(
            length(unique(unlist(cover_nuc))) /
                length(unique(unlist(chrom_nuc)))) + 0.03,
        ytop = sqrt(
            length(unique(unlist(cover_nuc))) /
                length(unique(unlist(chrom_nuc)))) + 0.03,
        border = "orange",
        lwd = 2)
    text(
        x = 0.045,
        y = sqrt(
            length(unique(unlist(cover_nuc))) /
                length(unique(unlist(chrom_nuc)))) + 0.05,
        labels = paste(
            "Detected peptide",
            round(length(unique(unlist(cover_nuc))) / 1000000, 1),
            "Mb", sep = " "),
        col = "orange", cex = 1.0, adj = 0)
    
    # Recored the plot
    pl <- recordPlot()
    
    # Return the plot result
    return(pl)
    
}


