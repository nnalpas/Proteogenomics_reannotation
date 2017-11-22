


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



### Packages import ------------------------------------------------------

# Create a function to load or install (then load) the required packages
load_package <- function(package) {
    
    if (
        suppressWarnings(require(
            package = eval(package),
            character.only = TRUE,
            quietly = TRUE))) {
        print(paste0(
            package, " is loaded correctly!"))
    }
    else {
        print(paste0(
            "Trying to install ",
            package))
        source(file = "http://bioconductor.org/biocLite.R", verbose = FALSE)
        biocLite(pkgs = eval(package), suppressUpdates = TRUE)
        if(
            require(
                package = eval(package),
                character.only = TRUE,
                quietly = TRUE)) {
            print(paste0(
                package,
                " is correctly installed and loaded!"))
        }
        else {
            stop(paste0(
                '"', "Could not install ",
                package, '"'))
        }
    }
    
    print(paste0(
        package, " version: ",
        try(utils::packageVersion(pkg = package))))
    
}



### Blast related functions ----------------------------------------------

# Function to read in any Blast file from specified path
blast_read <- function(
    file = NULL,
    blast_format = c(
        "0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11"),
    columns = NULL,
    ...) {
    
    # Safety: defined parameters
    if (is.null(file)) {
        stop("No defined Blast file name!")
    }
    format <- match.arg(blast_format)
    if (format == "6") {
        if (is.null(columns)) {
            columns <- c(
                "qseqid", "sseqid", "pident", "nident", "mismatch",
                "gaps", "length", "gapopen", "qstart", "qend", "qlen",
                "qframe", "qseq", "sstart", "send", "slen", "sframe",
                "sseq", "staxid", "ssciname", "sstrand", "evalue",
                "bitscore", "score")
            warning("Using default column names!")
        }
    } else {
        stop("The Blast format is not yet implemented, check with developper!")
    }
    
    # Read in the blast file
    blast_data <- fread(
        input = file, sep = "\t", header = FALSE,
        stringsAsFactors = FALSE, integer64 = "double",
        col.names = columns, quote = "", data.table = FALSE,
        ...)
    
    # Return the maxquant data
    return(blast_data)
    
}

# Function that select the best blast hit based on minimum e-value,
# maximum score and maximum identity sequentially
best_blast <- function(
    data = NULL,
    key = NULL,
    filter = NULL,
    multi_match = c("remove", "keep", "uniquify")) {
    
    # Check whether the data, key and filter have been submitted by users
    if (is.null(data)) {
        stop("Error: No data were provided!")
    }
    if (is.null(key)) {
        stop("Error: No key column was provided!")
    }
    multiple <- match.arg(multi_match)
    if (is.null(multi_match)) {
        multiple <- "keep"
    }
    
    # Turn the key column into symbol
    col.symb <- as.symbol(key)
    
    # Filter the data with specific fields in sequence
    if (is.null(filter)) {
        
        # Filter based on default
        data.filt <- data %>%
            unique(.) %>%
            dplyr::group_by_(., .dots = col.symb) %>%
            dplyr::filter(., evalue == min(evalue)) %>%
            dplyr::filter(., score == max(score)) %>%
            dplyr::filter(., pident == max(pident)) %>%
            dplyr::mutate(., best_count = n()) %>%
            base::as.data.frame(., stringsAsFactors = FALSE)
        
    } else {
        
        # Evaluate the filter parameter, this is required when parameters
        # contains space and was submitted from command line
        filter <- eval(filter)
        
        # Filter based on user provided criteria
        data.filt <- data %>%
            unique(.) %>%
            dplyr::group_by_(., .dots = col.symb)
        for (x in filter) {
            data.filt %<>%
                dplyr::filter_(., filter)
        }
        data.filt %<>%
            dplyr::mutate(., best_count = n()) %>%
            base::as.data.frame(., stringsAsFactors = FALSE)
        
    }
    
    # Check what should be done with entries matching multiple times
    if (length(unique(data[[key]])) == length(unique(data.filt[[key]]))) {
        warning("All blast hits matched uniquely!")
    } else if (multiple == "remove") {
        data.filt %<>%
            dplyr::filter(., best_count == 1)
        warning("All multiple blast hits were removed!")
    } else if (multiple == "uniquify") {
        orig <- data.filt$qseqid %>%
            unique(.) %>%
            length(.)
        data.filt %<>%
            dplyr::group_by(., qseqid) %>%
            dplyr::arrange(., sseqid) %>%
            dplyr::mutate(., rank_s = c(1:n())) %>%
            dplyr::ungroup() %>%
            dplyr::group_by(., sseqid) %>%
            dplyr::arrange(., qseqid) %>%
            dplyr::mutate(., rank_q = c(1:n())) %>%
            dplyr::ungroup() %>%
            dplyr::filter(., rank_q == rank_s) %>%
            dplyr::select(., -rank_q, -rank_s)
        retained <- data.filt$qseqid %>%
            unique(.) %>%
            length(.)
        warning(paste(
            "About", (orig - retained), "entries could not be uniquified!"))
    } else if (multiple == "keep") {
        warning("There may still be entries with multiple blast hits!")
    }
    
    # Return the best blast data
    return(base::as.data.frame(data.filt))
    
}



### Utility function for data reformatting -------------------------------

# Function to clean up the UniProt IDs (in case regular expression
# was missing in andromeda)
uni_id_clean <- function(x) {
    
    # Check submitted data is in the form of a character vector
    if (is.vector(x) & is.character(x)) {
        
        # Split ids and then use regular expression to get exact ids
        x %<>%
            strsplit(x = ., split = ";") %>%
            lapply(X = ., FUN = function(x) {
                val <- sub(
                    pattern = "^(REV__)?(gi.+|sp|tr)\\|(.+)\\|.*$",
                    replacement = "\\1\\3",
                    x = x) %>%
                    paste(., collapse = ";")
            }) %>%
            unlist(.)
        return(x)
        
    } else {
        
        stop("Data is not vector or not character!")
        
    }
    
}



### Computation functions ------------------------------------------------

# Function to determine the translation frame based on start coordinates
get_frame <- function(
    data = NULL,
    start = NULL,
    strand = NULL,
    genome_size = NULL) {
    
    # Check whether the correct variables have been submitted by users
    if (is.null(data)) {
        stop("Error: No dataframe of ORF coordinates provided!")
    }
    if (is.null(start)) {
        stop("Error: No start coordinate column provided!")
    }
    if (is.null(strand)) {
        stop("Error: No strand column provided!")
    }
    if (is.null(genome_size)) {
        stop("Error: No genome size provided!")
    }
    
    # Determine the frame of translation in a strand specific manner
    tmp <- as.integer(unlist(data[data[[strand]] == "+", start]))
    data[data[[strand]] == "+", "frame"] <- ((tmp + 2) / 3) %>%
        keep_decimals(.) %>%
        round(x = ((. + 0.33) * 3))
    tmp <- as.integer(unlist(data[data[[strand]] == "-", start]))
    data[data[[strand]] == "-", "frame"] <- ((genome_size - tmp + 3) / 3) %>%
        keep_decimals(.) %>%
        round(x = ((. + 0.33) * -3))
    
    # Return the dataframe with frame information
    return(data)
    
}


