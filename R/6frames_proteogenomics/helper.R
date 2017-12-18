


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

# Function to determine the position of a peptide within associated proteins
pept_locate <- function(
    data = NULL,
    peptide = "Sequence",
    proteins = "Proteins",
    fasta = NULL) {
    
    # Safety: defined parameters
    if (is.null(data)) {
        stop(paste(
            "No defined data, data should contain at least",
            " peptides and proteins columns!",
            sep = ""))
    }
    if (is.null(fasta)) {
        stop("No defined fasta file, recommended in 'SeqFastaAA' format!")
    }
    
    # Make sure the data is in character format
    data[[peptide]] <- as.character(data[[peptide]])
    data[[proteins]] <- as.character(data[[proteins]])
    
    # Check protein presence in database
    tmp <- nrow(data[!(data[[proteins]] %in% names(fasta)), ])
    if (tmp > 0) {
        warning(paste0(
            "There are ", tmp, " proteins from the peptide dataset that",
            " cannot be found in the database! These will be filtered out!"))
    }
    data <- data[data[[proteins]] %in% names(fasta), ]
    
    # Loop through all rows of peptide
    val <- apply(X = data, MARGIN = 1, FUN = function(x) {
        
        # Assign current entries to temporary variable
        pep <- x[[peptide]]
        prot <- x[[proteins]]
        prot.seq <- fasta[prot]
        
        # Locate all peptide pattern within protein sequence
        tmp <- str_locate_all(
            string = prot.seq %>%
                as.character(.),
            pattern = coll(pattern = pep, ignore_case = FALSE)) %>%
            set_names(prot)
        
        # Check whether locations were obtained
        if (
            length(tmp[[1]][, "start"]) == 0 |
            length(tmp[[1]][, "end"]) == 0) {
            
            # Warn user that there is no match for current peptide/protein
            warning(paste(
                "Current entry:", pep, "does not match in:", prot, sep = " "))
            
            # Add explicitely NA if there is no match
            tmp <- base::data.frame(
                pep, prot,
                start = NA_integer_, end = NA_integer_,
                stringsAsFactors = FALSE)
            
        } else {
            
            # Format as dataframe the current peptide positions
            tmp <- base::data.frame(
                pep, prot, tmp[[1]], stringsAsFactors = FALSE)
            
        }
        
        # Return the location data for current peptide
        return(tmp)
        
    }) %>%
        data.table::rbindlist(.)
    
    # Convert to numeric the peptidde positions
    val$start <- as.numeric(val$start)
    val$end <- as.numeric(val$end)
    
    # Return the overall location data
    return(val)
    
}

# Function to classify novel peptide into potential reasons for their novelty
novel_pep_classify <- function(
    x = NULL,
    coordinate = NULL,
    levenshtein = NULL,
    blast_ref = NULL,
    blast_all = NULL) {
    
    # Check whether the correct variables have been submitted by users
    if (is.null(x)) {
        stop("Error: No novel peptide was provided!")
    }
    if (is.null(coordinate)) {
        stop("Error: No peptide coordinates were provided!")
    }
    if (is.null(levenshtein)) {
        stop("Error: No levenshtein data were provided!")
    }
    if (is.null(blast_ref)) {
        stop("Error: No blast against reference data were provided!")
    }
    if (is.null(blast_all)) {
        stop("Error: No blast against all protein data were provided!")
    }
    
    # Get the peptide position within protein
    coordinate_tmp <- coordinate %>%
        dplyr::filter(., pep == x)
    
    # Check if a unique peptide location is available
    if (nrow(coordinate_tmp) != 1) {
        
        reason <- "Undetermined peptide location"
        
    } else {
        
        # Get the blast and levenshtein data
        blast_ref_tmp <- blast_ref[
            blast_ref$qseqid == coordinate_tmp$prot, ]
        blast_all_tmp <- blast_all[
            blast_all$qseqid == coordinate_tmp$prot, ]
        levenshtein_tmp <- levenshtein[
            levenshtein$Sequence == x &
                levenshtein$id %in% blast_ref_tmp$sseqid, ]
        
        # Check for blast data against reference protein
        if (nrow(blast_ref_tmp) == 0) {
            
            # Check for blast data against all uniprot or ncbi
            if (nrow(blast_all_tmp) == 0) {
                
                # Without blast data mark as potentially novel peptide
                reason <- "Potentially novel"
                
            } else {
                
                # With blast data against all organism mark as known
                # in other species
                reason <- "Known other species"
                
            }
            
        # Check for multiple reciprocal best hits
        } else if (nrow(blast_ref_tmp) > 1) {
            
            # With multiple best hits mark as so
            reason <- "Multiple blast hits"
            
        # Check for single reciprocal best hit
        } else {
            
            # Check whether the peptide is within the matching blast positions
            if (
                coordinate_tmp$start >= blast_ref_tmp$qstart_blast &
                coordinate_tmp$end <= blast_ref_tmp$qend_blast) {
                
                # If yes then mark as potential SAV
                reason <- "Potential SAV"
                
                # Check if the levenshtein data exist and
                # what is the distance score
                if (
                    nrow(levenshtein_tmp) > 0 &
                    unique(levenshtein_tmp$leven) == 1) {
                    
                    # If yes mark as confirmed SAV
                    reason <- "SAV"
                    
                } else if (
                    nrow(levenshtein_tmp) > 0 &
                    unique(levenshtein_tmp$leven) == 0) {
                    
                    # If yes not novel peptide
                    reason <- "Not novel"
                    
                } else if (
                    nrow(levenshtein_tmp) > 0 &
                    unique(levenshtein_tmp$leven) > 1) {
                    
                    # If yes multi SAVs or InDels
                    reason <- "SAVs/InDels"
                    
                }
                
            # Check whether the peptide overlap the reference protein start
            } else if (min(
                coordinate_tmp$start,
                blast_ref_tmp$qstart_blast) == coordinate_tmp$start) {
                
                # Check whether novel ORF and reference protein have
                # incompatible start site
                if (blast_ref_tmp$sstart_blast > blast_ref_tmp$qstart_blast) {
                    
                    reason <- "Incompatible start site"
                    
                } else {
                    
                    # Without blast data mark as potential alternate start
                    reason <- "Potential alternate start"
                    
                    # Check for blast data against all uniprot or ncbi
                    if (nrow(blast_all_tmp) != 0) {
                        
                        if (any(and(
                            coordinate_tmp$start >= blast_all_tmp$qstart_blast,
                            coordinate_tmp$end <= blast_all_tmp$qend_blast))) {
                            
                            # With blast data mark as alternate start
                            # characterised in other species
                            reason <- "Alternate start (known other species)"
                            
                        }
                        
                    }
                    
                }
                
            # Check whether the peptide overlap the reference protein end
            } else if (max(
                coordinate_tmp$end,
                blast_ref_tmp$qend_blast) == coordinate_tmp$end) {
                
                # Check whether novel ORF and reference protein have
                # incompatible end site
                if (blast_ref_tmp$send_blast > blast_ref_tmp$qend_blast) {
                    
                    reason <- "Incompatible end site"
                    
                } else {
                    
                    # Without blast data mark as potential alternate end
                    reason <- "Potential alternate end"
                    
                    # Check for blast data against all uniprot or ncbi
                    if (nrow(blast_all_tmp) != 0) {
                        
                        if (any(and(
                            coordinate_tmp$start >= blast_all_tmp$qstart_blast,
                            coordinate_tmp$end <= blast_all_tmp$qend_blast))) {
                            
                            # With blast data mark as alternate end
                            # characterised in other species
                            reason <- "Alternate end (known other species)"
                            
                        }
                        
                    }
                    
                }
                
            } else {
                
                reason <- "Undetermined reason"
                
            }
            
        }
        
    }
    
    # Return the reason for novelty
    return(reason)
    
}



### Number utilities functions -------------------------------------------

# Function to keep only decimals
keep_decimals <- function(x) {
    x - floor(x)
}



### MaxQuant related functions -------------------------------------------

# Function to read in any MaxQuant file from specified path
mq_read <- function(
    path = NULL,
    name = NULL,
    ...) {
    
    # Safety: defined parameters
    if (is.null(path)) {
        stop("No defined path to maxquant file!")
    }
    if (is.null(name)) {
        stop("No defined maxquant file name!")
    }
    
    # Locate the table file
    maxq.file <- list.files(
        path = path,
        pattern = name,
        full.names = TRUE)
    
    # Safety: located msScan table
    if (length(maxq.file) != 1) {
        stop("Number of located maxquant file not appropriate!")
    }
    
    # Read in the maxquant file
    maxquant <- fread(
        input = maxq.file, sep = "\t", header = TRUE,
        na.strings = c("NA", "N/A"), stringsAsFactors = FALSE,
        strip.white = TRUE, ...)
    
    # Return the maxquant data
    return(as.data.frame(x = maxquant, stringsAsFactors = FALSE))
    
}

# Define the function that filter Reverse and Contaminants hits from data
mq_rev_con_filt <- function (
    data = NULL) {
    
    # Check that the data variable have been provided
    if (is.null(data)) {
        stop("Error: User did not provide data!")
    }
    
    # Check that the data variable is of class dataframe
    if (!is.data.frame(data)) {
        stop("Error: User must provide data of class dataframe!")
    }
    
    # Define the column names to look for in data variable
    rev.col <- "Reverse"
    cont <- c("Potential.contaminant", "Potential contaminant")
    
    # Check whether filtering can be performed on Reverse column
    if (rev.col %in% colnames(data) & all(!is.na(data[[rev.col]]))) {
        data.noRev <- data[data[[rev.col]] != "+", ]
    } else {
        warning("Warning: No Reverse column in the provided data!")
        data.noRev <- data
    }
    
    # Check whether filtering can be performed on Contaminant column
    if (
        any(cont %in% colnames(data)) &
        all(!is.na(data[[cont[cont %in% colnames(data)]]]))) {
        cont.col <- cont[cont %in% colnames(data)]
        data.noCont <- data.noRev[data.noRev[[cont.col]] != "+", ]
    } else {
        warning("Warning: No Contaminant column in the provided data!")
        data.noCont <- data.noRev
    }
    
    # Output filtered data
    return(data.noCont)
}



### Plotting wrapper function --------------------------------------------

# Create function to draw boxplot (with-without zoom) acccording to
# certain key-values parameters
plots_box <- function(
    data = NULL,
    key = NULL,
    value = NULL,
    fill = "white",
    colour = "black",
    main = "Boxplot",
    xlabel = NULL,
    ylabel = NULL,
    textsize = 15,
    zoom = NULL,
    transf = "none",
    outlier_simplify = FALSE,
    label = NULL,
    bw = FALSE,
    xdir = "horizontal",
    ydir = "horizontal",
    auto_scale = 1) {
    
    # Check whether the data, key and value have been submitted by users
    if (is.null(data)) {
        stop("Error: No data were provided!")
    }
    if (is.null(key)) {
        stop("Error: No key column was provided!")
    }
    if (is.null(value)) {
        stop("Error: No value column was provided!")
    }
    xdir %<>%
        match.arg(arg = ., choices = c("horizontal", "vertical"))
    ydir %<>%
        match.arg(arg = ., choices = c("horizontal", "vertical"))
    
    # Create the list that will hold all generated plots
    myplot <- list()
    
    # Check whether the data values are integer64 or not numeric class
    if (is.integer64(data[[value]])) {
        
        # Give warning that values are integer64
        warning(paste(
            "Data values: '", value,
            "' are in integer 64, will attempt transformation to double!",
            sep = ""))
        
        # Change the data from integer64 to double for compatibility
        data[[value]] <- as.double(data[[value]])
        
    }
    if (!is.numeric(data[[value]])) {
        
        # Give warning that values are not numeric
        warning(paste(
            "Data values: '", value,
            "' are not numeric, will attempt transformation to numeric!",
            sep = ""))
        
        # Change the data to numeric
        data[[value]] <- as.numeric(as.character(data[[value]]))
        
    }
    
    # Reformat the key as factor
    data[[key]] <- factor(
        x = data[[key]],
        levels = as.character(unique(data[[key]])),
        labels = as.character(unique(data[[key]])),
        ordered = TRUE)
    
    # Format category name
    #cat.names <- data[[key]] %>%
    #    levels(.) %>%
    #    stringr::str_wrap(string = ., width = 15)
    
    # Check which direction to display x axis labels
    if (xdir == "horizontal") {
        x_custom <- element_text(
            angle = 0, vjust = 0.5, hjust = 0.5)
    } else if (xdir == "vertical") {
        x_custom <- element_text(
            angle = 90, vjust = 0.5, hjust = 1)
    }
    
    # Check which direction to display y axis labels
    if (ydir == "horizontal") {
        y_custom <- element_text(
            angle = 0, vjust = 0.5, hjust = 1)
    } else if (ydir == "vertical") {
        y_custom <- element_text(
            angle = 90, vjust = 0.5, hjust = 0.5)
    }
    
    # Check if fill value is a column or not and should be aes or not
    if (fill %in% colnames(data)) {
        
        # Reformat the fill as factor
        data[[fill]] <- factor(
            x = data[[fill]],
            levels = as.character(unique(data[[fill]])),
            labels = as.character(unique(data[[fill]])),
            ordered = TRUE)
        
        # Create the plot with fill as aes
        plot <- ggplot(
            data = data,
            mapping = aes_string(
                x = key, y = value, fill = fill)) +
            stat_boxplot(geom = "errorbar") +
            geom_boxplot(colour = colour)
        
    } else {
        
        # Create the plot with fill as not aes
        plot <- ggplot(
            data = data,
            mapping = aes_string(
                x = key, y = value)) +
            stat_boxplot(geom = "errorbar", width = 0.1) +
            geom_boxplot(fill = fill, colour = colour)
        
    }
    
    # Check whether axis labels were provided
    if (is.null(xlabel)) {
        xlabel <- key
    }
    if (is.null(ylabel)) {
        ylabel <- value
    }
    
    # Continue plot creation
    plot <- plot +
        xlab(label = xlabel) +
        ylab(label = ylabel) +
        ggtitle(label = main) +
        theme_bw() +
        theme(
            legend.position = "bottom",
            title = element_text(
                face = "bold",
                size = (textsize * 1.25)),
            text = element_text(size = textsize),
            plot.title = element_text(
                face = "bold",
                size = (textsize * 1.5)),
            axis.text.x = x_custom,
            axis.text.y = y_custom)# +
    #scale_x_discrete(labels = cat.names)
    
    # Check whether the plot should use black and white colors
    if (bw) {
        plot <- plot +
            scale_fill_grey(name = fill)
    } else {
        plot <- plot +
            scale_fill_discrete(name = fill)
    }
    
    # Check whether individual label should also be displayed
    if (!is.null(label)) {
        
        # Define whether fill should be used for grouping data
        group_var <- key
        if (fill %in% colnames(data)) {
            group_var <- c(group_var, fill)
        }
        
        # Define variable for grouping
        group_var %<>%
            lapply(., as.symbol)
        
        # Define variable to calculate on
        calc_var <- as.symbol(value)
        
        # Compute IQR, size, mean for the current data
        iqr <- data %>%
            as.data.frame(.) %>%
            dplyr::group_by_(., .dots = group_var) %>%
            dplyr::summarise_(
                .,
                N = ~n(),
                Median = lazyeval::interp(~ median(var), var = calc_var),
                Mean = lazyeval::interp(~ mean(var), var = calc_var),
                Min = lazyeval::interp(~ min(var), var = calc_var),
                Q1 = lazyeval::interp(~ quantile(var, probs = 0.25), var = calc_var),
                Q3 = lazyeval::interp(~ quantile(var, probs = 0.75), var = calc_var),
                Max = lazyeval::interp(~ max(var), var = calc_var)) %>%
            dplyr::mutate(
                .,
                IQR =  (Q3 - Q1),
                Q1_outlier = (Q1 - (1.5 * IQR)),
                Q3_outlier = (Q3 + (1.5 * IQR)),
                N_label = paste0("n = ", N),
                Md_label = paste0(
                    "md = ",
                    format(x = round(x = Median, digits = 1), nsmall = 1)),
                M_label = paste0(
                    "M = ",
                    format(x = round(x = Mean, digits = 1), nsmall = 1)))
        
        # Include all the labels requested by user
        v_adjust <- (0.05 * length(label)) + 0.05
        for (j in 1:length(label)) {
            
            iqr %<>%
                ungroup() %>%
                dplyr::mutate(
                    .,
                    y_val = (Max + (max(Max) * (v_adjust - (0.05 * j)))))
            y_val <- "y_val"
            plot <- plot +
                stat_summary(
                    data = iqr,
                    mapping = aes_string(
                        x = key,
                        y = y_val,
                        group = fill,
                        label = label[j]),
                    geom = "text",
                    position = position_dodge(width = 0.75))
            iqr$y_val <- NULL
            
        }
        
    }
    
    # Check whether key are of factor type
    if (is.factor(data[[key]])) {
        
        # Get the break value to display on x-axis
        breaks <- custom_scale(x = data[[key]] %>% levels(.), n = auto_scale)
        
        # Add a discrete x-axis with breaks to plot
        plot <- plot +
            scale_x_discrete(
                breaks = breaks,
                labels = breaks)
        
    }
    
    # Check whether the y axis need to be updated
    if (transf != "none") {
        
        plot <- conti_axis(
            plot = plot,
            axis = "y",
            transf = transf)
        
    }
    
    # Store the gTree plot into a list
    myplot["Boxplot"] <- list(plot)
    
    # Check whether the y axis need to be updated
    if (!is.null(zoom)) {
        
        plot <- conti_axis(
            plot = plot,
            axis = "y",
            zoom = zoom,
            transf = transf)
        
        # Store the gTree plot into a list
        myplot["Boxplot_zoom"] <- list(plot)
        
    }
    
    # Check whether outlier simplification was requested
    if (outlier_simplify) {
        
        myplot <- purrr::map(.x = myplot, .f = outlier_sampling)
        
    }
    
    # Return the plot list
    return(myplot)
    
}

# Create function to draw histogram (with-without zoom) acccording to
# certain key-values parameters
plots_hist <- function(
    data = NULL,
    key = NULL,
    value = NULL,
    group = NULL,
    fill = "grey",
    colour = "black",
    alpha = NULL,
    main = "Histogram",
    xlabel = NULL,
    ylabel = NULL,
    posit = "dodge",
    textsize = 15,
    zoom = NULL,
    transf = "none",
    label = NULL,
    legend = c("none", "bottom", "top", "left", "right"),
    bw = FALSE,
    xdir = "horizontal",
    ydir = "horizontal",
    auto_scale = 1) {
    
    # Check whether the data, key and value have been submitted by users
    if (is.null(data)) {
        stop("Error: No data were provided!")
    }
    if (is.null(key)) {
        stop("Error: No key column was provided!")
    }
    if (is.null(group)) {
        stop("Error: No group column was provided!")
    }
    if (is.null(fill)) {
        stop("Error: No fill column was provided!")
    }
    if (is.null(value)) {
        stop("Error: No value column was provided!")
    }
    xdir %<>%
        match.arg(arg = ., choices = c("horizontal", "vertical"))
    ydir %<>%
        match.arg(arg = ., choices = c("horizontal", "vertical"))
    
    # Create the list that will hold all generated plots
    myplot <- list()
    
    # Check whether the data values are integer64 or not numeric class
    if (is.integer64(data[[value]])) {
        
        # Give warning that values are integer64
        warning(paste(
            "Data values: '", value,
            "' are in integer 64, will attempt transformation to double!",
            sep = ""))
        
        # Change the data from integer64 to double for compatibility
        data[[value]] <- as.double(data[[value]])
        
    }
    if (!is.numeric(data[[value]])) {
        
        # Give warning that values are not numeric
        warning(paste(
            "Data values: '", value,
            "' are not numeric, will attempt transformation to numeric!",
            sep = ""))
        
        # Change the data to numeric
        data[[value]] <- as.numeric(as.character(data[[value]]))
        
    }
    
    # Check whether a legend is required
    leg <- match.arg(legend)
    
    # Reformat the key as factor
    if (!is.factor(data[[key]])) {
        
        data[[key]] <- factor(
            x = data[[key]],
            levels = as.character(unique(data[[key]])),
            labels = as.character(unique(data[[key]])),
            ordered = TRUE)
        
    }
    
    # Reformat the group as factor
    data[[group]] <- factor(
        x = data[[group]],
        levels = as.character(unique(data[[group]])),
        labels = as.character(unique(data[[group]])),
        ordered = TRUE)
    
    # Format category name
    cat.names <- data[[key]] %>%
        levels(.) %>%
        stringr::str_wrap(string = ., width = 15) %>%
        gsub("_", " ", .)
    
    # Check whether axis labels were provided
    if (is.null(xlabel)) {
        xlabel <- key
    }
    if (is.null(ylabel)) {
        ylabel <- value
    }
    
    # Check which direction to display x axis labels
    if (xdir == "horizontal") {
        x_custom <- element_text(
            angle = 0, vjust = 0.5, hjust = 0.5)
    } else if (xdir == "vertical") {
        x_custom <- element_text(
            angle = 90, vjust = 0.5, hjust = 1)
    }
    
    # Check which direction to display y axis labels
    if (ydir == "horizontal") {
        y_custom <- element_text(
            angle = 0, vjust = 0.5, hjust = 1)
    } else if (ydir == "vertical") {
        y_custom <- element_text(
            angle = 90, vjust = 0.5, hjust = 0.5)
    }
    
    # Check if fill value is a column or not and should be aes or not
    if (fill %in% colnames(data)) {
        
        # Reformat the fill as factor
        data[[fill]] <- factor(
            x = data[[fill]],
            levels = as.character(unique(data[[fill]])),
            labels = as.character(unique(data[[fill]])),
            ordered = TRUE)
        
        # Create the plot with fill as aes
        plot <- ggplot(
            data = data,
            mapping = aes_string(
                x = key,
                y = value,
                group = group,
                fill = fill)) +
            geom_bar(
                stat = "identity",
                position = posit,
                colour = colour)
        
    } else {
        
        # Create the plot with fill as not aes
        plot <- ggplot(
            data = data,
            mapping = aes_string(
                x = key,
                y = value,
                group = group)) +
            geom_bar(
                stat = "identity",
                position = posit,
                fill = fill,
                colour = colour)
        
    }
    
    # Continue plot creation
    plot <- plot +
        ggtitle(label = main) +
        xlab(label = xlabel) +
        ylab(label = ylabel) +
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
            axis.text.x = x_custom,
            axis.text.y = y_custom)# +
    #scale_x_discrete(labels = cat.names)
    
    # Check whether the plot should use black and white colors
    if (bw) {
        plot <- plot +
            scale_fill_grey(name = fill)
    } else {
        plot <- plot +
            scale_fill_discrete(name = fill)
    }
    
    # Check whether individual label should also be displayed
    if (!is.null(label)) {
        
        plot <- gg_labels(
            plot = plot,
            data = data,
            label = label,
            posit = posit,
            textsize = textsize)
        
    }
    
    # Check whether key are of factor type
    if (is.factor(data[[key]])) {
        
        # Get the break value to display on x-axis
        breaks <- custom_scale(x = data[[key]] %>% levels(.), n = auto_scale)
        
        # Add a discrete x-axis with breaks to plot
        plot <- plot +
            scale_x_discrete(
                breaks = breaks,
                labels = breaks)
        
    }
    
    # Check whether the y axis need to be updated
    if (transf != "none") {
        
        plot <- conti_axis(
            plot = plot,
            axis = "y",
            transf = transf)
        
    }
    
    # Store the gTree plot into a list
    myplot["histplot"] <- list(plot)
    
    # Check whether the y axis need to be updated
    if (!is.null(zoom)) {
        
        plot <- conti_axis(
            plot = plot,
            axis = "y",
            zoom = zoom,
            transf = transf)
        
        # Store the gTree plot into a list
        myplot["histplot_zoom"] <- list(plot)
        
    }
    
    # Return the plot list
    return(myplot)
    
}



### Utility function for plot reformatting -------------------------------

# Function to add labels on a ggplot in an optimised fashion, such as
# with vertical or horizontal display
gg_labels <- function(
    plot = NULL,
    data = NULL,
    label = NULL,
    posit = NULL,
    textsize = 15) {
    
    # Check whether the correct variables have been submitted by users
    if (is.null(plot)) {
        stop("Error: No plot was provided!")
    }
    if (is.null(data)) {
        stop("Error: No data were provided!")
    }
    if (is.null(label)) {
        stop("Error: No label values were provided!")
    }
    if (is.null(posit)) {
        stop("Error: No plot position value was provided!")
    }
    
    # Set the label display angle to 0
    angle <- 0
    vjust <- -0.3
    hjust <- 0.5
    
    # If more than 20 labels, use vertical label display
    if (length(data[[label]]) > 20) {
        angle <- 90
        vjust <- 0.5
        hjust <- -0.3
    }
    
    # Check the position parameter
    if (posit == "dodge") {
        posit <- position_dodge(width = 0.9)
    } else if (posit == "stack") {
        posit <- position_stack()
    } else {
        stop("This type of position was never tested with geom_text!")
    }
    
    # Add labels to the plot
    plot <- plot +
        geom_text(
            mapping = aes_string(
                label = deparse(
                    formatC(x = data[[label]], format = 'd', big.mark = ','))),
            position = posit,
            vjust = vjust,
            hjust = hjust,
            size = (textsize * 0.3),
            check_overlap = TRUE,
            angle = angle)
    
    # Return the updated plot
    return(plot)
    
}

# Function to simplify the x-axis by selecting only a few labels to display
custom_scale <- function(x = NULL, n = 2) {
    
    # Check whether the data has been submitted by users
    if (is.null(x)) {
        stop("Error: No vector of factor levels provided!")
    }
    
    # Keep only every n labels
    filt <- seq(from = 1, to = length(x), by = n)
    res <- x[filt]
    
    # Return the breaks results
    return(res)
    
}


