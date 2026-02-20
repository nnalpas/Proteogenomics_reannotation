


### Define bacterial proteogenomics helper functions ---------------------

# Several functions are used across several script for the bacterial
# proteogenomics pipeline



### Genomic coverage -----------------------------------------------------

# Function to draw rectangular venn diagram of the genomic coverage
plots_rectvenn <- function(
    ideo,
    ref,
    pep,
    set = c("Chromosome", "Protein-coding", "Expressed protein", "Detected reference peptide", "Detected novel peptide"),
    colour = NULL) {
    
    set <- match.arg(
        arg = set,
        choices = c("Chromosome", "Protein-coding", "Expressed protein", "Detected reference peptide", "Detected novel peptide"),
        several.ok = TRUE)
    
    if (is.null(colour)) {
        colour <- c("#a6cee3", "#1f78b4", "#b2df8a", "#33a02c", "#404040")
    }
    
    if (length(set) < length(colour)) {
        stop("Not enough colours defined!")
    }
    
    # Get all chromosome nucleotide position
    chrom_nuc <- gr_nucleotide_pos(
        grange = ideo)
    
    # Get all protein-coding associated nucleotide position
    coding_nuc <- gr_nucleotide_pos(
        grange = ref)
    
    # Get all expressed protein associated nucleotide position
    if (!"Expressed" %in% colnames(values(ref))) {
        expr_known <- pep$Proteins %>%
            strsplit(x = ., split = ", ", fixed = TRUE) %>%
            unlist(.) %>%
            unique(.) %>%
            data.frame(id = ., Expressed = TRUE) %>%
            dplyr::left_join(
                x = data.frame(id = names(ref)), y = ., by = "id") %>%
            dplyr::mutate(
                ., Expressed = ifelse(is.na(Expressed), FALSE, Expressed)) %>%
            set_rownames(.[["id"]]) %>%
            dplyr::select(., -id)
        values(ref) <- cbind(
            values(ref), expr_known)
    }
    exprs_nuc <- gr_nucleotide_pos(
        grange = ref, filter = "Expressed == TRUE")
    
    # Get all peptide associated nucleotide position
    cover_nuc_ref <- gr_nucleotide_pos(
        grange = pep, filter = 'grepl("Known|Target", Database)')
    cover_nuc_novel <- gr_nucleotide_pos(
        grange = pep, filter = 'grepl("^Novel$", Database)')
    
    set_list <- list(
        chrom_nuc,
        coding_nuc,
        exprs_nuc,
        cover_nuc_ref,
        cover_nuc_novel
    ) %>%
        set_names(set)
    
    # Plot a square venn diagram of chromosome coverage
    plot.new()
    for (i in 1:length(set_list)) {
        area <- sqrt(
            length(unique(unlist(set_list[i]))) /
                length(unique(unlist(set_list[["Chromosome"]]))))
        rect(
            xleft = (i/100),
            ybottom = (i/100),
            xright = area + (i/100),
            ytop = area + (i/100),
            border = colour[i],
            lwd = 2)
        text(
            x = (i/100),
            y = area + (i/100) + 0.02,
            labels = paste(
                names(set_list)[i],
                round(length(unique(unlist(set_list[i]))) / 1000000, 3),
                "Mb", sep = " "),
            col = colour[i], cex = 1.0, adj = 0)
    }
    
    # Recorded the plot
    pl <- recordPlot()
    
    # Return the plot result
    return(pl)
    
}

# Function to visualise nucleotide coverage depth for known and novel ORF
plot_nuc_coverage <- function(
    pep,
    count,
    set = c("Known|Target", "Novel"),
    colour = NULL) {
    
    set <- match.arg(
        arg = set,
        choices = c("Known|Target", "Novel"),
        several.ok = TRUE)
    
    if (is.null(colour)) {
        colour <- c("#387eb8", "#e21e25")
    }
    
    if (length(set) < length(colour)) {
        stop("Not enough colours defined!")
    }
    
    toplot <- data.frame(Count = as.character(c(0:50)))
    quantiles_toplot <- data.frame(
        Quartiles = c("0%", "25%", "50%", "75%", "100%", "Mean"))
    
    for (i in 1:length(set)) {
        
        # Get all peptide associated nucleotide position
        coverage_pep <- gr_nucleotide_pos(
            grange = pep, filter = paste0('grepl("', set[i], '", Database)'))
        
        # Convert to dataframe
        coverage_nucl <- coverage_pep %>%
            plyr::ldply(., "data.frame") %>%
            set_colnames(c("Sequence", "Nucl_pos")) %>%
            dplyr::filter(., !is.na(Nucl_pos))
        
        # Compute msms count for each peptide sequence
        coverage_nucl_count <- count %>%
            dplyr::select(., Sequence, `MS/MS count`) %>%
            dplyr::filter(., !is.na(`MS/MS count`)) %>%
            dplyr::group_by(., Sequence) %>%
            dplyr::summarise(., Count = sum(`MS/MS count`)) %>%
            dplyr::inner_join(coverage_nucl, ., by = "Sequence")
        coverage_nucl_count$Count <- factor(
            x = coverage_nucl_count$Count,
            levels = unique(coverage_nucl_count$Count),
            ordered = TRUE)
        
        # Calculate overall nucleotide coverage frequencies
        toplot <- data.frame(
            table(coverage_nucl_count$Count)) %>%
            set_colnames(c("Count", set[i])) %>%
            dplyr::full_join(x = toplot, y = ., by = "Count")
        
        # Calculate the quantiles of count frequencies
        quantiles_toplot_i <- quantile(
            x = as.integer(as.character(coverage_nucl_count$Count)),
            probs = c(0, 0.25, 0.5, 0.75, 1), na.rm = TRUE) %>%
            as.data.frame(.) %>%
            set_colnames("Values") %>%
            dplyr::mutate(., Quartiles = row.names(.))
        quantiles_toplot_i <- data.frame(
            Quartiles = "Mean",
            Values = base::mean(
                as.integer(as.character(coverage_nucl_count$Count)),
                na.rm = TRUE)) %>%
            dplyr::bind_rows(quantiles_toplot_i, .) %>%
            dplyr::select(., Quartiles, Values) %>%
            dplyr::mutate(., Values = round(x = Values, digits = 1)) %>%
            set_colnames(c("Quartiles", set[i]))
        
        quantiles_toplot <- dplyr::full_join(
            x = quantiles_toplot, y = quantiles_toplot_i, by = "Quartiles")
        
    }
    
    toplot %<>%
        dplyr::mutate_all(~tidyr::replace_na(data = ., replace = 0)) %>%
        tidyr::pivot_longer(
            data = ., cols = -Count, names_to = "Set", values_to = "Freq") %>%
        dplyr::mutate(
            ., Count = ifelse(
                as.integer(as.character(Count)) >= 50,
                "50+", as.character(Count))) %>%
        dplyr::group_by(., Count, Set) %>%
        dplyr::summarise(., Freq = sum(Freq)) %>%
        dplyr::ungroup(.)
    toplot$Count <- factor(
        x = toplot$Count, levels = c(0:49, "50+"), ordered = TRUE)
    
    ggplot(
        toplot,
        aes(x = Count, y = Freq, group = Set, fill = Set)) +
        geom_bar(stat = "identity", position = position_dodge(width = 0)) +
        ggpubr::theme_pubr() +
        xlab("MS/MS counts") +
        ylab("Nucleotide counts") +
        ggtitle("Coverage per nucleotide") +
        annotation_custom(
            grob = gridExtra::tableGrob(
                d = quantiles_toplot,
                theme = gridExtra::ttheme_minimal(), rows = NULL),
            xmin = 20, xmax = 50, ymin = 1E5, ymax = max(toplot$Freq)) +
        scale_x_discrete(
            breaks = c(seq(0, 49, 5), "50+"),
            limits = c(seq(0, 49, 1), "50+")) +
        scale_fill_manual(values = colour)
    
}

# Function to visualise location of expressed known and novel ORFs on circos
plot_circos <- function(
    ideo,
    ref,
    orf,
    pep,
    set = c("Known ORF (+)", "Known ORF (-)", "Reannotated ORF (+)", "Reannotated ORF (-)"),
    colour = NULL) {
    
    set <- match.arg(
        arg = set,
        choices = c("Known ORF (+)", "Known ORF (-)", "Reannotated ORF (+)", "Reannotated ORF (-)"),
        several.ok = TRUE)
    
    if (is.null(colour)) {
        colour <- c("#4682B4", "#BD5E5E", "#437A3C", "#F57C36")
    }
    
    if (length(set) < length(colour)) {
        stop("Not enough colours defined!")
    }
    
    # Get all expressed protein associated nucleotide position
    if (!"Expressed" %in% colnames(values(ref))) {
        expr_known <- subset(x = pep, grepl("Known|Target", Database)) %>%
            .$Proteins %>%
            strsplit(x = ., split = ", ", fixed = TRUE) %>%
            unlist(.) %>%
            unique(.) %>%
            data.frame(id = ., Expressed = TRUE) %>%
            dplyr::left_join(
                x = data.frame(id = names(ref)), y = ., by = "id") %>%
            dplyr::mutate(
                ., Expressed = ifelse(is.na(Expressed), FALSE, Expressed)) %>%
            set_rownames(.[["id"]]) %>%
            dplyr::select(., -id)
        values(ref) <- cbind(
            values(ref), expr_known)
    }
    if (!"Expressed" %in% colnames(values(orf))) {
        expr_known <- subset(x = pep, grepl("Novel", Database)) %>%
            .$Proteins %>%
            strsplit(x = ., split = ", ", fixed = TRUE) %>%
            unlist(.) %>%
            unique(.) %>%
            data.frame(id = ., Expressed = TRUE) %>%
            dplyr::left_join(
                x = data.frame(id = names(orf)), y = ., by = "id") %>%
            dplyr::mutate(
                ., Expressed = ifelse(is.na(Expressed), FALSE, Expressed)) %>%
            set_rownames(.[["id"]]) %>%
            dplyr::select(., -id)
        values(orf) <- cbind(
            values(orf), expr_known)
    }
    
    # Use ggbio extension to plot ORF location on genome as a circos graph
    ggplot() +
        #ggtitle(label = title) +
        ggbio::layout_circle(
            ideo, geom = "ideo", fill = "gray70",
            radius = 30, trackWidth = 2) +
        ggbio::layout_circle(
            ideo, geom = "scale", size = 4,
            radius = 33, trackWidth = 2) +
        ggbio::layout_circle(
            ideo, geom = "text", size = 7, aes(label = seqnames),
            angle = 0, radius = 37, trackWidth = 5) +
        ggbio::layout_circle(
            subset(x = ref, strand == "+" & Expressed),
            geom = "rect", color = colour[1],
            radius = 26, trackWidth = 4) +
        annotate(
            geom = "text", x = -6, y = 6, hjust = 0,
            label = set[1], colour = colour[1]) +
        ggbio::layout_circle(
            subset(x = ref, strand == "-" & Expressed),
            geom = "rect", color = colour[2],
            radius = 22, trackWidth = 4) +
        annotate(
            geom = "text", x = -6, y = 2, hjust = 0,
            label = set[2], colour = colour[2]) +
        ggbio::layout_circle(
            subset(x = orf, strand == "+" & Expressed),
            geom = "rect", color = colour[3],
            radius = 18, trackWidth = 4) +
        annotate(
            geom = "text", x = -6, y = -2, hjust = 0,
            label = set[3], colour = colour[3]) +
        ggbio::layout_circle(
            subset(x = orf, strand == "-" & Expressed),
            geom = "rect", color = colour[4],
            radius = 14, trackWidth = 4) +
        annotate(
            geom = "text", x = -6, y = -6, hjust = 0,
            label = set[4], colour = colour[4])
    
}



### Genomic ranges related functions -------------------------------------

# Function to get the nucleotide position for supplied GRange object
gr_nucleotide_pos <- function(grange = NULL, filter = NULL) {
    
    # Check whether the correct variables have been submitted by users
    if (is.null(grange)) {
        stop("Error: No GRange object was provided!")
    }
    
    # Apply the filter if it exist
    if (!is.null(filter)) {
        grange %<>%
            subset(x = ., subset = eval(parse(text = filter)))
    }
    
    # Get all nucleotide position for the GRange object
    #data <- lapply(X = ranges(grange), FUN = function(x) {
    #    x
    #})
    data <- data.frame(start = start(grange), end = end(grange))
    data <- apply(X = data, MARGIN = 1, FUN = function(x) {
        x[["start"]]:x[["end"]]
    })
    names(data) <- names(grange)
    
    # Return the nucleotide position for each feature
    return(data)
    
}

# Function to find preceding, following or nearest neighbour between
# genomic ranges and then compute their distance
gr_near_dist <- function(
    query = NULL,
    subject = NULL,
    near_type = c("nearest","precede", "follow"),
    select = NULL) {
    
    # Safety: defined parameters
    if (is.null(query)) {
        stop("No query genomic range data provided!")
    }
    if (is.null(subject)) {
        stop("No subject genomic range data provided!")
    }
    near_type %<>%
        match.arg(
            arg = .,
            choices = c("nearest","precede", "follow"),
            several.ok = TRUE)
    if (is.null(select)) {
        stop("No selection value for data provided!")
    }
    select %<>%
        match.arg(
            arg = .,
            choices = c("all","first", "last"),
            several.ok = FALSE)
    
    # List to compile data
    data_list <- list()
    
    # Loop through all function to apply
    for (fun in near_type) {
        
        # Call the specified function
        my_data <- do.call(
            what = fun,
            args = list(x = query, subject = subject, select = select))
        
        # Format data and compute distance between query and subject
        data_list[[fun]] <- my_data %>%
            as.data.frame(.) %>%
            dplyr::mutate(
                .,
                queryID = names(query)[queryHits],
                queryStrand = as.character(strand(query)[queryHits]),
                subjectID = names(subject)[subjectHits],
                subjectStrand = as.character(strand(subject)[subjectHits]),
                dist = GenomicRanges::distance(
                    query[queryID], subject[subjectID]))
        
    }
    
    # Return the results
    return(data_list)
    
}



### Packages import ------------------------------------------------------

# Create a function to load or install (then load) the required packages
load_package <- function(package) {
    
    stop("This function was discontinued to maintain package versionning!")
    
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
        requireNamespace("BiocManager", quietly = TRUE)
        BiocManager::install(pkgs = eval(package), update = FALSE)
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
    header = TRUE,
    columns = NULL,
    ...) {
    
    # Safety: defined parameters
    if (is.null(file)) {
        stop("No defined Blast file name!")
    }
    format <- match.arg(blast_format)
    if (format != "6") {
        warning(paste0(
            "The format is unsual: ", format,
            " data will be read as table, check results!"))
    }
    if (header == TRUE & is.null(columns)) {
        columns <- data.table::fread(
            input = file, sep = "\t", header = header, nrows = 0,
            stringsAsFactors = FALSE, integer64 = "double",
            quote = "", data.table = FALSE) %>%
            colnames(.)
    } else if (header == FALSE & is.null(columns)) {
        stop("The Blast column names must be specified from file header or user specification!")
    }
    
    # Read in the blast file
    blast_data <- fread(
        input = file, sep = "\t", header = header,
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
    bb_filter = NULL,
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
    key <- as.symbol(key)
    
    eps <- .Machine$double.xmin
    
    data.filt <- data %>%
        unique(.) %>%
        dplyr::mutate(
            ., ncover = qend - qstart + 1,
            pcover = ncover * 100 / qlen) %>%
        dplyr::group_by(., !!key) %>%
        dplyr::arrange(., evalue, dplyr::desc(score), dplyr::desc(pident)) %>%
        dplyr::mutate(
            ., 
            logevalue = -log((evalue+eps), 10),
            Devalue = (max(logevalue, na.rm = T)-logevalue)*100/max(logevalue, na.rm = T),
            Dbitscore = (max(bitscore, na.rm = T)-bitscore)*100/max(bitscore, na.rm = T),
            Dscore = (max(score, na.rm = T)-score)*100/max(score, na.rm = T),
            Dpident = (max(pident, na.rm = T)-pident)*100/max(pident, na.rm = T),
            Dnident = (max(nident, na.rm = T)-nident)*100/max(nident, na.rm = T),
            Dpcover = (max(pcover, na.rm = T)-pcover)*100/max(pcover, na.rm = T),
            Dncover = (max(ncover, na.rm = T)-ncover)*100/max(ncover, na.rm = T))
    
    # Filter the data with specific fields in sequence
    if (is.null(bb_filter)) {
        
        # Filter based on default
        data.filt %<>%
            dplyr::filter(., evalue == min(evalue)) %>%
            dplyr::filter(., score == max(score)) %>%
            dplyr::filter(., pident == max(pident)) %>%
            dplyr::mutate(., best_count = n()) %>%
            base::as.data.frame(., stringsAsFactors = FALSE)
        
    } else {
        
        # Evaluate the filter parameter, this is required when parameters
        # contains space and was submitted from command line
        bb_filter <- eval(bb_filter)
        
        # Filter based on user provided criteria
        for (x in bb_filter) {
            data.filt %<>%
                dplyr::filter_(., bb_filter)
        }
        data.filt %<>%
            dplyr::mutate(., best_count = n()) %>%
            base::as.data.frame(., stringsAsFactors = FALSE)
        
    }
    
    # Check what should be done with entries matching multiple times
    data.final <- data.filt
    if (length(unique(data[[key]])) == length(data.final[[key]])) {
        warning("All blast hits matched uniquely!")
    } else if (multiple == "remove") {
        data.final %<>%
            dplyr::filter(., best_count == 1)
        warning("All multiple blast hits were removed!")
    } else if (multiple == "uniquify") {
        orig <- data.final[[key]] %>%
            unique(.) %>%
            length(.)
        data.final %<>%
            dplyr::rowwise(.) %>%
            dplyr::filter(
                ., best_count == 1 |
                    grepl(sstart, Description) |
                    grepl(send, Description))
        retained <- data.final[[key]] %>%
            unique(.) %>%
            length(.)
        warning(paste(
            "About", (orig - retained), "entries could not be uniquified!"))
    } else if (multiple == "keep") {
        warning("There may still be entries with multiple blast hits!")
    }
    
    data.final %<>%
        dplyr::ungroup(.) %>%
        dplyr::arrange(., !!key)
    
    # Return the best blast data
    return(base::as.data.frame(data.final))
    
}



### Functions connected to generation of summarizedExperiment ------------

# Function to generate a phenodata file using exact column names
required_file_summExp <- function(my_data, my_col, input_path) {
    
    sample_name <- grep(my_col, colnames(my_data), value = TRUE) %>%
        grep("( peptides|___(1|2|3))$", ., value = TRUE, invert = TRUE)
    
    if (length(sample_name) == 0) {
        stop(paste0(
            "The columns '", j,
            "' are non-existing in the dataset: ",
            input_path))
    }
    
    data.table::data.table(
        SampleName = sample_name,
        Description = sample_name,
        SummExpSplit = my_col,
        stringsAsFactors = FALSE
    )
    
}

# Function to generate a single summarized experiment object
make_SummarizedExperiment <- function(name, assay, coldata) {
    
    # Import and format the phenodata
    my_pheno <- data.table::fread(
        input = coldata, sep = "\t", quote = "",
        stringsAsFactors = FALSE,
        colClasses = "character", header = TRUE)
    my_coldata <- my_pheno %>%
        dplyr::mutate(., id = SampleName) %>%
        tibble::column_to_rownames(.data = ., var = "id")
    
    # Collect expression from grange metadata or directly from file
    is_grange <- grepl("RDS$", assay)
    if (is_grange) {
        my_expr <- readRDS(assay)
        my_assay <- my_expr %>%
            GenomicRanges::values(.) %>%
            .[, c("id", my_pheno$SampleName)] %>%
            set_rownames(.[["id"]])
        my_assay$id <- NULL
        my_assay %<>%
            as.matrix(.) %>%
            list(.) %>%
            set_names(name)
        my_rowRanges <- my_expr
        GenomicRanges::values(my_rowRanges) <- my_rowRanges %>%
            GenomicRanges::values(.) %>%
            .[, !(colnames(.) %in% my_pheno$SampleName)] %>%
            .[, grep("___(1|2|3)$", colnames(.), invert = TRUE)]
        SummarizedExperiment(
            assays = my_assay, rowRanges = my_rowRanges, colData = my_coldata)
    } else {
        my_expr <- data.table::fread(
            input = assay, sep = "\t", quote = "",
            stringsAsFactors = FALSE,
            colClasses = "character", header = TRUE)
        my_assay <- my_expr %>%
            dplyr::select(., id, my_pheno$SampleName) %>%
            dplyr::mutate_at(., my_pheno$SampleName, as.double) %>%
            tibble::column_to_rownames(.data = ., var = "id") %>%
			as.matrix(.) %>%
            list(.) %>%
            set_names(name)
        SummarizedExperiment(
            assays = my_assay, colData = my_coldata)
    }
    
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
                    pattern = "^(REV__)?(gi.+|sp|tr|gb|ref)\\|(.+)\\|.*$",
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

# Function to warn user of NA values in specific columns and
# then filter them out
filter_na <- function(my_data, my_col) {
    if (any(is.na(my_data[[my_col]]))) {
        warning(paste(
            "There are ",
            sum(is.na(my_data[[my_col]])),
            "NA values in", my_col))
        my_data %<>%
            dplyr::filter(., !is.na(!!as.name(my_col)))
    }
    return(my_data)
}



### Computation functions ------------------------------------------------

# Function to identify the digestion rule,
# e.g. after which amino acid digestion occured
digestion_rule <- function(x, threshold = 0.05) {
    aa_freq <- table(x)
    my_res <- aa_freq[
        aa_freq > (sum(aa_freq) * threshold)] %>%
        names(.)
    my_res[!my_res %in% c("", "-")] %>%
        paste(., collapse = "|") %>%
        paste0("(", ., ")")
}

# Function that digest all protein sequences in fasta file and provide
# peptide position and sequence
prot_digest <- function(
    fasta = NULL,
    custom = "K|R",
    miss_cleav = c(0:2)) {
    
    # Check whether the fasta has been submitted by users
    if (is.null(fasta)) {
        stop("Error: No fasta was provided!")
    }
    
    # Check the class of entries in provided fasta
    if (lapply(fasta, class) %>% unlist(.) %>% unique(.) != "SeqFastaAA") {
        warning("Provided fasta entries are not of class 'SeqFastaAA'!")
    }
    
    # Digest all protein and store peptide position into list
    pep.pos <- cleavageRanges(
        x = fasta %>%
            as.character(.),
        custom = custom,
        missedCleavages = miss_cleav) %>% 
        set_names(names(fasta))
    
    # Replicate the proteins names per number of associated peptides
    tmp <- lapply(pep.pos, nrow) %>%
        rep(x = names(.), times = .)
    
    # Transform the list of peptide position to dataframe
    pep.pos %<>%
        do.call(rbind, .) %>%
        base::data.frame(Proteins = tmp, ., stringsAsFactors = FALSE)
    
    # Digest all protein and store peptide sequence into list
    pep.seq <- cleave(
        x = fasta %>%
            as.character(.),
        custom = custom,
        missedCleavages = miss_cleav,
        unique = FALSE) %>% 
        set_names(names(fasta))
    
    # Replicate the proteins names per number of associated peptides
    tmp <- lapply(pep.seq, length) %>%
        rep(x = names(.), times = .)
    
    # If the list of proteins names for pep.seq is identical to the one
    # for pep.pos then add peptide sequence to dataframe of peptide positions
    if (identical(x = tmp, y = pep.pos$Proteins)) {
        pep.list <- pep.seq %>%
            unlist(.) %>%
            base::cbind(pep.pos, Sequence = ., stringsAsFactors = FALSE)
    } else {
        stop("'pep.seq' and 'pep.pos' protein names list are not identical!")
    }
    
    # Check that the binding of peptide position and sequence was successfull
    if (any(nchar(pep.list$Sequence) != (pep.list$end - pep.list$start + 1))) {
        stop("Peptide seq length do not match peptide width (end - start)!")
    }
    
    # Return result
    return(pep.list)
    
}

# Function to parralelise the digestion of protein across different fasta files
prot_digest_foreach <- function(databases, custom, mc) {
    
    # In silico digest the fasta files
    all_pep <- foreach(
        x = databases, .inorder = T, .export = "prot_digest",
        .packages = c("magrittr", "cleaver")) %dopar%
        prot_digest(fasta = x, custom = custom, miss_cleav = mc)
    
    # Replicate the database names per entry
    tmp <- lapply(all_pep, nrow) %>%
        rep(x = names(databases), times = .)
    
    # Transform the list of peptide position to dataframe
    all_pep %<>%
        do.call(rbind, .) %>%
        base::data.frame(id = tmp, ., stringsAsFactors = FALSE)
    
    # Return the list of digested priteins
    return(all_pep)
    
}

# Function to identify which peptide are unique/common between databases
unique_to_database <- function(
    digest = NULL,
    pep = NULL) {
    
    # Check whether the data variables have been submitted by users
    if (is.null(digest)) {
        stop("Error: No digest was provided!")
    }
    if (length(unique(digest$id)) < 2) {
        stop("Error: Multiple named databases must be provided!")
    }
    if (is.null(pep)) {
        stop("Error: No peptide sequence vector was provided!")
    }
    
    # Filter out all peptides that were not identified
    digest %<>%
        dplyr::filter(., Sequence %in% pep)
    
    # Concatenate the database of origins for all identified peptides
    digest %>%
        dplyr::group_by(., Sequence) %>%
        dplyr::mutate(
            ., Dbuniqueness = paste(unique(id), collapse = ";")) %>%
        base::as.data.frame(., stringsAsFactors = FALSE)
    
}

# Locate start codon across genome
find_start_codon <- function(x) {
    res <- Biostrings::vmatchPattern(
        pattern = DNAString(x), subject = my_bsgeno)
    mcols(res) <- data.frame(Start_codon = rep(x = x, times = length(res)))
    res
}

# Function to determine the translation frame based on start coordinates
# for grange object (will need to create class)
calculate_frame <- function(gr) {
    gr_pos <- subset(gr, strand == "+")
    mcols(gr_pos) <- ((as.numeric(start(gr_pos)) + 2) / 3) %>%
        keep_decimals(.) %>%
        round(x = ((. + 0.33) * 3)) %>%
        cbind(mcols(gr_pos), Frame = .)
    gr_neg <- subset(gr, strand == "-")
    mcols(gr_neg) <- ((as.numeric(end(gr_neg)) + 2) / 3) %>%
        keep_decimals(.) %>%
        round(x = ((. + 0.33) * -3)) %>%
        cbind(mcols(gr_neg), Frame = .)
    c(gr_pos, gr_neg)
}

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
    data %<>%
        dplyr::rowwise(.) %>%
        dplyr::mutate(., frame = dplyr::case_when(
            strand == "+" ~ round(x = (
                (keep_decimals(((!!as.name(start) + 2) / 3)) + 0.33) * 3)),
            strand == "-" ~ round(x = (
                (keep_decimals(((genome_size[[sseqid]] - !!as.name(start) + 3) / 3)) + 0.33) * -3))
        )) %>%
        dplyr::ungroup(.)
    
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

# Function which compile the different novelty reason for any sequence that
# maps to multiple novel ORFs
novel_pep_classify <- function(
    x,
    coordinate,
    blast_ref,
    blast_all) {
    
    # Get the peptide position within protein
    coordinate_tmp <- coordinate %>%
        dplyr::filter(., Sequence == x)
    
    # Assess the sequence novelty for each combination of sequence to ORF
    my_rows <- c(min(1, nrow(coordinate_tmp)):nrow(coordinate_tmp))
    novel_final <- lapply(X = my_rows, FUN = function(i) {
        .novel_pep_classify(
            x = x, coord = coordinate_tmp[i, ],
            blast_ref = blast_ref, blast_all = blast_all)
    }) %>%
        set_names(coordinate_tmp$Proteins) %>%
        plyr::ldply(
            ., "data.frame", .id = "Proteins", stringsAsFactors = FALSE)
    
    # Define empty comment
    warn <- ""
    
    # Check if a unique peptide location is available
    if (nrow(coordinate_tmp) == 0) {
        
        warn <- "Peptide not located"
        
    } else if (
        nrow(coordinate_tmp) > 1 &
        length(unique(novel_final$NoveltyReason)) > 1) {
        
        warn <- "Mapping to several novel ORFs (different novelty reasons)"
        
    } else if (
        nrow(coordinate_tmp) > 1 &
        length(unique(novel_final$NoveltyReason)) == 1) {
        
        warn <- "Mapping to several novel ORFs (same novelty reasons)"
        
    }
    
    # Return the reason for novelty
    novel_final %<>%
        dplyr::bind_cols(., warning = rep(warn, times = nrow(coordinate_tmp)))
    novel_final
    
}

# Function which compile the different novelty reason obtained for a sequence
# and then defines the final reason (using reference and best blast)
.novel_pep_classify <- function(
    x,
    coord,
    blast_ref,
    blast_all) {
    
    # Get the blast and levenshtein data for current ORF to novel peptide map
    blast_ref_tmp <- blast_ref[
        blast_ref$qseqid == coord$Proteins, ]
    blast_all_tmp <- blast_all[
        blast_all$qseqid == coord$Proteins, ]
    
    # Get the best blast(s) for current ORF to novel peptide map
    best_blast_all_tmp <- blast_all_tmp %>%
        unique(.) %>%
        dplyr::group_by(., qseqid)
    
    if (nrow(best_blast_all_tmp) > 0) {
        best_blast_all_tmp %<>%
            dplyr::filter(., evalue_reciproc == min(evalue_reciproc)) %>%
            dplyr::filter(., score_reciproc == max(score_reciproc)) %>%
            dplyr::filter(., pident_reciproc == max(pident_reciproc))
    }

    # Perform novelty analysis on reference blast
    my_rows <- c(min(1, nrow(blast_ref_tmp)):nrow(blast_ref_tmp))
    ref_novel <- lapply(
        X = my_rows, FUN = function(i) {
            .novelty(x = x, coord = coord, blast = blast_ref_tmp[i, ])
        }) %>%
        plyr::ldply(., "data.frame", .id = NULL, stringsAsFactors = FALSE) %>%
        dplyr::bind_cols(Sequence = rep(x, times = nrow(.)), .)
    
    # Perform novelty analysis on best blast (only if best different from ref)
    if (identical(blast_ref_tmp$sseqid, best_blast_all_tmp$sseqid)) {
        best_novel <- ref_novel
    } else {
        my_rows <- c(min(1, nrow(best_blast_all_tmp)):nrow(best_blast_all_tmp))
        best_novel <- lapply(
            X = my_rows, FUN = function(i) {
                .novelty(x = x, coord = coord, blast = best_blast_all_tmp[i, ])
            }) %>%
            plyr::ldply(
                ., "data.frame", .id = NULL, stringsAsFactors = FALSE) %>%
            dplyr::bind_cols(Sequence = rep(x, times = nrow(.)), .)
    }
    
    # Combine reference and best novelty results
    novel_res <- dplyr::full_join(
        x = ref_novel, y = best_novel, by = "Sequence",
        suffix = c("_ref", "_best"))
    
    # Keep record of entries with multi blast results
    my_comment <- c()
    if (length(unique(novel_res$blast_best)) > 1) {
        my_comment <- c(my_comment, "Multiple best blast hits")
    }
    if (length(unique(novel_res$blast_ref)) > 1) {
        my_comment <- c(my_comment, "Multiple ref blast hits")
    }
    novel_res$comment <- paste(my_comment, collapse = ";")
    
    # Check sequence having multiple different reasons
    if (
        length(unique(novel_res$reason_best)) > 1 |
        length(unique(novel_res$reason_ref)) > 1) {
        novel_res$NoveltyReason <- "Multi incompatible reasons"
    } else {
        novel_res$NoveltyReason <- ""
    }
    
    # Summarise the sequence novelty info into a single row
    novel_res %<>%
        dplyr::group_by(., Sequence) %>%
        dplyr::summarise_all(
            .tbl = ., .funs = funs(paste(unique(.), collapse = ";")))
    
    # Finalise the novelty explanation by merging the explanation
    if (novel_res$NoveltyReason == "") {
        
        if (novel_res$reason_ref == novel_res$reason_best) {
            novel_res$NoveltyReason <- novel_res$reason_ref
        } else if (novel_res$reason_best == "Not novel") {
            novel_res$NoveltyReason <- paste(
                novel_res$reason_ref, "(exact match elsewhere)")
        } else if (novel_res$reason_best %in% c("SAVs/InDels", "SAV")) {
            novel_res$NoveltyReason <- paste(
                novel_res$reason_ref, "(SAV match elsewhere)")
        } else if (novel_res$reason_best %in% c("Potential alternate end", "Potential alternate start")) {
            novel_res$NoveltyReason <- paste(
                novel_res$reason_ref, "(partial match elsewhere)")
        } else if (novel_res$reason_best %in% c("Incompatible end site", "Incompatible start site")) {
            novel_res$NoveltyReason <- novel_res$reason_ref
        }
        
    }
    
    # Return the novelty reason
    novel_res
    
}

# Function to interprete peptide novelty based on blast and levenshtein data
.novelty <- function(
    x,
    coord,
    blast) {
    
    # Default reason is undetermined
    reason <- "Undetermined reason"
    
    # Dataframe compiling novelty reason for current entry
    my_res <- data.frame(
        reason = reason,
        blast = NA_character_,
        leven = NA_character_,
        stringsAsFactors = FALSE)
    
    # Check if blast results are available
    if (nrow(blast) == 0) {
        
        # Without blast data mark as potentially novel peptide
        reason <- "Potentially novel"
        
        # Return the results
        my_res$reason <- reason
        return(my_res)
        
    } else {
        
        # Concatenate the blast data
        my_res$blast <- paste(
            blast$sseqid, blast$evalue_reciproc,
            blast$score_reciproc, blast$pident_reciproc,
            blast$Description, blast$Taxon,
            sep = "||")
        
        # Check whether the peptide is within the matching blast positions
        if (
            coord$start >= blast$qstart_blast &
            coord$end <= blast$qend_blast) {
            
            # If yes then mark as potential SAV
            reason <- "Potential SAV"
            
            # Compute the levenshtein distance
            leven <- adist(
                x = x, y = blast$sseq_blast,
                partial = TRUE, ignore.case = TRUE) %>%
                as.data.frame(.) %>%
                cbind(blast$sseqid, .) %>%
                set_colnames(c("id", "leven"))
            
            # Concatenate the levenshtein data
            my_res$leven <- paste(leven$id, leven$leven, sep = "/")
            
            # Check if the levenshtein data exist and
            # what is the distance score
            if (
                nrow(leven) > 0 &
                unique(leven$leven) == 1) {
                
                # If yes mark as confirmed SAV
                reason <- "SAV"
                
            } else if (
                nrow(leven) > 0 &
                unique(leven$leven) == 0) {
                
                # If yes not novel peptide
                reason <- "Not novel"
                
            } else if (
                nrow(leven) > 0 &
                unique(leven$leven) > 1) {
                
                # If yes multi SAVs or InDels
                reason <- "SAVs/InDels"
                
            }
            
            # Return the results
            my_res$reason <- reason
            return(my_res)
            
            # Check whether the peptide overlap the reference protein start
        } else if (min(coord$start, blast$qstart_blast) == coord$start) {
            
            # Mark peptide as coming from an alternate start
            reason <- "Potential alternate start"
            
            # Check whether novel ORF and reference protein have
            # incompatible start site
            if (blast$sstart_blast > blast$qstart_blast) {
                
                reason <- "Incompatible start site"
                
            }
            
            # Return the results
            my_res$reason <- reason
            return(my_res)
            
            # Check whether the peptide overlap the reference protein end
        } else if (max(coord$end, blast$qend_blast) == coord$end) {
            
            # Mark peptide as coming from an alternate end
            reason <- "Potential alternate end"
            
            # Check whether novel ORF and reference protein have
            # incompatible end site
            if (blast$send_blast > blast$qend_blast) {
                
                reason <- "Incompatible end site"
                
            }
            
            # Return the results
            my_res$reason <- reason
            return(my_res)
            
        }
        
    }
    
    # Return the reason result
    my_res
    
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

# Function to describe a dataset via user-specified columns,
# it calculates a count and then plots the results
describe_count_plot <- function(
    my_data, main_id, x, separator = NULL, fill = NULL,
    xlabel = "Name", ylabel = "Count", gtitle = NULL,
    threshold = 5, flip = FALSE) {
    
    first_grouping <- c(x, fill)
    second_grouping <- c("Name", fill)
    
    if (is.null(separator)) {
        my_data_format <- my_data
    } else {
        my_data_format <- my_data %>%
            tidyr::separate_rows(data = ., !!as.name(x), sep = separator)
    }
    
    my_data_format %<>%
        dplyr::group_by_at(., vars(tidyselect::all_of(first_grouping))) %>%
        dplyr::summarise(
            ., Count = dplyr::n_distinct(!!as.name(main_id))) %>%
        dplyr::group_by(., !!as.name(x)) %>%
        dplyr::mutate(., Name = ifelse(
            any(Count >= threshold),
            !!as.name(x),
            paste0("Others (<=", threshold, ")"))) %>%
        dplyr::group_by_at(., vars(tidyselect::all_of(second_grouping))) %>%
        dplyr::summarise(., Count = sum(Count)) %>%
        dplyr::ungroup(.)
    
    if (!is.null(fill)) {
        my_data_format %<>%
            tidyr::complete(
                data = ., Name, !!as.name(fill),
                fill = list(Count = 0))
    }
    
    if (all(c("ordered", "factor") %in% class(my_data[[x]]))) {
        
        my_data_format$Name <- factor(
            x = my_data_format$Name,
            levels = c(
                levels(my_data[[x]]),
                paste0("Others (<=", threshold, ")")),
            ordered = TRUE)
        
    } else {
        
        my_data_format %<>%
            dplyr::arrange(., dplyr::desc(Count))
        
        my_data_format$Name <- factor(
            x = my_data_format$Name,
            levels = unique(my_data_format$Name),
            ordered = TRUE)
        
    }
    
    if (is.null(fill)) {
        pl <- ggplot(
            my_data_format,
            aes(x = Name, y = Count, label = Count))
    } else {
        pl <- ggplot(
            my_data_format,
            aes(x = Name, y = Count, fill = !!as.name(fill), label = Count))
    }
    
    pl <- pl +
        geom_bar(stat = "identity", position = "dodge") +
        ggpubr::theme_pubr() +
        ggtitle(gtitle) +
        xlab(xlabel) +
        ylab(ylabel)
    
    if (flip) {
        pl + 
            geom_text(position = position_dodge(width = 0.9), hjust = -0.4) +
            coord_flip()
    } else {
        pl + 
            geom_text(position = position_dodge(width = 0.9), vjust = -0.4)
    }
    
}

# Function to generate custom histogram for upset attribute plots
# the key change is the deletion of duplicated entries so that stack bars
# actually represent the data
custom_hist <- function(mydata, x, y) {
    
    # Remove duplicated entries (typically gray65) in case queries were used
    mydata_filt <- mydata %>%
        as.data.table(., key = x)
    mydata_format <- mydata_filt[
        , .(color = paste0(sort(color), collapse = ";")), by = c(x, y)]
    mydata_format$color %<>%
        sub("(;gray65|gray65;)", "", .)
    
    # Generate the histogram (stat bin)
    plot <- ggplot(
        data = mydata_format, aes(x = !!as.name(y), fill = color)) + 
        geom_histogram(stat = "bin", position = "stack", bins = 50) +
        scale_fill_identity() +
        theme_pubr()
    
}

# Function that computes the required values for boxplot and samples
# maximum 1000 outliers, the difference with the plots_box function
# comes from the value returned which is a ggplot (instead of ggtable)
# and is compatible with plotly
plots_box_plotly <- function(my_data, x, y, col_manual, txt_size) {
    
    # Format the data, keep required columns and remove NA
    toplot <- my_data %>%
        dplyr::select(., !!as.name(x), !!as.name(y))
    toplot %<>%
        filter_na(my_data = ., my_col = y)
    toplot[[x]] <- factor(
        x = toplot[[x]], levels = names(col_manual), ordered = TRUE)
    
    # Calculate the boxplot values (q25, median, q75...)
    toplot %<>%
        dplyr::group_by(., !!as.name(x)) %>%
        dplyr::mutate(
            .,
            y25 = quantile(!!as.name(y), 0.25),
            y50 = median(!!as.name(y)),
            y75 = quantile(!!as.name(y), 0.75),
            iqr = y75 - y25,
            ymin = ifelse(
                min(!!as.name(y)) < (y25 - (1.5 * iqr)),
                (y25 - (1.5 * iqr)),
                min(!!as.name(y))),
            ymax = ifelse(
                max(!!as.name(y)) > (y75 + (1.5 * iqr)),
                (y75 + (1.5 * iqr)),
                max(!!as.name(y)))
        )
    toplot_box <- toplot %>%
        dplyr::select(., !!as.name(x), ymin, y25, y50, y75, ymax) %>%
        unique(.)
    
    # Identify the outliers (outside 1.5 iqr) and select only a 1000 per group
    toplot_sort <- toplot %>%
        dplyr::group_by(., !!as.name(x)) %>%
        dplyr::filter(
            ., !!as.name(y) < unique(ymin) | !!as.name(y) > unique(ymax)) %>%
        dplyr::group_by(., !!as.name(x)) %>%
        dplyr::arrange(., !!as.name(y))
    toplot_outliers <- toplot_sort %>%
        dplyr::slice(., 1, dplyr::n())
    toplot_outliers <- toplot_sort %>%
        dplyr::slice(., -1, -dplyr::n()) %>%
        dplyr::sample_n(
            .,
            size = ifelse(
                length(!!as.name(y)) < 1000, length(!!as.name(y)), 1000),
            replace = FALSE) %>%
        dplyr::bind_rows(., toplot_outliers) %>%
        dplyr::select(., !!as.name(x), !!as.name(y))
    
    # Generate the boxplot together with jittered outliers
    ggplot() +
        geom_boxplot(
            data = toplot_box, mapping = aes(
                x = !!as.name(x),
                ymin = ymin,
                lower = y25,
                middle = y50,
                upper = y75,
                ymax = ymax,
                fill = !!as.name(x),
                colour = !!as.name(x)),
            stat = "identity",
            alpha = 0.6,
            size = 1) +
        geom_jitter(
            data = toplot_outliers, mapping = aes(
                x = !!as.name(x),
                y = !!as.name(y),
                fill = !!as.name(x),
                colour = !!as.name(x)),
            position = position_jitter(height = 0, width = 0.1),
            stroke = 0.5) +
        xlab(x) +
        ylab(y) +
        labs(fill = NULL) +
        ggpubr::theme_pubr() +
        theme(text = element_text(size = txt_size)) +
        scale_fill_manual(values = col_manual) +
        scale_colour_manual(guide = FALSE, values = col_manual)
    
}

# Function that plots reference ORF, novel ORF, associated peptide and sanger
# sequence in a genomic range context around a target
plots_orf_genomic <- function(
    x = NULL,
    bsgeno,
    ref_gr = NULL,
    orf_gr = NULL,
    pep_gr = NULL,
    sanger_gr = NULL,
    region_range = 1000,
    track_colour = NULL,
    ref_label = NULL) {
    
    # Check whether the correct variables have been submitted by users
    if (is.null(x)) {
        stop("Error: No target novel ORF ID was provided!")
    }
    if (is.null(ref_gr)) {
        stop("Error: No reference GRange object was provided!")
    }
    if (is.null(orf_gr)) {
        stop("Error: No novel ORF GRange object was provided!")
    }
    if (is.null(pep_gr)) {
        stop("Error: No peptide GRange object was provided!")
    }
    if (is.null(sanger_gr)) {
        warning("No Sanger sequence GRange object was provided!")
    }
    if (!is.numeric(region_range)) {
        stop("Error: The genomic region range should be numeric!")
    }
    if (is.null(track_colour)) {
        warning("No track colours were provided!")
        track_colour <- c(
            "#4682B4", "#BD5E5E", "#437A3C", "#F57C36", "#D58DEB", "#B2B83F")
    }
    if (is.null(ref_label)) {
        stop("Error: No reference column name to use for labelling provided!")
    }
    
    # Define genomic region of interest
    start_region <- start(orf_gr[x]) - region_range
    end_region <- end(orf_gr[x]) + region_range
    seqname_region <- orf_gr[x]@seqnames@values
    
    # The list containing all plots
    pl_list <- list()
    
    # The plot of known genes for that genomic region
    tmp <- ref_gr[
        (start(ref_gr) %in% c(start_region:end_region) |
             end(ref_gr) %in% c(start_region:end_region)) &
            seqnames(ref_gr) == seqname_region] %>%
        as.data.frame(.)
    tmp$value <- rep(x = 1, times = nrow(tmp))
    colnames(tmp)[colnames(tmp) == ref_label] <- "label_name"
    
    if (nrow(tmp) > 0) {
        pl <- autoplot(
            ref_gr[
                (start(ref_gr) %in% c(start_region:end_region) |
                     end(ref_gr) %in% c(start_region:end_region)) &
                    seqnames(ref_gr) == seqname_region],
            mapping = aes(fill = strand, alpha = Expressed),
            geom = "arrowrect", layout = "linear", colour = "black") +
            geom_text(
                data = tmp,
                mapping = aes(
                    x = start + ((end - start) / 2),
                    y = value,
                    label = label_name),
                nudge_y = 0.45,
                check_overlap = TRUE) +
            scale_fill_manual(
                values = c(`+` = track_colour[1], `-` = track_colour[2])) +
            scale_alpha_manual(
                values = c("TRUE" = 1, "FALSE" = 0.5))
    } else {
        pl <- ggplot()
    }
    pl_list["Known ORFs"] <- list(pl)
    
    # Loop through all frames from the ORF GRange
    for (y in unique(orf_gr$frame)) {
        
        # Store the current orf grange
        tmp_orf <- orf_gr[
            (start(orf_gr) %in% c(start_region:end_region) |
                 end(orf_gr) %in% c(start_region:end_region)) &
                orf_gr$frame == y & seqnames(orf_gr) == seqname_region]
        
        # The plot of novel ORFs for that genomic region if possible
        if (length(tmp_orf) > 0) {
            pl <- autoplot(
                tmp_orf,
                mapping = aes(fill = strand, alpha = Expressed),
                geom = "arrowrect", layout = "linear", colour = "black") +
                scale_fill_manual(
                    values = c(`+` = track_colour[3], `-` = track_colour[4])) +
                scale_alpha_manual(
                    values = c("TRUE" = 1, "FALSE" = 0.5))
        } else {
            pl <- ggplot()
        }
        pl_list[paste("Frame", y)] <- list(pl)
        
    }
    
    # Define zoomed in genomic region of interest
    start_region_zoom <- start(orf_gr[x])
    end_region_zoom <- end(orf_gr[x])
    
    # Format peptides sequences for genomic region visualisation
    tmp <- pep_gr[
        (start(pep_gr) %in% c(start_region_zoom:end_region_zoom) |
             end(pep_gr) %in% c(start_region_zoom:end_region_zoom)) &
            seqnames(pep_gr) == seqname_region] %>%
        as.data.frame(.)
    tmp$value <- rep(x = 1, times = nrow(tmp))
    tmp %<>%
        dplyr::arrange(., frame) %>%
        dplyr::mutate(., Database = sub("(Known).*", "\\1", Database))
    tmp$id <- factor(
        x = tmp$id,
        levels = as.character(tmp$id),
        labels = as.character(tmp$id),
        ordered = TRUE)
    
    if (nrow(tmp) > 0) {
        
        # Loop through all frames for the peptide GRange
        tmp <- cbind(tmp, y = NA)
        y_value <- 0
        for (aes_group in unique(tmp$Database)) {
            
            # Define the subset of peptides for the current database
            sub_tmp <- tmp %>%
                dplyr::filter(., Database == aes_group)
            
            # Define the y axis value to start with
            y_value <- (y_value + 1)
            
            # Loop through the subset of peptide as long as all y values
            # are not all defined
            while (any(is.na(sub_tmp$y))) {
                
                # Define the minimum x value threshold to add peptide
                # on current y axis value
                if (length(
                    sub_tmp[
                        !is.na(sub_tmp$y) & sub_tmp$y == y_value,
                        "end"]) == 0) {
                    min_thresh <- 0
                } else {
                    min_thresh <- max(
                        sub_tmp[!is.na(sub_tmp$y) & sub_tmp$y == y_value, "end"])
                }
                
                # Check for entry closest to the defined threshold
                min_start <- min(sub_tmp[
                    is.na(sub_tmp$y) & sub_tmp$start > min_thresh, "start"])
                min_end <- min(sub_tmp[
                    is.na(sub_tmp$y) & sub_tmp$start == min_start, "end"])
                
                # Get the closest entry ID
                pep_id <- as.character(sub_tmp[
                    is.na(sub_tmp$y) &
                        sub_tmp$start == min_start &
                        sub_tmp$end == min_end, "id"])
                
                # If there are no closest entry, increase y axis value by 1
                # if there is give that entry the current y axis value
                if (length(pep_id) == 0) {
                    y_value <- (y_value + 1)
                } else if (length(pep_id) == 1) {
                    sub_tmp[sub_tmp$id == pep_id, "y"] <- y_value
                    tmp[tmp$id == pep_id, "y"] <- y_value
                } else {
                    stop("Number of matching peptide ID not allowed!")
                }
                
            }
            
        }
        
        # The plot of peptides sequences for that genomic region
        pl <- ggplot(
            data = tmp,
            mapping = aes(
                xmin = start,
                xmax = end,
                ymin = as.integer(y),
                ymax = (as.integer(y) + 0.5),
                fill = factor(Database),
                label = id)) +
            geom_rect(colour = "black") +
            geom_text(mapping = aes(
                x = (start + ((end - start) / 2)),
                y = (as.integer(y) + 0.275)),
                colour = "white",
                size = 3,
                check_overlap = TRUE) +
            scale_fill_manual(
                name = "Database",
                values = c(
                    `Known` = track_colour[5],
                    `Target` = track_colour[5],
                    `Novel` = track_colour[6]))
        
    } else {
        pl <- ggplot()
    }
    pl_list["Peptide"] <- list(pl)
    
    # The plot of Sanger sequences for that genomic region
    tmp <- sanger_gr[
        (start(sanger_gr) %in% c(start_region_zoom:end_region_zoom) |
             end(sanger_gr) %in% c(start_region_zoom:end_region_zoom)) &
            seqnames(sanger_gr) == seqname_region] %>%
        as.data.frame(.)
    if (nrow(tmp) > 0) {
        pl <- ggplot(
            data = tmp,
            mapping = aes(
                xmin = start,
                xmax = end,
                ymin = as.integer(factor(id)),
                ymax = (as.integer(factor(id)) + 0.5),
                fill = strand)) +
            geom_rect() +
            scale_fill_manual(
                name = "Strand",
                values = c(`+` = "#CDCDC1", `-` = "#8B8B83")) +
            labs(fill = "Strand")
    } else {
        pl <- ggplot()
    }
    pl_list["Sequenced PCR"] <- list(pl)
    
    # The plot of genome sequence for that genomic region
    wh <- orf_gr[x] %>%
        range(.)
    
    # Renaming of sequences to be compatible with autoplot (not tested) 
    seq_rename <- paste0("chr", 1:length(seqnames(bsgeno))) %>%
        set_names(seqnames(bsgeno))
    seqlevels(wh) <- seq_rename
    seqlevels(bsgeno) <- seq_rename
    pl <- autoplot(bsgeno, which = wh, geom = "rect")
    pl_list["Genome"] <- list(pl)

    # Return all plots and the coordinate of the genomic region
    return(list(
        plots = pl_list,
        region_coordinates = c(start = start_region, end = end_region),
        region_zoom = c(start = start_region_zoom, end = end_region_zoom)))
    
}

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
    if (bit64::is.integer64(data[[value]])) {
        
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
                x = as.name(key),
                y = as.name(value),
                group = as.name(group),
                fill = as.name(fill))) +
            geom_bar(
                stat = "identity",
                position = posit,
                colour = colour)
        
    } else {
        
        # Create the plot with fill as not aes
        plot <- ggplot(
            data = data,
            mapping = aes_string(
                x = as.name(key),
                y = as.name(value),
                group = as.name(group))) +
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

# Function to modify a ggplot object by removing some outliers data
# so that in can easily be plotted
outlier_sampling <- function(
    plot = NULL) {
    
    # Check whether the correct variables have been submitted by users
    if (is.null(plot)) {
        stop("Error: No plot was provided!")
    }
    
    # Build the ggplot for rendering
    gg_build <- ggplot_build(plot)
    
    # Loop through the outlier data
    edit_outlier <- list()
    for (x in c(1:length(gg_build$data[[1]]$outliers))) {
        
        # Check if outlier vector is longer than 2500 values
        if (length(gg_build$data[[1]]$outliers[[x]]) > 2500) {
            
            tmp <- c(
                min(gg_build$data[[1]]$outliers[[x]]),
                max(gg_build$data[[1]]$outliers[[x]]),
                sample(gg_build$data[[1]]$outliers[[x]], 2498))
            
            warning("Outliers data were sampled to reduce number of points!")
            
        } else {
            tmp <- gg_build$data[[1]]$outliers[[x]]
        }
        
        # Compile the edited outlier into a new list
        edit_outlier[[x]] <- tmp
        
    }
    
    # Loop through all data in ggplot_build
    orig_outlier <- gg_build$data[[1]]$outliers
    for (x in c(1:length(gg_build$data))) {
        
        if (!identical(orig_outlier, gg_build$data[[x]]$outliers)) {
            stop("ggplot_build data are not identical!")
        } else {
            gg_build$data[[x]]$outliers <- edit_outlier
        }
        
    }
    
    # Return the edited ggplot_build object as a gtable
    return(ggplot_gtable(gg_build))
    
}

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
    
    # Format numeric with thousand separator
    #labels <- formatC(x = data[[label]], format = 'd', big.mark = ',')
    
    # Add labels to the plot
    plot <- plot +
        geom_text(
            aes_(label = as.name(label)),
            position = posit,
            vjust = vjust,
            hjust = hjust,
            size = (textsize * 0.3),
            check_overlap = TRUE,
            angle = angle,
            parse = TRUE)
    
    # Return the updated plot
    return(plot)
    
}

# Create a function to make colour coding per condition
color_grouping <- function(
    data,
    filter,
    reserved = NULL) {
    
    # Make a vector of the condition and count unique
    map <- as.character(data[, filter])
    
    if (!is.null(reserved)) {
        Nmap <- length(unique(map)) + length(reserved)
    } else {
        Nmap <- length(unique(map))
    }
    
    # Check colour brewer set to use depending on number of conditions
    if (Nmap <= 9) {
        
        colours <- brewer.pal(
            n = Nmap, name = "Set1")
        
    } else if (Nmap > 9 && Nmap <= 12) {
        
        colours <- brewer.pal(
            n = length(unique(map)), name = "Set3")
        
    } else {
        
        colours <- rainbow(n = ((Nmap + 1) / 2))
        colours <- c(colours, rainbow(n = ((Nmap + 1) / 2), alpha = 0.3))
        
    }
    
    # Remove colours that are reserved
    if (!is.null(reserved)) {
        colours %<>%
            .[!(. %in% reserved)]
    }
    
    # Include a color code column
    map <- data.frame(condition = map, color = map, stringsAsFactors = FALSE)
    
    # Loop through the number of unique conditions
    for(x in 1:Nmap) {
        
        # Replace the condition value by its associated colour
        map$color <- gsub(
            pattern = paste(
                "^", as.character(unique(map$color)[x]), "$", sep = ""),
            replacement = colours[x], x = map$color)
    }
    
    # Return group colouring
    legend <- map[!duplicated(map), ]
    map <- map$color
    names(map) <- rownames(data)
    return(list(map = map, legend = legend))
}

# Function to collect ggplot legend obtained from website: 
# 'http://stackoverflow.com/questions/11883844/'
# 'inserting-a-table-under-the-legend-in-a-ggplot2-histogram'
gg_legend <- function(a.gplot) {
    tmp <- ggplot_gtable(ggplot_build(a.gplot))
    leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
    legend <- tmp$grobs[[leg]]
    return(legend)
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

# Function to modify the axis of an existing ggplot
conti_axis <- function(
    plot = NULL,
    axis = "none",
    zoom = NULL,
    transf = "none") {
    
    # Check whether the plot has been submitted by users
    if (is.null(plot)) {
        stop("Error: No plot was provided!")
    }
    
    # Check arguments for axis against allowed list
    axis <- match.arg(
        arg = axis,
        choices = c("x", "y"))
    
    # Check the format of the provided zoom
    if (!is.null(zoom) & (!is.numeric(zoom) | length(zoom) != 2)) {
        stop("The axis limits parameter is not a numeric vector of length 2!")
    }
    
    # Check arguments for transform against allowed list
    transf <- match.arg(
        arg = transf,
        choices = c(
            "none", "asn", "atanh", "boxcox", "exp", "identity", "log",
            "log10", "log1p", "log2", "logit", "probability", "probit",
            "reciprocal", "reverse", "sqrt"))
    
    # Check that at least one parameters is valid
    if (is.null(zoom) & transf == "none") {
        stop("No axis settings are valid!")
    }
    
    # If a zoom parameter has been submitted, define new plot settings
    if (!is.null(zoom)) {
        
        zoom.param <- paste(zoom, collapse = ", ") %>%
            paste("limits = c(", ., ")", sep = "")
        plot$labels$title <- paste(zoom, collapse = " / ") %>%
            paste(plot$labels$title, " (", ., ")", sep = "")
        
    } else {
        
        zoom.param <- NA
        
    }
    
    # If a transform parameter has been submitted, define new plot settings
    if (transf != "none") {
        
        transf.param <- transf %>%
            paste("trans = '", ., "'", sep = "")
        
    } else {
        
        transf.param <- NA
        
    }
    
    # Remove the parameters that were not defined
    params <- c(zoom.param, transf.param)
    params %<>%
        .[!is.na(.)]
    
    # Combine all new axis parameters into a single command
    scaleConti <- paste(params, collapse = ", ") %>%
        paste("scale_", axis, "_continuous(", ., ")", sep = "")
    
    # Add the new axis settings to the plot
    plot <- plot +
        eval(parse(text = scaleConti))
    
    # Return the updated plot
    return(plot)
    
}

# Determine number of row and column required to plot all figures
# on same plot with grid
grid_dim <- function(x) {
    
    if (!is.null(x) | is.vector(x) | is.integer(x)) {
        
        col <- 1
        for (row in 1:x) {
            if ((row * col) >= x) {
                break
            }
            else {
                col <- row
                if ((row * col) >= x) {
                    break
                }
            }
        }
        return(list(x = row, y = col))
        
    } else {
        
        stop("Submitted data is null or not an integer vector!")
        
    }
    
}

# Function that determine the number of plot per page and the number
# of page required to display all plots
grids_display <- function(
    myplot = NULL,
    grid.x = NULL,
    grid.y = NULL) {
    
    # Check whether the plot list has been submitted by users
    if (is.null(myplot)) {
        stop("Error: No plot list was provided!")
    }
    
    # Check whether the grid number of rows and columns were submitted by user
    if (is.null(grid.x) | is.null(grid.y)) {
        
        # Determine required number of rows and columns for the plot
        grids <- grid_dim(x = length(myplot))
        grid.x <- grids$x
        grid.y <- grids$y
        
    }
    
    # Determine the number of page to plot all graphs
    nmb.page <- length(myplot) / (grid.x * grid.y)
    
    # Loop through the total number of page required for plotting
    for (x in 1:nmb.page) {
        
        # Determine the index of the plot that will go on current page
        plot.y <- (x * (grid.x * grid.y))
        plot.x <- (plot.y + 1 - (grid.x * grid.y))
        
        # Check whether there are less plot than spot on page
        if (plot.y > length(myplot)) {
            
            plot.y <- length(myplot)
            
        }
        
        # Create a graph
        do.call(
            what = "grid.arrange",
            args = c(myplot[plot.x:plot.y], list(
                nrow = grid.x, ncol = grid.y)))
        
    }
    
}



### Markdown report functions --------------------------------------------

# Function to create ioslides report given any 
report_markdown <- function(
    rmd_file = NULL,
    params = "ask",
    format = "ioslides_presentation",
    ext = "html") {
    
    # Check formatting of the params parameter
    if (!is.list(params)) {
        warning("'params' not a named list will try to help you with GUI!")
        params <- "ask"
    }
    
    # Check existence of markdown report file
    if (!file.exists(rmd_file)) {
        stop("Markdown not existing!")
    }
    
    # Scan the master markdown file for childs markdown and
    # if any find their absolute path
    rmd_childs <- scan(
        file = rmd_file, what = "character",
        sep = "", comment.char = "") %>%
        grep(".+\\.rmd", ., value = TRUE) %>%
        lapply(X = ., FUN = function(x) {
            list.files(
                path = dirname(rmd_file),
                pattern = x,
                full.names = TRUE,
                recursive = TRUE)
        }) %>%
        unlist(.)
    
    # Create the temporary directory
    tmp_report <- file.path(tempdir())
    
    # Copy markdown files to temporary location
    file.copy(
        from = c(rmd_file, rmd_childs),
        to = tmp_report,
        overwrite = TRUE)
    
    # Define output file name
    out_file <- paste(
        getwd(),
        "/",
        format(Sys.time(), '%Y%m%d_%H-%M'),
        "_",
        tools::file_path_sans_ext(basename(rmd_file)),
        ".",
        ext,
        sep = "")
    
    # Render the markdown report
    rmarkdown::render(
        input = file.path(tmp_report, basename(rmd_file)),
        output_format = format,
        output_file = out_file,
        params = params,
        envir = new.env(parent = globalenv()))
    
}



### Overrepresentation functions -----------------------------------------

# Function to perform overrepresentation analysis using fgsea package
fora_scored <- function(
    pathways, genes, universe,
    minSize = 1, maxSize = Inf, pval = 1, padj = 1) {
    
    my_oa <- fgsea::fora(
        pathways = pathways, genes = genes, universe = universe,
        minSize = opt$minsize, maxSize = opt$maxsize)
    
    if (is.numeric(opt$pval)) {
        my_oa %<>%
            dplyr::filter(., pval <= opt$pval)
    }
    
    if (is.numeric(opt$padj)) {
        my_oa %<>%
            dplyr::filter(., padj <= opt$padj)
    }
    
    uni_size <- length(universe)
    gene_size <- length(genes)
    
    my_oa %<>%
        dplyr::mutate(
            ., `-log10_pval` = -log10(pval), `-log10_padj` = -log10(padj),
            universe = uni_size, genes = gene_size,
            expected = genes * (size / universe),
            representation = ifelse(overlap >= expected, "+", "-"),
            score = overlap / expected,
            overlapGenes = sapply(overlapGenes, toString)) %>%
        dplyr::select(
            ., pathway, overlap, size,
            representation, score, pval, padj,
            `-log10_pval`, `-log10_padj`, universe,
            genes, expected, overlapGenes) %>%
        dplyr::arrange(., pval)
    
}

# Function to perform gene-set enrichment analysis using fgsea package
fgsea_scored <- function(
    pathways, stats, minSize = 1, maxSize = Inf, pval = 1, padj = 1) {
    
    my_gsea <- fgsea::fgsea(
        pathways = pathways, stats = stats,
        minSize = opt$minsize, maxSize = opt$maxsize)
    
    if (is.numeric(opt$pval)) {
        my_gsea %<>%
            dplyr::filter(., pval <= opt$pval)
    }
    
    if (is.numeric(opt$padj)) {
        my_gsea %<>%
            dplyr::filter(., padj <= opt$padj)
    }
    
    stats_size <- length(stats)
    
    my_gsea %<>%
        dplyr::mutate(
            ., `-log10_pval` = -log10(pval), `-log10_padj` = -log10(padj),
            genes = stats_size,
            leadingEdge = sapply(leadingEdge, toString)) %>%
        dplyr::select(
            ., pathway, pval, padj, log2err, ES, NES, size,
            `-log10_pval`, `-log10_padj`, genes, leadingEdge) %>%
        dplyr::arrange(., pval)
    
}

# Function to format pathways/function to gene association
# from dataframe to list format
fgsea_pathways <- function(annotation, gene, resource) {
    
    my_data <- annotation %>%
        dplyr::select(., dplyr::all_of(c(gene, resource))) %>%
        tidyr::separate_rows(
            data = ., as.name(resource), sep = ";", convert = FALSE)
    
    my_data %<>%
        dplyr::filter(
            ., !is.na(!!as.name(resource)) & !!as.name(resource) != "")
    
    my_path <- base::split(x = my_data[[gene]], f = my_data[[resource]])
    
}


