


# Peptide novelty
opt <- list(
    evidence = "C:\\Users\\kxmna01\\Dropbox\\Home_work_sync\\Work\\Srim_6frame\\Novel_res\\Group_evidence.RDS",
    reference_fasta = "C:\\Users\\kxmna01\\Dropbox\\Home_work_sync\\Work\\Srim_6frame\\Genome\\Streptomyces_rimosus_allStrains_2019-02-13.fasta",
    reciprocal_blast_ref = "C:\\Users\\kxmna01\\Dropbox\\Home_work_sync\\Work\\Srim_6frame\\ReciprocalBlast\\Best_Reciproc_Blast_Refprot_vs_ORFprot_annot",
    reciprocal_blast_uniprot = "C:\\Users\\kxmna01\\Dropbox\\Home_work_sync\\Work\\Srim_6frame\\ReciprocalBlast\\Best_Reciproc_Blast_Uniprot_vs_ORFprot_annot",
    reciprocal_blast_ncbi = "C:\\Users\\kxmna01\\Dropbox\\Home_work_sync\\Work\\Srim_6frame\\ReciprocalBlast\\Best_Reciproc_Blast_NCBIprot_vs_ORFprot_annot",
    peptide_location = "C:\\Users\\kxmna01\\Dropbox\\Home_work_sync\\Work\\Srim_6frame\\Novel_res\\Peptides_location.RDS",
    threads = 2,
    output = "NoveltyExplain"
)

# ORF novelty
opt <- list(
    novel_reason = "C:\\Users\\kxmna01\\Desktop\\NoveltyExplain\\Sequence_novelty_reason.RDS",
    ref_grange = "C:\\Users\\kxmna01\\Desktop/GRanges/Ref_prot_grange.RDS",
    orf_grange = "C:\\Users\\kxmna01\\Desktop/GRanges/Orf_prot_grange.RDS",
    operon_grange = character(0),
    pep_pos = "C:\\Users\\kxmna01\\Dropbox\\Home_work_sync\\Work\\Srim_6frame\\Novel_res\\Peptides_location.RDS",
    pep_grange = "C:\\Users\\kxmna01\\Dropbox\\Home_work_sync\\Work\\Srim_6frame\\GRanges\\Pept_seq_grange.RDS",
    sanger_grange = character(0),
    genome_grange = "C:\\Users\\kxmna01\\Dropbox\\Home_work_sync\\Work\\Srim_6frame\\GRanges\\Genome_grange.RDS",
    add_rbs = "",
    pep_class = c('class 1', 'class 2'),
    threads = 2,
    output = "NoveltyExplain"
)

# ref grange
opt <- list(
    coordinates = "C:\\Users\\kxmna01\\Dropbox\\Home_work_sync\\Work\\Srim_6frame\\ProtPosition\\Ref_prot_coordinates.txt",
    genome = "C:\\Users\\kxmna01\\Dropbox\\Home_work_sync\\Work\\Srim_6frame\\Genome\\Srim_G7_assembly.fasta",
    geno_name = "062019",
    circular = FALSE,
    annotations = "C:\\Users\\kxmna01\\Desktop\\ProtAnnotation\\Ref_prot_annotations.txt",
    output = "./GRanges/Ref_prot_grange.RDS"
)

# orf grange
opt <- list(
    coordinates = "C:\\Users\\kxmna01\\Dropbox\\Home_work_sync\\Work\\Srim_6frame\\ProtPosition\\Orf_prot_coordinates.txt",
    genome = "C:\\Users\\kxmna01\\Dropbox\\Home_work_sync\\Work\\Srim_6frame\\Genome\\Srim_G7_assembly.fasta",
    geno_name = "062019",
    circular = FALSE,
    annotations = "C:\\Users\\kxmna01\\Desktop\\ProtAnnotation\\Orf_prot_annotations.txt",
    output = "./GRanges/Orf_prot_grange.RDS"
)

# strathclyde grange
opt <- list(
    coordinates = "C:\\Users\\kxmna01\\Dropbox\\Home_work_sync\\Work\\Srim_6frame\\ProtPosition\\Strath_prot_coordinates.txt",
    genome = "C:\\Users\\kxmna01\\Dropbox\\Home_work_sync\\Work\\Srim_6frame\\Genome\\Srim_G7_assembly.fasta",
    geno_name = "062019",
    circular = FALSE,
    annotations = "",
    output = "C:\\Users\\kxmna01\\Dropbox\\Home_work_sync\\Work\\Srim_6frame/GRanges/Strath_prot_grange.RDS"
)


# Peptide novelty
opt <- list(
    evidence = "C:\\Users\\kxmna01\\Dropbox\\Home_work_sync\\Work\\Srim_6frame\\Novel_res\\Group_evidence.RDS",
    reference_fasta = "C:\\Users\\kxmna01\\Dropbox\\Home_work_sync\\Work\\Srim_6frame\\Genome\\Streptomyces_rimosus_allStrains_2019-02-13.fasta",
    reciprocal_blast_ref = "C:\\Users\\kxmna01\\Dropbox\\Home_work_sync\\Work\\Srim_6frame\\ReciprocalBlast\\Best_Reciproc_Blast_Refprot_vs_ORFprot_annot",
    reciprocal_blast_uniprot = "C:\\Users\\kxmna01\\Dropbox\\Home_work_sync\\Work\\Srim_6frame\\ReciprocalBlast\\Best_Reciproc_Blast_Uniprot_vs_ORFprot_annot",
    reciprocal_blast_ncbi = "C:\\Users\\kxmna01\\Dropbox\\Home_work_sync\\Work\\Srim_6frame\\ReciprocalBlast\\Best_Reciproc_Blast_NCBIprot_vs_ORFprot_annot",
    reciprocal_blast_strath = "C:\\Users\\kxmna01\\Dropbox\\Home_work_sync\\Work\\Srim_6frame\\ReciprocalBlast\\Best_Reciproc_Blast_Strathprot_vs_ORFprot_annot",
    peptide_location = "C:\\Users\\kxmna01\\Dropbox\\Home_work_sync\\Work\\Srim_6frame\\Novel_res\\Peptides_location.RDS",
    threads = 2,
    output = "NoveltyExplain"
)
# Import the reciprocal blast best hits between novel ORFs and
# all uniprot proteins
reciprocal_blast_strath <- opt$reciprocal_blast_strath %>%
    as.character(.) %>%
    read.table(
        file = ., header = TRUE, sep = "\t", quote = "",
        as.is = TRUE, comment.char = "") %>%
    dplyr::mutate(
        ., staxid_blast = as.integer(staxid_blast),
        staxid_reciproc = as.integer(staxid_reciproc))
# Compile all reciprocal best hits
reciprocal_blast_all <- dplyr::bind_rows(
    data.frame(DB = "Reference", reciprocal_blast_ref, stringsAsFactors = F),
    data.frame(DB = "UniProt", reciprocal_blast_uniprot, stringsAsFactors = F),
    data.frame(DB = "NCBI", reciprocal_blast_ncbi, stringsAsFactors = F),
    data.frame(DB = "Strathclyde", reciprocal_blast_strath, stringsAsFactors = F))

# ORF novelty
opt <- list(
    novel_reason = "C:\\Users\\kxmna01\\Dropbox\\Home_work_sync\\Work\\Srim_6frame\\NoveltyExplain\\Sequence_novelty_reason.RDS",
    ref_grange = "C:\\Users\\kxmna01\\Dropbox\\Home_work_sync\\Work\\Srim_6frame\\GRanges/Ref_prot_grange.RDS",
    strath_grange = "C:\\Users\\kxmna01\\Dropbox\\Home_work_sync\\Work\\Srim_6frame\\GRanges/Strath_prot_grange.RDS",
    orf_grange = "C:\\Users\\kxmna01\\Dropbox\\Home_work_sync\\Work\\Srim_6frame\\GRanges/Orf_prot_grange.RDS",
    operon_grange = character(0),
    pep_pos = "C:\\Users\\kxmna01\\Dropbox\\Home_work_sync\\Work\\Srim_6frame\\Novel_res\\Peptides_location.RDS",
    pep_grange = "C:\\Users\\kxmna01\\Dropbox\\Home_work_sync\\Work\\Srim_6frame\\GRanges\\Pept_seq_grange.RDS",
    sanger_grange = character(0),
    genome_grange = "C:\\Users\\kxmna01\\Dropbox\\Home_work_sync\\Work\\Srim_6frame\\GRanges\\Genome_grange.RDS",
    add_rbs = "",
    pep_class = c('class 1', 'class 2'),
    threads = 2,
    output = "NoveltyExplain"
)
# Import the reference genomic ranges file
strath_grange <- opt$strath_grange %>%
    as.character(.) %>%
    readRDS(file = .)
# Include the entries expression pattern into the GRange metadata
strath_grange_expr <- strath_grange
values(strath_grange_expr) <- cbind(
    values(strath_grange_expr), 
    Expressed = rep(FALSE, times = nrow(values(strath_grange_expr))))
# Function that plots reference ORF, novel ORF, associated peptide and sanger
# sequence in a genomic range context around a target
plots_orf_genomic <- function(
    x = NULL,
    bsgeno,
    ref_gr = NULL,
    strath_gr,
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
    pl_list["Known ORFs"] <- list(pl)
    
    # The plot of known genes for that genomic region
    tmp <- strath_gr[
        (start(strath_gr) %in% c(start_region:end_region) |
             end(strath_gr) %in% c(start_region:end_region)) &
            seqnames(strath_gr) == seqname_region] %>%
        as.data.frame(.)
    tmp$value <- rep(x = 1, times = nrow(tmp))
    colnames(tmp)[colnames(tmp) == ref_label] <- "label_name"
    pl <- autoplot(
        strath_gr[
            (start(strath_gr) %in% c(start_region:end_region) |
                 end(strath_gr) %in% c(start_region:end_region)) &
                seqnames(strath_gr) == seqname_region],
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
    pl_list["Strathclyde ORFs"] <- list(pl)
    
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
    
    # 
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
                values = c(`Known` = track_colour[5], `Novel` = track_colour[6]))
        
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
    pl <- autoplot(bsgeno, which = wh, geom = "rect")
    pl_list["Genome"] <- list(pl)
    
    # Return all plots and the coordinate of the genomic region
    return(list(
        plots = pl_list,
        region_coordinates = c(start = start_region, end = end_region),
        region_zoom = c(start = start_region_zoom, end = end_region_zoom)))
    
}


