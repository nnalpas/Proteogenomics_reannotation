


library(magrittr)

my_f <- "C:/Users/kxmna01/Desktop/CP023688_protein_FIXED.fasta"

my_fasta <- seqinr::read.fasta(file = my_f, seqtype = "AA", as.string = TRUE, whole.header = TRUE)

find_stop <- lapply(my_fasta, function(x) {
    grepl("\\*.", x)
}) %>%
    unlist(.)

my_fasta_filt <- my_fasta[!find_stop]

seqinr::write.fasta(
    sequences = my_fasta_filt, names = names(my_fasta_filt),
    file.out = sub(".fasta", "_FILT.fasta", my_f),
    open = "w", nbchar = 60, as.string = TRUE)

