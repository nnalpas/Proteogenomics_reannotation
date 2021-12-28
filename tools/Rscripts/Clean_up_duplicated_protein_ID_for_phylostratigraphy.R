


library(magrittr)

f <- "T:/User/Nicolas/Phylostratigraphy/Formatted/1232426.faa"
f <- "T:/User/Nicolas/Phylostratigraphy/Formatted/224911.faa"
f <- "T:/User/Nicolas/Phylostratigraphy/Formatted/243277.faa"

my_seqs <- seqinr::read.fasta(
    file = f, seqtype = "AA",
    as.string = TRUE, whole.header = TRUE)

table(duplicated(sub(" .+", "", names(my_seqs))))
table(duplicated(my_seqs))
table(duplicated(names(my_seqs)))
names(my_seqs)[duplicated(sub(" .+", "", names(my_seqs)))]

my_filt <- (!duplicated(sub(" .+", "", names(my_seqs))) & !duplicated(my_seqs))
my_seqs_filt <- my_seqs[my_filt]
names(my_seqs_filt)[duplicated(sub(" .+", "", names(my_seqs_filt)))] %<>%
    sub("\\.1 ", ".d ", .)

table(duplicated(sub(" .+", "", names(my_seqs_filt))))

#my_seqs_filt <- my_seqs[!duplicated(sub(" .+", "", names(my_seqs)))]

seqinr::write.fasta(
    sequences = my_seqs_filt, names = names(my_seqs_filt),
    file.out = f,
    open = "w", nbchar = 80, as.string = TRUE)


