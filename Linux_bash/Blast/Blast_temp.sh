
# Single ID retrieval and blasting against uniprot reference proteome for Bsu
blastdbcmd -db "/media/sf_F_DRIVE/Nicolas/Proteome/Bacillus_subtilis/Bsu_genome_assembly_GCA_000009045.1.out_FIXED_HEADER.fasta" -dbtype 'prot' -entry 'seq_149311' | blastp -query - -task 'blastp' -db "/media/sf_F_DRIVE/Nicolas/Proteome/Bacillus_subtilis/uniprot-proteome_Bacillus_subtilis_168_UP000001570_20150318.fasta" -out "blastp_test" -evalue 0.0001 -num_alignments 30 -num_threads 5 -outfmt '6 qseqid sseqid pident nident mismatch length gapopen qstart qend sstart send evalue bitscore score'

# Multi entries retrieval and blasting against uniprot reference proteome for Bsu
blastdbcmd -db "/media/sf_F_DRIVE/Nicolas/Proteome/Bacillus_subtilis/Bsu_genome_assembly_GCA_000009045.1.out_FIXED_HEADER.fasta" -dbtype 'prot' -entry_batch '/media/sf_E_DRIVE/processing/Nicolas/Data/Vaishnavi/combined - 6 frame translation/txt/Novel_pep_orfs.txt' | blastp -query - -task 'blastp' -db "/media/sf_F_DRIVE/Nicolas/Proteome/Bacillus_subtilis/uniprot-proteome_Bacillus_subtilis_168_UP000001570_20150318.fasta" -out "blastp_BsuRef_novel_candidates_21122016" -evalue 0.01 -num_alignments 30 -num_threads 5 -outfmt '6 qseqid sseqid pident nident mismatch length gapopen qstart qend sstart send evalue bitscore score'




./MakeBlastDb.sh '/media/sf_E_DRIVE/processing/Nicolas/Local_databases/Database_bacteria/uniprot_taxonomy_allBact/uniprot_bacteria.fasta'


# Multi entries retrieval and blasting against uniprot reference proteome for Bsu
blastdbcmd -db "/media/sf_F_DRIVE/Nicolas/Proteome/Bacillus_subtilis/Bsu_genome_assembly_GCA_000009045.1.out_FIXED_HEADER.fasta" -dbtype 'prot' -entry_batch '/media/sf_E_DRIVE/processing/Nicolas/Data/Vaishnavi/combined - 6 frame translation/txt/Novel_pep_orfs.txt' | blastp -query - -task 'blastp' -db "nr" -out "blastp_allprot_novel_candidates_21122016" -evalue 0.01 -num_alignments 30 -num_threads 5 -outfmt '6 qseqid sseqid pident nident mismatch length gapopen qstart qend sstart send evalue bitscore score'








blastp [-h] [-help] [-import_search_strategy filename]
    [-export_search_strategy filename] [-task task_name] [-db database_name]
    [-dbsize num_letters] [-gilist filename] [-seqidlist filename]
    [-negative_gilist filename] [-entrez_query entrez_query]
    [-db_soft_mask filtering_algorithm] [-db_hard_mask filtering_algorithm]
    [-subject subject_input_file] [-subject_loc range] [-query input_file]
    [-out output_file] [-evalue evalue] [-word_size int_value]
    [-gapopen open_penalty] [-gapextend extend_penalty]
    [-xdrop_ungap float_value] [-xdrop_gap float_value]
    [-xdrop_gap_final float_value] [-searchsp int_value]
    [-max_hsps_per_subject int_value] [-seg SEG_options]
    [-soft_masking soft_masking] [-matrix matrix_name]
    [-threshold float_value] [-culling_limit int_value]
    [-best_hit_overhang float_value] [-best_hit_score_edge float_value]
    [-window_size int_value] [-lcase_masking] [-query_loc range]
    [-parse_deflines] [-outfmt format] [-show_gis]
    [-num_descriptions int_value] [-num_alignments int_value] [-html]
    [-max_target_seqs num_sequences] [-num_threads int_value] [-ungapped]
    [-remote] [-comp_based_stats compo] [-use_sw_tback] [-version]

DESCRIPTION
   Protein-Protein BLAST 2.2.28+

OPTIONAL ARGUMENTS
 -h
   Print USAGE and DESCRIPTION;  ignore all other parameters
 -help
   Print USAGE, DESCRIPTION and ARGUMENTS; ignore all other parameters
 -version
   Print version number;  ignore other arguments

 *** Input query options
 -query <File_In>
   Input file name
   Default = `-'
 -query_loc <String>
   Location on the query sequence in 1-based offsets (Format: start-stop)

 *** General search options
 -task <String, Permissible values: 'blastp' 'blastp-short' 'deltablast' >
   Task to execute
   Default = `blastp'
 -db <String>
   BLAST database name
    * Incompatible with:  subject, subject_loc
 -out <File_Out>
   Output file name
   Default = `-'
 -evalue <Real>
   Expectation value (E) threshold for saving hits 
   Default = `10'
 -word_size <Integer, >=2>
   Word size for wordfinder algorithm
 -gapopen <Integer>
   Cost to open a gap
 -gapextend <Integer>
   Cost to extend a gap
 -matrix <String>
   Scoring matrix name (normally BLOSUM62)
 -threshold <Real, >=0>
   Minimum word score such that the word is added to the BLAST lookup table
 -comp_based_stats <String>
   Use composition-based statistics:
       D or d: default (equivalent to 2 )
       0 or F or f: No composition-based statistics
       1: Composition-based statistics as in NAR 29:2994-3005, 2001
       2 or T or t : Composition-based score adjustment as in Bioinformatics
   21:902-911,
       2005, conditioned on sequence properties
       3: Composition-based score adjustment as in Bioinformatics 21:902-911,
       2005, unconditionally
   Default = `2'

 *** BLAST-2-Sequences options
 -subject <File_In>
   Subject sequence(s) to search
    * Incompatible with:  db, gilist, seqidlist, negative_gilist,
   db_soft_mask, db_hard_mask
 -subject_loc <String>
   Location on the subject sequence in 1-based offsets (Format: start-stop)
    * Incompatible with:  db, gilist, seqidlist, negative_gilist,
   db_soft_mask, db_hard_mask, remote

 *** Formatting options
 -outfmt <String>
   alignment view options:
     0 = pairwise,
     1 = query-anchored showing identities,
     2 = query-anchored no identities,
     3 = flat query-anchored, show identities,
     4 = flat query-anchored, no identities,
     5 = XML Blast output,
     6 = tabular,
     7 = tabular with comment lines,
     8 = Text ASN.1,
     9 = Binary ASN.1,
    10 = Comma-separated values,
    11 = BLAST archive format (ASN.1) 
    
   
   Options 6, 7, and 10 can be additionally configured to produce
   a custom format specified by space delimited format specifiers.
   The supported format specifiers are:
   	    qseqid means Query Seq-id
   	       qgi means Query GI
   	      qacc means Query accesion
   	   qaccver means Query accesion.version
   	      qlen means Query sequence length
   	    sseqid means Subject Seq-id
   	 sallseqid means All subject Seq-id(s), separated by a ';'
   	       sgi means Subject GI
   	    sallgi means All subject GIs
   	      sacc means Subject accession
   	   saccver means Subject accession.version
   	   sallacc means All subject accessions
   	      slen means Subject sequence length
   	    qstart means Start of alignment in query
   	      qend means End of alignment in query
   	    sstart means Start of alignment in subject
   	      send means End of alignment in subject
   	      qseq means Aligned part of query sequence
   	      sseq means Aligned part of subject sequence
   	    evalue means Expect value
   	  bitscore means Bit score
   	     score means Raw score
   	    length means Alignment length
   	    pident means Percentage of identical matches
   	    nident means Number of identical matches
   	  mismatch means Number of mismatches
   	  positive means Number of positive-scoring matches
   	   gapopen means Number of gap openings
   	      gaps means Total number of gaps
   	      ppos means Percentage of positive-scoring matches
   	    frames means Query and subject frames separated by a '/'
   	    qframe means Query frame
   	    sframe means Subject frame
   	      btop means Blast traceback operations (BTOP)
   	   staxids means Subject Taxonomy ID(s), separated by a ';'
   	 sscinames means Subject Scientific Name(s), separated by a ';'
   	 scomnames means Subject Common Name(s), separated by a ';'
   	sblastnames means Subject Blast Name(s), separated by a ';'
   			 (in alphabetical order)
   	sskingdoms means Subject Super Kingdom(s), separated by a ';'
   			 (in alphabetical order) 
   	    stitle means Subject Title
   	salltitles means All Subject Title(s), separated by a '<>'
   	   sstrand means Subject Strand
   	     qcovs means Query Coverage Per Subject
   	   qcovhsp means Query Coverage Per HSP
   When not provided, the default value is:
   'qseqid sseqid pident length mismatch gapopen qstart qend sstart send
   evalue bitscore', which is equivalent to the keyword 'std'
   Default = `0'
 -show_gis
   Show NCBI GIs in deflines?
 -num_descriptions <Integer, >=0>
   Number of database sequences to show one-line descriptions for
   Not applicable for outfmt > 4
   Default = `500'
    * Incompatible with:  max_target_seqs
 -num_alignments <Integer, >=0>
   Number of database sequences to show alignments for
   Default = `250'
    * Incompatible with:  max_target_seqs
 -html
   Produce HTML output?

 *** Query filtering options
 -seg <String>
   Filter query sequence with SEG (Format: 'yes', 'window locut hicut', or
   'no' to disable)
   Default = `no'
 -soft_masking <Boolean>
   Apply filtering locations as soft masks
   Default = `false'
 -lcase_masking
   Use lower case filtering in query and subject sequence(s)?

 *** Restrict search or results
 -gilist <String>
   Restrict search of database to list of GI's
    * Incompatible with:  negative_gilist, seqidlist, remote, subject,
   subject_loc
 -seqidlist <String>
   Restrict search of database to list of SeqId's
    * Incompatible with:  gilist, negative_gilist, remote, subject,
   subject_loc
 -negative_gilist <String>
   Restrict search of database to everything except the listed GIs
    * Incompatible with:  gilist, seqidlist, remote, subject, subject_loc
 -entrez_query <String>
   Restrict search with the given Entrez query
    * Requires:  remote
 -db_soft_mask <String>
   Filtering algorithm ID to apply to the BLAST database as soft masking
    * Incompatible with:  db_hard_mask, subject, subject_loc
 -db_hard_mask <String>
   Filtering algorithm ID to apply to the BLAST database as hard masking
    * Incompatible with:  db_soft_mask, subject, subject_loc
 -culling_limit <Integer, >=0>
   If the query range of a hit is enveloped by that of at least this many
   higher-scoring hits, delete the hit
    * Incompatible with:  best_hit_overhang, best_hit_score_edge
 -best_hit_overhang <Real, (>=0 and =<0.5)>
   Best Hit algorithm overhang value (recommended value: 0.1)
    * Incompatible with:  culling_limit
 -best_hit_score_edge <Real, (>=0 and =<0.5)>
   Best Hit algorithm score edge value (recommended value: 0.1)
    * Incompatible with:  culling_limit
 -max_target_seqs <Integer, >=1>
   Maximum number of aligned sequences to keep 
   Not applicable for outfmt <= 4
   Default = `500'
    * Incompatible with:  num_descriptions, num_alignments

 *** Statistical options
 -dbsize <Int8>
   Effective length of the database 
 -searchsp <Int8, >=0>
   Effective length of the search space
 -max_hsps_per_subject <Integer, >=0>
   Override maximum number of HSPs per subject to save for ungapped searches
   (0 means do not override)
   Default = `0'

 *** Search strategy options
 -import_search_strategy <File_In>
   Search strategy to use
    * Incompatible with:  export_search_strategy
 -export_search_strategy <File_Out>
   File name to record the search strategy used
    * Incompatible with:  import_search_strategy

 *** Extension options
 -xdrop_ungap <Real>
   X-dropoff value (in bits) for ungapped extensions
 -xdrop_gap <Real>
   X-dropoff value (in bits) for preliminary gapped extensions
 -xdrop_gap_final <Real>
   X-dropoff value (in bits) for final gapped alignment
 -window_size <Integer, >=0>
   Multiple hits window size, use 0 to specify 1-hit algorithm
 -ungapped
   Perform ungapped alignment only?

 *** Miscellaneous options
 -parse_deflines
   Should the query and subject defline(s) be parsed?
 -num_threads <Integer, >=1>
   Number of threads (CPUs) to use in the BLAST search
   Default = `1'
    * Incompatible with:  remote
 -remote
   Execute search remotely?
    * Incompatible with:  gilist, seqidlist, negative_gilist, subject_loc,
   num_threads
 -use_sw_tback
   Compute locally optimal Smith-Waterman alignments?

