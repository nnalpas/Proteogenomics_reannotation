


# Retrieve a list of all available genome
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/genbank/bacteria/assembly_summary.txt

# Filter the list of geneome based on specific taxon ID
awk -F '\t' 'FNR==NR{k[$1]=1;next;} FNR==1 || FNR==2 || k[$6] || k[$7]' taxonomy_result.txt assembly_summary.txt > assembly_summary_taxon.txt

# Only keep the complete genome
awk -F '\t' '{if($12=="Complete Genome") print}' assembly_summary_taxon.txt > assembly_summary_taxon_complete_genomes.txt

# Download all genomes
for path in $(cut -f 20 ./work/Synechocystis_6frame/Synteny/assembly_summary_taxon_complete_genomes.txt); do
	wget -R "*_from_*" -P ./work/Synechocystis_6frame/Synteny/ "${path}/*genomic.fna.gz"
done

# Extract all data
cd ./work/Synechocystis_6frame/Synteny/; gunzip ./*.gz; cd -;

# Collect the info of each downloaded assembly
for infile in `find /home/tu/tu_tu/tu_kxmna01/work/Synechocystis_6frame/Synteny -name "GCA*"`; do
	mypattern=`basename $infile | cut -f1,2 -d_`
	grep $mypattern assembly_summary_taxon_complete_genomes.txt | awk -v f=$infile '{print f "\t" $0}'  >> downloaded_assembly.txt
done

# Use vim to clean up the organism name (remove parentheses '()' and make all names unique (i.e. 'Synechococcus sp. WH 8101', 'Synechocystis sp. PCC 6803'). This can be automatised as well.

# Replace all space by underscore in organisms name
awk -F "\t" '{OFS = "\t"; gsub(/ /,"_",$9); print}' downloaded_assembly.txt > downloaded_assembly_format.txt

# Create the Newick/map file required for Cactus
echo -n "(" > evolverCyanobacteria.txt
cut -f9 downloaded_assembly_format.txt | sed ':a;N;$!ba;s/\n/:1.0,/g' >> evolverCyanobacteria.txt
#awk -F "\t" '{OFS=":1.0,"; print}'
sed -i "s/$/:1.0);/" evolverCyanobacteria.txt
echo "" >> evolverCyanobacteria.txt
awk -F '\t' '{print $9 " " $1}' downloaded_assembly_format.txt >> evolverCyanobacteria.txt


