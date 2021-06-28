


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

# Use vim to clean up the organism name (remove parentheses '()' or '[]', slashes '/', equal '=' and make all names unique (i.e. 'Synechococcus sp. WH 8101', 'Synechocystis sp. PCC 6803'). This can be automatised as well.

# Replace all space by underscore in organisms name
awk -F "\t" '{OFS = "\t"; gsub(/ /,"_",$9); print}' downloaded_assembly.txt > downloaded_assembly_format.txt

# Create the Newick file required for Cactus
echo -n "(" > evolverCyanobacteria.txt
cut -f9 downloaded_assembly_format.txt | sed ':a;N;$!ba;s/\n/:1.0,/g' >> evolverCyanobacteria.txt
sed -i "s/$/:1.0);/" evolverCyanobacteria.txt
echo "" >> evolverCyanobacteria.txt

# Keep only the first sequence for each organism (make sure this is the chromosome sequence and not plasmid) and format the header within each fasta
cut -f1,9 downloaded_assembly_format.txt |
while IFS=$'\t' read -r -a myArray; do
	file=`basename ${myArray[0]} | sed "s/\\.gz$//"`
	outfile="`pwd`/${myArray[1]}.fasta"
	awk -v RS='>' 'NR>1 { gsub("\n", ";", $0); sub(";$", "", $0); print ">"$0 }' "${file}" | head -n 1 | tr ';' '\n' | grep -vE "^>" > ${myArray[1]}.fasta
	sed -i "1i >${myArray[1]}" "${myArray[1]}.fasta"
	echo "${myArray[1]} ${outfile}" >> evolverCyanobacteria.txt
done









# Attempt to run Cactus on my data
mkdir ./Synteny/TMP
nohup cactus ./Synteny/jobStore ./Synteny/evolverCyanobacteria.txt ./Synteny/evolverCyanobacteria.hal --stats --binariesMode local --logDebug --workDir ./Synteny/TMP --buildAvgs --defaultMemory 4Gi --defaultCores 2 --defaultDisk 100Gi --restart > ./Synteny/cactus.log 2>&1 &

nohup cactus ./Synteny_2021-03-02/jobStore ./Synteny_2021-03-02/evolverCyanobacteria.txt ./Synteny_2021-03-02/evolverCyanobacteria.hal --stats --binariesMode local --logDebug --workDir ./Synteny_2021-03-02/TMP --buildAvgs --defaultMemory 12Gi --defaultCores 2 --defaultDisk 100Gi > ./Synteny_2021-03-02/cactus.log 2>&1 &


nohup singularity exec ./Cactus_docker_1.3.0.sif cactus ./Synteny_2021-03-02/jobStore ./Synteny_2021-03-02/evolverCyanobacteria.txt ./Synteny_2021-03-02/evolverCyanobacteria.hal --stats --binariesMode local --logDebug --workDir ./Synteny_2021-03-02/TMP --buildAvgs > ./Synteny_2021-03-02/cactus.log 2>&1 &

nohup singularity exec ./Cactus_docker_1.3.0.sif cactus ./Synteny_2021-03-02/jobStore ./Synteny_2021-03-02/evolverMycobacterium.txt ./Synteny_2021-03-02/evolverMycobacterium.hal --stats --binariesMode local --logDebug --debugWorker --maxMemory 30Gi --maxDisk 600Gi --maxCores 6 --workDir ./Synteny_2021-03-02/TMP --buildAvgs > ./Synteny_2021-03-02/cactus.log 2>&1 &
nohup singularity exec ./Cactus_docker_1.3.0.sif cactus ./Synteny_2021-03-02/jobStore ./Synteny_2021-03-02/evolverMycobacterium.txt ./Synteny_2021-03-02/evolverMycobacterium.hal --stats --binariesMode local --logDebug --debugWorker --maxMemory 30Gi --maxDisk 600Gi --maxCores 6 --workDir ./Synteny_2021-03-02/TMP --buildAvgs --restart > ./Synteny_2021-03-02/cactus.log 2>&1 &


singularity exec /mnt/vol1000/Cactus_docker_1.3.0.sif cactus-prepare examples/evolverMammals.txt --outDir steps-output --outSeqFile steps-output/evolverMammals.txt --outHal steps-output/evolverMammals.hal --jobStore jobstore

