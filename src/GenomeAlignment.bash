


module load lastz/1.04.03

lastz ./Synechocystis_6frame/genbank/Chromosome.fasta ./Synechocystis_6frame/Conservation/GCA_000012525.1_ASM1252v1_genomic_chromosome.fna --notransition --step=20 --nogapped --format=maf --output=./Synechocystis_6frame/Conservation/Synechocystis_sp_PCC_6803_vs_Synechococcus_elongatus_PCC_7942.maf --rdotplot=./Synechocystis_6frame/Conservation/Synechocystis_sp_PCC_6803_vs_Synechococcus_elongatus_PCC_7942_for_plotting


lastz ./Synechocystis_6frame/genbank/Chromosome.fasta ./Synechocystis_6frame/Conservation/GCA_000012525.1_ASM1252v1_genomic_chromosome.fna --format=maf --output=./Synechocystis_6frame/Conservation/Synechocystis_sp_PCC_6803_vs_Synechococcus_elongatus_PCC_7942.maf --rdotplot=./Synechocystis_6frame/Conservation/Synechocystis_sp_PCC_6803_vs_Synechococcus_elongatus_PCC_7942_for_plotting


cp ./Synechocystis_6frame/genbank/Chromosome.fasta ./Synechocystis_6frame/genbank/Chromosome_duplicate.fasta
lastz ./Synechocystis_6frame/genbank/Chromosome.fasta ./Synechocystis_6frame/genbank/Chromosome_duplicate.fasta --format=maf --output=./Synechocystis_6frame/Conservation/Synechocystis_sp_PCC_6803_vs_itself.maf --rdotplot=./Synechocystis_6frame/Conservation/Synechocystis_sp_PCC_6803_vs_itself_for_plotting




lastz ./Synechocystis_6frame/Conservation/GCA_000010065.1_ASM1006v1_genomic.fna ./Synechocystis_6frame/Conservation/GCA_000012525.1_ASM1252v1_genomic_chromosome.fna --format=maf --output=./Synechocystis_6frame/Conservation/Synechococcus_elongatus_PCC_6301_vs_Synechococcus_elongatus_PCC_7942.maf --rdotplot=./Synechocystis_6frame/Conservation/Synechococcus_elongatus_PCC_6301_vs_Synechococcus_elongatus_PCC_7942_for_plotting





lastz ./lastz_test/bovis_3601.fasta ./lastz_test/bcg_1.fasta --notransition --step=20 --nogapped --format=maf --output=./lastz_test/bovis_vs_bcg.maf --rdotplot=./lastz_test/bovis_vs_bcg_for_plotting

lastz ./Conservation/Synechocystis_sp._PCC_6803_substr._GT-G.fasta ./Conservation/Synechocystis_sp._PCC_6803_ASM34078v1.fasta --notransition --step=20 --nogapped --format=maf --output=./Conservation/GT-G_vs_ASM34078v1.maf --rdotplot=./Conservation/GT-G_vs_ASM34078v1_for_plotting
lastz ./Conservation/Synechocystis_sp._IPPAS_B-1465.fasta ./Conservation/Synechocystis_sp._PCC_6803_ASM34078v1.fasta --notransition --step=20 --nogapped --format=maf --output=./Conservation/IPPAS_B-1465_vs_ASM34078v1.maf --rdotplot=./Conservation/IPPAS_B-1465_vs_ASM34078v1_for_plotting
lastz ./Conservation/GCA_000010065.1_ASM1006v1_genomic.fna ./Conservation/Synechocystis_sp._PCC_6803_ASM34078v1.fasta --notransition --step=20 --nogapped --format=maf --output=./Conservation/elongatus_vs_ASM34078v1.maf --rdotplot=./Conservation/elongatus_vs_ASM34078v1_for_plotting
