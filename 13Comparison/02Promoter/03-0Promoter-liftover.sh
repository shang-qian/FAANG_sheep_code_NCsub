#Comparison
workdir=~/FAANG/analysis
mkdir -p $workdir/12Comparison_sheep/02Promoter/03Tissue_Liftover
cd $workdir/12Comparison_sheep/02Promoter/03Tissue_Liftover

species="cattle"
mkdir $species
cd $species
for tissue in Cerebellum Cortex Lung Muscle
do
awk '{print $1"\t"$2"\t"$3}' $workdir/12Comparison/02Promoter/02Tissue_Promoter/02Prom_tissue/$species/02${species}_*${tissue}*.promoter \
> 01${species}_${tissue}_input_liftover.bed
liftOver -minMatch=0.8 01${species}_${tissue}_input_liftover.bed analysis/04Compare/chain/bosTau9ToHg38.over.chain.gz 02${species}2human_${tissue}_liftover_output.bed 03${species}_${tissue}_unMapped.bed
done
cd ..


workdir=~/FAANG/analysis
mkdir -p $workdir/12Comparison_sheep/02Promoter/03Tissue_Liftover
cd $workdir/12Comparison_sheep/02Promoter/03Tissue_Liftover

species="pig"
mkdir $species
cd $species
for tissue in Cerebellum Cortex Lung Muscle
do
awk '{print $1"\t"$2"\t"$3}' $workdir/12Comparison/02Promoter/02Tissue_Promoter/02Prom_tissue/$species/02${species}_*${tissue}*.promoter \
> 01${species}_${tissue}_input_liftover.bed
liftOver -minMatch=0.8 01${species}_${tissue}_input_liftover.bed analysis/04Compare/chain/susScr11ToHg38.over.chain.gz 02${species}2human_${tissue}_liftover_output.bed 03${species}_${tissue}_unMapped.bed
done
cd ..


workdir=~/FAANG/analysis
mkdir -p $workdir/12Comparison_sheep/02Promoter/03Tissue_Liftover
cd $workdir/12Comparison_sheep/02Promoter/03Tissue_Liftover

species="mouse"
mkdir $species
cd $species
for tissue in Cerebellum Cortex Lung Muscle
do
awk '{print $1"\t"$2"\t"$3}' $workdir/12Comparison/02Promoter/02Tissue_Promoter/02Prom_tissue/$species/02${species}_*${tissue}*.promoter \
> 01${species}_${tissue}_input_liftover.bed
liftOver -minMatch=0.8 01${species}_${tissue}_input_liftover.bed analysis/04Compare/chain/mm10ToHg38.over.chain.gz 02${species}2human_${tissue}_liftover_output.bed 03${species}_${tissue}_unMapped.bed
done
cd ..




######Sheep
workdir=~/FAANG/analysis
mkdir -p $workdir/12Comparison_sheep/02Promoter/03Tissue_Liftover
cd $workdir/12Comparison_sheep/02Promoter/03Tissue_Liftover

###sheep
species="sheep"
mkdir $species
cd $species
module load bowtie2 samtools bedtools
for tissue in Cerebellum Cortex Lung Muscle
do
echo $tissue
awk '{print $1"\t"$2"\t"$3}' $workdir/12Comparison_sheep/02Promoter/02Tissue_Promoter/$species/02${species}_*${tissue}*.promoter > 00${species}_${tissue}_promoter_Ramb2.bed
###bed to fastq
bedtools getfasta -fi /mnt/ceph/bmurdoch/ARS-UI_Ramb_v2.0/GCA_016772045.1_ARS-UI_Ramb_v2.0_genomic.fna -bed 00${species}_${tissue}_promoter_Ramb2.bed > 01${species}_${tissue}_promoter_Ramb2.seq
####OAR4
bowtie2 -p 32 -x ~/FAANG/analysis/12Comparison/01Enhancer/03Tissue_Liftover/sheep/genome/oar4_genome -f 01${species}_${tissue}_promoter_Ramb2.seq 2>> 02sheep_$tissue.promoter.mapping.log | samtools sort  - -o 03${species}_${tissue}_promoter.OAR4.sorted.bam
bedtools bamtobed -i 03${species}_${tissue}_promoter.OAR4.sorted.bam  > 04${species}_${tissue}_promoter_ovr4.bed
done


#
for tissue in Cerebellum Cortex Lung Muscle
do
sed -i "s/NC_019458.2/chr1/" 04${species}_${tissue}_promoter_ovr4.bed
sed -i "s/NC_019459.2/chr2/" 04${species}_${tissue}_promoter_ovr4.bed
sed -i "s/NC_019460.2/chr3/" 04${species}_${tissue}_promoter_ovr4.bed
sed -i "s/NC_019461.2/chr4/" 04${species}_${tissue}_promoter_ovr4.bed
sed -i "s/NC_019462.2/chr5/" 04${species}_${tissue}_promoter_ovr4.bed
sed -i "s/NC_019463.2/chr6/" 04${species}_${tissue}_promoter_ovr4.bed
sed -i "s/NC_019464.2/chr7/" 04${species}_${tissue}_promoter_ovr4.bed
sed -i "s/NC_019465.2/chr8/" 04${species}_${tissue}_promoter_ovr4.bed
sed -i "s/NC_019466.2/chr9/" 04${species}_${tissue}_promoter_ovr4.bed
sed -i "s/NC_019467.2/chr10/" 04${species}_${tissue}_promoter_ovr4.bed
sed -i "s/NC_019468.2/chr11/" 04${species}_${tissue}_promoter_ovr4.bed
sed -i "s/NC_019469.2/chr12/" 04${species}_${tissue}_promoter_ovr4.bed
sed -i "s/NC_019470.2/chr13/" 04${species}_${tissue}_promoter_ovr4.bed
sed -i "s/NC_019471.2/chr14/" 04${species}_${tissue}_promoter_ovr4.bed
sed -i "s/NC_019472.2/chr15/" 04${species}_${tissue}_promoter_ovr4.bed
sed -i "s/NC_019473.2/chr16/" 04${species}_${tissue}_promoter_ovr4.bed
sed -i "s/NC_019474.2/chr17/" 04${species}_${tissue}_promoter_ovr4.bed
sed -i "s/NC_019475.2/chr18/" 04${species}_${tissue}_promoter_ovr4.bed
sed -i "s/NC_019476.2/chr19/" 04${species}_${tissue}_promoter_ovr4.bed
sed -i "s/NC_019477.2/chr20/" 04${species}_${tissue}_promoter_ovr4.bed
sed -i "s/NC_019478.2/chr21/" 04${species}_${tissue}_promoter_ovr4.bed
sed -i "s/NC_019479.2/chr22/" 04${species}_${tissue}_promoter_ovr4.bed
sed -i "s/NC_019480.2/chr23/" 04${species}_${tissue}_promoter_ovr4.bed
sed -i "s/NC_019481.2/chr24/" 04${species}_${tissue}_promoter_ovr4.bed
sed -i "s/NC_019482.2/chr25/" 04${species}_${tissue}_promoter_ovr4.bed
sed -i "s/NC_019483.2/chr26/" 04${species}_${tissue}_promoter_ovr4.bed
sed -i "s/NC_019484.2/chrX/" 04${species}_${tissue}_promoter_ovr4.bed
sed -i "s/NC_001941.1/chrM/" 04${species}_${tissue}_promoter_ovr4.bed

liftOver -minMatch=0.8 04${species}_${tissue}_promoter_ovr4.bed analysis/04Compare/chain/oviAri4ToHg38.over.chain.gz 05${species}2human_${tissue}_liftover_output.bed 06${species}_${tissue}_unMapped.bed

done





workdir=~/FAANG/analysis
mkdir -p $workdir/12Comparison_sheep/02Promoter/03Tissue_Liftover
cd $workdir/12Comparison_sheep/02Promoter/03Tissue_Liftover

##human
species="human"
mkdir $species
cd $species
for tissue in Cerebellum Cortex Lung Muscle
do
awk '{print $1"\t"$2"\t"$3}' $workdir/12Comparison/02Promoter/02Tissue_Promoter/02Prom_tissue/$species/02${species}_*${tissue}*.promoter > 01${species}_${tissue}_input_liftover.bed
ln -sf 01${species}_${tissue}_input_liftover.bed 02human2human_${tissue}_liftover_output.bed 
done
cd ..






