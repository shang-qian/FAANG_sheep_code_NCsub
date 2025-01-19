#Comparison
workdir=~/FAANG/analysis
mkdir -p $workdir/12Comparison_sheep/01Enhancer/03Tissue_Liftover
cd $workdir/12Comparison_sheep/01Enhancer/03Tissue_Liftover

species="cattle"
mkdir $species
cd $species
for tissue in Cerebellum Cortex Lung Muscle
do
awk '{print $1"\t"$2"\t"$3}' $workdir/12Comparison_sheep/01Enhancer/02Tissue_Enhancer/02Enh_tissue/$species/02${species}_*${tissue}*.enhancer > 01${species}_${tissue}_input_liftover.bed
liftOver -minMatch=0.8 01${species}_${tissue}_input_liftover.bed 04Compare/chain/bosTau9ToHg38.over.chain.gz 02${species}2human_${tissue}_liftover_output8.bed 03${species}_${tissue}_unMapped.bed
done
cd ..


workdir=~/FAANG/analysis
mkdir -p $workdir/12Comparison_sheep/01Enhancer/03Tissue_Liftover
cd $workdir/12Comparison_sheep/01Enhancer/03Tissue_Liftover

species="pig"
mkdir $species
cd $species
for tissue in Cerebellum Cortex Lung Muscle
do
awk '{print $1"\t"$2"\t"$3}' $workdir/12Comparison_sheep/01Enhancer/02Tissue_Enhancer/02Enh_tissue/$species/02${species}_*${tissue}*.enhancer > 01${species}_${tissue}_input_liftover.bed
liftOver -minMatch=0.8 01${species}_${tissue}_input_liftover.bed 04Compare/chain/susScr11ToHg38.over.chain.gz 02${species}2human_${tissue}_liftover_output8.bed 03${species}_${tissue}_unMapped.bed
done
cd ..


workdir=~/FAANG/analysis
mkdir -p $workdir/12Comparison_sheep/01Enhancer/03Tissue_Liftover
cd $workdir/12Comparison_sheep/01Enhancer/03Tissue_Liftover

species="mouse"
mkdir $species
cd $species
for tissue in Cerebellum Cortex Lung Muscle
do
awk '{print $1"\t"$2"\t"$3}' $workdir/12Comparison_sheep/01Enhancer/02Tissue_Enhancer/02Enh_tissue/$species/02${species}_*${tissue}*.enhancer > 01${species}_${tissue}_input_liftover.bed
liftOver -minMatch=0.8 01${species}_${tissue}_input_liftover.bed 04Compare/chain/mm10ToHg38.over.chain.gz 02${species}2human_${tissue}_liftover_output8.bed 03${species}_${tissue}_unMapped.bed
done
cd ..



workdir=~/FAANG/analysis
mkdir -p $workdir/12Comparison_sheep/01Enhancer/03Tissue_Liftover
cd $workdir/12Comparison_sheep/01Enhancer/03Tissue_Liftover

###sheep
species="sheep"
mkdir $species
cd $species
for tissue in Cerebellum Cortex Lung Muscle
do
echo $tissue

liftOver -minMatch=0.8 $workdir/12Comparison/01Enhancer/03Tissue_Liftover/$species/04${species}_${tissue}_enhancer_ovr4.bed /mnt/ceph/bmurdoch/WGS_Dwarf_lambs/Shang/ATAC-seq20220203/comparison/analysis/04Compare/chain/oviAri4ToHg38.over.chain.gz 05${species}2human_${tissue}_liftover_output8.bed 06${species}_${tissue}_unMapped.bed
done
cd ..


workdir=~/FAANG/analysis
mkdir -p $workdir/12Comparison_sheep/01Enhancer/03Tissue_Liftover
cd $workdir/12Comparison_sheep/01Enhancer/03Tissue_Liftover

##human
species="human"
mkdir $species
cd $species
for tissue in Cerebellum Cortex Lung Muscle
do
awk '{print $1"\t"$2"\t"$3}' $workdir/12Comparison_sheep/01Enhancer/02Tissue_Enhancer/02Enh_tissue/$species/02${species}_*${tissue}*.enhancer > 01${species}_${tissue}_input_liftover.bed
ln -sf 01${species}_${tissue}_input_liftover.bed 02human2human_${tissue}_liftover_output.bed 
done
cd ..






