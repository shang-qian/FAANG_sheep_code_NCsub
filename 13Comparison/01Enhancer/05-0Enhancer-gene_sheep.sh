#Comparison
workdir=~/FAANG/analysis
mkdir -p $workdir/12Comparison_sheep/01Enhancer/05Enhancer_gene
cd $workdir/12Comparison_sheep/01Enhancer/05Enhancer_gene

#####sort liftover output
module load bedtools
for tissue in Cerebellum Cortex Lung Muscle
do
echo $tissue
mkdir $tissue
cd $tissue

sort -k1,1 -k2,2n $workdir/12Comparison_sheep/01Enhancer/03Tissue_Liftover/sheep/05sheep2human_${tissue}_liftover_output8.bed >01${tissue}_Liftover_sort.bed
awk '/^chr/ {print $1"\t"$2"\t"$3}' ../../04Venn_all4Tis/$tissue/06${tissue}_sheep_uniq_enhancer.bed |sort -k1,1 -k2,2n |uniq >02-1${tissue}_sheep_uniq_enhancer.bed
awk '/^chr/ {print $1"\t"$2"\t"$3}' ../../04Venn_all4Tis/$tissue/07${tissue}_sheep_ruminant_uniq_enhancer.bed |sort -k1,1 -k2,2n |uniq >02-2${tissue}_sheep_ruminant_uniq_enhancer.bed 
awk '/^chr/ {print $1"\t"$2"\t"$3}' ../../04Venn_all4Tis/$tissue/08${tissue}_sheep_common5_uniq_enhancer.bed |sort -k1,1 -k2,2n |uniq >02-3${tissue}_sheep_common5_uniq_enhancer.bed

bedtools intersect -loj -a 02-1${tissue}_sheep_uniq_enhancer.bed -b 01${tissue}_Liftover_sort.bed |awk '$1==$4&&$2==$5+1&&$3==$6 {split($7,a,":");split(a[2],b,"-");print a[1]"_"b[1]"_"b[2]}' |sort |uniq >03-1${tissue}_sheep_uniq_enhancer.name
bedtools intersect -loj -a 02-2${tissue}_sheep_ruminant_uniq_enhancer.bed -b 01${tissue}_Liftover_sort.bed |awk '$1==$4&&$2==$5+1&&$3==$6 {split($7,a,":");split(a[2],b,"-");print a[1]"_"b[1]"_"b[2]}' |sort |uniq >03-2${tissue}_sheep_ruminant_uniq_enhancer.name
bedtools intersect -loj -a 02-3${tissue}_sheep_common5_uniq_enhancer.bed -b 01${tissue}_Liftover_sort.bed |awk '$1==$4&&$2==$5+1&&$3==$6 {split($7,a,":");split(a[2],b,"-");print a[1]"_"b[1]"_"b[2]}' |sort |uniq >03-3${tissue}_sheep_common5_uniq_enhancer.name

cd ..
done


####Enhancer gene
mkdir 02Enh_gene
cd 02Enh_gene
###Venn 
k=0
for i in sheep_uniq sheep_ruminant_uniq sheep_common5_uniq
do
echo $i

mkdir $i

let k=$k+1
F1=~/FAANG/analysis/12Comparison/01Enhancer/05Enhancer_gene/Cerebellum/03-${k}Cerebellum_${i}_enhancer.name
F2=~/FAANG/analysis/12Comparison/01Enhancer/05Enhancer_gene/Cortex/03-${k}Cortex_${i}_enhancer.name
F3=~/FAANG/analysis/12Comparison/01Enhancer/05Enhancer_gene/Lung/03-${k}Lung_${i}_enhancer.name
F4=~/FAANG/analysis/12Comparison/01Enhancer/05Enhancer_gene/Muscle/03-${k}Muscle_${i}_enhancer.name

F5=~/FAANG/analysis/10interaction/06Each_Tissue/Cerebellum/Cerebellum_01Enh_prom_total.txt
F6=~/FAANG/analysis/10interaction/06Each_Tissue/CerebralCortex/CerebralCortex_01Enh_prom_total.txt
F7=~/FAANG/analysis/10interaction/06Each_Tissue/Lung/Lung_01Enh_prom_total.txt
F8=~/FAANG/analysis/10interaction/06Each_Tissue/MuscleSM/MuscleSM_01Enh_prom_total.txt


Rscript ~/FAANG/scripts/12Comparison/01Enhancer/05-1Enh_5species_comparison.r $F1 $F2 $F3 $F4 $i $F5 $F6 $F7 $F8

done















