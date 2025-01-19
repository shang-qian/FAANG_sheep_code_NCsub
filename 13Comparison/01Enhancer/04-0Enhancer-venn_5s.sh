#Comparison
workdir=~/FAANG/analysis
mkdir -p $workdir/12Comparison_sheep/01Enhancer/04Venn_all4Tis
cd $workdir/12Comparison_sheep/01Enhancer/04Venn_all4Tis

#####sort liftover output
for tissue in Cerebellum Cortex Lung Muscle
do
echo $tissue
mkdir $tissue
for species in cattle pig sheep mouse human
do
ls $workdir/12Comparison_sheep/01Enhancer/03Tissue_Liftover/$species/*${species}2human_${tissue}_liftover_output*.bed
sort -k1,1 -k2,2n $workdir/12Comparison_sheep/01Enhancer/03Tissue_Liftover/$species/*${species}2human_${tissue}_liftover_output*.bed >$tissue/01${species}_${tissue}_liftover8.enhancer

done
done


###barplot tissue cross-validate
workdir=~/FAANG/analysis
mkdir -p $workdir/12Comparison_sheep/01Enhancer/04Venn_all4Tis
cd $workdir/12Comparison_sheep/01Enhancer/04Venn_all4Tis

###barplot in each tissue cross validate by other four species 
for tissue in Lung Cerebellum Cortex Muscle
do
echo $tissue

sheep_enh=$tissue/01sheep_${tissue}_liftover8.enhancer
#cattle_enh1=$tissue/01cattle_${tissue}_liftover8.enhancer
#pig_enh=$tissue/01pig_${tissue}_liftover8.enhancer
#human_enh=$tissue/01human_${tissue}_liftover8.enhancer
#mouse_enh=$tissue/01mouse_${tissue}_liftover8.enhancer
ls $sheep_enh
Rscript ~/FAANG/scripts/12Comparison/01Enhancer/04-1Enh_peak_species_sheep.r $sheep_enh $tissue
done














###enhancer across five species Upset
####intervenn
export PATH=/mnt/ceph/bmurdoch/WGS_Dwarf_lambs/Shang/tenv/bin:$PATH
export PATH=/mnt/ceph/bmurdoch/WGS_Dwarf_lambs/Shang/Dwarf20211129/data/miniconda3/bin:$PATH

source activate /mnt/ceph/bmurdoch/WGS_Dwarf_lambs/Shang/tenv
workdir=~/20240401FAANG/analysis
mkdir -p $workdir/12Comparison_sheep/01Enhancer/04Venn_all4Tis
cd $workdir/12Comparison_sheep/01Enhancer/04Venn_all4Tis

for tissue in Cerebellum Cortex Lung Muscle
do
cd $tissue
for i in cattle pig human mouse 
do
awk '/^chr/ {print $1"\t"$2"\t"$3}' 02${tissue}${i}_peaks_Enh_sheep.bed > 05${tissue}_${i}.bed
done

awk '/^chr/ {print $0}' 01sheep_${tissue}_liftover8.enhancer >05${tissue}_sheep.bed

intervene upset -i 05${tissue}*.bed --filenames
intervene venn -i 05${tissue}*.bed --filenames
cd ..
done




