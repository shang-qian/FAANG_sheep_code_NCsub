#Comparison
workdir=~/FAANG/analysis
mkdir -p $workdir/12Comparison_sheep/02Promoter/04Venn_all4Tis
cd $workdir/12Comparison_sheep/02Promoter/04Venn_all4Tis

#####sort liftover output
for tissue in Cerebellum Cortex Lung Muscle
do
echo $tissue
mkdir $tissue
for species in cattle pig sheep mouse human
do
ls $workdir/12Comparison_sheep/02Promoter/03Tissue_Liftover/$species/*${species}2human_${tissue}_liftover_output*.bed
sort -k1,1 -k2,2n $workdir/12Comparison_sheep/02Promoter/03Tissue_Liftover/$species/*${species}2human_${tissue}_liftover_output*.bed >$tissue/01${species}_${tissue}_liftover8.promoter

done
done


###barplot tissue cross-validate
workdir=~/FAANG/analysis
mkdir -p $workdir/12Comparison_sheep/02Promoter/04Venn_all4Tis
cd $workdir/12Comparison_sheep/02Promoter/04Venn_all4Tis

###barplot in each tissue cross validate by other four species 
for tissue in Lung Cerebellum Cortex Muscle
do
echo $tissue

sheep_enh=$tissue/01sheep_${tissue}_liftover8.promoter
#cattle_enh1=$tissue/01cattle_${tissue}_liftover8.enhancer
#pig_enh=$tissue/01pig_${tissue}_liftover8.enhancer
#human_enh=$tissue/01human_${tissue}_liftover8.enhancer
#mouse_enh=$tissue/01mouse_${tissue}_liftover8.enhancer
ls $sheep_enh
Rscript ~/FAANG/scripts/12Comparison/02Promoter/04-1Prom_peak_species_sheep.r $sheep_enh $tissue
done


