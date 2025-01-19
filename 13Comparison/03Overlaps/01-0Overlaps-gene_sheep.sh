#Comparison
workdir=~/FAANG/analysis
mkdir -p $workdir/12Comparison_sheep/03Overlap_gene
cd $workdir/12Comparison_sheep/03Overlap_gene


####
$workdir/12Comparison/01Enhancer/05Enhancer_gene/09sheep_common5_uniq_4tissues_all_enhancer.bed
$workdir/12Comparison/02Promoter/05Promoter_gene/09sheep_common5_uniq_4tissues_all_promoter.bed

F1=$workdir/12Comparison_sheep/01Enhancer/05Enhancer_gene/02Enh_gene/sheep_ruminant_uniq/09sheep_ruminant_uniq_4tissues_any1_enhancer.bed
F2=$workdir/12Comparison_sheep/02Promoter/05Promoter_gene/02Prom_gene/sheep_ruminant_uniq/09sheep_ruminant_uniq_4tissues_any1_Promoter.bed

awk 'NR==FNR {a[$1]; next} $1 in a' $F1 $F2



for i in Cerebellum Cortex Lung Muscle Brain
do
for j in sheep_uniq sheep_ruminant_uniq sheep_common5_uniq
do

F1=$workdir/12Comparison_sheep/01Enhancer/05Enhancer_gene/02Enh_gene/$j/*${j}_${i}_uniq_enhancer.bed
F2=$workdir/12Comparison_sheep/02Promoter/05Promoter_gene/02Prom_gene/$j/*${j}_${i}_uniq_Promoter.bed
ls $F1
ls $F2
awk 'NR==FNR {a[$1]; next} $1 in a' $F1 $F2  >${j}_${i}_enh_prom_both.gene
done
done




F1=$workdir/12Comparison_sheep/01Enhancer/05Enhancer_gene/02Enh_gene/sheep_uniq/07sheep_uniq_Brain_uniq_enhancer.bed
F2=$workdir/12Comparison_sheep/02Promoter/05Promoter_gene/02Prom_gene/sheep_uniq/07sheep_uniq_Brain_uniq_Promoter.bed
awk 'NR==FNR {a[$1]; next} $1 in a' $F1 $F2  >Brain_enh_prom_both.gene



for i in sheep_uniq sheep_ruminant_uniq sheep_common5_uniq
do
echo $i

F1=$workdir/12Comparison/01Enhancer/05Enhancer_gene/09${i}_4tissues_all_enhancer.bed
F2=$workdir/12Comparison/02Promoter/05Promoter_gene/09${i}_4tissues_all_promoter.bed


Rscript ~/FAANG/scripts/12Comparison/03Overlaps/01-1Overlaps_gene.r $F1 $F2 $i

done













