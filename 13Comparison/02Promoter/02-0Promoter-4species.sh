#Comparison
workdir=~/FAANG/analysis
mkdir -p $workdir/12Comparison_sheep/02Promoter/02Tissue_Promoter
cd $workdir/12Comparison_sheep/02Promoter/02Tissue_Promoter

##1 merge peak
mkdir 01Prom_peak
cd 01Prom_peak

for species in "cattle" "pig"
do
for i in Cerebellum Cortex Lung Muscle
do
H3K4me3=$workdir/12Comparison/01Enhancer/01Public_RE/01Public_data/${species}-peak/*H3K4me3*${i}*
H3K27ac=$workdir/12Comparison/01Enhancer/01Public_RE/01Public_data/${species}-peak/*H3K27ac*${i}*
ls $H3K4me3
ls $H3K27ac
Rscript ~/FAANG/scripts/12Comparison/02Promoter/02-1Prom_peak_merge.r $H3K4me3 $H3K27ac $species $i
done
done


for species in "human" "mouse"
do
echo $species
for i in cerebellum cortex lung muscle
do
H3K4me3=$workdir/12Comparison/01Enhancer/01Public_RE/01Public_data/${species}-peak/*${i}*H3K4me3*
H3K27ac=$workdir/12Comparison/01Enhancer/01Public_RE/01Public_data/${species}-peak/*${i}*H3K27ac*
ls $H3K4me3
ls $H3K27ac
Rscript ~/FAANG/scripts/12Comparison/02Promoter/02-1Prom_peak_merge.r $H3K4me3 $H3K27ac $species $i
done
done



###2. generate bed file for four species and four tissues
workdir=~/FAANG/analysis
mkdir -p $workdir/12Comparison/02Promoter/02Tissue_Promoter/01Prom_peak
cd $workdir/12Comparison/02Promoter/02Tissue_Promoter/01Prom_peak

for species in "human" "mouse" "cattle" "pig"
do
echo $species
for tissue in "Cerebellum" "Cortex" "Lung" "Muscle"
do
tissue_name="${tissue:1}"
awk '/^chr/ {print $0}' $species/${species}_peaks_Promoter_*${tissue_name}.bed > $species/${species}_peaks_Prom_${tissue}_chr.bed

done
done



 

#####2 extract the tissue's enhancer using merge peak and potential 
workdir=~/FAANG/analysis
mkdir -p $workdir/12Comparison/02Promoter/02Tissue_Promoter
cd $workdir/12Comparison/02Promoter/02Tissue_Promoter

mkdir 02Prom_tissue
cd 02Prom_tissue

for species in human mouse pig cattle
do
mkdir $species
for tissue in Cerebellum Cortex Lung Muscle
do
echo $tissue
peakf=$workdir/12Comparison/02Promoter/02Tissue_Promoter/01Prom_peak/${species}/${species}_peaks_Prom_${tissue}_chr.bed
prom_region=$workdir/12Comparison/02Promoter/01Public_Promoter/01${species}_public_RE.bed
ls $peakf
ls $prom_region

Rscript ~/FAANG/scripts/12Comparison/02Promoter/02-2Prom_peak_RE_intersect.r $peakf $prom_region $species $tissue
sort -k1,1 -k2,2n $species/01${species}_promoter_$tissue.bed |awk '/^chr/ {print $1"\t"$2"\t"$3"\t"$4"\t"$5}' |uniq > $species/02${species}_${tissue}.promoter

done
done






###Sheep
workdir=~/FAANG/analysis
mkdir -p $workdir/12Comparison_sheep/02Promoter/02Tissue_Promoter
cd $workdir/12Comparison_sheep/02Promoter/02Tissue_Promoter

mkdir sheep

for tissue in Cerebellum CerebralCortex Lung MuscleSM
do
echo $tissue
awk '/^CM/ {print $1"\t"$9"\t"$10"\t"$6}' $workdir/08Promoter_peak/10Promoter_Tissue_Final/$tissue/${tissue}_01All_promoter_peaks.txt |sort -k1,1 -k2,2n |uniq > sheep/02sheep_${tissue}.promoter
done













####intervenn
export PATH=/mnt/ceph/bmurdoch/WGS_Dwarf_lambs/Shang/tenv/bin:$PATH
export PATH=/mnt/ceph/bmurdoch/WGS_Dwarf_lambs/Shang/Dwarf20211129/data/miniconda3/bin:$PATH

source activate /mnt/ceph/bmurdoch/WGS_Dwarf_lambs/Shang/tenv

for species in sheep cattle pig human mouse
do
cd $species
intervene upset -i 02*.promoter --filenames
intervene venn -i 02*.promoter --filenames
cd ..
done









