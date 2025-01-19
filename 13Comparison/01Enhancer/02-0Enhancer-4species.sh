#Comparison
workdir=~/FAANG/analysis
mkdir -p $workdir/12Comparison_sheep/01Enhancer/02Tissue_Enhancer
cd $workdir/12Comparison_sheep/01Enhancer/02Tissue_Enhancer


##1 merge peak
mkdir 01Enh_peak
cd 01Enh_peak

for species in "cattle" "pig"
do
for i in Cerebellum Cortex Lung Muscle
do
H3K4me1=$workdir/12Comparison/01Enhancer/01Public_RE/01Public_data/${species}-peak/*H3K4me1*${i}*
H3K27ac=$workdir/12Comparison/01Enhancer/01Public_RE/01Public_data/${species}-peak/*H3K27ac*${i}*
ls $H3K4me1
ls $H3K27ac
Rscript ~/FAANG/scripts/12Comparison/01Enhancer/02-1Enh_peak_merge.r $H3K4me1 $H3K27ac $species $i
done
done


for species in "human" "mouse"
do
echo $species
for i in cerebellum cortex lung muscle
do
H3K4me1=$workdir/12Comparison/01Enhancer/01Public_RE/01Public_data/${species}-peak/*${i}*H3K4me1*
H3K27ac=$workdir/12Comparison/01Enhancer/01Public_RE/01Public_data/${species}-peak/*${i}*H3K27ac*
ls $H3K4me1
ls $H3K27ac
Rscript ~/FAANG/scripts/12Comparison/01Enhancer/02-1Enh_peak_merge.r $H3K4me1 $H3K27ac $species $i
done
done





###2. generate bed file for four species and four tissues
workdir=~/FAANG/analysis
mkdir -p $workdir/12Comparison_sheep/01Enhancer/02Tissue_Enhancer/01Enh_peak
cd $workdir/12Comparison_sheep/01Enhancer/02Tissue_Enhancer/01Enh_peak

for species in "human" "mouse" "cattle" "pig"
do
echo $species
for tissue in "Cerebellum" "Cortex" "Lung" "Muscle"
do
modified_tissue=$(echo "$tissue" | sed 's/^.//')
ls $species/${species}_peaks_Enh_*${modified_tissue}.bed
awk '/^chr/ {print $0}' $species/${species}_peaks_Enh_*${modified_tissue}.bed > $species/${species}_peaks_Enh_${tissue}_chr.bed
done
done



 

#####2 extract the tissue's enhancer using merge peak and potential 
workdir=~/FAANG/analysis
mkdir -p $workdir/12Comparison_sheep/01Enhancer/02Tissue_Enhancer
cd $workdir/12Comparison_sheep/01Enhancer/02Tissue_Enhancer

mkdir 02Enh_tissue
cd 02Enh_tissue

for species in human mouse pig cattle
do
mkdir $species
for tissue in Cerebellum Cortex Lung Muscle
do
echo $tissue
peakf=$workdir/12Comparison/01Enhancer/02Tissue_Enhancer/01Enh_peak/${species}/${species}_peaks_Enh_${tissue}_chr.bed
ls $peakf
sort -k1,1 -k2,2n $peakf |awk '/^chr/ {print $1"\t"$2"\t"$3"\t"$4"\t"$5}' |uniq > $species/02${species}_${tissue}.enhancer
done
done


#####2 extract the tissue's enhancer using merge peak and potential 
workdir=~/FAANG/analysis
mkdir -p $workdir/12Comparison_sheep/01Enhancer/02Tissue_Enhancer
cd $workdir/12Comparison_sheep/01Enhancer/02Tissue_Enhancer

mkdir 02Enh_tissue
cd 02Enh_tissue

for species in human mouse pig cattle
do
mkdir $species
for tissue in Cerebellum Cortex Lung Muscle
do
echo $tissue
peakf=$workdir/12Comparison/01Enhancer/02Tissue_Enhancer/01Enh_peak/${species}/${species}_peaks_Enh_${tissue}_chr.bed
enh_region=$workdir/12Comparison/01Enhancer/01Public_RE/01${species}_public_RE.bed
ls $peakf
ls $enh_region

Rscript ~/FAANG/scripts/12Comparison/01Enhancer/02-2Enh_peak_RE_intersect.r $peakf $enh_region $species $tissue
sort -k1,1 -k2,2n $species/01${species}_enhancer_$tissue.bed |awk '/^chr/ {print $1"\t"$2"\t"$3"\t"$4"\t"$5}' |uniq > $species/02${species}_${tissue}.enhancer

done
done



###3.Sheep enhancer
mkdir sheep

for tissue in Cerebellum CerebralCortex Lung MuscleSM
do
echo $tissue
awk '/^CM/ {print $1"\t"$2"\t"$3"\t"$4}' $workdir/07Enhancer_peak/9Enhancer_Tissue_Final/$tissue/${tissue}_01All_enhancer.txt |sort -k1,1 -k2,2n |uniq > sheep/02sheep_${tissue}.enhancer
done




####intervenn
export PATH=/mnt/ceph/bmurdoch/WGS_Dwarf_lambs/Shang/tenv/bin:$PATH
export PATH=/mnt/ceph/bmurdoch/WGS_Dwarf_lambs/Shang/Dwarf20211129/data/miniconda3/bin:$PATH

source activate /mnt/ceph/bmurdoch/WGS_Dwarf_lambs/Shang/tenv

for species in sheep cattle pig human mouse
do
cd $species
intervene upset -i 02*.enhancer --filenames
intervene venn -i 02*.enhancer --filenames
cd ..
done









