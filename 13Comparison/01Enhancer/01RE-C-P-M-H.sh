#Comparison
workdir=~/FAANG/analysis
mkdir -p $workdir/12Comparison/01Enhancer/01Public_RE
cd $workdir/12Comparison/01Enhancer/01Public_RE

##cattle
awk 'NR>1 {print $1"\t"$2"\t"$3}' 01Public_data/cattle-bostau9-enhancer/Filtered_Cattle_Dataset.4bed.UCSC |sort -k1,1 -k2,2n > 01cattle_publih_RE.bed
##pig
awk 'NR>1 {print $1"\t"$2"\t"$3}' 01Public_data/pig-sus11-enhancer/Filtered_Pig_Dataset.4bed.UCSC |sort -k1,1 -k2,2n > 01pig_public_RE.bed
##mouse
awk 'NR>1 {print $1"\t"$2"\t"$3}' 01Public_data/mouse-mm10-enhancer/Allfiles_SortedMerged_Mm10.bed.UCSC |sort -k1,1 -k2,2n > 01mouse_public_RE.bed
##human
awk 'NR>1 {print $1"\t"$2"\t"$3}' 01Public_data/human-hg19-enhancer/human_hg38_enhancer.bed |sort -k1,1 -k2,2n > 01human_public_RE.bed
