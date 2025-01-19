#Comparison
workdir=~/FAANG/analysis
mkdir -p $workdir/12Comparison/02Promoter/01Public_Promoter
cd $workdir/12Comparison/02Promoter/01Public_Promoter


##cattle
awk 'NR>1 {print $1"\t"$2"\t"$3}' analysis/021Promoter/01PromoterRegion/promoter2k1k_cattle_chr_sort.bed |sort -k1,1 -k2,2n > 01cattle_publih_RE.bed

##pig
awk 'NR>1 {print $1"\t"$2"\t"$3}' analysis/021Promoter/01PromoterRegion/promoter2k1k_pig_chr_sort.bed |sort -k1,1 -k2,2n > 01pig_public_RE.bed

##mouse
awk 'NR>1 {print $1"\t"$2"\t"$3}' analysis/021Promoter/01PromoterRegion/promoter2k1k_mouse_chr_sort.bed |sort -k1,1 -k2,2n > 01mouse_public_RE.bed

##human
awk 'NR>1 {print $1"\t"$2"\t"$3}' analysis/021Promoter/01PromoterRegion/promoter2k1k_human_chr_sort.bed |sort -k1,1 -k2,2n > 01human_public_RE.bed


