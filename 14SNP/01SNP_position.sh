workdir=~/FAANG/analysis
mkdir -p $workdir/13SNP/02New
cd $workdir/13SNP/02New


###1. All snp position
awk '{OFS="\t"; if (NR > 1) print "chr"$1, $3, $3+1, $4, $5}' ../01data_SNP/HD_SNP_ramb2.txt |\
awk '{if(NF==5) {print $0}
      if(NF<5) {print $1"\t"$2"\t"$3"\t"$4"\t."}}' |sort -k1,1 -k2,2n |uniq > 01SNP_position_chr_all.bed
      

###2. genomic feature: exon intron promoter enhancer 
mkdir 02Genomic_region
cd 02Genomic_region
module load R/4.2.3
R


###3. QTL
workdir=~/FAANG/analysis
cd $workdir/13SNP/02New

awk -F "\t" '/^Chr/ {split($1,a,".");print "chr"a[2]"\t"$2"\t"$3"\t"$4}' ../01data_SNP/QTLdb_sheepOAR_rambo2.bed |sort -k1,1 -k2,2n |uniq >03QTL_position_chr.bed



###4. SNP intersect promoter_QTL
workdir=~/FAANG/analysis
cd $workdir/13SNP/02New
mkdir 10Intersect
cd 10Intersect

awk 'NR>1 {print $2"\t"$3"\t"$4"\t"$7"\t"$8"\t"$9"\t"$11"\t"$13}' ../04Promoter_QTL.txt >01Promoter_QTL.bed
awk 'NR>1 {print $2"\t"$3"\t"$4"\t"$7"\t"$8"\t"$9"\t"$11}' ../05Enhancer_QTL.txt >02Enhancer_QTL.bed



