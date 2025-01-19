# Comprehensive identification of cis-regulatory enhancers and promoters in sheep

Annotation of regulatory elements is essential for understanding mechanisms underlying gene regulation, particularly tissue-specific regulation in human and animals. We systematically characterized 274,682 enhancers and 25,975 promoters across 24 tissues in sheep using six high-resolution assays: ChIP-seq, ATAC-seq, CAGE-seq, RRBS, WGBS, and RNA-seq. This study provides a robust framework for exploring cis-regulatory mechanisms and tissue-specific regulation, advancing the functional annotation of the sheep reference genome. Here is the **code** for this study.


## Contents:
- 01 ChIP-seq
- 02 CAGE-seq
- 03 ATAC-seq
- 04 CTCF
- 05 ChromHMM
- 06 PeakAnnoation
- 07 Enahncer_peak
- 08 RNA-seq_TPM
- 09 Promoter_peak
- 10 Correlation
- 11 Tissue_specific
- 12 Methylation_geneExp
- 13 Comparison
- 14 SNP


## 01 ChIP-seq

Nextflow was used to identified peaks from ChIP-seq. 
> scripts for narrow marks:

```sh
workdir=~/20240401FAANG/analysis/01ChIP-seq/02nfcore_5e3/01Narrow1
genomedir=/mnt/ceph/bmurdoch/ARS-UI_Ramb_v2.0
for Tissue in '24tissues'
do
echo $Tissue
#narrow
output=$workdir/NF_$Tissue
mkdir -p $output 
cd $output
nextflow run nf-core/chipseq --input ~/20240401FAANG/analysis/00Sample_info/sample_information/${Tissue}_sample_narrow.csv \
--outdir $output -profile singularity --fasta $genomedir/GCA_016772045.1_ARS-UI_Ramb_v2.0_genomic.fa \
--gff $genomedir/GCF_016772045.1_ARS-UI_Ramb_v2.0_genomic_chrom.gff --macs_gsize 2600000000 --narrow_peak \
--aligner bowtie2 --bowtie2_index $genomedir --macs_pvalue 0.005 --skip_fastqc --skip_preseq --skip_plot_profile --skip_igv --skip_qc --skip_trimming --skip_picard_metrics --skip_multiqc --skip_deseq2_qc --skip_plot_finerprint 
cd ..
done
```
> scripts for broad marks:

```sh
workdir=~/20240401FAANG/analysis/01ChIP-seq/02nfcore_5e3/02Broad
genomedir=/mnt/ceph/bmurdoch/ARS-UI_Ramb_v2.0

mkdir -p $workdir/00Tools 
export NXF_SINGULARITY_CACHEDIR=$workdir/00Tools

for Tissue in '24tissues'
do
echo $Tissue
#broad
Boutput=$workdir/NF_$Tissue
mkdir -p $Boutput 
cd $Boutput
nextflow run nf-core/chipseq --input ~/20240401FAANG/analysis/00Sample_info/sample_information/${Tissue}_sample_broad.csv \
--outdir $Boutput -profile singularity --fasta $genomedir/GCA_016772045.1_ARS-UI_Ramb_v2.0_genomic.fa \
--gff $genomedir/GCF_016772045.1_ARS-UI_Ramb_v2.0_genomic_chrom.gff --macs_gsize 2600000000 --broad_cutoff 0.005 \
--aligner bowtie2 --bowtie2_index $genomedir --macs_pvalue 0.005 --skip_fastqc --skip_preseq --skip_plot_profile --skip_igv --skip_qc --skip_trimming --skip_picard_metrics --skip_multiqc --skip_deseq2_qc --skip_plot_fingerprint
rm -rf work genome
cd ..
done
```

## 02 CAGE-seq
Identify TSSs, Enhancers and Enhancer-TSS links from .

```r
library("CAGEfightR")
library("GenomicFeatures")
bw_plus_files <- list.files(path = "~/20240401FAANG/data/CAGE_data/bigWig", pattern = "\\.plus.bw$", full.names = TRUE)
bw_minus_files <- list.files(path = "~/20240401FAANG/data/CAGE_data/bigWig", pattern = "\\.minus.bw$", full.names = TRUE)
bw_plus_path=bw_plus_files
bw_minus_path=bw_minus_files
bw_plus <- BigWigFileList(bw_plus_path)
bw_minus <- BigWigFileList(bw_minus_path)
sampleDesign <- DataFrame(Name = samplename ,BigWigPlus = c(bw_plus),BigWigMinus = c(bw_minus))
names(bw_plus) <- sampleDesign$Name
names(bw_minus) <- sampleDesign$Name
CTSSs <- quantifyCTSSs(plusStrand=bw_plus,
                       minusStrand=bw_minus,
                       design=sampleDesign,
                       genome=genome_info)
supportedCTSSs <- subsetBySupport(CTSSs,inputAssay ="counts", outputColumn = "support", unexpressed=1, minSamples =0)
TSSs <- quickTSSs(supportedCTSSs)
enhancers <- quickEnhancers(supportedCTSSs)
TSSs <- assignTxType(TSSs, txModels=txdb)
enhancers <- assignTxType(enhancers, txModels=txdb)
enhancers <- subset(enhancers, txType %in% c("intergenic", "intron"))
rowRanges(TSSs)$clusterType <- "TSS"
rowRanges(enhancers)$clusterType <- "enhancer"
RSE <- combineClusters(TSSs, enhancers, removeIfOverlapping="object1")             
RSE <- assignGeneID(RSE, geneModels=txdb, outputColumn='geneID')
rowRanges(RSE)$clusterType <- factor(rowRanges(RSE)$clusterType,
                                     levels=c("TSS", "enhancer"))
for (i in 1:24)
{
print(paste0("The ",i," tissue: ",samplename[i]))
individual_RSE <- subset(RSE, RSEcount[,i] > 0)
ind_TSS <- subset(rowRanges(individual_RSE), clusterType == "TSS")
ind_Enh <- subset(rowRanges(individual_RSE), clusterType == "enhancer")
ind_links <- findLinks(individual_RSE, 
                   inputAssay="counts", 
                   maxDist=10000,                
                   method="kendall",   
                   directional="clusterType")
write.csv(as.data.frame(ind_TSS), paste0("56tissues/",samplename[i],"_TSSs_results.csv"), row.names = TRUE,quote=F)
write.csv(as.data.frame(ind_Enh), paste0("56tissues/",samplename[i],"_enhancers_results.csv"), row.names = TRUE,quote=F)
write.csv(as.data.frame(ind_links), paste0("56tissues/",samplename[i],"_links_results.csv"), row.names = TRUE,quote=F)
}
```

## 03 ATAC-seq
```sh
nextflow run nf-core/atacseq --input ~/20240401FAANG/analysis/03ATAC/00Sample/02tissue_format/${Tissue}.ATAC.sample.csv \
--outdir $output -profile singularity --fasta $genomedir/GCA_016772045.1_ARS-UI_Ramb_v2.0_genomic.fa \
--gff $genomedir/GCF_016772045.1_ARS-UI_Ramb_v2.0_genomic_chrom.gff --macs_gsize 2600000000 --narrow_peak \
--aligner bowtie2 --bowtie2_index $genomedir --macs_pvalue 0.005 --skip_fastqc --skip_preseq --skip_igv --skip_qc --skip_picard_metrics --skip_plot_fingerprint
```

## 04 CTCF
Identify peaks from CTCF.
```sh
workdir=~/20240401FAANG/analysis/04CTCF/01Peak/02${Tissue} 
genomedir=/mnt/ceph/bmurdoch/ARS-UI_Ramb_v2.0
output=$workdir
mkdir -p $output 
cd $output
nextflow run nf-core/chipseq --input ~/20240401FAANG/analysis/04CTCF/01Peak/01${Tissue}_CTCF_narrow.csv \
--outdir $output -profile singularity --fasta $genomedir/GCA_016772045.1_ARS-UI_Ramb_v2.0_genomic.fa \
--gff $genomedir/GCF_016772045.1_ARS-UI_Ramb_v2.0_genomic_chrom.gff --macs_gsize 2600000000 --narrow_peak \
--aligner bowtie2 --bowtie2_index $genomedir --macs_pvalue 0.01 --skip_fastqc --skip_preseq --skip_plot_profile --skip_igv --skip_qc --skip_trimming --skip_picard_metrics --skip_multiqc --skip_deseq2_qc --skip_plot_fingerprint -resume
```

## 05 ChromHMM

 <https://github.com/shang-qian/FAANG_sheep_code_NCsub/tree/main/05ChromHMM>
 
 

## 06 PeakAnnotation

 <https://github.com/shang-qian/FAANG_sheep_code_NCsub/tree/main/06PeakAnnotation>
 
 
 
## 07 Enhancer_peak

 <https://github.com/shang-qian/FAANG_sheep_code_NCsub/tree/main/07Enhancer_peak>
 
 
 
## 08 RNA-seq_TPM

 <https://github.com/shang-qian/FAANG_sheep_code_NCsub/tree/main/08RNA-seq_TPM)>
 
 
 
## 09 Promoter_peak

 <https://github.com/shang-qian/FAANG_sheep_code_NCsub/tree/main/09Promoter_peak>
 
 
 
## 10Correlation

 <https://github.com/shang-qian/FAANG_sheep_code_NCsub/tree/main/10Interaction>
 
 
 
## 11 Tissue_specific

 <https://github.com/shang-qian/FAANG_sheep_code_NCsub/tree/main/11Tissue_specific>
 
 
 
## 12 Methylation

 <https://github.com/shang-qian/FAANG_sheep_code_NCsub/tree/main/12Methylation_geneExp>
 
 

 
## 13 Comparison

 <https://github.com/shang-qian/FAANG_sheep_code_NCsub/tree/main/13Comparison>
 
 
## 14 SNP

 <https://github.com/shang-qian/FAANG_sheep_code_NCsub/tree/main/14SNP>
 
  
 
