#!/bin/bash
#SBATCH --job-name=nextflow
#SBATCH --partition=reg
#SBATCH --cpus-per-task=40
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH -t 168:00:00 # run time (hh:mm:ss)
#SBATCH -o sbatch_results-%j # output and error file name (%j expands to jobID)

export PATH=/mnt/ceph/sxie/FAANG/data/nextflow:$PATH
source activate nextflow

workdir=~/FAANG/analysis/01ChIP-seq/02nfcore_5e3/01Narrow
genomedir=/mnt/ceph/bmurdoch/ARS-UI_Ramb_v2.0

for Tissue in Adrenalcortex Adrenalmedulla Bladder Cerebellum Cerebralcortex Descendingcolon Duodenum Gallbladder Heartrightatrium Heartrightventricle 
Ileum Ileumpeyerspatch Lung Lymphnodemesenteric Skeletalmuscle Ovary Oviduct Reticulum Rumenatrium Spinalcord Spiralcolon Tongue Tonsil Uterus 
do
echo $Tissue
#narrow
output=$workdir/NF_$Tissue
mkdir -p $output 
cd $output
nextflow run nf-core/chipseq --input ~/FAANG/analysis/00Sample_info/sample_information/${Tissue}_sample_narrow.csv \
--outdir $output -profile singularity --fasta $genomedir/GCA_016772045.1_ARS-UI_Ramb_v2.0_genomic.fa \
--gff $genomedir/GCF_016772045.1_ARS-UI_Ramb_v2.0_genomic_chrom.gff --macs_gsize 2600000000 --narrow_peak \
--aligner bowtie2 --bowtie2_index $genomedir --macs_pvalue 0.005 --skip_fastqc --skip_preseq --skip_plot_profile --skip_igv --skip_qc --skip_trimming --skip_picard_metrics --skip_multiqc --skip_deseq2_qc --skip_plot_finerprint 
cd ..
done
