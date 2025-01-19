#!/bin/bash
#SBATCH --job-name=ATAC-seq
#SBATCH --partition=reg
#SBATCH --cpus-per-task=40
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH -t 168:00:00 # run time (hh:mm:ss)
#SBATCH -o sbatch_results-%j # output and error file name (%j expands to jobID)


#Tissue_name_ATAC
cd ~/FAANG/analysis/03ATAC/00Sample
for i in $(ls 02tissue_format/*.ATAC.sample.csv)
do
samplen=$(basename $i .ATAC.sample.csv)
echo $samplen >>ATAC_tissue.sample
done

workdir=~/20240401FAANG/analysis/03ATAC/02nfcore
genomedir=/mnt/ceph/bmurdoch/ARS-UI_Ramb_v2.0

for Tissue in "ATAC_tissue.sample"
do
echo $Tissue
output=$workdir/${Tissue}_ATAC
mkdir -p $output 
cd $output
nextflow run nf-core/atacseq --input ~/20240401FAANG/analysis/03ATAC/00Sample/02tissue_format/${Tissue}.ATAC.sample.csv \
--outdir $output -profile singularity --fasta $genomedir/GCA_016772045.1_ARS-UI_Ramb_v2.0_genomic.fa \
--gff $genomedir/GCF_016772045.1_ARS-UI_Ramb_v2.0_genomic_chrom.gff --macs_gsize 2600000000 --narrow_peak \
--aligner bowtie2 --bowtie2_index $genomedir --macs_pvalue 0.005 --skip_fastqc --skip_preseq --skip_igv --skip_qc --skip_picard_metrics --skip_plot_fingerprint

cd ..
done


for i in $(ls ~/FAANG/analysis/03ATAC/02nfcore/*_ATAC*/bowtie2/merged_library/macs2/narrow_peak/*_REP1.mLb.clN_peaks.narrowPeak)
do
echo $i
sample=$(basename $i _REP1.mLb.clN_peaks.narrowPeak)
echo $sample
awk '/^CM/&&$9>=3 {print $0}' $i > ${sample}.ATAC.peak
done
