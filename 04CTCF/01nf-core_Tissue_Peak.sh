#!/bin/bash
#SBATCH --job-name=CTCF
#SBATCH --partition=reg
#SBATCH --cpus-per-task=40
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH -t 168:00:00 # run time (hh:mm:ss)
#SBATCH -o sbatch_results-%j # output and error file name (%j expands to jobID)


for i in LIV CER SPL
do
for j in F1 F2 M1 M2 
do
echo $i
cat ~/FAANG/analysis/00Sample_info/00Sampleinfo/Head.info ${i}_${j}_CTCF.format ${i}_${j}_Input.format > 01${i}_${j}_CTCF_narrow.csv


Tissue=${i}_${j}
workdir=~/20240401FAANG/analysis/04CTCF/01Peak/02${Tissue} 
mkdir $workdir 
genomedir=/mnt/ceph/bmurdoch/ARS-UI_Ramb_v2.0

output=$workdir
mkdir -p $output 
cd $output
nextflow run nf-core/chipseq --input ~/20240401FAANG/analysis/04CTCF/01Peak/01${Tissue}_CTCF_narrow.csv \
--outdir $output -profile singularity --fasta $genomedir/GCA_016772045.1_ARS-UI_Ramb_v2.0_genomic.fa \
--gff $genomedir/GCF_016772045.1_ARS-UI_Ramb_v2.0_genomic_chrom.gff --macs_gsize 2600000000 --narrow_peak \
--aligner bowtie2 --bowtie2_index $genomedir --macs_pvalue 0.01 --skip_fastqc --skip_preseq --skip_plot_profile --skip_igv --skip_qc --skip_trimming --skip_picard_metrics --skip_multiqc --skip_deseq2_qc --skip_plot_fingerprint -resume
cd ..

done
done


