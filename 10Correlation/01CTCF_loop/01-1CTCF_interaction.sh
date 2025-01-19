#!/bin/bash

#CTCF
workdir=~/FAANG/analysis
mkdir -p $workdir/10interaction/01CTCF/01CTCF_motif
cd $workdir/10interaction/01CTCF/01CTCF_motif

module load bedtools samtools
Peak_BED=~/FAANG/analysis/10interaction/01CTCF/00Peak/01Peaks_CTCF.bed
awk 'NR>1 {print $0}' $Peak_BED > 00Input.bed
# get bed sequences
GENOME_FA="/mnt/ceph/bmurdoch/ARS-UI_Ramb_v2.0/GCA_016772045.1_ARS-UI_Ramb_v2.0_genomic.fa"
bedtools getfasta -fi $GENOME_FA -bed 00Input.bed -fo 01peak_sequences.fa

cp ~/FAANG/analysis/JASPAR2024_CORE_vertebrates_non-redundant_pfms.meme 02ctcf.meme

workdir=~/FAANG/analysis
mkdir -p $workdir/10interaction/01CTCF/01CTCF_motif/03fimo_output1
cd $workdir/10interaction/01CTCF/01CTCF_motif

#generate fimo.tsv
fimo --oc 03fimo_output --thresh 1e-4 --max-stored-scores 1000000 02ctcf.meme 01peak_sequences.fa





