#!/bin/bash

workdir=~/FAANG/analysis
mkdir -p $workdir/10interaction/01CTCF/00Peak
cd $workdir/10interaction/01CTCF/00Peak

module load R/4.2.3

C1=~/FAANG/analysis/04CTCF/01Peak/02CER_F1/bowtie2/mergedLibrary/macs2/narrowPeak/CER_F1_CTCF_peaks.narrowPeak
C2=~/FAANG/analysis/04CTCF/01Peak/02CER_F2/bowtie2/mergedLibrary/macs2/narrowPeak/CER_F2_CTCF_peaks.narrowPeak
C3=~/FAANG/analysis/04CTCF/01Peak/02CER_M1/bowtie2/mergedLibrary/macs2/narrowPeak/CER_M1_CTCF_peaks.narrowPeak
C4=~/FAANG/analysis/04CTCF/01Peak/02CER_M2/bowtie2/mergedLibrary/macs2/narrowPeak/CER_M2_CTCF_peaks.narrowPeak
L1=~/FAANG/analysis/04CTCF/01Peak/02LIV_F2/bowtie2/mergedLibrary/macs2/narrowPeak/LIV_F2_CTCF_peaks.narrowPeak
L2=~/FAANG/analysis/04CTCF/01Peak/02LIV_M1/bowtie2/mergedLibrary/macs2/narrowPeak/LIV_M1_CTCF_peaks.narrowPeak
L3=~/FAANG/analysis/04CTCF/01Peak/02LIV_M2/bowtie2/mergedLibrary/macs2/narrowPeak/LIV_M2_CTCF_peaks.narrowPeak
L4=~/FAANG/analysis/04CTCF/01Peak/02SPL_M2/bowtie2/mergedLibrary/macs2/narrowPeak/SPL_M2_CTCF_peaks.narrowPeak




