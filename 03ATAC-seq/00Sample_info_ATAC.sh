#!/bin/bash

###ATAC-seq sample
sampledir=~/FAANG/analysis/03ATAC/00Sample
mkdir -p $sampledir/01raw_format
cd $sampledir

for i in $(ls Benz2616_ATAC-seq_analyses/00-Raw/*_R1.fastq.gz)
do
R1=$i
sampleN=$(basename $i _R1.fastq.gz)
R2=Benz2616_ATAC-seq_analyses/00-Raw/${sampleN}_R2.fastq.gz
echo $sampleN |awk -v tis=${sampleN} -v fq1=$R1 -v fq2=$R2 'BEGIN {print "sample,fastq_1,fastq_2,replicate"}
                                     {print tis","fq1","fq2","1}' >01raw_format/01${sampleN}.format.raw
done

cd 01raw_format
mv 01X3034.Adn.C.format.raw ../02tissue_format/AdrenalCortex.ATAC.sample.csv
mv 01X3046.Adn.M.format.raw ../02tissue_format/AdrenalMedulla.ATAC.sample.csv
mv 01X4039.Bla.format.raw ../02tissue_format/Bladder.ATAC.sample.csv
mv 01X1010.Cor.format.raw ../02tissue_format/CerebralCortex.ATAC.sample.csv
mv 01X1016.Cer.format.raw ../02tissue_format/Cerebellum.ATAC.sample.csv
mv 01X2088.Duo.format.raw ../02tissue_format/Duodenum.ATAC.sample.csv
mv 01X3010.Gal.format.raw ../02tissue_format/Gallbladder.ATAC.sample.csv
mv 01X5074.Hrt.RA_MACS.NPB_2.format.raw ../02tissue_format/HeartRightAtrium.ATAC.sample.csv
mv 01X5086.Hrt.RV_MACS.NPB.format.raw ../02tissue_format/HeartRightVentricle.ATAC.sample.csv
mv 01X2114.Ile.format.raw ../02tissue_format/Ileum.ATAC.sample.csv
mv 01X2125.Ile.PP.format.raw ../02tissue_format/IleumPeyersPatch.ATAC.sample.csv
mv 01X5034.Lng.format.raw ../02tissue_format/Lung.ATAC.sample.csv
mv 01LN.mes.2108.format.raw ../02tissue_format/LymphNodeMesenteric.ATAC.sample.csv
mv 01X6042.Mus.SM_MACS.NPB.format.raw ../02tissue_format/MuscleSM.ATAC.sample.csv
mv 01X4057.Ovy.format.raw ../02tissue_format/Ovary.ATAC.sample.csv
mv 01Reticulum.2050.format.raw ../02tissue_format/Reticulum.ATAC.sample.csv
mv 01X2143.Scol.format.raw ../02tissue_format/SpiralColon.ATAC.sample.csv
mv 01X6133.SC.format.raw ../02tissue_format/SpiralCord.ATAC.sample.csv
mv 01X1151.Tng.format.raw ../02tissue_format/Tongue.ATAC.sample.csv
mv 01X1161.Ton_MACS.NPB.format.raw ../02tissue_format/Tonsil.ATAC.sample.csv
