#!/bin/bash

workdir=~FAANG/analysis/11_2Methy
mkdir -p $workdir/11_2Methy/01Methy_region
cd $workdir/11_2Methy/01Methy_region

##00 methylation data
cd $workdir/11_2Methy/01Methy_region/00Tissues

Kimdir=/mnt/ceph/bmurdoch/Kim/Benz2616_Methylation

#cp $Kimdir/CER_WGBS*_CpGs_10x.bedgraph ./Cerebellum.methy.mdr.WGBS
#cp $Kimdir/COR_WGBS_CpGs_10x.bedgraph ./CerebralCortex.methy.mdr.WGBS
#cp $Kimdir/LNG_WGBS*_CpGs_10x.bedgraph ./Lung.methy.mdr.WGBS
#cp $Kimdir/RUM_A*WGBS*_CpGs_10x.bedgraph ./RumenAtrium.methy.mdr.WGBS
#cp $Kimdir/OVY*WGBS*_CpGs_10x.bedgraph ./Ovary.methy.mdr.WGBS
#cp $Kimdir/ADN_C**RRBS*_CpGs_10x.bedgraph ./AdrenalCortex.methy.mdr
#cp $Kimdir/ADN_M**RRBS*_CpGs_10x.bedgraph ./AdrenalMedulla.methy.mdr
#cp $Kimdir/DCOL**RRBS*_CpGs_10x.bedgraph ./DescendingColon.methy.mdr
#cp $Kimdir/DUO**RRBS*_CpGs_10x.bedgraph ./Duodenum.methy.mdr
#cp $Kimdir/GAL**RRBS*_CpGs_10x.bedgraph ./Gallbladder.methy.mdr
#cp $Kimdir/HRT_RA**RRBS*_CpGs_10x.bedgraph ./HeartRightAtrium.methy.mdr
#cp $Kimdir/HRT_RV**RRBS*_CpGs_10x.bedgraph ./HeartRightVentricle.methy.mdr
#cp $Kimdir/ILE_RRBS*_CpGs_10x.bedgraph ./Ileum.methy.mdr
#cp $Kimdir/LN_MES**RRBS*_CpGs_10x.bedgraph ./LymphNodeMesenteric.methy.mdr
#cp $Kimdir/ILE_PP**RRBS*_CpGs_10x.bedgraph ./IleumPeyersPatch.methy.mdr
#cp $Kimdir/RET**RRBS*_CpGs_10x.bedgraph ./Reticulum.methy.mdr
#cp $Kimdir/SC_RRBS*_CpGs_10x.bedgraph ./SpinalCord.methy.mdr
#cp $Kimdir/TNG**RRBS*_CpGs_10x.bedgraph ./Tongue.methy.mdr
#cp $Kimdir/TON**RRBS*_CpGs_10x.bedgraph ./Tonsil.methy.mdr
#cp $Kimdir/BLA**RRBS*_CpGs_10x.bedgraph ./Bladder.methy.mdr
#cp $Kimdir/SCOL**RRBS*_CpGs_10x.bedgraph ./SpiralColon.methy.mdr
#cp $Kimdir/MUS_SM**RRBS*_CpGs_10x.bedgraph ./MuscleSM.methy.mdr

##01 methylation position
workdir=~/20240401FAANG/analysis/11_2Methy
mkdir -p $workdir/11_2Methy/01Methy_region/01Methy_loc
cd $workdir/11_2Methy/01Methy_region/01Methy_loc

for i in $(ls $workdir/11_2Methy/01Methy_region/00Tissues/*methy.mdr.WGBS)
do
echo $i
samplename=$(basename $i .methy.mdr.WGBS)
echo $samplename
awk '$4>0 {print "CM0"$1"\t"$2"\t"$3"\t"$4}' $i >$samplename.methy.bed.WGBS

done


