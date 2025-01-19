#!/bin/bash

workdir=~/FAANG/analysis/02CAGE
###enhancer
mkdir -p $workdir/01Enhancer
cd $workdir/01Enhancer
for enh in $(ls ../../24tissues/*enhancers_results.csv)
do
TN=$(basename $enh _enhancers_results.csv)
echo $TN
awk -F, -v Tis=$TN '/^NC/ {print $2"\t"$3"\t"$4"\t"$5}' $enh |uniq > ${TN}_enhancer
ln -sf ${TN}.CAGE.bed ${TN}_enhancer 

echo -e $TN $(wc -l ${TN}_enhancer) >>../Enhancer.stats
done

##TSS
mkdir -p $workdir/02TSS
cd $workdir/02TSS
for enh in $(ls ../../24tissues/*TSSs_results.csv)
do
TN=$(basename $enh _TSSs_results.csv)
echo $TN
awk -F, -v Tis=$TN '/^NC/ {print $2"\t"$3"\t"$4"\t"$5}' $enh |uniq > ${TN}_TSS
ln -sf ${TN}.CAGE.TSS.bed ${TN}_TSS
echo -e $TN $(wc -l ${TN}_TSS) >>../TSS.stats
done

###Link
workdir=~/20240401FAANG/analysis/02CAGE
mkdir -p $workdir/03Link
cd $workdir/03Link
for enh in $(ls ../../24tissues/*links_results.csv)
do
TN=$(basename $enh _links_results.csv)
echo $TN
awk -F, -v Tis=$TN '$2~/^NC/ {print $2"\t"$3"\t"$4"\t"$15"\t"$16"\t"$17}' $enh |uniq > ${TN}_link
TSSn=$(awk '{print $1"\t"$2"\t"$3}' ${TN}_link |sort -k1,1 -k2,2n |uniq |wc -l)
Enhn=$(awk '{print $4"\t"$5"\t"$6}' ${TN}_link |sort -k1,1 -k2,2n |uniq |wc -l)
echo -e $TN $(wc -l ${TN}_link) $Enhn $TSSn >>../Links.stats
done

