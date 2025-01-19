#!/bin/bash

workdir=~/FAANG/analysis
mkdir -p $workdir/10interaction/04Pairs
cd $workdir/10interaction/04Pairs

module load R/4.2.3

awk '{print FR}' $workdir/10interaction/03Cor8/06TPM_interaction_Effect_final.txt

awk '$10!="None"||$12!="NA" {print $0}' $workdir/10interaction/03Cor8/06TPM_interaction_Effect_final.txt >01Exp_gene_pair_final.txt
awk '{print $0}' $workdir/10interaction/03Cor8/07Unexp_interaction_final.txt >02Unexp_gene_pair_final.txt


###CTCF or CAGE #
awk '$10!="None" {print $1,$2}' 01Exp_gene_pair_final.txt |uniq |wc -l
awk '$10!="None" {print $1}' 01Exp_gene_pair_final.txt |sort |uniq |wc -l
awk '$10!="None" {print $2}' 01Exp_gene_pair_final.txt |sort |uniq |wc -l


