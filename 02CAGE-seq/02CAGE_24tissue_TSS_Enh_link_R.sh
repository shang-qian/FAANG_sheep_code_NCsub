#!/bin/bash

workdir=~/FAANG/analysis/02CAGE
mkdir -p $workdir
cd $workdir

###Rscripts
module load R/4.2.3
R

library("CAGEfightR")

bw_plus_files <- list.files(path = "~/FAANG/data/CAGE_data/bigWig", pattern = "\\.plus.bw$", full.names = TRUE)
bw_minus_files <- list.files(path = "~/FAANG/data/CAGE_data/bigWig", pattern = "\\.minus.bw$", full.names = TRUE)
bw_plus_path=bw_plus_files
bw_minus_path=bw_minus_files
bw_plus <- BigWigFileList(bw_plus_path)
bw_minus <- BigWigFileList(bw_minus_path)

samplename=c("AdrenalCortex","AdrenalMedulla","Bladder","Cerebellum","CerebralCortex","DescendingColon","Duodenum","Gallbladder","HeartRightAtrium",
"HeartRightVentricle","Ileum","IleumPeyersPatch","Lung","LymphNodeMesenteric","MuscleSM","Ovary","Oviduct","Reticulum","RumenAtrium","SpinalCord",
"SpiralColon","Tongue","Tonsil","Uterus")

sampleDesign <- DataFrame(Name = samplename ,BigWigPlus = c(bw_plus),BigWigMinus = c(bw_minus))
names(bw_plus) <- sampleDesign$Name
names(bw_minus) <- sampleDesign$Name

####rtracklayer extract FASTA seq info and name
#fasta_path <- "ram2_all.fa"
#fasta <- import(fasta_path, format = "fasta")
#seqnames <- names(fasta); seqlengths <- width(fasta)
#nc_values <- sub("^(NC_\\d+\\.\\d+).*", "\\1", seqnames)
#genome_info <- Seqinfo(seqnames = nc_values, seqlengths = seqlengths+1, genome = "Ramv2")
#save(genome_info, file = "genome_info.RData")
####annotation of TSSs and enhancers
#library(GenomicFeatures)
#gff_file <- "genomic.gff"
#txdb <- makeTxDbFromGFF(file = gff_file, format = "gff")
#saveDb(txdb, file = "Annotation_txdb.sqlite")

load("~/FAANG/analysis/02CAGE/genome_info.RData")
CTSSs <- quantifyCTSSs(plusStrand=bw_plus,
                       minusStrand=bw_minus,
                       design=sampleDesign,
                       genome=genome_info)

###total samples
supportedCTSSs <- subsetBySupport(CTSSs,inputAssay ="counts", outputColumn = "support", unexpressed=1, minSamples =0)
TSSs <- quickTSSs(supportedCTSSs)
enhancers <- quickEnhancers(supportedCTSSs)

library("GenomicFeatures")
txdb <- loadDb("~/FAANG/analysis/02CAGE/Annotation_txdb.sqlite")
TSSs <- assignTxType(TSSs, txModels=txdb)
enhancers <- assignTxType(enhancers, txModels=txdb)
enhancers <- subset(enhancers, txType %in% c("intergenic", "intron"))
rowRanges(TSSs)$clusterType <- "TSS"
rowRanges(enhancers)$clusterType <- "enhancer"

# Combine TSSs and enhancers, discarding TSSs if they overlap enhancers
RSE <- combineClusters(TSSs, enhancers, removeIfOverlapping="object1")             
#Add gene ID
RSE <- assignGeneID(RSE, geneModels=txdb, outputColumn='geneID')
rowRanges(RSE)$clusterType <- factor(rowRanges(RSE)$clusterType,
                                     levels=c("TSS", "enhancer"))

RSEcount=assay(RSE,"counts")
# write.csv(as.data.frame(RSEcount), "assays_counts_24tissues.csv", row.names = TRUE,quote=F)
row_ranges <- rowRanges(RSE)
# write.csv(as.data.frame(row_ranges), "row_ranges_24tissues.csv", row.names = TRUE,quote=F)

links <- findLinks(RSE, 
                   inputAssay="counts", 
                   maxDist=10000,                
                   method="kendall",   
                   directional="clusterType")
# write.csv(as.data.frame(links), "links_results_24tissues.csv", row.names = TRUE,quote=F)

###output TSSs, enhancers and Links
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
                   
print(length(ind_TSS))
print(length(ind_Enh))
print(length(ind_links))                 
                   
write.csv(as.data.frame(ind_TSS), paste0("24tissues/",samplename[i],"_TSSs_results.csv"), row.names = TRUE,quote=F)
write.csv(as.data.frame(ind_Enh), paste0("24tissues/",samplename[i],"_enhancers_results.csv"), row.names = TRUE,quote=F)
write.csv(as.data.frame(ind_links), paste0("24tissues/",samplename[i],"_links_results.csv"), row.names = TRUE,quote=F)

}

