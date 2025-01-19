# Rscript parameters
args <- commandArgs(trailingOnly = TRUE)

# load R packages
#
#library(rtracklayer)
library(data.table)
library(GenomicRanges)
#library(ChIPseeker)
#library(dplyr)

# read peak files
input_files<- list.files(path = "~/FAANG/analysis/08Promoter_peak/08Prom_border_TPM/02Prom_tis", pattern = "*_01final_promoter.txt", full.names = TRUE, recursive = F)
###Tissue name
tissue_names=unlist(lapply(strsplit(basename(input_files),"_"), function(x) x[1]))

read_peaks <- function(file) {
  prom_one <- fread(file, header = T)
  return(prom_one)
}

Prom_list <- lapply(input_files, read_peaks)    
All_prom=do.call(rbind, Prom_list)
gene_list=unique(All_prom$geneID)

library(parallel)

process_gene <- function(gene_id, All_pro) {
Gene_prom=All_prom[All_prom$geneID==gene_id,]

min_index=which.min(Gene_prom$distance)
# Step: select the closest region
prom_n=length(unique(Gene_prom$distance))
selected_promoter <- data.table(Gene_prom[min_index, ],prom_n)
#selected_promoter <- data.table(Gene_prom[which.min(abs(Gene_prom$start - selected_start) + abs(Gene_prom$end - selected_end)), ],prom_n)
return(selected_promoter)
}

result_list_gene <- mclapply(gene_list, function(gene_id) {process_gene(gene_id, All_prom)}, mc.cores = detectCores())  
Gene_promoter=do.call(rbind, result_list_gene)

# obatin All_prom  all promoters for one gene
new_promoter_gene=function(other_gene,All_prom){      
  gene_all_promoters <- All_prom[geneID == other_gene$geneID]
        
        know=other_gene$distance      
        filtered_promoters <- gene_all_promoters[gene_all_promoters$distance != know, ]
        # Step 2: calculate the frequency
        distance_table <- table(filtered_promoters$distance)
        mode_distance <- as.character(names(distance_table)[which.max(distance_table)])
        # Step 3: check the number
	if (length(mode_distance) > 0 && distance_table[mode_distance] > 1) {
	  # if existed，select the multiple 
          selected_promoters <- filtered_promoters[filtered_promoters$distance == mode_distance, ][1]
	} else {
	  # if no，average
	  mean_distance <- mean(filtered_promoters$distance)
	  closest_to_mean <- which.min(abs(filtered_promoters$distance - mean_distance))
  	  selected_promoters <- filtered_promoters[closest_to_mean, ]
	}

        prom_n=length(unique(gene_all_promoters$distance))
        new_promoter <- data.table(selected_promoters,prom_n)
return(new_promoter)
}

#### one-one gene-promoter
library(data.table)

# Step 1: Gene_promoter to data.table
Gene_pro <- as.data.table(Gene_promoter)

unique_combinations <- unique(Gene_pro[, .(start,end,width)])

# Step 2: Initializes results storage
results_list <- list()
count <- 1

# Step 3: run each one
for (i in 1:nrow(unique_combinations)) { 
  #
  current_combination <- unique_combinations[i]
  # Use all matching rows as a conditional subset
  matching_rows <- Gene_pro[start == current_combination$start & 
                            end == current_combination$end & 
                            width == current_combination$width]
  
  if (nrow(matching_rows) == 1) {
    # If the combination has only one line, add it directly to the result list
    results_list[[count]] <- matching_rows
    count <- count + 1
  } else {
    if(min(matching_rows$prom_n)>1)
    {    
    # Step 4.1: obtain the closest TSS distance promoter
    closest_promoter <- matching_rows[which.min(matching_rows$distance)]
    results_list[[count]] <- closest_promoter
    count <- count + 1
    # Step 4.2: rearrange promoter
    other_rows <- matching_rows[-which.min(matching_rows$distance)]    
    for (j in seq_len(nrow(other_rows))) {
      other_gene <- other_rows[j]    
      new_promoter=new_promoter_gene(other_gene,All_prom)  
      results_list[[count]] <- new_promoter
      count <- count + 1  
      }
    }else {
        if(length(matching_rows$prom_n[matching_rows$prom_n==1])==1)
        {
          results_list[[count]] <- matching_rows[matching_rows$prom_n==1,]
          count <- count + 1 
          other_rows <- matching_rows[-which(matching_rows$prom_n==1)]    
          for (j in seq_len(nrow(other_rows))) {
           other_gene <- other_rows[j]    
           new_promoter=new_promoter_gene(other_gene,All_prom)  
           results_list[[count]] <- new_promoter
           count <- count + 1  
           }
        }
        if(length(matching_rows$prom_n[matching_rows$prom_n==1])>1)
        {
          tmp_uniq <- matching_rows[matching_rows$prom_n==1,]          
          g_tmp=tmp_uniq[grepl("^gene",tmp_uniq$geneID)]
          if(nrow(g_tmp)!=0){
             closest_promoter <- g_tmp[which.min(g_tmp$distance)]
             results_list[[count]] <- closest_promoter
             count <- count + 1
             } else {
             closest_promoter <- tmp_uniq[which.min(tmp_uniq$distance)]
             results_list[[count]] <- closest_promoter
             count <- count + 1             
             }
          tmp_2 <- matching_rows[matching_rows$prom_n>1,]
          if(nrow(tmp_2)!=0) { 
            other_rows <- tmp_2    
            for (j in seq_len(nrow(other_rows))) {
              other_gene <- other_rows[j]    
              new_promoter=new_promoter_gene(other_gene,All_prom)  
              results_list[[count]] <- new_promoter
              count <- count + 1  }
           }
        }     
  }
 }
}

# Step 5: merge all results
final_df <- rbindlist(results_list)

Gene_pro=as.data.table(final_df)

#######Round2
unique_combinations <- unique(Gene_pro[, .(start,end,width)])

# Step 2: initialize
results_list1 <- list()
count <- 1

# Step 3:  each unique combination and process it
for (i in 1:nrow(unique_combinations)) { 
  current_combination <- unique_combinations[i]

  matching_rows <- Gene_pro[start == current_combination$start & 
                            end == current_combination$end & 
                            width == current_combination$width]
  
  if (nrow(matching_rows) == 1) {
    results_list1[[count]] <- matching_rows
    count <- count + 1
  } else {
          tmp_uniq <- matching_rows
          g_tmp=tmp_uniq[!grepl("^gene-LOC",tmp_uniq$geneID)&!grepl("^STRG",tmp_uniq$geneID)]
          if(nrow(g_tmp)!=0){
             closest_promoter <- g_tmp[which.min(g_tmp$distance)]
             results_list1[[count]] <- closest_promoter
             count <- count + 1
             } else {
             closest_promoter <- tmp_uniq[which.min(tmp_uniq$distance)]
             results_list1[[count]] <- closest_promoter
             count <- count + 1             
             }
  }
}


final_total_promoter <- rbindlist(results_list1)
write.table(as.data.frame(final_total_promoter), file = paste0("01final_promoter_gene_TSS.txt"), sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

final_anno_promoter=final_total_promoter[grepl("gene",final_total_promoter$geneID),]
write.table(as.data.frame(final_anno_promoter), file = paste0("02final_promoter_gene_TSS_annotated.txt"), sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

final_new_promoter=final_total_promoter[grepl("STRG",final_total_promoter$geneID),]
write.table(as.data.frame(final_new_promoter), file = paste0("03final_promoter_gene_TSS_new.txt"), sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

