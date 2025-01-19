####bash to workdir 
##cd ~/FAANG/analysis/05ChromHMM/04Chrom_fig

library(ggplot2)
library(reshape2)

#input data frame
data_raw <- read.table("~/FAANG/analysis/05ChromHMM/03Model/01Output/Cerebellum/10/emissions_10.txt", header = TRUE)
data=data_raw[c(5,8,9,7,6,10,4,2,3,1),c(1,5,4,6,3,2)]
data$state=c("TssA","TssAFlnk","TxFlnk","EnhA","EnhWk","EnhPois","EnhWkATAC","ReprPC","ReprPCWk","Quies")
data$state <- factor(data$state, levels = rev(unique(data$state)))

melted_data <- melt(data, id.vars = "state")
data$state <- factor(data$state, levels = unique(data$state))
# heapmap
heatmap_plot <- ggplot(melted_data, aes(x = variable, y = state, fill = value)) +
  geom_tile() +  
  scale_fill_gradient(low = "white", high = "#00008D", breaks = c(0, 1), labels = c("0", "1"), limits = c(0, 1)) + 
  coord_fixed() +  
  theme_minimal() + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, size = 16, hjust = 1),  
        axis.text.y = element_text(size = 16),  
        axis.title.x = element_blank(),  
        axis.title.y = element_blank(),  
        legend.position = "right",  
        legend.title = element_text(size = 10),  
        legend.text = element_text(size = 8)) + 
  labs(title = "Chromatin States Heatmap", fill = "Intensity") +  
  theme(plot.title = element_text(hjust = 0.5))
print(heatmap_plot)
ggsave("chromatin_states_heatmap.png", plot = heatmap_plot, width = 10, height = 10, dpi = 300)


####2.bash to workdir and calculated genome coverage in all tissues 
library(dplyr)

get_file_paths <- function(base_path, pattern) {
  file_paths <- list.files(base_path, pattern = pattern, full.names = TRUE, recursive = TRUE)
  return(file_paths)
}

base_path <- "~/FAANG/analysis/05ChromHMM/03Model/01Output"
file_suffix <- "_10_overlap.txt"

file_paths <- get_file_paths(base_path, file_suffix)

read_and_sort_data <- function(file_path) {
  sample_data <- read.table(file_path, header = F,skip = 1)
  sorted_values <- sort(sample_data[, 2], decreasing = TRUE) 
  return(sorted_values)
}

all_sorted_data <- lapply(file_paths, read_and_sort_data)

combined_sorted_data <- do.call(rbind, all_sorted_data)

data_values <- combined_sorted_data

mean_values <- colMeans(data_values)
sd_values <- apply(data_values, 2, sd)

result <- data.frame(State = 1:11, Mean = mean_values, SD = sd_values)

print(result)

write.table(result, file = "chromatin_states_statistics.txt", row.names = FALSE, sep = "\t")

