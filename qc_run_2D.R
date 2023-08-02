library(QuantumClone)
library(ggplot2)

args <- commandArgs(trailingOnly = TRUE)
INPUT_First = args[1]
INPUT_Second = args[2]
OUTPUT_DIR = args[3]

# "/data/project/Alzheimer/EM_cluster/QuantumClone/input/M1-3_M1-8_input_M1-3.txt"
# "/data/project/Alzheimer/EM_cluster/QuantumClone/input/M1-3_M1-8_input_M1-8.txt"

data_df_1 <- read.table(INPUT_First, header = TRUE)
data_df_2 <- read.table(INPUT_Second, header = TRUE)

data_df_1$Genotype <- as.factor(data_df_1$Genotype)
data_df_2$Genotype <- as.factor(data_df_2$Genotype)

list_df_data <- mget(ls(pattern = "data_df"))

Clustering_output <- QuantumClone(list_df_data, FREEC_list = NULL, contamination = c(0,0), nclone_range = 2:10,
                                  clone_priors = NULL, prior_weight = NULL, Initializations = 10, preclustering = "FLASH",
                                  simulated = FALSE, epsilon = NULL, save_plot = FALSE, ncores = 3,
                                  output_directory = NULL, model.selection = "BIC", optim = "default", keep.all.models = FALSE,
                                  force.single.copy = FALSE)
# Clustering_output
# length(Clustering_output$cluster)

first <- Clustering_output[["filtered.data"]][[1]]
second <- Clustering_output[["filtered.data"]][[2]]

first$cluster = Clustering_output$cluster
second$cluster = Clustering_output$cluster

write.table(first, paste0(OUTPUT_DIR, "/output0.txt"), append = FALSE, sep = "\t", dec = ".", row.names = FALSE, col.names = TRUE, quote = FALSE)
write.table(second, paste0(OUTPUT_DIR, "/output1.txt"), append = FALSE, sep = "\t", dec = ".", row.names = FALSE, col.names = TRUE, quote = FALSE)