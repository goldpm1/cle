library(QuantumClone)
library(ggplot2)

args <- commandArgs(trailingOnly = TRUE)
INPUT_First = args[1]
OUTPUT_DIR = args[2]

# "/data/project/Alzheimer/EM_cluster/QuantumClone/input/M1-3_M1-8_input_M1-3.txt"

data_df_1 <- read.table(INPUT_First, header = TRUE)
data_df_1$Genotype <- as.factor(data_df_1$Genotype)

list_df_data <- mget(ls(pattern = "data_df"))

Clustering_output <- QuantumClone(list_df_data, FREEC_list = NULL, contamination = c(0), nclone_range = 2:10,
                                  clone_priors = NULL, prior_weight = NULL, Initializations = 10, preclustering = "FLASH",
                                  simulated = FALSE, epsilon = NULL, save_plot = TRUE, ncores = 3,
                                  output_directory = NULL, model.selection = "BIC", optim = "default", keep.all.models = FALSE,
                                  force.single.copy = FALSE)
# Clustering_output
# length(Clustering_output$cluster)

first <- Clustering_output[["filtered.data"]][[1]]
first$cluster = Clustering_output$cluster

write.table(first, paste0(OUTPUT_DIR, "/output0.txt"), append = FALSE, sep = "\t", dec = ".", row.names = FALSE, col.names = TRUE, quote = FALSE)
