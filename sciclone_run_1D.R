library("IRanges")
library("ggplot2")
library(devtools)
library(sciClone)

args <- commandArgs(trailingOnly = TRUE)
INPUT_First = args[1]
OUTPUT_DIR = args[2]

v1 = read.table(INPUT_First, header = T)

names = c(INPUT_First)

sc = sciClone(vafs=list(v1), minimumDepth = 0,  sampleNames=names[1])

#create output
writeClusterTable(sc, paste0(OUTPUT_DIR , "/results.tsv"))
sc.plot1d(sc, paste0(OUTPUT_DIR, "/clusters.1d.pdf"))
