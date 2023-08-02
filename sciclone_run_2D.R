library("IRanges")
library("ggplot2")
library(devtools)
library(sciClone)

args <- commandArgs(trailingOnly = TRUE)
INPUT_First = args[1]
INPUT_Second = args[2]
OUTPUT_DIR = args[3]

v1 = read.table(INPUT_First, header = T)
v2 = read.table(INPUT_Second, header = T)

names = c(INPUT_First, INPUT_Second)

sc = sciClone(vafs=list(v1,v2), minimumDepth = 0,  sampleNames=names[1:2])

#create output
writeClusterTable(sc, paste0(OUTPUT_DIR , "/results.tsv"))
sc.plot1d(sc, paste0(OUTPUT_DIR, "/clusters.1d.pdf"))
sc.plot2d(sc, paste0(OUTPUT_DIR, "/clusters.2d.pdf"))
