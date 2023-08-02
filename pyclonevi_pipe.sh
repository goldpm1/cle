#!/bin/bash
#$ -cwd
#$ -S /bin/bash

# INPUT_TSV="/data/project/Alzheimer/EM_cluster/old/pilot/04.EM_input/pyclone_vi/pyclone_vi_220610.tsv"
# OUTPUT_H5="/data/project/Alzheimer/EM_cluster/old/pilot/04.EM_input/pyclone_vi/pyclone_vi_220610.h5"
# OUTPUT_TSV="/data/project/Alzheimer/EM_cluster/old/pilot/04.EM_input/pyclone_vi/pyclone_vi_output_220610.tsv"

INPUT_TSV=$1
OUTPUT_H5=$2
OUTPUT_TSV=$3


source /home/goldpm1/.bashrc
conda activate cnvpytor
module switch HDF5/1.10.7 HDF5/1.12.0

pyclone-vi fit -i ${INPUT_TSV} -o ${OUTPUT_H5} -c 10 -d beta-binomial -r 10   &> ${OUTPUT_TSV%/*}"/log"
pyclone-vi write-results-file -i ${OUTPUT_H5} -o ${OUTPUT_TSV}                       &>> ${OUTPUT_TSV%/*}"/log"

# pyclone vi는  python3, conda activate pyclone,  module add HDF5/1.12.0 을 해야 한다
