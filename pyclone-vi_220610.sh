#!/bin/bash
#$ -cwd
#$ -S /bin/bash

INPUT_TSV="/data/project/Alzheimer/EM_cluster/pilot/04.EM_input/pyclone_vi/pyclone_vi_220610.tsv"
OUTPUT_H5="/data/project/Alzheimer/EM_cluster/pilot/04.EM_input/pyclone_vi/pyclone_vi_220610.h5"
OUTPUT_TSV="/data/project/Alzheimer/EM_cluster/pilot/04.EM_input/pyclone_vi/pyclone_vi_output_220610.tsv"

pyclone-vi fit -i ${INPUT_TSV} -o ${OUTPUT_H5} -c 10 -d beta-binomial -r 10

echo -e "\n\nfit done\n\n"

#pyclone-vi write-results-file -i ${OUTPUT_H5} -o ${OUTPUT_TSV}
#python3 pyclone-vi_output_220610.py --in_file ${OUTPUT_H5} --out_file ${OUTPUT_TSV}
echo -e "\n\nresult print done"