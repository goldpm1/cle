#!/bin/bash
#$ -cwd
#$ -S /bin/bash



INPUT_DIR=$1
OUTPUT_FILENAME=$2
BENCHMARK_NO=$3




python3 "2.CellData_pipe3_benchmark.py" \
 --INPUT_DIR ${INPUT_DIR} \
 --OUTPUT_FILENAME ${OUTPUT_FILENAME} \
 --BENCHMARK_NO ${BENCHMARK_NO} 


 ###########$ -l h=!('compute15'|compute16')