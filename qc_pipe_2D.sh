#!/bin/bash
#$ -cwd
#$ -S /bin/bash
#$ -l h=!('compute15'|compute16')

INPUT_First=$1
INPUT_Second=$2
OUTPUT_DIR=$3

mkdir -m 774 -p $OUTPUT_DIR

#1. Rscript를 돌리되 log 파일 (stdout ,stderr)를 OUTPUT_DIR 안에 넣자

/opt/Yonsei/R/4.2.0/bin/Rscript /data/project/Alzheimer/YSscript/EM_MRS/qc_run_2D.R ${INPUT_First} ${INPUT_Second} $OUTPUT_DIR  &> $OUTPUT_DIR"/log"
