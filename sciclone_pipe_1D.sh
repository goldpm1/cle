#!/bin/bash
#$ -cwd
#$ -S /bin/bash
#$ -l h=!('compute15'|compute16')

INPUT_First=$1
OUTPUT_DIR=$2

mkdir -m 774 -p $OUTPUT_DIR

#1. Rscript를 돌리되 log 파일 (stdout ,stderr)를 OUTPUT_DIR 안에 넣자
Rscript /data/project/Alzheimer/YSscript/EM_MRS/sciclone_run_1D.R ${INPUT_First} $OUTPUT_DIR &> $OUTPUT_DIR"/log"
