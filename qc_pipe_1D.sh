#!/bin/bash
#$ -cwd
#$ -S /bin/bash
#$ -l h=!('compute15'|compute16')

SCRIPT_DIR=$1
INPUT_First=$2
OUTPUT_DIR=$3

mkdir -m 774 -p $OUTPUT_DIR

#1. Rscript를 돌리되 log 파일 (stdout ,stderr)를 OUTPUT_DIR 안에 넣자

/opt/Yonsei/R/4.2.0/bin/Rscript ${SCRIPT_DIR}"/qc_run_1D.R" ${INPUT_First} $OUTPUT_DIR  &> $OUTPUT_DIR"/log"
