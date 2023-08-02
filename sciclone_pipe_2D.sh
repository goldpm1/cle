#!/bin/bash
#$ -cwd
#$ -S /bin/bash
#$ -l h=!('compute15'|compute16')

SCRIPT_DIR=$1
INPUT_First=$2
INPUT_Second=$3
OUTPUT_DIR=$4

mkdir -m 774 -p $OUTPUT_DIR

#1. Rscript를 돌리되 log 파일 (stdout ,stderr)를 OUTPUT_DIR 안에 넣자

/opt/Yonsei/R/4.2.0/bin/Rscript ${SCRIPT_DIR}"/sciclone_run_2D.R" ${INPUT_First} ${INPUT_Second} $OUTPUT_DIR   &> $OUTPUT_DIR"/log"
