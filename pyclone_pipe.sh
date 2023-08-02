#!/bin/bash
#$ -cwd
#$ -S /bin/bash

INPUT_COMMAND=$1

COMMAND=${INPUT_COMMAND//"&"/" "}
echo $COMMAND

#COMMAND="/home/goldpm1/miniconda3/envs/pyclone/bin/PyClone  run_analysis_pipeline --in_files /data/project/Alzheimer/YSscript/EM_MRS/data/pyclone/M2-2_M2-12_input/block0.tsv /data/project/Alzheimer/YSscript/EM_MRS/data/pyclone/M2-2_M2-12_input/block1.tsv --samples "block0" "block1" --working_dir  /data/project/Alzheimer/YSscript/EM_MRS/data/pyclone/M2-2_M2-12_input  --tumour_contents 1.0 1.0  --num_iters 1000 --max_clusters 7 --min_cluster_size 10"

source /home/goldpm1/.bashrc
conda activate pyclone

$COMMAND



echo "Pyclone done"


# python3 /data/project/Alzheimer/YSscript/EM_MRS/pyclonesim.py \
# --OUTPUT_FILENAME ${OUTPUT_FILENAME} \
# --RANDOM_PICK ${RANDOM_PICK}