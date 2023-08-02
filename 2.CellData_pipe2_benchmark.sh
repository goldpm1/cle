#!/bin/bash
#$ -cwd
#$ -S /bin/bash



INPUT_TSV=$1
MODE=$2
NUM_CLONE_TRIAL_START=$3
NUM_CLONE_TRIAL_END=$4
NUM_CLONE_TRIAL_FORCE=$5
RANDOM_PICK=$6
AXIS_RATIO=$7
PARENT_RATIO=$8
NUM_PARENT=$9
FP_RATIO=${10}
FP_USEALL=${11}
TRIAL_NO=${12}
DEPTH_CUTOFF=${13}
MIN_CLUSTER_SIZE=${14}
VERBOSE=${15}
KMEANS_CLUSTERNO=${16}
RANDOM_SEED=${17}
SAMPLENAME=${18}
BENCHMARK_NO=${19}
NPVAF_DIR=${20}
CLEMENT_DIR=${21}
SCICLONE_DIR=${22}
PYCLONE_DIR=${23}
PYCLONEVI_DIR=${24}
QUANTUMCLONE_DIR=${25}
COMBINED_OUTPUT_DIR=${26}
SCORING=${27}
DIMENSION=${28}
MAKEONE_STRICT=${29}


#2. CLEMENT_hard, pyclonevi, sciclone_result 데이터를 읽은 후 시각화 및 저장

echo ${COMBINED_OUTPUT_DIR%/*}"/benchmark.jpg"
echo -e "FP_RATIO = "${FP_RATIO}

python3 2.CellData_pipe2_benchmark.py \
 --INPUT_DIR "/data/project/Alzheimer/YSscript/EM_MRS/data/combinedoutput/MRS_"${DIMENSION} \
 --SAMPLENAME ${SAMPLENAME} \
 --BENCHMARK_NO ${BENCHMARK_NO} \
 --RANDOM_PICK ${RANDOM_PICK} \
 --NUM_PARENT ${NUM_PARENT}  \
 --FP_RATIO ${FP_RATIO} \
 --AXIS_RATIO ${AXIS_RATIO}  \
 --OUTPUT_JPG ${COMBINED_OUTPUT_DIR%/*}"/benchmark.jpg"  --OUTPUT_TTEST ${COMBINED_OUTPUT_DIR%/*}"/ttest.txt"  \
 --MAKEONE_STRICT ${MAKEONE_STRICT}


 ###########$ -l h=!('compute15'|compute16')