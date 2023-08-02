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
MAKEONE_STRICT=${28}
MAXIMUM_NUM_PARENT=${29}



#1. 각 RANDOM_SEED에 대해서 tool을 돌린다
echo -e python3 EMhybrid.py --INPUT_TSV $INPUT_TSV --NUM_CLONE_TRIAL_START $NUM_CLONE_TRIAL_START --NUM_CLONE_TRIAL_END $NUM_CLONE_TRIAL_END \
    --RANDOM_PICK ${RANDOM_PICK} --AXIS_RATIO $AXIS_RATIO --PARENT_RATIO ${PARENT_RATIO} --NUM_PARENT ${NUM_PARENT} --FP_RATIO $FP_RATIO --FP_USEALL ${FP_USEALL}  \
    --DEPTH_CUTOFF ${DEPTH_CUTOFF} --MIN_CLUSTER_SIZE ${MIN_CLUSTER_SIZE}  --KMEANS_CLUSTERNO ${KMEANS_CLUSTERNO} --RANDOM_SEED ${RANDOM_SEED} \
    --NPVAF_DIR ${NPVAF_DIR} --CLEMENT_DIR ${CLEMENT_DIR} --SCICLONE_DIR ${SCICLONE_DIR} --PYCLONE_DIR ${PYCLONE_DIR} --PYCLONEVI_DIR ${PYCLONEVI_DIR} --QUANTUMCLONE_DIR ${QUANTUMCLONE_DIR} --COMBINED_OUTPUT_DIR ${COMBINED_OUTPUT_DIR} \
    --MODE ${MODE} --SCORING ${SCORING} --MAKEONE_STRICT ${MAKEONE_STRICT} --MAXIMUM_NUM_PARENT ${MAXIMUM_NUM_PARENT}  --VERBOSE $VERBOSE  --TRIAL_NO $TRIAL_NO



python3 EMhybrid.py --INPUT_TSV $INPUT_TSV --NUM_CLONE_TRIAL_START $NUM_CLONE_TRIAL_START --NUM_CLONE_TRIAL_END $NUM_CLONE_TRIAL_END --NUM_CLONE_TRIAL_FORCE $NUM_CLONE_TRIAL_FORCE \
    --RANDOM_PICK ${RANDOM_PICK} --AXIS_RATIO $AXIS_RATIO --PARENT_RATIO ${PARENT_RATIO} --NUM_PARENT ${NUM_PARENT} --FP_RATIO $FP_RATIO --FP_USEALL ${FP_USEALL}  \
    --DEPTH_CUTOFF ${DEPTH_CUTOFF} --MIN_CLUSTER_SIZE ${MIN_CLUSTER_SIZE}  --VERBOSE $VERBOSE --KMEANS_CLUSTERNO ${KMEANS_CLUSTERNO} --RANDOM_SEED ${RANDOM_SEED} \
    --NPVAF_DIR ${NPVAF_DIR} --CLEMENT_DIR ${CLEMENT_DIR} --SCICLONE_DIR ${SCICLONE_DIR} --PYCLONE_DIR ${PYCLONE_DIR} --PYCLONEVI_DIR ${PYCLONEVI_DIR} --QUANTUMCLONE_DIR ${QUANTUMCLONE_DIR} --COMBINED_OUTPUT_DIR ${COMBINED_OUTPUT_DIR} \
    --MODE ${MODE} --SCORING ${SCORING} --MAKEONE_STRICT ${MAKEONE_STRICT} --MAXIMUM_NUM_PARENT ${MAXIMUM_NUM_PARENT}  --TRIAL_NO $TRIAL_NO


###########  #$ -l h=!('compute15'|compute16')