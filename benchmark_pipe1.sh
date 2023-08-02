#!/bin/bash
#$ -cwd
#$ -S /bin/bash

# if ! options=$(getopt -o h --long INPUT_TSV:,MODE:,NUM_BLOCK:,NUM_CLONE_TRIAL_START:,NUM_CLONE_TRIAL_END:,NUM_CLONE_TRIAL_FORCE:,RANDOM_PICK:,AXIS_RATIO:,PARENT_RATIO:,NUM_PARENT:,FP_RATIO:,FP_2D:,\
#                     TRIAL_NO:,DEPTH_CUTOFF:,VERBOSE:,MIN_CLUSTER_SIZE:,KMEANS_CLUSTERNO:,RANDOM_SEED:, -- "$@")
# then
#     echo "ERROR: invalid options"
#     exit 1
# fi

# eval set -- $options

# while true; do
#     case "$1" in
#         -h|--help)
#             echo "Usage"
#         shift ;;
#         --INPUT_TSV)
#             INPUT_TSV=$2
#         shift 2 ;;
#         --MODE)
#             MODE=$2
#         shift 2 ;;
#         --NUM_CLONE_TRIAL_START)
#             NUM_CLONE_TRIAL_START=$2
#         shift 2 ;;
#         --NUM_CLONE_TRIAL_END)
#             NUM_CLONE_TRIAL_END=$2
#         shift 2 ;;
#         --NUM_CLONE_TRIAL_FORCE)
#             NUM_CLONE_TRIAL_FORCE=$2
#         shift 2 ;;
#         --RANDOM_PICK)
#             RANDOM_PICK=$2
#         shift 2 ;;
#         --AXIS_RATIO)
#             AXIS_RATIO=$2
#         shift 2 ;;
#         --PARENT_RATIO)
#             PARENT_RATIO=$2
#         shift 2 ;;
#         --NUM_PARENT)
#             NUM_PARENT=$2
#         shift 2 ;;
#         --FP_RATIO)
#             FP_RATIO=$2
#         shift 2 ;;
#         --FP_2D)
#             FP_2D=$2
#         shift 2 ;;
#         --TRIAL_NO)
#             TRIAL_NO=$2
#         shift 2 ;;
#         --TRIAL_NO)
#             TRIAL_NO=$2
#         shift 2 ;;
#         --DEPTH_CUTOFF)
#             DEPTH_CUTOFF=$2
#         shift 2 ;;
#         --VERBOSE)
#             VERBOSE=$2
#         shift 2 ;;
#         --MIN_CLUSTER_SIZE)
#             MIN_CLUSTER_SIZE=$2
#         shift 2 ;;
#         --KMEANS_CLUSTERNO)
#             KMEANS_CLUSTERNO=$2
#         shift 2 ;;
#         --RANDOM_SEED)
#             RANDOM_SEED=$2
#         shift 2 ;;
#         --)
#             shift
#             break
#     esac
# done

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
FP_2D=${11}
TRIAL_NO=${12}
DEPTH_CUTOFF=${13}
MIN_CLUSTER_SIZE=${14}
VERBOSE=${15}
KMEANS_CLUSTERNO=${16}
RANDOM_SEED=${17}
SAMPLENAME=${18}
BENCHMARK_NO=${19}
NPVAF_DIR=${20}
MYEM_DIR=${21}
SCICLONE_DIR=${22}
PYCLONE_DIR=${23}
PYCLONEVI_DIR=${24}
COMBINED_OUTPUT_DIR=${25}
MODE=${26}
SCORING=${27}



#1. 각 RANDOM_SEED에 대해서 tool을 돌린다
echo -e python3 EMhybrid.py --INPUT_TSV $INPUT_TSV --NUM_CLONE_TRIAL_START $NUM_CLONE_TRIAL_START --NUM_CLONE_TRIAL_END $NUM_CLONE_TRIAL_END --NUM_CLONE_TRIAL_FORCE $NUM_CLONE_TRIAL_FORCE \
    --RANDOM_PICK ${RANDOM_PICK} --AXIS_RATIO $AXIS_RATIO --PARENT_RATIO ${PARENT_RATIO} --NUM_PARENT ${NUM_PARENT} --FP_RATIO $FP_RATIO --FP_2D ${FP_2D}  \
    --TRIAL_NO $TRIAL_NO --DEPTH_CUTOFF ${DEPTH_CUTOFF} --MIN_CLUSTER_SIZE ${MIN_CLUSTER_SIZE}  --VERBOSE $VERBOSE --KMEANS_CLUSTERNO ${KMEANS_CLUSTERNO} --RANDOM_SEED ${RANDOM_SEED} \
    --NPVAF_DIR ${NPVAF_DIR} --MYEM_DIR ${MYEM_DIR} --SCICLONE_DIR ${SCICLONE_DIR} --PYCLONE_DIR ${PYCLONE_DIR} --PYCLONEVI_DIR ${PYCLONEVI_DIR} --COMBINED_OUTPUT_DIR ${COMBINED_OUTPUT_DIR} \
    --MODE ${MODE} --SCORING ${SCORING}



python3 EMhybrid.py --INPUT_TSV $INPUT_TSV --NUM_CLONE_TRIAL_START $NUM_CLONE_TRIAL_START --NUM_CLONE_TRIAL_END $NUM_CLONE_TRIAL_END --NUM_CLONE_TRIAL_FORCE $NUM_CLONE_TRIAL_FORCE \
    --RANDOM_PICK ${RANDOM_PICK} --AXIS_RATIO $AXIS_RATIO --PARENT_RATIO ${PARENT_RATIO} --NUM_PARENT ${NUM_PARENT} --FP_RATIO $FP_RATIO --FP_2D ${FP_2D}  \
    --TRIAL_NO $TRIAL_NO --DEPTH_CUTOFF ${DEPTH_CUTOFF} --MIN_CLUSTER_SIZE ${MIN_CLUSTER_SIZE}  --VERBOSE $VERBOSE --KMEANS_CLUSTERNO ${KMEANS_CLUSTERNO} --RANDOM_SEED ${RANDOM_SEED} \
    --NPVAF_DIR ${NPVAF_DIR} --MYEM_DIR ${MYEM_DIR} --SCICLONE_DIR ${SCICLONE_DIR} --PYCLONE_DIR ${PYCLONE_DIR} --PYCLONEVI_DIR ${PYCLONEVI_DIR} --COMBINED_OUTPUT_DIR ${COMBINED_OUTPUT_DIR} \
    --MODE ${MODE} --SCORING ${SCORING}
