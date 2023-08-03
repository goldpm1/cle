#!/bin/bash
#$ -cwd
#$ -S /bin/bash
#$ -l h=!('compute15'|compute16')

if ! options=$(getopt -o h --long SCRIPT_DIR:,INPUT_TSV:,NPVAF_DIR:,SIMPLE_KMEANS_DIR:,CLEMENT_DIR:,SCICLONE_DIR:,PYCLONEVI_DIR:,QUANTUMCLONE_DIR:,COMBINED_OUTPUT_DIR:,NUM_CLONE_TRIAL_START:,NUM_CLONE_TRIAL_END:,FP_RATIO:,MAXIMUM_NUM_PARENT:,DEPTH_CUTOFF:,VERBOSE:,TRIAL_NO:,RANDOM_SEED:,SCORING:,MAKEONE_STRICT:,MODE: -- "$@")
then
    echo "ERROR: invalid options"
    exit 1
fi

eval set -- $options

while true; do
    case "$1" in
        -h|--help)
            echo "Usage"
        shift ;;
        --SCRIPT_DIR)
            SCRIPT_DIR=$2
        shift 2 ;;
        --INPUT_TSV)
            INPUT_TSV=$2
        shift 2 ;;
        --NPVAF_DIR)
            NPVAF_DIR=$2
        shift 2 ;;
        --SIMPLE_KMEANS_DIR)
            SIMPLE_KMEANS_DIR=$2
        shift 2 ;;
        --CLEMENT_DIR)
            CLEMENT_DIR=$2
        shift 2 ;;
        --SCICLONE_DIR)
            SCICLONE_DIR=$2
        shift 2 ;;
        --PYCLONEVI_DIR)
            PYCLONEVI_DIR=$2
        shift 2 ;;
        --QUANTUMCLONE_DIR)
            QUANTUMCLONE_DIR=$2
        shift 2 ;;
        --COMBINED_OUTPUT_DIR)
            COMBINED_OUTPUT_DIR=$2
        shift 2 ;;
        --NUM_CLONE_TRIAL_START)
            NUM_CLONE_TRIAL_START=$2
        shift 2 ;;
        --NUM_CLONE_TRIAL_END)
            NUM_CLONE_TRIAL_END=$2
        shift 2 ;;
        --FP_RATIO)
            FP_RATIO=$2
        shift 2 ;;
        --MAXIMUM_NUM_PARENT)
            MAXIMUM_NUM_PARENT=$2
        shift 2 ;;
        --DEPTH_CUTOFF)
            DEPTH_CUTOFF=$2
        shift 2 ;;
        --VERBOSE)
            VERBOSE=$2
        shift 2 ;;
        --TRIAL_NO)
            TRIAL_NO=$2
        shift 2 ;;
        --RANDOM_SEED)
            RANDOM_SEED=$2
        shift 2 ;;
        --SCORING)
            SCORING=$2
        shift 2 ;;
        --MAKEONE_STRICT)
            MAKEONE_STRICT=$2
        shift 2 ;;
        --MODE)
            MODE=$2
        shift 2 ;;
        --)
            shift
            break
    esac
done



echo  python3 ${SCRIPT_DIR}"/CLEMENT_bm.py"    --INPUT_TSV ${INPUT_TSV}    --NPVAF_DIR ${NPVAF_DIR} --SIMPLE_KMEANS_DIR ${SIMPLE_KMEANS_DIR} --CLEMENT_DIR ${CLEMENT_DIR} --SCICLONE_DIR ${SCICLONE_DIR}    --PYCLONEVI_DIR ${PYCLONEVI_DIR} --QUANTUMCLONE_DIR ${QUANTUMCLONE_DIR}    --COMBINED_OUTPUT_DIR ${COMBINED_OUTPUT_DIR}      --NUM_CLONE_TRIAL_START ${NUM_CLONE_TRIAL_START}     --NUM_CLONE_TRIAL_END ${NUM_CLONE_TRIAL_END}  --FP_RATIO ${FP_RATIO}    --MAXIMUM_NUM_PARENT ${MAXIMUM_NUM_PARENT} --DEPTH_CUTOFF ${DEPTH_CUTOFF}     --VERBOSE ${VERBOSE}     --TRIAL_NO ${TRIAL_NO}   --RANDOM_SEED ${RANDOM_SEED}     --SCORING ${SCORING}   --MODE ${MODE}  --MAKEONE_STRICT ${MAKEONE_STRICT}

python3 ${SCRIPT_DIR}"/CLEMENT_bm.py"    --INPUT_TSV ${INPUT_TSV}    --NPVAF_DIR ${NPVAF_DIR} --SIMPLE_KMEANS_DIR ${SIMPLE_KMEANS_DIR} --CLEMENT_DIR ${CLEMENT_DIR} --SCICLONE_DIR ${SCICLONE_DIR}    --PYCLONEVI_DIR ${PYCLONEVI_DIR} --QUANTUMCLONE_DIR ${QUANTUMCLONE_DIR}    --COMBINED_OUTPUT_DIR ${COMBINED_OUTPUT_DIR}      --NUM_CLONE_TRIAL_START ${NUM_CLONE_TRIAL_START}     --NUM_CLONE_TRIAL_END ${NUM_CLONE_TRIAL_END}  --FP_RATIO ${FP_RATIO}    --MAXIMUM_NUM_PARENT ${MAXIMUM_NUM_PARENT} --DEPTH_CUTOFF ${DEPTH_CUTOFF}     --VERBOSE ${VERBOSE}     --TRIAL_NO ${TRIAL_NO}   --RANDOM_SEED ${RANDOM_SEED}     --SCORING ${SCORING}   --MODE ${MODE} --MAKEONE_STRICT ${MAKEONE_STRICT}

