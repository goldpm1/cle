#!/bin/bash
#$ -cwd
#$ -S /bin/bash

if ! options=$(getopt -o h --long SCRIPT_DIR:,NUM_CLONE:,NUM_BLOCK:,NUM_MUTATION:,FP_RATIO:,DEPTH_MEAN:,DEPTH_SD:,DEPTH_CUTOFF:,INPUT_TSV:,NPVAF_DIR:,BENCHMARK_I:,SIMDATA:, -- "$@")
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
        --NUM_CLONE)
            NUM_CLONE=$2
        shift 2 ;;
        --NUM_BLOCK)
            NUM_BLOCK=$2
        shift 2 ;;
        --NUM_MUTATION)
            NUM_MUTATION=$2
        shift 2 ;;
        --FP_RATIO)
            FP_RATIO=$2
        shift 2 ;;
        --DEPTH_MEAN)
            DEPTH_MEAN=$2
        shift 2 ;;
        --DEPTH_SD)
            DEPTH_SD=$2
        shift 2 ;;
        --DEPTH_CUTOFF)
            DEPTH_CUTOFF=$2
        shift 2 ;;
        --INPUT_TSV)
            INPUT_TSV=$2
        shift 2 ;;
        --NPVAF_DIR)
            NPVAF_DIR=$2
        shift 2 ;;
        --BENCHMARK_I)
            BENCHMARK_I=$2
        shift 2 ;;
        --SIMDATA)
            SIMDATA=$2
        shift 2 ;;
        --)
            shift
            break
    esac
done


echo -e python3 ${SCRIPT_DIR}"/1.SimData_pipe0_preparation.py" --NUM_CLONE ${NUM_CLONE} --NUM_BLOCK ${NUM_BLOCK} --NUM_MUTATION ${NUM_MUTATION} --FP_RATIO ${FP_RATIO} --DEPTH_MEAN ${DEPTH_MEAN} --DEPTH_SD ${DEPTH_SD} --DEPTH_CUTOFF ${DEPTH_CUTOFF}  --INPUT_TSV ${INPUT_TSV} --NPVAF_DIR ${NPVAF_DIR} --BENCHMARK_I ${BENCHMARK_I} --SimData ${SIMDATA}

python3 ${SCRIPT_DIR}"/1.SimData_pipe0_preparation.py" --NUM_CLONE ${NUM_CLONE} --NUM_BLOCK ${NUM_BLOCK} --NUM_MUTATION ${NUM_MUTATION} --FP_RATIO ${FP_RATIO} --DEPTH_MEAN ${DEPTH_MEAN} --DEPTH_SD ${DEPTH_SD} --DEPTH_CUTOFF ${DEPTH_CUTOFF} --INPUT_TSV ${INPUT_TSV} --NPVAF_DIR ${NPVAF_DIR} --BENCHMARK_I ${BENCHMARK_I} --SimData ${SIMDATA}
