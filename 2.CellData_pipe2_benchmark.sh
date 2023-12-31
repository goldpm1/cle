#!/bin/bash
#$ -cwd
#$ -S /bin/bash

if ! options=$(getopt -o h --long SCRIPT_DIR:,INPUT_DIR:,LOG_DIR:,SAMPLENAME:,CONDITIONNAME:,BENCHMARK_START:,BENCHMARK_END:,OUTPUT_JPG:,OUTPUT_TTEST:,OUTPUT_RESULT_TABLE:, -- "$@")
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
        --INPUT_DIR)
            INPUT_DIR=$2
        shift 2 ;;
        --LOG_DIR)
            LOG_DIR=$2
        shift 2 ;;
        --SAMPLENAME)
            SAMPLENAME=$2
        shift 2 ;;
        --CONDITIONNAME)
            CONDITIONNAME=$2
        shift 2 ;;
        --BENCHMARK_START)
            BENCHMARK_START=$2
        shift 2 ;;
        --BENCHMARK_END)
            BENCHMARK_END=$2
        shift 2 ;;
        --OUTPUT_JPG)
            OUTPUT_JPG=$2
        shift 2 ;;
        --OUTPUT_TTEST)
            OUTPUT_TTEST=$2
        shift 2 ;;
        --OUTPUT_RESULT_TABLE)
            OUTPUT_RESULT_TABLE=$2
        shift 2 ;;
        --)
            shift
            break
    esac
done

echo -e python3 ${SCRIPT_DIR}"/2.CellData_pipe2_benchmark.py" \
 --INPUT_DIR ${INPUT_DIR}  \
 --LOG_DIR ${LOG_DIR}  \
 --SAMPLENAME ${SAMPLENAME} \
 --CONDITIONNAME ${CONDITIONNAME} \
 --BENCHMARK_START ${BENCHMARK_START} \
 --BENCHMARK_END ${BENCHMARK_END} \
 --OUTPUT_JPG ${OUTPUT_JPG}  --OUTPUT_TTEST ${OUTPUT_TTEST}  --OUTPUT_RESULT_TABLE ${OUTPUT_RESULT_TABLE}


python3 ${SCRIPT_DIR}"/2.CellData_pipe2_benchmark.py" \
 --INPUT_DIR ${INPUT_DIR}  \
 --LOG_DIR ${LOG_DIR}  \
 --SAMPLENAME ${SAMPLENAME} \
 --CONDITIONNAME ${CONDITIONNAME} \
 --BENCHMARK_START ${BENCHMARK_START} \
 --BENCHMARK_END ${BENCHMARK_END} \
 --OUTPUT_JPG ${OUTPUT_JPG}  --OUTPUT_TTEST ${OUTPUT_TTEST}  --OUTPUT_RESULT_TABLE ${OUTPUT_RESULT_TABLE}


 ###########$ -l h=!('compute15'|compute16')