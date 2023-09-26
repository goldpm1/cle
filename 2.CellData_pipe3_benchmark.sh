#!/bin/bash
#$ -cwd
#$ -S /bin/bash


if ! options=$(getopt -o h --long SCRIPT_DIR:,INPUT_DIR:,CONDITIONNAME:,BENCHMARK_START:,BENCHMARK_END:,OUTPUT_MS_JPG:,OUTPUT_EC_JPG:,OUTPUT_FINAL_JPG:,OUTPUT_FINAL_TABLE:, -- "$@")
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
        --CONDITIONNAME)
            CONDITIONNAME=$2
        shift 2 ;;
        --BENCHMARK_START)
            BENCHMARK_START=$2
        shift 2 ;;
        --BENCHMARK_END)
            BENCHMARK_END=$2
        shift 2 ;;
        --OUTPUT_MS_JPG)
            OUTPUT_MS_JPG=$2
        shift 2 ;;
        --OUTPUT_EC_JPG)
            OUTPUT_EC_JPG=$2
        shift 2 ;;
        --OUTPUT_FINAL_JPG)
            OUTPUT_FINAL_JPG=$2
        shift 2 ;;
        --OUTPUT_FINAL_TABLE)
            OUTPUT_FINAL_TABLE=$2
        shift 2 ;;
        --)
            shift
            break
    esac
done

echo -e  python3 ${SCRIPT_DIR}"/2.CellData_pipe3_benchmark.py" \
 --INPUT_DIR ${INPUT_DIR} \
 --CONDITIONNAME ${CONDITIONNAME} \
 --BENCHMARK_START ${BENCHMARK_START} \
 --BENCHMARK_END ${BENCHMARK_END} \
 --OUTPUT_MS_JPG ${OUTPUT_MS_JPG} \
 --OUTPUT_EC_JPG ${OUTPUT_EC_JPG} \
 --OUTPUT_FINAL_JPG ${OUTPUT_FINAL_JPG} \
 --OUTPUT_FINAL_TABLE ${OUTPUT_FINAL_TABLE}


python3 ${SCRIPT_DIR}"/2.CellData_pipe3_benchmark.py" \
 --INPUT_DIR ${INPUT_DIR} \
 --CONDITIONNAME ${CONDITIONNAME} \
 --BENCHMARK_START ${BENCHMARK_START} \
 --BENCHMARK_END ${BENCHMARK_END} \
 --OUTPUT_MS_JPG ${OUTPUT_MS_JPG} \
 --OUTPUT_EC_JPG ${OUTPUT_EC_JPG} \
 --OUTPUT_FINAL_JPG ${OUTPUT_FINAL_JPG} \
 --OUTPUT_FINAL_TABLE ${OUTPUT_FINAL_TABLE}


 ###########$ -l h=!('compute15'|compute16')