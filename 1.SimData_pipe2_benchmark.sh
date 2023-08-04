#!/bin/bash
#$ -cwd
#$ -S /bin/bash
#$ -l h=!('compute15'|compute16')


if ! options=$(getopt -o h --long SCRIPT_DIR:,COMBINED_OUTPUT_DIR:,SAMPLENAME:,BENCHMARK_START:,BENCHMARK_END:,OUTPUT_TTEST:,OUTPUT_JPG:,FP_RATIO:, -- "$@")
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
        --COMBINED_OUTPUT_DIR)
            COMBINED_OUTPUT_DIR=$2
        shift 2 ;;
        --SAMPLENAME)
            SAMPLENAME=$2
        shift 2 ;;
        --BENCHMARK_START)
            BENCHMARK_START=$2
        shift 2 ;;
        --BENCHMARK_END)
            BENCHMARK_END=$2
        shift 2 ;;
        --OUTPUT_TTEST)
            OUTPUT_TTEST=$2
        shift 2 ;;
        --OUTPUT_JPG)
            OUTPUT_JPG=$2
        shift 2 ;;
        --FP_RATIO)
            FP_RATIO=$2
        shift 2 ;;
        --)
            shift
            break
    esac
done



echo -e python3 ${SCRIPT_DIR}"/1.SimData_pipe2_benchmark.py"  --COMBINED_OUTPUT_DIR ${COMBINED_OUTPUT_DIR} --SAMPLENAME ${SAMPLENAME}  --BENCHMARK_START ${BENCHMARK_START}  --BENCHMARK_END ${BENCHMARK_END}   --OUTPUT_TTEST ${OUTPUT_TTEST}     --OUTPUT_JPG ${OUTPUT_JPG} --FP_RATIO ${FP_RATIO}
python3 ${SCRIPT_DIR}"/1.SimData_pipe2_benchmark.py"  --COMBINED_OUTPUT_DIR ${COMBINED_OUTPUT_DIR} --SAMPLENAME ${SAMPLENAME}  --BENCHMARK_START ${BENCHMARK_START}  --BENCHMARK_END ${BENCHMARK_END}   --OUTPUT_TTEST ${OUTPUT_TTEST}     --OUTPUT_JPG ${OUTPUT_JPG} --FP_RATIO ${FP_RATIO}