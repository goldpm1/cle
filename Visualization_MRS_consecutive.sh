#!/bin/bash
#$ -cwd
#$ -S /bin/bash

if ! options=$(getopt -o h --long TRIAL_NUM:,RANDOM_PICK:,FP_RATIO:,OUTPUT_SCORE:,FIG_OUTPUT:, -- "$@")
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
        --TRIAL_NUM)
            TRIAL_NUM=$2
        shift 2 ;;
        --RANDOM_PICK)
            RANDOM_PICK=$2
        shift 2 ;;
        --FP_RATIO)
            FP_RATIO=$2
        shift 2 ;;
        --FIG_OUTPUT)
            FIG_OUTPUT=$2
        shift 2 ;;
        --OUTPUT_SCORE)
            OUTPUT_SCORE=$2
        shift 2 ;;
        --)
            shift
            break
    esac
done

echo $TRIAL_NUM
echo $RANDOM_PICK
echo $FP_RATIO
echo $FIG_OUTPUT
echo $OUTPUT_SCORE

python3 EM_MRS.py --TRIAL_NUM ${TRIAL_NUM} --RANDOM_PICK ${RANDOM_PICK} --FP_RATIO ${FP_RATIO} --OUTPUT_SCORE "./output/"${FP_RATIO}"_score.txt" --FIG_OUTPUT "./output/"${FP_RATIO}"_fig.png"