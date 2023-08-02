#!/bin/bash
#$ -cwd
#$ -S /bin/bash

if ! options=$(getopt -o h --long NUM_CLONE:,NUM_BLOCK:,NUM_MUTATION:,FP_RATIO:,LOWVAF_RATIO:,DEPTH_MEAN:,DEPTH_SD:,INPUT_TSV:,NPVAF_DIR:,BENCHMARK_I:,SIMDATA:, -- "$@")
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
        --LOWVAF_RATIO)
            LOWVAF_RATIO=$2
        shift 2 ;;
        --DEPTH_MEAN)
            DEPTH_MEAN=$2
        shift 2 ;;
        --DEPTH_SD)
            DEPTH_SD=$2
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



python3 1.SimData_pipe0_preparation.py --NUM_CLONE ${NUM_CLONE} --NUM_BLOCK ${NUM_BLOCK} --NUM_MUTATION ${NUM_MUTATION} --FP_RATIO ${FP_RATIO} --LOWVAF_RATIO ${LOWVAF_RATIO} --DEPTH_MEAN ${DEPTH_MEAN} --DEPTH_SD ${DEPTH_SD} --INPUT_TSV ${INPUT_TSV} --NPVAF_DIR ${NPVAF_DIR} --BENCHMARK_I ${BENCHMARK_I} --SimData ${SIMDATA}

                # command1 = " ".join(["python3 1.SimData_pipe0_preparation.py",
                #                     "--NUM_CLONE", str(NUM_CLONE),  "--NUM_BLOCK", str(NUM_BLOCK), "--NUM_MUTATION", str(NUM_MUTATION), "--FP_RATIO", str(kwargs["FP_RATIO"]),   "--LOWVAF_RATIO", str(kwargs["LOWVAF_RATIO"]), 
                #                     "--DEPTH_MEAN", str(kwargs["DEPTH_MEAN"]), "--DEPTH_SD", str(kwargs["DEPTH_SD"]),
                #                     "--INPUT_TSV", kwargs["INPUT_TSV"], "--NPVAF_DIR", kwargs["NPVAF_DIR"],
                #                     "--BENCHMARK_I", str(ii),
                #                     "--SimData", kwargs["SIMDATA"]
                #                     ])