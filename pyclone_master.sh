#!/bin/bash
#$ -cwd
#$ -S /bin/bash

# if ! options=$(getopt -o h --long output_filename:,RANDOM_PICK:, -- "$@")
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
#         --output_filename)
#             output_filename=$2
#         shift 2 ;;
#         --RANDOM_PICK)
#             RANDOM_PICK=$2
#         shift 2 ;;
#         --)
#             shift
#             break
#     esac
# done

# echo -e ${output_filename}"\n"${RANDOM_PICK}


source /home/goldpm1/.bashrc
conda activate pyclone

logPath="/data/project/Alzheimer/YSscript/EM_MRS/log"

qsub -pe smp 2 -e $logPath -o $logPath  -N "pyclone" /data/project/Alzheimer/YSscript/EM_MRS/pyclone_pipe.sh \
--output_filename "/data/project/Alzheimer/YSscript/EM_MRS/MRS_pyclone.jpg" --RANDOM_PICK 2000