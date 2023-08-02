#!/bin/bash

#PATH #

First='M1-1'
Second='M1-2'

mkdir -p /home/octo0410/script_organize/out/sci_clone
Log_PATH=/home/octo0410/script_organize/out/sci_clone

echo ${First}'_'${Second}
qsub -pe smp 2 -e $Log_PATH -o $Log_PATH -N 'sci.sub.'${First}'_'${Second} pipe_sciclone.sh ${First} ${Second}
