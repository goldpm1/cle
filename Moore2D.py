import datapreparationold, datapreparation220712, comb, extract, scoring, boundaryclone, graph, phylogeny
import pyclonesim, sciclonesim
import EMhard, EMsoft, EMhard
import scoring, filetype
import visualizationsingle, visualizationpair, visualizationsinglesoft
import numpy as np
import pandas as pd
from kneed import KneeLocator
import glob, subprocess, os


# 2D 가능한 모든 Moore dataset을 돈다


inputdf = pd.read_csv("/data/project/Alzheimer/EM_cluster/EM_input/summary/Moore_2_sample.txt", sep = "\t")
TISSUE_set = set([])

for k in range (inputdf.shape[0]):
    DONOR, TISSUE, SAMPLE = inputdf.iloc[k]["DONOR"], inputdf.iloc[k]["TISSUE"], inputdf.iloc[k]["SAMPLE"]
    INPUT_TSV = "/".join(["/data/project/Alzheimer/EM_cluster/EM_input/Moore_2_sample", DONOR, TISSUE, SAMPLE+"_input.txt"])
    
   
    kwargs = {"INPUT_TSV" : INPUT_TSV,  "NUM_CLONE_TRIAL_START" : 1, "NUM_CLONE_TRIAL_END" : 5, "NUM_CLONE_TRIAL_FORCE" : 1,
                "RANDOM_PICK":50, "AXIS_RATIO":0, "PARENT_RATIO": 0, "FP_RATIO":0, "FP_2D" : "False", "TRIAL_NO" : 4, "DEPTH_CUTOFF" : 10, "VERBOSE" : 1,  "MIN_CLUSTER_SIZE" : 5, "NUM_PARENT" : 0, 
                "NPVAF_DIR" : "", "MYEM_DIR" : "", "SCICLONE_DIR" : "", "PYCLONE_DIR" : "", "PYCLONEVI_DIR" : "",  "COMBINED_OUTPUT_DIR" : "", "BENCHMARK_NO" : 10, 
                "KMEANS_CLUSTERNO" : 10, "RANDOM_SEED" : 1, "MODE" : "Both", "SCORING" : False }

    INPUT_TSV = kwargs["INPUT_TSV"]
    INPUT_FILETYPE, NUM_BLOCK = filetype.main (INPUT_TSV)
    kwargs["NUM_BLOCK_INPUT"], kwargs["NUM_BLOCK"] = NUM_BLOCK, NUM_BLOCK
    SAMPLENAME = INPUT_TSV.split("/")[-1].split(".")[0]     
    kwargs["SAMPLENAME"] = SAMPLE      # 'C10_H9'

    kwargs["NPVAF_DIR"] = "/data/project/Alzheimer/YSscript/EM_MRS/data/npvaf/Moore2D/" + TISSUE + "/" + DONOR + "_" + SAMPLE
    kwargs["MYEM_DIR"] = "/data/project/Alzheimer/YSscript/EM_MRS/data/MyEM/Moore2D/" + TISSUE + "/" + DONOR + "_" + SAMPLE
    kwargs["SCICLONE_DIR"] = "/data/project/Alzheimer/YSscript/EM_MRS/data/sciclone/Moore2D/" + TISSUE + "/" + DONOR + "-" + SAMPLE
    kwargs["PYCLONE_DIR"] = "/data/project/Alzheimer/YSscript/EM_MRS/data/pyclone/Moore2D/" + TISSUE + "/" + DONOR + "-" + SAMPLE
    kwargs["PYCLONEVI_DIR"] = "/data/project/Alzheimer/YSscript/EM_MRS/data/pyclone-vi/Moore2D/" + TISSUE  + "/" + DONOR + "-" + SAMPLE
    kwargs["COMBINED_OUTPUT_DIR"] = "/data/project/Alzheimer/YSscript/EM_MRS/data/combinedoutput/Moore2D/" + TISSUE + "/" + DONOR + "-" + SAMPLE

    #print ("SAMPLENAME : {}\nDONOR : {}\nSAMPLE : {}\nTISSUE : {}".format (SAMPLENAME, DONOR, SAMPLE, TISSUE) )

    os.system ("mkdir -p " + kwargs["NPVAF_DIR"])   # 출력 디렉토리 만들기
    os.system ("mkdir -p " + kwargs["MYEM_DIR"])   # 출력 디렉토리 만들기
    os.system ("mkdir -p " + kwargs["SCICLONE_DIR"])   # 출력 디렉토리 만들기
    os.system ("mkdir -p " + kwargs["PYCLONE_DIR"])   # 출력 디렉토리 만들기
    os.system ("mkdir -p " + kwargs["PYCLONEVI_DIR"])   # 출력 디렉토리 만들기
    os.system ("mkdir -p " + kwargs["COMBINED_OUTPUT_DIR"])   # 출력 디렉토리 만들기


    logPath = "/data/project/Alzheimer/YSscript/EM_MRS/log/Moore2D/" + TISSUE + "/" + DONOR + "-" + SAMPLE
    os.system ("rm -rf " + logPath)
    os.system ("mkdir -p " + logPath)

    command = " ".join ( [ "qsub -pe smp 1", "-e", logPath, "-o", logPath, "-N", TISSUE + "_" + DONOR + "_" + SAMPLE, 
                                        "benchmark_pipe1.sh",  
                                        str(INPUT_TSV),  str(kwargs["MODE"]),  str(kwargs["NUM_CLONE_TRIAL_START"]),  str(kwargs["NUM_CLONE_TRIAL_END"]),  str(kwargs["NUM_CLONE_TRIAL_FORCE"]),
                                        str(kwargs["RANDOM_PICK"]), str(kwargs["AXIS_RATIO"]),  str(kwargs["PARENT_RATIO"]),  str(kwargs["NUM_PARENT"]),  str(kwargs["FP_RATIO"]),  str(kwargs["FP_2D"]),
                                        str(kwargs["TRIAL_NO"]), str(kwargs["DEPTH_CUTOFF"]),  str(kwargs["MIN_CLUSTER_SIZE"]),  str(kwargs["VERBOSE"]),
                                        str(kwargs["KMEANS_CLUSTERNO"]),  str(kwargs["RANDOM_SEED"]), str(kwargs["SAMPLENAME"]), str(kwargs["BENCHMARK_NO"]), 
                                        str(kwargs["NPVAF_DIR"]), str(kwargs["MYEM_DIR"]), str(kwargs["SCICLONE_DIR"]), str(kwargs["PYCLONE_DIR"]), str(kwargs["PYCLONEVI_DIR"]) ,str(kwargs["COMBINED_OUTPUT_DIR"]), str(kwargs["MODE"]), str(kwargs["SCORING"])    ] )
    
    os.system (command)

    # if k >= 5:
    #     break