import numpy as np
import pandas as pd
import glob, subprocess, os, filetype


def out(command): 
    result = subprocess.run(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True, shell=True) 
    return result.stdout.rstrip("\n")


DIR = "/data/project/Alzheimer/EM_cluster/Moore_data/Donor"
DONOR_LIST = glob.glob(DIR + "/*")
DONOR_LIST = [i.split("/")[-1] for i in DONOR_LIST]
k = 0

for DONOR in DONOR_LIST:               # PD42566
    TISSUE_LIST = glob.glob (DIR + "/" + DONOR + "/*")
    TISSUE_LIST = [i.split("/")[-1] for i in TISSUE_LIST]
    for TISSUE in TISSUE_LIST:               # colon_crypt
        SAMPLE_LIST = glob.glob (DIR + "/" + DONOR + "/" + TISSUE + "/*")
        SAMPLE_LIST = [i.split("/")[-1].split(".")[0] for i in SAMPLE_LIST]
        for SAMPLE in SAMPLE_LIST:       # PD42566b_lo00_A7.txt
            INPUT_TSV = DIR + "/" + DONOR + "/" + TISSUE + "/" + SAMPLE + ".txt"
            
            if int(out ("wc -l  " + INPUT_TSV).split(" ")[0]) < 200:       # 200줄이 안 되는 파일은 넘어간다
                continue 

            #print ( INPUT_TSV,  out ("wc -l  " + INPUT_TSV).split(" ")[0] ) 

            kwargs = {"INPUT_TSV" : INPUT_TSV,  "NUM_CLONE_TRIAL_START" : 1, "NUM_CLONE_TRIAL_END" : 5, "NUM_CLONE_TRIAL_FORCE" : 1,
                        "RANDOM_PICK":100, "AXIS_RATIO":0, "PARENT_RATIO": 0, "FP_RATIO":0,  "FP_2D" : "False", "TRIAL_NO" : 4, "DEPTH_CUTOFF" : 10, "VERBOSE" : 1,  "MIN_CLUSTER_SIZE" : 5, "NUM_PARENT" : 0, 
                        "NPVAF_DIR" : "", "MYEM_DIR" : "", "SCICLONE_DIR" : "", "PYCLONE_DIR" : "", "PYCLONEVI_DIR" : "",  "COMBINED_OUTPUT_DIR" : "",
                        "KMEANS_CLUSTERNO" : 10, "RANDOM_SEED" : 1, "BENCHMARK_NO" : 10, "MODE" : "Both", "SCORING" : False }

            INPUT_TSV = kwargs["INPUT_TSV"]
            INPUT_FILETYPE, NUM_BLOCK = filetype.main (INPUT_TSV)
            kwargs["NUM_BLOCK_INPUT"], kwargs["NUM_BLOCK"] = NUM_BLOCK, NUM_BLOCK
            SAMPLENAME = INPUT_TSV.split("/")[-1].split(".")[0]     
            kwargs["SAMPLENAME"] = SAMPLE      # 'C10_H9'

            kwargs["NPVAF_DIR"] = "/data/project/Alzheimer/YSscript/EM_MRS/data/npvaf/Moore1D/" + TISSUE + "/" + DONOR + "_" + SAMPLE
            kwargs["MYEM_DIR"] = "/data/project/Alzheimer/YSscript/EM_MRS/data/MyEM/Moore1D/" + TISSUE + "/" + DONOR + "_" + SAMPLE
            kwargs["SCICLONE_DIR"] = "/data/project/Alzheimer/YSscript/EM_MRS/data/sciclone/Moore1D/" + TISSUE + "/" + DONOR + "-" + SAMPLE
            kwargs["PYCLONE_DIR"] = "/data/project/Alzheimer/YSscript/EM_MRS/data/pyclone/Moore1D/" + TISSUE + "/" + DONOR + "-" + SAMPLE
            kwargs["PYCLONEVI_DIR"] = "/data/project/Alzheimer/YSscript/EM_MRS/data/pyclone-vi/Moore1D/" + TISSUE  + "/" + DONOR + "-" + SAMPLE
            kwargs["COMBINED_OUTPUT_DIR"] = "/data/project/Alzheimer/YSscript/EM_MRS/data/combinedoutput/Moore1D/" + TISSUE + "/" + DONOR + "-" + SAMPLE

            os.system ("mkdir -p " + kwargs["NPVAF_DIR"])   # 출력 디렉토리 만들기
            os.system ("mkdir -p " + kwargs["MYEM_DIR"])   # 출력 디렉토리 만들기
            os.system ("mkdir -p " + kwargs["SCICLONE_DIR"])   # 출력 디렉토리 만들기
            os.system ("mkdir -p " + kwargs["PYCLONE_DIR"])   # 출력 디렉토리 만들기
            os.system ("mkdir -p " + kwargs["PYCLONEVI_DIR"])   # 출력 디렉토리 만들기
            os.system ("mkdir -p " + kwargs["COMBINED_OUTPUT_DIR"])   # 출력 디렉토리 만들기
            

            logPath = "/data/project/Alzheimer/YSscript/EM_MRS/log/Moore1D/" + TISSUE + "/" + DONOR + "-" + SAMPLE
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
            #print (command)
            k = k  + 1

    #         if k >= 20:
    #             break
    #     if k >= 20:
    #         break
    # if k >= 20:
    #     break