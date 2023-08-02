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
    # if DONOR not in ["PD43850", "PD43851"]:
    #     continue

    for TISSUE in TISSUE_LIST:               # colon_crypt, pancreas_islet
        SAMPLE_LIST = glob.glob (DIR + "/" + DONOR + "/" + TISSUE + "/*")
        SAMPLE_LIST = [i.split("/")[-1].split(".")[0] for i in SAMPLE_LIST]

        for SAMPLE in SAMPLE_LIST:       # PD42566b_lo00_A7.txt
            INPUT_TSV = DIR + "/" + DONOR + "/" + TISSUE + "/" + SAMPLE + ".txt"
            
            kwargs = {"INPUT_TSV" : INPUT_TSV,  "MODE" : "Both",  "NUM_CLONE_TRIAL_START" : 1, "NUM_CLONE_TRIAL_END" : 5, "NUM_CLONE_TRIAL_FORCE" : 1,
                        "RANDOM_PICK":100, "AXIS_RATIO":0, "PARENT_RATIO": 0, "NUM_PARENT" : 0,  "FP_RATIO":0,  "FP_USEALL" : "False", "TRIAL_NO" : 5, "DEPTH_CUTOFF" : 10,  "MIN_CLUSTER_SIZE" : 5,  "VERBOSE" : 3, 
                        "KMEANS_CLUSTERNO" : 6, "RANDOM_SEED" : 1, "SAMPLENAME" : "", "BENCHMARK_NO" : 10, 
                        "NPVAF_DIR" : "", "MYEM_DIR" : "", "SCICLONE_DIR" : "", "PYCLONE_DIR" : "", "PYCLONEVI_DIR" : "",  "COMBINED_OUTPUT_DIR" : "",
                        "SCORING" : False, "MAKEONE_STRICT" :  2, "MAXIMUM_NUM_PARENT" : 1  }


            NUMBER_LINE = int(out ("wc -l  " + INPUT_TSV).split(" ")[0]) 
            print ("INPUT_TSV = {}\tNUMBER_LINE = {}".format (INPUT_TSV, NUMBER_LINE))
            if NUMBER_LINE < 100:       # 100줄이 안 되는 파일은 넘어간다
                continue 
            elif NUMBER_LINE < 300:
                kwargs ["RANDOM_PICK"] = 100
            else:
                kwargs ["RANDOM_PICK"] = 300

            INPUT_TSV = kwargs["INPUT_TSV"]
            INPUT_FILETYPE, NUM_BLOCK = filetype.main (INPUT_TSV)
            kwargs["NUM_BLOCK_INPUT"], kwargs["NUM_BLOCK"] = NUM_BLOCK, NUM_BLOCK
            SAMPLENAME = INPUT_TSV.split("/")[-1].split(".")[0]     
            kwargs["SAMPLENAME"] = SAMPLE      # 'C10_H9'

            kwargs["NPVAF_DIR"] = "/data/project/Alzheimer/YSscript/EM_MRS/data/npvaf/Moore_1D/" + TISSUE + "/" + DONOR + "_" + SAMPLE
            kwargs["MYEM_DIR"] = "/data/project/Alzheimer/YSscript/EM_MRS/data/CLEMENT/Moore_1D/" + TISSUE + "/" + DONOR + "_" + SAMPLE
            kwargs["SCICLONE_DIR"] = "/data/project/Alzheimer/YSscript/EM_MRS/data/sciclone/Moore_1D/" + TISSUE + "/" + DONOR + "-" + SAMPLE
            kwargs["PYCLONE_DIR"] = "/data/project/Alzheimer/YSscript/EM_MRS/data/pyclone/Moore_1D/" + TISSUE + "/" + DONOR + "-" + SAMPLE
            kwargs["PYCLONEVI_DIR"] = "/data/project/Alzheimer/YSscript/EM_MRS/data/pyclone-vi/Moore_1D/" + TISSUE  + "/" + DONOR + "-" + SAMPLE
            kwargs["QUANTUMCLONE_DIR"] = "/data/project/Alzheimer/YSscript/EM_MRS/data/quantumclone/Moore_1D/" + TISSUE  + "/" + DONOR + "-" + SAMPLE
            kwargs["COMBINED_OUTPUT_DIR"] = "/data/project/Alzheimer/YSscript/EM_MRS/data/combinedoutput/Moore_1D/" + TISSUE + "/" + DONOR + "-" + SAMPLE

            os.system ("mkdir -p " + kwargs["NPVAF_DIR"])   # 출력 디렉토리 만들기
            os.system ("mkdir -p " + kwargs["MYEM_DIR"])   # 출력 디렉토리 만들기
            os.system ("mkdir -p " + kwargs["SCICLONE_DIR"])   # 출력 디렉토리 만들기
            os.system ("mkdir -p " + kwargs["PYCLONE_DIR"])   # 출력 디렉토리 만들기
            os.system ("mkdir -p " + kwargs["PYCLONEVI_DIR"])   # 출력 디렉토리 만들기
            os.system ("mkdir -p " + kwargs["QUANTUMCLONE_DIR"])   # 출력 디렉토리 만들기
            os.system ("mkdir -p " + kwargs["COMBINED_OUTPUT_DIR"])   # 출력 디렉토리 만들기
            

            logPath = "/data/project/Alzheimer/YSscript/EM_MRS/log/Moore_1D/" + TISSUE + "/" + DONOR + "-" + SAMPLE
            os.system ("rm -rf " + logPath)
            os.system ("mkdir -p " + logPath)

            command = " ".join ( [ "qsub -pe smp 1", "-e", logPath, "-o", logPath, "-N", TISSUE + "_" + DONOR + "_" + SAMPLE, 
                                                "2.CellData_pipe1_EMhybrid.sh",  
                                                str(INPUT_TSV),  str(kwargs["MODE"]),  str(kwargs["NUM_CLONE_TRIAL_START"]),  str(kwargs["NUM_CLONE_TRIAL_END"]),  str(kwargs["NUM_CLONE_TRIAL_FORCE"]),
                                                str(kwargs["RANDOM_PICK"]), str(kwargs["AXIS_RATIO"]),  str(kwargs["PARENT_RATIO"]),  str(kwargs["NUM_PARENT"]),  str(kwargs["FP_RATIO"]),  str(kwargs["FP_USEALL"]),
                                                str(kwargs["TRIAL_NO"]), str(kwargs["DEPTH_CUTOFF"]),  str(kwargs["MIN_CLUSTER_SIZE"]),  str(kwargs["VERBOSE"]),
                                                str(kwargs["KMEANS_CLUSTERNO"]),  str(kwargs["RANDOM_SEED"]), str(kwargs["SAMPLENAME"]), str(kwargs["BENCHMARK_NO"]), 
                                                str(kwargs["NPVAF_DIR"]), str(kwargs["MYEM_DIR"]), str(kwargs["SCICLONE_DIR"]), str(kwargs["PYCLONE_DIR"]), str(kwargs["PYCLONEVI_DIR"]) , str(kwargs["QUANTUMCLONE_DIR"]),  str(kwargs["COMBINED_OUTPUT_DIR"]), 
                                                str(kwargs["SCORING"]), str(kwargs["MAKEONE_STRICT"]), str(kwargs["MAXIMUM_NUM_PARENT"])     ] )
    

            os.system (command)
            #print (command)
            k = k  + 1

    #         if k >= 5:
    #             break
    #     if k >= 5:
    #         break
    # if k >= 5:
    #     break