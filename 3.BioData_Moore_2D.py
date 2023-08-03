import numpy as np
import pandas as pd
import os, filetype


# 2D 가능한 모든 Moore dataset을 돈다
# 총 52개의 sample
# monoclonal (colon_crypt, small_bowel_crypt_, appendix_crypt) : 10개
# biclonal : (prostate_Acinus, testis_seminiferous_tubule) : 13개
# polyclonal (oesophagus_epithelium, bladder ,ureter, skin, thyroid_follicle, bronchus) : 29개


inputdf = pd.read_csv("/data/project/Alzheimer/EM_cluster/EM_input/summary/Moore_2_sample.txt", sep = "\t")
TISSUE_set = set([])

for k in range (inputdf.shape[0]):
    DONOR, TISSUE, SAMPLE = inputdf.iloc[k]["DONOR"], inputdf.iloc[k]["TISSUE"], inputdf.iloc[k]["SAMPLE"]
    INPUT_TSV = "/".join(["/data/project/Alzheimer/EM_cluster/EM_input/Moore_2_sample", DONOR, TISSUE, SAMPLE+"_input.txt"])
    
   
    kwargs = {"INPUT_TSV" : INPUT_TSV,  "MODE" : "Both",  "NUM_CLONE_TRIAL_START" : 1, "NUM_CLONE_TRIAL_END" : 5, "NUM_CLONE_TRIAL_FORCE" : 1,
                "RANDOM_PICK":50, "AXIS_RATIO":0, "PARENT_RATIO": 0, "NUM_PARENT" : 0, "FP_RATIO": 0, "FP_USEALL" : "False", "TRIAL_NO" : 5, "DEPTH_CUTOFF" : 10, "MIN_CLUSTER_SIZE" : 3,  "VERBOSE" : 1,  
                "KMEANS_CLUSTERNO" : 6, "RANDOM_SEED" : 1, "SAMPLENAME" : "", "BENCHMARK_NO" : 10, 
                "NPVAF_DIR" : "", "CLEMENT_DIR" : "", "SCICLONE_DIR" : "", "PYCLONE_DIR" : "", "PYCLONEVI_DIR" : "",  "QUANTUMCLONE_DIR" : "", "COMBINED_OUTPUT_DIR" : "", 
                "SCORING" : False, "MAKEONE_STRICT" : 2, "MAXIMUM_NUM_PARENT" : 2 }
    if int ( inputdf.iloc[k]["SHARED"] )  < 60:
        continue
    elif int ( inputdf.iloc[k]["SHARED"] )  < 300:
        kwargs ["RANDOM_PICK"] = 100
    else:
        kwargs ["RANDOM_PICK"] = 300

    INPUT_TSV = kwargs["INPUT_TSV"]
    INPUT_FILETYPE, NUM_BLOCK = filetype.main (INPUT_TSV)
    kwargs["NUM_BLOCK_INPUT"], kwargs["NUM_BLOCK"] = NUM_BLOCK, NUM_BLOCK
    SAMPLENAME = INPUT_TSV.split("/")[-1].split(".")[0]     
    kwargs["SAMPLENAME"] = SAMPLE      # 'C10_H9'

    kwargs["NPVAF_DIR"] = "/data/project/Alzheimer/YSscript/EM_MRS/data/npvaf/Moore_2D/" + TISSUE + "/" + DONOR + "_" + SAMPLE
    kwargs["CLEMENT_DIR"] = "/data/project/Alzheimer/YSscript/EM_MRS/data/CLEMENT/Moore_2D/" + TISSUE + "/" + DONOR + "_" + SAMPLE
    kwargs["SCICLONE_DIR"] = "/data/project/Alzheimer/YSscript/EM_MRS/data/sciclone/Moore_2D/" + TISSUE + "/" + DONOR + "-" + SAMPLE
    kwargs["PYCLONE_DIR"] = "/data/project/Alzheimer/YSscript/EM_MRS/data/pyclone/Moore_2D/" + TISSUE + "/" + DONOR + "-" + SAMPLE
    kwargs["PYCLONEVI_DIR"] = "/data/project/Alzheimer/YSscript/EM_MRS/data/pyclone-vi/Moore_2D/" + TISSUE  + "/" + DONOR + "-" + SAMPLE
    kwargs["QUANTUMCLONE_DIR"] = "/data/project/Alzheimer/YSscript/EM_MRS/data/quantumclone/Moore_2D/" + TISSUE  + "/" + DONOR + "-" + SAMPLE
    kwargs["COMBINED_OUTPUT_DIR"] = "/data/project/Alzheimer/YSscript/EM_MRS/data/combinedoutput/Moore_2D/" + TISSUE + "/" + DONOR + "-" + SAMPLE
    

    print ("k = {}\tSAMPLENAME : {}\tDONOR : {}\tSAMPLE : {}\tTISSUE : {}".format (k, SAMPLENAME, DONOR, SAMPLE, TISSUE) )

    os.system ("mkdir -p " + kwargs["NPVAF_DIR"])   # 출력 디렉토리 만들기
    os.system ("mkdir -p " + kwargs["CLEMENT_DIR"])   # 출력 디렉토리 만들기
    os.system ("mkdir -p " + kwargs["SCICLONE_DIR"])   # 출력 디렉토리 만들기
    os.system ("mkdir -p " + kwargs["PYCLONE_DIR"])   # 출력 디렉토리 만들기
    os.system ("mkdir -p " + kwargs["PYCLONEVI_DIR"])   # 출력 디렉토리 만들기
    os.system ("mkdir -p " + kwargs["QUANTUMCLONE_DIR"])   # 출력 디렉토리 만들기
    os.system ("mkdir -p " + kwargs["COMBINED_OUTPUT_DIR"])   # 출력 디렉토리 만들기


    logPath = "/data/project/Alzheimer/YSscript/EM_MRS/log/Moore_2D/" + TISSUE + "/" + DONOR + "-" + SAMPLE
    os.system ("rm -rf " + logPath)
    os.system ("mkdir -p " + logPath)

    command = " ".join ( [ "qsub -pe smp 1", "-e", logPath, "-o", logPath, "-N", TISSUE + "_" + DONOR + "_" + SAMPLE, 
                                        "2.CellData_pipe1_EMhybrid.sh",  
                                        str(INPUT_TSV),  str(kwargs["MODE"]),  str(kwargs["NUM_CLONE_TRIAL_START"]),  str(kwargs["NUM_CLONE_TRIAL_END"]),  str(kwargs["NUM_CLONE_TRIAL_FORCE"]),
                                        str(kwargs["RANDOM_PICK"]), str(kwargs["AXIS_RATIO"]),  str(kwargs["PARENT_RATIO"]),  str(kwargs["NUM_PARENT"]),  str(kwargs["FP_RATIO"]),  str(kwargs["FP_USEALL"]),
                                        str(kwargs["TRIAL_NO"]), str(kwargs["DEPTH_CUTOFF"]),  str(kwargs["MIN_CLUSTER_SIZE"]),  str(kwargs["VERBOSE"]),
                                        str(kwargs["KMEANS_CLUSTERNO"]),  str(kwargs["RANDOM_SEED"]), str(kwargs["SAMPLENAME"]), str(kwargs["BENCHMARK_NO"]), 
                                        str(kwargs["NPVAF_DIR"]), str(kwargs["CLEMENT_DIR"]), str(kwargs["SCICLONE_DIR"]), str(kwargs["PYCLONE_DIR"]), str(kwargs["PYCLONEVI_DIR"]) , str(kwargs["QUANTUMCLONE_DIR"]),  str(kwargs["COMBINED_OUTPUT_DIR"]), 
                                        str(kwargs["SCORING"])  , str(kwargs["MAKEONE_STRICT"]), str(kwargs["MAXIMUM_NUM_PARENT"])  ] )
    
    
    os.system (command)
    #print (command)

    # if k >= 1:
    #     break