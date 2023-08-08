import numpy as np
import pandas as pd
import os, subprocess, filetype, re

# python3 /data/project/Alzheimer/YSscript/cle/3.BioData_Moore_2D.py

# 2D 가능한 모든 Moore dataset을 돈다
# 총 52개의 sample
# monoclonal (colon_crypt, small_bowel_crypt_, appendix_crypt) : 10개
# biclonal : (prostate_Acinus, testis_seminiferous_tubule) : 13개
# polyclonal (oesophagus_epithelium, bladder ,ureter, skin, thyroid_follicle, bronchus) : 29개

def out(command): 
    result = subprocess.run(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True, shell=True) 
    return result.stdout.rstrip("\n")



####################################### Adrenal gland 서로 다른 층 제외 ##################################

if __name__ == "__main__":
    SCRIPT_DIR = os.path.dirname(__file__)
    print (SCRIPT_DIR, "\n")

    DIR = "/data/project/Alzheimer/CLEMENT/01.INPUT_TSV/3.BioData/Moore_2D/2.all_woMosaic"        # AG : 3.woMosaic_ver2

    inputdf = pd.read_csv("/data/project/Alzheimer/EM_cluster/EM_input/summary/Moore_2_sample.txt", sep = "\t")
    TISSUE_set = set([])


    for k in range (inputdf.shape[0]):
        DONOR, TISSUE, SAMPLE = inputdf.iloc[k]["DONOR"], inputdf.iloc[k]["TISSUE"], inputdf.iloc[k]["SAMPLE"]
        INPUT_TSV = "/".join(["/data/project/Alzheimer/EM_cluster/EM_input/Moore_2_sample", DONOR, TISSUE, SAMPLE+"_input.txt"])

        if TISSUE not in ["adrenal_gland_zona_glomerulosa", "visceral_fat", "bronchus_epithelium"]:
            continue
        

        kwargs = {"INPUT_TSV" : INPUT_TSV,  "MODE" : "Both",  "NUM_CLONE_TRIAL_START" : 1, "NUM_CLONE_TRIAL_END" : 5, 
                        "TRIAL_NO" : 5, "DEPTH_CUTOFF" : 10,  "KMEANS_CLUSTERNO" : 6, "MIN_CLUSTER_SIZE" : 10,  "MAKEONE_STRICT" :  2,
                        "RANDOM_PICK" : 0, "AXIS_RATIO":0.4, "PARENT_RATIO": 0, "NUM_PARENT" : 0,  "FP_RATIO":0,  "FP_USEALL" : "False", 
                        "RANDOM_SEED" : 0, "SAMPLENAME" : "", "BENCHMARK_NO" : 0, 
                        "NPVAF_DIR" : "", "SIMPLE_KMEANS_DIR" : "", "CLEMENT_DIR" : "", "SCICLONE_DIR" : "", "PYCLONEVI_DIR" : "",  "COMBINED_OUTPUT_DIR" : "",
                        "SCORING" : False,  "MAXIMUM_NUM_PARENT" : 1, "VERBOSE" : 1 }
        
        if int ( inputdf.iloc[k]["SHARED"] )  < 60:
            continue
        elif int ( inputdf.iloc[k]["TOTAL"] )  < 350:
            continue
        else:
            kwargs ["RANDOM_PICK"] = 300

        INPUT_TSV = kwargs["INPUT_TSV"]
        INPUT_FILETYPE, NUM_BLOCK = filetype.main (INPUT_TSV)
        kwargs["NUM_BLOCK_INPUT"] = kwargs["NUM_BLOCK"] = NUM_BLOCK
        SAMPLENAME = kwargs["SAMPLENAME"] = re.split(r'[_ .]', INPUT_TSV.split("/")[-1])[0]


        kwargs["NPVAF_DIR"] = "/data/project/Alzheimer/CLEMENT/02.npvaf/3.BioData/Moore_2D/" + TISSUE + "/" + DONOR + "-" + SAMPLENAME
        kwargs["COMBINED_OUTPUT_DIR"] = "/data/project/Alzheimer/CLEMENT/03.combinedoutput/3.BioData/Moore_2D/" + TISSUE + "/" + DONOR + "-" + SAMPLENAME
        kwargs["SIMPLE_KMEANS_DIR"] = "/data/project/Alzheimer/YSscript/cle/data/SIMPLE_KMEANS/3.BioData/Moore_2D/" + TISSUE + "/" + DONOR + "-" + SAMPLENAME
        kwargs["CLEMENT_DIR"] = "/data/project/Alzheimer/YSscript/cle/data/CLEMENT/3.BioData/Moore_2D/" + TISSUE + "/" + DONOR + "-" + SAMPLENAME
        kwargs["PYCLONEVI_DIR"] = "/data/project/Alzheimer/YSscript/cle/data/pyclone-vi/3.BioData/Moore_2D/" + TISSUE  + "/" + DONOR + "-" + SAMPLENAME
        kwargs["SCICLONE_DIR"] = "/data/project/Alzheimer/YSscript/cle/data/sciclone/3.BioData/Moore_2D/" + TISSUE + "/" + DONOR + "-" + SAMPLENAME
        kwargs["QUANTUMCLONE_DIR"] = "/data/project/Alzheimer/YSscript/cle/data/quantumclone/3.BioData/Moore_2D/" + TISSUE  + "/" + DONOR + "-" + SAMPLENAME
        

        print ("k = {}\t{}/{}-{} (SAMPLE = {})\tTOTAL : {}\tSHARED : {}".format (k, TISSUE, DONOR, SAMPLENAME, SAMPLE, int ( inputdf.iloc[k]["TOTAL"] ), int ( inputdf.iloc[k]["SHARED"] )) )

        os.system ("mkdir -p " + kwargs["NPVAF_DIR"])   # 출력 디렉토리 만들기
        os.system ("mkdir -p " + kwargs["SIMPLE_KMEANS_DIR"])   # 출력 디렉토리 만들기
        os.system ("mkdir -p " + kwargs["CLEMENT_DIR"])   # 출력 디렉토리 만들기
        os.system ("mkdir -p " + kwargs["SCICLONE_DIR"])   # 출력 디렉토리 만들기
        os.system ("mkdir -p " + kwargs["PYCLONEVI_DIR"])   # 출력 디렉토리 만들기
        os.system ("mkdir -p " + kwargs["QUANTUMCLONE_DIR"])   # 출력 디렉토리 만들기
        os.system ("mkdir -p " + kwargs["COMBINED_OUTPUT_DIR"])   # 출력 디렉토리 만들기


        logPath = "/data/project/Alzheimer/YSscript/cle/log/3.BioData/Moore_2D/" + TISSUE + "/" + DONOR + "-" + SAMPLE
        os.system ("rm -rf " + logPath)
        os.system ("mkdir -p " + logPath)
        hold_j = TISSUE + "_" + DONOR + "_" + SAMPLE
        command = " ".join ( [ "qsub -pe smp 1", "-e", logPath, "-o", logPath, 
                                            "-N", TISSUE + "_" + DONOR + "_" + SAMPLE, 
                                            SCRIPT_DIR  + "/2.CellData_pipe1_CLEMENT_bm.sh",  
                                            str(SCRIPT_DIR), str(INPUT_TSV),  str(kwargs["MODE"]),  str(kwargs["NUM_CLONE_TRIAL_START"]),  str(kwargs["NUM_CLONE_TRIAL_END"]), 
                                            str(kwargs["RANDOM_PICK"]), str(kwargs["AXIS_RATIO"]),  str(kwargs["PARENT_RATIO"]),  str(kwargs["NUM_PARENT"]),  str(kwargs["FP_RATIO"]),  str(kwargs["FP_USEALL"]),
                                            str(kwargs["TRIAL_NO"]), str(kwargs["DEPTH_CUTOFF"]),  str(kwargs["MIN_CLUSTER_SIZE"]),  str(kwargs["VERBOSE"]),
                                            str(kwargs["KMEANS_CLUSTERNO"]),  str(kwargs["RANDOM_SEED"]), str(kwargs["SAMPLENAME"]), str(kwargs["BENCHMARK_NO"]), 
                                            str(kwargs["NPVAF_DIR"]), str(kwargs["SIMPLE_KMEANS_DIR"]), str(kwargs["CLEMENT_DIR"]), str(kwargs["SCICLONE_DIR"]), str(kwargs["PYCLONEVI_DIR"]) , str(kwargs["QUANTUMCLONE_DIR"]),  str(kwargs["COMBINED_OUTPUT_DIR"]), 
                                            str(kwargs["SCORING"]), str(kwargs["MAKEONE_STRICT"]), str(kwargs["MAXIMUM_NUM_PARENT"])     ] )
        
        os.system (command)
        #print (command)

        # if k >= 1:
        #     break