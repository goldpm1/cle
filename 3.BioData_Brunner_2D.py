import numpy as np
import pandas as pd
import os, subprocess, filetype, re, glob

# python3 /data/project/Alzheimer/YSscript/cle/3.BioData_Brunner_2D.py

def out(command): 
    result = subprocess.run(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True, shell=True) 
    return result.stdout.rstrip("\n")



####################################### Adrenal gland 서로 다른 층 제외 ##################################

if __name__ == "__main__":
    INPUT_TSV_DIR = "/data/project/Alzheimer/CLEMENT/01.INPUT_TSV/3.BioData/Brunner_2D/PASS"

    SCRIPT_DIR = os.path.dirname(__file__)
    print (SCRIPT_DIR, "\n")


    from natsort import natsorted
    INPUT_DIR_LIST = natsorted (  glob.glob(INPUT_TSV_DIR + "/*") ) 


    toollist = ["CLEMENT_decision", "pyclonevi", "sciclone", "quantumclone", "simpleK_elbow", "simpleK_silhouette", "simpleK_gap"]
    toollist_concise = ["CLEMENT", "pyclonevi", "sciclone", "qc", "simpleK_elb", "simpleK_sil", "simpleK_gap*"]


    INPUT_DIR_LIST_SELECTED = []
    for k in INPUT_DIR_LIST:
        if os.path.isfile (k) == True:
            if "input" in k.split("/")[-1]:
                SAMPLENAME = k.split("/")[-1].split("_input")[0]                
                INPUT_DIR_LIST_SELECTED.append (k)

    INPUT_DIR_LIST = natsorted (  INPUT_DIR_LIST_SELECTED ) 
    
    

    for INPUT_TSV in INPUT_DIR_LIST:
        kwargs = {"INPUT_TSV" : INPUT_TSV,  "MODE" : "Both",  "NUM_CLONE_TRIAL_START" : 1, "NUM_CLONE_TRIAL_END" : 5, 
                        "TRIAL_NO" : 5, "DEPTH_CUTOFF" : 10,  "KMEANS_CLUSTERNO" : 6, "MIN_CLUSTER_SIZE" : 5,  "MAKEONE_STRICT" :  3,
                        "RANDOM_PICK" : -1, "AXIS_RATIO": -1, "PARENT_RATIO": 0, "NUM_PARENT" : 0,  "FP_RATIO":0,  "FP_USEALL" : "False", 
                        "RANDOM_SEED" : 0, "SAMPLENAME" : "", "BENCHMARK_NO" : 0, 
                        "NPVAF_DIR" : "", "SIMPLE_KMEANS_DIR" : "", "CLEMENT_DIR" : "", "SCICLONE_DIR" : "", "PYCLONEVI_DIR" : "",  "COMBINED_OUTPUT_DIR" : "",
                        "SCORING" : False,  "MAXIMUM_NUM_PARENT" : 0, "VERBOSE" : 1 }
            
        INPUT_FILETYPE, NUM_BLOCK = filetype.main (INPUT_TSV)
        kwargs["NUM_BLOCK_INPUT"] = kwargs["NUM_BLOCK"] = NUM_BLOCK
        SAMPLENAME = kwargs["SAMPLENAME"] = INPUT_TSV.split("/")[-1].split("_input")[0]

        kwargs["NPVAF_DIR"] = "/data/project/Alzheimer/CLEMENT/02.npvaf/3.BioData/Brunner_2D/" + SAMPLENAME
        kwargs["COMBINED_OUTPUT_DIR"] = "/data/project/Alzheimer/CLEMENT/03.combinedoutput/3.BioData/Brunner_2D/" + SAMPLENAME
        kwargs["SIMPLE_KMEANS_DIR"] = "/data/project/Alzheimer/YSscript/cle/data/SIMPLE_KMEANS/3.BioData/Brunner_2D/" + SAMPLENAME
        kwargs["CLEMENT_DIR"] = "/data/project/Alzheimer/YSscript/cle/data/CLEMENT/3.BioData/Brunner_2D/"  + SAMPLENAME
        kwargs["PYCLONEVI_DIR"] = "/data/project/Alzheimer/YSscript/cle/data/pyclone-vi/3.BioData/Brunner_2D/" + SAMPLENAME
        kwargs["SCICLONE_DIR"] = "/data/project/Alzheimer/YSscript/cle/data/sciclone/3.BioData/Brunner_2D/" + SAMPLENAME
        kwargs["QUANTUMCLONE_DIR"] = "/data/project/Alzheimer/YSscript/cle/data/quantumclone/3.BioData/Brunner_2D/" + SAMPLENAME
        

        os.system ("mkdir -p " + kwargs["NPVAF_DIR"])   # 출력 디렉토리 만들기
        os.system ("mkdir -p " + kwargs["SIMPLE_KMEANS_DIR"])   # 출력 디렉토리 만들기
        os.system ("mkdir -p " + kwargs["CLEMENT_DIR"])   # 출력 디렉토리 만들기
        os.system ("mkdir -p " + kwargs["SCICLONE_DIR"])   # 출력 디렉토리 만들기
        os.system ("mkdir -p " + kwargs["PYCLONEVI_DIR"])   # 출력 디렉토리 만들기
        os.system ("mkdir -p " + kwargs["QUANTUMCLONE_DIR"])   # 출력 디렉토리 만들기
        os.system ("mkdir -p " + kwargs["COMBINED_OUTPUT_DIR"])   # 출력 디렉토리 만들기


        logPath = "/data/project/Alzheimer/YSscript/cle/log/3.BioData/Brunner_2D/" + SAMPLENAME
        os.system ("rm -rf " + logPath)
        os.system ("mkdir -p " + logPath)
        
        command = " ".join ( [ "qsub -pe smp 1", "-e", logPath, "-o", logPath, 
                                                "-N Brunner_2D_" + SAMPLENAME, 
                                                SCRIPT_DIR  + "/2.CellData_pipe1_CLEMENT_bm.sh",  
                                                str(SCRIPT_DIR), str(INPUT_TSV),  str(kwargs["MODE"]),  str(kwargs["NUM_CLONE_TRIAL_START"]),  str(kwargs["NUM_CLONE_TRIAL_END"]), 
                                                str(kwargs["RANDOM_PICK"]), str(kwargs["AXIS_RATIO"]),  str(kwargs["PARENT_RATIO"]),  str(kwargs["NUM_PARENT"]),  str(kwargs["FP_RATIO"]),  str(kwargs["FP_USEALL"]),
                                                str(kwargs["TRIAL_NO"]), str(kwargs["DEPTH_CUTOFF"]),  str(kwargs["MIN_CLUSTER_SIZE"]),  str(kwargs["VERBOSE"]),
                                                str(kwargs["KMEANS_CLUSTERNO"]),  str(kwargs["RANDOM_SEED"]), str(kwargs["SAMPLENAME"]), str(kwargs["BENCHMARK_NO"]), 
                                                str(kwargs["NPVAF_DIR"]), str(kwargs["SIMPLE_KMEANS_DIR"]), str(kwargs["CLEMENT_DIR"]), str(kwargs["SCICLONE_DIR"]), str(kwargs["PYCLONEVI_DIR"]) , str(kwargs["QUANTUMCLONE_DIR"]),  str(kwargs["COMBINED_OUTPUT_DIR"]), 
                                                str(kwargs["SCORING"]), str(kwargs["MAKEONE_STRICT"]), str(kwargs["MAXIMUM_NUM_PARENT"])     ] )
            
        os.system (command)
        #print (command)