import filetype
import argparse
import os
import glob
import numpy as np
import pandas as pd
import random

# [1D]  10 x (BENCHMARK_NO) x (BENCHMARK_NUM_PARENT + 1)
# python3 2.CellData_master.py --DIMENSION 1D --DATA_DATE 230111  --BENCHMARK_NO 10 --FP_RATIO 0  --FP_USEALL "False" --AXIS_RATIO 0   --BENCHMARK_NUM_PARENT  1   --MAXIMUM_NUM_PARENT 1
# python3 2.CellData_master.py --DIMENSION 1D --DATA_DATE 230111  --BENCHMARK_NO 10 --FP_RATIO 0.05 --FP_USEALL "False" --AXIS_RATIO 0  --BENCHMARK_NUM_PARENT  1   --MAXIMUM_NUM_PARENT 1
# python3 2.CellData_master.py --DIMENSION 1D --DATA_DATE 230111  --BENCHMARK_NO 10 --FP_RATIO 0.1 --FP_USEALL "False" --AXIS_RATIO 0    --BENCHMARK_NUM_PARENT  1    --MAXIMUM_NUM_PARENT 1
# python3 2.CellData_master.py --DIMENSION 1D --DATA_DATE 230111  --BENCHMARK_NO 10 --FP_RATIO 0  --FP_USEALL "False" --AXIS_RATIO 0.05   --BENCHMARK_NUM_PARENT  1   --MAXIMUM_NUM_PARENT 1
# python3 2.CellData_master.py --DIMENSION 1D --DATA_DATE 230111  --BENCHMARK_NO 10 --FP_RATIO 0.05 --FP_USEALL "False" --AXIS_RATIO 0.05  --BENCHMARK_NUM_PARENT  1   --MAXIMUM_NUM_PARENT 1
# python3 2.CellData_master.py --DIMENSION 1D --DATA_DATE 230111  --BENCHMARK_NO 10 --FP_RATIO 0.1 --FP_USEALL "False" --AXIS_RATIO 0.05    --BENCHMARK_NUM_PARENT  1   --MAXIMUM_NUM_PARENT 1
# python3 2.CellData_master.py --DIMENSION 1D --DATA_DATE 230111  --BENCHMARK_NO 10 --FP_RATIO 0  --FP_USEALL "False" --AXIS_RATIO 0.1  --BENCHMARK_NUM_PARENT  1   --MAXIMUM_NUM_PARENT 1
# python3 2.CellData_master.py --DIMENSION 1D --DATA_DATE 230111  --BENCHMARK_NO 10 --FP_RATIO 0.05 --FP_USEALL "False" --AXIS_RATIO 0.1  --BENCHMARK_NUM_PARENT  1   --MAXIMUM_NUM_PARENT 1
# python3 2.CellData_master.py --DIMENSION 1D --DATA_DATE 230111  --BENCHMARK_NO 10 --FP_RATIO 0.1 --FP_USEALL "False" --AXIS_RATIO 0.1    --BENCHMARK_NUM_PARENT  1   --MAXIMUM_NUM_PARENT 1
# -MAKEONE_STRICT old --MAKEONE_STRICT sil 
# [2D] 21 x (BENCHMARK_NO + 1) x (BENCHMARK_NUM_PARENT + 1)
# python3 2.CellData_master.py --DIMENSION 2D --DATA_DATE 230111  --BENCHMARK_NO 10 --FP_RATIO 0 --FP_USEALL "False" --AXIS_RATIO 0  --BENCHMARK_NUM_PARENT  1
# python3 2.CellData_master.py --DIMENSION 2D --DATA_DATE 230111  --BENCHMARK_NO 10 --FP_RATIO 0.05 --FP_USEALL "False" --AXIS_RATIO 0  --BENCHMARK_NUM_PARENT  1
# python3 2.CellData_master.py --DIMENSION 2D --DATA_DATE 230111  --BENCHMARK_NO 10 --FP_RATIO 0.1 --FP_USEALL "False" --AXIS_RATIO 0  --BENCHMARK_NUM_PARENT  1
# python3 2.CellData_master.py --DIMENSION 2D --DATA_DATE 230111  --BENCHMARK_NO 10 --FP_RATIO 0 --FP_USEALL "False" --AXIS_RATIO 0.05  --BENCHMARK_NUM_PARENT  1
# python3 2.CellData_master.py --DIMENSION 2D --DATA_DATE 230111  --BENCHMARK_NO 10 --FP_RATIO 0.05 --FP_USEALL "False" --AXIS_RATIO 0.05  --BENCHMARK_NUM_PARENT  1
# python3 2.CellData_master.py --DIMENSION 2D --DATA_DATE 230111  --BENCHMARK_NO 10 --FP_RATIO 0.1 --FP_USEALL "False" --AXIS_RATIO 0.05  --BENCHMARK_NUM_PARENT  1
# python3 2.CellData_master.py --DIMENSION 2D --DATA_DATE 230111  --BENCHMARK_NO 10 --FP_RATIO 0 --FP_USEALL "False" --AXIS_RATIO 0.1  --BENCHMARK_NUM_PARENT  1
# python3 2.CellData_master.py --DIMENSION 2D --DATA_DATE 230111  --BENCHMARK_NO 10 --FP_RATIO 0.05 --FP_USEALL "False" --AXIS_RATIO 0.1  --BENCHMARK_NUM_PARENT  1
# python3 2.CellData_master.py --DIMENSION 2D --DATA_DATE 230111  --BENCHMARK_NO 10 --FP_RATIO 0.1  --FP_USEALL "False" --AXIS_RATIO 0.1  --BENCHMARK_NUM_PARENT  1 

# [3D] 24 x (BENCHMARK_NO + 1) x (BENCHMARK_NUM_PARENT + 1)
# python3 2.CellData_master.py --DIMENSION 3D --DATA_DATE 230111  --BENCHMARK_NO 10 --FP_RATIO 0 --FP_USEALL "False" --AXIS_RATIO 0  --BENCHMARK_NUM_PARENT  1
# python3 2.CellData_master.py --DIMENSION 3D --DATA_DATE 230111  --BENCHMARK_NO 10 --FP_RATIO 0.05 --FP_USEALL "False" --AXIS_RATIO 0  --BENCHMARK_NUM_PARENT  1
# python3 2.CellData_master.py --DIMENSION 3D --DATA_DATE 230111  --BENCHMARK_NO 10 --FP_RATIO 0.1 --FP_USEALL "False" --AXIS_RATIO 0  --BENCHMARK_NUM_PARENT  1
# python3 2.CellData_master.py --DIMENSION 3D --DATA_DATE 230111  --BENCHMARK_NO 10 --FP_RATIO 0 --FP_USEALL "False" --AXIS_RATIO 0.05  --BENCHMARK_NUM_PARENT  1
# python3 2.CellData_master.py --DIMENSION 3D --DATA_DATE 230111  --BENCHMARK_NO 10 --FP_RATIO 0.05 --FP_USEALL "False" --AXIS_RATIO 0.05  --BENCHMARK_NUM_PARENT  1
# python3 2.CellData_master.py --DIMENSION 3D --DATA_DATE 230111  --BENCHMARK_NO 10 --FP_RATIO 0.1 --FP_USEALL "False" --AXIS_RATIO 0.05  --BENCHMARK_NUM_PARENT  1
# python3 2.CellData_master.py --DIMENSION 3D --DATA_DATE 230111  --BENCHMARK_NO 10 --FP_RATIO 0 --FP_USEALL "False" --AXIS_RATIO 0.1 --BENCHMARK_NUM_PARENT  1
# python3 2.CellData_master.py --DIMENSION 3D --DATA_DATE 230111  --BENCHMARK_NO 10 --FP_RATIO 0.05 --FP_USEALL "False" --AXIS_RATIO 0.1 --BENCHMARK_NUM_PARENT  1
# python3 2.CellData_master.py --DIMENSION 3D --DATA_DATE 230111  --BENCHMARK_NO 10 --FP_RATIO 0.1 --FP_USEALL "False" --AXIS_RATIO 0.1  --BENCHMARK_NUM_PARENT  1

BENCHMARK_NUM_PARENT = 0

parser = argparse.ArgumentParser(description='The below is usage direction.')
parser.add_argument('--DIMENSION', type=str, default="2D")
parser.add_argument('--STOP_K', type=int, default=5000)
parser.add_argument('--BENCHMARK_NUM_PARENT', type=int, default=0)
parser.add_argument('--MAXIMUM_NUM_PARENT', type=int, default=3)
parser.add_argument('--DATA_DATE', type=str, default="230101")
parser.add_argument('--AXIS_RATIO', type=float, default=0.1)
parser.add_argument('--FP_RATIO', type=float, default=0.1)
parser.add_argument('--FP_USEALL', type=str, default="False")
parser.add_argument('--BENCHMARK_NO', type=int, default=30)
args = parser.parse_args()

kwargs = {}
kwargs["DIMENSION"] = args.DIMENSION
kwargs["STOP_K"] = args.STOP_K
kwargs["DATA_DATE"] = args.DATA_DATE
kwargs["FP_RATIO"] = float(args.FP_RATIO)
kwargs["FP_USEALL"] = str(args.FP_USEALL)
kwargs["AXIS_RATIO"] = float(args.AXIS_RATIO)
kwargs["BENCHMARK_NO"] = int(args.BENCHMARK_NO)
kwargs["BENCHMARK_NUM_PARENT"] = int(args.BENCHMARK_NUM_PARENT)
kwargs["MAXIMUM_NUM_PARENT"] = int(args.MAXIMUM_NUM_PARENT)


kwargs["MAKEONE_STRICT"] = 1
#kwargs["MODE"] = "Both" if kwargs["DIMENSION"] == "1D" else "Hard"
kwargs["MODE"] = "Both"
kwargs["NUM_CLONE_TRIAL_START"], kwargs["NUM_CLONE_TRIAL_END"], kwargs["NUM_CLONE_TRIAL_FORCE"] = 2, 7, 4
kwargs["RANDOM_PICK"], kwargs["PARENT_RATIO"], kwargs["TRIAL_NO"], kwargs["DEPTH_CUTOFF"], kwargs["MIN_CLUSTER_SIZE"], kwargs["KMEANS_CLUSTERNO"], kwargs["SCORING"]  \
    = 500, 0, 20, 60, 15, 8, "True"
kwargs["VERBOSE"]  = 2

if kwargs["DATA_DATE"] == "original":   # Original MRS dataset
    DIR = "/data/project/Alzheimer/EM_cluster/EM_input/MRS_" +    kwargs["DIMENSION"] + "_original"
    DIR_LIST = sorted (  glob.glob(DIR + "/*_input.txt") ) 
    kwargs["FP_USEALL"] = "False"
elif kwargs["DATA_DATE"] == "230111":      # Revision dataset (2023.01.11)
    DIR = "/data/project/Alzheimer/EM_cluster/EM_input/Revision/MRS_" +  kwargs["DIMENSION"] + "_original_dp_filter_off"
    DIR_LIST = sorted (  glob.glob(DIR + "/*_input.txt") )
    kwargs["FP_USEALL"] = "False"
    #kwargs["FP_RATIO"], kwargs["FP_USEALL"] = 0, "False"
SAMPLENAME_LIST = [i.split("/")[-1] for i in DIR_LIST]

global MRS_set
MRS_set = [['M1-2', 'M1-5', 'M1-6', 'M1-7', 'M1-8'],
           ['M2-2', 'M2-4', 'M2-6', 'M2-8', 'M2-10']]


def M1only(SAMPLENAME, DIMENSION):
    if DIMENSION == "1D":
        for i in range(len(MRS_set)):
            if (SAMPLENAME.split("_")[0] in MRS_set[i]): 
                return True
        # if 'M3' in SAMPLENAME.split("_")[0]:   # M1, M2를 다 돌려보자
        #     return False
        # return True
    elif DIMENSION == "2D":
        for i in range(len(MRS_set)):
            if (SAMPLENAME.split("_")[0] in MRS_set[i]) & (SAMPLENAME.split("_")[1] in MRS_set[i]):
                return True
    elif DIMENSION == "3D":
        for i in range(len(MRS_set)):
            if (SAMPLENAME.split("_")[0] in MRS_set[i]) & (SAMPLENAME.split("_")[1] in MRS_set[i]) & (SAMPLENAME.split("_")[2] in MRS_set[i]):
                return True

    return False


n = 0


hold_jj = []
for INPUT_TSV_k, INPUT_TSV in enumerate(DIR_LIST):
    kwargs["INPUT_TSV"] = INPUT_TSV
    INPUT_FILETYPE, NUM_BLOCK = filetype.main(INPUT_TSV)
    kwargs["NUM_BLOCK_INPUT"], kwargs["NUM_BLOCK"] = NUM_BLOCK, NUM_BLOCK
    SAMPLENAME = INPUT_TSV.split("/")[-1].split(".")[0]     # 'M1-5_M1-8_input'
    kwargs["SAMPLENAME"] = SAMPLENAME

    if M1only(SAMPLENAME, kwargs["DIMENSION"]) == True:
        for NUM_PARENT in range(0, kwargs["BENCHMARK_NUM_PARENT"] + 1)  :
            kwargs["NUM_PARENT"] = NUM_PARENT
            print("# {}:  {}, NUM_PARENT = {}".format(INPUT_TSV_k, SAMPLENAME, NUM_PARENT))

            # 1.
            hold_j = []
            for i in range(0, kwargs["BENCHMARK_NO"]):
                kwargs["NPVAF_DIR"] = "/data/project/Alzheimer/YSscript/EM_MRS/data/npvaf/MRS_" + kwargs["DIMENSION"] + "/" + str( kwargs["RANDOM_PICK"]) + "_" + str(kwargs["NUM_PARENT"]) + "_" + str(kwargs["FP_RATIO"]) + "_" + str(kwargs["AXIS_RATIO"]) + "/" +  SAMPLENAME + "/" + str(i)
                kwargs["CLEMENT_DIR"] = "/data/project/Alzheimer/YSscript/EM_MRS/data/CLEMENT/MRS_" + kwargs["DIMENSION"] + "/" +  str(kwargs["RANDOM_PICK"]) + "_" + str(kwargs["NUM_PARENT"]) + "_" + str(kwargs["FP_RATIO"]) + "_" + str(kwargs["AXIS_RATIO"]) + "/" + SAMPLENAME  + "/" + str(i)
                kwargs["SCICLONE_DIR"] = "/data/project/Alzheimer/YSscript/EM_MRS/data/sciclone/MRS_" + kwargs["DIMENSION"] + "/" + str( kwargs["MAKEONE_STRICT"]) + "/" + str(kwargs["RANDOM_PICK"]) + "_" + str(kwargs["NUM_PARENT"]) + "_" + str(kwargs["FP_RATIO"]) + "_" + str(kwargs["AXIS_RATIO"]) + "/" + SAMPLENAME + "/" + str(i)
                kwargs["PYCLONE_DIR"] = "/data/project/Alzheimer/YSscript/EM_MRS/data/pyclone/MRS_" + kwargs["DIMENSION"] + "/" +  str(kwargs["RANDOM_PICK"]) + "_" + str(kwargs["NUM_PARENT"]) + "_" + str(kwargs["FP_RATIO"]) + "_" + str(kwargs["AXIS_RATIO"]) + "/" + SAMPLENAME + "/" + str(i)
                kwargs["PYCLONEVI_DIR"] = "/data/project/Alzheimer/YSscript/EM_MRS/data/pyclone-vi/MRS_" + kwargs["DIMENSION"] + "/"  + str(kwargs["RANDOM_PICK"]) + "_" + str(kwargs["NUM_PARENT"]) + "_" + str(kwargs["FP_RATIO"]) + "_" + str(kwargs["AXIS_RATIO"]) + "/" + SAMPLENAME +   "/" + str(i)
                kwargs["QUANTUMCLONE_DIR"] = "/data/project/Alzheimer/YSscript/EM_MRS/data/quantumclone/MRS_" + kwargs["DIMENSION"] + "/"  + str(kwargs["RANDOM_PICK"]) + "_" + str(kwargs["NUM_PARENT"]) + "_" + str(kwargs["FP_RATIO"]) + "_" + str(kwargs["AXIS_RATIO"]) + "/" + SAMPLENAME +  "/" + str(i)
                kwargs["COMBINED_OUTPUT_DIR"] = "/data/project/Alzheimer/YSscript/EM_MRS/data/combinedoutput/MRS_" + kwargs["DIMENSION"] + "/" + str(kwargs["RANDOM_PICK"]) + "_" + str(kwargs["NUM_PARENT"]) + "_" + str(kwargs["FP_RATIO"]) + "_" + str(kwargs["AXIS_RATIO"]) + "/" + SAMPLENAME +  "/" + str(i)

                for DIR in [ kwargs["NPVAF_DIR"], kwargs["CLEMENT_DIR"], kwargs["SCICLONE_DIR"], kwargs["PYCLONE_DIR"], kwargs["PYCLONEVI_DIR"], kwargs["QUANTUMCLONE_DIR"], kwargs["COMBINED_OUTPUT_DIR"]  ]:
                    if os.path.isdir( DIR   ) == False:
                        os.system("mkdir -p " + DIR  )

                logPath = "/data/project/Alzheimer/YSscript/EM_MRS/log/MRS_" + kwargs["DIMENSION"] + "/"  + str(kwargs["RANDOM_PICK"]) + "_" + str(kwargs["NUM_PARENT"]) + "_" + str(kwargs["FP_RATIO"]) + "_" + str(kwargs["AXIS_RATIO"]) + "/" + SAMPLENAME  + "/" + str(i)
                os.system("rm -rf " + logPath)
                os.system("mkdir -p " + logPath)

                                    #"-q ", COMPUTE_RANDOM, 
                #COMPUTE_RANDOM = "cpu.q@compute" + str( random.randint (1,14) ).zfill(2)
                hold_j.append( "_" + kwargs["DIMENSION"] + "_" + str(kwargs["RANDOM_PICK"]) + "_" + str( kwargs["NUM_PARENT"]) + "_" + str(kwargs["FP_RATIO"]) + "_" + str(kwargs["AXIS_RATIO"]) + "_" + SAMPLENAME + "_" + str(i) )
                command = " ".join(["qsub -pe smp 1", "-e", logPath, "-o", logPath, "-N", "_" + kwargs["DIMENSION"] + "_"  + str(kwargs["RANDOM_PICK"]) + "_" + str( kwargs["NUM_PARENT"]) + "_" + str(kwargs["FP_RATIO"]) + "_" + str(kwargs["AXIS_RATIO"]) + "_" + SAMPLENAME + "_" + str(i),
                                    "2.CellData_pipe1_EMhybrid.sh",
                                    str(INPUT_TSV),  str(kwargs["MODE"]),  str(kwargs["NUM_CLONE_TRIAL_START"]),  str( kwargs["NUM_CLONE_TRIAL_END"]),  str(kwargs["NUM_CLONE_TRIAL_FORCE"]),
                                    str(kwargs["RANDOM_PICK"]), str(kwargs["AXIS_RATIO"]),  str(kwargs["PARENT_RATIO"]),  str( kwargs["NUM_PARENT"]),  str(kwargs["FP_RATIO"]),  str(kwargs["FP_USEALL"]),
                                    str(kwargs["TRIAL_NO"]), str(kwargs["DEPTH_CUTOFF"]),  str( kwargs["MIN_CLUSTER_SIZE"]),  str(kwargs["VERBOSE"]),
                                    str(kwargs["KMEANS_CLUSTERNO"]),  str(i), str( kwargs["SAMPLENAME"]), str(kwargs["BENCHMARK_NO"]),
                                    str(kwargs["NPVAF_DIR"]), str(kwargs["CLEMENT_DIR"]), str(kwargs["SCICLONE_DIR"]), str(kwargs["PYCLONE_DIR"]), str( kwargs["PYCLONEVI_DIR"]), str(kwargs["QUANTUMCLONE_DIR"]), str(kwargs["COMBINED_OUTPUT_DIR"]),
                                    str(kwargs["SCORING"]), str(kwargs["MAKEONE_STRICT"]), str(kwargs["MAXIMUM_NUM_PARENT"])])
                
                #print (command)
                os.system(command)
                n = n+1

            # 2.
            logPath = "/data/project/Alzheimer/YSscript/EM_MRS/log/MRS_" + kwargs["DIMENSION"] +  "/" + str(kwargs["RANDOM_PICK"]) + "_" + str(kwargs["NUM_PARENT"]) + "_" + str(kwargs["FP_RATIO"]) + "_" + str(kwargs["AXIS_RATIO"]) + "/" + SAMPLENAME + "/visualization"
            os.system("rm -rf " + logPath)
            os.system("mkdir -p " + logPath)

            #print (",".join (hold_j))

            command = " ".join(["qsub -pe smp 1", "-e", logPath, "-o", logPath, "-N", "_" + kwargs["DIMENSION"] + "_" + str(kwargs["RANDOM_PICK"]) + "_" + str(kwargs["NUM_PARENT"]) + "_" + str(kwargs["FP_RATIO"]) + "_" + str(kwargs["AXIS_RATIO"]) + "_" + SAMPLENAME +  "_visualization",
                                "-hold_jid",  str(",".join(hold_j)), 
                                "2.CellData_pipe2_benchmark.sh",
                                str(INPUT_TSV),  str(kwargs["MODE"]),  str(kwargs["NUM_CLONE_TRIAL_START"]),  str( kwargs["NUM_CLONE_TRIAL_END"]),  str(kwargs["NUM_CLONE_TRIAL_FORCE"]),
                                str(kwargs["RANDOM_PICK"]), str(kwargs["AXIS_RATIO"]),  str(kwargs["PARENT_RATIO"]),  str(kwargs["NUM_PARENT"]),  str(kwargs["FP_RATIO"]),  str(kwargs["FP_USEALL"]),
                                str(kwargs["TRIAL_NO"]), str(kwargs["DEPTH_CUTOFF"]),  str( kwargs["MIN_CLUSTER_SIZE"]),  str(kwargs["VERBOSE"]),
                                str(kwargs["KMEANS_CLUSTERNO"]),  str(i), str(kwargs["SAMPLENAME"]), str(kwargs["BENCHMARK_NO"]),
                                str(kwargs["NPVAF_DIR"]), str(kwargs["CLEMENT_DIR"]), str(kwargs["SCICLONE_DIR"]), str(kwargs["PYCLONE_DIR"]), str(kwargs["PYCLONEVI_DIR"]), str(kwargs["QUANTUMCLONE_DIR"]), str(kwargs["COMBINED_OUTPUT_DIR"]),
                                str(kwargs["SCORING"]), str(kwargs["DIMENSION"]), str(kwargs["MAKEONE_STRICT"])])

            #print (command)
            #os.system(command)
            n = n+1
            hold_jj.append(  "_" + kwargs["DIMENSION"] + "_" + str(kwargs["RANDOM_PICK"]) + "_" + str(kwargs["NUM_PARENT"]) + "_" + str(kwargs["FP_RATIO"]) + "_" + str(kwargs["AXIS_RATIO"]) + "_" + SAMPLENAME +  "_visualization"  )
            


#3. 다 끝나야 최종 그림을 그림

for NUM_PARENT in range(0, kwargs["BENCHMARK_NUM_PARENT"] + 1)  :   
    kwargs["NUM_PARENT"] = NUM_PARENT
    INPUT_DIR="/data/project/Alzheimer/YSscript/EM_MRS/data/combinedoutput/MRS_" + kwargs["DIMENSION"] + "/" + str(kwargs["RANDOM_PICK"]) + "_" + str(kwargs["NUM_PARENT"]) + "_" + str(kwargs["FP_RATIO"]) + "_" + str(kwargs["AXIS_RATIO"])
    os.system ("rm -rf " +  INPUT_DIR + "/" + "benchmark.jpg" )
    os.system ("rm -rf " +  INPUT_DIR + "/" + str(kwargs["RANDOM_PICK"]) + "_" + str(kwargs["NUM_PARENT"]) + "_" + str(kwargs["FP_RATIO"]) + "_" + str(kwargs["AXIS_RATIO"]) + "benchmarktotal.jpg" )
    OUTPUT_FILENAME = INPUT_DIR + "/" + str(kwargs["RANDOM_PICK"]) + "_" + str(kwargs["NUM_PARENT"]) + "_" + str(kwargs["FP_RATIO"]) + "_" + str(kwargs["AXIS_RATIO"]) + "_total.jpg"
    logPath = "/data/project/Alzheimer/YSscript/EM_MRS/log/MRS_" + kwargs["DIMENSION"] +  "/" + str(kwargs["RANDOM_PICK"]) + "_" + str(kwargs["NUM_PARENT"]) + "_" + str(kwargs["FP_RATIO"]) + "_" + str(kwargs["AXIS_RATIO"]) + "/visualizationtotal"
    os.system("rm -rf " + logPath)
    os.system("mkdir -p " + logPath)
    
    command = " ".join(  ["qsub -pe smp 1", "-e", logPath, "-o", logPath, "-N", "_" + kwargs["DIMENSION"] + "_" + str(kwargs["RANDOM_PICK"]) + "_" + str(kwargs["NUM_PARENT"]) + "_" + str(kwargs["FP_RATIO"]) + "_" + str(kwargs["AXIS_RATIO"]) + "_visualizationtotal",
                        "-hold_jid",  str(",".join(hold_jj)), 
                        "2.CellData_pipe3_benchmark.sh",
                        str(INPUT_DIR),  str(OUTPUT_FILENAME),  str(kwargs["BENCHMARK_NO"])]  )
    
    print (command)
    os.system(command)
    n = n+1


print ("Total job = {}".format(n))