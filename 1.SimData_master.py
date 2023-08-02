import argparse
import numpy as np
import pandas as pd
from random import *
import os

# [FP 없는 경우]  python3 1.SimData_master.py --NUM_CLONE 4 --NUM_BLOCK 2 --BENCHMARK_START 0 --BENCHMARK_END 99  --SIMDATA decoy
# python3 1.SimData_master.py --NUM_CLONE 2 --NUM_BLOCK 1 --BENCHMARK_START 0 --BENCHMARK_END 99 --FP_RATIO 0.1  --SIMDATA decoy
# python3 1.SimData_master.py --NUM_CLONE 2 --NUM_BLOCK 2 --BENCHMARK_START 0 --BENCHMARK_END 99 --FP_RATIO 0.1  --SIMDATA lump
# python3 1.SimData_master.py --NUM_CLONE 2 --NUM_BLOCK 3 --BENCHMARK_START 0 --BENCHMARK_END 99 --FP_RATIO 0.1  --SIMDATA decoy

if __name__ == "__main__":
    kwargs = {}

    parser = argparse.ArgumentParser(description='The below is usage direction.')
    parser.add_argument('--NUM_CLONE', type=int, default=4)
    parser.add_argument('--NUM_BLOCK', type=int, default=3)
    parser.add_argument('--NUM_MUTATION', type=int, default=500)
    parser.add_argument('--FP_RATIO', type = float, default = 0)
    parser.add_argument('--LOWVAF_RATIO', type = float, default = 0)
    parser.add_argument('--DEPTH_MEAN', type=int, default=100)
    parser.add_argument('--DEPTH_SD', type=int, default=8)
    parser.add_argument('--BENCHMARK_START', type=int, default=0)
    parser.add_argument('--BENCHMARK_END', type=int, default=99)
    parser.add_argument('--SIMDATA', type=str, default="decoy")

    args = parser.parse_args()

    global NUM_CLONE, NUM_BLOCK, NUM_MUTATION
    kwargs["NUM_CLONE"], kwargs["NUM_BLOCK"], kwargs["NUM_MUTATION"] = args.NUM_CLONE, args.NUM_BLOCK, args.NUM_MUTATION
    kwargs["FP_RATIO"] = args.FP_RATIO
    FP_EXISTANCE = "noFP" if kwargs["FP_RATIO"] == 0 else "FP"
    kwargs["LOWVAF_RATIO"] = args.LOWVAF_RATIO
    NUM_CLONE, NUM_BLOCK, NUM_MUTATION = args.NUM_CLONE, args.NUM_BLOCK, args.NUM_MUTATION
    kwargs["DEPTH_MEAN"] = int(args.DEPTH_MEAN)
    kwargs["DEPTH_SD"] = int(args.DEPTH_SD)
    kwargs["BENCHMARK_START"] = int(args.BENCHMARK_START)
    kwargs["BENCHMARK_END"] = int(args.BENCHMARK_END)
    kwargs["SIMDATA"] = str(args.SIMDATA)

    hold_j = []
    print("NUM_CLONE = {}\tNUM_BLOCK = {}".format(NUM_CLONE, NUM_BLOCK))
    
    for ii in range ( kwargs["BENCHMARK_START"],  kwargs["BENCHMARK_END"] + 1 ):
        kwargs["INPUT_TSV"] = "/data/project/Alzheimer/EM_cluster/EM_input/simulation_" + str(NUM_BLOCK) + "D/" + kwargs ["SIMDATA"] + "/" + str(FP_EXISTANCE) + "/clone_" + str(NUM_CLONE) + "/" + str(NUM_BLOCK) + "D_clone" + str(NUM_CLONE) + "_" + str(ii) + ".txt"
        kwargs["NPVAF_DIR"] = "/data/project/Alzheimer/YSscript/EM_MRS/data/npvaf/simulation_" + str(NUM_BLOCK) + "D/" +  kwargs ["SIMDATA"] + "/" + str(FP_EXISTANCE) + "/clone_" + str(NUM_CLONE) + "/" + str(ii)
        kwargs["CLEMENT_DIR"] = "/data/project/Alzheimer/YSscript/EM_MRS/data/CLEMENT/simulation_" + str(NUM_BLOCK) + "D/" +  kwargs ["SIMDATA"] + "/" + str(FP_EXISTANCE) + "/clone_" + str(NUM_CLONE) + "/" + str(ii)
        kwargs["SCICLONE_DIR"] = "/data/project/Alzheimer/YSscript/EM_MRS/data/sciclone/simulation_" + str(NUM_BLOCK) + "D/" +  kwargs ["SIMDATA"] + "/" + str(FP_EXISTANCE) + "/clone_" + str(NUM_CLONE) + "/" + str(ii)
        kwargs["PYCLONEVI_DIR"] = "/data/project/Alzheimer/YSscript/EM_MRS/data/pyclone-vi/simulation_" + str(NUM_BLOCK) + "D/" +  kwargs ["SIMDATA"] + "/" + str(FP_EXISTANCE) + "/clone_" + str(NUM_CLONE) + "/" + str(ii)
        kwargs["QUANTUMCLONE_DIR"] = "/data/project/Alzheimer/YSscript/EM_MRS/data/quantumclone/simulation_" + str(NUM_BLOCK) + "D/" +  kwargs ["SIMDATA"] + "/" + str(FP_EXISTANCE) + "/clone_" + str(NUM_CLONE) + "/" + str(ii)
        kwargs["COMBINED_OUTPUT_DIR"] = "/data/project/Alzheimer/YSscript/EM_MRS/data/combinedoutput/simulation_" + str(NUM_BLOCK) + "D/" +  kwargs ["SIMDATA"] + "/" + str(FP_EXISTANCE) + "/clone_" + str(NUM_CLONE) + "/" + str(ii)

        # print("\n\t{}번째 trial : ".format(ii))

        # print ( "INPUT_TSV = {}\nNPVAF_DIR = {}\nCLEMENT_DIR = {}\nSCICLONE_DIR = {}\nPYCLOENVI_DIR = {}\nCOMBINED_OUTPUT_DIR = {}"
        #                         .format ( kwargs["INPUT_TSV"], kwargs["NPVAF_DIR"], kwargs["CLEMENT_DIR"], kwargs["SCICLONE_DIR"], kwargs["PYCLONEVI_DIR"], kwargs["COMBINED_OUTPUT_DIR"] )  )

        # 1. Simulation dataset 생성 (decoy or sparse)
        command1 = " ".join(["python3 1.SimData_pipe0_preparation.py",
                             "--NUM_CLONE", str(NUM_CLONE),  "--NUM_BLOCK", str(NUM_BLOCK), "--NUM_MUTATION", str(NUM_MUTATION), "--FP_RATIO", str(kwargs["FP_RATIO"]),   "--LOWVAF_RATIO", str(kwargs["LOWVAF_RATIO"]), 
                             "--DEPTH_MEAN", str(kwargs["DEPTH_MEAN"]), "--DEPTH_SD", str(kwargs["DEPTH_SD"]),
                             "--INPUT_TSV", kwargs["INPUT_TSV"], "--NPVAF_DIR", kwargs["NPVAF_DIR"],
                             "--BENCHMARK_I", str(ii),
                             "--SimData", kwargs["SIMDATA"]
                             ])
        #print (command1)
        os.system(command1)

        # 2. EM 돌리기 (2-7, Hard)
        logPath = "/data/project/Alzheimer/YSscript/EM_MRS/log/simulation_" + str(NUM_BLOCK) + "D/" +  kwargs ["SIMDATA"] + "/" + str(FP_EXISTANCE) + "/clone_" + str(NUM_CLONE) + "/" + str(ii)
        os.system("rm -rf " + logPath)
        os.system("mkdir -p " + logPath)

        hold_j.append("simulation_" + str(NUM_BLOCK) + "D_clone_" + str(NUM_CLONE) + "_" + kwargs["SIMDATA"] + "_" + str(ii))
        command2 = " ".join(["qsub -pe smp 1 -e", logPath, "-o", logPath, "-N simulation_" + str(NUM_BLOCK) + "D_clone_" + str(NUM_CLONE) + "_" + kwargs["SIMDATA"] + "_" + str(ii),
                             "1.SimData_pipe1_EMhybrid.sh",
                             "--INPUT_TSV", str(kwargs["INPUT_TSV"]),
                             "--NPVAF_DIR", str(kwargs["NPVAF_DIR"]),
                             "--CLEMENT_DIR", str(kwargs["CLEMENT_DIR"]),
                             "--SCICLONE_DIR", str(kwargs["SCICLONE_DIR"]),
                             "--PYCLONEVI_DIR", str(kwargs["PYCLONEVI_DIR"]),
                             "--QUANTUMCLONE_DIR", str(kwargs["QUANTUMCLONE_DIR"]),
                            "--COMBINED_OUTPUT_DIR", str(kwargs["COMBINED_OUTPUT_DIR"]),
                            "--NUM_CLONE_TRIAL_START", str(2), "--NUM_CLONE_TRIAL_END", str(7),
                            "--FP_RATIO", str(kwargs["FP_RATIO"]),
                            "--DEPTH_CUTOFF", str(10),  "--VERBOSE", str(1), "--TRIAL_NO", str(8), "--RANDOM_SEED", str(ii), "--SCORING True", "--MODE Both"
        ])
        #print(command2)
        print ("python3 EMhybrid.py --INPUT_TSV {} --NPVAF_DIR {} --CLEMENT_DIR {} --SCICLONE_DIR {} --PYCLONEVI_DIR {} --QUANTUMCLONE_DIR {} --COMBINED_OUTPUT_DIR {} --NUM_CLONE_TRIAL_START {} --NUM_CLONE_TRIAL_END {} --FP_RATIO {} --DEPTH_CUTOFF {} --VERBOSE {} --TRIAL_NO {} --RANDOM_SEED {}  --SCORING {} --MODE Both".format(
                str(kwargs["INPUT_TSV"]), str(kwargs["NPVAF_DIR"]), str(kwargs["CLEMENT_DIR"]), str(kwargs["SCICLONE_DIR"]), str(kwargs["PYCLONEVI_DIR"]), str(kwargs["QUANTUMCLONE_DIR"]), str(kwargs["COMBINED_OUTPUT_DIR"]), 2, 7, str(kwargs["FP_RATIO"]), 10, 1, 8, str(ii), "True"  ))
        #os.system(command2)

    # 채점하기
    kwargs["COMBINED_OUTPUT_DIR"] = "/data/project/Alzheimer/YSscript/EM_MRS/data/combinedoutput/simulation_" +  str(NUM_BLOCK) + "D/" +  kwargs ["SIMDATA"] + "/" + str(FP_EXISTANCE) + "/clone_" + str(NUM_CLONE)
    kwargs["SAMPLENAME"] = "simulation_" +  str(NUM_BLOCK) + "D_clone_" + str(NUM_CLONE)
    kwargs["OUTPUT_TTEST"] = kwargs["COMBINED_OUTPUT_DIR"] + "/ttest.txt"
    kwargs["OUTPUT_JPG"] = kwargs["COMBINED_OUTPUT_DIR"] + "/benchmark.jpg"
    logPath = "/data/project/Alzheimer/YSscript/EM_MRS/log/simulation_" +  str(NUM_BLOCK) + "D/" + kwargs ["SIMDATA"] + "/" + str(FP_EXISTANCE) + "/clone_" + str(NUM_CLONE) + "/visualization"
    
    os.system("rm -rf " + logPath)
    os.system("mkdir -p " + logPath)

    command3 = " ".join(["qsub -pe smp 1 -e", logPath, "-o", logPath, "-N visualization",  "-hold_jid",  str(",".join(hold_j)),
                         "1.SimData_pipe2_benchmark.sh",
                         "--COMBINED_OUTPUT_DIR", str(kwargs["COMBINED_OUTPUT_DIR"]),
                        "--SAMPLENAME", str(kwargs["SAMPLENAME"]),
                        "--BENCHMARK_START", str(kwargs["BENCHMARK_START"]),
                        "--BENCHMARK_END", str(kwargs["BENCHMARK_END"]),
                        "--OUTPUT_TTEST", str(kwargs["OUTPUT_TTEST"]),
                        "--OUTPUT_JPG", str(kwargs["OUTPUT_JPG"])
                        ])
    ##print (command3)
    #os.system(command3)



# python3 1.SimData_master.py --NUM_BLOCK 2 --NUM_CLONE 4 --BENCHMARK_NO 1

# python3 EMhybrid.py --INPUT_TSV /data/project/Alzheimer/EM_cluster/EM_input/simulation_1D/clone_4/1D_clone4_0.txt --NPVAF_DIR /data/project/Alzheimer/YSscript/EM_MRS/data/npvaf/simulation_1D/clone_4/0 --CLEMENT_DIR /data/project/Alzheimer/YSscript/EM_MRS/data/CLEMENT/simulation_1D/clone_4/0 --SCICLONE_DIR /data/project/Alzheimer/YSscript/EM_MRS/data/sciclone/simulation_1D/clone_4/0 --PYCLONEVI_DIR /data/project/Alzheimer/YSscript/EM_MRS/data/pyclone-vi/simulation_1D/clone_4/0 --COMBINED_OUTPUT_DIR /data/project/Alzheimer/YSscript/EM_MRS/data/combinedoutput/simulation_1D/clone_4/0 --NUM_CLONE_TRIAL_START 2 --NUM_CLONE_TRIAL_END 5 --DEPTH_CUTOFF 10 --VERBOSE 1 --TRIAL_NO 3 --RANDOM_SEED 0 --SCORING True --MODE Hard


# bash 1.SimData_pipe1_EMhybrid.sh --INPUT_TSV /data/project/Alzheimer/EM_cluster/EM_input/simulation_2D/clone_4/2D_clone4_5.txt \
# --NPVAF_DIR /data/project/Alzheimer/YSscript/EM_MRS/data/npvaf/simulation_2D/clone_4/5 
# --CLEMENT_DIR /data/project/Alzheimer/YSscript/EM_MRS/data/CLEMENT/simulation_2D/clone_4/5 
# --SCICLONE_DIR /data/project/Alzheimer/YSscript/EM_MRS/data/sciclone/simulation_2D/clone_4/5 
# --PYCLONEVI_DIR /data/project/Alzheimer/YSscript/EM_MRS/data/pyclone-vi/simulation_2D/clone_4/5 
# --QUANTUMCLONE_DIR /data/project/Alzheimer/YSscript/EM_MRS/data/quantumclone/simulation_2D/clone_4/5 
# --COMBINED_OUTPUT_DIR /data/project/Alzheimer/YSscript/EM_MRS/data/combinedoutput/simulation_2D/clone_4/5 
# --NUM_CLONE_TRIAL_START 2 --NUM_CLONE_TRIAL_END 7 --FP_RATIO 0.1 --DEPTH_CUTOFF 10 --VERBOSE 1 --TRIAL_NO 3 --RANDOM_SEED 5 --SCORING True --MODE Hard