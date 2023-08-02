import argparse, random, os
import numpy as np
import pandas as pd

SCRIPT_DIR = os.getcwd()

# python3 1.SimData_master_all.py --SIMDATA decoy  --FP_RATIO 0 --NUM_CLONE_START 2 --NUM_CLONE_END 7  --NUM_BLOCK_START 1 --NUM_BLOCK_END 3 --BENCHMARK_START 0 --BENCHMARK_END 9
# python3 1.SimData_master_all.py --SIMDATA decoy  --FP_RATIO 0.01 --NUM_CLONE_START 2 --NUM_CLONE_END 7  --NUM_BLOCK_START 1 --NUM_BLOCK_END 3 --BENCHMARK_START 0 --BENCHMARK_END 9
# python3 1.SimData_master_all.py --SIMDATA decoy  --FP_RATIO 0.03 --NUM_CLONE_START 2 --NUM_CLONE_END 7  --NUM_BLOCK_START 1 --NUM_BLOCK_END 3 --BENCHMARK_START 0 --BENCHMARK_END 9
# python3 1.SimData_master_all.py --SIMDATA decoy  --FP_RATIO 0.05 --NUM_CLONE_START 2 --NUM_CLONE_END 7  --NUM_BLOCK_START 1 --NUM_BLOCK_END 3 --BENCHMARK_START 0 --BENCHMARK_END 9
# python3 1.SimData_master_all.py --SIMDATA decoy  --FP_RATIO 0.1 --NUM_CLONE_START 2 --NUM_CLONE_END 7  --NUM_BLOCK_START 1 --NUM_BLOCK_END 3 --BENCHMARK_START 0 --BENCHMARK_END 9
# python3 1.SimData_master_all.py --SIMDATA lump  --FP_RATIO 0 --NUM_CLONE_START 2 --NUM_CLONE_END 7  --NUM_BLOCK_START 1 --NUM_BLOCK_END 3 --BENCHMARK_START 0 --BENCHMARK_END 9
# python3 1.SimData_master_all.py --SIMDATA lump  --FP_RATIO 0.01 --NUM_CLONE_START 2 --NUM_CLONE_END 7  --NUM_BLOCK_START 1 --NUM_BLOCK_END 3 --BENCHMARK_START 0 --BENCHMARK_END 9
# python3 1.SimData_master_all.py --SIMDATA lump  --FP_RATIO 0.03 --NUM_CLONE_START 2 --NUM_CLONE_END 7  --NUM_BLOCK_START 1 --NUM_BLOCK_END 3 --BENCHMARK_START 0 --BENCHMARK_END 9
# python3 1.SimData_master_all.py --SIMDATA lump  --FP_RATIO 0.05 --NUM_CLONE_START 2 --NUM_CLONE_END 7  --NUM_BLOCK_START 1 --NUM_BLOCK_END 3 --BENCHMARK_START 0 --BENCHMARK_END 9
# python3 1.SimData_master_all.py --SIMDATA lump  --FP_RATIO 0.1 --NUM_CLONE_START 2 --NUM_CLONE_END 7  --NUM_BLOCK_START 1 --NUM_BLOCK_END 3 --BENCHMARK_START 0 --BENCHMARK_END 9

if __name__ == "__main__":
    kwargs = {}

    parser = argparse.ArgumentParser( description='The below is usage direction.')
    parser.add_argument('--NUM_CLONE_START', type=int, default=2)
    parser.add_argument('--NUM_CLONE_END', type=int, default=4)
    parser.add_argument('--NUM_BLOCK_START', type=int, default=1)
    parser.add_argument('--NUM_BLOCK_END', type=int, default=3)
    parser.add_argument('--NUM_MUTATION', type=int, default=500)
    parser.add_argument('--FP_RATIO', type=float, default=0)
    parser.add_argument('--LOWVAF_RATIO', type=float, default=0)
    parser.add_argument('--DEPTH_MEAN', type=int, default=100)
    parser.add_argument('--DEPTH_SD', type=int, default=8)
    parser.add_argument('--BENCHMARK_START', type=int, default=0)
    parser.add_argument('--BENCHMARK_END', type=int, default=19)
    parser.add_argument('--SIMDATA', type=str, default="decoy")
    parser.add_argument('--MAXIMUM_NUM_PARENT', type=int, default = 1)
    
    args = parser.parse_args()

    kwargs["FP_RATIO"] = args.FP_RATIO
    FP_EXISTANCE = "noFP" if kwargs["FP_RATIO"] == 0 else "FP"
    kwargs["LOWVAF_RATIO"] = args.LOWVAF_RATIO
    kwargs["NUM_MUTATION"] = args.NUM_MUTATION
    NUM_MUTATION = args.NUM_MUTATION
    kwargs["DEPTH_MEAN"] = int(args.DEPTH_MEAN)
    kwargs["DEPTH_SD"] = int(args.DEPTH_SD)
    kwargs["BENCHMARK_START"] = int(args.BENCHMARK_START)
    kwargs["BENCHMARK_END"] = int(args.BENCHMARK_END)
    kwargs["SIMDATA"] = str(args.SIMDATA)
    kwargs["MAXIMUM_NUM_PARENT"] = int(args.MAXIMUM_NUM_PARENT)

    n  = 0 
    for NUM_CLONE in range(args.NUM_CLONE_START, args.NUM_CLONE_END + 1):
        kwargs["NUM_CLONE"] = NUM_CLONE
        for NUM_BLOCK in range(args.NUM_BLOCK_START, args.NUM_BLOCK_END + 1):
            kwargs["NUM_BLOCK"] = NUM_BLOCK

            hold_j = []

            for ii in range(kwargs["BENCHMARK_START"],  kwargs["BENCHMARK_END"] + 1):
                print("NUM_CLONE = {}\tNUM_BLOCK = {}\tBENCHMARK_NO = {}".format(NUM_CLONE, NUM_BLOCK, ii))
                kwargs["INPUT_TSV"] = "/data/project/Alzheimer/EM_cluster/EM_input/simulation_" + str(NUM_BLOCK) + "D/" + kwargs["SIMDATA"] + "/" + str(kwargs["FP_RATIO"]) + "/clone_" + str(NUM_CLONE) + "/" + str(NUM_BLOCK) + "D_clone" + str(NUM_CLONE) + "_" + str(ii) + ".txt"
                kwargs["NPVAF_DIR"] = "/data/project/Alzheimer/YSscript/EM_MRS/data/npvaf/simulation_" + str(NUM_BLOCK) + "D/" + kwargs["SIMDATA"] + "/" + str(kwargs["FP_RATIO"]) + "/clone_" + str(NUM_CLONE) + "/" + str(ii)
                kwargs["CLEMENT_DIR"] = "/data/project/Alzheimer/YSscript/EM_MRS/data/CLEMENT/simulation_" + str(NUM_BLOCK) + "D/" + kwargs["SIMDATA"] + "/" + str(kwargs["FP_RATIO"]) + "/clone_" + str(NUM_CLONE) + "/" + str(ii)
                kwargs["SCICLONE_DIR"] = "/data/project/Alzheimer/YSscript/EM_MRS/data/sciclone/simulation_" + str(NUM_BLOCK) + "D/" + kwargs["SIMDATA"] + "/" + str(kwargs["FP_RATIO"]) + "/clone_" + str(NUM_CLONE) + "/" + str(ii)
                kwargs["PYCLONEVI_DIR"] = "/data/project/Alzheimer/YSscript/EM_MRS/data/pyclone-vi/simulation_" + str(NUM_BLOCK) + "D/" + kwargs["SIMDATA"] + "/" + str(kwargs["FP_RATIO"]) + "/clone_" + str(NUM_CLONE) + "/" + str(ii)
                kwargs["QUANTUMCLONE_DIR"] = "/data/project/Alzheimer/YSscript/EM_MRS/data/quantumclone/simulation_" + str(NUM_BLOCK) + "D/" + kwargs["SIMDATA"] + "/" + str(kwargs["FP_RATIO"]) + "/clone_" + str(NUM_CLONE) + "/" + str(ii)
                kwargs["COMBINED_OUTPUT_DIR"] = "/data/project/Alzheimer/YSscript/EM_MRS/data/combinedoutput/simulation_" + str(NUM_BLOCK) + "D/" + kwargs["SIMDATA"] + "/" + str(kwargs["FP_RATIO"]) + "/clone_" + str(NUM_CLONE) + "/" + str(ii)

                for DIR in ["/".join(kwargs["INPUT_TSV"].split("/")[:-1]), kwargs["NPVAF_DIR"], kwargs["CLEMENT_DIR"], kwargs["SCICLONE_DIR"], kwargs["PYCLONEVI_DIR"], kwargs["QUANTUMCLONE_DIR"], kwargs["COMBINED_OUTPUT_DIR"]]:
                    # 얘는 공통분모 디렉토리를 지워버리면 안된다
                    if DIR == "/".join(kwargs["INPUT_TSV"].split("/")[:-1]):
                        if os.path.exists(kwargs["INPUT_TSV"]) == True:
                            os.system("rm -rf  " + kwargs["INPUT_TSV"])
                        if os.path.exists(DIR) == False:
                            os.system("mkdir -p " + DIR)
                    else:
                        if os.path.exists(DIR) == True:
                            os.system("rm -rf " + DIR)
                        if os.path.exists(DIR) == False:
                            os.system("mkdir -p " + DIR)



                # 0. Simulation dataset 생성 (decoy or sparse)
                COMPUTE_RANDOM = "cpu.q@compute" + str( random.randint (1,14) ).zfill(2)
                logPath = "/data/project/Alzheimer/YSscript/EM_MRS/log/simulation_" + str(NUM_BLOCK) + "D/" + kwargs["SIMDATA"] + "/" + str(kwargs["FP_RATIO"]) + "/clone_" + str(NUM_CLONE) + "/" + str(ii)
                os.system("rm -rf " + logPath)
                os.system("mkdir -p " + logPath)
                command1 = " ".join(["qsub -pe smp 1 -e", logPath, "-o", logPath, "-N SimDataFormation_" + str(NUM_BLOCK) + "D_clone_" + str(NUM_CLONE) + "_" + kwargs["SIMDATA"] + "_" + str(ii),
                                     "-q ", COMPUTE_RANDOM, 
                                     "1.SimData_pipe0_preparation.sh",
                                    "--SCRIPT_DIR", str(SCRIPT_DIR),  "--NUM_CLONE", str(NUM_CLONE),  "--NUM_BLOCK", str(NUM_BLOCK), "--NUM_MUTATION", str(NUM_MUTATION), "--FP_RATIO", str(kwargs["FP_RATIO"]),   "--LOWVAF_RATIO", str(kwargs["LOWVAF_RATIO"]),
                                     "--DEPTH_MEAN", str(kwargs["DEPTH_MEAN"]), "--DEPTH_SD", str(kwargs["DEPTH_SD"]),
                                     "--INPUT_TSV", kwargs["INPUT_TSV"], "--NPVAF_DIR", kwargs["NPVAF_DIR"],
                                     "--BENCHMARK_I", str(ii),
                                     "--SIMDATA", kwargs["SIMDATA"]
                                     ])
                #print (command1)
                n = n + 1
                os.system(command1)

                # 2. EM 돌리기 (2-7, Hard)
                logPath = "/data/project/Alzheimer/YSscript/EM_MRS/log/simulation_" + str(NUM_BLOCK) + "D/" + kwargs["SIMDATA"] + "/" + str(kwargs["FP_RATIO"]) + "/clone_" + str(NUM_CLONE) + "/" + str(ii)

                hold_j.append("simulation_" + str(NUM_BLOCK) + "D_clone_" + str(NUM_CLONE) + "_" + kwargs["SIMDATA"] + "_" + str(ii))
                command2 = " ".join(["qsub -pe smp 1 -e", logPath, "-o", logPath, "-N simulation_" + str(NUM_BLOCK) + "D_clone_" + str(NUM_CLONE) + "_" + kwargs["SIMDATA"] + "_" + str(ii),
                                     "-hold_jid SimDataFormation_" + str(NUM_BLOCK) + "D_clone_" + str(NUM_CLONE) + "_" + kwargs["SIMDATA"] + "_" + str(ii),
                                     "-q ", COMPUTE_RANDOM, 
                                    "1.SimData_pipe1_CLEMENT_bm.sh",
                                     "--INPUT_TSV", str(kwargs["INPUT_TSV"]), "--NPVAF_DIR", str(kwargs["NPVAF_DIR"]),
                                     "--CLEMENT_DIR", str(kwargs["CLEMENT_DIR"]),
                                     "--SCICLONE_DIR", str(kwargs["SCICLONE_DIR"]),
                                     "--PYCLONEVI_DIR", str(kwargs["PYCLONEVI_DIR"]),
                                     "--QUANTUMCLONE_DIR", str(kwargs["QUANTUMCLONE_DIR"]),
                                     "--COMBINED_OUTPUT_DIR", str(kwargs["COMBINED_OUTPUT_DIR"]),
                                     "--NUM_CLONE_TRIAL_START", str(2), "--NUM_CLONE_TRIAL_END", str(7),
                                     "--FP_RATIO", str(kwargs["FP_RATIO"]),
                                     "--MAXIMUM_NUM_PARENT", str( kwargs["MAXIMUM_NUM_PARENT"]) ,
                                     "--DEPTH_CUTOFF", str(10),  "--VERBOSE", str(1), "--TRIAL_NO", str(8), "--RANDOM_SEED", str(ii), "--SCORING True", "--MODE Both"
                                     ])
                # print(command2)
                os.system(command2)

            # 채점하기
            kwargs["COMBINED_OUTPUT_DIR"] = "/data/project/Alzheimer/YSscript/EM_MRS/data/combinedoutput/simulation_" + str(NUM_BLOCK) + "D/" + kwargs["SIMDATA"] + "/" +  str(kwargs["FP_RATIO"]) + "/clone_" + str(NUM_CLONE)
            kwargs["SAMPLENAME"] = "simulation_" + str(NUM_BLOCK) + "D_clone_" + str(NUM_CLONE)
            kwargs["OUTPUT_TTEST"] = kwargs["COMBINED_OUTPUT_DIR"] + "/ttest.txt"
            kwargs["OUTPUT_JPG"] = kwargs["COMBINED_OUTPUT_DIR"] +  "/benchmark.jpg"
            logPath = "/data/project/Alzheimer/YSscript/EM_MRS/log/simulation_" + str(NUM_BLOCK) + "D/" + kwargs["SIMDATA"] + "/" + str(kwargs["FP_RATIO"]) + "/clone_" + str(NUM_CLONE) + "/visualization"

            os.system("rm -rf " + logPath)
            os.system("mkdir -p " + logPath)

            command3 = " ".join(["qsub -pe smp 1 -e", logPath, "-o", logPath, "-N visualization",  "-hold_jid",  str(",".join(hold_j)),
                                 "-q ", COMPUTE_RANDOM, 
                                "1.SimData_pipe2_benchmark.sh",
                                 "--COMBINED_OUTPUT_DIR", str(kwargs["COMBINED_OUTPUT_DIR"]),
                                 "--SAMPLENAME", str(kwargs["SAMPLENAME"]),
                                 "--BENCHMARK_START", str(kwargs["BENCHMARK_START"]),
                                 "--BENCHMARK_END", str(kwargs["BENCHMARK_END"]),
                                 "--OUTPUT_TTEST", str(kwargs["OUTPUT_TTEST"]),
                                 "--OUTPUT_JPG", str(kwargs["OUTPUT_JPG"])
                                 ])
            ##print (command3)
            os.system(command3)
    print ("Total job = {}".format( n ))