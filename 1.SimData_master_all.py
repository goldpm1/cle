import argparse, random, os
import numpy as np
import pandas as pd

SCRIPT_DIR = "/".join ( os.path.abspath(__file__).split("/")[:-1] )
print (SCRIPT_DIR, "\n")

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

    NUM_BLOCK_LIST = [1, 2]
    NUM_MUTATION_LIST = [500, 100]
    DEPTH_MEAN_LIST = [100, 30]
    FP_RATIO_LIST = [0.0, 0.1]
    SIMDATA_LIST = ["decoy", "lump"]
    NUM_CLONE_LIST = [2, 4]
    BENCHMARK_LIST = [0, 4]; kwargs["BENCHMARK_START"] = BENCHMARK_LIST[0];  kwargs["BENCHMARK_END"] = BENCHMARK_LIST[1]

    kwargs["MAXIMUM_NUM_PARENT"] = 0

    n  = 0 


    for NUM_BLOCK in NUM_BLOCK_LIST:
        kwargs["NUM_BLOCK"] = NUM_BLOCK
        for NUM_MUTATION in NUM_MUTATION_LIST:
            kwargs["NUM_MUTATION"] = NUM_MUTATION
            for DEPTH_MEAN in DEPTH_MEAN_LIST:
                kwargs["DEPTH_MEAN"] = DEPTH_MEAN
                kwargs["DEPTH_SD"] = 8 if DEPTH_MEAN >= 100 else 5
                for SIMDATA in SIMDATA_LIST:
                    kwargs["SIMDATA"] = SIMDATA
                    for FP_RATIO in FP_RATIO_LIST:
                        kwargs["FP_RATIO"] = FP_RATIO

                        print("====================== 1.SimData_{}D/n{}_{}x/{}/{} ===============================".format( kwargs["NUM_BLOCK"], kwargs["NUM_MUTATION"], kwargs["DEPTH_MEAN"], kwargs["SIMDATA"], kwargs["FP_RATIO"]))
                            
                        for NUM_CLONE in range( NUM_CLONE_LIST[0], NUM_CLONE_LIST[1] ):
                            kwargs["NUM_CLONE"] = NUM_CLONE

                            hold_j = []
                            for ii in range(kwargs["BENCHMARK_START"],  kwargs["BENCHMARK_END"] + 1):
            
                                kwargs["INPUT_TSV"] = "/data/project/Alzheimer/CLEMENT/01.INPUT_TSV/1.SimData/SimData_" + str(NUM_BLOCK) + "D/n" + str(NUM_MUTATION) + "_" + str(DEPTH_MEAN)  + "x/" + str("SIMDATA") + "/" + str(FP_RATIO) + "/clone_" + str(NUM_CLONE) + "/" +  str(ii) + ".txt"
                                kwargs["NPVAF_DIR"] = "/data/project/Alzheimer/CLEMENT/02.npvaf/1.SimData/SimData_" + str(NUM_BLOCK) + "D/n" + str(NUM_MUTATION) + "_" + str(DEPTH_MEAN)  + "x/" + str("SIMDATA") + "/" + str(FP_RATIO) + "/clone_" + str(NUM_CLONE) + "/" +  str(ii) 
                                kwargs["COMBINED_OUTPUT_DIR"] = "/data/project/Alzheimer/CLEMENT/03.combinedoutput/1.SimData/SimData_" + str(NUM_BLOCK) + "D/n" + str(NUM_MUTATION) + "_" + str(DEPTH_MEAN)  + "x/" + str("SIMDATA") + "/" + str(FP_RATIO) + "/clone_" + str(NUM_CLONE) + "/" +  str(ii) 
                                kwargs["SIMPLE_KMEANS_DIR"] = "/data/project/Alzheimer/YSscript/cle/data/SIMPLE_KMEANS/1.SimData/SimData_" + str(NUM_BLOCK) + "D/n" + str(NUM_MUTATION) + "_" + str(DEPTH_MEAN)  + "x/" + str("SIMDATA") + "/" + str(FP_RATIO) + "/clone_" + str(NUM_CLONE) + "/" +  str(ii) 
                                kwargs["CLEMENT_DIR"] = "/data/project/Alzheimer/YSscript/cle/data/CLEMENT/1.SimData/SimData_" + str(NUM_BLOCK) + "D/n" + str(NUM_MUTATION) + "_" + str(DEPTH_MEAN)  + "x/" + str("SIMDATA") + "/" + str(FP_RATIO) + "/clone_" + str(NUM_CLONE) + "/" +  str(ii) 
                                kwargs["SCICLONE_DIR"] = "/data/project/Alzheimer/YSscript/cle/data/sciclone/1.SimData/SimData_" + str(NUM_BLOCK) + "D/n" + str(NUM_MUTATION) + "_" + str(DEPTH_MEAN)  + "x/" + str("SIMDATA") + "/" + str(FP_RATIO) + "/clone_" + str(NUM_CLONE) + "/" +  str(ii) 
                                kwargs["PYCLONEVI_DIR"] = "/data/project/Alzheimer/YSscript/cle/data/pyclone-vi/1.SimData/SimData_" + str(NUM_BLOCK) + "D/n" + str(NUM_MUTATION) + "_" + str(DEPTH_MEAN)  + "x/" + str("SIMDATA") + "/" + str(FP_RATIO) + "/clone_" + str(NUM_CLONE) + "/" +  str(ii) 
                                kwargs["QUANTUMCLONE_DIR"] = "/data/project/Alzheimer/YSscript/cle/data/quantumclone/1.SimData/SimData_" + str(NUM_BLOCK) + "D/n" + str(NUM_MUTATION) + "_" + str(DEPTH_MEAN)  + "x/" + str("SIMDATA") + "/" + str(FP_RATIO) + "/clone_" + str(NUM_CLONE) + "/" +  str(ii) 

                                for DIR in ["/".join(kwargs["INPUT_TSV"].split("/")[:-1]), kwargs["NPVAF_DIR"], kwargs["SIMPLE_KMEANS_DIR"],  kwargs["CLEMENT_DIR"], kwargs["SCICLONE_DIR"], kwargs["PYCLONEVI_DIR"], kwargs["QUANTUMCLONE_DIR"], kwargs["COMBINED_OUTPUT_DIR"]]:
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


                                # 0. Simulation dataset 생성 
                                #COMPUTE_RANDOM = "cpu.q@compute" + str( random.randint (1,14) ).zfill(2)
                                logPath = "/data/project/Alzheimer/YSscript/cle/log/1.SimData/SimData_" + str(NUM_BLOCK) + "D/n" + str(NUM_MUTATION) + "_" + str(DEPTH_MEAN)  + "x/" + str("SIMDATA") + "/" + str(FP_RATIO) + "/clone_" + str(NUM_CLONE) + "/" +  str(ii) 
                                os.system("rm -rf " + logPath)
                                os.system("mkdir -p " + logPath)
                                command1 = " ".join([ "qsub -pe smp 1 -e", logPath, "-o", logPath, "-N SimData_Formation_" + str(NUM_BLOCK) + "D_n" + str(NUM_MUTATION) + "_" + str(DEPTH_MEAN)  + "x_" + str("SIMDATA") + "_" + str(FP_RATIO) + "_clone_" + str(NUM_CLONE) + "_" +  str(ii) ,
                                                    #"-q ", COMPUTE_RANDOM, 
                                                    SCRIPT_DIR + "/1.SimData_pipe0_preparation.sh",
                                                    "--SCRIPT_DIR", str(SCRIPT_DIR),  "--NUM_CLONE", str(NUM_CLONE),  "--NUM_BLOCK", str(NUM_BLOCK), "--NUM_MUTATION", str(NUM_MUTATION), "--FP_RATIO", str(kwargs["FP_RATIO"]),  
                                                    "--DEPTH_MEAN", str(kwargs["DEPTH_MEAN"]), "--DEPTH_SD", str(kwargs["DEPTH_SD"]),
                                                    "--INPUT_TSV", kwargs["INPUT_TSV"], "--NPVAF_DIR", kwargs["NPVAF_DIR"],
                                                    "--BENCHMARK_I", str(ii),
                                                    "--SIMDATA", kwargs["SIMDATA"]
                                                    ])
                                #print (command1)
                                n = n + 1
                                os.system(command1)

    #             # 2. EM 돌리기 (2-7, Hard)
    #             logPath = "/data/project/Alzheimer/YSscript/EM_MRS/log/simulation_" + str(NUM_BLOCK) + "D/" + kwargs["SIMDATA"] + "/" + str(kwargs["FP_RATIO"]) + "/clone_" + str(NUM_CLONE) + "/" + str(ii)

    #             hold_j.append("simulation_" + str(NUM_BLOCK) + "D_clone_" + str(NUM_CLONE) + "_" + kwargs["SIMDATA"] + "_" + str(ii))
    #             command2 = " ".join(["qsub -pe smp 1 -e", logPath, "-o", logPath, "-N simulation_" + str(NUM_BLOCK) + "D_clone_" + str(NUM_CLONE) + "_" + kwargs["SIMDATA"] + "_" + str(ii),
    #                                  "-hold_jid SimDataFormation_" + str(NUM_BLOCK) + "D_clone_" + str(NUM_CLONE) + "_" + kwargs["SIMDATA"] + "_" + str(ii),
    #                                  "-q ", COMPUTE_RANDOM, 
    #                                 SCRIPT_DIR + "/1.SimData_pipe1_CLEMENT_bm.sh",
    #                                  "--INPUT_TSV", str(kwargs["INPUT_TSV"]), "--NPVAF_DIR", str(kwargs["NPVAF_DIR"]),
    #                                  "--CLEMENT_DIR", str(kwargs["CLEMENT_DIR"]),
    #                                  "--SCICLONE_DIR", str(kwargs["SCICLONE_DIR"]),
    #                                  "--PYCLONEVI_DIR", str(kwargs["PYCLONEVI_DIR"]),
    #                                  "--QUANTUMCLONE_DIR", str(kwargs["QUANTUMCLONE_DIR"]),
    #                                  "--COMBINED_OUTPUT_DIR", str(kwargs["COMBINED_OUTPUT_DIR"]),
    #                                  "--NUM_CLONE_TRIAL_START", str(2), "--NUM_CLONE_TRIAL_END", str(7),
    #                                  "--FP_RATIO", str(kwargs["FP_RATIO"]),
    #                                  "--MAXIMUM_NUM_PARENT", str( kwargs["MAXIMUM_NUM_PARENT"]) ,
    #                                  "--DEPTH_CUTOFF", str(10),  "--VERBOSE", str(1), "--TRIAL_NO", str(8), "--RANDOM_SEED", str(ii), "--SCORING True", "--MODE Both"
    #                                  ])
    #             # print(command2)
    #             os.system(command2)

    #         # 채점하기
    #         kwargs["COMBINED_OUTPUT_DIR"] = "/data/project/Alzheimer/YSscript/EM_MRS/data/combinedoutput/simulation_" + str(NUM_BLOCK) + "D/" + kwargs["SIMDATA"] + "/" +  str(kwargs["FP_RATIO"]) + "/clone_" + str(NUM_CLONE)
    #         kwargs["SAMPLENAME"] = "simulation_" + str(NUM_BLOCK) + "D_clone_" + str(NUM_CLONE)
    #         kwargs["OUTPUT_TTEST"] = kwargs["COMBINED_OUTPUT_DIR"] + "/ttest.txt"
    #         kwargs["OUTPUT_JPG"] = kwargs["COMBINED_OUTPUT_DIR"] +  "/benchmark.jpg"
    #         logPath = "/data/project/Alzheimer/YSscript/EM_MRS/log/simulation_" + str(NUM_BLOCK) + "D/" + kwargs["SIMDATA"] + "/" + str(kwargs["FP_RATIO"]) + "/clone_" + str(NUM_CLONE) + "/visualization"

    #         os.system("rm -rf " + logPath)
    #         os.system("mkdir -p " + logPath)

    #         command3 = " ".join(["qsub -pe smp 1 -e", logPath, "-o", logPath, "-N visualization",  "-hold_jid",  str(",".join(hold_j)),
    #                              "-q ", COMPUTE_RANDOM, 
    #                             "1.SimData_pipe2_benchmark.sh",
    #                              "--COMBINED_OUTPUT_DIR", str(kwargs["COMBINED_OUTPUT_DIR"]),
    #                              "--SAMPLENAME", str(kwargs["SAMPLENAME"]),
    #                              "--BENCHMARK_START", str(kwargs["BENCHMARK_START"]),
    #                              "--BENCHMARK_END", str(kwargs["BENCHMARK_END"]),
    #                              "--OUTPUT_TTEST", str(kwargs["OUTPUT_TTEST"]),
    #                              "--OUTPUT_JPG", str(kwargs["OUTPUT_JPG"])
    #                              ])
    #         ##print (command3)
    #         os.system(command3)
    # print ("Total job = {}".format( n ))