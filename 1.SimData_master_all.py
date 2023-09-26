import os

SCRIPT_DIR = SCRIPT_DIR = os.path.dirname(__file__)
print (SCRIPT_DIR, "\n")


# python3 /data/project/Alzheimer/YSscript/cle/1.SimData_master_all.py

if __name__ == "__main__":
    kwargs = {}

    NUM_BLOCK_LIST = [ 1, 2, 3 ]             # 1, 2, 3
    NUM_MUTATION_LIST = [ 1000, 100 ]    # 100, 500, 100
    DEPTH_MEAN_LIST = [ 250]       # 250, 125, 30
    FP_RATIO_LIST = [ 0.0, 0.1  ]        # 0.0, 0.1
    SIMDATA_LIST = [ "decoy", "lump"] # "decoy", "lump"
    NUM_CLONE_LIST = [2, 3, 4, 5, 6, 7]      # 2, 3, 4, 5, 6, 7
    BENCHMARK_LIST = [0, 3]; kwargs["BENCHMARK_START"] = BENCHMARK_LIST[0];  kwargs["BENCHMARK_END"] = BENCHMARK_LIST[1]

    kwargs["NUM_CLONE_TRIAL_START"], kwargs["NUM_CLONE_TRIAL_END"] = 2, 7
    kwargs["MAXIMUM_NUM_PARENT"] = 0
    kwargs["TRIAL_NO"] = 8
    kwargs["MAKEONE_STRICT"] = 1
    kwargs["SCORING"] = "True"
    kwargs["MODE"] = "Both"
    kwargs["VERBOSE"] = 1
                           

    n  = 0 

    for NUM_BLOCK in NUM_BLOCK_LIST:
        kwargs["NUM_BLOCK"] = NUM_BLOCK
        for NUM_MUTATION in NUM_MUTATION_LIST:
            kwargs["NUM_MUTATION"] = NUM_MUTATION
            kwargs["MIN_CLUSTER_SIZE"] = 15 if NUM_MUTATION >= 100 else 5
            for DEPTH_MEAN in DEPTH_MEAN_LIST:
                kwargs["DEPTH_MEAN"] = DEPTH_MEAN
                kwargs["DEPTH_SD"] = 8 if DEPTH_MEAN >= 100 else 5
                kwargs["DEPTH_CUTOFF"] = 30 if DEPTH_MEAN >= 100 else 10
                kwargs["MAKEONE_STRICT"] = 1 if DEPTH_MEAN > 100 else 2

                for SIMDATA in SIMDATA_LIST:
                    kwargs["SIMDATA"] = SIMDATA
                    for FP_RATIO in FP_RATIO_LIST:
                        kwargs["FP_RATIO"] = FP_RATIO

                        print("\n======================\t1.SimData_{}D/n{}_{}x/{}/{}\t===============================".format( kwargs["NUM_BLOCK"], kwargs["NUM_MUTATION"], kwargs["DEPTH_MEAN"], kwargs["SIMDATA"], kwargs["FP_RATIO"]))
                        
                        hold_jj = []
                        for NUM_CLONE in NUM_CLONE_LIST:
                            kwargs["NUM_CLONE"] = NUM_CLONE
                            if NUM_BLOCK == 1:
                                kwargs["KMEANS_CLUSTERNO"] = 8 if NUM_CLONE >= 6 else 6
                                kwargs["TRIAL_NO"] = 5
                            if NUM_BLOCK >= 2:
                                kwargs["KMEANS_CLUSTERNO"] = 9 if NUM_CLONE >= 6 else 7
                                kwargs["TRIAL_NO"] = 8

                            hold_j = []
                            for ii in range(kwargs["BENCHMARK_START"],  kwargs["BENCHMARK_END"] + 1):
            
                                kwargs["INPUT_TSV"] = "/data/project/Alzheimer/CLEMENT/01.INPUT_TSV/1.SimData/SimData_" + str(NUM_BLOCK) + "D/n" + str(NUM_MUTATION) + "_" + str(DEPTH_MEAN)  + "x/" + str(SIMDATA) + "/" + str(FP_RATIO) + "/clone_" + str(NUM_CLONE) + "/" +  str(ii) + ".txt"
                                kwargs["NPVAF_DIR"] = "/data/project/Alzheimer/CLEMENT/02.npvaf/1.SimData/SimData_" + str(NUM_BLOCK) + "D/n" + str(NUM_MUTATION) + "_" + str(DEPTH_MEAN)  + "x/" + str(SIMDATA) + "/" + str(FP_RATIO) + "/clone_" + str(NUM_CLONE) + "/" +  str(ii) 
                                kwargs["COMBINED_OUTPUT_DIR"] = "/data/project/Alzheimer/CLEMENT/03.combinedoutput/1.SimData/SimData_" + str(NUM_BLOCK) + "D/n" + str(NUM_MUTATION) + "_" + str(DEPTH_MEAN)  + "x/" + str(SIMDATA) + "/" + str(FP_RATIO) + "/clone_" + str(NUM_CLONE) + "/" +  str(ii) 
                                kwargs["SIMPLE_KMEANS_DIR"] = "/data/project/Alzheimer/YSscript/cle/data/SIMPLE_KMEANS/1.SimData/SimData_" + str(NUM_BLOCK) + "D/n" + str(NUM_MUTATION) + "_" + str(DEPTH_MEAN)  + "x/" + str(SIMDATA) + "/" + str(FP_RATIO) + "/clone_" + str(NUM_CLONE) + "/" +  str(ii) 
                                kwargs["CLEMENT_DIR"] = "/data/project/Alzheimer/YSscript/cle/data/CLEMENT/1.SimData/SimData_" + str(NUM_BLOCK) + "D/n" + str(NUM_MUTATION) + "_" + str(DEPTH_MEAN)  + "x/" + str(SIMDATA) + "/" + str(FP_RATIO) + "/clone_" + str(NUM_CLONE) + "/" +  str(ii) 
                                kwargs["PYCLONEVI_DIR"] = "/data/project/Alzheimer/YSscript/cle/data/pyclone-vi/1.SimData/SimData_" + str(NUM_BLOCK) + "D/n" + str(NUM_MUTATION) + "_" + str(DEPTH_MEAN)  + "x/" + str(SIMDATA) + "/" + str(FP_RATIO) + "/clone_" + str(NUM_CLONE) + "/" +  str(ii) 
                                kwargs["SCICLONE_DIR"] = "/data/project/Alzheimer/YSscript/cle/data/sciclone/1.SimData/SimData_" + str(NUM_BLOCK) + "D/n" + str(NUM_MUTATION) + "_" + str(DEPTH_MEAN)  + "x/" + str(SIMDATA) + "/" + str(FP_RATIO) + "/clone_" + str(NUM_CLONE) + "/" +  str(ii) 
                                kwargs["QUANTUMCLONE_DIR"] = "/data/project/Alzheimer/YSscript/cle/data/quantumclone/1.SimData/SimData_" + str(NUM_BLOCK) + "D/n" + str(NUM_MUTATION) + "_" + str(DEPTH_MEAN)  + "x/" + str(SIMDATA) + "/" + str(FP_RATIO) + "/clone_" + str(NUM_CLONE) + "/" +  str(ii) 

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
                                logPath = "/data/project/Alzheimer/YSscript/cle/log/1.SimData/SimData_" + str(NUM_BLOCK) + "D/n" + str(NUM_MUTATION) + "_" + str(DEPTH_MEAN)  + "x/" + str(SIMDATA) + "/" + str(FP_RATIO) + "/clone_" + str(NUM_CLONE) + "/" +  str(ii) 
                                os.system("rm -rf " + logPath)
                                os.system("mkdir -p " + logPath)
                                command1 = " ".join([ "qsub -pe smp 1 -e", logPath, "-o", logPath, "-N SimData_Formation_" + str(NUM_BLOCK) + "D_n" + str(NUM_MUTATION) + "_" + str(DEPTH_MEAN)  + "x_" + str(SIMDATA) + "_" + str(FP_RATIO) + "_clone_" + str(NUM_CLONE) + "_" +  str(ii) ,
                                                    #"-q ", COMPUTE_RANDOM, 
                                                    SCRIPT_DIR + "/1.SimData_pipe0_preparation.sh",
                                                    "--SCRIPT_DIR", str(SCRIPT_DIR),  "--NUM_CLONE", str(NUM_CLONE),  "--NUM_BLOCK", str(NUM_BLOCK), "--NUM_MUTATION", str(NUM_MUTATION), "--FP_RATIO", str(kwargs["FP_RATIO"]),  
                                                    "--DEPTH_MEAN", str(kwargs["DEPTH_MEAN"]), "--DEPTH_SD", str(kwargs["DEPTH_SD"]), "--DEPTH_CUTOFF", str(kwargs["DEPTH_CUTOFF"]),
                                                    "--INPUT_TSV", kwargs["INPUT_TSV"], "--NPVAF_DIR", kwargs["NPVAF_DIR"],
                                                    "--BENCHMARK_I", str(ii),
                                                    "--SIMDATA", kwargs["SIMDATA"]
                                                    ])
                                n = n + 1
                                os.system(command1)

                                #1. EM 돌리기 
                                hold_j.append( "SimData_EM_" + str(NUM_BLOCK) + "D_n" + str(NUM_MUTATION) + "_" + str(DEPTH_MEAN)  + "x_" + str(SIMDATA) + "_" + str(FP_RATIO) + "_clone_" + str(NUM_CLONE) + "_" +  str(ii) )
                                command2 = " ".join(["qsub -pe smp 1 -e", logPath, "-o", logPath, "-N SimData_EM_" + str(NUM_BLOCK) + "D_n" + str(NUM_MUTATION) + "_" + str(DEPTH_MEAN)  + "x_" + str(SIMDATA) + "_" + str(FP_RATIO) + "_clone_" + str(NUM_CLONE) + "_" +  str(ii) ,
                                                     "-hold_jid SimData_Formation_" + str(NUM_BLOCK) + "D_n" + str(NUM_MUTATION) + "_" + str(DEPTH_MEAN)  + "x_" + str(SIMDATA) + "_" + str(FP_RATIO) + "_clone_" + str(NUM_CLONE) + "_" +  str(ii),
                                                     #"-q ", COMPUTE_RANDOM, 
                                                    SCRIPT_DIR + "/1.SimData_pipe1_CLEMENT_bm.sh",
                                                    "--SCRIPT_DIR", str(SCRIPT_DIR),
                                                     "--INPUT_TSV", str(kwargs["INPUT_TSV"]), 
                                                     "--NPVAF_DIR", str(kwargs["NPVAF_DIR"]),
                                                     "--COMBINED_OUTPUT_DIR", str(kwargs["COMBINED_OUTPUT_DIR"]),
                                                     "--SIMPLE_KMEANS_DIR", str(kwargs["SIMPLE_KMEANS_DIR"]),
                                                     "--CLEMENT_DIR", str(kwargs["CLEMENT_DIR"]),
                                                     "--SCICLONE_DIR", str(kwargs["SCICLONE_DIR"]),
                                                     "--PYCLONEVI_DIR", str(kwargs["PYCLONEVI_DIR"]),
                                                     "--QUANTUMCLONE_DIR", str(kwargs["QUANTUMCLONE_DIR"]),
                                                     "--NUM_CLONE_TRIAL_START", str(kwargs["NUM_CLONE_TRIAL_START"]), "--NUM_CLONE_TRIAL_END", str(kwargs["NUM_CLONE_TRIAL_END"]),
                                                     "--FP_RATIO", str(kwargs["FP_RATIO"]),
                                                     "--MAXIMUM_NUM_PARENT", str( kwargs["MAXIMUM_NUM_PARENT"]) ,
                                                     "--DEPTH_CUTOFF", str(kwargs["DEPTH_CUTOFF"]),  
                                                     "--VERBOSE", str( kwargs["VERBOSE"]), 
                                                     "--TRIAL_NO", str( kwargs["TRIAL_NO"]), 
                                                     "--MAKEONE_STRICT", str ( kwargs["MAKEONE_STRICT"] ),
                                                     "--RANDOM_SEED", str(ii), 
                                                     "--SCORING", str(kwargs["SCORING"]), 
                                                     "--KMEANS_CLUSTERNO", str(kwargs["KMEANS_CLUSTERNO"]), 
                                                     "--MIN_CLUSTER_SIZE", str(kwargs["MIN_CLUSTER_SIZE"]), 
                                                     "--MODE", str(kwargs["MODE"]),
                                                     "--VERBOSE", str(kwargs["VERBOSE"])
                                                     ])
                                os.system(command2)

                            #2. 각 iteratio 단위에서 채점하고 benchmark 하기
                            kwargs["COMBINED_OUTPUT_DIR"] = "/data/project/Alzheimer/CLEMENT/03.combinedoutput/1.SimData/SimData_" + str(NUM_BLOCK) + "D/n" + str(NUM_MUTATION) + "_" + str(DEPTH_MEAN)  + "x/" + str(SIMDATA) + "/" + str(FP_RATIO) + "/clone_" + str(NUM_CLONE) 
                            kwargs["SAMPLENAME"] = "SimData_" + str(NUM_BLOCK) + "D/n" + str(NUM_MUTATION) + "_" + str(DEPTH_MEAN)  + "x/" + str(SIMDATA) + "/" + str(FP_RATIO) + "/clone_" + str(NUM_CLONE) 
                            kwargs["OUTPUT_TTEST"] = kwargs["COMBINED_OUTPUT_DIR"] + "/ttest.txt"
                            kwargs["OUTPUT_JPG"] = kwargs["COMBINED_OUTPUT_DIR"] +  "/bm.jpg"

                            logPath = "/data/project/Alzheimer/YSscript/cle/log/1.SimData/SimData_" + str(NUM_BLOCK) + "D/n" + str(NUM_MUTATION) + "_" + str(DEPTH_MEAN)  + "x/" + str(SIMDATA) + "/" + str(FP_RATIO) + "/clone_" + str(NUM_CLONE) + "/benchmark"
                            os.system("rm -rf " + logPath)
                            os.system("mkdir -p " + logPath)

                            hold_jj.append ( "bm_SimData_" + str(NUM_BLOCK) + "D_n" + str(NUM_MUTATION) + "_" + str(DEPTH_MEAN)  + "x_" + str(SIMDATA) + "_" + str(FP_RATIO) + "_clone_" + str(NUM_CLONE) )
                            command3 = " ".join(["qsub -pe smp 1 -e", logPath, "-o", logPath, 
                                                 "-N bm_" + "SimData_" + str(NUM_BLOCK) + "D_n" + str(NUM_MUTATION) + "_" + str(DEPTH_MEAN)  + "x_" + str(SIMDATA) + "_" + str(FP_RATIO) + "_clone_" + str(NUM_CLONE) ,  
                                                 "-hold_jid",  str(",".join(hold_j)),
                                                str(SCRIPT_DIR) + "/1.SimData_pipe2_benchmark.sh",
                                                "--SCRIPT_DIR", str(SCRIPT_DIR),
                                                "--COMBINED_OUTPUT_DIR", str(kwargs["COMBINED_OUTPUT_DIR"]),
                                                "--SAMPLENAME", str(kwargs["SAMPLENAME"]),
                                                "--BENCHMARK_START", str(kwargs["BENCHMARK_START"]),
                                                "--BENCHMARK_END", str(kwargs["BENCHMARK_END"]),
                                                "--OUTPUT_TTEST", str(kwargs["OUTPUT_TTEST"]),
                                                "--OUTPUT_JPG", str(kwargs["OUTPUT_JPG"]),
                                                "--FP_RATIO", str(kwargs["FP_RATIO"])
                                                ])
                            #print (command3)
                            os.system(command3)
                            n += 1

                        #3. Clone 2 ~7  단위에서의 benchmark visualization
                        INPUT_DIR = "/".join(kwargs["COMBINED_OUTPUT_DIR"].split("/")[ : -1])
                        logPath = "/data/project/Alzheimer/YSscript/cle/log/1.SimData/SimData_" + str(NUM_BLOCK) + "D/n" + str(NUM_MUTATION) + "_" + str(DEPTH_MEAN)  + "x/" + str(SIMDATA) + "/" + str(FP_RATIO) + "/BM"
        
                        os.system("rm -rf " + logPath)
                        os.system("mkdir -p " + logPath)
                        
                        CONDITIONNAME = "SimData_" + str(NUM_BLOCK) + "D_n" + str(NUM_MUTATION) + "_" + str(DEPTH_MEAN)  + "x/" + str(SIMDATA) + "/" + str(FP_RATIO) 

                        command4 = " ".join(  ["qsub -pe smp 1", "-e", logPath, "-o", logPath, 
                                            "-N", "BM_SimData_" + str(NUM_BLOCK) + "D_" + CONDITIONNAME.replace("/", "_"),
                                            "-hold_jid",  str(",".join(hold_jj)), " ",
                                            SCRIPT_DIR + "/1.SimData_pipe3_benchmark.sh",
                                                                "--SCRIPT_DIR", str(SCRIPT_DIR), 
                                                                "--INPUT_DIR", str(INPUT_DIR) , 
                                                                "--CONDITIONNAME", str(CONDITIONNAME), 
                                                                "--BENCHMARK_START", str( kwargs["BENCHMARK_START"]), 
                                                                "--BENCHMARK_END", str( kwargs["BENCHMARK_END"]), 
                                                                "--OUTPUT_MS_JPG", str(INPUT_DIR) + "/BM_MS.jpg", 
                                                                "--OUTPUT_EC_JPG", str(INPUT_DIR) + "/BM_EC.jpg",
                                                                "--OUTPUT_FINAL_JPG", str(INPUT_DIR) + "/BM_FINAL.jpg",
                                                                "--OUTPUT_FINAL_TABLE", str(INPUT_DIR) + "/BM_FINAL.tsv"
                                                                    ]  )
                        
                        # print (command4)
                        os.system(command4)
                        n += 1
    print ("Total job = {}".format( n ))