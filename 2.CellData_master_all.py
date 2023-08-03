import os, glob,re

SCRIPT_DIR = SCRIPT_DIR = os.path.dirname(__file__)
print (SCRIPT_DIR, "\n")


# python3 /data/project/Alzheimer/YSscript/cle/2.CellData_master_all.py

def M1only(SAMPLENAME, NUM_BLOCK):
    global MRS_set
    #MRS_set = [['M1-2', 'M1-5', 'M1-6', 'M1-7', 'M1-8'],  ['M2-2', 'M2-4', 'M2-6', 'M2-8', 'M2-10']]
    MRS_set = [['M1-2', 'M1-5', 'M1-6', 'M1-7', 'M1-8'] ]

    if NUM_BLOCK == 1:
        for i in range(len(MRS_set)):
            if (SAMPLENAME.split("_")[0] in MRS_set[i]): 
                return True
        # if 'M3' in SAMPLENAME.split("_")[0]:   # M1, M2를 다 돌려보자
        #     return False
        # return True
    elif NUM_BLOCK == 2:
        for i in range(len(MRS_set)):
            if (SAMPLENAME.split("_")[0] in MRS_set[i]) & (SAMPLENAME.split("_")[1] in MRS_set[i]):
                return True
    elif NUM_BLOCK == 3:
        for i in range(len(MRS_set)):
            if (SAMPLENAME.split("_")[0] in MRS_set[i]) & (SAMPLENAME.split("_")[1] in MRS_set[i]) & (SAMPLENAME.split("_")[2] in MRS_set[i]):
                return True

    return False



if __name__ == "__main__":
    kwargs = {}
    
    NUM_BLOCK_LIST = [2]             # 1, 2, 3
    NUM_MUTATION_LIST = [500]    # 500, 100
    DEPTH_MEAN_LIST = [125]       # 125, 30
    NUM_PARENT_LIST = [  0 ]       # 0 , 1
    FP_RATIO_LIST = [0.0 ]        # 0.0, 0.1
    AXIS_RATIO_LIST = [0.0, 0.1 ]        # 0.0, 0.1
    BENCHMARK_LIST = [0, 4]; kwargs["BENCHMARK_START"] = BENCHMARK_LIST[0];  kwargs["BENCHMARK_END"] = BENCHMARK_LIST[1]

    kwargs["NUM_CLONE_TRIAL_START"], kwargs["NUM_CLONE_TRIAL_END"] = 2, 7
    kwargs["MAXIMUM_NUM_PARENT"] = 1
    kwargs["TRIAL_NO"] = 8
    kwargs["MAKEONE_STRICT"] = 1
    kwargs["VERBOSE"] = 1
    kwargs["FP_USEALL"] = "False"
    kwargs["SCORING"] = "True"
    kwargs["MODE"] = "Both"
    kwargs["MIN_CLUSTER_SIZE"] = 15
    kwargs["KMEANS_CLUSTERNO"] = 8
    kwargs["PARENT_RATIO"] = 0



    n = 0


    for NUM_BLOCK in NUM_BLOCK_LIST:
        kwargs["NUM_BLOCK"] = NUM_BLOCK
        for NUM_MUTATION in NUM_MUTATION_LIST:
            kwargs["NUM_MUTATION"] = NUM_MUTATION
            for DEPTH_MEAN in DEPTH_MEAN_LIST:
                kwargs["DEPTH_MEAN"] = DEPTH_MEAN
                kwargs["DEPTH_CUTOFF"] = 30 if DEPTH_MEAN >= 100 else 10
                for NUM_PARENT in NUM_PARENT_LIST:
                    kwargs["NUM_PARENT"] = NUM_PARENT
                    

                    for FP_RATIO in FP_RATIO_LIST:
                        kwargs["FP_RATIO"] = FP_RATIO
                        for AXIS_RATIO in AXIS_RATIO_LIST:
                            kwargs["AXIS_RATIO"] = AXIS_RATIO
                            print("\n======================\t2.CellData_{}D/n{}_{}x/parent_{}/fp_{}/axis_{}\t===============================".format( kwargs["NUM_BLOCK"], kwargs["NUM_MUTATION"], kwargs["DEPTH_MEAN"], kwargs["NUM_PARENT"], kwargs["FP_RATIO"], kwargs["AXIS_RATIO"] ))

                            DIR = "/data/project/Alzheimer/CLEMENT/01.INPUT_TSV/2.CellData/CellData_" +  str(kwargs["NUM_BLOCK"]) + "D/" + str(kwargs["DEPTH_MEAN"]) + "x"
                            DIR_LIST = sorted (  glob.glob(DIR + "/*.txt") ) 

                            hold_jj = []
                            for INPUT_TSV in DIR_LIST:
                                kwargs["INPUT_TSV"] = INPUT_TSV
                                SAMPLENAME = INPUT_TSV.split("/")[-1].split("_input.txt") [0]     # 'M1-5_M1-8_input'
                                kwargs["SAMPLENAME"] = SAMPLENAME

                                if M1only(SAMPLENAME, kwargs["NUM_BLOCK"]) == True:
                                    print("# n = {},  {}".format(n + 1, SAMPLENAME, NUM_PARENT))


                                    # 1. EM 돌리기 
                                    hold_j = []
                                    for ii in range(kwargs["BENCHMARK_START"], kwargs["BENCHMARK_END"] + 1):
                                        #print("\fp_{}/axis_{}/{}/ii = {}".format( kwargs["FP_RATIO"], kwargs["AXIS_RATIO"], kwargs["SAMPLENAME"], ii ))

                                        kwargs["NPVAF_DIR"] = "/data/project/Alzheimer/CLEMENT/02.npvaf/2.CellData/CellData_" + str(NUM_BLOCK) + "D/n" + str(NUM_MUTATION) + "_" + str(DEPTH_MEAN)  + "x/parent_" + str(NUM_PARENT) + "/fp_" + str(FP_RATIO) + "/axis_" + str(AXIS_RATIO) + "/" + SAMPLENAME + "/" + str(ii) 
                                        kwargs["COMBINED_OUTPUT_DIR"] = "/data/project/Alzheimer/CLEMENT/03.combinedoutput/2.CellData/CellData_" + str(NUM_BLOCK) + "D/n" + str(NUM_MUTATION) + "_" + str(DEPTH_MEAN)  + "x/parent_" + str(NUM_PARENT) + "/fp_" + str(FP_RATIO) + "/axis_" + str(AXIS_RATIO) + "/" + SAMPLENAME + "/" + str(ii) 
                                        
                                        kwargs["SIMPLE_KMEANS_DIR"] = "/data/project/Alzheimer/YSscript/cle/data/SIMPLE_KMEANS/2.CellData/CellData_" +  str(NUM_BLOCK) + "D/n" + str(NUM_MUTATION) + "_" + str(DEPTH_MEAN)  + "x/parent_" + str(NUM_PARENT) + "/fp_" + str(FP_RATIO) + "/axis_" + str(AXIS_RATIO) + "/" + SAMPLENAME + "/" + str(ii) 
                                        kwargs["CLEMENT_DIR"] = "/data/project/Alzheimer/YSscript/cle/data/CLEMENT/2.CellData/CellData_" +  str(NUM_BLOCK) + "D/n" + str(NUM_MUTATION) + "_" + str(DEPTH_MEAN)  + "x/parent_" + str(NUM_PARENT) + "/fp_" + str(FP_RATIO) + "/axis_" + str(AXIS_RATIO) + "/" + SAMPLENAME + "/" + str(ii) 
                                        kwargs["PYCLONEVI_DIR"] = "/data/project/Alzheimer/YSscript/cle/data/pyclone-vi/2.CellData/CellData_" +  str(NUM_BLOCK) + "D/n" + str(NUM_MUTATION) + "_" + str(DEPTH_MEAN)  + "x/parent_" + str(NUM_PARENT) + "/fp_" + str(FP_RATIO) + "/axis_" + str(AXIS_RATIO) + "/" + SAMPLENAME + "/" + str(ii) 
                                        kwargs["SCICLONE_DIR"] = "/data/project/Alzheimer/YSscript/cle/data/sciclone/2.CellData/CellData_" +  str(NUM_BLOCK) + "D/n" + str(NUM_MUTATION) + "_" + str(DEPTH_MEAN)  + "x/parent_" + str(NUM_PARENT) + "/fp_" + str(FP_RATIO) + "/axis_" + str(AXIS_RATIO) + "/" + SAMPLENAME + "/" + str(ii) 
                                        kwargs["QUANTUMCLONE_DIR"] = "/data/project/Alzheimer/YSscript/cle/data/quantumclone/2.CellData/CellData_" +  str(NUM_BLOCK) + "D/n" + str(NUM_MUTATION) + "_" + str(DEPTH_MEAN)  + "x/parent_" + str(NUM_PARENT) + "/fp_" + str(FP_RATIO) + "/axis_" + str(AXIS_RATIO) + "/" + SAMPLENAME + "/" + str(ii) 

                                        for DIR in [ kwargs["NPVAF_DIR"], kwargs["SIMPLE_KMEANS_DIR"], kwargs["CLEMENT_DIR"], kwargs["SCICLONE_DIR"], kwargs["PYCLONEVI_DIR"], kwargs["QUANTUMCLONE_DIR"], kwargs["COMBINED_OUTPUT_DIR"]  ]:
                                            if os.path.exists ( DIR ) == True:
                                                os.system("rm -rf " + DIR  )
                                            if os.path.exists ( DIR ) == False:
                                                os.system("mkdir -p " + DIR  )
                                            
                                        logPath = "/data/project/Alzheimer/YSscript/cle/log/2.CellData/CellData_" +  str(NUM_BLOCK) + "D/n" + str(NUM_MUTATION) + "_" + str(DEPTH_MEAN)  + "x/parent_" + str(NUM_PARENT) + "/fp_" + str(FP_RATIO) + "/axis_" + str(AXIS_RATIO) + "/" + SAMPLENAME + "/" + str(ii) 
                                        os.system("rm -rf " + logPath)
                                        os.system("mkdir -p " + logPath)

                                        hold_j.append( "CellData_EM_" + str(NUM_BLOCK) + "D_n" + str(NUM_MUTATION) + "_" + str(DEPTH_MEAN)  + "x_parent_" + str(NUM_PARENT) + "_fp_" + str(FP_RATIO) + "_axis_" + str(AXIS_RATIO) + "_" + SAMPLENAME + "_" + str(ii)  )
                                        command = " ".join(["qsub -pe smp 1", "-e", logPath, "-o", logPath, 
                                                                        "-N", "CellData_EM_" + str(NUM_BLOCK) + "D_n" + str(NUM_MUTATION) + "_" + str(DEPTH_MEAN)  + "x_parent_" + str(NUM_PARENT) + "_fp_" + str(FP_RATIO) + "_axis_" + str(AXIS_RATIO) + "_" + SAMPLENAME + "_" + str(ii) ,
                                                                        SCRIPT_DIR + "/2.CellData_pipe1_CLEMENT_bm.sh",
                                                                        str(SCRIPT_DIR), str(INPUT_TSV),  str(kwargs["MODE"]),  str(kwargs["NUM_CLONE_TRIAL_START"]),  str( kwargs["NUM_CLONE_TRIAL_END"]),  
                                                                        str(kwargs["NUM_MUTATION"]), str(kwargs["AXIS_RATIO"]),  str(kwargs["PARENT_RATIO"]),  str( kwargs["NUM_PARENT"]),  str(kwargs["FP_RATIO"]),  str(kwargs["FP_USEALL"]),
                                                                        str(kwargs["TRIAL_NO"]), str(kwargs["DEPTH_CUTOFF"]),  str( kwargs["MIN_CLUSTER_SIZE"]),  str(kwargs["VERBOSE"]),
                                                                        str(kwargs["KMEANS_CLUSTERNO"]),  str(ii), str( kwargs["SAMPLENAME"]), str(kwargs["BENCHMARK_END"]),
                                                                        str(kwargs["NPVAF_DIR"]), str(kwargs["SIMPLE_KMEANS_DIR"]),  str(kwargs["CLEMENT_DIR"]), str(kwargs["SCICLONE_DIR"]), str( kwargs["PYCLONEVI_DIR"]), str(kwargs["QUANTUMCLONE_DIR"]), str(kwargs["COMBINED_OUTPUT_DIR"]),
                                                                        str(kwargs["SCORING"]), str(kwargs["MAKEONE_STRICT"]), str(kwargs["MAXIMUM_NUM_PARENT"])])
                                        os.system(command)
                                        n += 1

                                    # 2. i 개 trial의 benchmark
                                    logPath = "/data/project/Alzheimer/YSscript/cle/log/2.CellData/CellData_" + str(NUM_BLOCK) + "D/n" + str(NUM_MUTATION) + "_" + str(DEPTH_MEAN)  + "x/parent_" + str(NUM_PARENT) + "/fp_" + str(FP_RATIO) + "/axis_" + str(AXIS_RATIO) + "/" + str(SAMPLENAME) + "/bm"
                                    os.system("rm -rf " + logPath)
                                    os.system("mkdir -p " + logPath)

                                    INPUT_DIR = "/".join(kwargs["COMBINED_OUTPUT_DIR"].split("/")[ : -1])
                                    CONDITIONNAME = "n" + str(NUM_MUTATION) + "_" + str(DEPTH_MEAN)  + "x/parent_" + str(NUM_PARENT) + "/fp_" + str(FP_RATIO) + "/axis_" + str(AXIS_RATIO)
                                    hold_jj.append ( "bm_CellData_" + str(NUM_BLOCK) + "D_n" + str(NUM_MUTATION) + "_" + str(DEPTH_MEAN)  + "x_parent_" + str(NUM_PARENT) + "_fp_" + str(FP_RATIO) + "_axis_" + str(AXIS_RATIO) + "_" + str(SAMPLENAME) )
                                    command = " ".join(["qsub -pe smp 1", "-e", logPath, "-o", logPath, 
                                                                    "-N", "bm_CellData_" + str(NUM_BLOCK) + "D_n" + str(NUM_MUTATION) + "_" + str(DEPTH_MEAN)  + "x_parent_" + str(NUM_PARENT) + "_fp_" + str(FP_RATIO) + "_axis_" + str(AXIS_RATIO) + "_" + str(SAMPLENAME),
                                                                    "-hold_jid",  str(",".join(hold_j)), 
                                                                    SCRIPT_DIR + "/2.CellData_pipe2_benchmark.sh",
                                                                    "--SCRIPT_DIR", str(SCRIPT_DIR), 
                                                                    "--INPUT_DIR", str(INPUT_DIR) , 
                                                                    "--SAMPLENAME", str(SAMPLENAME), 
                                                                    "--CONDITIONNAME", str(CONDITIONNAME), 
                                                                    "--BENCHMARK_START", str( kwargs["BENCHMARK_START"]), 
                                                                    "--BENCHMARK_END", str( kwargs["BENCHMARK_END"]), 
                                                                    "--OUTPUT_JPG", str(INPUT_DIR) + "/bm.jpg"
                                                                    "--OUTPUT_TTEST", str(INPUT_DIR) + "/ttest.txt" ])
                                    os.system(command)
                                    n += 1

            


                            #3. 큰 단위에서의 benchmark visualization
                            OUTPUT_FILENAME = INPUT_DIR + "/" + str(kwargs["NUM_MUTATION"]) + "_" + str(kwargs["NUM_PARENT"]) + "_" + str(kwargs["FP_RATIO"]) + "_" + str(kwargs["AXIS_RATIO"]) + "_total.jpg"
                            logPath = "/data/project/Alzheimer/YSscript/cle/log/2.CellData/CellData_" + str(NUM_BLOCK) + "D/n" + str(NUM_MUTATION) + "_" + str(DEPTH_MEAN)  + "x/parent_" + str(NUM_PARENT) + "/fp_" + str(FP_RATIO) + "/axis_" + str(AXIS_RATIO) + "/BM"
                            os.system("rm -rf " + logPath)
                            os.system("mkdir -p " + logPath)
                            
                            INPUT_DIR = "/".join(kwargs["COMBINED_OUTPUT_DIR"].split("/")[ : -2])
                            CONDITIONNAME = "n" + str(NUM_MUTATION) + "_" + str(DEPTH_MEAN)  + "x/parent_" + str(NUM_PARENT) + "/fp_" + str(FP_RATIO) + "/axis_" + str(AXIS_RATIO)

                            command = " ".join(  ["qsub -pe smp 1", "-e", logPath, "-o", logPath, 
                                                "-N", "BM_CellData_" + str(NUM_BLOCK) + "D_n" + str(NUM_MUTATION) + "_" + str(DEPTH_MEAN)  + "x_parent_" + str(NUM_PARENT) + "_fp_" + str(FP_RATIO) + "_axis_" + str(AXIS_RATIO) + "_" + str(SAMPLENAME),
                                                "-hold_jid",  str(",".join(hold_jj)), 
                                                SCRIPT_DIR + "/2.CellData_pipe3_benchmark.sh",
                                                                    "--SCRIPT_DIR", str(SCRIPT_DIR), 
                                                                    "--INPUT_DIR", str(INPUT_DIR) , 
                                                                    "--SAMPLENAME", str(SAMPLENAME), 
                                                                    "--CONDITIONNAME", str(CONDITIONNAME), 
                                                                    "--BENCHMARK_START", str( kwargs["BENCHMARK_START"]), 
                                                                    "--BENCHMARK_END", str( kwargs["BENCHMARK_END"]), 
                                                                    "--OUTPUT_JPG", str(INPUT_DIR) + "/BM.jpg"]  )
                            
                            #print (command)
                            os.system(command)
                            n += 1


    print ("Total job = {}".format(n))