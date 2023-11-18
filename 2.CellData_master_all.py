import os, glob,re, random

SCRIPT_DIR = SCRIPT_DIR = os.path.dirname(__file__)
print (SCRIPT_DIR, "\n")


# python3 /data/project/Alzheimer/YSscript/cle/2.CellData_master_all.py

def Mingle_intra(SAMPLENAME, NUM_BLOCK):
    global MRS_set
    
    #MRS_set = [['M1-2', 'M1-5', 'M1-6', 'M1-7', 'M1-8'],  ['M2-2', 'M2-4', 'M2-6', 'M2-8', 'M2-10']]
    MRS_set = [['M1-2', 'M1-4', 'M1-6', 'M1-8'],  ['M2-2', 'M2-4', 'M2-6', 'M2-8'] ]


    if NUM_BLOCK == 1:
        for i in range(len(MRS_set)):
            if (SAMPLENAME.split("_")[0] in MRS_set[i]): 
                return True
    elif NUM_BLOCK == 2:
        for i in range(len(MRS_set)):
            if (SAMPLENAME.split("_")[0] in MRS_set[i]) & (SAMPLENAME.split("_")[1] in MRS_set[i]):
                return True
    elif NUM_BLOCK == 3:
        for i in range(len(MRS_set)):
            if (SAMPLENAME.split("_")[0] in MRS_set[i]) & (SAMPLENAME.split("_")[1] in MRS_set[i]) & (SAMPLENAME.split("_")[2] in MRS_set[i]):
                return True

    return False

def Mingle_inter (SAMPLENAME, NUM_BLOCK):
    global MRS_set
    #MRS_set = [['M1-1', 'M1-4', 'M1-5', 'M1-6', 'M1-8'],  ['M2-2', 'M2-4', 'M2-6', 'M2-10', 'M2-12'], ['M3-1', 'M3-2', 'M3-4', 'M3-5', 'M3-11']]
    MRS_set = [['M1-2', 'M1-4', 'M1-6', 'M1-8'],  ['M2-2', 'M2-4', 'M2-6', 'M2-8'] ]

    if NUM_BLOCK == 1:
        for i in range(len(MRS_set)):
            if (SAMPLENAME.split("_")[0] in MRS_set[i]): 
                return True
        return False
    elif NUM_BLOCK == 2:
        if (SAMPLENAME.split("_")[0] in MRS_set[0]) &  (SAMPLENAME.split("_")[1] in MRS_set[1]):
            return True
        else:
            return False
    elif NUM_BLOCK == 3:
        if (SAMPLENAME.split("_")[0] in MRS_set[0]) &  (SAMPLENAME.split("_")[1] in MRS_set[0]) &   (SAMPLENAME.split("_")[2] in MRS_set[1]):
            return True
        elif (SAMPLENAME.split("_")[0] in MRS_set[0]) &  (SAMPLENAME.split("_")[1] in MRS_set[1]) &   (SAMPLENAME.split("_")[2] in MRS_set[1]):
            return True
        else:
            return False


if __name__ == "__main__":
    kwargs = {}
    
    NUM_BLOCK_LIST = [  2  ]             # 1, 2, 3
    NUM_MUTATION_LIST = [  500  ]    # 1000, 500, 100
    DEPTH_MEAN_LIST = [ 125  ]       # 250, 125, 30
    NUM_PARENT_LIST = [ 0  ]       # 0 , 1
    FP_RATIO_LIST = [  0.0 , 0.1  ]        # 0.0, 0.025, 0.05, 0.075, 0.1
    AXIS_RATIO_LIST = [ -1 ]        # -1, 0.0, 0.2
    BENCHMARK_LIST = [ 0, 3 ]; kwargs["BENCHMARK_START"] = BENCHMARK_LIST[0];  kwargs["BENCHMARK_END"] = BENCHMARK_LIST[1]

    kwargs["NUM_CLONE_TRIAL_START"], kwargs["NUM_CLONE_TRIAL_END"] = 2, 7
    kwargs["PARENT_RATIO"] = 0
    kwargs["FP_USEALL"] = "False"
    kwargs["SCORING"] = "True"
    kwargs["MODE"] = "Both"
    kwargs["IMAGE_FORMAT"] = "jpg"
    kwargs["VISUALIZATION"] = "True"
    kwargs["VERBOSE"] = 1

    #COMPUTE_RANDOM = "cpu.q@compute" + str( random.choice ( list ( range (1,15) ) + list ( range (16, 19) ) ) ).zfill(2)            

    n = 0

    DELETE_PREVIOUS_HISTORY = input("DELETE_PREVIOUS_HISTORY?: ")
    DELETE_PREVIOUS_HISTORY = True if DELETE_PREVIOUS_HISTORY in ["True", "Treu", "true", "T", "Y", "Yes"] else False
    print ( DELETE_PREVIOUS_HISTORY )

    for NUM_BLOCK in NUM_BLOCK_LIST:
        kwargs["NUM_BLOCK"] = NUM_BLOCK
        for NUM_MUTATION in NUM_MUTATION_LIST:
            kwargs["NUM_MUTATION"] = NUM_MUTATION
            kwargs["MIN_CLUSTER_SIZE"] = int (NUM_MUTATION / 33)
            for DEPTH_MEAN in DEPTH_MEAN_LIST:
                kwargs["DEPTH_MEAN"] = DEPTH_MEAN
                kwargs["DEPTH_CUTOFF"] = 30 if DEPTH_MEAN > 100 else 10
                kwargs["MAKEONE_STRICT"] = 1 if DEPTH_MEAN > 100 else 2
                for NUM_PARENT in NUM_PARENT_LIST:
                    kwargs["NUM_PARENT"] =  kwargs["MAXIMUM_NUM_PARENT"] = NUM_PARENT
                    for FP_RATIO in FP_RATIO_LIST:
                        kwargs["FP_RATIO"] = FP_RATIO
                        kwargs["FP_PRIOR"] = 0.01 + FP_RATIO
                        AXIS_RATIO_LIST_TEMP = [-1] if NUM_BLOCK == 1 else AXIS_RATIO_LIST      # 1D에서는 AXIS의 개념이 딱히 없다
                        if NUM_BLOCK == 1:
                            kwargs["KMEANS_CLUSTERNO"] = 6
                            kwargs["NUM_CLONE_TRIAL_END"] = 5
                            kwargs["TRIAL_NO"] = 5
                        elif NUM_BLOCK == 2:
                            kwargs["KMEANS_CLUSTERNO"] = 7 + NUM_PARENT
                            kwargs["NUM_CLONE_TRIAL_END"] = 7
                            kwargs["TRIAL_NO"] = 12
                        elif NUM_BLOCK == 3:
                            kwargs["KMEANS_CLUSTERNO"] = 7 + NUM_PARENT
                            kwargs["NUM_CLONE_TRIAL_END"] = 7
                            kwargs["TRIAL_NO"] = 15
                        for AXIS_RATIO in AXIS_RATIO_LIST_TEMP:
                            kwargs["AXIS_RATIO"] = AXIS_RATIO
                            #print("\n======================\t2.CellData_{}D/n{}_{}x/parent_{}/fp_{}/axis_{}\t===============================".format( kwargs["NUM_BLOCK"], kwargs["NUM_MUTATION"], kwargs["DEPTH_MEAN"], kwargs["NUM_PARENT"], kwargs["FP_RATIO"], kwargs["AXIS_RATIO"] ))

                            DIR = "/data/project/Alzheimer/CLEMENT/01.INPUT_TSV/2.CellData/CellData_" +  str(kwargs["NUM_BLOCK"]) + "D/" + str(kwargs["DEPTH_MEAN"]) + "x"
                            DIR_LIST = sorted (  glob.glob(DIR + "/*.txt") ) 

                            if DELETE_PREVIOUS_HISTORY == True:
                                if os.path.exists ( "/data/project/Alzheimer/CLEMENT/03.combinedoutput/2.CellData/CellData_" + str(NUM_BLOCK) + "D/n" + str(NUM_MUTATION) + "_" + str(DEPTH_MEAN)  + "x/parent_" + str(NUM_PARENT) + "/fp_" + str(FP_RATIO) + "/axis_" + str(AXIS_RATIO) + "/BM_EC.jpg" ) == True:
                                    os.system("rm -rf " + "/data/project/Alzheimer/CLEMENT/03.combinedoutput/2.CellData/CellData_" + str(NUM_BLOCK) + "D/n" + str(NUM_MUTATION) + "_" + str(DEPTH_MEAN)  + "x/parent_" + str(NUM_PARENT) + "/fp_" + str(FP_RATIO) + "/axis_" + str(AXIS_RATIO) + "/BM_EC.jpg"  )
                                    os.system("rm -rf " + "/data/project/Alzheimer/CLEMENT/03.combinedoutput/2.CellData/CellData_" + str(NUM_BLOCK) + "D/n" + str(NUM_MUTATION) + "_" + str(DEPTH_MEAN)  + "x/parent_" + str(NUM_PARENT) + "/fp_" + str(FP_RATIO) + "/axis_" + str(AXIS_RATIO) + "/BM_FINAL.jpg"  )
                                    os.system("rm -rf " + "/data/project/Alzheimer/CLEMENT/03.combinedoutput/2.CellData/CellData_" + str(NUM_BLOCK) + "D/n" + str(NUM_MUTATION) + "_" + str(DEPTH_MEAN)  + "x/parent_" + str(NUM_PARENT) + "/fp_" + str(FP_RATIO) + "/axis_" + str(AXIS_RATIO) + "/BM_FINAL.tsv"  )
                                    os.system("rm -rf " + "/data/project/Alzheimer/CLEMENT/03.combinedoutput/2.CellData/CellData_" + str(NUM_BLOCK) + "D/n" + str(NUM_MUTATION) + "_" + str(DEPTH_MEAN)  + "x/parent_" + str(NUM_PARENT) + "/fp_" + str(FP_RATIO) + "/axis_" + str(AXIS_RATIO) + "/BM_MS.jpg"  )

                            hold_jj = []
                            for INPUT_TSV in DIR_LIST:
                                kwargs["INPUT_TSV"] = INPUT_TSV
                                SAMPLENAME = INPUT_TSV.split("/")[-1].split("_input.txt") [0]     # 'M1-5_M1-8_input'
                                kwargs["SAMPLENAME"] = SAMPLENAME

                                if (NUM_BLOCK == 2) & ("M1-6" in SAMPLENAME):
                                    kwargs["KMEANS_CLUSTERNO"] = 8 + NUM_PARENT
                                
                                if (Mingle_inter(SAMPLENAME, kwargs["NUM_BLOCK"]) == True) | (Mingle_intra(SAMPLENAME, kwargs["NUM_BLOCK"]) == True):

                                    # 1. EM 돌리기 
                                    hold_j = []
                                    for ii in range(kwargs["BENCHMARK_START"], kwargs["BENCHMARK_END"] + 1):
                                        kwargs["NPVAF_DIR"] = "/data/project/Alzheimer/CLEMENT/02.npvaf/2.CellData/CellData_" + str(NUM_BLOCK) + "D/n" + str(NUM_MUTATION) + "_" + str(DEPTH_MEAN)  + "x/parent_" + str(NUM_PARENT) + "/fp_" + str(FP_RATIO) + "/axis_" + str(AXIS_RATIO) + "/" + SAMPLENAME + "/" + str(ii) 
                                        kwargs["COMBINED_OUTPUT_DIR"] = "/data/project/Alzheimer/CLEMENT/03.combinedoutput/2.CellData/CellData_" + str(NUM_BLOCK) + "D/n" + str(NUM_MUTATION) + "_" + str(DEPTH_MEAN)  + "x/parent_" + str(NUM_PARENT) + "/fp_" + str(FP_RATIO) + "/axis_" + str(AXIS_RATIO) + "/" + SAMPLENAME + "/" + str(ii) 
                                        
                                        kwargs["CLEMENT_DIR"] = "/data/project/Alzheimer/YSscript/cle/data/CLEMENT/2.CellData/CellData_" +  str(NUM_BLOCK) + "D/n" + str(NUM_MUTATION) + "_" + str(DEPTH_MEAN)  + "x/parent_" + str(NUM_PARENT) + "/fp_" + str(FP_RATIO) + "/axis_" + str(AXIS_RATIO) + "/" + SAMPLENAME + "/" + str(ii) 
                                        kwargs["SIMPLE_KMEANS_DIR"] = "/data/project/Alzheimer/YSscript/cle/data/SIMPLE_KMEANS/2.CellData/CellData_" +  str(NUM_BLOCK) + "D/n" + str(NUM_MUTATION) + "_" + str(DEPTH_MEAN)  + "x/parent_" + str(NUM_PARENT) + "/fp_" + str(FP_RATIO) + "/axis_" + str(AXIS_RATIO) + "/" + SAMPLENAME + "/" + str(ii) 
                                        kwargs["PYCLONEVI_DIR"] = "/data/project/Alzheimer/YSscript/cle/data/pyclone-vi/2.CellData/CellData_" +  str(NUM_BLOCK) + "D/n" + str(NUM_MUTATION) + "_" + str(DEPTH_MEAN)  + "x/parent_" + str(NUM_PARENT) + "/fp_" + str(FP_RATIO) + "/axis_" + str(AXIS_RATIO) + "/" + SAMPLENAME + "/" + str(ii) 
                                        kwargs["SCICLONE_DIR"] = "/data/project/Alzheimer/YSscript/cle/data/sciclone/2.CellData/CellData_" +  str(NUM_BLOCK) + "D/n" + str(NUM_MUTATION) + "_" + str(DEPTH_MEAN)  + "x/parent_" + str(NUM_PARENT) + "/fp_" + str(FP_RATIO) + "/axis_" + str(AXIS_RATIO) + "/" + SAMPLENAME + "/" + str(ii) 
                                        kwargs["QUANTUMCLONE_DIR"] = "/data/project/Alzheimer/YSscript/cle/data/quantumclone/2.CellData/CellData_" +  str(NUM_BLOCK) + "D/n" + str(NUM_MUTATION) + "_" + str(DEPTH_MEAN)  + "x/parent_" + str(NUM_PARENT) + "/fp_" + str(FP_RATIO) + "/axis_" + str(AXIS_RATIO) + "/" + SAMPLENAME + "/" + str(ii) 

                                        if DELETE_PREVIOUS_HISTORY == True:
                                            print ("Delete previous directory")
                                            for subdir in [kwargs["NPVAF_DIR"] , kwargs["COMBINED_OUTPUT_DIR"], kwargs["CLEMENT_DIR"], kwargs["SIMPLE_KMEANS_DIR"], kwargs["PYCLONEVI_DIR"], kwargs["SCICLONE_DIR"], kwargs["QUANTUMCLONE_DIR"] ] :
                                                if os.path.exists ( subdir ) == True:
                                                    os.system("rm -rf " + subdir  )
                                                if os.path.exists ( subdir ) == False:
                                                    os.system("mkdir -p " + subdir  )
                                            
                                        logPath = "/data/project/Alzheimer/YSscript/cle/log/2.CellData/CellData_" +  str(NUM_BLOCK) + "D/n" + str(NUM_MUTATION) + "_" + str(DEPTH_MEAN)  + "x/parent_" + str(NUM_PARENT) + "/fp_" + str(FP_RATIO) + "/axis_" + str(AXIS_RATIO) + "/" + SAMPLENAME + "/" + str(ii) 
                                        os.system("rm -rf " + logPath)
                                        os.system("mkdir -p " + logPath)
                                        logPath = "/data/project/Alzheimer/YSscript/cle/log/2.CellData/CellData_" +  str(NUM_BLOCK) + "D/n" + str(NUM_MUTATION) + "_" + str(DEPTH_MEAN)  + "x/parent_" + str(NUM_PARENT) + "/fp_" + str(FP_RATIO) + "/axis_" + str(AXIS_RATIO) + "/" + SAMPLENAME + "/" + str(ii) 
                                        logPath_copy = kwargs["COMBINED_OUTPUT_DIR"] 

                                        hold_j.append( "CellData_EM_" + str(NUM_BLOCK) + "D_n" + str(NUM_MUTATION) + "_" + str(DEPTH_MEAN)  + "x_parent_" + str(NUM_PARENT) + "_fp_" + str(FP_RATIO) + "_axis_" + str(AXIS_RATIO) + "_" + SAMPLENAME + "_" + str(ii)  )
                                        COMPUTE_RANDOM = "cpu.q@compute" + str( random.choice ( list ( range (1,14) ) + list ( range (16, 19) ) ) ).zfill(2)
                                
                                        command = " ".join(["qsub -pe smp 1", "-e", logPath, "-o", logPath, 
                                                                        "-N", "CellData_EM_" + str(NUM_BLOCK) + "D_n" + str(NUM_MUTATION) + "_" + str(DEPTH_MEAN)  + "x_parent_" + str(NUM_PARENT) + "_fp_" + str(FP_RATIO) + "_axis_" + str(AXIS_RATIO) + "_" + SAMPLENAME + "_" + str(ii) ,
                                                                        #"-q", COMPUTE_RANDOM,
                                                                        SCRIPT_DIR + "/2.CellData_pipe1_CLEMENT_bm.sh",
                                                                        str(SCRIPT_DIR), str(INPUT_TSV),  str(kwargs["MODE"]),  str(kwargs["NUM_CLONE_TRIAL_START"]),  str( kwargs["NUM_CLONE_TRIAL_END"]),  
                                                                        str(kwargs["NUM_MUTATION"]), str(kwargs["AXIS_RATIO"]),  str(kwargs["PARENT_RATIO"]),  str( kwargs["NUM_PARENT"]),  str(kwargs["FP_RATIO"]),  str(kwargs["FP_USEALL"]), str(kwargs["FP_PRIOR"]),
                                                                        str(kwargs["TRIAL_NO"]), str(kwargs["DEPTH_CUTOFF"]),  str( kwargs["MIN_CLUSTER_SIZE"]),  str(kwargs["VERBOSE"]),
                                                                        str(kwargs["KMEANS_CLUSTERNO"]),  str(ii), str( kwargs["SAMPLENAME"]), str(kwargs["BENCHMARK_END"]),
                                                                        str(kwargs["NPVAF_DIR"]), str(kwargs["SIMPLE_KMEANS_DIR"]),  str(kwargs["CLEMENT_DIR"]), str(kwargs["SCICLONE_DIR"]), str( kwargs["PYCLONEVI_DIR"]), str(kwargs["QUANTUMCLONE_DIR"]), str(kwargs["COMBINED_OUTPUT_DIR"]),
                                                                        str(kwargs["SCORING"]), str(kwargs["MAKEONE_STRICT"]), str(kwargs["MAXIMUM_NUM_PARENT"]), str(kwargs["IMAGE_FORMAT"]), str (kwargs["VISUALIZATION"]),
                                                                         str(logPath), str(logPath_copy) ])
                                        os.system(command)
                                        n += 1

                                    # 2. i 개 trial의 benchmark
                                    # logPath = "/data/project/Alzheimer/YSscript/cle/log/2.CellData/CellData_" + str(NUM_BLOCK) + "D/n" + str(NUM_MUTATION) + "_" + str(DEPTH_MEAN)  + "x/parent_" + str(NUM_PARENT) + "/fp_" + str(FP_RATIO) + "/axis_" + str(AXIS_RATIO) + "/" + str(SAMPLENAME) + "/bm"
                                    # os.system("rm -rf " + logPath)
                                    # os.system("mkdir -p " + logPath)

                                    if DELETE_PREVIOUS_HISTORY == True:
                                        if os.path.exists ( "/data/project/Alzheimer/CLEMENT/03.combinedoutput/2.CellData/CellData_" + str(NUM_BLOCK) + "D/n" + str(NUM_MUTATION) + "_" + str(DEPTH_MEAN)  + "x/parent_" + str(NUM_PARENT) + "/fp_" + str(FP_RATIO) + "/axis_" + str(AXIS_RATIO) + "/" + SAMPLENAME + "/bm.jpg" ) == True:
                                            os.system("rm -rf " + "/data/project/Alzheimer/CLEMENT/03.combinedoutput/2.CellData/CellData_" + str(NUM_BLOCK) + "D/n" + str(NUM_MUTATION) + "_" + str(DEPTH_MEAN)  + "x/parent_" + str(NUM_PARENT) + "/fp_" + str(FP_RATIO) + "/axis_" + str(AXIS_RATIO) + "/" + SAMPLENAME + "/bm.jpg"  )


                                    INPUT_DIR = "/".join(kwargs["COMBINED_OUTPUT_DIR"].split("/")[ : -1])
                                    logPath = INPUT_DIR
                                    os.system("rm -rf " + logPath + "/bm_CellData_" + str(NUM_BLOCK) + "D_n" + str(NUM_MUTATION) + "_" + str(DEPTH_MEAN)  + "x_parent_" + str(NUM_PARENT) + "_fp_" + str(FP_RATIO) + "_axis_" + str(AXIS_RATIO) + "_" + str(SAMPLENAME) + "*" )
                                    CONDITIONNAME = "n" + str(NUM_MUTATION) + "_" + str(DEPTH_MEAN)  + "x/parent_" + str(NUM_PARENT) + "/fp_" + str(FP_RATIO) + "/axis_" + str(AXIS_RATIO)
                                    hold_jj.append ( "bm_CellData_" + str(NUM_BLOCK) + "D_n" + str(NUM_MUTATION) + "_" + str(DEPTH_MEAN)  + "x_parent_" + str(NUM_PARENT) + "_fp_" + str(FP_RATIO) + "_axis_" + str(AXIS_RATIO) + "_" + str(SAMPLENAME) )
                                    command = " ".join(["qsub -pe smp 1", "-e", logPath, "-o", logPath, 
                                                                    "-N", "bm_CellData_" + str(NUM_BLOCK) + "D_n" + str(NUM_MUTATION) + "_" + str(DEPTH_MEAN)  + "x_parent_" + str(NUM_PARENT) + "_fp_" + str(FP_RATIO) + "_axis_" + str(AXIS_RATIO) + "_" + str(SAMPLENAME),
                                                                    "-hold_jid",  str(",".join(hold_j)), 
                                                                    #"-q", COMPUTE_RANDOM,
                                                                    SCRIPT_DIR + "/2.CellData_pipe2_benchmark.sh",
                                                                    "--SCRIPT_DIR", str(SCRIPT_DIR), 
                                                                    "--INPUT_DIR", str(INPUT_DIR) , 
                                                                    "--LOG_DIR", "/data/project/Alzheimer/YSscript/cle/log/2.CellData/CellData_" +  str(NUM_BLOCK) + "D/n" + str(NUM_MUTATION) + "_" + str(DEPTH_MEAN)  + "x/parent_" + str(NUM_PARENT) + "/fp_" + str(FP_RATIO) + "/axis_" + str(AXIS_RATIO) + "/" + SAMPLENAME , 
                                                                    "--SAMPLENAME", str(SAMPLENAME), 
                                                                    "--CONDITIONNAME", str(CONDITIONNAME), 
                                                                    "--BENCHMARK_START",  str(0),   #str( kwargs["BENCHMARK_START"]), 
                                                                    "--BENCHMARK_END", str( kwargs["BENCHMARK_END"]), 
                                                                    "--OUTPUT_JPG", str(INPUT_DIR) + "/bm.jpg",
                                                                    "--OUTPUT_TTEST", str(INPUT_DIR) + "/ttest.txt",
                                                                    "--OUTPUT_RESULT_TABLE",  str(INPUT_DIR) + "/bm.tsv" ])
                                    os.system(command)
                                    print ("")
                                    #n += 1

            


                            #3. 큰 단위에서의 benchmark visualization
                            # logPath = "/data/project/Alzheimer/YSscript/cle/log/2.CellData/CellData_" + str(NUM_BLOCK) + "D/n" + str(NUM_MUTATION) + "_" + str(DEPTH_MEAN)  + "x/parent_" + str(NUM_PARENT) + "/fp_" + str(FP_RATIO) + "/axis_" + str(AXIS_RATIO) + "/BM"
                            # os.system("rm -rf " + logPath)
                            # os.system("mkdir -p " + logPath)
                            
                            INPUT_DIR = "/".join(kwargs["COMBINED_OUTPUT_DIR"].split("/")[ : -2])
                            logPath = INPUT_DIR
                            CONDITIONNAME = "CellData_" + str(NUM_BLOCK) + "D_n" + str(NUM_MUTATION) + "_" + str(DEPTH_MEAN)  + "x/parent_" + str(NUM_PARENT) + "/fp_" + str(FP_RATIO) + "/axis_" + str(AXIS_RATIO)
                            os.system("rm -rf " + logPath + "/BM_CellData_" + str(NUM_BLOCK) + "D_" + CONDITIONNAME.replace("/", "_") + "*" )

                            command = " ".join(  ["qsub -pe smp 1", "-e", logPath, "-o", logPath, 
                                                "-N", "BM_CellData_" + str(NUM_BLOCK) + "D_" + CONDITIONNAME.replace("/", "_"),
                                                "-hold_jid",  str(",".join(hold_jj)), " ",
                                                #"-q", COMPUTE_RANDOM,
                                                SCRIPT_DIR + "/2.CellData_pipe3_benchmark.sh",
                                                                    "--SCRIPT_DIR", str(SCRIPT_DIR), 
                                                                    "--INPUT_DIR", str(INPUT_DIR) , 
                                                                    "--CONDITIONNAME", str(CONDITIONNAME), 
                                                                    "--BENCHMARK_START", str(0),  # str( kwargs["BENCHMARK_START"]), 
                                                                    "--BENCHMARK_END", str( kwargs["BENCHMARK_END"]), 
                                                                    "--OUTPUT_MS_JPG", str(INPUT_DIR) + "/BM_MS.jpg", 
                                                                    "--OUTPUT_EC_JPG", str(INPUT_DIR) + "/BM_EC.jpg",
                                                                    "--OUTPUT_FINAL_JPG", str(INPUT_DIR) + "/BM_FINAL.jpg",
                                                                    "--OUTPUT_FINAL_TABLE", str(INPUT_DIR) + "/BM_FINAL.tsv"
                                                                      ]  )
                            
                            os.system(command)
                            print ("\n\n")
                            #n += 1


    print ("Total job = {}".format(n))