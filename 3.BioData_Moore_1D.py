import numpy as np
import pandas as pd
import glob, subprocess, os, filetype
import os, re

# python3 /data/project/Alzheimer/YSscript/cle/3.BioData_Moore_1D.py


def out(command): 
    result = subprocess.run(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True, shell=True) 
    return result.stdout.rstrip("\n")


########################################################################################################################################################################

def Other_tissues (SCRIPT_DIR ) :
    DIR = "/data/project/Alzheimer/CLEMENT/01.INPUT_TSV/3.BioData/Moore_1D/2.woMosaic_ver1"        # AG : 3.woMosaic_ver2
    DONOR_LIST = sorted(glob.glob(DIR + "/*"))
    DONOR_LIST = [i.split("/")[-1] for i in DONOR_LIST]
    n = 0

    for DONOR in DONOR_LIST:               # PD42566
        TISSUE_LIST = sorted( glob.glob (DIR + "/" + DONOR + "/*") ) 
        TISSUE_LIST = [i.split("/")[-1] for i in TISSUE_LIST]

        for TISSUE in TISSUE_LIST:               # colon_crypt, pancreas_islet
        #for TISSUE in [ "bronchus_epithelium", "adrenal_gland_zona_glomerulosa", "adrenal_gland_zona_fasciculata",  "adrenal_gland_zona_reticularis" ]:

            SAMPLE_LIST = sorted ( glob.glob (DIR + "/" + DONOR + "/" + TISSUE + "/*") ) 
            SAMPLE_LIST = [i.split("/")[-1].split(".")[0] for i in SAMPLE_LIST]

            for SAMPLE in SAMPLE_LIST:       # PD42566b_lo00_A7.txt
                INPUT_TSV = DIR + "/" + DONOR + "/" + TISSUE + "/" + SAMPLE + ".txt"

                if os.path.exists(INPUT_TSV) == False:
                    continue
                
                kwargs = {"INPUT_TSV" : INPUT_TSV,  "MODE" : "Both",  "NUM_CLONE_TRIAL_START" : 1, "NUM_CLONE_TRIAL_END" : 5, 
                                "TRIAL_NO" : 10, "DEPTH_CUTOFF" : 16,  "KMEANS_CLUSTERNO" : 7, "MIN_CLUSTER_SIZE" : 5,  "MAKEONE_STRICT" :  3,
                                "RANDOM_PICK" : 0, "AXIS_RATIO":0, "PARENT_RATIO": 0, "NUM_PARENT" : 0,  "FP_RATIO":0,  "FP_USEALL" : "False", "FP_PRIOR" : 0.01,
                                "RANDOM_SEED" : 0, "SAMPLENAME" : "", "BENCHMARK_NO" : 0, 
                                "NPVAF_DIR" : "", "SIMPLE_KMEANS_DIR" : "", "CLEMENT_DIR" : "", "SCICLONE_DIR" : "", "PYCLONEVI_DIR" : "",  "COMBINED_OUTPUT_DIR" : "",
                                "SCORING" : False,  "MAXIMUM_NUM_PARENT" : 1, "VERBOSE" : 1,
                                 "IMAGE_FORMAT" : "jpg", "VISUALIZATION" : "True" }

                NUMBER_LINE = int(out ("wc -l  " + INPUT_TSV).split(" ")[0]) 

                if (NUMBER_LINE < 350) & ("adrenal_gland" not in TISSUE):       # 350줄이 안 되는 파일은 넘어간다
                    continue 

                kwargs ["RANDOM_PICK"] = -1    # Select all

                INPUT_FILETYPE, NUM_BLOCK = filetype.main (INPUT_TSV)
                kwargs["NUM_BLOCK_INPUT"] = kwargs["NUM_BLOCK"] = NUM_BLOCK
                SAMPLENAME = kwargs["SAMPLENAME"] = re.split(r'[_ .]', INPUT_TSV.split("/")[-1])[0]
                
        

                kwargs["NPVAF_DIR"] = "/data/project/Alzheimer/CLEMENT/02.npvaf/3.BioData/Moore_1D/" + TISSUE + "/" + DONOR + "-" + SAMPLENAME
                kwargs["COMBINED_OUTPUT_DIR"] = "/data/project/Alzheimer/CLEMENT/03.combinedoutput/3.BioData/Moore_1D/" + TISSUE + "/" + DONOR + "-" + SAMPLENAME
                kwargs["SIMPLE_KMEANS_DIR"] = "/data/project/Alzheimer/YSscript/cle/data/SIMPLE_KMEANS/3.BioData/Moore_1D/" + TISSUE + "/" + DONOR + "-" + SAMPLENAME
                kwargs["CLEMENT_DIR"] = "/data/project/Alzheimer/YSscript/cle/data/CLEMENT/3.BioData/Moore_1D/" + TISSUE + "/" + DONOR + "-" + SAMPLENAME
                kwargs["PYCLONEVI_DIR"] = "/data/project/Alzheimer/YSscript/cle/data/pyclone-vi/3.BioData/Moore_1D/" + TISSUE  + "/" + DONOR + "-" + SAMPLENAME
                kwargs["SCICLONE_DIR"] = "/data/project/Alzheimer/YSscript/cle/data/sciclone/3.BioData/Moore_1D/" + TISSUE + "/" + DONOR + "-" + SAMPLENAME
                kwargs["QUANTUMCLONE_DIR"] = "/data/project/Alzheimer/YSscript/cle/data/quantumclone/3.BioData/Moore_1D/" + TISSUE  + "/" + DONOR + "-" + SAMPLENAME

                # for SUBDIR in [ kwargs["NPVAF_DIR"],  kwargs["COMBINED_OUTPUT_DIR"] , kwargs["SIMPLE_KMEANS_DIR"],  kwargs["CLEMENT_DIR"], kwargs["PYCLONEVI_DIR"], kwargs["SCICLONE_DIR"],  kwargs["QUANTUMCLONE_DIR"] ]:
                #     if os.path.exists(SUBDIR) == True:
                #         os.system("rm -rf " + SUBDIR)
                #     if os.path.exists(SUBDIR) == False:
                #         os.system("mkdir -p " + SUBDIR)

                print ("# n = {}\tNUMBER_LINE = {} {}".format (n,  NUMBER_LINE, SAMPLENAME))

                #1. EM 돌리기
                logPath = "/data/project/Alzheimer/YSscript/cle/log/3.BioData/Moore_1D/" + TISSUE + "/" + DONOR + "-" + SAMPLENAME
                logPath_copy = kwargs["COMBINED_OUTPUT_DIR"] 
                os.system ("rm -rf " + logPath)
                os.system ("mkdir -p " + logPath)
                hold_j = TISSUE + "_" + DONOR + "_" + SAMPLE
                command = " ".join ( [ "qsub -pe smp 1", "-e", logPath, "-o", logPath, 
                                                    "-N", TISSUE + "_" + DONOR + "_" + SAMPLE, 
                                                    SCRIPT_DIR  + "/2.CellData_pipe1_CLEMENT_bm.sh",  
                                                    str(SCRIPT_DIR), str(INPUT_TSV),  str(kwargs["MODE"]),  str(kwargs["NUM_CLONE_TRIAL_START"]),  str(kwargs["NUM_CLONE_TRIAL_END"]), 
                                                    str(kwargs["RANDOM_PICK"]), str(kwargs["AXIS_RATIO"]),  str(kwargs["PARENT_RATIO"]),  str(kwargs["NUM_PARENT"]),  str(kwargs["FP_RATIO"]),  str(kwargs["FP_USEALL"]), str(kwargs ["FP_PRIOR"]),
                                                    str(kwargs["TRIAL_NO"]), str(kwargs["DEPTH_CUTOFF"]),  str(kwargs["MIN_CLUSTER_SIZE"]),  str(kwargs["VERBOSE"]),
                                                    str(kwargs["KMEANS_CLUSTERNO"]),  str(kwargs["RANDOM_SEED"]), str(kwargs["SAMPLENAME"]), str(kwargs["BENCHMARK_NO"]), 
                                                    str(kwargs["NPVAF_DIR"]), str(kwargs["SIMPLE_KMEANS_DIR"]), str(kwargs["CLEMENT_DIR"]), str(kwargs["SCICLONE_DIR"]), str(kwargs["PYCLONEVI_DIR"]) , str(kwargs["QUANTUMCLONE_DIR"]),  str(kwargs["COMBINED_OUTPUT_DIR"]), 
                                                    str(kwargs["SCORING"]), str(kwargs["MAKEONE_STRICT"]), str(kwargs["MAXIMUM_NUM_PARENT"]), str (kwargs["IMAGE_FORMAT"] ) , str (kwargs["VISUALIZATION"]),
                                                    str(logPath),  str(logPath_copy)     ] )
                n = n  + 1
                os.system (command)


                # #2. MatrixFormation + SigProfiler
                # toollist = ["CLEMENT_decision", "pyclonevi", "sciclone", "quantumclone"]
                # for tool in toollist:
                #     logPath = "/data/project/Alzheimer/YSscript/cle/log/3.BioData/Moore_1D/" + TISSUE + "/" + DONOR + "-" + SAMPLENAME + "/" + tool
                #     os.system ("rm -rf " + logPath)
                #     os.system ("mkdir -p " + logPath)

                #     if tool == "CLEMENT_decision":
                #         kwargs["OUTPUT_DIR"] = kwargs["CLEMENT_DIR"] + "/SigProfiler"
                #         kwargs["DECISION_MEMBERSHIP_PATH"] = kwargs["COMBINED_OUTPUT_DIR"] + "/result/CLEMENT_decision.membership.txt" 
                #     elif tool == "pyclonevi":
                #         kwargs["OUTPUT_DIR"] = kwargs["PYCLONEVI_DIR"] + "/SigProfiler"
                #         kwargs["DECISION_MEMBERSHIP_PATH"] = kwargs["COMBINED_OUTPUT_DIR"] + "/result/pyclonevi.membership.txt" 
                #     elif tool == "sciclone":
                #         kwargs["OUTPUT_DIR"] = kwargs["SCICLONE_DIR"] + "/SigProfiler"
                #         kwargs["DECISION_MEMBERSHIP_PATH"] = kwargs["COMBINED_OUTPUT_DIR"] + "/result/sciclone.membership.txt" 
                #     elif tool == "quantumclone":
                #         kwargs["OUTPUT_DIR"] = kwargs["QUANTUMCLONE_DIR"] + "/SigProfiler"
                #         kwargs["DECISION_MEMBERSHIP_PATH"] = kwargs["COMBINED_OUTPUT_DIR"] + "/result/quantumclone.membership.txt" 

                #     os.system ("rm -rf " + kwargs["OUTPUT_DIR"])
                #     os.system ("rm -rf " + kwargs["OUTPUT_DIR"] + "MatrixGenerator")
                #     os.system ("mkdir -p " + kwargs["OUTPUT_DIR"])
                #     command = " ".join(["qsub -pe smp 4", "-e", logPath, "-o", logPath, 
                #                     "-N", "Sig_" + TISSUE + "_" + DONOR + "_" + SAMPLE + "_" + tool,
                #                     "-hold_jid",  str( hold_j ), 
                #                     SCRIPT_DIR + "/3.BioData_pipe1_Signature.sh",
                #                     "--SCRIPT_DIR", str(SCRIPT_DIR), 
                #                     "--DECISION_MEMBERSHIP_PATH", kwargs["DECISION_MEMBERSHIP_PATH"] , 
                #                     "--PLOT_OUTPUT_PATH", kwargs["COMBINED_OUTPUT_DIR"] + "/AssignmentPlot_" + tool + ".pdf",
                #                     "--NPVAF_PATH", kwargs["NPVAF_DIR"] + "/npvaf.txt", 
                #                     "--DONOR", DONOR,
                #                     "--TISSUE", TISSUE,
                #                     "--OUTPUT_DIR", str( kwargs["OUTPUT_DIR"] ) ])

                #     # if os.path.exists( kwargs["OUTPUT_DIR"] + "/output/Assignment/Assignment_Solution/Activities/Assignment_Solution_Activities.txt" ) == False:
                #     #     print (kwargs["OUTPUT_DIR"])
                #     os.system (command)


########################################################################################################################################################################

def adrenal_gland_continuous ( SCRIPT_DIR  ):
    import glob

    DIR = "/data/project/Alzheimer/CLEMENT/01.INPUT_TSV/3.BioData/Moore_1D/3.woMosaic_ver2"        # AG : 3.woMosaic_ver2
    DONOR_LIST = sorted(glob.glob(DIR + "/*"))
    DONOR_LIST = [i.split("/")[-1] for i in DONOR_LIST]
    n = 0

    for DONOR in DONOR_LIST:               # PD28690
        AG_TISSUE_LIST = sorted(glob.glob(DIR + "/" + DONOR + "/adrenal_gland_zona/*"  ))
        AG_TISSUE_LIST = [i.split("/")[-1].split("_input.txt")[0] for i in AG_TISSUE_LIST]

        for AG_TISSUE in AG_TISSUE_LIST:            
            if "reticularis" in AG_TISSUE:
                continue

            INPUT_TSV = DIR + "/" + DONOR + "/adrenal_gland_zona/" + AG_TISSUE + "_input.txt"

            kwargs = {"INPUT_TSV" : INPUT_TSV,  "MODE" : "Both",  "NUM_CLONE_TRIAL_START" : 1, "NUM_CLONE_TRIAL_END" : 5, 
                    "TRIAL_NO" : 10, "DEPTH_CUTOFF" : 10,  "KMEANS_CLUSTERNO" : 8, "MIN_CLUSTER_SIZE" : 10,  "MAKEONE_STRICT" :  3,
                    "RANDOM_PICK" : 0, "AXIS_RATIO":0.4, "PARENT_RATIO": 0, "NUM_PARENT" : 0,  "FP_RATIO":0,  "FP_USEALL" : "False", 
                    "RANDOM_SEED" : 0, "SAMPLENAME" : "", "BENCHMARK_NO" : 0, 
                    "NPVAF_DIR" : "", "SIMPLE_KMEANS_DIR" : "", "CLEMENT_DIR" : "", "SCICLONE_DIR" : "", "PYCLONEVI_DIR" : "",  "COMBINED_OUTPUT_DIR" : "",
                    "SCORING" : False,  "MAXIMUM_NUM_PARENT" : 1, "VERBOSE" : 1 ,
                    "IMAGE_FORMAT" : "jpg", "VISUALIZATION" : "True" }
            kwargs["RANDOM_PICK"] = -1


            INPUT_FILETYPE, NUM_BLOCK = filetype.main (INPUT_TSV)
            kwargs["NUM_BLOCK_INPUT"] = kwargs["NUM_BLOCK"] = NUM_BLOCK
            SAMPLENAME = kwargs["SAMPLENAME"] = AG_TISSUE
            TISSUE = "adrenal_gland_zona"

            kwargs["NPVAF_DIR"] = "/data/project/Alzheimer/CLEMENT/02.npvaf/3.BioData/Moore_1D_AG/" +  AG_TISSUE
            kwargs["COMBINED_OUTPUT_DIR"] = "/data/project/Alzheimer/CLEMENT/03.combinedoutput/3.BioData/Moore_1D_AG/" +  AG_TISSUE
            kwargs["SIMPLE_KMEANS_DIR"] = "/data/project/Alzheimer/YSscript/cle/data/SIMPLE_KMEANS/3.BioData/Moore_1D_AG/" +  AG_TISSUE
            kwargs["CLEMENT_DIR"] = "/data/project/Alzheimer/YSscript/cle/data/CLEMENT/3.BioData/Moore_1D_AG/" +  AG_TISSUE
            kwargs["PYCLONEVI_DIR"] = "/data/project/Alzheimer/YSscript/cle/data/pyclone-vi/3.BioData/Moore_1D_AG/" + TISSUE  + "/" + DONOR + "-" + AG_TISSUE
            kwargs["SCICLONE_DIR"] = "/data/project/Alzheimer/YSscript/cle/data/sciclone/3.BioData/Moore_1D_AG/" +  AG_TISSUE
            kwargs["QUANTUMCLONE_DIR"] = "/data/project/Alzheimer/YSscript/cle/data/quantumclone/3.BioData/Moore_1D_AG/" +  AG_TISSUE
            
            NUMBER_LINE = int(out ("wc -l  " + INPUT_TSV).split(" ")[0]) 

            print ("# n = {}\t{} (wc - l = {})".format ( n, AG_TISSUE, NUMBER_LINE ) )

            os.system ("mkdir -p " + kwargs["NPVAF_DIR"])   # 출력 디렉토리 만들기
            os.system ("mkdir -p " + kwargs["SIMPLE_KMEANS_DIR"])   # 출력 디렉토리 만들기
            os.system ("mkdir -p " + kwargs["CLEMENT_DIR"])   # 출력 디렉토리 만들기
            os.system ("mkdir -p " + kwargs["SCICLONE_DIR"])   # 출력 디렉토리 만들기
            os.system ("mkdir -p " + kwargs["PYCLONEVI_DIR"])   # 출력 디렉토리 만들기
            os.system ("mkdir -p " + kwargs["QUANTUMCLONE_DIR"])   # 출력 디렉토리 만들기
            os.system ("mkdir -p " + kwargs["COMBINED_OUTPUT_DIR"])   # 출력 디렉토리 만들기

            #1. EM 돌리기
            logPath = "/data/project/Alzheimer/YSscript/cle/log/3.BioData/Moore_1D_AG/" +  AG_TISSUE
            os.system ("rm -rf " + logPath)
            os.system ("mkdir -p " + logPath)
            logPath_copy = kwargs["COMBINED_OUTPUT_DIR"] 

            hold_j = "Moore_1D_AG_" + AG_TISSUE
            command = " ".join ( [ "qsub -pe smp 1", "-e", logPath, "-o", logPath, 
                                                "-N", "Moore_1D_AG_" + AG_TISSUE, 
                                                SCRIPT_DIR  + "/2.CellData_pipe1_CLEMENT_bm.sh",  
                                                str(SCRIPT_DIR), str(INPUT_TSV),  str(kwargs["MODE"]),  str(kwargs["NUM_CLONE_TRIAL_START"]),  str(kwargs["NUM_CLONE_TRIAL_END"]), 
                                                str(kwargs["RANDOM_PICK"]), str(kwargs["AXIS_RATIO"]),  str(kwargs["PARENT_RATIO"]),  str(kwargs["NUM_PARENT"]),  str(kwargs["FP_RATIO"]),  str(kwargs["FP_USEALL"]), str(kwargs ["FP_PRIOR"]),
                                                str(kwargs["TRIAL_NO"]), str(kwargs["DEPTH_CUTOFF"]),  str(kwargs["MIN_CLUSTER_SIZE"]),  str(kwargs["VERBOSE"]),
                                                str(kwargs["KMEANS_CLUSTERNO"]),  str(kwargs["RANDOM_SEED"]), str(kwargs["SAMPLENAME"]), str(kwargs["BENCHMARK_NO"]), 
                                                str(kwargs["NPVAF_DIR"]), str(kwargs["SIMPLE_KMEANS_DIR"]), str(kwargs["CLEMENT_DIR"]), str(kwargs["SCICLONE_DIR"]), str(kwargs["PYCLONEVI_DIR"]) , str(kwargs["QUANTUMCLONE_DIR"]),  str(kwargs["COMBINED_OUTPUT_DIR"]), 
                                                str(kwargs["SCORING"]), str(kwargs["MAKEONE_STRICT"]), str(kwargs["MAXIMUM_NUM_PARENT"]), str (kwargs["IMAGE_FORMAT"] ) , str (kwargs["VISUALIZATION"]),
                                                str(logPath),  str(logPath_copy)        ] )
            os.system (command)
            n = n + 1


            #2. MatrixFormation + SigProfiler
            toollist = ["CLEMENT_decision", "pyclonevi", "sciclone", "quantumclone"]
            for tool in toollist:
                logPath = "/data/project/Alzheimer/YSscript/cle/log/3.BioData/Moore_1D_AG/" + AG_TISSUE + "/" + tool
                os.system ("rm -rf " + logPath)
                os.system ("mkdir -p " + logPath)

                if tool == "CLEMENT_decision":
                    kwargs["OUTPUT_DIR"] = kwargs["CLEMENT_DIR"] + "/SigProfiler"
                    kwargs["DECISION_MEMBERSHIP_PATH"] = kwargs["COMBINED_OUTPUT_DIR"] + "/result/CLEMENT_decision.membership.txt" 
                elif tool == "pyclonevi":
                    kwargs["OUTPUT_DIR"] = kwargs["PYCLONEVI_DIR"] + "/SigProfiler"
                    kwargs["DECISION_MEMBERSHIP_PATH"] = kwargs["COMBINED_OUTPUT_DIR"] + "/result/pyclonevi.membership.txt" 
                elif tool == "sciclone":
                    kwargs["OUTPUT_DIR"] = kwargs["SCICLONE_DIR"] + "/SigProfiler"
                    kwargs["DECISION_MEMBERSHIP_PATH"] = kwargs["COMBINED_OUTPUT_DIR"] + "/result/sciclone.membership.txt" 
                elif tool == "quantumclone":
                    kwargs["OUTPUT_DIR"] = kwargs["QUANTUMCLONE_DIR"] + "/SigProfiler"
                    kwargs["DECISION_MEMBERSHIP_PATH"] = kwargs["COMBINED_OUTPUT_DIR"] + "/result/quantumclone.membership.txt" 

                os.system ("rm -rf " + kwargs["OUTPUT_DIR"])
                os.system ("rm -rf " + kwargs["OUTPUT_DIR"] + "MatrixGenerator")
                os.system ("mkdir -p " + kwargs["OUTPUT_DIR"])

                command = " ".join(["qsub -pe smp 2", "-e", logPath, "-o", logPath, 
                                "-N", "Sig_Moore_1D_AG_" + AG_TISSUE,
                                #"-hold_jid",  str( hold_j ), 
                                SCRIPT_DIR + "/3.BioData_pipe1_Signature.sh",
                                "--SCRIPT_DIR", str(SCRIPT_DIR), 
                                "--DECISION_MEMBERSHIP_PATH",  kwargs["DECISION_MEMBERSHIP_PATH"], 
                                "--PLOT_OUTPUT_PATH", kwargs["COMBINED_OUTPUT_DIR"] + "/AssignmentPlot_" + tool + ".pdf",
                                "--NPVAF_PATH", kwargs["NPVAF_DIR"] + "/npvaf.txt", 
                                "--DONOR", DONOR,
                                "--TISSUE", AG_TISSUE,
                                "--OUTPUT_DIR", str( kwargs["OUTPUT_DIR"] ) ])
                os.system (command)




if __name__ == "__main__":
    SCRIPT_DIR = os.path.dirname(__file__)
    print (SCRIPT_DIR, "\n")

    Other_tissues (SCRIPT_DIR)
    print ("\n\n")
    #adrenal_gland_continuous ( SCRIPT_DIR )
