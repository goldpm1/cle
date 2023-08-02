import  filetype, argparse, os
import numpy as np
import pandas as pd

# python3 benchmark_master.py --INPUT_TSV "/data/project/Alzheimer/EM_cluster/EM_input/MRS_2_sample/M1-5_M1-6_input.txt" --FP_RATIO 0.1 --AXIS_RATIO 0.1  --NUM_PARENT 2 --VERBOSE 1 --NUM_CLONE_TRIAL_START 2 --NUM_CLONE_TRIAL_END 7 --TRIAL_NO 6

kwargs = {}

parser = argparse.ArgumentParser(description='The below is usage direction.')
parser.add_argument('--INPUT_TSV', type = str, default="/data/project/Alzheimer/EM_cluster/EM_input/MRS_2_sample/M1-3_M1-8_input.txt", help = "Input data whether TSV or VCF. The tool automatically detects the number of samples")
parser.add_argument('--MODE', type = str, choices = ["Hard", "Soft", "Both"], default="Hard")
parser.add_argument('--NUM_CLONE_TRIAL_START', type = int, default=2, help = "Minimum number of expected cluster_hards (initation of K)")
parser.add_argument('--NUM_CLONE_TRIAL_END', type = int, default=7, choices=range(1, 11), help = "Maximum number of expected cluster_hards (termination of K)")
parser.add_argument('--NUM_CLONE_TRIAL_FORCE', type = int, default=4, help = "Designate the number of expected cluster_hards by force")
parser.add_argument('--NPVAF_DIR', default = None, help = "Directory where selected datasets are")
parser.add_argument('--MYEM_DIR', default = None, help = "Directory where input and output of MY EM deposits")
parser.add_argument('--SCICLONE_DIR', default = None, help = "Directory where input and output of SCICLONE deposits")
parser.add_argument('--PYCLONE_DIR', default = None, help = "Directory where input and output of PYCLONE deposits")
parser.add_argument('--PYCLONEVI_DIR', default = None, help = "Directory where input and output of PYCLONEVI deposits")
parser.add_argument('--COMBINED_OUTPUT_DIR', default = None, help = "Directory where input and output of MYEM, PYCLONEVI, SCICLONE deposits")
parser.add_argument('--RANDOM_PICK', type = int, default=500, help = "The number of mutations that are randomly selected in each trials")
parser.add_argument('--AXIS_RATIO', default=0,  type = float, help = "The fraction of the mutations not shared at least one sample")
parser.add_argument('--PARENT_RATIO', default=0,  type = float, help = "The fraction of parent clone mutations. If this values is designated, do not set NUM_PARENT")
parser.add_argument('--NUM_PARENT',  default=0, type = int, help = "The fraction of parent clones being inserted by large order. If this values is designated, do not set PARENT_RATIO")
parser.add_argument('--FP_RATIO', default=0,  type = float, help = "The fraction of false positive mutations regardless of all samples")
parser.add_argument('--FP_2D', default="False", choices = ["True", "False"], help = "True : extract ALL FPs,   False : Do not extact FPs")
parser.add_argument('--TRIAL_NO', default=3, type = int, choices=range(1, 21), help = "Trial number in each candidate cluster_hard number. DO NOT recommend over 15")
parser.add_argument('--DEPTH_CUTOFF', default=100, type = int, help = "The mutation of which depth below this values is abandoned")
parser.add_argument('--MIN_CLUSTER_SIZE', type = int, default=5)
parser.add_argument('--VERBOSE', type = int, choices = [0, 1, 2], default=1)
parser.add_argument('--KMEANS_CLUSTERNO',  type = int , default=15, choices = range(8,20), help = "Number of initial K-means cluster_hard")
parser.add_argument('--RANDOM_SEED', type = int, default=1, help = "random_seed for regular random sampling")
parser.add_argument('--BENCHMARK_NO', type = int, default=10, help = "the number of random_seed sampling")



args = parser.parse_args()
kwargs["INPUT_TSV"] = args.INPUT_TSV

INPUT_TSV = kwargs["INPUT_TSV"]
INPUT_FILETYPE, NUM_BLOCK = filetype.main (INPUT_TSV)
kwargs["NUM_BLOCK_INPUT"], kwargs["NUM_BLOCK"] = NUM_BLOCK, NUM_BLOCK
SAMPLENAME = INPUT_TSV.split("/")[-1].split(".")[0]     # 'M1-5_M1-8_input'


kwargs["MODE"] = args.MODE
kwargs["NUM_CLONE_TRIAL_START"], kwargs["NUM_CLONE_TRIAL_END"], kwargs["NUM_CLONE_TRIAL_FORCE"] = args.NUM_CLONE_TRIAL_START, args.NUM_CLONE_TRIAL_END, args.NUM_CLONE_TRIAL_FORCE
kwargs["RANDOM_PICK"] = int(args.RANDOM_PICK)
kwargs["AXIS_RATIO"] = float(args.AXIS_RATIO)
kwargs["PARENT_RATIO"] = float(args.PARENT_RATIO)
kwargs["NUM_PARENT"] = int(args.NUM_PARENT)
kwargs["FP_RATIO"] = float(args.FP_RATIO)
kwargs["FP_2D"] = args.FP_2D
kwargs["TRIAL_NO"] = int(args.TRIAL_NO)
kwargs["DEPTH_CUTOFF"] = int(args.DEPTH_CUTOFF)
kwargs["VERBOSE"] = int(args.VERBOSE)
kwargs["MIN_CLUSTER_SIZE"] = int(args.MIN_CLUSTER_SIZE)
kwargs["KMEANS_CLUSTERNO"] = args.KMEANS_CLUSTERNO
kwargs["RANDOM_SEED"] = int(args.RANDOM_SEED)
kwargs["BENCHMARK_NO"] = int(args.BENCHMARK_NO)
kwargs["SAMPLENAME"] = SAMPLENAME
kwargs["SCORING"] = True



#1.

hold_j = []
for i in range (1, kwargs["BENCHMARK_NO"] + 1):
    if args.NPVAF_DIR == None:
        kwargs["NPVAF_DIR"] = "/data/project/Alzheimer/YSscript/EM_MRS/data/npvaf/MRS/" + SAMPLENAME + "/" + str(kwargs["RANDOM_PICK"]) + "_" + str(kwargs["NUM_PARENT"]) + "_" + str(kwargs["FP_RATIO"]) + "_" + str(kwargs["AXIS_RATIO"]) + "/"  + str(i)
    else:
        kwargs["NPVAF_DIR"] = args.NPVAF_DIR
    if args.MYEM_DIR == None:
        kwargs["MYEM_DIR"] = "/data/project/Alzheimer/YSscript/EM_MRS/data/MyEM/MRS/" + SAMPLENAME + "/" + str(kwargs["RANDOM_PICK"]) + "_" + str(kwargs["NUM_PARENT"]) + "_" + str(kwargs["FP_RATIO"]) + "_" + str(kwargs["AXIS_RATIO"]) + "/" + str(i)
    else:
        kwargs["MYEM_DIR"] = args.MYEM_DIR
    if args.SCICLONE_DIR == None:
        kwargs["SCICLONE_DIR"] = "/data/project/Alzheimer/YSscript/EM_MRS/data/sciclone/MRS/" + SAMPLENAME + "/" + str(kwargs["RANDOM_PICK"]) + "_" + str(kwargs["NUM_PARENT"]) + "_" + str(kwargs["FP_RATIO"]) + "_" + str(kwargs["AXIS_RATIO"]) + "/" + str(i)
    else:
        kwargs["SCICLONE_DIR"] = args.SCICLONE_DIR
    if args.PYCLONE_DIR == None:
        kwargs["PYCLONE_DIR"] = "/data/project/Alzheimer/YSscript/EM_MRS/data/pyclone/MRS/" + SAMPLENAME + "/" + str(kwargs["RANDOM_PICK"]) + "_" + str(kwargs["NUM_PARENT"]) + "_" + str(kwargs["FP_RATIO"]) + "_" + str(kwargs["AXIS_RATIO"]) + "/" + str(i)
    else:
        kwargs["PYCLONE_DIR"] = args.PYCLONE_DIR
    if args.PYCLONEVI_DIR == None:
        kwargs["PYCLONEVI_DIR"] = "/data/project/Alzheimer/YSscript/EM_MRS/data/pyclone-vi/MRS/" + SAMPLENAME + "/" + str(kwargs["RANDOM_PICK"]) + "_" + str(kwargs["NUM_PARENT"]) + "_" + str(kwargs["FP_RATIO"]) + "_" + str(kwargs["AXIS_RATIO"]) + "/" + str(i)
    else:
        kwargs["PYCLONEVI_DIR"] = args.PYCLONEVI_DIR
    if args.COMBINED_OUTPUT_DIR == None:
        kwargs["COMBINED_OUTPUT_DIR"] = "/data/project/Alzheimer/YSscript/EM_MRS/data/combinedoutput/MRS/" + SAMPLENAME + "/" + str(kwargs["RANDOM_PICK"]) + "_" + str(kwargs["NUM_PARENT"]) + "_" + str(kwargs["FP_RATIO"]) + "_" + str(kwargs["AXIS_RATIO"]) + "/" + str(i)
    else:
        kwargs["COMBINED_OUTPUT_DIR"] = args.COMBINED_OUTPUT_DIR



    logPath = "/data/project/Alzheimer/YSscript/EM_MRS/log/" + SAMPLENAME + "/" + str(kwargs["RANDOM_PICK"]) + "_" + str(kwargs["NUM_PARENT"]) + "_" + str(kwargs["FP_RATIO"]) + "_" + str(kwargs["AXIS_RATIO"]) + "/"  + str(i)
    os.system ("rm -rf " + logPath)
    os.system ("mkdir -p " + logPath)

    hold_j.append (  SAMPLENAME + "_" + str(kwargs["RANDOM_PICK"]) + "_" + str(kwargs["NUM_PARENT"]) + "_" + str(kwargs["FP_RATIO"]) + "_" + str(kwargs["AXIS_RATIO"]) + "_"  + str(i)  )
    command = " ".join ( [ "qsub -pe smp 1", "-e", logPath, "-o", logPath, "-N", SAMPLENAME + "_" + str(kwargs["RANDOM_PICK"]) + "_" + str(kwargs["NUM_PARENT"]) + "_" + str(kwargs["FP_RATIO"]) + "_" + str(kwargs["AXIS_RATIO"]) + "_"  + str(i), 
                                        "benchmark_pipe1.sh",  
                                        str(INPUT_TSV),  str(kwargs["MODE"]),  str(kwargs["NUM_CLONE_TRIAL_START"]),  str(kwargs["NUM_CLONE_TRIAL_END"]),  str(kwargs["NUM_CLONE_TRIAL_FORCE"]),
                                        str(kwargs["RANDOM_PICK"]), str(kwargs["AXIS_RATIO"]),  str(kwargs["PARENT_RATIO"]),  str(kwargs["NUM_PARENT"]),  str(kwargs["FP_RATIO"]),  str(kwargs["FP_2D"]),
                                        str(kwargs["TRIAL_NO"]), str(kwargs["DEPTH_CUTOFF"]),  str(kwargs["MIN_CLUSTER_SIZE"]),  str(kwargs["VERBOSE"]),
                                        str(kwargs["KMEANS_CLUSTERNO"]),  str(i), str(kwargs["SAMPLENAME"]), str(kwargs["BENCHMARK_NO"]), 
                                        str(kwargs["NPVAF_DIR"]), str(kwargs["MYEM_DIR"]), str(kwargs["SCICLONE_DIR"]), str(kwargs["PYCLONE_DIR"]), str(kwargs["PYCLONEVI_DIR"]) ,str(kwargs["COMBINED_OUTPUT_DIR"]) ,
                                        str(kwargs["MODE"]), str(kwargs["SCORING"])   ] )
    os.system (command)



#2.

logPath = "/data/project/Alzheimer/YSscript/EM_MRS/log/MRS/" + SAMPLENAME + "/" + str(kwargs["RANDOM_PICK"]) + "_" + str(kwargs["NUM_PARENT"]) + "_" + str(kwargs["FP_RATIO"]) + "_" + str(kwargs["AXIS_RATIO"]) + "/visualization"
os.system ("rm -rf " + logPath)
os.system ("mkdir -p " + logPath)

#print (",".join (hold_j))

command = " ".join ( [ "qsub -pe smp 1", "-e", logPath, "-o", logPath, "-N", SAMPLENAME + "_" + str(kwargs["RANDOM_PICK"]) + "_" + str(kwargs["NUM_PARENT"]) + "_" + str(kwargs["FP_RATIO"]) + "_" + str(kwargs["AXIS_RATIO"]) + "_visualization" ,
                                     "-hold_jid",  str( ",".join (hold_j)) , "benchmark_pipe2.sh",  
                                    str(INPUT_TSV),  str(kwargs["MODE"]),  str(kwargs["NUM_CLONE_TRIAL_START"]),  str(kwargs["NUM_CLONE_TRIAL_END"]),  str(kwargs["NUM_CLONE_TRIAL_FORCE"]),
                                    str(kwargs["RANDOM_PICK"]), str(kwargs["AXIS_RATIO"]),  str(kwargs["PARENT_RATIO"]),  str(kwargs["NUM_PARENT"]),  str(kwargs["FP_RATIO"]),  str(kwargs["FP_2D"]),
                                    str(kwargs["TRIAL_NO"]), str(kwargs["DEPTH_CUTOFF"]),  str(kwargs["MIN_CLUSTER_SIZE"]),  str(kwargs["VERBOSE"]),
                                    str(kwargs["KMEANS_CLUSTERNO"]),  str(i), str(kwargs["SAMPLENAME"]), str(kwargs["BENCHMARK_NO"]), 
                                    str(kwargs["NPVAF_DIR"]), str(kwargs["MYEM_DIR"]), str(kwargs["SCICLONE_DIR"]), str(kwargs["PYCLONE_DIR"]), str(kwargs["PYCLONEVI_DIR"]) ,str(kwargs["COMBINED_OUTPUT_DIR"]),
                                    str(kwargs["MODE"]), str(kwargs["SCORING"])    ] )
os.system (command)