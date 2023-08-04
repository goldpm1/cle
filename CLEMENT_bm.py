import os, re, subprocess, sys, datetime, time, io, contextlib, argparse
import numpy as np
import pandas as pd

print ( "Package directory : {}\n".format (  os.path.dirname(__file__) ) )
SCRIPT_DIR = os.path.dirname(__file__)
if os.path.dirname(__file__) not in sys.path:
    sys.path.append  ( os.path.dirname(__file__) )
    
import EMhard, EMsoft, Estep, Mstep, Bunch, miscellaneous, datapreparation, phylogeny, visualizationsingle, visualizationpair, visualizationsinglesoft, filetype, result, scoring, simplekmeans, pyclonevisim, sciclonesim, quantumclonesim

pd.options.mode.chained_assignment = None

kwargs = {}

parser = argparse.ArgumentParser(description='The below is usage direction.')
parser.add_argument('--INPUT_TSV', type=str, default="/data/project/Alzheimer/EM_cluster/EM_input/MRS_2_sample/M1-3_M1-8_input.txt",  help="Input data whether TSV or VCF. The tool automatically detects the number of samples")
parser.add_argument('--CLEMENT_DIR', default=None,   help="Directory where input and output of CLEMENT deposits")
parser.add_argument('--MODE', type=str, choices=["Hard", "Soft", "Both"], default="Hard")
parser.add_argument('--RANDOM_PICK', type=int, default=500,  help="The number of mutations that are randomly selected in each trials")
parser.add_argument('--KMEANS_CLUSTERNO',  type=int, default=8,  choices=range(5, 20), help="Number of initial K-means cluster_hard")
parser.add_argument('--NUM_CLONE_TRIAL_START', type=int, default=5,  help="Minimum number of expected cluster_hards (initation of K)")
parser.add_argument('--NUM_CLONE_TRIAL_END', type=int, default=6, choices=range(1, 11), help="Maximum number of expected cluster_hards (termination of K)")
parser.add_argument('--TRIAL_NO', default=3, type=int, choices=range(1, 21),  help="Trial number in each candidate cluster_hard number. DO NOT recommend over 15")
parser.add_argument('--MAXIMUM_NUM_PARENT',  default=1, type=int,  help="The maximum number of parents in the given samples.")
parser.add_argument('--MIN_CLUSTER_SIZE', type=int, default = 9)
parser.add_argument('--RANDOM_SEED', type=int, default=1,  help="random_seed for regular random sampling")
parser.add_argument('--MAKEONE_STRICT', type=int,  choices=[1, 2, 3], default = 1)
parser.add_argument('--TN_CONFIDENTIALITY', default=0.995, type=float, help="Confidentiality that negative being negative (TN). Recommendation : > 0.99. Default : 0.995")
parser.add_argument('--FONT_FAMILY', type=str, default="arial", help="Font family that displayed in the plots. Default : arial")
parser.add_argument('--IMAGE_FORMAT', type=str, default="jpg", choices = ["jpg", "pdf"], help="Image format that displayed in the plots. Default : jpg")
parser.add_argument('--VERBOSE', type=int, choices=[0, 1, 2, 3], default=2, help="0: Verbose, 3: Concise. Default : 2")


parser.add_argument('--NPVAF_DIR', default=None,  help="Directory where selected datasets are")
parser.add_argument('--COMBINED_OUTPUT_DIR', default=None,    help="Directory where input and output of CLEMENT, PYCLONEVI, SCICLONE deposits")
parser.add_argument('--SIMPLE_KMEANS_DIR', default=None,   help="Directory where input and output of SIMPLE_KMEANS deposits")
parser.add_argument('--SCICLONE_DIR', default=None,   help="Directory where input and output of SCICLONE deposits")
parser.add_argument('--PYCLONEVI_DIR', default=None,   help="Directory where input and output of PYCLONEVI deposits")
parser.add_argument('--QUANTUMCLONE_DIR', default=None,   help="Directory where input and output of QUANTUMCLONE deposits")

parser.add_argument('--AXIS_RATIO', default=0,  type=float,  help="The fraction of the mutations not shared at least one sample")
parser.add_argument('--PARENT_RATIO', default=0,  type=float,  help="The fraction of parent clone mutations. If this values is designated, do not set NUM_PARENT")
parser.add_argument('--NUM_PARENT',  default=0, type=int,  help="The fraction of parent clones being inserted by large order. If this values is designated, do not set PARENT_RATIO")
parser.add_argument('--FP_RATIO', default=0,  type=float, help="The fraction of false positive mutations regardless of all samples")
parser.add_argument('--FP_USEALL', default="False", choices=["True", "False"], help="True : extract ALL FPs,   False : Do not extact FPs")
parser.add_argument('--DEPTH_CUTOFF', default=100, type=int, help="The mutation of which depth below this values is abandoned")
parser.add_argument('--SCORING', type=str, choices=["True", "False"], default="True", help="True : comparing with the answer set,  False : just visualization")

args = parser.parse_args()


kwargs["INPUT_TSV"] = args.INPUT_TSV
INPUT_TSV = kwargs["INPUT_TSV"]
INPUT_FILETYPE, NUM_BLOCK = filetype.main(INPUT_TSV)
kwargs["NUM_BLOCK_INPUT"], kwargs["NUM_BLOCK"] = NUM_BLOCK, NUM_BLOCK
SAMPLENAME = re.split(r'[_ input]', INPUT_TSV.split("/")[-1])      # 'M1-5_M1-8_input' -> 'M1-5_M1-8'
kwargs["CLEMENT_DIR"] = args.CLEMENT_DIR
kwargs["MODE"] = args.MODE
kwargs["RANDOM_PICK"] = int(args.RANDOM_PICK)
kwargs["NUM_CLONE_TRIAL_START"], kwargs["NUM_CLONE_TRIAL_END"] = args.NUM_CLONE_TRIAL_START, args.NUM_CLONE_TRIAL_END
kwargs["TRIAL_NO"] = int(args.TRIAL_NO)
kwargs["MAXIMUM_NUM_PARENT"] = int(args.MAXIMUM_NUM_PARENT)
kwargs["KMEANS_CLUSTERNO"] = args.KMEANS_CLUSTERNO
kwargs["MIN_CLUSTER_SIZE"] = int(args.MIN_CLUSTER_SIZE)
kwargs["RANDOM_SEED"] = int(args.RANDOM_SEED)
kwargs["MAKEONE_STRICT"] = int(args.MAKEONE_STRICT)
kwargs["TN_CONFIDENTIALITY"] = float(args.TN_CONFIDENTIALITY)
kwargs["FONT_FAMILY"] = str(args.FONT_FAMILY)
kwargs["IMAGE_FORMAT"] = str(args.IMAGE_FORMAT)
kwargs["VERBOSE"] = int(args.VERBOSE)


kwargs["NPVAF_DIR"] = args.NPVAF_DIR
kwargs["SIMPLE_KMEANS_DIR"] = args.SIMPLE_KMEANS_DIR
kwargs["SCICLONE_DIR"] = args.SCICLONE_DIR
kwargs["PYCLONEVI_DIR"] = args.PYCLONEVI_DIR
kwargs["QUANTUMCLONE_DIR"] = args.QUANTUMCLONE_DIR
kwargs["COMBINED_OUTPUT_DIR"] = args.COMBINED_OUTPUT_DIR

kwargs["AXIS_RATIO"] = float(args.AXIS_RATIO)
kwargs["PARENT_RATIO"] = float(args.PARENT_RATIO)
kwargs["NUM_PARENT"] = int(args.NUM_PARENT)
kwargs["FP_RATIO"] = float(args.FP_RATIO)
kwargs["FP_USEALL"] = args.FP_USEALL
kwargs["DEPTH_CUTOFF"] = int(args.DEPTH_CUTOFF)
if args.SCORING in ["False", "false"]:
    kwargs["SCORING"] = False
elif args.SCORING in ["True", "true"]:
    kwargs["SCORING"] = True


kwargs["method"] = "gap+normal"
kwargs["adjustment"] = "half"
kwargs["STEP_NO"] = 30
kwargs["DECISION_STANDARD"]= 0.8 + (0.03 * kwargs ["NUM_BLOCK"])




print("\n\n\n\nNOW RUNNING IS STARTED  :  {}h:{}m:{}s\n".format(time.localtime().tm_hour, time.localtime().tm_min, round(time.localtime().tm_sec)))
print("NUMBER OF INPUT SAMPLES = {}\n\n\n".format(NUM_BLOCK))



print("============================== STEP #1.   DATA EXTRACTION FROM THE ANSWER SET  ==============================")

for DIR in [kwargs["NPVAF_DIR"], kwargs["SIMPLE_KMEANS_DIR"],  kwargs["SIMPLE_KMEANS_DIR"] + "/result", kwargs["SIMPLE_KMEANS_DIR"] + "/elbow", kwargs["SIMPLE_KMEANS_DIR"] + "/silhouette", kwargs["SIMPLE_KMEANS_DIR"] + "/gap", kwargs["CLEMENT_DIR"], kwargs["SCICLONE_DIR"], kwargs["PYCLONEVI_DIR"], kwargs["QUANTUMCLONE_DIR"], kwargs["COMBINED_OUTPUT_DIR"], kwargs["CLEMENT_DIR"] + "/trial",  kwargs["CLEMENT_DIR"] + "/Kmeans",   kwargs["CLEMENT_DIR"] + "/candidate",   kwargs["CLEMENT_DIR"]  + "/result", kwargs["COMBINED_OUTPUT_DIR"] + "/result"]:
    if os.path.exists(DIR) == True:
        os.system("rm -rf  " + DIR)
    if os.path.exists(DIR) == False:
        os.system("mkdir -p " + DIR)

output_logfile = open (kwargs["COMBINED_OUTPUT_DIR"] + "/0.commandline.txt", "w" , encoding="utf8" )
print ("python3 /data/project/Alzheimer/YSscript/cle/CLEMENT_bm.py  --INPUT_TSV {} --NUM_CLONE_TRIAL_START {} --NUM_CLONE_TRIAL_END {}  --RANDOM_PICK {} --AXIS_RATIO {} --PARENT_RATIO {} --NUM_PARENT {} --FP_RATIO {} --FP_USEALL {}   --DEPTH_CUTOFF {} --MIN_CLUSTER_SIZE {}   --KMEANS_CLUSTERNO {} --RANDOM_SEED {}  --NPVAF_DIR {} --SIMPLE_KMEANS_DIR {} --CLEMENT_DIR {} --SCICLONE_DIR {} --PYCLONEVI_DIR {} --QUANTUMCLONE_DIR {} --COMBINED_OUTPUT_DIR {}   --MODE {} --SCORING {} --MAKEONE_STRICT {} --MAXIMUM_NUM_PARENT {} --TRIAL_NO {}   --VERBOSE {}".format( kwargs["INPUT_TSV"], kwargs["NUM_CLONE_TRIAL_START"], kwargs["NUM_CLONE_TRIAL_END"],  kwargs["RANDOM_PICK"], kwargs["AXIS_RATIO"],  kwargs["PARENT_RATIO"],  kwargs["NUM_PARENT"],  kwargs["FP_RATIO"],  kwargs["FP_USEALL"],  kwargs["DEPTH_CUTOFF"],  kwargs["MIN_CLUSTER_SIZE"], kwargs["KMEANS_CLUSTERNO"],  kwargs["RANDOM_SEED"],  kwargs["NPVAF_DIR"], kwargs["SIMPLE_KMEANS_DIR"], kwargs["CLEMENT_DIR"], kwargs["SCICLONE_DIR"], kwargs["PYCLONEVI_DIR"], kwargs["QUANTUMCLONE_DIR"], kwargs["COMBINED_OUTPUT_DIR"], kwargs["MODE"], kwargs["SCORING"], kwargs["MAKEONE_STRICT"], kwargs["MAXIMUM_NUM_PARENT"], kwargs["TRIAL_NO"],   kwargs["VERBOSE"]  ),         file = output_logfile)
output_logfile.close()



inputdf, df, np_vaf, np_BQ, membership_answer, mixture_answer, mutation_id, kwargs = datapreparation.main( **kwargs)
NUM_MUTATION = kwargs["NUM_MUTATION"] = kwargs["RANDOM_PICK"]
membership_answer_numerical = np.zeros(kwargs["NUM_MUTATION"], dtype="int")
membership_answer_numerical_nofp_index = []

# membership_answer: ['V1', 'V2', 'FP', 'S0', 'V2', 'S0', 'V2', ....
# membership_answer_numerical : [0 1 2 3 1 3 1 ...

if type(inputdf) != type(False):
    # {0: 'FP', 1: 'V2', 2: 'S0', 3: 'V1'}
    kwargs["samplename_dict_NumToCharacter"] = {v: k for k, v in kwargs["samplename_dict_CharacterToNum"].items()}

    print ("samplename_dict_CharacterToNum = {}\nsamplename_dict_NumToCharacter = {}".format( kwargs["samplename_dict_CharacterToNum"], kwargs["samplename_dict_NumToCharacter"] ))

    with open(kwargs["COMBINED_OUTPUT_DIR"] + "/0.input_membership_numerical.txt", "w", encoding="utf8") as result_answer:
        for k in range(kwargs["NUM_MUTATION"]):
            membership_answer_numerical[k] = kwargs["samplename_dict_CharacterToNum"][membership_answer[k]]
            if (membership_answer[k] != "FP"):
                membership_answer_numerical_nofp_index.append(k)
            print(membership_answer_numerical[k], file=result_answer)

    with open(kwargs["COMBINED_OUTPUT_DIR"] + "/0.input_membership_letter.txt", "w", encoding="utf8") as result_answer:
        for k in range(kwargs["NUM_MUTATION"]):
            print(membership_answer[k], file=result_answer)

    with open(kwargs["COMBINED_OUTPUT_DIR"] + "/0.input_mixture.txt", "w", encoding="utf8") as result_answer:
        print(mixture_answer, file=result_answer)
        for i in range(NUM_BLOCK):
            sum_mixture = 0
            for j in range(mixture_answer.shape[1]):
                if "FP" not in kwargs["samplename_dict_CharacterToNum"]:
                    sum_mixture = sum_mixture + mixture_answer[i][j]
                else:
                    if kwargs["samplename_dict_CharacterToNum"]["FP"] != j:
                        sum_mixture = sum_mixture + mixture_answer[i][j]
            print("Sum of mixture (except FP) in sample {} : {}".format(i, sum_mixture), file=result_answer)
        print(kwargs["samplename_dict_NumToCharacter"],  file=result_answer)
        print(np.unique(membership_answer_numerical, return_counts=True)[1], file=result_answer)

    if kwargs["NUM_BLOCK"] == 1:
        x_median = miscellaneous.VAFdensitogram(np_vaf, "INPUT DATA", kwargs["COMBINED_OUTPUT_DIR"] + "/0.inputdata." + kwargs["IMAGE_FORMAT"], **kwargs)
        visualizationsingle.drawfigure_1d(membership_answer_numerical, "ANSWER_SET (n={})".format(kwargs["NUM_MUTATION"]), kwargs["COMBINED_OUTPUT_DIR"] + "/0.inputdata." + kwargs["IMAGE_FORMAT"], np_vaf, kwargs["samplename_dict_NumToCharacter"], False, -1, [], **kwargs )
    elif kwargs["NUM_BLOCK"] == 2:
        visualizationsingle.drawfigure_2d(membership_answer, "ANSWER_SET (n={})".format(kwargs["NUM_MUTATION"]), kwargs["COMBINED_OUTPUT_DIR"] + "/0.inputdata." + kwargs["IMAGE_FORMAT"], np_vaf, kwargs["samplename_dict_CharacterToNum"], False,  -1, "None", **kwargs)
    elif kwargs["NUM_BLOCK"] >= 3:
        visualizationsingle.drawfigure_2d(membership_answer, "ANSWER_SET (n={})".format(kwargs["NUM_MUTATION"]), kwargs["COMBINED_OUTPUT_DIR"] + "/0.inputdata." + kwargs["IMAGE_FORMAT"], np_vaf, kwargs["samplename_dict_CharacterToNum"], False,  -1, "SVD", **kwargs)
    subprocess.run(["cp " + kwargs["COMBINED_OUTPUT_DIR"] + "/0.inputdata." + kwargs["IMAGE_FORMAT"] + " "  +  kwargs["CLEMENT_DIR"] + "/candidate/0.inputdata." + kwargs["IMAGE_FORMAT"]], shell=True)
    subprocess.run(["cp " + kwargs["COMBINED_OUTPUT_DIR"] + "/0.inputdata." + kwargs["IMAGE_FORMAT"] + " "  +  kwargs["CLEMENT_DIR"] + "/trial/0.inputdata." + kwargs["IMAGE_FORMAT"]], shell=True)





START_TIME = datetime.datetime.now()


np_vaf = miscellaneous.np_vaf_extract(df)
mixture_kmeans, kwargs["KMEANS_CLUSTERNO"] = miscellaneous.initial_kmeans (np_vaf, kwargs["KMEANS_CLUSTERNO"], kwargs["CLEMENT_DIR"] + "/trial/0.inqitial_kmeans." + kwargs["IMAGE_FORMAT"])

cluster_hard = Bunch.Bunch2(**kwargs)
cluster_soft = Bunch.Bunch2(**kwargs)



print ("\n\n\n\n=============================================== STEP #2. SIMPLE_KMEANS  =======================================")

kwargs, simpleK = simplekmeans.clustering ( np_vaf, **kwargs )
simpleK = simplekmeans.scoring ( membership_answer, membership_answer_numerical, simpleK, **kwargs )
simplekmeans.visualization ( simpleK, np_vaf, **kwargs )
simplekmeans.save (simpleK, round((datetime.datetime.now() - START_TIME).total_seconds()) , **kwargs)

print ( "\t▸ Elbow method : k = {}, score = {}".format (simpleK.elbow_K, simpleK.elbow_K_score))
print ( "\t▸ Silhouette method : k = {}, score = {}".format (simpleK.silhouette_K, simpleK.silhouette_K_score ))
print ( "\t▸ Gap* method : k = {}, score = {}".format (simpleK.gap_K, simpleK.gap_K_score ))






######################################################### Step 4. Hard clustering ########################################################
print ("\n\n\n======================================== STEP #4.   EM HARD  ========================================")

cluster_hard = EMhard.main (df, np_vaf, np_BQ, mixture_kmeans, **kwargs)
subprocess.run (["cp -r " + kwargs["CLEMENT_DIR"] + "/candidate  " + kwargs["COMBINED_OUTPUT_DIR"]  ], shell = True)
subprocess.run (["cp -r " + kwargs["CLEMENT_DIR"] + "/trial  " + kwargs["COMBINED_OUTPUT_DIR"]  ], shell = True)





######################################################### Step 5. Hard -> Soft clustering ########################################################
if kwargs["MODE"] in ["Soft", "Both"]:
    print ("\n\n\n======================================== STEP #5.   EM SOFT  ========================================")

    for NUM_CLONE in range(kwargs["NUM_CLONE_TRIAL_START"], kwargs["NUM_CLONE_TRIAL_END"] + 1):
        print("\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n\nNUM_CLONE = {0}".format(NUM_CLONE))
        kwargs["NUM_CLONE"], kwargs["NUM_CLONE_NOMINAL"] = NUM_CLONE, NUM_CLONE
        kwargs["OPTION"] = "soft"

        if cluster_hard.likelihood_record[ NUM_CLONE ] !=  float("-inf"):
            print("\n\n\tSequential Soft clustering (TRIAL_NO = {}, STEP_NO = {})".format ( cluster_hard.trialindex_record[ NUM_CLONE ], cluster_hard.stepindex_record [ NUM_CLONE ] ))
            step_soft = Bunch.Bunch1(kwargs["NUM_MUTATION"] , NUM_BLOCK, NUM_CLONE, cluster_hard.stepindex_record [ NUM_CLONE ] + kwargs["STEP_NO"])
            step_soft.copy (cluster_hard, 0, NUM_CLONE)  # 0번 step에 cluster_hard를 복사한다

            for step_index in range(1, kwargs["STEP_NO"]):   # 0번은 채웠으니 1번부터 시작
                kwargs["STEP"], kwargs["TRIAL"] = step_index, cluster_hard.trialindex_record[ NUM_CLONE ]
                kwargs["STEP_TOTAL"] = step_index + cluster_hard.stepindex_record [ NUM_CLONE ] - 1
                
                print("\t\t#{}번째 step ( = TOTAL {}번째 step)".format(kwargs["STEP"], kwargs["STEP_TOTAL"]) )

                step_soft = Estep.main(df, np_vaf, np_BQ, step_soft, **kwargs)                   # 주어진 mixture 내에서 새 membership 정하기
                print ( "\t\t\tE step 이후 : {}\tmakeone_index : {}".format( np.unique(step_soft.membership  , return_counts=True), step_soft.makeone_index ) )
                step_soft = Mstep.main(df, np_vaf, np_BQ, step_soft, "Soft", **kwargs)     # 새 memberhsip에서 새 mixture구하기
                print("\t\t\tM step 이후 : fp_index {}\tmakeone_index {}\tlikelihood : {}".format( step_soft.fp_index, step_soft.makeone_index, round (step_soft.likelihood, 1), step_soft.mixture ))


                step_soft.acc(step_soft.mixture, step_soft.membership, step_soft.likelihood, step_soft.membership_p, step_soft.membership_p_normalize, step_soft.makeone_index, step_soft.fp_index, step_index + 1, step_soft.fp_member_index, step_soft.includefp, step_soft.fp_involuntary, step_soft.makeone_prenormalization, step_index, step_index)

                if (miscellaneous.GoStop(step_soft, **kwargs) == "Stop")  :
                    break
                if ( EMsoft.iszerocolumn (step_soft, **kwargs) == True) :
                    print ("\t\t\t\t→ 빈 mixture가 있어서 종료\t{}".format(step_soft.mixture))
                    break
                if ( len ( set (step_soft.membership) ) < NUM_CLONE ) :
                    print ("\t\t\t\t→ 빈 clone이 있어서 종료")
                    break


            step_soft.max_step_index =  step_soft.find_max_likelihood(1, step_soft.stepindex - 2 )   # 합쳐서 무조건 1이 되게 한다면 현실과 안 맞을수도 있음...
            i = step_soft.max_step_index

            # soft clustering에서 아예 답을 못 찾을 경우
            if i == -1:
                print ("\t\t\t1번째 step부터 망해서 이번 clone은 망함")
            elif  (step_soft.likelihood_record [i]  <= -9999999) :
                print ("\t\t\t모든 step에서 망해서 (-9999999) 이번 clone은 망함")

            else:  # 대부분의경우:  Soft clustering에서 답을 찾은 경우
                os.system ("cp " + kwargs["CLEMENT_DIR"] + "/trial/clone" + str (kwargs["NUM_CLONE"]) + "." + str( kwargs["TRIAL"] ) + "-"  + str(step_soft.max_step_index  + cluster_hard.stepindex_record [ NUM_CLONE ] - 1) + "\(soft\)." + kwargs["IMAGE_FORMAT"] + "  " + kwargs["CLEMENT_DIR"] + "/candidate/clone" + str (kwargs["NUM_CLONE"])  + ".\(soft\)." + kwargs["IMAGE_FORMAT"]  )
                cluster_soft.acc ( step_soft.mixture_record [i], step_soft.membership_record [i], step_soft.likelihood_record [i], step_soft.membership_p_record [i], step_soft.membership_p_normalize_record [i], step_soft.stepindex_record[i], cluster_hard.trialindex, step_soft.max_step_index_record[i], step_soft.makeone_index_record[i], step_soft.fp_index_record[i], step_soft.includefp_record[i], step_soft.fp_involuntary_record[i], step_soft.fp_member_index_record[i]   ,**kwargs )


        else:   # hard clustering에서 아예 답을 못 찾은 경우
            print ("Hard clustering에서조차 모두 망해서 이번 clone은 더 돌리기 싫다")


subprocess.run (["cp -r " + kwargs["CLEMENT_DIR"] + "/candidate  " + kwargs["COMBINED_OUTPUT_DIR"]  ], shell = True)
subprocess.run (["cp -r " + kwargs["CLEMENT_DIR"] + "/trial  " + kwargs["COMBINED_OUTPUT_DIR"]  ], shell = True)






print ("\n\n\n\n==================================== STEP #6.  OPTIMAL K DETERMINATION  =======================================")

NUM_CLONE_hard , NUM_CLONE_soft = [], []    # Hard clustering에서의 order, Soft clustering에서의 order

print ("\n\n★★★ Gap Statistics method (Hard clustering)\n")

f = io.StringIO()
with contextlib.redirect_stdout(f):
    NUM_CLONE_hard = miscellaneous.decision_gapstatistics (cluster_hard, np_vaf, **kwargs)
print ( f.getvalue() )
with open (kwargs["CLEMENT_DIR"] + "/CLEMENT_hard.gapstatistics.txt", "w", encoding = "utf8") as gap_myEM:
    print (f.getvalue(), file = gap_myEM)
subprocess.run (["cp " + kwargs["CLEMENT_DIR"] + "/CLEMENT_hard.gapstatistics.txt  " + kwargs["COMBINED_OUTPUT_DIR"]  + "/CLEMENT_hard.gapstatistics.txt"], shell = True)

if kwargs["MODE"] in ["Soft", "Both"]:
    if NUM_BLOCK >= 1:
        print ("\n\n\n★★★ XieBeni index method (2D, 3D Soft clustering)\n")

        f = io.StringIO()
        with contextlib.redirect_stdout(f):
            NUM_CLONE_soft = miscellaneous.decision_XieBeni (cluster_soft, np_vaf, **kwargs)

        print ( f.getvalue() )
        with open (kwargs["CLEMENT_DIR"] + "/CLEMENT_soft.xiebeni.txt", "w", encoding = "utf8") as xiebeni_myEM:
            print (f.getvalue(), file = xiebeni_myEM)
        subprocess.run (["cp " + kwargs["CLEMENT_DIR"] + "/CLEMENT_soft.xiebeni.txt  " + kwargs["COMBINED_OUTPUT_DIR"]  + "/CLEMENT_soft.xiebeni.txt"], shell = True)

    if NUM_BLOCK == 1:
        print ("\n\n\n★★★ Max likelihood method (1D Soft clustering)\n")

        f = io.StringIO()
        with contextlib.redirect_stdout(f):
            NUM_CLONE_soft = miscellaneous.decision_max (cluster_soft, np_vaf, **kwargs)

        print ( f.getvalue() )
        with open (kwargs["CLEMENT_DIR"] + "/CLEMENT_soft.maxlikelihood.txt", "w", encoding = "utf8") as maxlikelihood_myEM:
            print (f.getvalue(), file = maxlikelihood_myEM)
        subprocess.run (["cp " + kwargs["CLEMENT_DIR"] + "/CLEMENT_soft.maxlikelihood.txt  " + kwargs["COMBINED_OUTPUT_DIR"]  + "/CLEMENT_soft.maxlikelihood.txt"], shell = True)


print ("\n현재 시각 : {}h:{}m:{}s    (걸린 시간 : {})\n".format(time.localtime().tm_hour, time.localtime().tm_min, round(time.localtime().tm_sec), datetime.datetime.now() - START_TIME ))





print ("\n\n\n\n=============================================== STEP #7.  SCORING :  EM HARD  =======================================")

subprocess.run (["cp -r " +  kwargs["CLEMENT_DIR"]+ "/candidate  " + kwargs["COMBINED_OUTPUT_DIR"] ], shell = True)
DECISION = "hard_1st"
max_score_CLEMENT = 0

if kwargs["MODE"] in ["Hard", "Both"]:
    if kwargs["SCORING"] == True:
        print ("\n\nOrder : {}".format(NUM_CLONE_hard))

        for i, priority in enumerate(["1st"]):   # "2nd"
            if ( i >= len (NUM_CLONE_hard) ):
                break
            if ( cluster_hard.mixture_record [NUM_CLONE_hard[i]] == []):
                break
            
            if cluster_soft.mixture_record [ NUM_CLONE_hard[i] ] == []:       # Soft clustering이 완전히 망했을 때
                if (priority == "1st") & (kwargs["MODE"] in ["Both"]):
                    DECISION = "hard_1st"
                    print ("DECISION\t{}".format(DECISION))
                    with open (kwargs["CLEMENT_DIR"] + "/result/CLEMENT_hard_vs_fuzzy.txt", "w", encoding = "utf8") as output_hard_vs_fuzzy:
                        print ("DECISION\t{}".format(DECISION) , file = output_hard_vs_fuzzy )   
            else:
                if (priority == "1st") & (kwargs["MODE"] in ["Both"]):
                    moved_col_list = miscellaneous.movedcolumn ( cluster_hard, cluster_soft,  NUM_CLONE_hard[i]  )
                    hard_std = np.std(  cluster_hard.mixture_record [ NUM_CLONE_hard[i] ] [ : , moved_col_list ]   )
                    soft_std = np.std(  cluster_soft.mixture_record [ NUM_CLONE_hard[i] ] [ : , moved_col_list ]   )
                    if ( len (moved_col_list) == 1 )  & (  kwargs["NUM_BLOCK"] == 1 ) :  # Incalculable std in this condition
                        not_moved_col_list = [i for i in list(range ( NUM_CLONE_hard[i] )) if i not in moved_col_list]
                        print ( "moved_col_list = {}\nnot_moved_col_list = {}\ncluster_hard.mixture_record = {}".format (moved_col_list, not_moved_col_list, cluster_hard.mixture_record [ NUM_CLONE_hard[i] ] ))
                        not_moved_col_mean  =  np.mean ( cluster_hard.mixture_record [ NUM_CLONE_hard[i] ]  [0] [ not_moved_col_list ] ) 
                        moved_col_mean = np.mean ( cluster_hard.mixture_record [ NUM_CLONE_hard[i] ] [ 0 , moved_col_list ]  ) 
                        print ( "not_moved_col_mean = {}\nmoved_col_mean = {}".format (not_moved_col_mean, moved_col_mean) )

                        if abs (moved_col_mean - not_moved_col_mean) > 0.1:  
                            DECISION = "hard_1st"
                        else:
                            DECISION = "soft_1st"
                    
                    if  (soft_std < hard_std * kwargs["DECISION_STANDARD"]):
                        DECISION = "soft_1st"
                    
                    with open (kwargs["COMBINED_OUTPUT_DIR"] + "/result/CLEMENT_hard_vs_fuzzy.txt", "w", encoding = "utf8") as output_hard_vs_fuzzy:
                        print ( "Moved column : {}".format(moved_col_list), file = output_hard_vs_fuzzy)
                        print ( "Hard (n = {}) : std = {}\tstep = {}\nhard_mixture = {}\n".format( cluster_hard.mixture_record [NUM_CLONE_hard[i]].shape[1],  round( hard_std, 3) ,  cluster_hard.stepindex_record [ NUM_CLONE_hard[i] ],  cluster_hard.mixture_record [ NUM_CLONE_hard[i] ]   ) , file = output_hard_vs_fuzzy  )
                        print ( "Soft (n = {}) : std = {}\tstep = {}\nhard_mixture = {}".format( cluster_soft.mixture_record [NUM_CLONE_hard[i]].shape[1],  round( soft_std, 3) ,  cluster_soft.stepindex_record [ NUM_CLONE_hard[i] ],  cluster_soft.mixture_record [ NUM_CLONE_hard[i] ]   )  , file = output_hard_vs_fuzzy )
                        print ( "ratio : {}".format ( round(soft_std, 3) / round(hard_std, 3) ), file = output_hard_vs_fuzzy )
                        print ("\nsoft 선택 기준 :  < {}\nDECISION\t{}".format(kwargs["DECISION_STANDARD"], DECISION) , file = output_hard_vs_fuzzy )
                

            #정답set과 점수 맞춰보고 색깔 맞추기  (Hard clustering)
            
            score_df, score = scoring.mixturebased(mixture_answer, cluster_hard.mixture_record [NUM_CLONE_hard[i]], membership_answer, cluster_hard.membership_record [NUM_CLONE_hard[i]], kwargs["samplename_dict_CharacterToNum"], kwargs["samplename_dict_NumToCharacter"], cluster_hard.includefp_record [NUM_CLONE_hard[i]], cluster_hard.fp_index_record [NUM_CLONE_hard[i]] , "CLEMENT", **kwargs)
            
            max_score, sample_dict_PtoA, sample_dict_AtoP = scoring.Scoring ( membership_answer, membership_answer_numerical ,  cluster_hard.membership_record [NUM_CLONE_hard[i]], cluster_hard.fp_index_record [NUM_CLONE_hard[i]], 
                                                                                                                            set( list (range(0, NUM_CLONE )) ) - set( cluster_hard.makeone_index_record [NUM_CLONE_hard[i]] ) - set ( [ cluster_hard.fp_index_record [NUM_CLONE_hard[i]]  ] ), **kwargs  )
            if (max_score_CLEMENT < max_score) & ("hard" in DECISION) & (priority == "1st"):
                max_score_CLEMENT = max_score
            
            print ("max_score_CLEMENT = {}".format (max_score_CLEMENT))
           
            ARI_CLEMENT = result.ARI ( np.array ( [ membership_answer_numerical [j] for j in membership_answer_numerical_nofp_index  ] ) , 
                                                   np.array ( [ cluster_hard.membership_record [NUM_CLONE_hard[i]][j] for j in membership_answer_numerical_nofp_index] ) )
            print ("\n[■ {} BEST RESULTS]\n\nCLEMENT\t{}/{}\nNUM_CLONE\t{}\nARI\t{}\nmixture\t{}\nmakeone_index_record\t{}\n".format(priority, max_score,  kwargs["NUM_MUTATION"], cluster_hard.mixture_record [NUM_CLONE_hard[i]].shape[1], round (ARI_CLEMENT, 2) , list(cluster_hard.mixture_record [NUM_CLONE_hard[i]],), cluster_hard.makeone_index_record [NUM_CLONE_hard[i]]))
            #print ("\n(Greedy 방식) score : {}점 / {}점".format(score, kwargs["NUM_MUTATION"]))
            print ("(모두 다 돌려본 결과) score : {}점 / {}점,  변환표 (A to P)= {}\n".format ( max_score, kwargs["NUM_MUTATION"], sample_dict_AtoP))
            print ("{}".format(score_df))
            
            if (kwargs["FP_RATIO"] != 0) | (kwargs["FP_USEALL"] == "True"):
                if cluster_hard.fp_index_record [ NUM_CLONE_hard[i] ] == -1:        # FP가 분명히 없는 경우
                    answeronly, intersection, myonly, sensitivity, PPV, F1 = 0, 0, 0, None, None, None
                else:
                    #print ("[{} BEST FP ANALYSIS]".format(priority))
                    answeronly, intersection, myonly, sensitivity, PPV, F1 = result.FPmatrix(score_df)
                    #print ("answer FP {}개 중에 {}개 일치함".format( answeronly + intersection ,  intersection ) )
                    #print ("\tanswerFP only : {}\n\tintersection : {}\n\tmyFP only : {}\n\n".format( answeronly, intersection, myonly )
            else:
                answeronly, intersection, myonly, sensitivity, PPV, F1 = 0, 0, 0, None, None, None
                
            NUM_CLONE_CLEMENT = cluster_hard.mixture_record [NUM_CLONE_hard[i]].shape[1] - int (cluster_hard.includefp_record [ NUM_CLONE_hard[i] ])
            NUM_CHILD_CLEMENT = len (cluster_hard.makeone_index_record [NUM_CLONE_hard[i]])

            with open (kwargs["CLEMENT_DIR"] + "/result/CLEMENT_hard_" + priority + ".results.txt", "w", encoding = "utf8") as output_myEM:
                print ("NUM_CLONE\t{}\nNUM_CHILD\t{}\nscore\t{}/{}\nARI\t{}\nrunningtime\t{}\nFPexistence\t{}\nmakeone_index\t{}".
                    format(NUM_CLONE_CLEMENT, NUM_CHILD_CLEMENT, max_score, kwargs["NUM_MUTATION"], ARI_CLEMENT,  round((datetime.datetime.now() - START_TIME).total_seconds()),  cluster_hard.includefp_record [NUM_CLONE_hard[i]], cluster_hard.makeone_index_record [NUM_CLONE_hard[i]] ), file = output_myEM)
                
            subprocess.run (["cp " + kwargs["CLEMENT_DIR"] + "/result/CLEMENT_hard_" + priority + ".results.txt  " + kwargs["COMBINED_OUTPUT_DIR"]  + "/CLEMENT_hard_" + priority + ".results.txt"], shell = True)


            pd.DataFrame(cluster_hard.membership_record [NUM_CLONE_hard[i]]).to_csv (kwargs["CLEMENT_DIR"] + "/result/CLEMENT_hard_" + priority + ".membership.txt", index = False, header= False,  sep = "\t" )
            pd.DataFrame(cluster_hard.membership_record [NUM_CLONE_hard[i]]).to_csv (kwargs["COMBINED_OUTPUT_DIR"] + "/CLEMENT_hard_" + priority + ".membership.txt", index = False, header= False,  sep = "\t" )
            pd.DataFrame(cluster_hard.mixture_record [NUM_CLONE_hard[i]],).to_csv (kwargs["CLEMENT_DIR"] + "/result/CLEMENT_hard_" + priority + ".mixture.txt", index = False, header= False,  sep = "\t" )
            pd.DataFrame(cluster_hard.mixture_record [NUM_CLONE_hard[i]],).to_csv (kwargs["COMBINED_OUTPUT_DIR"] + "/CLEMENT_hard_" + priority + ".mixture.txt", index = False, header= False,  sep = "\t" )
            pd.DataFrame( np.unique( cluster_hard.membership_record [NUM_CLONE_hard[i]], return_counts = True ) ).to_csv (kwargs["CLEMENT_DIR"] + "/result/CLEMENT_hard_" + priority + ".membership_count.txt", index = False, header= False,  sep = "\t" )
            pd.DataFrame( np.unique( cluster_hard.membership_record [NUM_CLONE_hard[i]], return_counts = True ) ).to_csv (kwargs["COMBINED_OUTPUT_DIR"] + "/CLEMENT_hard_" + priority + ".membership_count.txt", index = False, header= False,  sep = "\t" )
            pd.DataFrame(score_df).to_csv (kwargs["CLEMENT_DIR"] + "/result/CLEMENT_hard_" + priority + ".scoredf.txt", index = False, header= True,  sep = "\t")
            pd.DataFrame(score_df).to_csv (kwargs["COMBINED_OUTPUT_DIR"] + "/CLEMENT_hard_" + priority + ".scoredf.txt", index = False, header= True,  sep = "\t")

            #visualizationpair
            if kwargs["NUM_BLOCK"] >= 2:
                visualizationpair.drawfigure_2d (membership_answer, mixture_answer, cluster_hard.membership_record [NUM_CLONE_hard[i]], np.round(cluster_hard.mixture_record [NUM_CLONE_hard[i]], 2),
                            score_df, kwargs["CLEMENT_DIR"]+ "/CLEMENT_hard_" + priority + "." + kwargs["IMAGE_FORMAT"] ,  "ANSWER\n", "CLEMENT\n{}/{}, ARI={}".format(score, kwargs["NUM_MUTATION"], round (ARI_CLEMENT, 2) ), np_vaf, cluster_hard.includefp_record [NUM_CLONE_hard[i]],  cluster_hard.makeone_index_record [NUM_CLONE_hard[i]], dimensionreduction="None", **kwargs)
            else:  # 어차피 pair로 안 그릴거면 그냥 복사해오자
                subprocess.run (["cp " +  kwargs["CLEMENT_DIR"]+ "/candidate/clone" + str( NUM_CLONE_hard[i] ) + ".\(hard\)." + kwargs["IMAGE_FORMAT"] + " "  + kwargs["CLEMENT_DIR"]+ "/CLEMENT_hard_" + priority + "." + kwargs["IMAGE_FORMAT"]], shell = True)    
                # samplename_dict = {k:"clone {}".format(k) for k in range(0, np.max( cluster_hard.membership_record [NUM_CLONE_hard[i]] )+ 1)}
                # visualizationsingle.drawfigure_1d(cluster_hard.membership_record [NUM_CLONE_hard[i]], "CLEMENT", kwargs["CLEMENT_DIR"]+ "/CLEMENT_hard_" + priority + "." + kwargs["IMAGE_FORMAT"], np_vaf,
                #                                                         samplename_dict, cluster_hard.includefp_record [NUM_CLONE_hard[i]], cluster_hard.fp_index_record [NUM_CLONE_hard[i]],  cluster_hard.makeone_index_record [NUM_CLONE_hard[i]])
            subprocess.run (["cp " +  kwargs["CLEMENT_DIR"]+ "/CLEMENT_hard_" + priority + "." + kwargs["IMAGE_FORMAT"] + " "  + kwargs["COMBINED_OUTPUT_DIR"]  + "/CLEMENT_hard_" + priority + "." + kwargs["IMAGE_FORMAT"]], shell = True)

            if len (cluster_hard.makeone_index_record [NUM_CLONE_hard[i]]) + int ( cluster_hard.includefp_record [NUM_CLONE_hard[i]])  < NUM_CLONE_hard[i]:    # FP
                print ("\n\n\t▶▶▶  PHYLOGENY RECONSTRUCTION IN HARD CLUSTERING ◀◀◀")
                kwargs["PHYLOGENY_DIR"] = kwargs["CLEMENT_DIR"] + "/CLEMENT_hard_" + priority + ".phylogeny.txt"
                f = io.StringIO()
                with contextlib.redirect_stdout(f):
                    membership_child = set ( cluster_hard.makeone_index_record[ NUM_CLONE_hard[i] ] )            # child의 번호만 뽑아준다 ( e.g.  0, 1, 3)
                    membership_outside = set (range (0, NUM_CLONE_hard [i] )) - membership_child - set ( [cluster_hard.fp_index_record[ NUM_CLONE_hard[i] ] ] )   # outlier의 번호만 뽑아준다 (e.g. 2)

                    g = phylogeny.main(membership_child, membership_outside, cluster_hard.mixture_record[ NUM_CLONE_hard[i] ],  **kwargs)

                print ( f.getvalue() )
                with open ( kwargs["PHYLOGENY_DIR"] , "w", encoding = "utf8") as phylogeny_file:
                    print (f.getvalue(), file = phylogeny_file)
                subprocess.run (["cp " + kwargs["PHYLOGENY_DIR"] + "  " + kwargs["COMBINED_OUTPUT_DIR"]  + "/CLEMENT_hard_" + priority + ".phylogeny.txt"], shell = True)
                


    else:  # SCORING 없는 단독 결과일 때
        for i, priority in enumerate(["1st"]):
            if i >= len (NUM_CLONE_hard):
                break
            
            if cluster_soft.mixture_record [ NUM_CLONE_hard[i] ] == []:
                if (priority == "1st") & (kwargs["MODE"] in ["Both"]):
                    print ("DECISION\t{}".format(DECISION))
                    with open (kwargs["COMBINED_OUTPUT_DIR"] + "/CLEMENT_hard_vs_fuzzy.txt", "w", encoding = "utf8") as output_hard_vs_fuzzy:
                        print ("DECISION\t{}".format(DECISION) , file = output_hard_vs_fuzzy )            
                
            else:
                if (priority == "1st") & (kwargs["MODE"] in ["Both"]):
                    moved_col_list = miscellaneous.movedcolumn ( cluster_hard, cluster_soft,  NUM_CLONE_hard[i]  )
                    hard_std = np.std(  cluster_hard.mixture_record [ NUM_CLONE_hard[i] ] [ : , moved_col_list ]   )
                    soft_std = np.std(  cluster_soft.mixture_record [ NUM_CLONE_hard[i] ] [ : , moved_col_list ]   )
                    if ( len (moved_col_list) == 1 )  & (  kwargs["NUM_BLOCK"] == 1 ) :  # Incalculable std in this condition
                        not_moved_col_list = [i for i in list(range ( NUM_CLONE_hard[i] )) if i not in moved_col_list]
                        print ( "moved_col_list = {}\nnot_moved_col_list = {}\ncluster_hard.mixture_record = {}".format (moved_col_list, not_moved_col_list, cluster_hard.mixture_record [ NUM_CLONE_hard[i] ] ))
                        not_moved_col_mean  =  np.mean ( cluster_hard.mixture_record [ NUM_CLONE_hard[i] ]  [0] [ not_moved_col_list ] ) 
                        moved_col_mean = np.mean ( cluster_hard.mixture_record [ NUM_CLONE_hard[i] ] [ 0 , moved_col_list ]  ) 
                        print ( "not_moved_col_mean = {}\nmoved_col_mean = {}".format (not_moved_col_mean, moved_col_mean) )

                        if abs (moved_col_mean - not_moved_col_mean) > 0.1:  
                            DECISION = "hard_1st"
                        else:
                            DECISION = "soft_1st"
                    
                    if  (soft_std < hard_std * kwargs["DECISION_STANDARD"]):
                        DECISION = "soft_1st"
                    
                    with open (kwargs["COMBINED_OUTPUT_DIR"] + "/result/CLEMENT_hard_vs_fuzzy.txt", "w", encoding = "utf8") as output_hard_vs_fuzzy:
                        print ( "Moved column : {}".format(moved_col_list), file = output_hard_vs_fuzzy)
                        print ( "Hard (n = {}) : std = {}\tstep = {}\nhard_mixture = {}\n".format( cluster_hard.mixture_record [NUM_CLONE_hard[i]].shape[1],  round( hard_std, 3) ,  cluster_hard.stepindex_record [ NUM_CLONE_hard[i] ],  cluster_hard.mixture_record [ NUM_CLONE_hard[i] ]   ) , file = output_hard_vs_fuzzy  )
                        print ( "Soft (n = {}) : std = {}\tstep = {}\nhard_mixture = {}".format( cluster_soft.mixture_record [NUM_CLONE_hard[i]].shape[1],  round( soft_std, 3) ,  cluster_soft.stepindex_record [ NUM_CLONE_hard[i] ],  cluster_soft.mixture_record [ NUM_CLONE_hard[i] ]   )  , file = output_hard_vs_fuzzy )
                        print ( "ratio : {}".format ( round(soft_std, 3) / round(hard_std, 3) ), file = output_hard_vs_fuzzy )
                        print ("\nsoft 선택 기준 :  < {}\nDECISION\t{}".format(kwargs["DECISION_STANDARD"], DECISION) , file = output_hard_vs_fuzzy )


            with open (kwargs["CLEMENT_DIR"] + "/result/CLEMENT_hard_" + priority + ".results.txt", "w", encoding = "utf8") as output_myEM:
                print ("NUM_CLONE\t{}\nNUM_CHILD\t{}\nrunningtime\t{}\nFPexistence\t{}\nFPindex\t{}\nmakeone_index\t{}".
                        format(cluster_hard.mixture_record [NUM_CLONE_hard[i]].shape[1] - int (cluster_hard.includefp_record [ NUM_CLONE_hard[i] ]) , len (cluster_hard.makeone_index_record [NUM_CLONE_hard[i]]),   round((datetime.datetime.now() - START_TIME).total_seconds()), cluster_hard.includefp_record [NUM_CLONE_hard[i]], cluster_hard.fp_index_record [NUM_CLONE_hard[i]] , cluster_hard.makeone_index_record [NUM_CLONE_hard[i]]  ), file = output_myEM)
            subprocess.run (["cp " + kwargs["CLEMENT_DIR"] + "/result/CLEMENT_hard_" + priority + ".results.txt  " + kwargs["COMBINED_OUTPUT_DIR"]  + "/CLEMENT_hard_" + priority + ".results.txt"], shell = True)

            pd.DataFrame(cluster_hard.membership_record [NUM_CLONE_hard[i]]).to_csv (kwargs["CLEMENT_DIR"] + "/result/CLEMENT_hard_" + priority + ".membership.txt", index = False, header= False,  sep = "\t" )
            pd.DataFrame(cluster_hard.membership_record [NUM_CLONE_hard[i]]).to_csv (kwargs["COMBINED_OUTPUT_DIR"] + "/CLEMENT_hard_" + priority + ".membership.txt", index = False, header= False,  sep = "\t" )
            pd.DataFrame(cluster_hard.mixture_record [NUM_CLONE_hard[i]],).to_csv (kwargs["CLEMENT_DIR"] + "/result/CLEMENT_hard_" + priority + ".mixture.txt", index = False, header= False,  sep = "\t" )
            pd.DataFrame(cluster_hard.mixture_record [NUM_CLONE_hard[i]],).to_csv (kwargs["COMBINED_OUTPUT_DIR"] + "/CLEMENT_hard_" + priority + ".mixture.txt", index = False, header= False,  sep = "\t" )
            pd.DataFrame( np.unique( cluster_hard.membership_record [NUM_CLONE_hard[i]], return_counts = True ) ).to_csv (kwargs["CLEMENT_DIR"] + "/result/CLEMENT_hard_" + priority + ".membership_count.txt", index = False, header= False,  sep = "\t" )
            pd.DataFrame( np.unique( cluster_hard.membership_record [NUM_CLONE_hard[i]], return_counts = True ) ).to_csv (kwargs["COMBINED_OUTPUT_DIR"] + "/CLEMENT_hard_" + priority + ".membership_count.txt", index = False, header= False,  sep = "\t" )

            samplename_dict = {k:k for k in range(0, np.max(cluster_hard.membership_record [NUM_CLONE_hard[i]])+ 1)}

            if kwargs["NUM_BLOCK"] == 1:
                visualizationsingle.drawfigure_1d (cluster_hard.membership_record [NUM_CLONE_hard[i]], "CLEMENT_hard : {}".format( round ( cluster_hard.likelihood_record [NUM_CLONE_hard[i]] ) ), kwargs["CLEMENT_DIR"]+ "/CLEMENT_hard_" + priority + "." + kwargs["IMAGE_FORMAT"], np_vaf, samplename_dict, cluster_hard.includefp_record [NUM_CLONE_hard[i]] , cluster_hard.fp_index_record[NUM_CLONE_hard[i]], cluster_hard.makeone_index_record [NUM_CLONE_hard[i]], **kwargs )
            elif kwargs["NUM_BLOCK"] == 2:
                visualizationsingle.drawfigure_2d (cluster_hard.membership_record [NUM_CLONE_hard[i]], "CLEMENT_hard : {}".format( round ( cluster_hard.likelihood_record [NUM_CLONE_hard[i]] ) ), kwargs["CLEMENT_DIR"]+ "/CLEMENT_hard_" + priority + "." + kwargs["IMAGE_FORMAT"], np_vaf, samplename_dict, cluster_hard.includefp_record [NUM_CLONE_hard[i]] , cluster_hard.fp_index_record[NUM_CLONE_hard[i]], "None", **kwargs )
            else:
                visualizationsingle.drawfigure_2d (cluster_hard.membership_record [NUM_CLONE_hard[i]], "CLEMENT_hard : {}".format( round ( cluster_hard.likelihood_record [NUM_CLONE_hard[i]] ) ), kwargs["CLEMENT_DIR"]+ "/CLEMENT_hard_" + priority + "." + kwargs["IMAGE_FORMAT"], np_vaf, samplename_dict, cluster_hard.includefp_record [NUM_CLONE_hard[i]] , cluster_hard.fp_index_record[NUM_CLONE_hard[i]], "SVD", **kwargs )


            subprocess.run (["cp -rf " + kwargs["CLEMENT_DIR"] + "/candidate/clone" + str(NUM_CLONE_hard[i]) +".\(hard\)." + kwargs["IMAGE_FORMAT"] + " " + kwargs["CLEMENT_DIR"]+ "/CLEMENT_hard_" + priority + "." + kwargs["IMAGE_FORMAT"]], shell = True)
            subprocess.run (["cp -rf " +  kwargs["CLEMENT_DIR"]+ "/CLEMENT_hard_" + priority + "." + kwargs["IMAGE_FORMAT"] + "  " + kwargs["COMBINED_OUTPUT_DIR"]  + "/CLEMENT_hard_" + priority + "." + kwargs["IMAGE_FORMAT"]], shell = True)
            print ("\tHard {} results printed".format(priority))
        









if kwargs["MODE"] in ["Soft", "Both"]:
    print ("\n\n\n======================================= STEP #8.  SCORING :  EM SOFT  =======================================")
    if kwargs["SCORING"] == True:
        print ("\n\nOrder : {}".format(NUM_CLONE_soft))
        for i, priority in enumerate(["1st"]):        
            if i >= len (NUM_CLONE_soft):
                break
            if ( cluster_hard.mixture_record [NUM_CLONE_soft[i]] == []): 
                break
            if len (cluster_soft.makeone_index_record [NUM_CLONE_soft[i]]) == 0:
                print ( "NUM_CHILD = 0" )
                print ( cluster_soft.mixture_record [NUM_CLONE_soft[i]] )
                break


            #정답set과 점수 맞춰보고 색깔 맞추기  (Soft clustering)
            score_df, score = scoring.mixturebased(mixture_answer, cluster_soft.mixture_record [NUM_CLONE_soft[i]], membership_answer, cluster_soft.membership_record [NUM_CLONE_soft[i]], kwargs["samplename_dict_CharacterToNum"], kwargs["samplename_dict_NumToCharacter"], cluster_soft.includefp_record [NUM_CLONE_soft[i]], cluster_soft.fp_index_record [NUM_CLONE_soft[i]], "CLEMENT", **kwargs)
        
            max_score, sample_dict_PtoA, sample_dict_AtoP = scoring.Scoring ( membership_answer, membership_answer_numerical ,  cluster_soft.membership_record [NUM_CLONE_soft[i]], cluster_soft.fp_index_record [NUM_CLONE_soft[i]], 
                                                                                                                            set( list (range(0, NUM_CLONE )) ) - set( cluster_soft.makeone_index_record [NUM_CLONE_soft[i]] ) - set ( [ cluster_soft.fp_index_record [NUM_CLONE_soft[i]]  ] ), **kwargs  )
            if (max_score_CLEMENT < max_score) & ("soft" in DECISION) & (priority == "1st"):
                max_score_CLEMENT = max_score

            ARI_CLEMENT = result.ARI ( np.array ( [ membership_answer_numerical [j] for j in membership_answer_numerical_nofp_index  ] ) , 
                                                   np.array ( [ cluster_soft.membership_record [NUM_CLONE_soft[i]] [j] for j in membership_answer_numerical_nofp_index] ) )

            #visualization
            if kwargs["NUM_BLOCK"] >= 3:
                visualizationpair.drawfigure_2d (membership_answer, mixture_answer, cluster_soft.membership_record [NUM_CLONE_soft[i]], np.round(cluster_soft.mixture_record [NUM_CLONE_soft[i]], 2),
                            score_df, kwargs["CLEMENT_DIR"]+ "/CLEMENT_soft_" + priority + "." + kwargs["IMAGE_FORMAT"] ,  "ANSWER\n", "CLEMENT ", np_vaf, cluster_soft.includefp_record [NUM_CLONE_soft[i]],  cluster_soft.makeone_index_record [NUM_CLONE_soft[i]], dimensionreduction="SVD", **kwargs)
            elif kwargs["NUM_BLOCK"] == 2:
                visualizationpair.drawfigure_2d (membership_answer, mixture_answer, cluster_soft.membership_record [NUM_CLONE_soft[i]], np.round(cluster_soft.mixture_record [NUM_CLONE_soft[i]], 2),
                            score_df, kwargs["CLEMENT_DIR"]+ "/CLEMENT_soft_" + priority + "." + kwargs["IMAGE_FORMAT"] ,  "ANSWER\n", "CLEMENT ", np_vaf, cluster_soft.includefp_record [NUM_CLONE_soft[i]],  cluster_soft.makeone_index_record [NUM_CLONE_soft[i]], dimensionreduction="None", **kwargs)
                samplename_dict = {k:k for k in range(0, np.max(cluster_soft.membership_record [NUM_CLONE_soft[i]])+ 1)}
                # visualizationsinglesoft.drawfigure_2d( cluster_soft.membership_record [NUM_CLONE_soft[i]], cluster_soft.mixture_record [NUM_CLONE_soft[i]], cluster_soft.membership_p_normalize_record [NUM_CLONE_soft[i]], 
                #                                                             "CLEMENT", kwargs["CLEMENT_DIR"]+ "/CLEMENT_soft_soft_" + priority + "." + kwargs["IMAGE_FORMAT"], np_vaf, samplename_dict, cluster_soft.includefp_record [NUM_CLONE_soft[i]], cluster_soft.fp_index_record[NUM_CLONE_soft[i]], cluster_soft.makeone_index_record[NUM_CLONE_soft[i]], dimensionreduction="None")
            elif kwargs["NUM_BLOCK"] == 1:
                samplename_dict = {k:k for k in range(0, np.max(cluster_soft.membership_record [NUM_CLONE_soft[i]])+ 1)}
                visualizationsinglesoft.drawfigure_1d (cluster_soft.membership_record [NUM_CLONE_soft[i]], cluster_soft.mixture_record [NUM_CLONE_soft[i]], cluster_soft.membership_p_normalize_record [NUM_CLONE_soft[i]],
                                                                        "CLEMENT_soft : {}".format( round (cluster_soft.likelihood_record [NUM_CLONE_soft[i]] ) ), kwargs["CLEMENT_DIR"]+ "/CLEMENT_soft_" + priority + "." + kwargs["IMAGE_FORMAT"], np_vaf, samplename_dict, cluster_soft.includefp_record [NUM_CLONE_soft[i]], cluster_soft.fp_index_record[NUM_CLONE_soft[i]], cluster_soft.makeone_index_record[NUM_CLONE_soft[i]] )
            
            subprocess.run (["cp " +  kwargs["CLEMENT_DIR"]+ "/CLEMENT_soft_" + priority + "." + kwargs["IMAGE_FORMAT"]  + " " + kwargs["COMBINED_OUTPUT_DIR"]  + "/CLEMENT_soft_" + priority + "." + kwargs["IMAGE_FORMAT"]], shell = True)


            print ("\n[■ {} BEST RESULTS]\n\nCLEMENT\t{}/{}\nNUM_CLONE\t{}\n\nARI\t{}\nmakeone_index\t{}\n{}\n".format(priority, max_score,  kwargs["NUM_MUTATION"], cluster_soft.mixture_record [NUM_CLONE_soft[i]].shape[1], ARI_CLEMENT, cluster_soft.makeone_index_record [NUM_CLONE_soft[i]], list(cluster_soft.mixture_record [NUM_CLONE_soft[i]], )))
            
            #print ("\n(Greedy 방식) score : {}점 / {}점".format(score, kwargs["NUM_MUTATION"]))
            print ("(모두 다 돌려본 결과) score : {}점 / {}점,  변환표 (A to P)= {}\n".format ( max_score, kwargs["NUM_MUTATION"], sample_dict_AtoP))
            print ("{}".format(score_df))
            
            if (kwargs["FP_RATIO"] != 0) | (kwargs["FP_USEALL"] == "True"):
                if cluster_soft.fp_index_record [ NUM_CLONE_soft[i] ] == -1:        # FP가 분명히 없는 경우
                    answeronly, intersection, myonly, sensitivity, PPV, F1 = 0, 0, 0, None, None, None
                else:
                    #print ("[FP ANALYSIS]")
                    answeronly, intersection, myonly, sensitivity, PPV, F1 = result.FPmatrix(score_df)
                    #print ("answer FP {}개 중에 {}개 일치함".format( answeronly + intersection ,  intersection ) )
                    #print ("\tanswerFP only : {}\n\tintersection : {}\n\tmyFP only : {}\n".format( answeronly, intersection, myonly ))
            else:
                answeronly, intersection, myonly, sensitivity, PPV, F1 = 0, 0, 0, None, None, None

            NUM_CLONE_CLEMENT = cluster_soft.mixture_record [NUM_CLONE_soft[i]].shape[1] - int (cluster_soft.includefp_record [ NUM_CLONE_soft[i] ]) 
            NUM_CHILD_CLEMENT =  len (cluster_soft.makeone_index_record [NUM_CLONE_soft[i]])

            with open (kwargs["CLEMENT_DIR"] + "/result/CLEMENT_soft_" + priority + ".results.txt", "w", encoding = "utf8") as output_myEM:
                print ("NUM_CLONE\t{}\nNUM_CHILD\t{}\nscore\t{}/{}\nARI\t{}\nrunningtime\t{}\nFPexistence\t{}\nmakeone_index\t{}\n".
                    format(NUM_CLONE_CLEMENT, NUM_CHILD_CLEMENT, max_score, kwargs["NUM_MUTATION"],  ARI_CLEMENT, round((datetime.datetime.now() - START_TIME).total_seconds()), cluster_soft.includefp_record [NUM_CLONE_soft[i]], cluster_soft.makeone_index_record [NUM_CLONE_soft[i]]  ), file = output_myEM)
            subprocess.run (["cp -rf " + kwargs["CLEMENT_DIR"] + "/result/CLEMENT_soft_" + priority + ".results.txt  " + kwargs["COMBINED_OUTPUT_DIR"]  + "/CLEMENT_soft_" + priority + ".results.txt"], shell = True)


            pd.DataFrame(cluster_soft.membership_record [NUM_CLONE_soft[i]]).to_csv (kwargs["CLEMENT_DIR"] + "/result/CLEMENT_soft_" + priority + ".membership.txt", index = False, header= False,  sep = "\t" )
            pd.DataFrame(cluster_soft.membership_record [NUM_CLONE_soft[i]]).to_csv (kwargs["COMBINED_OUTPUT_DIR"] + "/CLEMENT_soft_" + priority + ".membership.txt", index = False, header= False,  sep = "\t" )
            pd.DataFrame(cluster_soft.mixture_record [NUM_CLONE_soft[i]],).to_csv (kwargs["CLEMENT_DIR"] + "/result/CLEMENT_soft_" + priority + ".mixture.txt", index = False, header= False,  sep = "\t" )
            pd.DataFrame(cluster_soft.mixture_record [NUM_CLONE_soft[i]],).to_csv (kwargs["COMBINED_OUTPUT_DIR"] + "/CLEMENT_soft_" + priority + ".mixture.txt", index = False, header= False,  sep = "\t" )
            pd.DataFrame( np.unique( cluster_soft.membership_record [NUM_CLONE_soft[i]], return_counts = True ) ).to_csv (kwargs["CLEMENT_DIR"] + "/result/CLEMENT_soft_" + priority + ".membership_count.txt", index = False, header= False,  sep = "\t" )
            pd.DataFrame( np.unique( cluster_soft.membership_record [NUM_CLONE_soft[i]], return_counts = True ) ).to_csv (kwargs["COMBINED_OUTPUT_DIR"] + "/CLEMENT_soft_" + priority + ".membership_count.txt", index = False, header= False,  sep = "\t" )
            pd.DataFrame(score_df).to_csv (kwargs["CLEMENT_DIR"] + "/result/CLEMENT_soft_" + priority + ".scoredf.txt", index = False, header= True,  sep = "\t")
            pd.DataFrame(score_df).to_csv (kwargs["COMBINED_OUTPUT_DIR"] + "/CLEMENT_soft_" + priority + ".scoredf.txt", index = False, header= True,  sep = "\t")

            if len (cluster_soft.makeone_index_record [NUM_CLONE_soft[i]]) + int ( cluster_soft.includefp_record [NUM_CLONE_soft[i]])  < NUM_CLONE_soft[i]:    # FP
                ISPARENT = True
            else:
                ISPARENT = False

            if ISPARENT == True:
                print ("\n\n\n\n▶▶▶  PHYLOGENY ANALYSIS ◀◀◀")
                kwargs["PHYLOGENY_DIR"] = kwargs["CLEMENT_DIR"] + "/CLEMENT_soft_" + priority + ".phylogeny.txt"
                f = io.StringIO()
                with contextlib.redirect_stdout(f):
                    membership_child = set ( cluster_soft.makeone_index_record[ NUM_CLONE_soft[i] ] )            # child의 번호만 뽑아준다 ( e.g.  0, 1, 3)
                    membership_outside = set (range (0, NUM_CLONE_soft [i] )) - membership_child - set ( [cluster_soft.fp_index_record[ NUM_CLONE_soft[i] ] ] )   # outlier의 번호만 뽑아준다 (e.g. 2)

                    g = phylogeny.main(membership_child, membership_outside, cluster_soft.mixture_record[ NUM_CLONE_soft[i] ],  **kwargs)

                print ( f.getvalue() )
                with open ( kwargs["PHYLOGENY_DIR"] , "w", encoding = "utf8") as phylogeny_file:
                    print (f.getvalue(), file = phylogeny_file)
                subprocess.run (["cp -rf " + kwargs["PHYLOGENY_DIR"] + "  " + kwargs["COMBINED_OUTPUT_DIR"]  + "/CLEMENT_hard_" + priority + ".phylogeny.txt"], shell = True)

    else:     # Scoring 안 하는 경우
        for i, priority in enumerate(["1st"]):
            if i >= len (NUM_CLONE_soft):
                break

            if cluster_soft.mixture_record [NUM_CLONE_soft[i]] == []:
                print ("Soft : Empty")
                break

            print ( "NUM_CLONE_soft (by order): {}".format(NUM_CLONE_soft))

            with open (kwargs["CLEMENT_DIR"] + "/result/CLEMENT_soft_" + priority + ".results.txt", "w", encoding = "utf8") as output_myEM:
                try:
                    print ("NUM_CLONE\t{}\nNUM_CHILD\t{}\nrunningtime\t{}\nFPexistence\t{}\nFPindex\t{}".
                        format(cluster_soft.mixture_record [NUM_CLONE_soft[i]].shape[1] - int (cluster_soft.includefp_record [ NUM_CLONE_soft[i] ]) , len (cluster_soft.makeone_index_record [NUM_CLONE_soft[i]]),   round((datetime.datetime.now() - START_TIME).total_seconds()),  cluster_soft.includefp_record [NUM_CLONE_soft[i]], cluster_soft.fp_index_record [NUM_CLONE_soft[i]]  ), file = output_myEM)
                except:
                    print ("NUM_CLONE\t-1\nNUM_CHILD\t{}\nrunningtime\t{}\nFPexistence\t{}\nFPindex\t{}", file = output_myEM)
            subprocess.run (["cp -rf " + kwargs["CLEMENT_DIR"] + "/result/CLEMENT_soft_" + priority + ".results.txt  " + kwargs["COMBINED_OUTPUT_DIR"]  + "/CLEMENT_soft_" + priority + ".results.txt"], shell = True)


            pd.DataFrame(cluster_soft.membership_record [NUM_CLONE_soft[i]]).to_csv (kwargs["CLEMENT_DIR"] + "/CLEMENT_soft_" + priority + ".membership.txt", index = False, header= False,  sep = "\t" )
            pd.DataFrame(cluster_soft.membership_record [NUM_CLONE_soft[i]]).to_csv (kwargs["COMBINED_OUTPUT_DIR"] + "/CLEMENT_soft_" + priority + ".membership.txt", index = False, header= False,  sep = "\t" )
            pd.DataFrame(cluster_soft.mixture_record [NUM_CLONE_soft[i]],).to_csv (kwargs["CLEMENT_DIR"] + "/CLEMENT_soft_" + priority + ".mixture.txt", index = False, header= False,  sep = "\t" )
            pd.DataFrame(cluster_soft.mixture_record [NUM_CLONE_soft[i]],).to_csv (kwargs["COMBINED_OUTPUT_DIR"] + "/CLEMENT_soft_" + priority + ".mixture.txt", index = False, header= False,  sep = "\t" )
            pd.DataFrame( np.unique( cluster_soft.membership_record [NUM_CLONE_soft[i]], return_counts = True ) ).to_csv (kwargs["CLEMENT_DIR"] + "/CLEMENT_soft_" + priority + ".membership_count.txt", index = False, header= False,  sep = "\t" )
            pd.DataFrame( np.unique( cluster_soft.membership_record [NUM_CLONE_soft[i]], return_counts = True ) ).to_csv (kwargs["COMBINED_OUTPUT_DIR"] + "/CLEMENT_soft_" + priority + ".membership_count.txt", index = False, header= False,  sep = "\t" )

            samplename_dict = {k:k for k in range(0, np.max(cluster_soft.membership_record [NUM_CLONE_soft[i]])+ 1)}

            if kwargs["NUM_BLOCK"] == 1:
                visualizationsinglesoft.drawfigure_1d (cluster_soft.membership_record [NUM_CLONE_soft[i]], cluster_soft.mixture_record [NUM_CLONE_soft[i]], cluster_soft.membership_p_normalize_record [NUM_CLONE_soft[i]],
                                                                        "CLEMENT_soft : {}".format( round (cluster_soft.likelihood_record [NUM_CLONE_soft[i]] ) ), kwargs["CLEMENT_DIR"]+ "/CLEMENT_soft_" + priority + "." + kwargs["IMAGE_FORMAT"], np_vaf, samplename_dict, cluster_soft.includefp_record [NUM_CLONE_soft[i]] , cluster_soft.fp_index_record[NUM_CLONE_soft[i]] , cluster_soft.makeone_index_record[NUM_CLONE_soft[i]] )
            elif kwargs["NUM_BLOCK"] == 2:
                visualizationsinglesoft.drawfigure_2d (cluster_soft.membership_record [NUM_CLONE_soft[i]], cluster_soft.mixture_record [NUM_CLONE_soft[i]], cluster_soft.membership_p_normalize_record [NUM_CLONE_soft[i]],
                                                                        "CLEMENT_soft : {}".format( round (cluster_soft.likelihood_record [NUM_CLONE_soft[i]] ) ), kwargs["CLEMENT_DIR"]+ "/CLEMENT_soft_" + priority + "." + kwargs["IMAGE_FORMAT"], np_vaf, samplename_dict, cluster_soft.includefp_record [NUM_CLONE_soft[i]] , cluster_soft.fp_index_record[NUM_CLONE_soft[i]] , cluster_soft.makeone_index_record[NUM_CLONE_soft[i]], "None" )
            else:
                visualizationsinglesoft.drawfigure_2d (cluster_soft.membership_record [NUM_CLONE_soft[i]], cluster_soft.mixture_record [NUM_CLONE_soft[i]], cluster_soft.membership_p_normalize_record [NUM_CLONE_soft[i]],
                                                                        "CLEMENT_soft : {}".format( round (cluster_soft.likelihood_record [NUM_CLONE_soft[i]]) ), kwargs["CLEMENT_DIR"]+ "/CLEMENT_soft_" + priority + "." + kwargs["IMAGE_FORMAT"], np_vaf, samplename_dict, cluster_soft.includefp_record [NUM_CLONE_soft[i]] , cluster_soft.fp_index_record[NUM_CLONE_soft[i]] , cluster_soft.makeone_index_record[NUM_CLONE_soft[i]], "SVD" )

            subprocess.run (["cp -rf " + kwargs["CLEMENT_DIR"] + "/candidate/clone" + str(NUM_CLONE_soft[i]) +".\(soft\)." + kwargs["IMAGE_FORMAT"]  +" " + kwargs["CLEMENT_DIR"]+ "/CLEMENT_soft_" + priority + "." + kwargs["IMAGE_FORMAT"]], shell = True)
            subprocess.run (["cp -rf " +  kwargs["CLEMENT_DIR"]+ "/CLEMENT_soft_" + priority + "." + kwargs["IMAGE_FORMAT"] + "  " + kwargs["COMBINED_OUTPUT_DIR"]  + "/CLEMENT_soft_" + priority + "." + kwargs["IMAGE_FORMAT"]], shell = True)
        
            print ("\tSoft {} results printed".format(priority))



if kwargs["MODE"] in ["Soft", "Both"]:
    #DECISION : Hard 인지 Soft인지
    print ("\n\n\n======================================= STEP #9. DECISION:  HARD VS SOFT  =======================================")
    if soft_std != 0 :
        print ( "DECISION : {}\t\thard_std : {}\tsoft_std : {}\tratio : {}\tcriteria : < {}".format( DECISION, round(hard_std, 3), round(soft_std, 3), round ( round(soft_std, 3) / round(hard_std, 3), 2) , round (kwargs["DECISION_STANDARD"], 2) ))
    else:
        try:
            print ( "DECISION : {}\t\tmoved_col_mean : {}\tnot_moved_col_mean : {}\tcriteria : > 0.1".format( DECISION, moved_col_mean, not_moved_col_mean ) )
        except:
            print ( "DECISION : {}\t\tmoved_col_list : {}".format( DECISION, moved_col_list) )
    
    if "hard" in DECISION:
        NUM_CLONE_CLEMENT = cluster_hard.mixture_record [NUM_CLONE_hard[0]].shape[1] - int (cluster_hard.includefp_record [ NUM_CLONE_hard[0] ])
        NUM_CHILD_CLEMENT = len (cluster_hard.makeone_index_record [NUM_CLONE_hard[0]])
    elif "soft" in DECISION:
        NUM_CLONE_CLEMENT = cluster_soft.mixture_record [NUM_CLONE_soft[0]].shape[1] - int (cluster_soft.includefp_record [ NUM_CLONE_soft[0] ])
        NUM_CHILD_CLEMENT = len (cluster_soft.makeone_index_record [NUM_CLONE_soft[0]])
        

    subprocess.run ([ "cp -rf " +  kwargs["CLEMENT_DIR"]+ "/CLEMENT_" + DECISION + "." + kwargs["IMAGE_FORMAT"]  + " " + kwargs["COMBINED_OUTPUT_DIR"]  + "/result/CLEMENT_decision." + kwargs["IMAGE_FORMAT"] ], shell = True)
    subprocess.run ([ "cp -rf " +  kwargs["CLEMENT_DIR"]+ "/result/CLEMENT_" + DECISION + ".results.txt  " + kwargs["COMBINED_OUTPUT_DIR"]  + "/result/CLEMENT_decision.results.txt" ], shell = True)
    subprocess.run ([ "cp -rf " +  kwargs["CLEMENT_DIR"]+ "/result/CLEMENT_" + DECISION + ".membership.txt  " + kwargs["COMBINED_OUTPUT_DIR"]  + "/result/CLEMENT_decision.membership.txt" ], shell = True)
    subprocess.run ([ "cp -rf " +  kwargs["CLEMENT_DIR"]+ "/result/CLEMENT_" + DECISION + ".membership_count.txt  " + kwargs["COMBINED_OUTPUT_DIR"]  + "/result/CLEMENT_decision.membership_count.txt" ], shell = True)
    subprocess.run ([ "cp -rf " +  kwargs["CLEMENT_DIR"]+ "/result/CLEMENT_" + DECISION + ".mixture.txt  " + kwargs["COMBINED_OUTPUT_DIR"]  + "/result/CLEMENT_decision.mixture.txt" ], shell = True)
    with open (kwargs["CLEMENT_DIR"]  + "/result/CLEMENT_decision.results.txt", "a", encoding = "utf8") as output_myEM:
        print ( "DECISION\t{}\nhard_std\t{}\nsoft_std\t{}".format( DECISION, round(hard_std, 3), round(soft_std, 3) ), file = output_myEM)
        

        


print ("\n\n\n\n============================== STEP #11.   PYCLONEVI RUNNING ==================================")

PYCLONEVI_START_TIME = datetime.datetime.now()

INPUT_TSV=kwargs["PYCLONEVI_DIR"] + "/input.tsv"
OUTPUT_H5=kwargs["PYCLONEVI_DIR"]  + "/output.h5"
OUTPUT_TSV=kwargs["PYCLONEVI_DIR"]  + "/output.tsv"

subprocess.run (["bash " + SCRIPT_DIR + "/pyclonevi_pipe.sh " + INPUT_TSV + " " + OUTPUT_H5 + " " + OUTPUT_TSV], shell = True)

INPUT_PYCLONEVI_RESULT = OUTPUT_TSV
INPUT_NPVAF = kwargs["NPVAF_DIR"] + "/npvaf.txt"
OUTPUT_FILENAME = kwargs["PYCLONEVI_DIR"]  + "/pyclonevi." + kwargs["IMAGE_FORMAT"]

score_df_pyclonevi, score_pyclonevi, max_score_pyclonevi, membership_pyclonevi, mixture_pyclonevi, sample_dict_PtoA_pyclonevi, sample_dict_AtoP_pyclonevi  = \
    pyclonevisim.main(INPUT_PYCLONEVI_RESULT,  INPUT_NPVAF, OUTPUT_FILENAME, mixture_answer, membership_answer, membership_answer_numerical, kwargs["samplename_dict_CharacterToNum"], kwargs["samplename_dict_NumToCharacter"], **kwargs)


if kwargs["SCORING"] == True:
    Y_index_pyclonevi = result.Yindex ( score_df_pyclonevi )
    print ( "FP 빼고 {} - {}개만 다룬다".format ( kwargs["NUM_MUTATION"], len ( membership_answer_numerical_nofp_index ) ))
    ARI_pyclonevi = result.ARI ( np.array ( [ membership_answer_numerical [i] for i in membership_answer_numerical_nofp_index  ] ) , 
                                                    np.array ( [ membership_pyclonevi [i] for i in membership_answer_numerical_nofp_index] ) )
    print ("\n[ ■  PYCLONEVI RESULTS]\n\npyclone_vi\t{}/{}\nNUM_CLONE\t{}\nARI\t{}\n".format(max_score_pyclonevi, kwargs["NUM_MUTATION"], mixture_pyclonevi.shape[1],  round (ARI_pyclonevi, 2)  ))
    
    #print ("\n(Greedy 방식) score : {}점 / {}점".format(score_pyclonevi, kwargs["NUM_MUTATION"]))
    print ("(모두 다 돌려본 결과) score : {}점 / {}점,  변환표 (A to P)= {}\n".format ( max_score_pyclonevi, kwargs["NUM_MUTATION"], sample_dict_AtoP_pyclonevi))
    print (score_df_pyclonevi)
    

    if (kwargs["FP_RATIO"] != 0) | (kwargs["FP_USEALL"] == "True"):
        # print ("[FP ANALYSIS]")
        answeronly_pyclonevi, intersection_pyclonevi, pyclonevi_only, sensitivity_pyclonevi, PPV_pyclonevi, F1_pyclonevi = result.FPmatrix(score_df_pyclonevi)
        # print ("answer FP {}개 중에 {}개 일치함".format( answeronly_pyclonevi + intersection_pyclonevi ,  intersection_pyclonevi ) )
        # print ("\tanswerFP only : {}\n\tintersection : {}\n\tpycloneviFP only : {}".format( answeronly_pyclonevi, intersection_pyclonevi, pyclonevi_only ))
    else:
        answeronly_pyclonevi, intersection_pyclonevi, pyclonevi_only, sensitivity_pyclonevi, PPV_pyclonevi, F1_pyclonevi = 0, 0, 0, None, None, None

    NUM_CLONE_pyclonevi, NUM_CHILD_pyclonevi = mixture_pyclonevi.shape[1], mixture_pyclonevi.shape[1]
    with open (kwargs["PYCLONEVI_DIR"]  + "/results.txt", "w", encoding = "utf8") as output_pyclonevi:
        print ("NUM_CLONE\t{}\nNUM_CHILD\t{}\nscore\t{}/{}\nARI\t{}\nrunningtime\t{}".
                format(NUM_CLONE_pyclonevi, NUM_CLONE_pyclonevi, max_score_pyclonevi, kwargs["NUM_MUTATION"], round (ARI_pyclonevi, 2) ,  round((datetime.datetime.now() - PYCLONEVI_START_TIME).total_seconds()) ), file = output_pyclonevi)
    subprocess.run (["cp -rf " + kwargs["PYCLONEVI_DIR"] + "/results.txt  " + kwargs["COMBINED_OUTPUT_DIR"]  + "/result/pyclonevi.results.txt"], shell = True)


    pd.DataFrame(membership_pyclonevi).to_csv (kwargs["PYCLONEVI_DIR"] + "/membership.txt", index = False, header= False,  sep = "\t" )
    pd.DataFrame(membership_pyclonevi).to_csv (kwargs["COMBINED_OUTPUT_DIR"]  + "/result/pyclonevi.membership.txt", index = False, header= False,  sep = "\t" )
    pd.DataFrame(mixture_pyclonevi).to_csv (kwargs["PYCLONEVI_DIR"] + "/results.mixture.txt", index = False, header= False,  sep = "\t" )
    pd.DataFrame(mixture_pyclonevi).to_csv (kwargs["COMBINED_OUTPUT_DIR"]  + "/result/pyclonevi.mixture.txt", index = False, header= False,  sep = "\t" )
    pd.DataFrame(score_df_pyclonevi).to_csv (kwargs["PYCLONEVI_DIR"] + "/results.scoredf.txt", index = False, header= True,  sep = "\t")
    pd.DataFrame(score_df_pyclonevi).to_csv (kwargs["COMBINED_OUTPUT_DIR"]  + "/result/pyclonevi.scoredf.txt", index = False, header= True,  sep = "\t" )

    samplename_dict = {k:k for k in range(0, np.max(membership_pyclonevi)+ 1)}
    if kwargs["NUM_BLOCK"] == 1:
        visualizationsingle.drawfigure_1d (membership_pyclonevi, "pyclone_vi", kwargs["PYCLONEVI_DIR"]+ "/pyclonevi." + kwargs["IMAGE_FORMAT"], np_vaf, samplename_dict, False, -1, list( set (membership_pyclonevi)), **kwargs )
    elif kwargs["NUM_BLOCK"] == 2:
        visualizationpair.drawfigure_2d (membership_answer, mixture_answer, membership_pyclonevi, mixture_pyclonevi, score_df_pyclonevi, OUTPUT_FILENAME,  "ANSWER\n", "pyclone_vi\n{}/{}, ARI={}".format(score_pyclonevi, kwargs["NUM_MUTATION"], round (ARI_pyclonevi, 2)), np_vaf, "No",  [], dimensionreduction="None", **kwargs)
    else:
        visualizationsingle.drawfigure_2d (membership_pyclonevi, "pyclone_vi", kwargs["PYCLONEVI_DIR"]+ "/pyclonevi." + kwargs["IMAGE_FORMAT"], np_vaf, samplename_dict, False, -1, "SVD", **kwargs  )

    subprocess.run (["cp -rf " + OUTPUT_FILENAME + "  " + kwargs["COMBINED_OUTPUT_DIR"]  + "/result/pyclonevi." + kwargs["IMAGE_FORMAT"]], shell = True)


else:
    with open (kwargs["PYCLONEVI_DIR"] + "/pyclonevi.results.txt", "w", encoding = "utf8") as output_PYCLONEVI:
        print ("NUM_CLONE\t{}\nNUM_CHILD\t{}\nrunningtime\t{}".
                format(mixture_pyclonevi.shape[1], mixture_pyclonevi.shape[1],   round((datetime.datetime.now() - START_TIME).total_seconds()) ), file = output_PYCLONEVI)

    subprocess.run (["cp -rf " + kwargs["PYCLONEVI_DIR"] + "/pyclonevi.results.txt  " + kwargs["COMBINED_OUTPUT_DIR"]  + "/result/pyclonevi.results.txt"], shell = True)

    pd.DataFrame(membership_pyclonevi).to_csv (kwargs["PYCLONEVI_DIR"] + "/membership.txt", index = False, header= False,  sep = "\t" )
    pd.DataFrame(membership_pyclonevi).to_csv (kwargs["COMBINED_OUTPUT_DIR"]  + "/result/pyclonevi.membership.txt", index = False, header= False,  sep = "\t" )
    pd.DataFrame(mixture_pyclonevi).to_csv (kwargs["PYCLONEVI_DIR"] + "/results.mixture.txt", index = False, header= False,  sep = "\t" )
    pd.DataFrame(mixture_pyclonevi).to_csv (kwargs["COMBINED_OUTPUT_DIR"]  + "/result/pyclonevi.mixture.txt", index = False, header= False,  sep = "\t" )

    samplename_dict = {k:k for k in range(0, np.max(membership_pyclonevi)+ 1)}
    if kwargs["NUM_BLOCK"] == 1:
        samplename_dict = {k:"clone {}".format(k) for k in range(0, np.max(membership_pyclonevi)+ 1)}
        visualizationsingle.drawfigure_1d (membership_pyclonevi, "pyclone_vi", kwargs["PYCLONEVI_DIR"]+ "/pyclonevi." + kwargs["IMAGE_FORMAT"], np_vaf, samplename_dict, False, -1 , list( set (membership_pyclonevi)), **kwargs  )
    elif kwargs["NUM_BLOCK"] == 2:
        visualizationsingle.drawfigure_2d (membership_pyclonevi, "pyclone_vi", kwargs["PYCLONEVI_DIR"]+ "/pyclonevi." + kwargs["IMAGE_FORMAT"], np_vaf, samplename_dict, False, -1, "None", **kwargs  )
    else:
        visualizationsingle.drawfigure_2d (membership_pyclonevi, "pyclone_vi", kwargs["PYCLONEVI_DIR"]+ "/pyclonevi." + kwargs["IMAGE_FORMAT"], np_vaf, samplename_dict, False, -1, "SVD", **kwargs  )

    subprocess.run (["cp -rf " +  kwargs["PYCLONEVI_DIR"]+ "/pyclonevi." + kwargs["IMAGE_FORMAT"] + " "  + kwargs["COMBINED_OUTPUT_DIR"]  + "/result/pyclonevi." + kwargs["IMAGE_FORMAT"]], shell = True)

print ("\n현재 시각 : {}h:{}m:{}s     (걸린 시간 : {})".format(time.localtime().tm_hour, time.localtime().tm_min, round(time.localtime().tm_sec), datetime.datetime.now() - PYCLONEVI_START_TIME ))











print("\n\n\n\n================================ STEP #12.   SCICLONE RUNNING ==================================")

SCICLONE_START_TIME = datetime.datetime.now()
print("\nNOW SCICLONE IS STARTED  :  {}h:{}m:{}s\n\n".format( time.localtime().tm_hour, time.localtime().tm_min, round(time.localtime().tm_sec)))

INPUT_First = kwargs["SCICLONE_DIR"] + "/block0.dat"
INPUT_Second = kwargs["SCICLONE_DIR"] + "/block1.dat"
INPUT_Third = kwargs["SCICLONE_DIR"] + "/block2.dat"

if NUM_BLOCK == 3:
    subprocess.run(["bash " + SCRIPT_DIR + "/sciclone_pipe_3D.sh " + SCRIPT_DIR + " " + INPUT_First + " " + INPUT_Second + " " + INPUT_Third + " " + kwargs["SCICLONE_DIR"]], shell=True)
elif NUM_BLOCK == 2:
    subprocess.run(["bash " + SCRIPT_DIR + "/sciclone_pipe_2D.sh " + SCRIPT_DIR + " " + INPUT_First + " " + INPUT_Second + " " + kwargs["SCICLONE_DIR"]], shell=True)
elif NUM_BLOCK == 1:
    subprocess.run(["bash " + SCRIPT_DIR + "/sciclone_pipe_1D.sh " + SCRIPT_DIR + " " + INPUT_First + " " + kwargs["SCICLONE_DIR"]], shell=True)


INPUT_SCICLONE_RESULT = kwargs["SCICLONE_DIR"] + "/results.tsv"
INPUT_NPVAF = kwargs["NPVAF_DIR"] + "/npvaf.txt"
OUTPUT_FILENAME = kwargs["SCICLONE_DIR"] + "/sciclone." + kwargs["IMAGE_FORMAT"]

score_df_sciclone, score_sciclone, max_score_sciclone, membership_sciclone, mixture_sciclone, sample_dict_PtoA_sciclone, sample_dict_AtoP_sciclone = \
    sciclonesim.main(INPUT_SCICLONE_RESULT,  INPUT_NPVAF, OUTPUT_FILENAME, mixture_answer, membership_answer,  membership_answer_numerical, kwargs["samplename_dict_CharacterToNum"], kwargs["samplename_dict_NumToCharacter"], **kwargs)


if kwargs["SCORING"] == True:
    Y_index_sciclone = result.Yindex(score_df_sciclone)
    ARI_sciclone = result.ARI(np.array([membership_answer_numerical[i] for i in membership_answer_numerical_nofp_index]),
                              np.array([membership_sciclone[i] for i in membership_answer_numerical_nofp_index]))

    print("\n[■  SCICLONE RESULTS]\n\nSciClone\t{}/{}\nNUM_CLONE\t{}\nARI\t{}\n".format(max_score_sciclone,
          kwargs["NUM_MUTATION"], mixture_sciclone.shape[1],  round (ARI_sciclone, 2)  ))

    #print("\n(Greedy 방식) score : {}점 / {}점".format(score_sciclone,kwargs["NUM_MUTATION"]))
    print("(모두 다 돌려본 결과) score : {}점 / {}점,  변환표 (A to P)= {}\n".format(max_score_sciclone, kwargs["NUM_MUTATION"], sample_dict_AtoP_sciclone))
    print(score_df_sciclone)

    if (kwargs["FP_RATIO"] != 0) | (kwargs["FP_USEALL"] == "True"):
        # print("[FP ANALYSIS]")
        answeronly_sciclone, intersection_sciclone, sciclone_only, sensitivity_sciclone, PPV_sciclone, F1_sciclone = result.FPmatrix(score_df_sciclone)
        # print("answer FP {}개 중에 {}개 일치함".format(answeronly_sciclone + intersection_sciclone,  intersection_sciclone))
        # print("\tanswerFP only : {}\n\tintersection : {}\n\tscicloneFP only : {}".format(answeronly_sciclone, intersection_sciclone, sciclone_only))
        print("")
    else:
        answeronly_sciclone, intersection_sciclone, sciclone_only, sensitivity_sciclone, PPV_sciclone, F1_sciclone = 0, 0, 0, None, None, None

    NUM_CLONE_sciclone, NUM_CHILD_sciclone = mixture_sciclone.shape[1], mixture_sciclone.shape[1]
    with open(kwargs["SCICLONE_DIR"] + "/results.txt", "w", encoding="utf8") as output_sciclone:
        print("NUM_CLONE\t{}\nNUM_CHILD\t{}\nscore\t{}/{}\nARI\t{}\nrunningtime\t{}".
              format(NUM_CLONE_sciclone, NUM_CHILD_sciclone, max_score_sciclone, kwargs["NUM_MUTATION"],  ARI_sciclone,  round((datetime.datetime.now() - SCICLONE_START_TIME).total_seconds())), file=output_sciclone)
    subprocess.run(["cp -rf " + kwargs["SCICLONE_DIR"] + "/results.txt  " + kwargs["COMBINED_OUTPUT_DIR"] + "/result/sciclone.results.txt"], shell=True)

    pd.DataFrame(membership_sciclone).to_csv(kwargs["SCICLONE_DIR"] + "/membership.txt", index=False, header=False,  sep="\t")
    pd.DataFrame(membership_sciclone).to_csv(kwargs["COMBINED_OUTPUT_DIR"] + "/result/sciclone.membership.txt", index=False, header=False,  sep="\t")
    pd.DataFrame(mixture_sciclone).to_csv(kwargs["SCICLONE_DIR"] + "/results.mixture.txt", index=False, header=False,  sep="\t")
    pd.DataFrame(mixture_sciclone).to_csv(kwargs["COMBINED_OUTPUT_DIR"] + "/result/sciclone.mixture.txt", index=False, header=False,  sep="\t")
    pd.DataFrame(score_df_sciclone).to_csv(kwargs["SCICLONE_DIR"] + "/results.scoredf.txt", index=False, header=True,  sep="\t")
    pd.DataFrame(score_df_sciclone).to_csv(kwargs["COMBINED_OUTPUT_DIR"] + "/result/sciclone.scoredf.txt", index=False, header=True,  sep="\t")

    samplename_dict = {k: k for k in range(0, np.max(membership_sciclone) + 1)}
    if kwargs["NUM_BLOCK"] == 1:
        visualizationsingle.drawfigure_1d(membership_sciclone, "SciClone", kwargs["SCICLONE_DIR"] + "/sciclone." + kwargs["IMAGE_FORMAT"], np_vaf, samplename_dict, False, -1 , list (set (membership_sciclone)), **kwargs   )
    elif kwargs["NUM_BLOCK"] == 2:
        visualizationpair.drawfigure_2d(membership_answer, mixture_answer, membership_sciclone, mixture_sciclone, score_df_sciclone, OUTPUT_FILENAME,  "ANSWER\n",
                                            "SciClone\n{}/{}, ARI={}".format(max_score_sciclone, kwargs["NUM_MUTATION"], round(ARI_sciclone, 2)), np_vaf, "No",  [],  dimensionreduction="None", **kwargs )
    elif kwargs["NUM_BLOCK"] == 2:
        visualizationpair.drawfigure_2d(membership_answer, mixture_answer, membership_sciclone, mixture_sciclone, score_df_sciclone, OUTPUT_FILENAME,  "ANSWER\n",
                                            "SciClone\n{}/{}, ARI={}".format(max_score_sciclone, kwargs["NUM_MUTATION"], round(ARI_sciclone, 2)), np_vaf, "No",  [],  dimensionreduction="SVD", **kwargs )

    subprocess.run(["cp -rf " + OUTPUT_FILENAME + "  " +  kwargs["COMBINED_OUTPUT_DIR"] + "/result/sciclone." + kwargs["IMAGE_FORMAT"]], shell=True)


else:
    with open(kwargs["SCICLONE_DIR"] + "/sciclone.results.txt", "w", encoding="utf8") as output_sciclone:
        print("NUM_CLONE\t{}\nNUM_CHILD\t{}\nrunningtime\t{}".
              format(mixture_sciclone.shape[1], mixture_sciclone.shape[1],   round((datetime.datetime.now() - START_TIME).total_seconds())), file=output_sciclone)
    subprocess.run(["cp -rf " + kwargs["SCICLONE_DIR"] + "/sciclone.results.txt  " + kwargs["COMBINED_OUTPUT_DIR"] + "/result/sciclone.results.txt"], shell=True)

    pd.DataFrame(membership_sciclone).to_csv(  kwargs["SCICLONE_DIR"] + "/membership.txt", index=False, header=False,  sep="\t" )
    pd.DataFrame(membership_sciclone).to_csv(  kwargs["COMBINED_OUTPUT_DIR"] + "/result/sciclone.membership.txt", index=False, header=False,  sep="\t" )
    pd.DataFrame(mixture_sciclone).to_csv(  kwargs["SCICLONE_DIR"] + "/results.mixture.txt", index=False, header=False,  sep="\t" )
    pd.DataFrame(mixture_sciclone).to_csv(  kwargs["COMBINED_OUTPUT_DIR"] + "/result/sciclone.mixture.txt", index=False, header=False,  sep="\t" )

    samplename_dict = {k: k for k in range(0, np.max(membership_sciclone) + 1)}
    if kwargs["NUM_BLOCK"] == 1:
        samplename_dict = {k: "clone {}".format(k) for k in range(0, np.max(membership_sciclone) + 1)}
        visualizationsingle.drawfigure_1d( membership_sciclone, "SciClone", kwargs["SCICLONE_DIR"] + "/sciclone." + kwargs["IMAGE_FORMAT"], np_vaf, samplename_dict, False, -1, list (set (membership_sciclone)), **kwargs  )
    elif kwargs["NUM_BLOCK"] == 2:
        visualizationsingle.drawfigure_2d( membership_sciclone, "SciClone", kwargs["SCICLONE_DIR"] + "/sciclone." + kwargs["IMAGE_FORMAT"], np_vaf, samplename_dict, False, -1, "None", **kwargs )
    else:
        visualizationsingle.drawfigure_2d( membership_sciclone, "SciClone", kwargs["SCICLONE_DIR"] + "/sciclone." + kwargs["IMAGE_FORMAT"], np_vaf, samplename_dict, False, -1, "SVD", **kwargs )

    subprocess.run(["cp -rf " + kwargs["SCICLONE_DIR"] + "/sciclone." + kwargs["IMAGE_FORMAT"] + " "  +  kwargs["COMBINED_OUTPUT_DIR"] + "/result/sciclone." + kwargs["IMAGE_FORMAT"]], shell=True)


print("\n현재 시각 : {}h:{}m:{}s     (걸린 시간 : {})".format(time.localtime().tm_hour, time.localtime().tm_min, round(time.localtime().tm_sec), datetime.datetime.now() - SCICLONE_START_TIME))










print("\n\n\n\n================================ STEP #13.   QUANTUMCLONE RUNNING ==================================")

QUANTUMCLONE_START_TIME = datetime.datetime.now()
print("\nNOW QUANTUMCLONE IS STARTED  :  {}h:{}m:{}s\n\n".format(time.localtime().tm_hour, time.localtime().tm_min, round(time.localtime().tm_sec)))

INPUT_First = kwargs["QUANTUMCLONE_DIR"] + "/block0.dat"
INPUT_Second = kwargs["QUANTUMCLONE_DIR"] + "/block1.dat"
INPUT_Third = kwargs["QUANTUMCLONE_DIR"] + "/block2.dat"

if NUM_BLOCK == 3:
    subprocess.run(["bash " + SCRIPT_DIR + "/qc_pipe_3D.sh " + SCRIPT_DIR + " " + INPUT_First + " " + INPUT_Second + " " + INPUT_Third + " " + kwargs["QUANTUMCLONE_DIR"]], shell=True)
elif NUM_BLOCK == 2:
    subprocess.run(["bash " + SCRIPT_DIR + "/qc_pipe_2D.sh " + SCRIPT_DIR + " " + INPUT_First + " " + INPUT_Second + " " + kwargs["QUANTUMCLONE_DIR"]], shell=True)
elif NUM_BLOCK == 1:
    subprocess.run(["bash " + SCRIPT_DIR + "/qc_pipe_1D.sh " + SCRIPT_DIR + " "+ INPUT_First + " " + kwargs["QUANTUMCLONE_DIR"]], shell=True)

INPUT_QUANTUMCLONE_RESULT = kwargs["QUANTUMCLONE_DIR"]
INPUT_NPVAF = kwargs["NPVAF_DIR"] + "/npvaf.txt"
OUTPUT_FILENAME = kwargs["QUANTUMCLONE_DIR"] + "/quantumclone." + kwargs["IMAGE_FORMAT"]

score_df_quantumclone, score_quantumclone, max_score_quantumclone, membership_quantumclone, mixture_quantumclone, sample_dict_PtoA_quantumclone, sample_dict_AtoP_quantumclone = \
    quantumclonesim.main(INPUT_QUANTUMCLONE_RESULT,  INPUT_NPVAF, OUTPUT_FILENAME, mixture_answer, membership_answer, membership_answer_numerical, kwargs["samplename_dict_CharacterToNum"], kwargs["samplename_dict_NumToCharacter"], **kwargs)

if kwargs["SCORING"] == True:
    Y_index_quantumclone = result.Yindex(score_df_quantumclone)
    ARI_quantumclone = result.ARI(np.array([membership_answer_numerical[i] for i in membership_answer_numerical_nofp_index]),
                                  np.array([membership_quantumclone[i] for i in membership_answer_numerical_nofp_index]))

    print("\n[ ■ QUANTUMCLONE RESULTS]\n\nQuantumClone\t{}/{}\nNUM_CLONE\t{}\nARI\t{}\n{}\n"
          .format(score_quantumclone, kwargs["NUM_MUTATION"], mixture_quantumclone.shape[1],  ARI_quantumclone, list(mixture_quantumclone)))

    #print("\n(Greedy 방식) score : {}점 / {}점".format(score_quantumclone, kwargs["NUM_MUTATION"]))
    print("(모두 다 돌려본 결과) score : {}점 / {}점,  변환표 (A to P)= {}\n".format(max_score_quantumclone, kwargs["NUM_MUTATION"], sample_dict_AtoP_quantumclone))
    print(score_df_quantumclone)

    if (kwargs["FP_RATIO"] != 0) | (kwargs["FP_USEALL"] == "True"):
        #print("[FP ANALYSIS]")
        answeronly_quantumclone, intersection_quantumclone, quantumclone_only, sensitivity_quantumclone, PPV_quantumclone, F1_quantumclone = result.FPmatrix(
            score_df_quantumclone)
        # print("answer FP {}개 중에 {}개 일치함".format(answeronly_quantumclone + intersection_quantumclone,  intersection_quantumclone))
        # print("\tanswerFP only : {}\n\tintersection : {}\n\tquantumcloneFP only : {}".format(answeronly_quantumclone, intersection_quantumclone, quantumclone_only))
    else:
        answeronly_quantumclone, intersection_quantumclone, quantumclone_only, sensitivity_quantumclone, PPV_quantumclone, F1_quantumclone = 0, 0, 0, None, None, None

    NUM_CLONE_quantumclone, NUM_CHILD_quantumclone = mixture_quantumclone.shape[1], mixture_quantumclone.shape[1]
    with open(kwargs["QUANTUMCLONE_DIR"] + "/results.txt", "w", encoding="utf8") as output_quantumclone:
        print("NUM_CLONE\t{}\nNUM_CHILD\t{}\nscore\t{}/{}\nARI\t{}\nrunningtime\t{}".
              format(NUM_CLONE_quantumclone, NUM_CHILD_quantumclone, score_quantumclone, kwargs["NUM_MUTATION"], ARI_quantumclone,   round((datetime.datetime.now() - QUANTUMCLONE_START_TIME).total_seconds())), file=output_quantumclone)

    subprocess.run(["cp -rf " + kwargs["QUANTUMCLONE_DIR"] + "/results.txt  " + kwargs["COMBINED_OUTPUT_DIR"] + "/result/quantumclone.results.txt"], shell=True)

    pd.DataFrame(membership_quantumclone).to_csv(kwargs["QUANTUMCLONE_DIR"] + "/membership.txt", index=False, header=False,  sep="\t")
    pd.DataFrame(membership_quantumclone).to_csv(kwargs["COMBINED_OUTPUT_DIR"] + "/result/quantumclone.membership.txt", index=False, header=False,  sep="\t")
    pd.DataFrame(mixture_quantumclone).to_csv(kwargs["QUANTUMCLONE_DIR"] + "/results.mixture.txt", index=False, header=False,  sep="\t")
    pd.DataFrame(mixture_quantumclone).to_csv(kwargs["COMBINED_OUTPUT_DIR"] + "/result/quantumclone.mixture.txt", index=False, header=False,  sep="\t")
    pd.DataFrame(score_df_quantumclone).to_csv(kwargs["QUANTUMCLONE_DIR"] + "/results.scoredf.txt", index=False, header=True,  sep="\t")
    pd.DataFrame(score_df_quantumclone).to_csv(kwargs["COMBINED_OUTPUT_DIR"] + "/result/quantumclone.scoredf.txt", index=False, header=True,  sep="\t")

    samplename_dict = {k: k for k in range( 0, np.max(membership_quantumclone) + 1)}
    if kwargs["NUM_BLOCK"] == 1:
        visualizationsingle.drawfigure_1d(membership_quantumclone, "QuantumClone", kwargs["QUANTUMCLONE_DIR"] + "/quantumclone." + kwargs["IMAGE_FORMAT"], np_vaf, samplename_dict, False, -1, list ( set (membership_quantumclone)), **kwargs ) 
    elif kwargs["NUM_BLOCK"] == 2:
        visualizationpair.drawfigure_2d(membership_answer, mixture_answer, membership_quantumclone, mixture_quantumclone, score_df_quantumclone, OUTPUT_FILENAME,  "ANSWER\n",
                                        "QuantumClone\n{}/{}, ARI={}".format(score_quantumclone, kwargs["NUM_MUTATION"], round(ARI_quantumclone, 2)), np_vaf, "No",  [],  dimensionreduction="None", **kwargs )
    elif kwargs["NUM_BLOCK"] >= 3:
        visualizationpair.drawfigure_2d(membership_answer, mixture_answer, membership_quantumclone, mixture_quantumclone, score_df_quantumclone, OUTPUT_FILENAME,  "ANSWER\n",
                                        "QuantumClone\n{}/{}, ARI={}".format(score_quantumclone, kwargs["NUM_MUTATION"], round(ARI_quantumclone, 2)), np_vaf, "No",  [],  dimensionreduction="SVD", **kwargs )

    subprocess.run(["cp -rf " + OUTPUT_FILENAME + "  " + kwargs["COMBINED_OUTPUT_DIR"] + "/result/quantumclone." + kwargs["IMAGE_FORMAT"]], shell=True)


else:
    with open(kwargs["QUANTUMCLONE_DIR"] + "/quantumclone.results.txt", "w", encoding="utf8") as output_quantumclone:
        print("NUM_CLONE\t{}\nNUM_CHILD\t{}\nrunningtime\t{}".
              format(mixture_quantumclone.shape[1], mixture_quantumclone.shape[1],   round((datetime.datetime.now() - START_TIME).total_seconds())), file=output_quantumclone)
    subprocess.run(["cp -rf " + kwargs["QUANTUMCLONE_DIR"] + "/quantumclone.results.txt  " +  kwargs["COMBINED_OUTPUT_DIR"] + "/result/quantumclone.results.txt"], shell=True)

    pd.DataFrame(membership_quantumclone).to_csv( kwargs["QUANTUMCLONE_DIR"] + "/membership.txt", index=False, header=False,  sep="\t")
    pd.DataFrame(membership_quantumclone).to_csv( kwargs["COMBINED_OUTPUT_DIR"] + "/result/quantumclone.membership.txt", index=False, header=False,  sep="\t")
    pd.DataFrame(mixture_quantumclone).to_csv( kwargs["QUANTUMCLONE_DIR"] + "/results.mixture.txt", index=False, header=False,  sep="\t")
    pd.DataFrame(mixture_quantumclone).to_csv( kwargs["COMBINED_OUTPUT_DIR"] + "/result/quantumclone.mixture.txt", index=False, header=False,  sep="\t")

    samplename_dict = {k: k for k in range( 0, np.max(membership_quantumclone) + 1)}
    if kwargs["NUM_BLOCK"] == 1:
        samplename_dict = {k: "clone {}".format(k) for k in range( 0, np.max(membership_quantumclone) + 1)}
        visualizationsingle.drawfigure_1d(membership_quantumclone, "QuantumClone", kwargs["QUANTUMCLONE_DIR"] + "/quantumclone." + kwargs["IMAGE_FORMAT"], np_vaf, samplename_dict, False, -1, list ( set (membership_quantumclone)), **kwargs )
    elif kwargs["NUM_BLOCK"] == 2:
        visualizationsingle.drawfigure_2d(membership_quantumclone, "QuantumClone", kwargs["QUANTUMCLONE_DIR"] + "/quantumclone." + kwargs["IMAGE_FORMAT"], np_vaf, samplename_dict, False, -1, "None", **kwargs )
    else:
        visualizationsingle.drawfigure_2d(membership_quantumclone, "QuantumClone", kwargs["QUANTUMCLONE_DIR"] + "/quantumclone." + kwargs["IMAGE_FORMAT"], np_vaf, samplename_dict, False, -1, "SVD", **kwargs )

    subprocess.run(["cp -rf " + kwargs["QUANTUMCLONE_DIR"] + "/quantumclone." + kwargs["IMAGE_FORMAT"] + " "  + kwargs["COMBINED_OUTPUT_DIR"] + "/result/quantumclone." + kwargs["IMAGE_FORMAT"]], shell=True)



if kwargs["SCORING"] == True:    
    with open (kwargs["COMBINED_OUTPUT_DIR"] + "/★results_score.txt", "w", encoding = "utf8") as results_score:
        print ( "ANSWER_NUM_CLONE = {}  (FP : {})\n".format ( len (kwargs["samplename_dict_NumToCharacter"].keys()), bool("FP" in kwargs["samplename_dict_NumToCharacter"].values()) ) , file = results_score)
        print ("CLEMENT: {}/{}\tNUM_CLONE = {}\tNUM_CHILD = {}".format( max_score_CLEMENT, kwargs["NUM_MUTATION"], NUM_CLONE_CLEMENT, NUM_CHILD_CLEMENT), file = results_score)
        print ("simpleKmeans_elbow: {}/{}\tNUM_CLONE = {}".format(simpleK.elbow_K_score, kwargs["NUM_MUTATION"], simpleK.elbow_K), file = results_score)
        print ("simpleKmeans_silhouette: {}/{}\tNUM_CLONE = {}".format(simpleK.silhouette_K_score, kwargs["NUM_MUTATION"], simpleK.silhouette_K), file = results_score)
        print ("simpleKmeans_gap*: {}/{}\tNUM_CLONE = {}".format(simpleK.gap_K_score, kwargs["NUM_MUTATION"], simpleK.gap_K), file = results_score)
        print ("pyclone_vi : {}/{}\tNUM_CLONE = {}".format(max_score_pyclonevi, kwargs["NUM_MUTATION"], NUM_CLONE_pyclonevi ), file = results_score)
        print ("sciclone: {}/{}\tNUM_CLONE = {}".format(max_score_sciclone, kwargs["NUM_MUTATION"], NUM_CLONE_sciclone ), file = results_score)
        print ("quantumclone : {}/{}\tNUM_CLONE = {}".format(max_score_quantumclone, kwargs["NUM_MUTATION"], NUM_CLONE_quantumclone ), file = results_score)

print("\n현재 시각 : {}h:{}m:{}s     (걸린 시간 : {})".format(time.localtime().tm_hour, time.localtime().tm_min, round(time.localtime().tm_sec), datetime.datetime.now() - QUANTUMCLONE_START_TIME))
