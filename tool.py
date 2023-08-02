import argparse, datetime, time, subprocess
from collections import Counter
import filetype, datapreparation220712, comb, extract, scoring, boundaryclone, graph, phylogeny, scoring, result, isparent, fppick, onemore, os
import pyclonesim, pyclonevisim, sciclonesim
import EMhard, EMsoft
import visualizationsingle, visualizationpair, visualizationsinglesoft, visualizationpairsoft
import numpy as np
import pandas as pd
import scipy
from kneed import KneeLocator
pd.options.mode.chained_assignment = None

kwargs = {}

parser = argparse.ArgumentParser(description='The below is usage direction.')
parser.add_argument('--INPUT_TSV', type = str, default="/data/project/Alzheimer/EM_cluster/EM_input/MRS_2_sample/M1-5_M1-6_input.txt", help = "Input data whether TSV or VCF. The tool automatically detects the number of samples")
parser.add_argument('--NPVAF_DIR', default = None, help = "Directory where selected datasets are")
parser.add_argument('--MYEM_DIR', default = None, help = "Directory where input and output of MY EM deposits")
parser.add_argument('--SCICLONE_DIR', default = None, help = "Directory where input and output of SCICLONE deposits")
parser.add_argument('--PYCLONE_DIR', default = None, help = "Directory where input and output of PYCLONE deposits")
parser.add_argument('--PYCLONEVI_DIR', default = None, help = "Directory where input and output of PYCLONEVI deposits")
parser.add_argument('--COMBINED_OUTPUT_DIR', default = None, help = "Directory where input and output of MYEM, PYCLONEVI, SCICLONE deposits")
parser.add_argument('--NUM_CLONE_TRIAL_START', type = int, default=2, help = "Minimum number of expected clusters (initation of K)")
parser.add_argument('--NUM_CLONE_TRIAL_END', type = int, default=6, choices=range(1, 11), help = "Maximum number of expected clusters (termination of K)")
parser.add_argument('--NUM_CLONE_TRIAL_FORCE', type = int, default=4, help = "Designate the number of expected clusters by force")
parser.add_argument('--RANDOM_PICK', type = int, default=500, help = "The number of mutations that are randomly selected in each trials")
parser.add_argument('--AXIS_RATIO', default=0,  type = float, help = "The fraction of the mutations not shared at least one sample")
parser.add_argument('--PARENT_RATIO', default=0,  type = float, help = "The fraction of parent clone mutations. If this values is designated, do not set NUM_PARENT")
parser.add_argument('--NUM_PARENT',  default=0, type = int, help = "The fraction of parent clones being inserted by large order. If this values is designated, do not set PARENT_RATIO")
parser.add_argument('--FP_RATIO', default=0,  type = float, help = "The fraction of false positive mutations regardless of all samples")
parser.add_argument('--FP_2D', default="False", choices = ["True", "False"], help = "True : extract ALL FPs,   False : Do not extact FPs")
parser.add_argument('--TRIAL_NO', default=10, type = int, choices=range(1, 21), help = "Trial number in each candidate cluster number. DO NOT recommend over 15")
parser.add_argument('--DEPTH_CUTOFF', default=100, type = int, help = "The mutation of which depth below this values is abandoned")
parser.add_argument('--VERBOSE', type = int, choices = [0, 1, 2], default=0)
parser.add_argument('--ELBOW_S', type = float, default=3)
parser.add_argument('--GAUSSIAN_SD', type = float, default=1)
parser.add_argument('--MIN_CLUSTER_SIZE', type = int, default=10)
parser.add_argument('--OUTLIER_STANDARD', choices = ["looser", "stricter"], default="looser")
parser.add_argument('--OUTLIER_METHOD', default="Distance", choices = ["Distance", "ML"], help = "The method of delineating outliers after 1st soft clustering (Distance or ML)")
parser.add_argument('--KMEANS_CLUSTERNO',  type = int , default=15, choices = range(8,20), help = "Number of initial K-means cluster")
parser.add_argument('--RANDOM_SEED', type = int, default=1, help = "random_seed for regular random sampling")


args = parser.parse_args()
kwargs["INPUT_TSV"] = args.INPUT_TSV

INPUT_TSV = kwargs["INPUT_TSV"]
INPUT_FILETYPE, NUM_BLOCK = filetype.main (INPUT_TSV)
kwargs["NUM_BLOCK_INPUT"], kwargs["NUM_BLOCK"] = NUM_BLOCK, NUM_BLOCK
SAMPLENAME = INPUT_TSV.split("/")[-1].split(".")[0]     # 'M1-5_M1-8_input'

if args.NPVAF_DIR == None:
    kwargs["NPVAF_DIR"] = "/data/project/Alzheimer/YSscript/EM_MRS/data/npvaf/" + SAMPLENAME
if args.MYEM_DIR == None:
    kwargs["MYEM_DIR"] = "/data/project/Alzheimer/YSscript/EM_MRS/data/MyEM/" + SAMPLENAME
if args.SCICLONE_DIR == None:
    kwargs["SCICLONE_DIR"] = "/data/project/Alzheimer/YSscript/EM_MRS/data/sciclone/" + SAMPLENAME 
if args.PYCLONE_DIR == None:
    kwargs["PYCLONE_DIR"] = "/data/project/Alzheimer/YSscript/EM_MRS/data/pyclone/" + SAMPLENAME 
if args.PYCLONEVI_DIR == None:
    kwargs["PYCLONEVI_DIR"] = "/data/project/Alzheimer/YSscript/EM_MRS/data/pyclone-vi/" + SAMPLENAME 
if args.COMBINED_OUTPUT_DIR == None:
    kwargs["COMBINED_OUTPUT_DIR"] = "/data/project/Alzheimer/YSscript/EM_MRS/data/combinedoutput/" + SAMPLENAME 


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
kwargs["ELBOW_S"] = float(args.ELBOW_S)
kwargs["GAUSSIAN_SD"] = float(args.GAUSSIAN_SD)
kwargs["MIN_CLUSTER_SIZE"] = int(args.MIN_CLUSTER_SIZE)
kwargs["OUTLIER_STANDARD"] = args.OUTLIER_STANDARD
kwargs["OUTLIER_METHOD"] = args.OUTLIER_METHOD
kwargs["KMEANS_CLUSTERNO"] = args.KMEANS_CLUSTERNO
kwargs["RANDOM_SEED"] = int(args.RANDOM_SEED)

NUM_MUTATION = kwargs["RANDOM_PICK"]


START_TIME = datetime.datetime.now()
#print ("\nNOW RUNNING IS STARTED  {}  : \n\n\t{}\n\n".format(START_TIME  , kwargs))
print ("\nNOW RUNNING IS STARTED  :  {}h:{}m:{}s\n\n".format( time.localtime().tm_hour, time.localtime().tm_min, round(time.localtime().tm_sec) ))


print ("=============== STEP #1.   DATA EXTRACTION FROM THE ANSWER SET  ===============")
os.system ("mkdir -p " +  kwargs["NPVAF_DIR"])
os.system ("mkdir -p " +  kwargs["MYEM_DIR"])
os.system ("rm -rf " +  kwargs["MYEM_DIR"] + "/trial")
os.system ("mkdir -p " +  kwargs["MYEM_DIR"] + "/trial")
os.system ("mkdir -p " + kwargs["SCICLONE_DIR"] )
os.system ("mkdir -p " + kwargs["PYCLONE_DIR"] )
os.system ("mkdir -p " + kwargs["PYCLONEVI_DIR"] )
os.system ("mkdir -p " + kwargs["COMBINED_OUTPUT_DIR"] )

inputdf, df, np_vaf, membership_answer, mixture_answer,  mutation_id, samplename_dict_input  = datapreparation220712.main(**kwargs)
if type(inputdf) != type(False):
    samplename_dict_input_rev = {v: k for k, v in samplename_dict_input.items()}
    visualizationsingle.drawfigure_2d (membership_answer, "NotSave", np_vaf, samplename_dict_input, "")


print ("\n\n=============== STEP #2.   EXPLORATIVE CLUSTERING (DETECT WHETHER THE PARENT CLONE IS)  ===============")
#2-1. EM 돌리기  (원하는 df만 넣어주면 가능)     (elbow,max)+(normal, outlier)    # outlier라면 뭐든지 +1 해서 출력된다
kwargs["ELBOW_S"] = 1; kwargs["GAUSSIAN_SD"] = 2; kwargs["TRIAL_NO"] = 5;  kwargs["MIN_CLUSTER_SIZE"] = 15; kwargs["OUTLIER_STANDARD"] = "looser"
method = "gap+normal"; adjustment = "Half"
NUM_CLONE_NOADJUST, membership, mixture, includeoutlier_inside = EMhard.main(df, method, adjustment, **kwargs)    


#2-2. 정답set과 점수 맞춰보고 색깔 맞추기  (Hard clustering)
score_df, score = scoring.mixturebased(mixture_answer, mixture, membership_answer, membership, samplename_dict_input, samplename_dict_input_rev, includeoutlier_inside, **kwargs)
print ("\nscore : {}점 / {}점".format(score, kwargs["RANDOM_PICK"]))

#2-3. visualizationpair
if kwargs["NUM_BLOCK"] >= 2:
    visualizationpair.drawfigure_2d (membership_answer, mixture_answer, membership, np.round(mixture, 2), score_df, kwargs["MYEM_DIR"]+ "/2.hard.jpg" ,  "ANSWER : MRS", "MY ALGORITHM ", np_vaf, includeoutlier_inside,  dimensionreduction="None")


#2-4. Parent가 있는지 확인해서 ISPARENT 값을 얻기
ISPARENT, CANDIDATE_CHILD_index, CANDIDATE_CHILD_mixture, CANDIDATE_parent_index = isparent.main(membership, mixture, includeoutlier_inside, **kwargs)
print ("\n가장 그럴듯한 CHILD 조합 ({}) : sum = {}".format(CANDIDATE_CHILD_index, CANDIDATE_CHILD_mixture))
if ISPARENT == True:
    print ("{}번째 clone ( {} )이 2개 이상의 child보다 크기에 parent 가능성이 있다".format(CANDIDATE_parent_index, mixture[:, CANDIDATE_parent_index]) )
else:
    print ("Parent clone은 없어 보인다")

print ("\n현재 시각 : {}h:{}m:{}s    (걸린 시간 : {})".format(time.localtime().tm_hour, time.localtime().tm_min, round(time.localtime().tm_sec), datetime.datetime.now() - START_TIME ))


ISPARENT = True
NUM_CLONE_NOADJUST=6


print ("\n\n============ STEP #3.   SOFT CLUSTERING  (OUTLIER CUTOFF DIFFERS BY EXISTENCE OF THE PARENT CLONE) ============")
#3-1. EMsoft 돌리기  (원하는 df만 넣어주면 가능)     (elbow,max,force)+(normal, outlier)    # includeoutlier == Yes라면 +1 해서 출력된다
kwargs["NUM_CLONE_TRIAL_END"] = np.max([6, NUM_CLONE_NOADJUST]); kwargs["MIN_CLUSTER_SIZE"] = int(NUM_MUTATION / 50); 
method = "gap+outlier"; adjustment = "Half"
if ISPARENT == True:       # Outlier 를 매우 많이 생성되게 유도한다
    kwargs["ISPARENT"] = True; kwargs["ELBOW_S"] = 3; kwargs["GAUSSIAN_SD"] = 1.0; kwargs["OUTLIER_STANDARD"] = "stricter"
else:       # Outlier를 매우 적게 생성되게 유도한다
    kwargs["ISPARENT"] = False;  kwargs["ELBOW_S"] = 1; kwargs["GAUSSIAN_SD"] = 1.5; kwargs["OUTLIER_STANDARD"] = "looser"
    
NUM_CLONE, membership, membership_p, membership_p_normalize, mixture, includeoutlier_inside = EMsoft.main(df, method, adjustment, **kwargs)   
membership_p_normalize_new = membership_p_normalize


#3-2. visualizationsingle
samplename_dict = {i:i for i in range(0,NUM_CLONE)}         # index_no :  clone_no = color_no   (clone 번호와 color 번호는 일치시킨다)
        ## Hard visualization
if kwargs["NUM_BLOCK"] >= 2:
    visualizationsingle.drawfigure_2d (membership, "NotSave", np_vaf, samplename_dict, includeoutlier_inside, "")
if kwargs["NUM_BLOCK"] == 1:
    visualizationsingle.drawfigure_1d (membership, "NotSave", np_vaf, samplename_dict, includeoutlier_inside)

        ## Soft visualization
if kwargs["NUM_BLOCK"] >= 2:
    visualizationsinglesoft.drawfigure_2d (membership, np.round(mixture, 3), membership_p_normalize,kwargs["MYEM_DIR"]+ "/3.soft.jpg" , np_vaf, samplename_dict, includeoutlier_inside, "")
if kwargs["NUM_BLOCK"] == 1:
    visualizationsinglesoft.drawfigure_1d (membership, np.round(mixture, 3), membership_p_normalize, kwargs["MYEM_DIR"]+ "/3.soft.jpg" , np_vaf, samplename_dict, includeoutlier_inside)


print ("\n현재 시각 : {}h:{}m:{}s    (걸린 시간 : {})".format(time.localtime().tm_hour, time.localtime().tm_min, round(time.localtime().tm_sec), datetime.datetime.now() - START_TIME ))








print ("\n\n=============== STEP #4.   OUTLIER HARD CLUSTERING  AND ADD ===============")
# 5. Outlier 뽑기 (Insdie, Outside구분)         # hard membership을 바탕으로 분류해도 충분하다
Outlier_method = kwargs["OUTLIER_METHOD"]
#Outlier_method = "Distance"  # "Distance"
OUTLIER_NO = NUM_CLONE - 1

if Outlier_method in ["ML", "ml"]:
    try :
        boundary_mixture = boundaryclone.main(mixture[:,:-1])    # 맨 마지막 outlier group은 빼고 돌려야지
        df_inside, df_inside_index, np_vaf_inside, df_outside, df_outside_index, np_vaf_outside = \
            extract.classifier_lightgbm (df, np_vaf, membership, boundary_mixture, OUTLIER_NO)
    except:
        df_inside, df_inside_index, np_vaf_inside, df_outside, df_outside_index, np_vaf_outside = \
            extract.classifier_distance (df, np_vaf, membership, mixture [:,:-1], OUTLIER_NO)       

elif Outlier_method in ["Distance", "distance"]:
    df_inside, df_inside_index, np_vaf_inside, df_outside, df_outside_index, np_vaf_outside = \
        extract.classifier_distance (df, np_vaf, membership, mixture [:,:-1], OUTLIER_NO)   
        # 맨 마지막 outlier group은 거리계산에서 뺴줘야지



if ISPARENT == False:   # 더 돌릴 것도 없이 다 outlier 처리한다
    membership_total = [i for i in membership]   # 얕은 사본을 만들어주고
    for i in df_outside_index:
        membership_total[i] = OUTLIER_NO
    mixture_total = mixture.copy()

    if includeoutlier_inside == "Yes":
        includeoutlier_total = "Yes"
        child_clone = set(membership) -  set([OUTLIER_NO])
    elif includeoutlier_inside == "No":
        includeoutlier_total = "No"
        child_clone = set(membership)

    if includeoutlier_total == "Yes":             # 맨 끝 column에 outlier mixture를 재정비해주기 (M step)
        for i in range (kwargs["NUM_BLOCK"]):
            sum_depth, sum_alt = 0, 0
            for k in range(kwargs["RANDOM_PICK"]):
                if membership_total[k] == np.max(membership_total):    # outlier 찾기
                    sum_depth = sum_depth + df[k][i]["depth"]
                    sum_alt = sum_alt + df[k][i]["alt"]
            mixture_total[i][-1] = round((sum_alt * 2) / sum_depth, 2)

if ISPARENT == True:
    ############### INSIDE OUTLIER ##########################################################
    membership_inside = [0] * len(df_inside)
    mixture_inside  = EMhard.Mstep (membership_inside, df_inside, 1, kwargs["NUM_BLOCK"], len(df_inside), "False", "Hard")    # inside outlier 의 mixture를 구해주기
    membership_inside = [OUTLIER_NO for i in membership_inside ]            # 번호를 뒤로 밀어 준다


    ############### OUTSIDE OUTIER :  df_outside 들끼리만 EM 돌리기   (hard clustering을 돌려도 될 듯) ############
    kwargs["NUM_CLONE_TRIAL_START"] = 1; kwargs["NUM_CLONE_TRIAL_END"] = 5;  kwargs["NUM_CLONE_FORCE"] = 1;  kwargs["MIN_CLUSTER_SIZE"] = int ( len(df_outside) / 10);  kwargs["ELBOW_S"] = 0.99; kwargs["GAUSSIAN_SD"] = 1.5; kwargs["VERBOSE"] = 0
    if len(df_outside) < kwargs["NUM_CLONE_TRIAL_END"]:
        kwargs["NUM_CLONE_TRIAL_END"] = len(df_outside)

    if len(df_outside) == 0:       # 아예 outlier가 하나도 발견 안 됐을 경우
        includeoutlier_outside  = "No"
        membership_outside = []

    else:
        method = "sil+outlier"; adjustment = "False"
        NUM_CLONE_outside, membership_outside, mixture_outside, includeoutlier_outside = EMhard.main(df_outside, method, adjustment, **kwargs)
        membership_outside = [i + OUTLIER_NO for i in membership_outside ]            # 번호를 뒤로 밀어 준다
        #includeoutlier_outside = "Yes"

        if includeoutlier_outside == "Yes":
            OUTLIER_OUTLIER_NO = np.max (membership_outside)
            membership_inside = [OUTLIER_OUTLIER_NO for i in membership_inside ]            # 번호를 맨 뒤로 한번 더 밀어 준다   (outlier 끼리는 색깔을 통일하려고)
            # print ("previous membership : {} + {} \nmembership_inside : {} -> {} \nmembership_outside : {} + {}"\
            # .format(set(membership) - set([OUTLIER_NO]), OUTLIER_NO, OUTLIER_NO, OUTLIER_OUTLIER_NO, set(membership_outside) - set([OUTLIER_OUTLIER_NO]), OUTLIER_OUTLIER_NO))
        else:
            OUTLIER_OUTLIER_NO = np.max (membership_outside) + 1
            membership_inside = [OUTLIER_OUTLIER_NO for i in membership_inside ]            # 번호를 맨 뒤로 한번 더 밀어 준다   
            # print ("previous membership : {} + {} \nmembership_inside : {} -> {} \nmembership_outside : {}"\
            # .format(set(membership) - set([OUTLIER_NO]), OUTLIER_NO, OUTLIER_NO, OUTLIER_OUTLIER_NO, set(membership_outside) ))


        #5-2. visualizationsingle
        #samplename_dict = {i:i  for i in range(OUTLIER_NO , OUTLIER_OUTLIER_NO + 1)}
        #visualizationsingle.drawfigure_2d (membership_outside, "NotSave", np_vaf_outside, samplename_dict, includeoutlier_outside, "")
        #visualizationsingle.drawfigure_1d (membership_outside, "NotSave", np_vaf_outside, samplename_dict, includeoutlier_outside)



    # 6. Total 한번에 그림그리기
    if (includeoutlier_inside == "Yes") | (includeoutlier_outside == "Yes"):
        includeoutlier_total = "Yes"
    else:
        includeoutlier_total = "No"

    membership_total = [i for i in membership]   # 얕은 사본을 만들어주고
    for i, index in enumerate(df_inside_index):
        membership_total[index] = membership_inside[i]    
    for i, index in enumerate(df_outside_index):
        membership_total[index] = membership_outside[i]


    # Mixture total 만들기
    if includeoutlier_inside == "Yes":
        child_clone = set(membership) -  set([OUTLIER_NO])
        if len(df_outside) > 0:     # 아예 outlier가 없으면
            mixture_total = np.concatenate (((mixture[: , : -1]), (mixture_outside)), axis = 1)
            if includeoutlier_outside == "No":     # 맨 마지막에 insideoutlier 를 넣어줘야 하니까
                    mixture_total = np.concatenate ( (mixture_total, np.zeros((kwargs["NUM_BLOCK"], 1), dtype = "float")) , axis = 1)
        else:    # outlier_outside가 하나도 없을 경우
            mixture_total = mixture
    else:
        child_clone = set(membership)

    if includeoutlier_total == "Yes":             # 맨 끝 column에 outlier mixture를 재정비해주기 (M step)
        parent_clone = set(membership_outside)  - set ( [np.max(membership_total)] )
        for i in range (kwargs["NUM_BLOCK"]):
            sum_depth, sum_alt = 0, 0
            for k in range(kwargs["RANDOM_PICK"]):
                if membership_total[k] == np.max(membership_total):    # outlier 찾기
                    sum_depth = sum_depth + df[k][i]["depth"]
                    sum_alt = sum_alt + df[k][i]["alt"]
            mixture_total[i][-1] = round((sum_alt * 2) / sum_depth, 2)
    else:
        parent_clone = set(membership_outside)
    
    print ("child clone : {}\ninside outlier : {}\nparent clone : {}\noutside outlier : {}\n\ntotal mixture : \n{}\n".format(child_clone, set(membership_inside), set(membership_outside)  - set ( [np.max(membership_total)] ) , set([np.max(membership_total)]),  mixture_total))



# 채점하고 색깔 다시 맞춰서 그리기
score_df, score = scoring.mixturebased(mixture_answer, mixture_total, membership_answer, membership_total, samplename_dict_input, samplename_dict_input_rev, includeoutlier_total, **kwargs)
visualizationpair.drawfigure_2d (membership_answer, mixture_answer, membership_total, mixture_total, score_df, kwargs["MYEM_DIR"] + "/4.total.jpg", "ANSWER", "MY ALGORITHM", np_vaf, includeoutlier_total,  dimensionreduction="None")


print ("\n\n=============== STEP #5.   ONE MORE EM  ===============")

mixture_onemore_soft,  mixture_onemore_hard, membership_onemore, membership_p_onemore = \
    onemore.main(NUM_MUTATION, df, mixture_total, child_clone,  **kwargs)

# 채점하고 색깔 다시 맞춰서 그리기
score_df_onemore, score_onemore = scoring.mixturebased(mixture_answer, mixture_onemore_hard, membership_answer, membership_onemore, samplename_dict_input, samplename_dict_input_rev, includeoutlier_total, **kwargs)
visualizationpair.drawfigure_2d (membership_answer, mixture_answer, membership_onemore, mixture_onemore_hard, score_df_onemore, kwargs["MYEM_DIR"]  + "/5.onemore_hard.jpg" , "ANSWER : MRS", "MY ALGORITHM", np_vaf, includeoutlier_total,  dimensionreduction="None")
visualizationpair.drawfigure_2d (membership_answer, mixture_answer, membership_onemore, mixture_onemore_soft, score_df_onemore, kwargs["MYEM_DIR"]  + "/5.onemore_soft.jpg" , "ANSWER : MRS", "MY ALGORITHM", np_vaf, includeoutlier_total,  dimensionreduction="None")
subprocess.run (["cp " + kwargs["MYEM_DIR"] + "/5.onemore_hard.jpg  " + kwargs["COMBINED_OUTPUT_DIR"]  + "/TGIL_hard.jpg"], shell = True)

score, score_df = score_onemore, score_df_onemore
mixture_total, membership_total = mixture_onemore_hard, membership_onemore




print ("\n\n=============== STEP #6.   RESULTS & EVALUATION ===============")


Y_index = result.Yindex(score_df, 5)
print ("\n[RESULTS]\n\nMyEM\t{}/{}\nNUM_CLONE\t{}\nY-index\t{}\n{}\n".format(score,  kwargs["RANDOM_PICK"], mixture_total.shape[1], Y_index, list(mixture_total)))
if (kwargs["FP_RATIO"] != 0) | (kwargs["FP_2D"] == "True"):
    print ("[FP ANALYSIS]")
    answeronly, intersection, myonly, sensitivity, PPV, F1 = result.FPmatrix(score_df)
    print ("answer FP {}개 중에 {}개 일치함".format( answeronly + intersection ,  intersection ) )
    print ("\tanswerFP only : {}\n\tintersection : {}\n\tmyFP only : {}\n".format( answeronly, intersection, myonly ))
else:
    answeronly, intersection, myonly, sensitivity, PPV, F1 = 0, 0, 0, None, None, None

with open (kwargs["MYEM_DIR"] + "/results.txt", "w", encoding = "utf8") as output_myEM:
    print ("NUM_CLONE\t{}\nscore\t{}/{}\nY-index\t{}\nansweronly\t{}\nintersection\t{}\nthistoolonly\t{}\nsensitivity\t{}\nPPV\t{}\nF1\t{}".
        format(mixture_total.shape[1], score, kwargs["RANDOM_PICK"], Y_index, answeronly, intersection, myonly, sensitivity, PPV, F1 ), file = output_myEM)
subprocess.run (["cp " + kwargs["MYEM_DIR"] + "/results.txt  " + kwargs["COMBINED_OUTPUT_DIR"]  + "/TGIL.results.txt"], shell = True)


pd.DataFrame(membership_total).to_csv (kwargs["MYEM_DIR"] + "/membership.txt", index = False, header= False,  sep = "\t" )
pd.DataFrame(membership_total).to_csv (kwargs["COMBINED_OUTPUT_DIR"] + "/TGIL.membership.txt", index = False, header= False,  sep = "\t" )
pd.DataFrame(mixture_total).to_csv (kwargs["MYEM_DIR"] + "/mixture.txt", index = False, header= False,  sep = "\t" )
pd.DataFrame(mixture_total).to_csv (kwargs["COMBINED_OUTPUT_DIR"] + "/TGIL.mixture.txt", index = False, header= False,  sep = "\t" )
pd.DataFrame(score_df).to_csv (kwargs["MYEM_DIR"] + "/scoredf.txt", index = False, header= True,  sep = "\t")
pd.DataFrame(score_df).to_csv (kwargs["COMBINED_OUTPUT_DIR"] + "/TGIL.scoredf.txt", index = False, header= True,  sep = "\t")

if ISPARENT == True:
    kwargs["PHYLOGENY_DIR"] = kwargs["MYEM_DIR"] + "/TGIL.phylogeny.txt"
    g = phylogeny.main(membership, membership_outside, mixture_total, includeoutlier_outside, **kwargs)


print ("\n현재 시각 : {}h:{}m:{}s     (걸린 시간 : {})".format(time.localtime().tm_hour, time.localtime().tm_min, round(time.localtime().tm_sec), datetime.datetime.now() - START_TIME ))





# # print ("\n\n=============== STEP #6.   PYCLONE RUNNING ===============")

# # COMMAND = " ".join( ["\"PyClone", "run_analysis_pipeline", "--in_files"] + [kwargs["PYCLONE_DIR"] + "/block{0}.tsv".format(i) for i in range(kwargs["NUM_BLOCK"])] + \
# # ["--samples \"block0\" \"block1\" --working_dir ", kwargs["PYCLONE_DIR"], " --tumour_contents 1.0 1.0  --num_iters 1000  --min_cluster_size 10 --max_clusters 7\""] )
# # INPUT_COMMAND = COMMAND.replace(" ", "&")
# # subprocess.run (["bash pyclone_pipe.sh " + INPUT_COMMAND], shell = True)


# # INPUT_NPVAF =  "/data/project/Alzheimer/YSscript/EM_MRS/data/npvaf/{0}.npvaf".format( ".".join( kwargs["INPUT_TSV"].split("/")[-1].split(".")[:-1]) )
# # INPUT_FILENAME = kwargs["PYCLONE_DIR"] + "/tables/loci.tsv"
# # OUTPUT_FILENAME = kwargs["PYCLONE_DIR"]  + "/4.total.jpg" 

# # score_df_pyclone, score_pyclone, membership_pyclone, mixture_pyclone = \
# #     pyclonesim.main(INPUT_NPVAF, INPUT_FILENAME, OUTPUT_FILENAME, mixture_answer, membership_answer, samplename_dict_input, samplename_dict_input_rev, **kwargs)


# # visualizationpair.drawfigure_2d (membership_answer, mixture_answer, membership_pyclone, mixture_pyclone, score_df_pyclone, OUTPUT_FILENAME,  "ANSWER : MRS", "PyClone", np_vaf, "No",  dimensionreduction="None")

# # Y_index_pyclone = result.Yindex(score_df_pyclone, 5)
# # print ("PyClone : {}개 / {}개\nY-index : {}\n\n".format( score_pyclone, kwargs["RANDOM_PICK"], Y_index_pyclone))

# #if kwargs["FP_RATIO"] != 0:
# # print ("[FP ANALYSIS]")
# # answeronly_pyclone, intersection_pyclone, myonly_pyclone, sensitivity_pyclone, PPV_pyclone, F1_pyclone = result.FPmatrix(score_df_pyclone)
# # print ("answer FP {}개 중에 {}개 일치함".format( answeronly_pyclone + intersection_pyclone ,  intersection_pyclone ) )
# # print ("\tanswerFP only : {}\n\tintersection : {}\n\tpycloneFP only : {}".format( answeronly_pyclone, intersection_pyclone , myonly_pyclone ))
# # else:
# #     answeronly_pyclone, intersection_pyclone, pyclone_only, sensitivity_pyclone, PPV_pyclone, F1_pyclone = 0, 0, 0, None, None, None


# # with open (kwargs["PYCLONE_DIR"]  + "/results.txt", "w", encoding = "utf8") as output_pyclone:
# #     print ("{}\n{}\n{}/{}\n{}\n{}\n{}\n{}\n{}\n{}\n{}".format(mixture_pyclone, mixture_pyclone.shape[1], score_pyclone, kwargs["RANDOM_PICK"], Y_index_pyclone, answeronly_pyclone, intersection_pyclone, myonly_pyclone, sensitivity_pyclone, PPV_pyclone, F1_pyclone ), file = output_pyclone)
# #     print (score_df_pyclone, file = output_pyclone)

# # print ("\n현재 시각 : {}h:{}m:{}s     (걸린 시간 : {})".format(time.localtime().tm_hour, time.localtime().tm_min, round(time.localtime().tm_sec), datetime.datetime.now() - START_TIME ))



print ("\n\n=============== STEP #7.   PYCLONEVI RUNNING ===============")

PYCLONEVI_START_TIME = datetime.datetime.now()

INPUT_TSV=kwargs["PYCLONEVI_DIR"] + "/input.tsv"
OUTPUT_H5=kwargs["PYCLONEVI_DIR"]  + "/output.h5"
OUTPUT_TSV=kwargs["PYCLONEVI_DIR"]  + "/output.tsv"

subprocess.run (["bash pyclonevi_pipe.sh " + INPUT_TSV + " " + OUTPUT_H5 + " " + OUTPUT_TSV], shell = True)

INPUT_PYCLONEVI_RESULT = OUTPUT_TSV
INPUT_NPVAF = kwargs["NPVAF_DIR"] + ".npvaf"
OUTPUT_FILENAME = kwargs["PYCLONEVI_DIR"]  + "/4.total.jpg"  

score_df_pyclonevi, score_pyclonevi, membership_pyclonevi, mixture_pyclonevi = \
    pyclonevisim.main(INPUT_PYCLONEVI_RESULT,  INPUT_NPVAF, OUTPUT_FILENAME, mixture_answer, membership_answer, samplename_dict_input, samplename_dict_input_rev, **kwargs)

visualizationpair.drawfigure_2d (membership_answer, mixture_answer, membership_pyclonevi, mixture_pyclonevi, score_df_pyclonevi, OUTPUT_FILENAME,  "ANSWER", "pyclone_vi", np_vaf, "No",  dimensionreduction="None")
subprocess.run (["cp " + OUTPUT_FILENAME + "  " + kwargs["COMBINED_OUTPUT_DIR"]  + "/pyclonevi.jpg"], shell = True)

Y_index_pyclonevi = result.Yindex(score_df_pyclonevi)
print ("\n[RESULTS]\n\npyclone_vi\t{}/{}\nNUM_CLONE\t{}\nY-index\t{}\n{}\n".format(score_pyclonevi, kwargs["RANDOM_PICK"], mixture_pyclonevi.shape[1],  Y_index_pyclonevi, list(mixture_pyclonevi)))

if (kwargs["FP_RATIO"] != 0) | (kwargs["FP_2D"] == "True"):
    print ("[FP ANALYSIS]")
    answeronly_pyclonevi, intersection_pyclonevi, pyclonevi_only, sensitivity_pyclonevi, PPV_pyclonevi, F1_pyclonevi = result.FPmatrix(score_df_pyclonevi)
    print ("answer FP {}개 중에 {}개 일치함".format( answeronly_pyclonevi + intersection_pyclonevi ,  intersection_pyclonevi ) )
    print ("\tanswerFP only : {}\n\tintersection : {}\n\tpycloneviFP only : {}".format( answeronly_pyclonevi, intersection_pyclonevi, pyclonevi_only ))
else:
    answeronly_pyclonevi, intersection_pyclonevi, pyclonevi_only, sensitivity_pyclonevi, PPV_pyclonevi, F1_pyclonevi = 0, 0, 0, None, None, None


with open (kwargs["PYCLONEVI_DIR"]  + "/results.txt", "w", encoding = "utf8") as output_pyclonevi:
    print ("NUM_CLONE\t{}\nscore\t{}/{}\nY-index\t{}\nansweronly\t{}\nintersection\t{}\nthistoolonly\t{}\nsensitivity\t{}\nPPV\t{}\nF1\t{}".
            format(mixture_pyclonevi.shape[1], score_pyclonevi, kwargs["RANDOM_PICK"], Y_index_pyclonevi, answeronly_pyclonevi, intersection_pyclonevi, pyclonevi_only, sensitivity_pyclonevi, PPV_pyclonevi, F1_pyclonevi ), file = output_pyclonevi)
subprocess.run (["cp " + kwargs["PYCLONEVI_DIR"] + "/results.txt  " + kwargs["COMBINED_OUTPUT_DIR"]  + "/pyclonevi.results.txt"], shell = True)


pd.DataFrame(membership_pyclonevi).to_csv (kwargs["PYCLONEVI_DIR"] + "/membership.txt", index = False, header= False,  sep = "\t" )
pd.DataFrame(membership_pyclonevi).to_csv (kwargs["COMBINED_OUTPUT_DIR"]  + "/pyclonevi.membership.txt", index = False, header= False,  sep = "\t" )
pd.DataFrame(mixture_pyclonevi).to_csv (kwargs["PYCLONEVI_DIR"] + "/results.mixture.txt", index = False, header= False,  sep = "\t" )
pd.DataFrame(mixture_pyclonevi).to_csv (kwargs["COMBINED_OUTPUT_DIR"]  + "/pyclonevi.mixture.txt", index = False, header= False,  sep = "\t" )
pd.DataFrame(score_df_pyclonevi).to_csv (kwargs["PYCLONEVI_DIR"] + "/results.scoredf.txt", index = False, header= True,  sep = "\t")
pd.DataFrame(score_df_pyclonevi).to_csv (kwargs["COMBINED_OUTPUT_DIR"]  + "/pyclonevi.scoredf.txt", index = False, header= True,  sep = "\t" )


print ("\n현재 시각 : {}h:{}m:{}s     (걸린 시간 : {})".format(time.localtime().tm_hour, time.localtime().tm_min, round(time.localtime().tm_sec), datetime.datetime.now() - PYCLONEVI_START_TIME ))





print ("\n\n=============== STEP #8.   SCICLONE RUNNING ===============")

SCICLONE_START_TIME = datetime.datetime.now()
print ("\nNOW SCICLONE IS STARTED  :  {}h:{}m:{}s\n\n".format( time.localtime().tm_hour, time.localtime().tm_min, round(time.localtime().tm_sec) ))

INPUT_First=kwargs["SCICLONE_DIR"] + "/block0.dat"
INPUT_Second=kwargs["SCICLONE_DIR"] + "/block1.dat"

#os.system("bash sciclone_pipe.sh " + INPUT_First + " " + INPUT_Second + " " + kwargs["SCICLONE_DIR"])
subprocess.run (["bash sciclone_pipe.sh " + INPUT_First + " " + INPUT_Second + " " + kwargs["SCICLONE_DIR"]], shell = True)


INPUT_SCICLONE_RESULT = kwargs["SCICLONE_DIR"]+ "/results.tsv"
INPUT_NPVAF = kwargs["NPVAF_DIR"] + ".npvaf"
OUTPUT_FILENAME = kwargs["SCICLONE_DIR"]  + "/4.total.jpg"  

score_df_sciclone, score_sciclone, membership_sciclone, mixture_sciclone = \
    sciclonesim.main(INPUT_SCICLONE_RESULT,  INPUT_NPVAF, OUTPUT_FILENAME, mixture_answer, membership_answer, samplename_dict_input, samplename_dict_input_rev, **kwargs)

visualizationpair.drawfigure_2d (membership_answer, mixture_answer, membership_sciclone, mixture_sciclone, score_df_sciclone, OUTPUT_FILENAME,  "ANSWER", "SciClone", np_vaf, "No",  dimensionreduction="None")
subprocess.run (["cp " + OUTPUT_FILENAME + "  " + kwargs["COMBINED_OUTPUT_DIR"]  + "/sciclone.jpg"], shell = True)

Y_index_sciclone = result.Yindex(score_df_sciclone)
print ("\n[RESULTS]\n\nSciClone_vi\t{}/{}\nNUM_CLONE\t{}\nY-index\t{}\n{}\n".format(score_sciclone, kwargs["RANDOM_PICK"], mixture_sciclone.shape[1],  Y_index_sciclone, list(mixture_sciclone)))

if (kwargs["FP_RATIO"] != 0) | (kwargs["FP_2D"] == "True"):
    print ("[FP ANALYSIS]")
    answeronly_sciclone, intersection_sciclone, sciclone_only, sensitivity_sciclone, PPV_sciclone, F1_sciclone = result.FPmatrix(score_df_sciclone)
    print ("answer FP {}개 중에 {}개 일치함".format( answeronly_sciclone + intersection_sciclone ,  intersection_sciclone ) )
    print ("\tanswerFP only : {}\n\tintersection : {}\n\tscicloneFP only : {}".format( answeronly_sciclone, intersection_sciclone, sciclone_only ))
else:
    answeronly_sciclone, intersection_sciclone, sciclone_only, sensitivity_sciclone, PPV_sciclone, F1_sciclone = 0, 0, 0, None, None, None


with open (kwargs["SCICLONE_DIR"]  + "/results.txt", "w", encoding = "utf8") as output_sciclone:
    print ("NUM_CLONE\t{}\nscore\t{}/{}\nY-index\t{}\nansweronly\t{}\nintersection\t{}\nthistoolonly\t{}\nsensitivity\t{}\nPPV\t{}\nF1\t{}".
        format(mixture_sciclone.shape[1], score_sciclone, kwargs["RANDOM_PICK"], Y_index_sciclone, answeronly_sciclone, intersection_sciclone, sciclone_only, sensitivity_sciclone, PPV_sciclone, F1_sciclone ), file = output_sciclone)
subprocess.run (["cp " + kwargs["SCICLONE_DIR"] + "/results.txt  " + kwargs["COMBINED_OUTPUT_DIR"]  + "/sciclone.results.txt"], shell = True)

pd.DataFrame(membership_sciclone).to_csv (kwargs["SCICLONE_DIR"] + "/membership.txt", index = False, header= False,  sep = "\t" )
pd.DataFrame(membership_sciclone).to_csv (kwargs["COMBINED_OUTPUT_DIR"]  + "/sciclone.membership.txt", index = False, header= False,  sep = "\t" )
pd.DataFrame(mixture_sciclone).to_csv (kwargs["SCICLONE_DIR"] + "/results.mixture.txt", index = False, header= False,  sep = "\t" )
pd.DataFrame(mixture_sciclone).to_csv (kwargs["COMBINED_OUTPUT_DIR"]  + "/sciclone.mixture.txt", index = False, header= False,  sep = "\t" )
pd.DataFrame(score_df_sciclone).to_csv (kwargs["SCICLONE_DIR"] + "/results.scoredf.txt", index = False, header= True,  sep = "\t" )
pd.DataFrame(score_df_sciclone).to_csv (kwargs["COMBINED_OUTPUT_DIR"]  + "/sciclone.scoredf.txt", index = False, header= True,  sep = "\t" )


print ("\n현재 시각 : {}h:{}m:{}s     (걸린 시간 : {})".format(time.localtime().tm_hour, time.localtime().tm_min, round(time.localtime().tm_sec), datetime.datetime.now() - SCICLONE_START_TIME ))
