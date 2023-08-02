import datapreparationold, datapreparation220712, comb, extract, scoring, boundaryclone, graph, phylogeny, scoring, os
import fppick
import pyclonesim, sciclonesim
import EMhard, EMsoft
import visualizationsingle, visualizationpair, visualizationsinglesoft, visualizationpairsoft
import numpy as np
from kneed import KneeLocator


#1. MRS :  Data preparationold  (1000개 random으로 고르고 + Pyclone dataset 준비)
# kwargs = {"NUM_BLOCK_INPUT": 3, "NUM_BLOCK": 2, "NUM_CLONE_TRIAL_START" : 3, "NUM_CLONE_TRIAL_END" : 5, "NUM_CLONE_FORCE" : 4,
# "RANDOM_PICK":2000, "FP_RATIO":0.3, "TRIAL_NO" : 10, "VERBOSE" : 1 }
# inputdf, df, np_vaf, membership_answer, mixture_answer, mutation_id, samplename_dict = datapreparationold.main(**kwargs)

# visualizationsingle.drawfigure_2d (membership_answer, "./output/0.MRS_old_answer.jpg", np_vaf, samplename_dict, "")





#2. MRS  ; Data preparation_new

# M1-5_M1-8 : S0, V1, V2 이고 가장 쉬운 문제 (pyclone n = 4 : 89.6점,   n = 10 : 40점)

kwargs = {"INPUT_TSV" : "/data/project/Alzheimer/EM_cluster/EM_input/MRS_2_sample/M1-5_M1-8_input.txt",  "NUM_BLOCK_INPUT": 2, "NUM_BLOCK": 2, "NUM_CLONE_TRIAL_START" : 2, "NUM_CLONE_TRIAL_END" : 6, "NUM_CLONE_FORCE" : 4, 
                "RANDOM_PICK":300, "AXIS_RATIO":0.1, "PARENT_RATIO": 0.2, "FP_RATIO":0.2, "TRIAL_NO" : 5, "DEPTH_CUTOFF" : 100, "VERBOSE" : 1, "ELBOW_S" : 3 , "GAUSSIAN_SD" : 1, "MIN_CLUSTER_SIZE" : 5, "OUTLIER_STANDARD" : "looser"}
SAMPLENAME = kwargs["INPUT_TSV"].split("/")[-1].split(".")[0]     # 'M1-5_M1-8_input'
kwargs["SCICLONE_DIR"] = "/data/project/Alzheimer/YSscript/EM_MRS/data/sciclone/" + SAMPLENAME 
kwargs["PYCLONE_DIR"] = "/data/project/Alzheimer/YSscript/EM_MRS/data/pyclone/" + SAMPLENAME 
kwargs["PYCLONEVI_DIR"] = "/data/project/Alzheimer/YSscript/EM_MRS/data/pyclone-vi/" + SAMPLENAME 


inputdf, df, np_vaf, membership_answer, mixture_answer,  mutation_id, samplename_dict  = datapreparation220712.main(**kwargs)

visualizationsingle.drawfigure_2d (membership_answer, "./output/0.MRS_new_answer.jpg", np_vaf, samplename_dict, "")





#3. Moore :  Data preparation_new 
# "/data/project/Alzheimer/EM_cluster/EM_input/Moore_2_sample/PD43851/colon_crypt/F5_G6_input.txt"       # 2차원 monolonal
# "/data/project/Alzheimer/EM_cluster/Moore_data/Donor/PD43851/stomach_gland/PD43851x_P52_STM_B10.txt"    # 1차원 monoclonal인데 biclonal로 계속 오류남.   soft clustering으로 바꾸고 나서 monoclonal로 줌
# "/data/project/Alzheimer/EM_cluster/Moore_data/Donor/PD43851/oesophagus_epithelium/PD43851k_P52_OSPHG_E12.txt"  # 1차원 binclonal
# "/data/project/Alzheimer/EM_cluster/Moore_data/Donor/PD43851/oesophagus_epithelium/PD43851k_P53_OSPHG_B2.txt"  # 1차원  biclonal이 맞는 듯.  soft clustering에서도 biclonal로 줌
# "/data/project/Alzheimer/EM_cluster/EM_input/Moore_2_sample/PD28690/visceral_fat/L3_L4_input.txt"      # 2차원 monoclonal.  soft clustering으로 바꾸고 해결함
#  "/data/project/Alzheimer/EM_cluster/Moore_data/Donor/PD28690/oesophagus_epithelium/PD28690bl_OES2_CU2.txt"      # 1차원 polyclonal.  soft clustering + max로 해야 답을 4로 줌

# kwargs = {"INPUT_TSV" : "/data/project/Alzheimer/EM_cluster/Moore_data/Donor/PD43851/oesophagus_epithelium/PD43851k_P53_OSPHG_C2.txt",  "NUM_BLOCK_INPUT": 1, "NUM_BLOCK": 1, "NUM_CLONE_TRIAL_START" : 1, "NUM_CLONE_TRIAL_END" : 3, "NUM_CLONE_FORCE" : 1,
#                 "RANDOM_PICK":1000, "AXIS_RATIO":0, "PARENT_RATIO": 0, "FP_RATIO":0, "TRIAL_NO" : 5, "DEPTH_CUTOFF" : 5, "VERBOSE" : 0, "ELBOW_S" : 3 , "GAUSSIAN_SD" : 1, "MIN_CLUSTER_SIZE" : 5}
# OUTPUT_DIR = "/data/project/Alzheimer/EM_cluster/Moore_data/Donor/YSChung/Moore_1_sample/" + kwargs["INPUT_TSV"].split("/")[-2] + "/" + kwargs["INPUT_TSV"].split("/")[-1].replace(".txt", "")
# kwargs["SCICLONE_TSV"], kwargs["PYCLONE_TSV"], kwargs["PYCLONEVI_TSV"] = OUTPUT_DIR, OUTPUT_DIR, OUTPUT_DIR

# inputdf, df, np_vaf, membership_answer, mixture_answer,  mutation_id, samplename_dict  = datapreparation220712.main(**kwargs)
# membership_answer = [0] * kwargs["RANDOM_PICK"]
# samplename_dict = {0 : 0}

# if kwargs["NUM_BLOCK"] >= 2:         # 2차원이라면 정답 그림 그려준다
#     visualizationsingle.drawfigure_2d (membership_answer, "./output/0." + str(kwargs["INPUT_TSV"].split("/")[-2:]) + ".jpg", np_vaf, samplename_dict, "")



import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import kde
import palettable
tabl = palettable.tableau.Tableau_20.mpl_colors
Gr_10 = palettable.scientific.sequential.GrayC_20.mpl_colors
colorlist = [i for i in tabl]
plt.figure (figsize = (3,3))
try:
    plt.suptitle ("{0}".format(' : '. join(kwargs["INPUT_TSV"].split("/")[-2:])))
except:
    plt.suptitle ("MRS dataset")
plt.xlabel ("VAF")
plt.ylabel ("Density")
max_y = 0

for i in range(kwargs["NUM_BLOCK"]):
    np_vaf[:,i]
    x = np.linspace(0,1,101)
    kde_np_vaf_new = kde.gaussian_kde(np_vaf[:,i])
    y = kde_np_vaf_new(x)
    if max_y < np.max(y):
        max_y = np.max(y)
    plt.plot(x, y, color = colorlist[i], label = "sample {0}".format(i))

plt.axis ([0,  1.01,  0,  max_y * 1.3])
plt.legend()

try:
    plt.savefig("./output/1." + kwargs["INPUT_TSV"].split("/")[-2:] + ".jpg")
except:
    plt.savefig("./output/1.MRS dataset.jpg")





#2-1. EM 돌리기  (원하는 df만 넣어주면 가능)     (elbow,max,force)+(normal, outlier)    # includeoutlier == Yes라면 +1 해서 출력된다
kwargs["ELBOW_S"] = 1; kwargs["GAUSSIAN_SD"] = 1; kwargs["TRIAL_NO"] = 10; kwargs["NUM_CLONE_TRIAL_START"] = 3;  kwargs["NUM_CLONE_TRIAL_END"] = 6; kwargs["NUM_CLONE_FORCE"] = 1; kwargs["MIN_CLUSTER_SIZE"] = 15; kwargs["VERBOSE"] =   0; kwargs["OUTLIER_STANDARD"] = "looser"
method = "best+outlier"; adjustment = "True"
NUM_CLONE, membership, membership_p, membership_p_normalize, mixture, includeoutlier = EMsoft.main(df, method, adjustment, **kwargs)   

#2-2. visualizationsingle
samplename_dict = {i:i for i in range(0,NUM_CLONE)}         # index_no :  clone_no = color_no   (clone 번호와 color 번호는 일치시킨다)
## Hard visualization
if kwargs["NUM_BLOCK"] >= 2:
    visualizationsingle.drawfigure_2d (membership, "./output/2.EMsoft_h_2D_" + method + ".jpg", np_vaf, samplename_dict, includeoutlier, "")
if kwargs["NUM_BLOCK"] == 1:
    visualizationsingle.drawfigure_1d (membership, "./output/2.EMsoft_h_1D_" + method + ".jpg", np_vaf, samplename_dict, includeoutlier)


## Soft visualization
if kwargs["NUM_BLOCK"] >= 2:
    visualizationsinglesoft.drawfigure_2d (membership, mixture, membership_p_normalize,"./output/2.EMsoft_s_2D_" + method + ".jpg", np_vaf, samplename_dict, includeoutlier, "")
if kwargs["NUM_BLOCK"] >= 1:
    visualizationsinglesoft.drawfigure_1d (membership, mixture, membership_p_normalize,"./output/2.EMsoft_s_1D_" + method + ".jpg", np_vaf, samplename_dict, includeoutlier)

# if "outlier" in method:
#     print ("NUM_CLONE : {0} (outlier clone 포함)\nmixture : {1}\n각 clone 당 개수 : {2}\n".format(NUM_CLONE, mixture, np.unique(membership, return_counts = True)))
# else:
#     print ("NUM_CLONE : {0} \nmixture : {1}\n각 clone 당 개수 : {2}\n".format(NUM_CLONE, mixture, np.unique(membership, return_counts = True)))




mixture, membership, membership_p, membership_p_normalize_new, axis_index, NUM_CLONE = fppick.main(mixture, membership, membership_p, df, **kwargs)


# 4. Outlier 뽑기 (Insdie, Outside구분)         # hard membership을 바탕으로 분류해도 충분하다
Outlier_method = "Distance"  # "Distance"
OUTLIER_NO = NUM_CLONE - 1

if Outlier_method in ["ML", "ml"]:
    boundary_mixture = boundaryclone.main(mixture[:,:-1])    # 맨 마지막 outlier group은 빼고 돌려야지
    df_inside, df_inside_index, np_vaf_inside, df_outside, df_outside_index, np_vaf_outside = \
        extract.classifier_lightgbm (df, np_vaf, membership, boundary_mixture, OUTLIER_NO)

elif Outlier_method in ["Distance", "distance"]:
    df_inside, df_inside_index, np_vaf_inside, df_outside, df_outside_index, np_vaf_outside = \
        extract.classifier_distance (df, np_vaf, membership, mixture [:,:-1], OUTLIER_NO)   
        # 맨 마지막 outlier group은 거리계산에서 뺴줘야지


############### INSIDE OUTLIER ##########################################################
membership_inside = [0] * len(df_inside)
mixture_inside  = EMhard.Mstep (membership_inside, df_inside, 1, kwargs["NUM_BLOCK"], len(df_inside), "False", "")    # inside outlier 의 mixture를 구해주기
membership_inside = [OUTLIER_NO for i in membership_inside ]            # 번호를 뒤로 밀어 준다


############### OUTSIDE OUTIER :  df_outside 들끼리만 EM 돌리기   (hard clustering을 돌려도 될 듯) ############

kwargs["NUM_CLONE_TRIAL_START"] = 1; kwargs["NUM_CLONE_TRIAL_END"] = 4;  kwargs["NUM_CLONE_FORCE"] = 1;  kwargs["ELBOW_S"] = 0.99; kwargs["GAUSSIAN_SD"] = 1.5; kwargs["VERBOSE"] = 0

method = "best+outlier"; adjustment = "False"
NUM_CLONE_outside, membership_outside, mixture_outside, includeoutlier = EMhard.main(df_outside, method, adjustment, **kwargs)
membership_outside = [i + OUTLIER_NO for i in membership_outside ]            # 번호를 뒤로 밀어 준다
if includeoutlier == "Yes":
    OUTLIER_OUTLIER_NO = np.max (membership_outside)
    membership_inside = [OUTLIER_OUTLIER_NO for i in membership_inside ]            # 번호를 맨 뒤로 한번 더 밀어 준다   (outlier 끼리는 색깔을 통일하려고)
    print ("previous membership : {} + {} \nmembership_inside : {} -> {} \n:membership_outside : {} + {}"\
    .format(set(membership) - set([OUTLIER_NO]), OUTLIER_NO, OUTLIER_NO, OUTLIER_OUTLIER_NO, set(membership_outside) - set([OUTLIER_OUTLIER_NO]), OUTLIER_OUTLIER_NO))
else:
    OUTLIER_OUTLIER_NO = np.max (membership_outside) + 1
    membership_inside = [OUTLIER_OUTLIER_NO for i in membership_inside ]            # 번호를 맨 뒤로 한번 더 밀어 준다   
    print ("previous membership : {} + {} \nmembership_inside : {} -> {} \n:membership_outside : {}"\
    .format(set(membership) - set([OUTLIER_NO]), OUTLIER_NO, OUTLIER_NO, OUTLIER_OUTLIER_NO, set(membership_outside) ))


#5-2. visualizationsingle
samplename_dict = {i:i  for i in range(OUTLIER_NO , OUTLIER_OUTLIER_NO + 1)}
visualizationsingle.drawfigure_2d (membership_outside, "./output/5.EM_outlieronly_2D_" + method + ".jpg", np_vaf_outside, samplename_dict, includeoutlier, "")
#visualizationsingle.drawfigure_1d (membership_outside, "./output/5.EM_outlieronly_1D_" + method + ".jpg", np_vaf_outside, samplename_dict, includeoutlier)


# 6. Total 한번에 그림그리기
membership_total = [i for i in membership]   # 얕은 사본을 만들어주고
for i, index in enumerate(df_inside_index):
    membership_total[index] = membership_inside[i]    
for i, index in enumerate(df_outside_index):
    membership_total[index] = membership_outside[i]

mixture_total = np.concatenate (((mixture[: , : -1]), (mixture_outside)), axis = 1)


print ("previous membership : {0}\ninside membership : {1}\noutside membership : {2}\n\ntotal mixture : \n{3}\n".format(set(membership), set(membership_inside), set(membership_outside), mixture_total))
print (np.unique(membership_total, return_counts = True))

samplename_dict = {i:i for i in range(np.max(membership_total) + 1)}

# 채점하고 색깔 다시 맞춰서 그리기
membership_answer_max, score, sample_dict_rev = scoring.main(membership_answer, list(membership_total), **kwargs)
print ("score : {0}점 / 100점".format(score))

visualizationpair.drawfigure_2d (membership_answer_max, sample_dict_rev, membership_total, mixture_total ,"./output/6.EM_total_" + method + "_paired.jpg", np_vaf, "")
visualizationpairsoft.drawfigure_2d (membership_answer_max, sample_dict_rev, membership_total, membership_outside, membership_p_normalize, membership_p_normalize_new, mixture_total , 
                        axis_index, df_inside_index, df_outside_index, "./output/6.EM_total_" + method + "_paired.jpg", np_vaf, includeoutlier, "")



from sklearn.metrics import confusion_matrix
from sklearn.metrics import precision_score, recall_score, f1_score

for key, value in sample_dict_rev.items():
    if value == "FP":
        my_fp_index =  list ( np.array (np.where( np.array(membership_total) == key ) ) [0] )
        answer_fp_index = list ( np.array (np.where ( np.array(membership_answer) == "FP")) [0] )

        print ("answer FP {}개 중에 {}개 일치함".format( len(set(answer_fp_index)) , len( set(my_fp_index) & set(answer_fp_index) ) ) )
        print ("\tanswerFP only : {}\n\tintersection : {}\n\tmyFP only : {}".format( len( set(answer_fp_index) - set(my_fp_index) ) , len( set(my_fp_index) & set(answer_fp_index) ) ,  len( set(my_fp_index) - set(answer_fp_index) )))

        sensitivity =  round( len( set(my_fp_index) & set(answer_fp_index) ) / len(set(answer_fp_index)) , 2)
        PPV = round( len( set(my_fp_index) & set(answer_fp_index) ) / len(set(my_fp_index)) , 2)
        F1 = round( 2 * (sensitivity*PPV) / (sensitivity + PPV), 2)

        print ("\t\tsensitivity : {}\tPPV : {}\tF1 : {}".format(sensitivity, PPV, F1))

        break



# 7. phylogeny 관계를 유추한다
g = phylogeny.main(membership, membership_outside, mixture_total, includeoutlier, **kwargs)