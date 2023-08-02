import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import palettable
import scoring
import copy

def main (INPUT_QUANTUMCLONE_RESULT, INPUT_NPVAF, OUTPUT_FILENAME,  mixture_answer, membership_answer, membership_answer_numerical, samplename_dict_input, samplename_dict_input_rev, **kwargs):   
    #  INPUT_QUANTUMCLONE_RESULT =  "/data/project/Alzheimer/EM_cluster/old/pilot/04.EM_input/quantumclone/result/results.tsv"
    #  OUTPUT_FILENAME = "./output/MRS_quantumclone.jpg" 

    
    df = pd.read_csv (INPUT_QUANTUMCLONE_RESULT + "/output0.txt", sep = "\t")
    df = df.sort_values (by = "id", axis = 0, ascending = True)
    df["id"] = df["id"] - 1
    df = df.reset_index (drop = True)

    
    # 각 mutation_id의 vaf 정보를 가져오기 위함
    np_vaf = pd.read_csv(INPUT_NPVAF, sep = "\t")
    np_vaf.rename(columns = {"Unnamed: 0":"id"}, inplace = True) 
    for k in range (kwargs["NUM_MUTATION"]):
        #np_vaf["id"].iloc[k] = np_vaf["id"].iloc[k].split("_")[1]
        np_vaf["id"].iloc[k] = k
    np_vaf = np_vaf.astype ( { 'id':'int' } )

    b = pd.merge (np_vaf,  df, left_on = "id", right_on = "id")

    # mixture_quantumclone, membership_quantumclone 정해주기
    print ( "원래 clustering number : {}".format ( set (b["cluster"])) )
    
    membership_quantumclone =  list(b["cluster"] )        # quantumclone은 cluster 가 1부터 시작하니까.. 
    min_membership_quantumclone = np.min( membership_quantumclone )  
    for k in range( len( membership_quantumclone) ):
        membership_quantumclone[k] = membership_quantumclone[k] - min_membership_quantumclone
    
    dict = {}
    dict_rev = {}
    #print ( "중간 clustering number : {}".format ( set ( membership_quantumclone ) ) )
    
    if len (set(membership_quantumclone))  != (np.max (list (set (membership_quantumclone))) + 1) :  #중간에 붕 떠버리는 경우가 잇따  {0, 1,2, 3, 4, 6, 7, 8, 9}
        print ("Hi")
        for j1, j2 in enumerate ( list(set(membership_quantumclone)) ):
            dict [j2] = j1
            dict_rev [j1] = j2

        new_membership_quantumclone = []
        for k in range ( len(membership_quantumclone) ) :
            new_membership_quantumclone.append ( dict [  membership_quantumclone[k]  ] )
        membership_quantumclone = copy.deepcopy ( new_membership_quantumclone )
    else:
        for j1, j2 in enumerate ( list(set(membership_quantumclone)) ):
            dict[j2] = j1
            dict_rev[j1] = j2

    print ( "조절한 다음 clustering number : {}".format( set (membership_quantumclone)) )


    mixture_quantumclone = np.zeros ((mixture_answer.shape[0], len(set(membership_quantumclone))), dtype = "float")

    if mixture_answer.shape[0] == 3:
        for j in set(membership_quantumclone):
            x_mean = np.mean (b[b["cluster"] == dict_rev[j] + min_membership_quantumclone]["block0"])
            y_mean = np.mean (b[b["cluster"] == dict_rev[j] + min_membership_quantumclone]["block1"])
            z_mean = np.mean (b[b["cluster"] == dict_rev[j] + min_membership_quantumclone]["block2"])
            mixture_quantumclone[:, j ] = [x_mean * 2, y_mean * 2, z_mean * 2]
    elif mixture_answer.shape[0] == 2:
        for j in set(membership_quantumclone):
            x_mean = np.mean (b[b["cluster"] == dict_rev[j] + min_membership_quantumclone]["block0"])
            y_mean = np.mean (b[b["cluster"] == dict_rev[j] + min_membership_quantumclone]["block1"])
            mixture_quantumclone[:, j ] = [x_mean * 2, y_mean * 2]
    elif mixture_answer.shape[0] == 1:
        for j in set(membership_quantumclone):
            x_mean = np.mean (b[b["cluster"] == dict_rev[j] + min_membership_quantumclone]["block0"])
            mixture_quantumclone[:, j ] = [x_mean * 2]


    # 채점하기
    if kwargs["SCORING"] == True:
        score_df, score = \
            scoring.mixturebased(mixture_answer, mixture_quantumclone, membership_answer, membership_quantumclone, samplename_dict_input, samplename_dict_input_rev, "No", -1, "QuantumClone", **kwargs) # FP는 무조건 못 잡을테니 -1로 넘겨준다
        max_score, sample_dict_PtoA, sample_dict_AtoP = scoring.Scoring ( membership_answer, membership_answer_numerical, membership_quantumclone, -1, []  ) # fp를 designate 하지 못하니까 무조건 fp_index는 -1  , parent_index는 []

        mixture_quantumclone = np.round(mixture_quantumclone, 2)
        return score_df, score, max_score, membership_quantumclone, mixture_quantumclone, sample_dict_PtoA, sample_dict_AtoP 

    else:
        return pd.DataFrame(), -1, -1, membership_quantumclone, mixture_quantumclone, {}, {}



