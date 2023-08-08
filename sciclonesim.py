import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import palettable
import scoring

def main (INPUT_SCICLONE_RESULT, INPUT_NPVAF, OUTPUT_FILENAME,  mixture_answer, membership_answer, membership_answer_numerical, samplename_dict_input, samplename_dict_input_rev, **kwargs):   
    #  INPUT_SCICLONE_RESULT =  "/data/project/Alzheimer/EM_cluster/old/pilot/04.EM_input/sciclone/result/results.tsv"
    #  OUTPUT_FILENAME = "./output/MRS_sciclone.jpg" 

    df = pd.read_csv (INPUT_SCICLONE_RESULT, sep = "\t")

    vaf_col = [i for i in df.columns if "vaf" in i]

    df.loc[:,vaf_col] =  df.loc[:,vaf_col] / 100
    df["mutation_id"] = df["chr"] + "_" + df.astype({"st":"str"})["st"]

    # 각 mutation_id의 vaf 정보를 가져오기 위함
    np_vaf = pd.read_csv(INPUT_NPVAF, sep = "\t")
    np_vaf.rename(columns = {"Unnamed: 0":"mutation_id"}, inplace = True) 
    np_vaf.head()

    b = pd.merge (np_vaf,  df, left_on = "mutation_id", right_on = "mutation_id")

    # mixture_sciclone, membership_sciclone 정해주기
    membership_sciclone =  list(b["cluster"])        # sciclone은 cluster 가 1부터 시작하니까.. 
    
    print ( "원래 clustering number : {}".format ( set (membership_sciclone)) )
    min_membership_sciclone = np.min( membership_sciclone )   # 원래는 이 값이 1인게 정상. 그런데 가끔 0일 수 있다
        
    for k in range( len( membership_sciclone) ):
        membership_sciclone[k] = membership_sciclone[k] - min_membership_sciclone
    mixture_sciclone = np.zeros ((mixture_answer.shape[0], len(set(membership_sciclone))), dtype = "float")
    
    print ( "조절한 다음 clustering number : {}".format( set (membership_sciclone)) )

    if mixture_answer.shape[0] == 3:
        for j in set(membership_sciclone):
            x_mean = np.mean (b[b["cluster"] == j + min_membership_sciclone]["block0"])
            y_mean = np.mean (b[b["cluster"] == j + min_membership_sciclone]["block1"])
            z_mean = np.mean (b[b["cluster"] == j + min_membership_sciclone]["block2"])
            mixture_sciclone[:, j ] = [x_mean * 2, y_mean * 2, z_mean * 2]
    elif mixture_answer.shape[0] == 2:
        for j in set(membership_sciclone):
            x_mean = np.mean (b[b["cluster"] == j + min_membership_sciclone]["block0"])
            y_mean = np.mean (b[b["cluster"] == j + min_membership_sciclone]["block1"])
            mixture_sciclone[:, j ] = [x_mean * 2, y_mean * 2]
    elif mixture_answer.shape[0] == 1:
        for j in set(membership_sciclone):
            x_mean = np.mean (b[b["cluster"] == j + min_membership_sciclone]["block0"])
            mixture_sciclone[:, j ] = [x_mean * 2]


    print (mixture_sciclone)

    # 채점하기
    if kwargs["SCORING"] == True:
        
        score_df, score = \
            scoring.mixturebased(mixture_answer, mixture_sciclone, membership_answer, membership_sciclone, samplename_dict_input, samplename_dict_input_rev, "No", -1,  "SciClone", **kwargs)  # FP는 무조건 못 잡을테니 -1로 넘겨준다
        try:
            max_score, sample_dict_PtoA, sample_dict_AtoP  = scoring.Scoring ( membership_answer, membership_answer_numerical, membership_sciclone, -1 , [] ) # fp를 designate 하지 못하니까 무조건 fp_index는 -1, parent_index는 []
        except:
            print ("SciClone 뭔가 이상하다 membership_answer = {}\nmembership_sciclone = {}".format(membership_answer, membership_sciclone))
            return score_df, score, score, membership_sciclone, mixture_sciclone, {}, {}

        mixture_sciclone = np.round(mixture_sciclone, 2)
        return score_df, score, max_score, membership_sciclone, mixture_sciclone,  sample_dict_PtoA, sample_dict_AtoP         # 100점 만점 score 반환

    else:
        return pd.DataFrame(), -1, -1, membership_sciclone, mixture_sciclone, {}, {}





if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description='Here is usage direction.')
    parser.add_argument('--INPUT_SCICLONE_RESULT')
    parser.add_argument('--INPUT_NPVAF')
    parser.add_argument('--OUTPUT_FILENAME')
    args = parser.parse_args()

    main (args.INPUT_SCICLONE_RESULT, args.INPUT_NPVAF, args.OUTPUT_FILENAME )
