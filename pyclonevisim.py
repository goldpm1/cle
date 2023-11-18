def main (INPUT_PYCLONEVI_RESULT, INPUT_NPVAF, OUTPUT_FILENAME,  mixture_answer, membership_answer, membership_answer_numerical, **kwargs):   
    import pandas as pd
    import numpy as np
    import scoring

    samplename_dict_input, samplename_dict_input_rev = kwargs["samplename_dict_CharacterToNum"], kwargs["samplename_dict_NumToCharacter"]

    df = pd.read_csv (INPUT_PYCLONEVI_RESULT, sep = "\t")
    df = df.drop_duplicates(['mutation_id'], keep = 'first')  
    

    # 각 mutation_id의 vaf 정보를 가져오기 위함
    np_vaf = pd.read_csv(INPUT_NPVAF, sep = "\t")
    np_vaf.rename(columns = {"Unnamed: 0":"mutation_id"}, inplace = True) 
    np_vaf.head()

    b = pd.merge (np_vaf,  df, left_on = "mutation_id", right_on = "mutation_id")


    # mixture_pyclonevi, membership_pyclonevi 정해주기
    membership_pyclonevi =  list(b["cluster_id"])        # pyclonevi은 cluster 가 1부터 시작하니까.. 
    mixture_pyclonevi = np.zeros ((mixture_answer.shape[0], len(set(membership_pyclonevi))), dtype = "float")

    if mixture_answer.shape[0] == 3:
        for j in set(membership_pyclonevi):
            x_mean = np.mean (b[b["cluster_id"] == j ]["block0"])
            y_mean = np.mean (b[b["cluster_id"] == j ]["block1"])
            z_mean = np.mean (b[b["cluster_id"] == j ]["block2"])
            mixture_pyclonevi[:, j ] = [x_mean * 2, y_mean * 2, z_mean * 2]
    elif mixture_answer.shape[0] == 2:
        for j in set(membership_pyclonevi):
            x_mean = np.mean (b[b["cluster_id"] == j ]["block0"])
            y_mean = np.mean (b[b["cluster_id"] == j ]["block1"])
            mixture_pyclonevi[:, j ] = [x_mean * 2, y_mean * 2]
    elif mixture_answer.shape[0] == 1:
        for j in set(membership_pyclonevi):
            x_mean = np.mean (b[b["cluster_id"] == j ]["block0"])
            mixture_pyclonevi[:, j ] = [x_mean * 2]



    # 채점하기
    if kwargs["SCORING"] == True:
        #score_df, score = scoring.mixturebased(mixture_answer, mixture_pyclonevi, membership_answer, membership_pyclonevi, samplename_dict_input, samplename_dict_input_rev, "No", -1,  "PyCloneVI", **kwargs)   # FP는 무조건 못 잡을테니 -1로 넘겨준다
        max_score, sample_dict_rev, sample_dict, score_df = scoring.Scoring ( membership_answer, membership_answer_numerical, membership_pyclonevi, -1, [] , **kwargs  ) # fp를 designate 하지 못하니까 무조건 fp_index는 -1, parent_index는 []

        mixture_pyclonevi = np.round(mixture_pyclonevi, 2)
        return score_df, max_score, max_score, membership_pyclonevi, mixture_pyclonevi, sample_dict_rev, sample_dict 

    else:
        return pd.DataFrame(), -1, -1, membership_pyclonevi, mixture_pyclonevi, {}, {}