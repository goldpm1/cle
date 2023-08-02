import EMsoft, EMhard
import numpy as np
import math

def main (NUM_MUTATION, df, mixture_total, child_clone, **kwargs):

    #1. 다시 E step을 돌려서 likelihood를 평가하고 membership을 새로 줌
    likelihood_onemore, membership_onemore, membership_p_onemore = EMsoft.Estep (mixture_total [:, :-1], df)
    NUM_BLOCK = mixture_total.shape[0]
    NUM_CLONE = mixture_total.shape[1]

    kwargs["ELBOW_S"] = 3; kwargs["GAUSSIAN_SD"] = 1.2


    #2. Outlier를 추출
    #개별 cluster 단위로 outlier 추출
    outlier_index1 = []
    for j in range(NUM_CLONE):
        if j in child_clone:        # Child clone에서만 FP를 뽑도록 조작하는 것
            sel_index = np.where (np.array (membership_onemore ) == j )[0]
            x, y = np.arange(0, len(sel_index), dtype="int32"),  np.zeros(len(sel_index), dtype="float64")
            for k_i, k in enumerate(sel_index):  # p 중에 그래도 가장 큰 것을 골라서 준다
                y[k_i] = np.max(membership_p_onemore[k])         #  그건 너무 가혹하니까...    (y[k] = np.power(10, np.max(membership_p_onemore[k])) )
            y_sorted = sorted(y, reverse=True)

            threshold_x, threshold_y = EMsoft.Outlier_detection(x, y_sorted, "decreasing", kwargs["ELBOW_S"],  3, 0, "NUM_CLONE = {}  CLONE = {}  GAUSSIAN_SD = {}".format(NUM_CLONE, j, kwargs["GAUSSIAN_SD"]), "looser")
            for k_i, k in enumerate(sel_index):
                if y[k_i] < threshold_y:
                    outlier_index1.append (k)
    outlier_index1 = sorted(outlier_index1)

    # 전체 단위로 outlier 추출
    x, y = np.arange(0, NUM_MUTATION, dtype="int32"),  np.zeros(NUM_MUTATION, dtype="float64")
    for k in range(NUM_MUTATION):  # p 중에 그래도 가장 큰 것을 골라서 준다
        y[k] = np.max(membership_p_onemore[k])         #  그건 너무 가혹하니까...    (y[k] = np.power(10, np.max(membership_p_onemore[k])) )
    y_sorted = sorted(y, reverse=True)
    threshold_x, threshold_y = EMsoft.Outlier_detection(x, y_sorted, "decreasing", kwargs["ELBOW_S"], kwargs["GAUSSIAN_SD"] , 0, "NUM_CLONE = {}   GAUSSIAN_SD = {}".format(NUM_CLONE, kwargs["GAUSSIAN_SD"]), "looser")
    outlier_index2 = sorted((list ( np.where (np.array (y) < threshold_y) [0] )))


    outlier_index = list (set(outlier_index1) & set(outlier_index2))
    print ("threshold_x : {}\noutlier_index : {}".format(threshold_x, outlier_index))

    if len(outlier_index) == 0:
        includeoutlier_total = "No"
    else:
        includeoutlier_total = "Yes"
        # 다시 뽑은 outlier를 mixture, membership에 반영하기
        mixture_total, membership_onemore, membership_p_onemore = EMsoft.Outlier_resetting(outlier_index, y, df, mixture_total [:, :-1], membership_onemore, membership_p_onemore)
        np.unique(membership_onemore, return_counts = True)



    # Soft clustring으로 끝내기

    mixture_onemore_soft = EMsoft.Mstep (membership_onemore, membership_p_onemore, df,  mixture_total.shape[1] - (includeoutlier_total == "Yes"), mixture_total.shape[0], len(df) , list(child_clone), "Soft")  #

    mixture_onemore_soft = np.append(mixture_onemore_soft, np.zeros((NUM_BLOCK, 1), dtype="float64"), axis=1)
    if includeoutlier_total == "Yes":             # 맨 끝 column에 outlier mixture를 재정비해주기 (M step)
        for i in range (kwargs["NUM_BLOCK"]):
            sum_depth, sum_alt = 0, 0
            for k in range(kwargs["RANDOM_PICK"]):
                if membership_onemore[k] == np.max(membership_onemore):    # outlier 찾기
                    sum_depth = sum_depth + df[k][i]["depth"]
                    sum_alt = sum_alt + df[k][i]["alt"]
            mixture_onemore_soft[i][-1] = round((sum_alt * 2) / sum_depth, 2)



    # Hard clustring으로 끝내기

    mixture_onemore_hard = EMsoft.Mstep (membership_onemore, membership_p_onemore, df,  mixture_total.shape[1] - (includeoutlier_total == "Yes"), mixture_total.shape[0], len(df) , list(child_clone), "Hard")  

    mixture_onemore_hard = np.append(mixture_onemore_hard, np.zeros((NUM_BLOCK, 1), dtype="float64"), axis=1)
    if includeoutlier_total == "Yes":             # 맨 끝 column에 outlier mixture를 재정비해주기 (M step)
        for i in range (kwargs["NUM_BLOCK"]):
            sum_depth, sum_alt = 0, 0
            for k in range(kwargs["RANDOM_PICK"]):
                if membership_onemore[k] == np.max(membership_onemore):    # outlier 찾기
                    sum_depth = sum_depth + df[k][i]["depth"]
                    sum_alt = sum_alt + df[k][i]["alt"]
            mixture_onemore_hard[i][-1] = round((sum_alt * 2) / sum_depth, 2)



    return mixture_onemore_soft,  mixture_onemore_hard, membership_onemore, membership_p_onemore 