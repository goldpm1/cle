from collections import Counter
import scipy as sp
import scipy.stats
from scipy.special import beta, gamma
import math, visualizationsingle, visualizationsinglesoft, EMhard
import numpy as np
import pandas as pd
import random, palettable, isparent, isparent
import EMhard, EMsoft, Estep, Mstep, Bunch
from kneed import KneeLocator
from scipy.ndimage import gaussian_filter
import matplotlib.pyplot as plt
from sklearn.metrics import silhouette_samples, silhouette_score
from sklearn.cluster import KMeans
import warnings
warnings.simplefilter (action = 'ignore', category = FutureWarning)
warnings.filterwarnings("ignore")



def Outlier_detection(x, y, **kwargs):
    warnings.simplefilter (action = 'ignore', category = FutureWarning)
        
    if len(x) <= 2:            # 해당 clone에서의 원소의 개수가 1, 2개면 에러 난다
        threshold_x = len(x)
        threshold_y = float("-inf")

    else:
        y_prime = np.array(y)[1:] - np.array(y)[:-1]                       # 1차 미분값
        y_prime = np.append(y_prime, y_prime[-1])
        smoothed = scipy.ndimage.gaussian_filter(y_prime, 3.)
        smoothed_sd = scipy.ndimage.standard_deviation(y_prime[int(len(y_prime)*0.1):int(len(y_prime)*0.99)])
        y_doubleprime = np.array(smoothed)[1:] - np.array(smoothed)[:-1]
        y_doubleprime = np.append(y_doubleprime, y_doubleprime[-1])
        Median = np.percentile(y_prime, 50,  interpolation='midpoint')   # median
    
        
        if kwargs["VERBOSE"] >= 1:
            fig, ax = plt.subplots(ncols=2, figsize=(15, 6))
            fig.suptitle( kwargs["FIG_TITLE"] )
            ax[0].set_title ("y = Posterior ")
            ax[1].set_title ("y' = prime of y")
            ax[0].plot(x[:-1], y[:-1], '.', c="g", alpha=0.3)
            ax[1].plot(x[:int(len(smoothed) * 0.99 )], smoothed[:int(len(smoothed) * 0.99 )] , '.', c="r", alpha=0.3)



        #print ("Clone No : {0}\nprob_list : {1},\ny : {2}".format(j,max_prob_list, y))

        # 후보1.  Elbow 법
        threshold_x_elbow, threshold_y_elbow = len(x), float("-inf")
        threshold_x_gaussian, threshold_y_gaussian = len(x), float("-inf")
        for S in list(np.arange( kwargs["ELBOW_S"] , 0, -0.4)):
            kneedle = KneeLocator(x, y, S=S, curve="concave", direction= kwargs["DIRECTION"] )                 # " increasing", "decreasing"
            if kneedle.knee != None:
                threshold_x_elbow = round(kneedle.knee, 1)
                threshold_y_elbow = round(kneedle.knee_y, 1)        # elbow 값을 구함
                if kwargs["VERBOSE"] >= 2:
                    print("S = {0}".format(S))
                if kwargs["VERBOSE"] >= 1:
                    ax[0].axvline(threshold_x_elbow, label="knee/elbow", color="b")
                    ax[1].axvline(threshold_x_elbow, label="knee/elbow", color="b")
                    #kneedle.plot_knee()
                break

        # 후보2.  Gaussian filter를 통해 찾아보자
        #print ("y_prime : {}\nsmoothed : {}\ny_doubleprime : {}".format(y_prime, smoothed,  y_doubleprime) )
        for i in range(int(len(smoothed) / 3), len(smoothed)):
            if (smoothed[i] < Median - ( kwargs["GAUSSIAN_SD"] * smoothed_sd)) & (y_doubleprime[i] < 0):
                convexdetector = 0           # y'' 를 주위로 1%씩 분석해봐서 대부분 concave여야 통과해준다

                width = 0
                for j in range( np.max( [ 0, int(i - len(y) / 100) ] ), np.min ( [ len(y_doubleprime) - 1, int(i + len(y) / 100) ] )):
                    width = width + 1
                    if y_doubleprime[j] > 0:                  # convex의 개수 세기
                        convexdetector = convexdetector + 1
            
                    
                if convexdetector < width / 4:           # convex < 25%
                    threshold_x_gaussian = i
                    threshold_y_gaussian = round(y[i], 1)
                    if kwargs["VERBOSE"] >= 1:
                        #print ("i : {0} , smoothed[i] : {1},  median : {2},  smoothed_sd : {3}, y_doubleprime : {4}". format(i, smoothed[i], Median, smoothed_sd, y_doubleprime))
                        ax[0].axvline(threshold_x_gaussian, label="Gaussian filter", color="g")
                        ax[1].axvline(threshold_x_gaussian, label="Gaussian filter", color="g")
                        ax[1].axhline(Median , label="y' median", color="#90a4ae")
                        ax[1].axhline(Median - (kwargs["GAUSSIAN_SD"] * smoothed_sd), label="y' cutoff", color="#90a4ae")
                        #ax[1].axis([0,  len(smoothed), np.min(smoothed) * 1.1,  0])
                        ax[0].legend()
                        ax[1].legend()
                break

        # 둘 중 threshold가 더 엄격한 것을 골라주자
        if kwargs ["OUTLIER_STANDARD"] == "stricter":
            if threshold_x_elbow < threshold_x_gaussian:
                threshold_x, threshold_y = threshold_x_elbow, threshold_y_elbow
            else:
                threshold_x, threshold_y = threshold_x_gaussian, threshold_y_gaussian
        else:
            if (threshold_x_elbow == len(x))  & (threshold_x_gaussian != len(x)) :        # elbow법으로 못 찾았을 떄에는 그냥 gaussian을 주자
                threshold_x, threshold_y = threshold_x_gaussian, threshold_y_gaussian
            elif (threshold_x_elbow != len(x))  & (threshold_x_gaussian == len(x)) :        # gaussian법으로 못 찾았을 떄에는 그냥 elbow를 주자
                threshold_x, threshold_y = threshold_x_elbow, threshold_y_elbow
            else:
                if threshold_x_elbow > threshold_x_gaussian:
                    threshold_x, threshold_y = threshold_x_elbow, threshold_y_elbow
                else:
                    threshold_x, threshold_y = threshold_x_gaussian, threshold_y_gaussian

        
        # # 조금 더 느슨하게 해 주고 싶어서...
        # if ISPARENT == False:
        #     if threshold_x < len(x):
        #         if kwargs["VERBOSE"] >= 2:
        #             ax[0].axhline(y[threshold_x - int(len(x) / 50)]  , color="#90a4ae")
        #             ax[0].axhline(threshold_y  , color="g")
        #         threshold_y = threshold_y - (y[threshold_x - int(len(x) / 50)] - threshold_y)
        #         for k in range (threshold_x, len(x)):
        #             if y[k] < threshold_y:                  
        #                 if kwargs["VERBOSE"] >= 2:
        #                     ax[0].axvline(threshold_x, label="Final", color="#ff7043")
        #                     ax[1].axvline(threshold_x, label="Final", color="#ff7043")
        #                     ax[0].axhline(y[k] , color="#ff7043")
        #                     ax[0].legend()
        #                     ax[1].legend()
        #                 threshold_x = k
        #                 threshold_y = y[k]
        #                 break

        #         plt.show()

        if kwargs["VERBOSE"] >= 2:
            print ("\t\t SUPTITLE : {0} \t\t threshold_x : {1} \t\t threshold_y : {2} \t\t  len(y) = {3}".format( kwargs["FIG_TITLE"], threshold_x, threshold_y, len(y)))
    
    if kwargs["FIG_DIRECTORY"] != "NotSave":
        fig.savefig( kwargs["FIG_DIRECTORY"] )

    return threshold_x, threshold_y




def elaboration_outlierclone (cluster, df, np_vaf, **kwargs):
    NUM_MUTATION = kwargs["NUM_MUTATION"]
    NUM_BLOCK = kwargs["NUM_BLOCK"]

    # 이미 kwargs["NUM_CLONE"]은 하나 줄이고 온 상황

    nuts = Bunch.Bunch1(NUM_MUTATION , NUM_BLOCK, kwargs["NUM_CLONE"], kwargs["STEP_NO"])

    # 1. outlier clone을 제외한 NUM_CLONE  - 1로 E step 돌고 그림 그리기

    kwargs["STEP"] = 0
    nuts.mixture =  np.delete (cluster.mixture_record[ kwargs["NUM_CLONE_NOMINAL"] ], cluster.fp_index_record[ kwargs["NUM_CLONE_NOMINAL"] ], axis = 1)
    nuts = Estep.main(df, np_vaf, nuts, **kwargs)                   # 주어진 mixture 내에서 새 membership 정하기
    kwargs["OPTION"] , kwargs["STEP"], nuts.likelihood, nuts.includeoutlier  = "start", "nuts", 0, False        # 일단은 n-1개로 그림 그리는 용도니까
    Mstep.drawfigure_2d (nuts, np_vaf, kwargs["MYEM_DIR"] + "/trial/clone" + str (kwargs["NUM_CLONE_NOMINAL"]) + ".nuts.pre.jpg" , **kwargs)



    # 2. 위 조건에서 Outlier 고르기  (전체에서 고르고, 원래 outlier clone과 교집합을 구해줌)
    x = np.arange(0, NUM_MUTATION, dtype="int32")
    y = np.amax(nuts.membership_p, axis=1)
    y_sorted = sorted ( y , reverse = True)

    kwargs ["DIRECTION"] = "decreasing"
    kwargs ["OUTLIER_STANDARD"] = "looser"
    kwargs ["GAUSSIAN_SD"] = kwargs["GAUSSIAN_SD"] + 0.15 * ( kwargs["NUM_CLONE"])
    kwargs ["FIG_TITLE"] = "NUM_CLONE = {}    GAUSSIAN_SD = {}    {}".format(kwargs["NUM_CLONE_NOMINAL"], kwargs["GAUSSIAN_SD"], kwargs["OUTLIER_STANDARD"])
    kwargs ["FIG_DIRECTORY"] = kwargs["MYEM_DIR"] + "/trial/clone" + str (kwargs["NUM_CLONE_NOMINAL"]) + ".nuts.outlier.cutoff.jpg"
     
    threshold_x, threshold_y = Outlier_detection(x, y_sorted, **kwargs)
    outlier_index_1 = sorted((list ( np.where (np.array (y) < threshold_y) [0] )))
    outlier_index_2 = np.where ( np.array( cluster.membership_record [kwargs["NUM_CLONE_NOMINAL"] ] ) == cluster.fp_index_record [kwargs["NUM_CLONE_NOMINAL"]]  ) [0]
    outlier_index_common = sorted( list ( set (outlier_index_1) & set(outlier_index_2) ) )
        
    #print ("공통 : {}  ({}개)".format (set (outlier_index_1) & set(outlier_index_2),  len( set (outlier_index_1) & set(outlier_index_2) )))



    nuts.includeoutlier = True if  len(outlier_index_common) > 0 else False

    #3. Mixture 재조절하고 정리해서 cluster로 넘기기
    if nuts.includeoutlier == True:         #  M step + Outlier로 분류된 애들을 재정비해서 mixture 만들어주기
        kwargs ["STEP"], kwargs["STEP_TOTAL"], kwargs["TRIAL"], nuts.includeoutlier  = 0, cluster.stepindex, cluster.trialindex, False
        nuts = Mstep.main(df, np_vaf, nuts, "Hard", **kwargs)   # 새 memberhsip에서 새 mixture구하기
        nuts.includeoutlier = True

        # outlier mixture 만들어주기
        nuts.mixture = np.append (nuts.mixture, np.zeros ( (kwargs["NUM_BLOCK"], 1) ), axis = 1)
        for i in range (kwargs["NUM_BLOCK"]):
            sum_depth, sum_alt = 0, 0
            for k in outlier_index_common:
                sum_depth = sum_depth + df[k][i]["depth"]
                sum_alt = sum_alt + df[k][i]["alt"]
            nuts.mixture [i][-1] = round((sum_alt * 2) / sum_depth, 2)

        nuts.outlier_index = outlier_index_common           # mutation 단위의 index
        nuts.fp_index = np.max(nuts.membership) + 1
        nuts.membership [outlier_index_common] = nuts.fp_index

        #print ("Post      FP_index : {}	Child_index : {}".format(nuts.fp_index, nuts.makeone_index  ))

        nuts = Estep.main(df, np_vaf, nuts, **kwargs)                   # 주어진 mixture 내에서 새 membership 정하기

        # outlier membership_p만들어주기  (1칸만 추가해주면 됨)
        nuts.membership_p = np.append (nuts.membership_p, np.full((kwargs["NUM_MUTATION"], 1), -999),  axis = 1)


    print ("nuts : .fp_index : {}\t.makeone_index : {}\t.includeoutlier : {}\t.outlier_index : {}\t.likelihood : {}".format( nuts.fp_index, nuts.makeone_index, nuts.includeoutlier, nuts.outlier_index, nuts.likelihood  ))


    kwargs["STEP"], kwargs["STEP_TOTAL"] = cluster.stepindex, cluster.stepindex
    kwargs ["OPTION"] = "Outlier"
    Mstep.drawfigure_2d (nuts, np_vaf, kwargs["MYEM_DIR"] + "/trial/clone" + str (kwargs["NUM_CLONE_NOMINAL"]) + ".nuts.outlier.post.jpg" , **kwargs)
    Mstep.drawfigure_2d (nuts, np_vaf, kwargs["MYEM_DIR"] + "/trial/clone" + str (kwargs["NUM_CLONE_NOMINAL"]) + "." + str(cluster.trialindex) + "-" + str(cluster.stepindex) + "(hard).jpg" , **kwargs)


    cluster.stepindex = cluster.stepindex + 1
    cluster.acc ( nuts.mixture, nuts.membership, nuts.likelihood, nuts.membership_p, nuts.membership_p_normalize, cluster.stepindex, 
                        cluster.trialindex, nuts.makeone_index, nuts.fp_index, True, nuts.outlier_index, **kwargs )      # 원래의 NUM_CLONE에 nuts 정보를 집어넣기 (IncludeOutlier : True)

    return nuts
