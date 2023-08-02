from collections import Counter
import scipy
import scipy.stats
from scipy.special import beta, gamma
import math
import numpy as np
import pandas as pd
import random
from kneed import KneeLocator
from scipy.ndimage import gaussian_filter
import matplotlib.pyplot as plt


def set_initial_parameter(NUM_CLONE, trial, mixture):
    while True:
        for i in range (NUM_BLOCK):
            li = list(range(NUM_CLONE))   # 분율이 0인 경우도 포함
            choiceLIst = Counter([random.choice(li) for i in range(10)])
            
            for j in range (NUM_CLONE):
                mixture[i][j] = choiceLIst[j] / 10

        break
        # if 0 not in mixture:   # 처음부터 분율이 0인게 있으면 다시 돌려라
        #     break
      

    # print ("INITIAL : trial " + str(trial + 1) )
    # for i in range (NUM_BLOCK):
    #    print ("\t{0}".format(mixture[i]))

    return mixture


def calc_likelihood (mixture, df,  k, NUM_CLONE, NUM_BLOCK):
    max_prob = -99999; max_clone = -1

    for j in range (NUM_CLONE):
        prob = 0; prob2 = 0
        for i in range (NUM_BLOCK):
            depth_calc = int(df[k][i]["depth"] * mixture[i][j])
            alt_calc = int(df[k][i]["depth"] * mixture[i][j] * 0.5)
            depth_obs = int(df[k][i]["depth"])
            alt_obs = int(df[k][i]["alt"])
            
            # Binomial probability
            # if depth_calc >= alt_obs:
            #     try:
            #         prob = prob + math.log10(scipy.stats.binom.pmf (n = depth_calc,k = alt_obs, p = 0.5))            # 이렇게 하면 안되고 beta binomial로 계산해야 함  ( Beta  (alt_expected, depth - alt_expected , vaf_observed))
            #     except:
            #         prob = prob - 400
            #     #math.log10(sys.float_info.min*sys.float_info.epsilon)
            # else :
            #     prob = prob +  math.log10(SEQ_ERROR) * (alt_obs-depth_calc)


        
            # # Beta binomial distribution
            a = df[k][i]["depth"] * mixture[i][j] * 0.5              # alt_expected
            b = depth_obs - a            # ref_expected
            try:
                prob = prob + math.log10(scipy.stats.betabinom.pmf(alt_obs, depth_obs, a+1, b+1))
            except:
                prob = prob - 400


            #print ("{0}번째 mutation : {1} 번째 clone, {2}번째 block → alt_expected : {3},  alt_observed : {4}, depth_observed : {5}, likelihood : {6}"\
            #         .format(k, j, i, round(a,1), alt_obs, depth_obs,  scipy.stats.betabinom.pmf(alt_obs, depth_obs, a+1, b+1) ))
        #print ("{0}번째 mutation : {1} 번째 clone →  likelihood : {2}".format(k, j,   prob ))
                  
        if prob > max_prob:
            max_prob = prob
            max_prob_clone_candidate = [j]
        elif prob == max_prob:
            max_prob_clone_candidate.append(j)
    
    max_clone = random.choice(max_prob_clone_candidate)
    #print ("mutation_{0} : {1}번 clone이 가장 적절 (log p = {2})\n".format(k, max_clone, max_prob))
    return max_clone, max_prob


def Outlier_detection (membership, max_prob_list, NUM_CLONE):
    
    for j in range(NUM_CLONE):
        membership_j_index = list(np.where(membership == j)[0])         # j 번째 clone의 mutation만 (i)을 구함
        if len (membership_j_index) <= 2 :     # 해당 membership의 원소가 너무 작으면 y_prime, elbow 계산도 안 된다
            continue

        max_prob_j = []
        for i in membership_j_index:
            max_prob_j.append (max_prob_list[i])                                       # i번째 mutation의 probabilyt를 list화시킴
        x= np.arange(start = 1, stop = len(max_prob_j) + 1 )
        y = sorted(max_prob_j, reverse = True) 
        y_prime = np.array(y)[1:] - np.array(y)[:-1]
        smoothed = scipy.ndimage.gaussian_filter(y_prime, 3.)
        smoothed_sd = scipy.ndimage.standard_deviation(y_prime)
        Median = np.percentile(y_prime, 50,  interpolation = 'midpoint')   # median

        #print ("CLONE = {0} \t\t {1} \t {2}".format(j, list(x), y))

        #1. 먼저 Gaussian filter를 통해 찾아보자
        threshold_x = -1
        for i in range(len(smoothed)):
            if y_prime[i] < Median - (1 * smoothed_sd):
                threshold_x = i
                threshold_y = round (y[i], 1)
                plt.plot(x[:-1] , y[:-1], '.', label = "CLONE : " + str(j))
                plt.vlines (threshold_x, ymin = plt.axis()[2],  ymax = plt.axis()[3])
                break

        # 2. 안 되면  Elbow 법이라도..
        if threshold_x == -1:
            for S in list(np.arange(3, 0, -0.4)):
                kneedle = KneeLocator(x, y, S=S, curve="concave", direction="decreasing")
                if kneedle.knee != None:
                    print ("S = {0}".format(S))
                    threshold_x = round(kneedle.knee, 1)
                    threshold_y = round(kneedle.knee_y, 1)        # elbow 값을 구함
                    kneedle.plot_knee()
                    break

        print ("CLONE = {0}\t\t threshold_x : {1} \t\t threshold_y : {2}".format(j, threshold_x, threshold_y))

        for i in membership_j_index:
            if (max_prob_list[i] < threshold_y):           # elbow 값보다 낮으면
                membership[i] = NUM_CLONE         # 그 i번째 membership은 outlier clone으로 준다
                #print ("{0}번째 mutation : {1} clone -> {2} outlier로 변경   (prob = {3}, threshold = {4})". format(i, membership[i], NUM_CLONE, max_prob_list[i], threshold_y ))
    plt.legend()

    return membership


def Estep( mixture , df, NUM_CLONE, NUM_BLOCK):  ## MEMBERSHIP 정하는 과정
    total_prob = 0
    max_prob_list = []
    membership = np.zeros (NUM_MUTATION, dtype = "int32")

    for k in range (NUM_MUTATION):
        #print ("{0}번째 mutation : {1}, {2}".format(k, round(df[k][0]["alt"] / df[k][0]["depth"], 2) , round(df[k][1]["alt"] / df[k][1]["depth"], 2)))
        j, max_prob = calc_likelihood (mixture, df,  k, NUM_CLONE, NUM_BLOCK)  # k번째 mutation은 membership을 누구에게 주는 게 좋을까?
        membership[k] = j
        total_prob = total_prob + max_prob
        max_prob_list.append(max_prob)
        

    #print ("\tTotal likelihood : {0}\n\tMembership now : {1}".format(round(total_prob,2), membership))
    membership = Outlier_detection(membership, max_prob_list, NUM_CLONE)

    return round(total_prob,2), membership



def main(df, initial_parameter, **kwargs):
    global NUM_BLOCK, RANDOM_PICK, NUM_MUTATION, NUM_CLONE_TRIAL_START, NUM_CLONE_TRIAL_END, TRIAL_NO
    global max_likelihood, max_membership, max_mixture;


    NUM_BLOCK = len(df[0])
    RANDOM_PICK = len(df)
    NUM_MUTATION = RANDOM_PICK
    print (NUM_MUTATION)
    
    NUM_CLONE_TRIAL_START = initial_parameter.shape[1]
    NUM_CLONE_TRIAL_END = initial_parameter.shape[1]
    TRIAL_NO = 1


    # df : 2d list + dictionary  (k번째 mutation, i번째 block에서의 depth, ref, alt)
    # mixture : 2d ndarray       (i번째 block에서 j번째 clone의 분율)
    # membership : 1d ndarray    (k번째 mutation이 어떤 clone인지)


    for NUM_CLONE in range (NUM_CLONE_TRIAL_START ,NUM_CLONE_TRIAL_END + 1 ):
        print ("NUM_CLONE = {0}".format(NUM_CLONE))
        mixture = initial_parameter
        membership = np.zeros (NUM_MUTATION, dtype = "int32")




        total_prob, membership = Estep(mixture, df,  NUM_CLONE, NUM_BLOCK )            # 주어진 mixture 내에서 새 membership 정하기. Outlier 도 elbow법으로 정해준다 

        return NUM_CLONE + 1, membership, mixture, total_prob
    
    

