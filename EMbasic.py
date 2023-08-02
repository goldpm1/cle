from collections import Counter
import scipy as sp
import scipy.stats
from scipy.special import beta, gamma
import math
import numpy as np
import pandas as pd
import random
from kneed import KneeLocator
from scipy.ndimage import gaussian_filter
import matplotlib.pyplot as plt


def set_initial_parameter(NUM_CLONE, NUM_MUTATION, df , trial, mixture, adjustment):
    if adjustment in ["True", "T"]:
        while True:
            for i in range (NUM_BLOCK):
                li = list(range(NUM_CLONE))   # 분율이 0인 경우도 포함
                choiceLIst = Counter([random.choice(li) for i in range(10)])
                
                for j in range (NUM_CLONE):
                    mixture[i][j] = choiceLIst[j] / 10

            break
            # if 0 not in mixture:   # 처음부터 분율이 0인게 있으면 다시 돌려라
            #     break
    else:
        for j, k in enumerate(random.sample(range(NUM_MUTATION), NUM_CLONE)):
            for i in range(NUM_BLOCK):
                mixture[i][j] = round(int(df[k][i]["alt"]) / int(df[k][i]["depth"]) * 2, 2)


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


            # print ("{0}번째 mutation : {1} 번째 clone, {2}번째 block → alt_expected : {3},  alt_observed : {4}, depth_observed : {5}, likelihood : {6}"\
            #         .format(k, j, i, round(a,1), alt_obs, depth_obs,  scipy.stats.betabinom.pmf(alt_obs, depth_obs, a+1, b+1) ))
                  
        if prob > max_prob:
            max_prob = prob
            max_prob_clone_candidate = [j]
        elif prob == max_prob:
            max_prob_clone_candidate.append(j)
    
    max_clone = random.choice(max_prob_clone_candidate)
    # print ("mutation_{0} : {1}번 clone이 가장 적절 (log p = {2})\n".format(k, max_clone, max_prob))
    return max_clone, max_prob


def Estep( mixture , df, NUM_CLONE, NUM_BLOCK):  ## MEMBERSHIP 정하는 과정
    total_prob = 0
    max_prob_list = []
    membership = np.zeros (NUM_MUTATION, dtype = "int32")

    for k in range (NUM_MUTATION):
        j, max_prob = calc_likelihood (mixture, df,  k, NUM_CLONE, NUM_BLOCK)  # k번째 mutation은 membership을 누구에게 주는 게 좋을까?
        membership[k] = j
        total_prob = total_prob + max_prob
        max_prob_list.append(max_prob)
        

    #print ("\tTotal likelihood : {0}".format(round(total_prob,2)))

    return round(total_prob,2), membership


def Mstep(membership, df, NUM_CLONE, NUM_BLOCK, adjustment , Last):  ## 새로운 MIXTURE 정하는 과정
    mixture = np.zeros ((NUM_BLOCK, NUM_CLONE), dtype = 'float')     #mixture 값을 일단 초기화
    NUM_MUTATION = len(membership)

    for j in range (NUM_CLONE):       
        ind_list = []          # membership == j 인 index를 구하기
        for k in range(NUM_MUTATION):
            if membership[k] == j:
                ind_list.append(k)

        
        for i in range (NUM_BLOCK):
            sum_depth = 0; sum_alt = 0
            for ind in ind_list:       # depth, alt를 다 더하기
                sum_depth = sum_depth + df[ind][i]["depth"]
                sum_alt = sum_alt + df[ind][i]["alt"]
            
            #print (i, j, sum_depth, sum_alt)
            if sum_depth == 0:      # i번째 block에서 j번째 clone이 아예 없으면 0을 준다
                mixture[i][j] = 0
            else:                   # j번째 clone만을 생각한 이상적인 분율을 일단 assign
                mixture[i][j] = round((sum_alt * 2) / sum_depth,2)
           
     # Block당 Mixture의 합이 1이 되도록 재조정       
    if adjustment in  ["True", "T"]:
        for i in range (NUM_BLOCK):
            sum = 0

            if Last == "Last":
                for j in range (NUM_CLONE - 1):
                    sum = sum + mixture[i][j]
                for j in range (NUM_CLONE - 1):
                    mixture[i][j] = np.round(mixture[i][j] / sum, 2)
            else:
                for j in range (NUM_CLONE):
                    sum = sum + mixture[i][j]
                mixture[i] = np.round(mixture[i] / sum, 2)

    return mixture


def Sorting ( membership, mixture , NUM_CLONE):
    # membership, mixture를 앞에 나오는 순서대로 numbering을 다시 해줌
    mixture_sort = np.zeros ((NUM_BLOCK, NUM_CLONE ), dtype = 'float') 
    membership_sort = []
    dict={}
    num = 0

    # membership을 정렬해서 보내줌
    for k in range(NUM_MUTATION):
        if membership[k] not in dict:
            membership_sort.append (num)
            dict[membership[k]] = num
            num = num + 1
        else:
            membership_sort.append (dict[membership[k]])
    membership_sort = np.array(membership_sort)

    # 그 순서대로 mixture도 정렬해서 보내줌
    for i in range (0, NUM_BLOCK):
        for j in range (0, NUM_CLONE):
            if j in dict:
                mixture_sort[i][dict[j]] = mixture[i][j]

    return membership_sort, mixture_sort


def GoStop(t, total_prob, membership, mixture, previous_No, previous_membership, previous_likelihood ): # 이전 membership과 겹치거나 likelihood가 나아질 기미가 보이지 않으면 stop
    global max_likelihood; global max_membership; global max_mixture; global NUM_MUTATION

    # 방금 step에서 구한게 최적값을 경신했는지 판다
    if total_prob > max_likelihood:
        #print ("\t\t\t\t !!!!!! Renewal : ", total_prob)
        max_likelihood = total_prob
        max_membership = membership
        max_mixture = mixture


    # 첫 5번 정도는 판단 유예
    if t >= previous_No:
        #print ("Concordance score : ", end = " " )
        for p in range(previous_No):       # 지난 5개의 membership을 보고 비교해본다
            #print (str(np.sum(np.equal(membership, previous_membership[p]))) + "점", end = "\t")
            if np.sum(np.equal(membership, previous_membership[p])) >= int(NUM_MUTATION * 0.9):          # membership이 90% 이상 겹치면 stop 해주자
                #print ("\n\n지난 {0} 번째와 membership 변화가 없어서 stop : {1}개 겹침".format(p, np.sum(np.equal(membership, previous_membership[p]))))
                return "Stop", previous_membership, previous_likelihood
        # if (total_prob) < 0.99 * np.max(previous_likelihood):
        #     print ("\n\nLikelihood 기준으로 stop")
        #     return "Stop", previous_membership, previous_likelihood


    # 차곡차곡 채워넣어준다
    previous_membership[t % previous_No] = membership
    previous_likelihood[t % previous_No] = total_prob


    return "Go", previous_membership, previous_likelihood


def Outlier_detection_IQR (prob_list):
    Q1 = np.percentile(prob_list, 25,  interpolation = 'midpoint')
    Q3 = np.percentile(prob_list, 75,  interpolation = 'midpoint')
    IQR = Q3 - Q1
    
    # Upper bound
    upper = list( np.where(prob_list>= (Q3+3*IQR)) [0])
    # Lower bound
    lower = list( np.where(prob_list<= (Q1-3*IQR)) [0])

    return upper + lower, (Q1-3*IQR),  (Q3+3*IQR)


def Outlier_detection (max_prob_list, j):

    x= np.arange(start = 1, stop = len(max_prob_list) + 1 )
    y = sorted(max_prob_list, reverse = True) 

    if len (x) <= 2:            # 해당 clone에서의 원소의 개수가 1, 2개면 에러 난다
        threshold_y  = float("-inf")

    else:
        y_prime = np.array(y)[1:] - np.array(y)[:-1]
        smoothed = scipy.ndimage.gaussian_filter(y_prime, 3.)
        smoothed_sd = scipy.ndimage.standard_deviation(y_prime[0:int(len(y_prime)*0.99)])
        Median = np.percentile(y_prime, 50,  interpolation = 'midpoint')   # median

        #print ("Clone No : {0}\nprob_list : {1},\ny : {2}".format(j,max_prob_list, y))


        #후보1.  Elbow 법
        threshold_x_elbow = len(max_prob_list); threshold_y_elbow  = float("-inf")
        threshold_x_gaussian= len(max_prob_list); threshold_y_gaussian  = float("-inf")
        for S in list(np.arange(3, 0, -0.4)):
            kneedle = KneeLocator(x, y, S=S, curve="concave", direction="decreasing")
            if kneedle.knee != None:
                #print ("S = {0}".format(S))
                threshold_x_elbow = round(kneedle.knee, 1)
                threshold_y_elbow = round(kneedle.knee_y, 1)        # elbow 값을 구함
                #kneedle.plot_knee()
                break

        #후보2.  Gaussian filter를 통해 찾아보자
        for i in range(len(smoothed)):
            if smoothed[i] < Median - (1 * smoothed_sd):
                threshold_x_gaussian = i
                threshold_y_gaussian = round (y[i], 1)
                # plt.plot(x[:-1] , y[:-1], '.', label = "CLONE : " + str(j), c = "g", alpha = 0.3)
                # plt.axvline (threshold_x_gaussian, label = "Gaussian filter", color = "g")
                # plt.legend()
                break

        # 둘 중 threshold가 더 엄격한 것을 골라주자
        if threshold_x_elbow < threshold_x_gaussian:
            threshold_x = threshold_x_elbow
            threshold_y = threshold_y_elbow
        else:
            threshold_x = threshold_x_gaussian
            threshold_y = threshold_y_gaussian


        #print ("\t\t threshold_x : {0} \t\t threshold_y : {1}".format(threshold_x, threshold_y))

    upperlower=[]
    for i in range(len(max_prob_list)):
        if (max_prob_list[i] < threshold_y):           # elbow 값보다 낮으면
            upperlower.append(i)

    return upperlower


def Estep_final ( df, maxmaxmax_mixture , maxmaxmax_membership, j ):  ## MEMBERSHIP 정하는 과정
    prob_list = []
    prob_list_index = []
    
    for k in range (NUM_MUTATION):
        if maxmaxmax_membership[k] != j:
            continue

        prob = 0
        for i in range (NUM_BLOCK):
            # # Beta binomial distribution
            depth_obs = int(df[k][i]["depth"])
            alt_obs = int(df[k][i]["alt"])
            a = df[k][i]["depth"] * maxmaxmax_mixture[i][j] * 0.5              # alt_expected
            b = depth_obs - a            # ref_expected
            try:
                prob = prob + math.log10(scipy.stats.betabinom.pmf(alt_obs, depth_obs, a+1, b+1))
            except:
                prob = prob - 400
            
        prob_list.append(prob)
        prob_list_index.append(k)

    
        # print ("{0}번째 mutation (block1 : {1},  block2 : {2})의 p → {3}    (alt_obs = {4}, alt_exp = {5}, depth = {6})" \
        #     .format(k, np_vaf[k][0] * 2, np_vaf[k][1] * 2, round(prob,2), alt_obs, a, depth_obs))        
    #fp_candidate_prob_list, lowerbound, upperbound = Outlier_detection_IQR (prob_list)
    fp_candidate_prob_list = Outlier_detection(prob_list, j)
    fp_candidate_index = [prob_list_index[i] for i in fp_candidate_prob_list]


    return fp_candidate_index




def main(df, method, adjustment, **kwargs):
    global NUM_BLOCK, RANDOM_PICK, NUM_MUTATION, NUM_CLONE_TRIAL_START, NUM_CLONE_TRIAL_END, TRIAL_NO
    global max_likelihood, max_membership, max_mixture;


    NUM_BLOCK = len(df[0])
    RANDOM_PICK = len(df)
    NUM_MUTATION = RANDOM_PICK
    print ("NUM_MUTATION : {}".format(NUM_MUTATION))

    NUM_CLONE_TRIAL_START = kwargs["NUM_CLONE_TRIAL_START"]
    NUM_CLONE_TRIAL_END = kwargs["NUM_CLONE_TRIAL_END"]
    TRIAL_NO = kwargs["TRIAL_NO"]

    gradient_max = float("-inf")
    previous_maxmax_likelihood = 0


    # df : 2d list + dictionary  (k번째 mutation, i번째 block에서의 depth, ref, alt)
    # mixture : 2d ndarray       (i번째 block에서 j번째 clone의 분율)
    # membership : 1d ndarray    (k번째 mutation이 어떤 clone인지)


    maxmaxmax_likelihood = float("-inf")    # 모든 NUM_CLONE에서 최대값을 갖는 것을 저장
    maxmaxmax_membership = np.zeros (NUM_MUTATION, dtype = "int32") 
    maxmaxmax_NUM_CLONE = 0

    maxmaxmax_mixture_record = []
    maxmaxmax_membership_record = []
    maxmaxmax_likelihood_record = []



    for NUM_CLONE in range (NUM_CLONE_TRIAL_START ,NUM_CLONE_TRIAL_END + 1 ):
        # print ("NUM_CLONE = {0}".format(NUM_CLONE))
        mixture = np.zeros ((NUM_BLOCK, NUM_CLONE), dtype = 'float') 
        membership = np.zeros (NUM_MUTATION, dtype = "int32")

        previous_No = 5
        previous_membership = np.zeros ((previous_No, NUM_MUTATION), dtype = "int32")
        previous_likelihood = np.zeros (previous_No, dtype = 'float'); previous_likelihood.fill(float("-inf"))

        max_likelihood = float("-inf")          # 한 trial 내에서 최대값을 갖는 것을 저장
        max_membership = np.zeros (NUM_MUTATION, dtype = "int32")
        max_mixture = np.zeros ((NUM_BLOCK, NUM_CLONE + 1), dtype = 'float') 

        maxmax_likelihood = float("-inf")       # 한 NUM_CLONE에서 최대값을 갖는 것을 저장
        maxmax_membership = np.zeros (NUM_MUTATION, dtype = "int32")
        maxmax_mixture = np.zeros ((NUM_BLOCK, NUM_CLONE + 1), dtype = 'float') 

        mixture_record = []
        membership_record = []
        likelihood_record = []



        for trial in range(TRIAL_NO):
            mixture = set_initial_parameter(NUM_CLONE, NUM_MUTATION, df, trial, mixture, adjustment)
            max_mixture = mixture
            max_likelihood = float("-inf")
            max_membership = np.zeros (NUM_MUTATION, dtype = "int32")

            # mixture = np.array([[0.3,0.6,0.1],[0.7,0.1,0.2]])

            for t in range(1,30):
                #print ("\n\n #{0}번째 step\n".format(t))


                total_prob, membership = Estep(mixture, df,  NUM_CLONE, NUM_BLOCK )            # 주어진 mixture 내에서 새 membership 정하기
                mixture = Mstep(membership , df, NUM_CLONE, NUM_BLOCK, adjustment, "")        # 새 memberhsip에서 새 mixture구하기
                membership_sort, mixture_sort = Sorting ( membership, mixture, NUM_CLONE )
            
                #print ("E-step >\n\tTotal prob : {0} \t  membership : {1} ".format(total_prob, membership_sort))
                #print ("M-step >\n{0}". format(mixture_sort))

                message, previous_membership, previous_likelihood = GoStop(t, total_prob, membership_sort, mixture_sort, previous_No, previous_membership, previous_likelihood)
                if message == "Stop":
                    break
                
                #print ("지금 step까지 max likelihood : {0}\n지금 step까지 max_membership : {1}".format(max_likelihood, max_membership))


            check = 0
            for j in range(NUM_CLONE):
                if j not in np.unique(max_membership):
                    check = 1
            if check == 1:   # 설정한 NUM_CLONE과 정답의 NUM_CLONE 개수가 맞지 않으면 넘어가자
                maxmax_likelihood = max_likelihood
                maxmax_mixture = max_mixture
                maxmax_membership = max_membership
                continue

            #print ("이번 trial의 Total Likelihood = {0}\n{1}".format(max_likelihood, max_membership))
            #print (max_mixture,"\n")

            if maxmax_likelihood < max_likelihood:
                maxmax_likelihood = max_likelihood
                maxmax_mixture = max_mixture
                maxmax_membership = max_membership


        #print ("NUM_CLONE = {0} : {1}\n{2}\n\n".format(NUM_CLONE,maxmax_likelihood, maxmax_mixture))
        #print (likelihood_record)
        
        maxmaxmax_likelihood_record.append (maxmax_likelihood)
        maxmaxmax_mixture_record.append (maxmax_mixture)
        maxmaxmax_membership_record.append (maxmax_membership)

        if method.split("+")[0] == "max":
            if maxmaxmax_likelihood < maxmax_likelihood :   # 최고의 값만 구하고 싶을 때
                previous_maxmax_likelihood = maxmax_likelihood
                maxmaxmax_likelihood = maxmax_likelihood
                maxmaxmax_mixture = maxmax_mixture
                maxmaxmax_membership = maxmax_membership
                maxmaxmax_NUM_CLONE = NUM_CLONE


    #print ("X : ", np.arange(start = NUM_CLONE_TRIAL_START , stop = NUM_CLONE_TRIAL_END + 1 ))
    #print ("Likelihood record : " , maxmaxmax_likelihood_record)


    x= np.arange(start = NUM_CLONE_TRIAL_START , stop = NUM_CLONE_TRIAL_END + 1 )
    y = maxmaxmax_likelihood_record
    if NUM_CLONE_TRIAL_END > NUM_CLONE_TRIAL_START:
        for S in list(np.arange(0.1, 3.0, 0.5)):
            kneedle = KneeLocator(x, y, S=S, curve="concave", direction="increasing")
            if kneedle.knee != None:
                kneedle.plot_knee()
                break
    
    print (list(x), y)

    if method.split("+")[0] == "elbow":
        if NUM_CLONE_TRIAL_START == NUM_CLONE_TRIAL_END:       # 1개밖에 없으면  kneedle이 오류난다
            threshold_x = NUM_CLONE_TRIAL_END
        else:
            if kneedle.knee != None:             # 가장 정상적인 경우
                threshold_x = round(kneedle.knee, 1)
            else:                                             # elbow 법에서 제대로 못 찾았을 경우  단순히 가장 높은 값으로 한다
                threshold_x = np.argmax(maxmaxmax_likelihood_record) + NUM_CLONE_TRIAL_START
        
        maxmaxmax_NUM_CLONE = int (threshold_x)
        maxmaxmax_likelihood = maxmaxmax_likelihood_record[int (threshold_x) - NUM_CLONE_TRIAL_START]
        maxmaxmax_mixture = maxmaxmax_mixture_record[int (threshold_x) - NUM_CLONE_TRIAL_START]
        maxmaxmax_membership = maxmaxmax_membership_record[int (threshold_x) - NUM_CLONE_TRIAL_START]





    #print (method.split("+"))

    if method.split("+")[1] == "outlier":
        fp_candidate_index_total = []
        for j in range(maxmaxmax_NUM_CLONE):       # 각 cluster마다 돌면서 outlier 찾아보기
            fp_candidate_index = Estep_final ( df, maxmaxmax_mixture, maxmaxmax_membership, j )            # outlier 구하러 떠나가자
            fp_candidate_index_total = fp_candidate_index_total  + fp_candidate_index                                   # 각 clone의 outlier index를 모아주기
            #print ("cluster {0} : {1}".format(j, fp_candidate_index))
        #print ("\n")
    
        for i in fp_candidate_index_total:
            maxmaxmax_membership[i] = maxmaxmax_NUM_CLONE
        maxmaxmax_NUM_CLONE = maxmaxmax_NUM_CLONE + 1         # fp cluster가 생기는 것이니까 + 1

        # 새로 mixture를 정해준다
        maxmaxmax_mixture = Mstep(maxmaxmax_membership , df, maxmaxmax_NUM_CLONE, NUM_BLOCK, adjustment, "Last") 

        #print ("Outlier removal")
        #print ("최종 정답\n\tclone={0}\n{1}".format(maxmaxmax_NUM_CLONE, maxmaxmax_likelihood))
        return maxmaxmax_NUM_CLONE, maxmaxmax_membership, maxmaxmax_mixture    

    else:
        #print ("최종 정답\n\tclone={0}\n{1}".format(maxmaxmax_NUM_CLONE, maxmaxmax_likelihood))
        return maxmaxmax_NUM_CLONE, maxmaxmax_membership, maxmaxmax_mixture
