import numpy as np
import scipy
import copy
from scipy.special import expit
import random
import math

def phred_to_percentile (phred_score):
    return round ( 10 ** (-phred_score / 10) , 5)  #  Converts a Phred score to a percentile.


def expected_calculator (  i, j, k, mixture, df, input_containpos, **kwargs ):
    import re

    if (kwargs["SEX"] == "M") & ( bool(re.search(r'X|Y', input_containpos.iloc[k]["pos"]))  == True  ) :
        depth_calc, alt_calc = int(df[k][i]["depth"] ), int(df[k][i]["depth"] * mixture[i][j])
    else:
        depth_calc, alt_calc = int(df[k][i]["depth"] ), int(df[k][i]["depth"] * mixture[i][j] * 0.5)

    depth_obs, alt_obs = int(df[k][i]["depth"]), int(df[k][i]["alt"])
    a, b = alt_calc, depth_obs - alt_calc              # alt_expected, ref_expected

    return (depth_calc, alt_calc, depth_obs, alt_obs, a, b)


def calc_likelihood(input_containpos, df,  np_vaf, np_BQ, step, k, **kwargs):
    import re
    mixture = step.mixture

    max_prob = float("-inf")
    max_clone = -1

    prob = np.zeros(kwargs["NUM_CLONE"], dtype="float64")            # 합치면 1 되게
    prob_abs = np.zeros(kwargs["NUM_CLONE"], dtype="float64")    # 순수한 beta binomal만..

    check = 0


    if kwargs["DEBUG"] == True:
        #debug_k = np.where(  ( (np_vaf[:, 0] > 0.3 ) &  (np_vaf[:, 0] < 0.5 ) & ( np_vaf[:, 1] < 0.3 ) )  )[0]
        debug_k = [1]
        if (k in debug_k):
            print ("\t\t\tk = {}\tNUM_CLONE = {}, NUM_BLOCK = {}, df = [{},{}]".format(k, kwargs["NUM_CLONE"], kwargs["NUM_BLOCK"], len(df), len(df[0]), mixture))
    else:
        debug_k = []


    for j in range(kwargs["NUM_CLONE"]): 
        for i in range(kwargs["NUM_BLOCK"]):
            alt_obs = int(df[k][i]["alt"])

            if alt_obs == 0:    # TN or FN?   (합치면 1이 되어야 한다))
                if mixture [i][j] == 0 :  # TN
                    p = p_abs =  math.log10 ( kwargs ["TN_CONFIDENTIALITY"] )            # nearly 0

                else: # FN
                    p1 = math.log10 ( 1 - kwargs ["TN_CONFIDENTIALITY"] )

                    depth_calc, alt_calc, depth_obs, alt_obs, a, b = expected_calculator (  i, j, k, mixture, df, input_containpos, **kwargs )
                    try:
                        p2_numerator = scipy.stats.betabinom.pmf(alt_obs, depth_obs, a+1, b+1)   # 분자
                        #p2_numerator = scipy.spatial.distance.euclidean ( np_vaf[k] * 2 ,  np.array ( mixture [: , j] ) )
                        
                        p2_denominator = 0   # 분모
                        for whole_j in range(kwargs["NUM_CLONE"]): 
                            if (whole_j != step.fp_index):
                                depth_calc, alt_calc, depth_obs, alt_obs, a, b = expected_calculator (  i, whole_j, k, mixture, df, input_containpos, **kwargs )
                                try:
                                    p2_denominator += scipy.stats.betabinom.pmf( alt_obs, depth_obs, a+1, b+1 )
                                    #p2_denominator += scipy.spatial.distance.euclidean(  np_vaf[k] * 2  ,  np.array ( mixture [:, whole_j] ) )
                                except:
                                    p2_denominator += 0
                        p2 = p2_numerator / p2_denominator
                        try:
                            p_abs = math.log10 (p2_numerator)
                        except:
                            p_abs = -400

                    except:
                        p2 = p_abs = float("-inf")
                
                    try:
                        p = p1 + math.log10 ( p2 )
                    except:
                        p = -400

                # if(kwargs ["DEBUG"] == True ) :
                #     if  ( k in debug_k)  :            
                #         np.set_printoptions(suppress=True)   # Scientific expression이 싫어요
                #         print ( "\t\t\t\t\tj = {}\ti = {}\tp = {}".format( j, i,  round(p, 2 )) )


            else:   # TP or FP?
                SEQ_ERROR = phred_to_percentile ( np_BQ[k][i] )
                depth_calc, alt_calc, depth_obs, alt_obs, a, b = expected_calculator (  i, j, k, mixture, df, input_containpos, **kwargs )

                if mixture [i][j] == 0: # FP
                    try:
                        p = p_abs = math.log10(scipy.stats.binom.pmf(n = depth_obs, p = SEQ_ERROR, k = alt_obs))
                    except:
                        p = p_abs = -400
                else:  # TP
                    try:
                        p1 = math.log10(1 - scipy.stats.binom.pmf(n = depth_obs, p = SEQ_ERROR, k = alt_obs))   # Not FP
                    except:
                        p1 = -400

                    try:
                        p2_numerator = scipy.stats.betabinom.pmf(alt_obs, depth_obs, a+1, b+1)   # 분자
                        #p2_numerator = scipy.spatial.distance.euclidean ( np_vaf[k] * 2 ,  np.array ( mixture [: , j] ) )
                        p2_denominator = 0   # 분모
                        for whole_j in range(kwargs["NUM_CLONE"]): 
                            if (whole_j != step.fp_index):
                                depth_calc, alt_calc, depth_obs, alt_obs, a, b = expected_calculator (  i, whole_j, k, mixture, df, input_containpos, **kwargs )
                                try:
                                    p2_denominator += scipy.stats.betabinom.pmf( alt_obs, depth_obs, a+1, b+1 )
                                    #p2_denominator += scipy.spatial.distance.euclidean(  np_vaf[k] * 2  ,  np.array ( mixture [:, whole_j] ) )
                                except:
                                    p2_denominator += 0
                        
                        p2 = p2_numerator / p2_denominator
                        try:
                            p_abs = math.log10 (p2_numerator)
                        except:
                            p_abs = -400

                        
                    except:
                        p2 = p_abs = float ("-inf")
                    

                    try:
                        p = p1 + math.log10 ( p2 )
                    except:
                        p = p1 - 400


                # if(kwargs ["DEBUG"] == True ) :
                #     if  ( k in debug_k)  :            
                #         np.set_printoptions(suppress=True)   # Scientific expression이 싫어요
                #         print ( "\t\t\t\t\tj = {}\ti = {}\tp = {}".format( j, i,  round(p, 2 )) )


            prob[j] += p
            prob_abs[j] += p_abs

        if prob[j] > max_prob:
            max_prob = prob[j]
            max_prob_clone_candidate = [j]
        elif prob[j] == max_prob:
            max_prob_clone_candidate.append(j)
    


    max_clone = random.choice(max_prob_clone_candidate)

    # if(kwargs ["DEBUG"] == True ) :
    #     if  ( k in debug_k)  :
    #         print ("prob = {}".format (prob))
    #         np.set_printoptions(suppress=False)   
    #         print ( "\t\t\t\t\tk = {} ({})\tmax_clone = {}\tprob = {}\tprob_abs = {}\tmax_prob_abs = {}".format( k,  np.round( np_vaf[k] * 2, 2 ) , max_clone, np.round ( np.power (10, prob) , 3) , np.round ( np.power (10, prob_abs) , 3) , round ( prob_abs[max_clone], 3)  ))
    #         np.set_printoptions(suppress=True)   # Scientific expression이 싫어요
    #         print ( "\t\t\t\t\tk = {} ({})\tmax_clone = {}\tprob = {}\tprob_abs = {}\tmax_prob_abs = {}".format( k,  np.round( np_vaf[k] * 2, 2 ) , max_clone, np.round ( np.power (10, prob) , 3) , np.round ( np.power (10, prob_abs) , 3) , round ( prob_abs[max_clone], 3)  ))


    if kwargs["OPTION"] in ["Hard", "hard"]:
        return list(prob), max_prob, prob_abs [max_clone], max_clone

    elif kwargs["OPTION"] in ["Soft", "soft"]:
        weight = np.power (10, prob_abs)   # prob? prob_abs?
        if  ( k in debug_k)  :
            np.set_printoptions(suppress=True)   # Scientific expression이 싫어요
            print (weight )
        new_likelihood = round(np.average(prob_abs, weights = weight), 3)       # Likelihood in Soft clustering


        return list(prob), new_likelihood, new_likelihood, max_clone 



def main (input_containpos, df, np_vaf, np_BQ, step, **kwargs):
    total_prob = total_prob_abs = 0

    # if kwargs["DEBUG"] == True:
    #     debug_k = np.where(  ( (np_vaf[:, 1]  == 0 ) &  (np_vaf[:, 0] < 0.03 )) )  [0]
        #print (debug_k)

    max_prob_abs_list = []
    for k in range(kwargs["NUM_MUTATION"]):
        step.membership_p[k], max_prob, max_prob_abs, step.membership[k] = calc_likelihood(input_containpos, df,  np_vaf, np_BQ, step, k, **kwargs)
        total_prob += max_prob
        total_prob_abs += max_prob_abs
        max_prob_abs_list.append (max_prob_abs)

    max_prob_abs_list = np.array (max_prob_abs_list )
    
    # if (kwargs["DEBUG"] == True):
    #     #print ("total_prob = {}\ttotal_prob_abs = {}".format( round( total_prob, 3) , round(total_prob_abs, 3) ))
    #     for j in range ( kwargs["NUM_CLONE"]):
    #         if j in step.membership:
    #             print ( "\t\t\t\t\tCLONE = {}\tsum (max_prob_abs_list) = {}".format (j , round ( np.sum (max_prob_abs_list [ np.where (step.membership == j)[0] ]) , 1) , np.round ( np.mean (np_vaf [ np.where (step.membership == j)[0]  ] * 2,  axis = 0 ), 2)  ) )

    step.likelihood = total_prob_abs
    step.likelihood_record[kwargs["STEP"]] = total_prob_abs


    if  step.fp_index in set(step.membership):   # FP를 선택한 variant가 하나라도 있어야 비로소 includefp = True로 변경
        step.includefp = True
        step.fp_member_index = list(np.where(step.membership == step.fp_index)[0])
    else:
        step.includefp = False
        step.fp_member_index = []


    if (kwargs["VERBOSE"] >= 1):
        print ("\t\t\tEstep.py : set(step.membership) = {}\tcounts = {}\tfp_index = {}\tincludefp = {}\tstep.likelihood = {}".format ( set(step.membership), np.unique(step.membership  , return_counts=True)[1], step.fp_index, step.includefp, round(step.likelihood)) )


    # membership_p : extermely low value if the variant is  fp
    for k in range(kwargs["NUM_MUTATION"]):
        if step.membership[k] == step.fp_index:
            step.membership_p[k] = [-999] * len(step.membership_p[k])

    return step
