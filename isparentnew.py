import scipy.stats
import numpy as np
import itertools, comb, math
from sklearn.decomposition import TruncatedSVD, PCA
from matplotlib import pyplot as plt

def greaterall (a, b, boole):
    if boole == "<":
        for i in range(len(a)):
            if a[i] > b[i]:
                return False
        return True
    if boole == ">":
        for i in range(len(a)):
            if a[i] < b[i]:
                return False
        return True

def checkall (sum_mixture):
    for i in range ( len(sum_mixture) ):
        if ( sum_mixture[i] < 0.77 ) | ( sum_mixture[i] > 1.25 ):       # 한 block이라도 1에 걸맞지 않은 게 있으면 False를 return함
            return False
    return True


def makeone (step, **kwargs):
    membership = step.membership
    mixture = step.mixture

    global subset_list_acc, subset_mixture_acc, sum_mixture_acc

    
    if step.fp_index != -1:       # Outlier가 있으면 그건 빼자
        subset_list_acc, subset_mixture_acc, sum_mixture_acc = comb.comball(list(set(membership) - set( [step.fp_index ])), mixture)   # 모든 덧셈 조합을 구하기
    elif step.fp_index == -1:
        subset_list_acc, subset_mixture_acc, sum_mixture_acc = comb.comball(list(set(membership))[: ], mixture)   # 모든 덧셈 조합을 구하기


    if kwargs["NUM_CLONE_NOMINAL"] == 1:   # 이상하게 NUM_CLONE == 1일 때에는 comb.comball이 작동하지 않는다
        return [0], [[1, 0]], -1          # makeone_index , p_list, fp_index

    
    p_max = float("-inf")
    p_list, j2_list = [], []
    fp_possible = []
    fp_possible_clone = np.zeros (len(subset_mixture_acc), dtype = "int")

    for j2 in range(len(subset_mixture_acc)):      #여러 조합을 돈다   # 맨 뒤는 outlier ( = inside membership + outside membership이니까 제외)
        subset_list, subset_mixture, sum_mixture = subset_list_acc[j2], subset_mixture_acc[j2], sum_mixture_acc[j2]

        if checkall(sum_mixture) == False:         # 한 block이라도 1에 걸맞지 않은 게 있으면 False를 return함
            continue

        p = 0
        for i in range (kwargs["NUM_BLOCK"]):
            depth = 1000
            a = int(sum_mixture[i] * 1000 / 2) 
            b = depth - a

            target_a = 500
            try:
                p = p + math.log10(scipy.stats.betabinom.pmf(target_a, depth, a + 1, b+1))
            except:
                p = p - 400
        if p > -400:
            ISSMALLER_cnt, SMALLER_ind = 0, []

            for j3 in (set(range(0, mixture.shape[1]  )) - set(subset_list)) :    # 나머지 clone
                if (j3 == step.fp_index) | ( len (set (membership)) <=  j3 )  :
                    continue

                ISLARGER_cnt = 0
                for j4 in subset_list :    #선정된 boundary 후보          나머지 clone중에 child clone보다 작으면 안됨.  그러나 딱 1개만 있고 FP clone이면 용서해준다
                    if greaterall ( mixture[:,j3] , mixture[:,j4], "<" ) == True:
                        ISLARGER_cnt = ISLARGER_cnt + 1
                        #break
                if ISLARGER_cnt >= 1:       # 적어도 1개의 boundry 후보 (child 후보)들 보다 작다는게 말이 안됨.  다만 FP clone이 있을 경우에는 다르다
                    ISSMALLER_cnt = ISSMALLER_cnt + 1
                    SMALLER_ind.append (j3)


            if (step.includeoutlier == False) &  (ISSMALLER_cnt == 1):       # 그동안 FP clone 없었는데 유일하게 안쪽에서 발견할 경우
                
                #print ("0 :  {} 번째 clone은 {}보다 안쪽에 있는 유일한 clone이라 살리는걸 고려.  개수는 {}개. ".format(SMALLER_ind[0], subset_list,   int(np.unique(membership, return_counts=True)[1] [SMALLER_ind[0]] )))
                check = 0

                for j3 in set(range(0, mixture.shape[1] )) - set(subset_list) - set([SMALLER_ind[0]]):   # 나머지 clone (outlier, FP도 빼고) 을 돈다
                    if j3 == step.fp_index:
                        continue

                    tt = []
                    for j4 in subset_list:   # boundary clone (putative child clone)을 돈다
                        if greaterall ( mixture[:,j3] , mixture[:,j4], ">" ) == True:
                            tt.append (mixture[:,j4])

                    if len(tt) < 2:        # 나머지 clone은 2개 이상의 child 조합으로 만들어지는 parent여야 한다. 그게 만족 안하면 이 조합이 오류임을 알 수 있다
                        #print ("{} ( {} ) 는 2개 이상의 child clone의 합으로 만들어지지지 않아서 이 조합 (FP 포함, {}) 은 아예 기각".format(j3 , mixture[:, j3], subset_mixture))
                        check = 1
                        break

                if check == 0:
                    if SMALLER_ind[0] >= np.max(membership) + 1:
                        print ("mixture : {}".format(list(mixture)))
                        print ("SMALLER_ind : {}".format(SMALLER_ind))
                        print ("membership : {}".format(np.unique(membership, return_counts=True)[1]))

                    # if len(subset_list) == 3:
                    #         print ("조합 : {},   FP (  {} ( {} )  {}  )  :   {}  MIN_CLUSTER_SIZE : {}".format(subset_list, j3 , mixture[:, j3], SMALLER_ind[0], int(np.unique(membership, return_counts=True)[1] [SMALLER_ind[0]] ) , kwargs["MIN_CLUSTER_SIZE"]  )   )

                    try:
                        if ( int(np.unique(membership, return_counts=True)[1] [SMALLER_ind[0]] )  > kwargs["MIN_CLUSTER_SIZE"]  ):
                            #if ( int(np.unique(membership, return_counts=True)[1] [SMALLER_ind[0]] )  < 30  ):
                            fp_possible.append (j2)
                            fp_possible_clone[j2] =  SMALLER_ind[0]
                            p_list.append ( [ p, j2] )
                    except:
                        print ("뭔가 문제가 있다")
                        print ("FP로 생각되는 clone number : {}".format(SMALLER_ind[0]))
                        print ( int(np.unique(membership, return_counts=True) ) )



            if ISSMALLER_cnt == 0:      # FP clone 없이 순전하게 child clone이 구성됐을 경우
                check = 0
                for j3 in set(range(0, mixture.shape[1] )) - set(subset_list):   # 나머지 clone (outlier를 빼고) 을 돈다
                    if j3 == step.fp_index:
                        continue
                    tt = []
                    for j4 in subset_list:   # boundary clone (putative child clone)을 돈다
                        if greaterall ( mixture[:,j3] , mixture[:,j4], ">" ) == True:
                            tt.append (mixture[:,j4])

                    if len(tt) < 2:        # 나머지 clone은 2개 이상의 child 조합으로 만들어지는 parent여야 한다. 그게 만족 안하면 이 조합이 오류임을 알 수 있다
                        #print ("{} ( {} ) 는 2개 이상의 child clone의 합으로 만들어지지지 않아서 이 조합은 (FP 미포함, {})은 아예 기각".format(j3 , mixture[:, j3], subset_mixture ))
                        check = 1
                        break
                if check == 0:
                    p_list.append ( [ p, j2] )        # boundary clone (child clone) 조합을 넣어주자.  그래도 최고 후보만 아니라 몇 개는 보자



    if p_list == []:
        return [], [], -1

    p_list = np.array(p_list).transpose()
    p_list = p_list[ :, p_list[0].argsort()]
    p_list = np.flip(p_list , 1)

    if kwargs["VERBOSE"] >= 2:
        if step.fp_index != -1:
            for i in range(len(p_list[1])):
                print ("step={}\tclone = {}, sum = {} : p = {}".format (kwargs["STEP_TOTAL"], subset_list_acc[ int(p_list[1][i]) ] ,  sum_mixture_acc[ int(p_list[1][i]) ],  p_list[0][i]))
            print ("")

    best_j2 = int( p_list[1,0] )

    if best_j2 in fp_possible:
        if kwargs["VERBOSE"] >= 2:
            print ("\t\t{} 번째 clone은 {}보다 안쪽에 있는 유일한 clone이라 살려줌. 그런게 그게 ({}) 가장 확률 높은 makeone이 되버렸네 ".format(fp_possible_clone[best_j2], subset_list_acc[best_j2], subset_list_acc[best_j2]))
        return subset_list_acc[best_j2], p_list, fp_possible_clone[best_j2]

        # try:
        #     print ("NUM_CLONE : {}\t{}-{}\tmakeone_index\t{}\tsum\t{}".format(mixture.shape[1], kwargs["TRIAL"], kwargs["STEP"], subset_list_acc[best_j2], sum_mixture_acc[best_j2] ) )
        # except:
        #     print ("NUM_CLONE : {}\tmakeone_index\t{}\tsum\t{}".format(mixture.shape[1], subset_list_acc[best_j2], sum_mixture_acc[best_j2] ) )

    return subset_list_acc[best_j2], p_list, step.fp_index

