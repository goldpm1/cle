import  filetype, argparse, os, scipy, math
import numpy as np
import pandas as pd
import random
import copy
import collections
from random import *


# python3  1.SimData_pipe0_preparation.py --NUM_CLONE 4 --NUM_BLOCK 2 --NUM_MUTATION 500 --FP_RATIO 0.1
def TN_prior_cal(x):
    from scipy.special import expit
    return (1 - expit( 200*x - 5)) * 0.15

def fp_sampling ( **kwargs ):
    import random 

    y_acc, x_acc = [], np.zeros (100, dtype = float)

    x = np.linspace(0, 1.6, 101)

    t = 0
    for i, j in zip(x , TN_prior_cal(x)):
        t = t + j
        y_acc.append (t)

    for y_index, y in enumerate ( np.random.uniform(0, np.max(y_acc), 100) ) :
        if y < y_acc [0]:
            x_acc [y_index]  = 0
        for k in range ( len(x_acc) - 1 ):
            if ( y_acc[k] < y ) &  ( y < y_acc[k + 1] ) :
                x_acc [y_index] = k / 100
                break

    fp_pool =  sorted (x_acc) 

    return (random.sample (fp_pool, kwargs["NUM_BLOCK"]))

###################################################################################################################################################################################################

def dirichlet_sampling ( trial, **kwargs ):
    import random, math
    import numpy as np

    random.seed( trial )
    np.random.seed( trial )

    df = [[None] * NUM_BLOCK for i in range(NUM_MUTATION)]
    np_vaf = np.zeros((NUM_MUTATION, NUM_BLOCK), dtype = 'float')
    mutation_id = np.zeros((NUM_MUTATION), dtype = 'int')
    membership = np.zeros (NUM_MUTATION, dtype = "str")
    mixture = np.zeros ((NUM_BLOCK, NUM_CLONE), dtype = 'float') 

    for row in range (NUM_MUTATION):
        for col in range (NUM_BLOCK):
            df[row][col] = {"depth":0, "ref":0, "alt":0}


    # Dirichlet sampling : Clone별 개수 정하고, VAF 정하고 , Depth/Alt 정해주기
    print ("NUM_CLONE : {}".format (NUM_CLONE))
    print ("FP_RATIO : {}\n".format (kwargs["FP_RATIO"], type(kwargs["FP_RATIO"])))

   
    NUM_CLONE_ACTIVE = NUM_CLONE


    #0. 각 cluster마다 membership을 몇 개 줄지 정함   (clone 0은 무조건 low vaf로 주자)   :  equal distribution
    li = list(range(0 , NUM_CLONE_ACTIVE))   # 분율이 0인 경우도 포함
    membership_count = collections.Counter([random.choice(li) for i in range( int (NUM_MUTATION * ( 1 - kwargs["FP_RATIO"] )) )])        # membership 개수 정해주기
    membership_count = collections.OrderedDict ( sorted (membership_count.items()) )                        # key 순서대로 sort 해주기
    print ( "►membership_count = {}".format (membership_count ))
    print ( "\t→ space = {}\tFP = {}\n".format (np.sum ( [values for key, values in membership_count.items()] ), NUM_MUTATION - np.sum ( [values for key, values in membership_count.items()] )))


    

    for i in range (NUM_BLOCK):
        print ("\n#######################################   {}TH BLOCK  #########################################".format (i))

        #1. 일단 Dirichlet sampling에 들어갈 alpha 계수를 정해준다
        li = np.zeros ( NUM_CLONE_ACTIVE, dtype = "int")         # 기준이 되는 지점
        li_2 = np.zeros ( NUM_CLONE_ACTIVE, dtype = "int")      # 진짜 계수

        if kwargs["SimData"] == "decoy":
            num_zero = int ( np.random.uniform (0, math.ceil (NUM_CLONE_ACTIVE / 2) ) ) if NUM_BLOCK >= 2 else 0      # 맘에 안들면 무조건 0으로 하면 됨
            if NUM_BLOCK == 1:         # 1차원에서는 0을 빼주자
                num_zero = 0
            num_zero = 0
            #li = np.array( [0] * num_zero +  list ( np.round ( np.random.uniform (1, 100, NUM_CLONE_ACTIVE - num_zero) ) ) )    # 0 ~100 중에 j - num_zero 개를 뽑음
            li = np.array(  list ( np.round ( np.random.uniform (1, 100, NUM_CLONE_ACTIVE - num_zero) ) ) )    # 0 ~100 중에 j - num_zero 개를 뽑음
        elif kwargs["SimData"] == "lump":
            num_zero = int ( np.random.uniform (0, math.ceil (NUM_CLONE_ACTIVE / 3) ) ) if NUM_BLOCK >= 2 else 0      # 맘에 안들면 무조건 0으로 하면 됨. 좀더 확률을 낮추자
            if NUM_BLOCK == 1:         # 1차원에서는 0을 빼주자
                num_zero = 0
            for j in range (num_zero, NUM_CLONE_ACTIVE):  
                li[j] = (NUM_CLONE) * 10 - (j * 6)


        li = np.round (li)
        random.shuffle (li)       # 섞어주자
        for j in range (NUM_CLONE_ACTIVE):      # 기준값(li)에서 binomial로 난수 추출
            if kwargs["SimData"] == "decoy":
                if li[j] < 25:
                    p = 0.3
                elif li[j] < 40:
                    p = 0.4
                elif li[j]  < 70:
                    p = 0.5
                else:
                    p = 0.6
            elif kwargs["SimData"] == "lump":
                p = 0.5

            li_2[j] = np.random.binomial  ( li[j], 0.5 )
        li_2 = np.array (li_2)
        print ("\tli = {}\n\tli_2 = {}\n".format (li, li_2))
        
                
        #2. Dirichlet sampling 수행
        s = np.zeros ( (int (NUM_MUTATION * ( 1 - kwargs["FP_RATIO"] )),  NUM_CLONE_ACTIVE), dtype = "float" )
        ss = np.random.dirichlet ( [k for k in li_2 if k != 0] , int (NUM_MUTATION * ( 1 - kwargs["FP_RATIO"] )) )    # li_2가 0인 것은 뺴고 수행
        jj = 0
        for j in range( s.shape[1] ):       # 0때문에 발생한 k * 3  -> k * 4로 늘려주기
            if li_2[j] != 0:
                s [:, j] = ss [:, jj]
                jj += 1
        #print ( s, s.shape  )
        

        #3. membership 개수만큼 np_vaf에 집어넣어주기
        k_now = 0
        for j, j_count in membership_count.items(): 
            if li_2[j] == 0:
                for k in range (k_now, k_now + j_count ):
                    np_vaf [k][i] = 0
                    membership[k] = j
                    mutation_id[k] = k                    
                k_now += j_count
            else:  # 정상인 경우
                for k in range (k_now, k_now + j_count ):
                    np_vaf [k][i] = s [k][j] / 2        
                    membership[k] = j
                    mutation_id[k] = k
                k_now = k_now + j_count
        


        # 결과 출력
        print ( "\t► alpha coefficient (Dirichlet distribution) : {}".format(li_2) )
        ttt = 0
        for j in range ( s.shape[1] ):
            mixture[i][j] = np.mean ( s [: membership_count[j], j] )
            print ("\t\tclone {} (0 ~ {})의 mean : {}".format ( j, membership_count[j], round( mixture[i][j], 2)  ) )
            ttt = ttt + np.mean ( s [: membership_count[j], j] )
        print ("\t\ttotal mean : {}\n".format ( round ( ttt, 2)  ))


        #4.  Depth/Alt 정해주기
        for k in range (int (NUM_MUTATION * ( 1 - kwargs["FP_RATIO"] ))):
            while True:
                df[k][i]["depth"] = int(np.random.normal( kwargs["DEPTH_MEAN"], kwargs["DEPTH_SD"], size = 1))
                if df[k][i]["depth"] > kwargs["DEPTH_CUTOFF"]:
                    break
            df[k][i]["alt"] =  round( df[k][i]["depth"] * np_vaf [k][i] )
            df[k][i]["ref"] = df[k][i]["depth"]  -  df[k][i]["alt"] 
            df[k][i]["membership_answer"] = membership[k]


    # FP sampling
    k = int (NUM_MUTATION * ( 1 - kwargs["FP_RATIO"] ))  # 맨 뒤에부터 하나씩 추가
    while k < NUM_MUTATION:
        pp = fp_sampling (**kwargs)

        if np.sum (pp) != 0:
            np_vaf [k] = pp
            membership[k] = "FP"
            
            for i in range (NUM_BLOCK):
                while True:
                    df[k][i]["depth"] = int(np.random.normal( kwargs["DEPTH_MEAN"], kwargs["DEPTH_SD"], size = 1))
                    if df[k][i]["depth"] > kwargs["DEPTH_CUTOFF"]:
                        break
                df[k][i]["alt"] =  round( df[k][i]["depth"] * np_vaf [k][i] )
                df[k][i]["ref"] = df[k][i]["depth"]  -  df[k][i]["alt"] 
                df[k][i]["membership_answer"] = "FP"

            k = k + 1


    return df, np_vaf, mixture, membership, mutation_id
    

#####################################################################################################################

def printresult():
    global NUM_CLONE, NUM_BLOCK, NUM_MUTATION, mixture, membership, df, np_vaf, mutation_id

    os.system ("rm -rf " + kwargs ["INPUT_TSV"]  )
    os.system ("rm -rf " + kwargs ["NPVAF_DIR"] + "/npvaf.txt" )
    
    output_file_inputtsv = open( kwargs ["INPUT_TSV"], "w" )
    output_file_npvaf = open( kwargs ["NPVAF_DIR"] + "/npvaf.txt" , "w" )

    # print 1st row
    print ("", end = "\t", file = output_file_npvaf)
    for i in range(NUM_BLOCK):
        print ("block{0}".format(i), end = "\t", file = output_file_npvaf)
    print ("membership_answer", file = output_file_npvaf)

    # print 2 ~ k+1th row
    for k in range(NUM_MUTATION):
        print ("mut_{}".format(k), end = "\t", file = output_file_npvaf)
        for i in range (0, NUM_BLOCK):
            print ( np_vaf[k][i] , end = "\t", file = output_file_npvaf)
            #print (df[k][i]["depth"], ".", df[k][i]["ref"], ".", df[k][i]["alt"], sep = "", end = "\t", file = output_file_npvaf)
        #print ( "{} {}".format(df[k][i]["membership_answer"], membership[k]), file = output_file_npvaf)
        print ( "{}".format(df[k][i]["membership_answer"]), file = output_file_npvaf)
    output_file_npvaf.close()

    # print df (INPUT_TSV)
    for k in range(0, NUM_MUTATION):
        print ("mut_{}\t{}\t".format( k, df[k][i]["membership_answer"] ) , end = "", file = output_file_inputtsv)
        for i in range (NUM_BLOCK - 1):
            print ( "{},{}".format( df[k][i]["depth"], df[k][i]["alt"] ), end = ",", file = output_file_inputtsv)
        print ( "{},{}".format( df[k][NUM_BLOCK - 1]["depth"], df[k][NUM_BLOCK - 1]["alt"] ), file = output_file_inputtsv)
    output_file_inputtsv.close()



#####################################################################################################################


if __name__ == "__main__":
    kwargs = {}

    parser = argparse.ArgumentParser(description='The below is usage direction.')
    parser.add_argument('--NUM_CLONE', type = int, default=4)
    parser.add_argument('--NUM_BLOCK', type = int, default=3)
    parser.add_argument('--NUM_MUTATION', type = int, default = 500)
    parser.add_argument('--FP_RATIO', type = float, default = 0)
    parser.add_argument('--LOWVAF_RATIO', type = float, default = 0)
    parser.add_argument('--DEPTH_MEAN', type = int, default = 100)
    parser.add_argument('--DEPTH_SD', type = int, default = 8)
    parser.add_argument('--DEPTH_CUTOFF', type = int, default = 30)
    parser.add_argument('--INPUT_TSV', default = None)
    parser.add_argument('--NPVAF_DIR', default = None)
    parser.add_argument('--BENCHMARK_I', default = 0)
    parser.add_argument('--SimData' , default = "decoy")       # "decoy or sparse"
    
    args = parser.parse_args()

    global NUM_CLONE, NUM_BLOCK, NUM_MUTATION, mixture, membership, df, np_vaf, mutation_id

    NUM_CLONE, NUM_BLOCK, NUM_MUTATION = args.NUM_CLONE, args.NUM_BLOCK, args.NUM_MUTATION
    kwargs["NUM_CLONE"], kwargs["NUM_BLOCK"], kwargs["NUM_MUTATION"] = args.NUM_CLONE, args.NUM_BLOCK, args.NUM_MUTATION
    kwargs ["FP_RATIO"] = float (args.FP_RATIO)
    kwargs ["LOWVAF_RATIO"] = float (args.LOWVAF_RATIO)
    kwargs["DEPTH_MEAN"] = int(args.DEPTH_MEAN)
    kwargs["DEPTH_SD"] = float(args.DEPTH_SD)
    kwargs["DEPTH_CUTOFF"] = float(args.DEPTH_CUTOFF)
    kwargs["NPVAF_DIR"] = args.NPVAF_DIR
    kwargs["BENCHMARK_I"] = int(args.BENCHMARK_I)
    ii = kwargs["BENCHMARK_I"]
    kwargs["SimData"] = args.SimData
    kwargs["INPUT_TSV"] = args.INPUT_TSV
    kwargs["NPVAF_DIR"] = args.NPVAF_DIR
    kwargs["RANDOM_SEED"] = int(args.BENCHMARK_I)
    

    # 원점이 centroid로 들어가있으면 안되는데..
    trial = kwargs["BENCHMARK_I"]
    while True:
        df, np_vaf, mixture, membership, mutation_id = dirichlet_sampling (trial, **kwargs )
        if np.any(np.all(mixture == 0, axis=0)) == False:      # 원점은 없어야지 허락해주고 끝내준다
            break
        else:
            trial += 100
    

    # 파일 출력
    printresult()

