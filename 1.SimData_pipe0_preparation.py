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
    return (1 - expit( 100*x - 5)) * 0.5

def fp_sampling ( **kwargs ):
    import random 

    y_acc, x_acc = [], np.zeros (100, dtype = float)

    x = np.linspace(0, 1, 101)

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

def dirichlet_sampling ( **kwargs ):
    import random
    df = [[None] * NUM_BLOCK for i in range(NUM_MUTATION)]
    np_vaf = np.zeros((NUM_MUTATION, NUM_BLOCK), dtype = 'float')
    mutation_id = np.zeros((NUM_MUTATION), dtype = 'int')
    membership = np.zeros (NUM_MUTATION, dtype = "str")
    mixture = np.zeros ((NUM_BLOCK, NUM_CLONE), dtype = 'float') 

    for row in range (NUM_MUTATION):
        for col in range (NUM_BLOCK):
            df[row][col] = {"depth":0, "ref":0, "alt":0}


    # Dirichlet sampling : Clone별 개수 정하고, VAF 정하고 , Depth/Alt 정해주기

    if kwargs["FP_RATIO"]  != 0:
        NUM_CLONE_ACTIVE = NUM_CLONE
    else:
        NUM_CLONE_ACTIVE = NUM_CLONE + 1

    #0. 각 cluster마다 membership을 몇 개 줄지 정함   (clone 0은 무조건 low vaf로 주자)
    if (kwargs ["LOWVAF_RATIO"] != 0):
        li = list(range(1 , NUM_CLONE_ACTIVE))   # 분율이 0인 경우도 포함
    else:
        li = list(range(0 , NUM_CLONE_ACTIVE))   # 분율이 0인 경우도 포함

    membership_count = collections.Counter([random.choice(li) for i in range( int (NUM_MUTATION * ( 1 - kwargs["FP_RATIO"] - kwargs["LOWVAF_RATIO"])) )])        # membership 개수 정해주기
    if (kwargs ["LOWVAF_RATIO"] != 0):
        membership_count [0] = int (NUM_MUTATION *  ( kwargs["LOWVAF_RATIO"] ) )
        
    membership_count = collections.OrderedDict ( sorted (membership_count.items()) )                        # key 순서대로 sort 해주기
    #print ( membership_count  )
    
    for i in range (NUM_BLOCK):
        #print ("\n[{}번때 sample(block)]".format (i))

        #1. 일단 Dirichlet sampling에 들어갈 alpha 계수를 정해준다
        li = np.zeros ( NUM_CLONE_ACTIVE, dtype = "int")
        li_2 = np.zeros ( NUM_CLONE_ACTIVE, dtype = "int")
        for j in range (NUM_CLONE_ACTIVE):
            if kwargs["SimData"] == "decoy":
                if (kwargs["LOWVAF_RATIO"] != 0) & (j == 0) & ( random.randrange (0, NUM_BLOCK) != 0):
                    li[j] = 6
                else:
                    li[j] = 15 + (j * 10)
            elif kwargs["SimData"] == "lump":
                if (kwargs["LOWVAF_RATIO"] != 0) & (j == 0) & ( random.randrange (0, NUM_BLOCK) != 0) :
                    li[j] = 6
                else:
                    li[j] = (NUM_CLONE) * 10 - (j * 6)

        random.shuffle (li)       # 섞어주자
        for j in range (NUM_CLONE_ACTIVE):     
            while True:
                li_2[j] = np.random.binomial  ( li[j], 0.5 )
                if li_2[j] >= 1:
                    break
                else:
                    continue
        li = copy.deepcopy ( li_2 ) 
                
        #2. Dirichlet sampling 수행
        s  = np.random.dirichlet ( li, int (NUM_MUTATION * ( 1 - kwargs["FP_RATIO"] )) )

        #3. membership 개수만큼 np_vaf에 집어넣어주기
        k_now = 0
        for j, j_count in membership_count.items():
            for k in range (k_now, k_now + j_count):
                np_vaf [k][i] = s [k - k_now][j] / 2                    # j번째 clone의 vaf 정보를 넣어주기  (vaf니까 나누기 2)
                membership[k] = j
                mutation_id[k] = k
            k_now = k_now + j_count

        # print ( "alpha coefficient (Dirichlet distribution) : {}".format(li) )
        # ttt = 0
        # for j in range ( s.shape[1] ):
        #     mixture[i][j] = np.mean ( s [: membership_count[j], j] )
        #     print ("clone {} (0 ~ {})의 mean : {}".format ( j, membership_count[j], mixture[i][j] ) )
        #     ttt = ttt + np.mean ( s [: membership_count[j], j] )
        # print ("total mean : {}".format ( ttt ))


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
    k = int (NUM_MUTATION * ( 1 - kwargs["FP_RATIO"] ))
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

    #print ("INPUT_TSV\t{}\nNPVAF_DIR\t{}".format (kwargs["INPUT_TSV"], kwargs["NPVAF_DIR"]))

    return df, np_vaf, mixture, membership, mutation_id
    

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

    global NUM_CLONE, NUM_BLOCK, NUM_MUTATION, mixture, membership, df, np_vaf, mutation_id, samplename_dict

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
    
    os.system ("rm -rf " + kwargs ["NPVAF_DIR"] + "/preparation.npvaf.txt" )
    #os.system ("rm -rf " + kwargs ["NPVAF_DIR"] + "/" + str(NUM_BLOCK) + "D_clone" + str(NUM_CLONE) + "_" + str(ii) + ".npvaf" )
    #os.system ("mkdir -p " + kwargs ["NPVAF_DIR"] )

    output_file_inputtsv = open( kwargs ["INPUT_TSV"], "w" )
    #output_file_npvaf = open( kwargs ["NPVAF_DIR"] + "/" + str(NUM_BLOCK) + "D_clone" + str(NUM_CLONE) + "_" + str(ii) + ".npvaf", "w" )
    output_file_npvaf = open( kwargs ["NPVAF_DIR"] + "/preparation.npvaf.txt" , "w" )


    df, np_vaf, mixture, membership, mutation_id = dirichlet_sampling ( **kwargs )
    
    samplename_dict = { i:i for i in range (NUM_BLOCK) }


    def printresult():
        mixture_sort = np.zeros ((NUM_BLOCK, NUM_CLONE), dtype = 'float') 
        membership_sort = []
        dict={}
        num = 0

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
            print ( "{}".format(df[k][i]["membership_answer"]), file = output_file_npvaf)
        output_file_npvaf.close()


        # print df (INPUT_TSV)
        for k in range(0, NUM_MUTATION):
            print ("mut_{}\t{}\t".format( k, df[k][i]["membership_answer"] ) , end = "", file = output_file_inputtsv)
            for i in range (NUM_BLOCK - 1):
                print ( "{},{}".format( df[k][i]["depth"], df[k][i]["alt"] ), end = ",", file = output_file_inputtsv)
            print ( "{},{}".format( df[k][NUM_BLOCK - 1]["depth"], df[k][NUM_BLOCK - 1]["alt"] ), file = output_file_inputtsv)
        output_file_inputtsv.close()


    printresult()

    #return df, np_vaf, mixture, membership, mutation_id, samplename_dict
