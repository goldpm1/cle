import numpy as np
import pandas as pd 
import re
import random
import scipy

def makedf ():
    global  df, inputdf, np_vaf, membership, mutation_id, depth_list, NUM_MUTATION, samplename_dict

    input_containpos = pd.read_csv(INPUT_DIR + "M3_MT_MT2_all_input.txt",  header = None, names =["pos", "sample", "info"], sep = "\t") 
    samplename_dict = {}
    NUM_MUTATION = input_containpos.shape[0]

    np_vaf = np.zeros((NUM_MUTATION, NUM_BLOCK_INPUT), dtype = 'float')
    inputdf = pd.DataFrame (np.zeros((NUM_MUTATION, NUM_BLOCK_INPUT), dtype = 'object'), columns = ['block' + str(i + 1) for i in range(NUM_BLOCK_INPUT)])
    mutation_id = []
    membership = []
    depth_list = []
    
    # input 형식은 n * 3 으로 구성 :   ID (chr_pos), membmership(정답 set 일 경우),  NUM_BLOCK_INPUT(3)만큼의 depth, alt 정보

    for row in range(NUM_MUTATION):
        mutation_id.append( str(input_containpos.iloc[row][0]) )
        membership.append( str(input_containpos.iloc[row][1]) )
        if str(input_containpos.iloc[row][1]) not in samplename_dict.keys():
            samplename_dict[str(input_containpos.iloc[row][1])] = int (len(samplename_dict))      # {'other': 0, 'V5': 1, 'V3': 2, 'V1': 3}

        rmv_bracket = re.sub("[\[\] ]", '', str(input_containpos.iloc[row][2])).split(",")            # [194, 25, 193, 66, 0, 0] 라고 되어 있는데 bracket과 한 칸 공백을 지움
        depth_row=[]
        for i in range(0, len(rmv_bracket), 2 ):
            depth = int(rmv_bracket[i])
            alt = int(rmv_bracket[i+1])
            ref = depth - alt

            col = int(i / 2)

            if depth == 0:
                np_vaf[row][col] = 0
                inputdf.iloc[row][col] = "0:0:0"
            else:    
                np_vaf[row][col] = round (alt / depth , 2)
                inputdf.iloc[row][col] = str(depth) + ":" + str(ref) + ":" + str(alt)
                depth_row.append(depth)

        # "0.0.0"을 그대로 놔둘 수 없다.  평균 depth로 갈음해서 바꿔 넣는다
        for  i in range(0, len(rmv_bracket), 2 ):
            col = int(i / 2)
            if inputdf.iloc[row][col] == "0:0:0":
                inputdf.iloc[row][col] = str(round(np.mean(depth_row))) + ":" + str(round(np.mean(depth_row))) + ":0"
        
        depth_list.append(np.mean(depth_row))

    df = [[None] * NUM_BLOCK for i in range(inputdf.shape[0])]
    for row in range (inputdf.shape[0]):
        for col in range (NUM_BLOCK):
            df[row][col] = {"depth":int(inputdf.iloc[row][col].split(":")[0]), "ref":int(inputdf.iloc[row][col].split(":")[1]), "alt":int(inputdf.iloc[row][col].split(":")[2])}
            if df[row][col]["depth"] == 0:
                print (df[row][col], row, col)



def random_pick_fun():
    global  df, inputdf, np_vaf, membership, mutation_id, depth_list, NUM_MUTATION, NUM_CLONE, NUM_BLOCK
    global membership_answer, mutation_id, mixture_answer

    # RANDOM하게 1000개만 뽑기. FP 의 비율도 조절해서...

    TP_index = [i for i in range (0, NUM_MUTATION) if ( ( 'het' in membership[i] ) & (depth_list[i] > 800) ) ]
    FP_index = [i for i in range (0, NUM_MUTATION) if ( (membership[i] == 'FP') & (depth_list[i] > 800) ) ]
    random_index = sorted(random.sample(TP_index, int(RANDOM_PICK * (1-FP_RATIO))) + random.sample(FP_index, int(RANDOM_PICK * (FP_RATIO))))
    inputdf  = inputdf.iloc[random_index]
    df = [df[i] for i in random_index]
    np_vaf = np_vaf[random_index]
    membership_answer = [membership[i] for i in random_index]
    mutation_id = [mutation_id[i] for i in random_index]

    #np_vaf + membership를 df 형식으로 하고 1000개만 출력 
    t = pd.DataFrame(np_vaf, columns = ["block{0}".format(i) for i in range(NUM_BLOCK_INPUT)], index = mutation_id)
    t["membership_answer"] = pd.Series(membership_answer, index = mutation_id)
    t.to_csv (INPUT_DIR + "sampling_{0}.txt".format(RANDOM_PICK), index = True, header=True, sep = "\t")


    NUM_CLONE = len(set(membership_answer))
    mixture_answer = np.zeros ((NUM_BLOCK, NUM_CLONE), dtype = 'float')     #mixture 값을 일단 초기화
    for i in range(NUM_BLOCK):
        for j in range(NUM_CLONE):
            #print (i, j, list(set(membership_answer))[j])
            mixture_answer[i][j] = round(np.mean(np_vaf[[x  for x in range(len(membership_answer)) if membership_answer[x] == list(set(membership_answer))[j]]][:,i] * 2), 2)

def pyclone_dataset():
    for col in range (NUM_BLOCK):
        PYCLONE_OUTPUT="/data/project/Alzheimer/EM_cluster/old/pilot/04.EM_input/pyclone/block" + str(col) + ".tsv"
        with open (PYCLONE_OUTPUT, "w", encoding = "utf8") as output_pyclone:
            print ("\t".join(["mutation_id", "ref_counts", "var_counts", "normal_cn", "minor_cn", "major_cn"]), file = output_pyclone)
            for row in range(inputdf.shape[0]):
                pyclone_row = [mutation_id[row] , str(df[row][col]["ref"]), str(df[row][col]["alt"]), str(2), str(1), str(1)]
                print("\t".join(pyclone_row), file = output_pyclone)

def pyclonevi_dataset():
    PYCLONE_VI_OUTPUT="/data/project/Alzheimer/EM_cluster/old/pilot/04.EM_input/pyclone_vi/pyclone_vi_220610.tsv"
    with open (PYCLONE_VI_OUTPUT, "w", encoding = "utf8") as output_pyclone_vi:
        print ("\t".join(["mutation_id", "sample_id", "ref_counts", "alt_counts", "normal_cn", "major_cn", "minor_cn", "tumour_content"]), file = output_pyclone_vi)
        for row in range(inputdf.shape[0]):
            for col in range(NUM_BLOCK):
                pyclone_vi_row = [mutation_id[row] , "block" + str(col), str(df[row][col]["ref"]), str(df[row][col]["alt"]), str(2), str(1), str(1), str(1.0)]
                print("\t".join(pyclone_vi_row), file = output_pyclone_vi)




def main (**kwargs):
    global NUM_BLOCK_INPUT, NUM_BLOCK, RANDOM_PICK, FP_RATIO, INPUT_DIR, OUTPUT_DIR, mixture_answer, membership_answer

    NUM_BLOCK_INPUT = kwargs["NUM_BLOCK_INPUT"] 
    NUM_BLOCK = kwargs["NUM_BLOCK"] 
    RANDOM_PICK = kwargs["RANDOM_PICK"]
    FP_RATIO = kwargs["FP_RATIO"]
    INPUT_DIR = "/data/project/Alzheimer/EM_cluster/old/pilot/04.EM_input/"
    OUTPUT_DIR ="./output/"

    makedf()
    random_pick_fun()
    pyclone_dataset()
    pyclonevi_dataset()

    return (inputdf, df, np_vaf, membership_answer, mixture_answer,  mutation_id, samplename_dict)