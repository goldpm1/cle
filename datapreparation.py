import numpy as np
import pandas as pd 
import random

def makedf ( **kwargs ):
    global  df, inputdf, input_containpos, np_vaf, np_BQ, membership, mutation_id, depth_list, membership_answer, mutation_id, mixture_answer, parent_type, parent_type_selected

    input_containpos = pd.read_csv(kwargs["INPUT_TSV"],  header = None,  sep = "\t") 

    if input_containpos.shape[1] == 3: #  If 4th column (BQ) is absent
        input_containpos.columns = ["pos", "sample", "info"]
        input_containpos.astype ({"pos":"str", "sample":"str", "info":"str"})
    elif input_containpos.shape[1] == 4: #  If 4th column (BQ) is present
        input_containpos.columns = ["pos", "sample", "info", "BQ"]
        input_containpos.astype ({"pos":"str", "sample":"str", "info":"str", "BQ":"str"})


    input_containpos ["cha1"] = "exclusive"  # monoclone이면 exclusive, 둘 이상의 clone이 합쳐진거면 parent
    input_containpos ["cha2"] = "space"       # 축 상에 있으면 axis, 공간 상에 있으면 space
    kwargs["samplename_dict_CharacterToNum"] = {}
    kwargs["samplename_dict_NumToCharacter"] = {}
    kwargs["NUM_MUTATION"] = input_containpos.shape[0]

    np_vaf = np.zeros( (kwargs["NUM_MUTATION"], kwargs["NUM_BLOCK"] ), dtype = 'float')
    np_BQ = np.zeros( (kwargs["NUM_MUTATION"], kwargs["NUM_BLOCK"] ), dtype = 'float')
    np_BQ.fill(20)  # default BQ is 20
    inputdf = pd.DataFrame (np.zeros((kwargs["NUM_MUTATION"], kwargs["NUM_BLOCK"]), dtype = 'object'), columns = ['block' + str(i + 1) for i in range(kwargs["NUM_BLOCK"] )])
    mutation_id = []
    membership = []
    depth_list = []
    
    # input 형식은 n * 3 으로 구성 :   ID (chr_pos), membmership(정답 set 일 경우),  kwargs["kwargs["NUM_BLOCK"]"](3)만큼의 depth, alt 정보

    depth_col = [[]] * int(len(input_containpos.iloc[0][2].split(","))/2)
    depth_row = []
    for row in range(kwargs["NUM_MUTATION"]):
        depth_row_mini = []
        mutation_id.append( str(input_containpos.iloc[row][0]) )            # "pos"
        membership.append( str(input_containpos.iloc[row][1]) )           # "sample"
        if "," in str(input_containpos.iloc[row][1]) :
            input_containpos.loc[row,"cha1"] = "parent"

        if str(input_containpos.iloc[row][1]) not in kwargs["samplename_dict_CharacterToNum"].keys():
            kwargs["samplename_dict_CharacterToNum"][str(input_containpos.iloc[row][1])] = int (len(kwargs["samplename_dict_CharacterToNum"]))      # {'other': 0, 'V5': 1, 'V3': 2, 'V1': 3}           # 각 sample name을 숫자화시킴

        #rmv_bracket = re.sub("[\[\] ]", '', str(input_containpos.iloc[row][2])).split(",")            # [194, 25, 193, 66, 0, 0] 라고 되어 있는데 bracket과 한 칸 공백을 지움
        rmv_bracket=input_containpos.iloc[row][2].split(",")
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
                depth_row_mini.append(depth)
                depth_col[col].append(depth)

            if "BQ" in input_containpos.columns: #  If 4th column (BQ) is present
                if kwargs["NUM_BLOCK"] == 1:
                    BQ_input = [input_containpos.iloc[row][3]]     
                else:
                    BQ_input = input_containpos.iloc[row][3].split(",")   # 4th column : 20,20
                for i in range (0, len(BQ_input) ):
                    np_BQ [row][i] = int ( BQ_input[i] )
                    if int ( BQ_input[i] ) == 0:     # Set BQ=20 in axis (FN) variant
                        np_BQ [row][i] = 20 

        depth_row.append (depth_row_mini)

    # "0.0.0"을 그대로 놔둘 수 없다.  평균 depth로 갈음해서 바꿔 넣는다  (alt는 0으로 유지)

    for row in range( kwargs["NUM_MUTATION"] ):
        for  i in range(0, len(rmv_bracket), 2 ):
            col = int(i / 2)
            if inputdf.iloc[row][col] == "0:0:0":
                inputdf.iloc[row][col] = str(round(np.mean(depth_col[col]))) + ":" + str(round(np.mean(depth_col[col]))) + ":0"
                input_containpos.loc[row,"cha2"] = "axis"
        depth_list.append(np.mean(depth_row[row]))

    df = [[None] * kwargs["NUM_BLOCK"] for i in range(kwargs["NUM_MUTATION"])]
    for row in range (kwargs["NUM_MUTATION"]):
        for col in range ( kwargs["NUM_BLOCK"] ):
            df[row][col] = {"depth":int(inputdf.iloc[row][col].split(":")[0]), "ref":int(inputdf.iloc[row][col].split(":")[1]), "alt":int(inputdf.iloc[row][col].split(":")[2])}
            if df[row][col]["depth"] == 0:
                print (df[row][col], row, col)




    # 일단 기존 data에서 answer별로 count가 얼마나 되는지 출력하기
    df_counts = pd.DataFrame (np.unique ( [membership[i] for i in [i for i in range (0, kwargs["NUM_MUTATION"]) if (  (depth_list[i] > kwargs["DEPTH_CUTOFF"]) ) ] ], return_counts = True ))
    #df_counts = df_counts.loc[1:, list(kwargs["samplename_dict_CharacterToNum"].keys())]
    for s_index, s in enumerate(df_counts.iloc[0]):
        kwargs["samplename_dict_CharacterToNum"][s] = s_index
        kwargs["samplename_dict_NumToCharacter"][s_index] = s

    print ("조정 전")
    for i in range(df_counts.shape[0]):
        print ("\t", end = "")
        for j in range(df_counts.shape[1]) :
            print (df_counts.iloc[i][j], end ="\t")
        print ("")


    # PARENT_ RATIO가 아닌 kwargs["NUM_PARENT"]만큼 뽑고 싶을 때
    lowvaf_clone = []
    NUM_CLONE = len(set(membership))
    mixture_pre = np.zeros ((kwargs["NUM_BLOCK"], NUM_CLONE), dtype = 'float')     #mixture 값을 일단 초기화
    for i in range(kwargs["NUM_BLOCK"]):
        for j in range(NUM_CLONE):
            mixture_pre[i][j] = round(np.mean(np_vaf[[x  for x in range(len(membership)) if membership[x] == list(kwargs["samplename_dict_CharacterToNum"].keys())[j]   ]] [: , i] * 2), 3)        
        if np.all (mixture_pre[i] < 0.06):
            lowvaf_clone.append ( kwargs["samplename_dict_CharacterToNum"][i] )
            

    if kwargs["NUM_PARENT"] != 0:
        parent_type, count = np.unique ( [membership[i] for i in [i for i in range (0, kwargs["NUM_MUTATION"]) if ( ( "parent" in input_containpos.loc[i,"cha1"] ) & (depth_list[i] > kwargs["DEPTH_CUTOFF"]) ) ] ], return_counts = True )
        count_sort_ind = np.argsort(-count)
        parent_type, count = parent_type[count_sort_ind] , count[count_sort_ind]    # 숫자대로 내림차순 정렬

        parent_type_selected = []
        for k in range ( len ( [i for i in kwargs["samplename_dict_CharacterToNum"].keys()  if  ("," in i)   ]  ) ):
            if k == kwargs["NUM_PARENT"]:
                break
            check = 0
            for j in lowvaf_clone:
                if ( j in parent_type[k] ):
                    print ("{}는 lowvaf_clone ({})이 더해진 parent이므로 탈락".format (parent_type[k], j))
                    check = 1
                    break
                
            if check == 0:  # 위에서 문제가 없어야 통과
                parent_type_selected.append ( parent_type[k] )     # 개수대로 상위 kwargs["NUM_PARENT"]개의 cluster 이름을 뽑는다  (char name)
        
        print ("parent_type_selected = {}".format( parent_type_selected))  
        if len (parent_type_selected) == 0:  #위에서 전혀 못 뽑는다면...
            parent_type_selected = parent_type [0 : kwargs["NUM_PARENT"]]


        if len ([i for i in range (0, kwargs["NUM_MUTATION"]) if ( ( membership[i] in parent_type_selected  ) & (depth_list[i] > kwargs["DEPTH_CUTOFF"]) ) ]) / kwargs["RANDOM_PICK"] < kwargs["PARENT_RATIO"]:             # depth_list : 해당 mutation의 평균 depth        # input 받은 kwargs["PARENT_RATIO"]보다 내 database에서 적게 보유하고 있으면
            kwargs["PARENT_RATIO"] = round( len ([i for i in range (0, kwargs["NUM_MUTATION"]) if ( ( "parent" in input_containpos.loc[i,"cha1"] ) & (depth_list[i] > kwargs["DEPTH_CUTOFF"]) ) ])  / kwargs["RANDOM_PICK"], 3)
    else:
        if len ([i for i in range (0, kwargs["NUM_MUTATION"]) if ( ( "parent" in input_containpos.loc[i,"cha1"]  ) & (depth_list[i] > kwargs["DEPTH_CUTOFF"]) ) ]) / kwargs["RANDOM_PICK"] < kwargs["PARENT_RATIO"]:             # depth_list : 해당 mutation의 평균 depth        # input 받은 kwargs["PARENT_RATIO"]보다 내 database에서 적게 보유하고 있으면
            kwargs["PARENT_RATIO"] = round( len ([i for i in range (0, kwargs["NUM_MUTATION"]) if ( ( "parent" in input_containpos.loc[i,"cha1"] ) & (depth_list[i] > kwargs["DEPTH_CUTOFF"]) ) ])  / kwargs["RANDOM_PICK"], 3)



    print ("FP_RATIO : {}\tFP_USEALL :{}\tDEPTH_CUTOFF : {}".format (kwargs["FP_RATIO"], kwargs["FP_USEALL"], kwargs["DEPTH_CUTOFF"]))

    try:
        p = [i for i in range (0, kwargs["NUM_MUTATION"]) if ( ( "axis" in input_containpos.loc[i,"cha2"] ) & (depth_list[i] > kwargs["DEPTH_CUTOFF"]) ) ]
        if len (p) / kwargs["RANDOM_PICK"] < kwargs["AXIS_RATIO"]:      # input 받은 kwargs["AXIS_RATIO"]보다 내 database에서 적게 보유하고 있으면
            kwargs["AXIS_RATIO"] = round( len ([i for i in range (0, kwargs["NUM_MUTATION"]) if ( ( "axis" in input_containpos.loc[i,"cha2"] ) & (depth_list[i] > kwargs["DEPTH_CUTOFF"]) ) ]) / kwargs["RANDOM_PICK"], 3)
    except:
        print ("AXIS가 없음")
        kwargs["AXIS_RATIO"] = 0

    try:
        p = [i for i in range (0, kwargs["NUM_MUTATION"]) if ( ( "FP" in input_containpos.loc[i,"sample"] ) & (depth_list[i] > kwargs["DEPTH_CUTOFF"]) ) ]
        if len (p) / kwargs["RANDOM_PICK"] < kwargs["FP_RATIO"]:      # input 받은 kwargs["FP_RATIO"]보다 내 database에서 적게 보유하고 있으면
            kwargs["FP_RATIO"] = round( len ([i for i in range (0, kwargs["NUM_MUTATION"]) if ( ( "FP" in input_containpos.loc[i,"sample"] ) & (depth_list[i] > kwargs["DEPTH_CUTOFF"]) ) ]) / kwargs["RANDOM_PICK"], 3)
    except:
        print ("FP가 없음")
        kwargs["FP_RATIO"] = 0

    return kwargs
    


def RANDOM_PICK_fun(**kwargs):
    global  df, inputdf, input_containpos, np_vaf, np_BQ, membership, mutation_id, depth_list, membership_answer, mixture_answer, parent_type, parent_type_selected

    # RANDOM하게 1000개만 뽑기. FP 의 비율도 조절해서...

    random.seed(kwargs["RANDOM_SEED"])

    if ((kwargs["FP_USEALL"] != "True") & (kwargs["FP_RATIO"] != 0)):
        FP_index = [i for i in range (0, kwargs["NUM_MUTATION"]) if ( (membership[i] == 'FP') & (depth_list[i] > kwargs["DEPTH_CUTOFF"]) ) ]        # Depth 너무 낮으면 곤란하니 input을 받자
    elif kwargs["FP_USEALL"] == "True":
        #FP_index = [i for i in range (0, kwargs["NUM_MUTATION"]) if ( (membership[i] == 'FP') & ( input_containpos.loc[i,"cha2"] == "space" ) &  (depth_list[i] > kwargs["DEPTH_CUTOFF"]) ) ]        # FP는 space에 있는 것만 밭자
        FP_index = [i for i in range (0, kwargs["NUM_MUTATION"]) if ( (membership[i] == 'FP') & (depth_list[i] > kwargs["DEPTH_CUTOFF"]) ) ]        # 굳이 space에만 있는 것을 뽑을 필요는 없다


    if kwargs["NUM_PARENT"] != 0:
        PARENT_index = [i for i in range (0, kwargs["NUM_MUTATION"]) if ( ( membership[i] in parent_type_selected  ) & (depth_list[i] > kwargs["DEPTH_CUTOFF"]) ) ]        #    "parent" in input_containpos.loc[i,"cha1"] 
    else:
        PARENT_index = [i for i in range (0, kwargs["NUM_MUTATION"]) if ( ( "parent" in input_containpos.loc[i,"cha1"]  ) & (depth_list[i] > kwargs["DEPTH_CUTOFF"]) ) ]        #    "parent" in input_containpos.loc[i,"cha1"] 

    CHILD_index = [i for i in range (0, kwargs["NUM_MUTATION"]) if ( (membership[i] != 'FP')  &  ( "exclusive" in input_containpos.loc[i,"cha1"] )  & (depth_list[i] > kwargs["DEPTH_CUTOFF"]) ) ]               #
    AXIS_index = [i for i in range (0, kwargs["NUM_MUTATION"]) if ( ( "axis" in input_containpos.loc[i,"cha2"] ) & (depth_list[i] > kwargs["DEPTH_CUTOFF"]) ) ]        #


    CHILD_SPACE_index = list(set(CHILD_index) - set(AXIS_index))      # CHILD 이면서 SPACE에 위치해 있는 것들
    CHILD_AXIS_index = list(set(CHILD_index) & set(AXIS_index))       # CHILD 이면서  AXIS에 위치해 있는 것들
    
    print ( "len (CHILD_index) : {}\tlen (CHILD_SPACE_index) : {}".format (len (CHILD_index), len (CHILD_SPACE_index)))
    
    if (kwargs["FP_USEALL"] == "True"):        # kwargs["FP"] == True면 full로 뽑는다
        FP_randomsample = FP_index
        kwargs["FP_RATIO"] = len(FP_index) / kwargs["RANDOM_PICK"]
    elif (kwargs["FP_USEALL"] == "False") & (kwargs["FP_RATIO"] == 0):      # 둘다 False, 0이면 하나도 안 뽑느다
        FP_randomsample = FP_index = []
        kwargs["FP_RATIO"] = 0
    elif (kwargs["FP_USEALL"] == "False") & (kwargs["FP_RATIO"] != 0):
        try:
            FP_randomsample = random.sample(FP_index, int(kwargs["RANDOM_PICK"] * (kwargs["FP_RATIO"])))
        except:
            print (kwargs["INPUT_TSV"] + " - Can't extract FP data as requested")
            FP_randomsample = FP_index
            kwargs["FP_RATIO"] = len(FP_index) / kwargs["RANDOM_PICK"]
            return False
    
    if kwargs["NUM_PARENT"] == 0:     # 확률로 뽑고싶을 때
        try:
            PARENT_randomsample = random.sample(PARENT_index, int(kwargs["RANDOM_PICK"] * (kwargs["PARENT_RATIO"])))
        except:
            print (kwargs["INPUT_TSV"] + " - Can't extract Parent data as requested")
            PARENT_randomsample = PARENT_index
            kwargs["PARENT_RATIO"] = len(PARENT_index) / kwargs["RANDOM_PICK"]
            return False
    else:   # 상위 몇 개 cluster를 뽑고싶을 때
        PARENT_index = [i for i in range (0, kwargs["NUM_MUTATION"]) if ( ( membership[i] in parent_type_selected  ) & (depth_list[i] > kwargs["DEPTH_CUTOFF"]) ) ]        #    "parent" in input_containpos.loc[i,"cha1"] 
        PARENT_randomsample = PARENT_index
        kwargs["PARENT_RATIO"] = len(PARENT_index) / kwargs["RANDOM_PICK"]

    try:
        CHILD_AXIS_randomsample = random.sample(CHILD_AXIS_index, int(kwargs["RANDOM_PICK"] * (kwargs["AXIS_RATIO"])))
    except:
        print (kwargs["INPUT_TSV"] + " - Can't extract Child_axis data as requested")
        return False, kwargs


    print ( "RANDOM_PICK = {}\tFP_randomsample = {}\tPARENT_randomsample = {}\tCHILD_AXIS_randomsample = {}".format ( kwargs["RANDOM_PICK"], len(FP_randomsample), len(PARENT_randomsample), len (CHILD_AXIS_randomsample) ) )

    try:
        CHILD_SPACE_randomsample = random.sample(CHILD_SPACE_index, kwargs["RANDOM_PICK"] - (len(FP_randomsample) + len(PARENT_randomsample) + len(CHILD_AXIS_randomsample)))  
    except:
        print (kwargs["INPUT_TSV"] + " - Can't extract Child_space data as requested")
        return False, kwargs
    

    print ("조정 후\n\tFP 개수 : {}\tPARENT 개수 : {}\tCHILD_SPACE 개수 : {}\tCHILD_AXIS 개수 : {}".format(len(FP_randomsample), len (PARENT_randomsample), len (CHILD_SPACE_randomsample), len (CHILD_AXIS_randomsample)))



    random_index = sorted( FP_randomsample + PARENT_randomsample + CHILD_SPACE_randomsample + CHILD_AXIS_randomsample )             # 다 합치면 RADOM_PICK 개수가 되겠지

    # kwargs["RANDOM_PICK"] 개만으로 재정비
    input_containpos =  input_containpos.iloc[random_index]
    #input_containpos.to_csv (OUTPUT_DIR + "sampling_{0}.df.tsv".format(kwargs["RANDOM_PICK"]), index = True, header=True, sep = "\t")

    inputdf  = inputdf.iloc[random_index]
    df = [df[i] for i in random_index]
    np_vaf = np_vaf[random_index]
    membership_answer = [membership[i] for i in random_index]
    mutation_id = [mutation_id[i] for i in random_index]

    #np_vaf + membership를 df 형식으로 하고 kwargs["RANDOM_PICK"]개만 출력 
    t = pd.DataFrame(np_vaf, columns = ["block{0}".format(i) for i in range(kwargs["NUM_BLOCK"])], index = mutation_id)
    t["membership_answer"] = pd.Series(membership_answer, index = mutation_id)
    t.to_csv ("{0}/npvaf.txt".format( kwargs["NPVAF_DIR"] ), index = True, header=True, sep = "\t")


    kwargs["samplename_dict_CharacterToNum"], cnt = {}, 0
    for k in membership_answer:
        if k not in kwargs["samplename_dict_CharacterToNum"].keys():
            kwargs["samplename_dict_CharacterToNum"][k] = cnt
            cnt = cnt + 1
    
    NUM_CLONE = len(set(membership_answer))
    mixture_answer = np.zeros ((kwargs["NUM_BLOCK"], NUM_CLONE), dtype = 'float')     #mixture 값을 일단 초기화
    for i in range(kwargs["NUM_BLOCK"]):
        for j in range(NUM_CLONE):
            mixture_answer[i][j] = round(np.mean(np_vaf[[x  for x in range(len(membership_answer)) if membership_answer[x] == list(kwargs["samplename_dict_CharacterToNum"].keys())[j]   ]] [: , i] * 2), 3)


    return True, kwargs


def pyclone_dataset( DIR, **kwargs ):
    for col in range (kwargs["NUM_BLOCK"]):
        PYCLONE_OUTPUT=DIR + "/block" + str(col) + ".tsv"
        with open (PYCLONE_OUTPUT, "w", encoding = "utf8") as output_pyclone:
            print ("\t".join(["mutation_id", "ref_counts", "var_counts", "normal_cn", "minor_cn", "major_cn"]), file = output_pyclone)
            for row in range( kwargs["RANDOM_PICK"] ):
                pyclone_row = [mutation_id[row] , str(df[row][col]["ref"]), str(df[row][col]["alt"]), str(2), str(1), str(1)]
                print("\t".join(pyclone_row), file = output_pyclone)

def pyclonevi_dataset( DIR , **kwargs):
    PYCLONE_VI_OUTPUT=DIR + "/input.tsv"
    with open (PYCLONE_VI_OUTPUT, "w", encoding = "utf8") as output_pyclone_vi:
        print ("\t".join(["mutation_id", "sample_id", "ref_counts", "alt_counts", "normal_cn", "major_cn", "minor_cn", "tumour_content"]), file = output_pyclone_vi)
        for row in range(kwargs["RANDOM_PICK"]):
            for col in range(kwargs["NUM_BLOCK"]):
                # print (mutation_id [row])
                # print (df[row][col] )
                pyclone_vi_row = [mutation_id[row] , "block" + str(col), str(df[row][col]["ref"]), str(df[row][col]["alt"]), str(2), str(1), str(1), str(1.0)]
                print("\t".join(pyclone_vi_row), file = output_pyclone_vi)

def sciclone_dataset( DIR, **kwargs ):
    import re

    for col in range (kwargs["NUM_BLOCK"]):
        SCICLONE_OUTPUT=DIR + "/block" + str(col) + ".dat"
        with open (SCICLONE_OUTPUT, "w", encoding = "utf8") as output_sciclone:
            print ("\t".join(["chr", "pos", "ref_reads", "var_reads", "vaf"]), file = output_sciclone)
            for row in range(kwargs["RANDOM_PICK"]):
                sciclone_row = [re.split(r'[_ :]', mutation_id[row])[0], re.split(r'[_ :]', mutation_id[row])[1]  , str(df[row][col]["ref"]), str(df[row][col]["alt"]), str(round((df[row][col]["alt"] / ( df[row][col]["alt"] + df[row][col]["ref"] )), 3) * 100)  ]
                print("\t".join(sciclone_row), file = output_sciclone)


def quantumclone_dataset( DIR, **kwargs ):
    import re

    for col in range (kwargs["NUM_BLOCK"]):
        QUANTUMCLONE_OUTPUT=DIR + "/block" + str(col) + ".dat"
        with open (QUANTUMCLONE_OUTPUT, "w", encoding = "utf8") as output_quantumclone:
            print ("\t".join(["Sample", "SampleName", "Chr", "Start", "Alt", "Depth", "Genotype"]), file = output_quantumclone)
            for row in range(kwargs["RANDOM_PICK"]):
                quantumclone_row = ["Sample" + str(col), "Sample" + str(col) , re.split(r'[_ :]', mutation_id[row])[0], re.split(r'[_ :]', mutation_id[row])[1]  , str(df[row][col]["alt"]), str(df[row][col]["ref"] + df[row][col]["alt"]),  "AB"  ]
                print("\t".join(quantumclone_row), file = output_quantumclone)



def main (**kwargs):
    global  df, inputdf, input_containpos, np_vaf, np_BQ, membership, mutation_id, depth_list, membership_answer,  mixture_answer, parent_type, parent_type_selected


    kwargs = makedf( **kwargs)
    check, kwargs = RANDOM_PICK_fun(**kwargs)

    if check == True:
        if type(kwargs["PYCLONEVI_DIR"]) == type("string"):
            pyclonevi_dataset( kwargs["PYCLONEVI_DIR"], **kwargs )
        if type(kwargs["SCICLONE_DIR"]) == type("string"):
            sciclone_dataset( kwargs["SCICLONE_DIR"], **kwargs )
        if type(kwargs["QUANTUMCLONE_DIR"]) == type("string"):
            quantumclone_dataset( kwargs["QUANTUMCLONE_DIR"], **kwargs )

            
    return (input_containpos, inputdf, df, np_vaf, np_BQ, membership_answer, mixture_answer, mutation_id, kwargs)
