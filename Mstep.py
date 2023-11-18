def main(input_containpos, df, np_vaf, np_BQ, step, option, **kwargs):  
    import math, re
    import isparent, miscellaneous, Estep
    import visualizationsinglesoft,  visualizationeachstep
    import numpy as np

    NUM_BLOCK, NUM_CLONE, NUM_MUTATION = kwargs["NUM_BLOCK"], kwargs["NUM_CLONE"], kwargs["NUM_MUTATION"]
    NUM_MUTATION = kwargs["RANDOM_PICK"]

    kwargs["OPTION"] = option

    if option in ["Hard", "hard"]:
        ############################### HARD CLUSTERING ##############################
        for j in range(NUM_CLONE):
            ind_list = np.where(step.membership == j)[0]   # Find the index  where membership == j
            #check = 0
            for i in range(NUM_BLOCK):
                sum_depth, sum_alt = 0, 0
                for ind in ind_list:       # Summing depth and alt
                    if j in step.makeone_index:  
                        if df[ind][i]["alt"] == 0:  # TP cluster 에서는 FN을 평균치 계산에 넣지 않는다
                            continue
                    if j == step.fp_index:  
                        step.mixture[i][j] = 0   # FP cluster는 무조건 원점으로 박아 넣는다
                        continue
                    if j in step.tn_index:
                        zero_dimension = np.where( np.round ( step.mixture [:, j], 2 ) == 0 )[0]   # 0인 축 (sample)
                        # if check == 0:
                        #     check = 1
                        if ( i in zero_dimension )  &  (df[ind][i]["alt"] != 0):    # TN : 0이어야 하는 축에서 0이 아니면 넘기자
                            continue
                    
                    if (kwargs["SEX"] == "M") & ( bool(re.search(r'X|Y', input_containpos.iloc[ind]["pos"]))  == True  ) :
                        sum_depth = sum_depth + df[ind][i]["depth"]
                        sum_alt = sum_alt + int( df[ind][i]["alt"] / 2 )
                    else: #Most of the cases
                        sum_depth = sum_depth + df[ind][i]["depth"]
                        sum_alt = sum_alt + df[ind][i]["alt"]
                    
                step.mixture[i][j] = round((sum_alt * 2) / sum_depth, 2) if sum_depth != 0 else 0   # Ideal centroid allocation
        
        # if kwargs["VERBOSE"] >= 1:
        #     print ("\t\t\t", end = " " )
        #     print( ",".join( str(row) for row in step.mixture ))

        p_list, step, condition = isparent.makeone(input_containpos, df, np_vaf, np_BQ, step, **kwargs)


        if step.makeone_index == []:   # if failed  (첫 3개는 웬만하면 봐주려고 하지만, 그래도 잘 안될 때)
            step.likelihood = step.likelihood_record[kwargs["STEP"]] = -9999999
            step.checkall_lenient =  step.checkall_strict = False
            sum_mixture = np.sum (step.mixture, axis = 1)
            if kwargs["VERBOSE"] >= 1:
                print ("\t\t\tMstep.py : condition = {}\tcheckall({}) = False\tUnable to make 1 or appropriate superclone-clone structure ({}) ".format ( condition, condition, "\t".join(str(np.round(row, 2)) for row in step.mixture )  ) )

        else:    # most of the cases
            # step.checkall_strict, sum_mixture =  miscellaneous.checkall (step, "strict", np_vaf, **kwargs) 
            kwargs["CHECKALLFROM"] = "Mstep.py"
            if kwargs["MAKEONE_STRICT"] == 3:       # BioData에서는 극단적으로 lenient하게 잡아줘야 하니까
                step.checkall_lenient, sum_mixture, p =  miscellaneous.checkall (step, "lenient", np_vaf, **kwargs) 
                step.checkall_strict, sum_mixture, p =  miscellaneous.checkall (step, "strict", np_vaf, **kwargs) 
            else:
                step.checkall_lenient, sum_mixture, p =  miscellaneous.checkallttest (step, "lenient", np_vaf, **kwargs) 
                step.checkall_strict, sum_mixture, p =  miscellaneous.checkallttest (step, "strict", np_vaf, **kwargs) 

            if kwargs["STEP"] <= (kwargs["COMPULSORY_NORMALIZATION"] - 1):
                if kwargs["VERBOSE"] >= 1:
                    print ("\t\t\tMstep.py : condition = lenient, checkall (lenient) = {}, checkall (strict) = {}".format(step.checkall_lenient, step.checkall_strict), end = "\t") 
                    print( "(sum = {})".format  ( " ".join(str( np.round(row, 2) ) for row in sum_mixture )) )

                for i in range(NUM_BLOCK):     # Normalization 
                    sum = 0
                    for j in range(NUM_CLONE):
                        if j in step.makeone_index:   
                            sum = sum + step.mixture[i][j]
                    step.mixture[i] = np.round( step.mixture[i] / sum, 2) if sum != 0 else 0   # If sum = 0, let mixture = 0

            else:
                if kwargs["VERBOSE"] >= 1:
                    print ("\t\t\tMstep.py : condition = strict, checkall (strict) = {}".format(step.checkall_strict), end = "\t") 
                    print( "(sum = {})".format  ( " ".join(str( np.round(row, 2) ) for row in sum_mixture )) )

                                
        if  kwargs["VISUALIZATION"] == True:
            if (kwargs["NUM_BLOCK"] == 1):
                visualizationeachstep.drawfigure_1d_hard(step, np_vaf, kwargs["CLEMENT_DIR"] + "/trial/clone" + str(kwargs["NUM_CLONE_ITER"]) + "." + str(kwargs["TRIAL"]) + "-" + str(kwargs["STEP_TOTAL"]) + "(hard)." + kwargs["IMAGE_FORMAT"], sum_mixture,**kwargs)
            elif (kwargs["NUM_BLOCK"] == 2):
                visualizationeachstep.drawfigure_2d(step, np_vaf, kwargs["CLEMENT_DIR"] + "/trial/clone" + str(kwargs["NUM_CLONE_ITER"]) + "." + str(kwargs["TRIAL"]) + "-" + str(kwargs["STEP_TOTAL"]) + "(hard)." + kwargs["IMAGE_FORMAT"], sum_mixture, **kwargs)
            elif (kwargs["NUM_BLOCK"] >= 3):
                #visualizationeachstep.drawfigure_3d(step, np_vaf, kwargs["CLEMENT_DIR"] + "/trial/clone" + str(kwargs["NUM_CLONE_ITER"]) + "." + str(kwargs["TRIAL"]) + "-" + str(kwargs["STEP_TOTAL"]) + "(hard)." + kwargs["IMAGE_FORMAT"], sum_mixture, **kwargs)
                visualizationeachstep.drawfigure_3d_SVD(step, np_vaf, kwargs["CLEMENT_DIR"] + "/trial/clone" + str(kwargs["NUM_CLONE_ITER"]) + "." + str(kwargs["TRIAL"]) + "-" + str(kwargs["STEP_TOTAL"]) + "(hard)." + kwargs["IMAGE_FORMAT"], sum_mixture, **kwargs)
    ###############################################################################

    ################################ SOFT CLUSTERING ##############################

    if option in ["Soft", "soft"]:
        #print ("\t\ta. Mixture (before soft clustering) : {}". format(list(step.mixture)))

        makeone_index_i = []
        for k in range(NUM_MUTATION):
            if step.membership[k] in step.makeone_index:
                makeone_index_i.append(k)

        for j in range(NUM_CLONE):
            if j not in step.makeone_index:
                for i in range(NUM_BLOCK):
                    # Summing all depth and alt
                    sum_depth, sum_alt = 0, 0
                    for ind in np.where(step.membership == j)[0]:
                        if df[ind][i]["alt"] != 0:  # TP cluster 라면 FN을 제거하기 위해 이래야 하고,  TN cluster라면 이럴 필요가 없다
                            if (kwargs["SEX"] == "M") & ( bool(re.search(r'X|Y', input_containpos.iloc[ind]["pos"]))  == True  ) :
                                sum_depth = sum_depth + df[ind][i]["depth"] * 2
                                sum_alt = sum_alt + df[ind][i]["alt"]
                            else: #Most of the cases
                                sum_depth = sum_depth + df[ind][i]["depth"]
                                sum_alt = sum_alt + df[ind][i]["alt"]
                    
                    step.mixture[i][j] = round((sum_alt * 2) / sum_depth, 2) if sum_depth != 0 else 0 

            elif j in step.makeone_index:  
                for i in range(NUM_BLOCK):   # Calculate the weighted mean
                    vaf, weight = np.zeros(NUM_MUTATION, dtype="float"), np.zeros(NUM_MUTATION, dtype="float")
                    for k in range(NUM_MUTATION):
                        if (kwargs["SEX"] == "M") & ( bool(re.search(r'X|Y', input_containpos.iloc[k]["pos"]))  == True  ) :
                            vaf[k] = int(df[k][i]["alt"])  / ( int(df[k][i]["depth"]) * 2)
                        else:
                            vaf[k] = int(df[k][i]["alt"]) / int(df[k][i]["depth"])

                        if step.membership[k] in step.makeone_index:
                            weight[k] = math.pow(10, step.membership_p[k][j])
                        #print ("{} : (clone {})   weight = {},  vaf = {}".format(k, step.membership[k], weight[k], vaf[k]))

                    step.mixture[i][j] = round(np.average(vaf[makeone_index_i], weights = weight [ makeone_index_i ]), 2) * 2

        

 
        if NUM_CLONE == 1:
            step.mixture = np.array([[1.0]] * kwargs["NUM_BLOCK"])

        p_list, step, condition = isparent.makeone(input_containpos, df, np_vaf, np_BQ, step, **kwargs)     # 첫 3개는 lenient로, 그 다음은 strict로 거른다

        if step.makeone_index == []:   # if failed  (첫 1개는 웬만하면 봐주려고 하지만, 그래도 잘 안될 때)
            step.likelihood = step.likelihood_record[kwargs["STEP"]] = -9999999
            step.checkall_lenient = step.checkall_strict = False
            if kwargs["VERBOSE"] >= 1:
                print ("\t\t\tMstep.py : condition = {}\tcheckall({}) = False\tUnable to make 1 or appropriate superclone-clone structure ({})".format ( condition, condition, step.mixture.flatten()  ) )

        else:    # if succeedeed
            if kwargs["STEP"] <= (kwargs["COMPULSORY_NORMALIZATION"] ):     # 첫 3개는 강제 normalization 해주자
                for i in range(NUM_BLOCK):
                    sum = 0
                    for j in range(NUM_CLONE):
                        if j in step.makeone_index:            
                            sum = sum + step.mixture[i][j]           

                    for j in range(NUM_CLONE):
                        if j in step.makeone_index:      
                            step.mixture[i][j] = np.round(step.mixture[i][j] / sum, 2) if sum != 0 else 0

            if kwargs["VERBOSE"] >= 1:
                if kwargs["STEP"] <= (kwargs["COMPULSORY_NORMALIZATION"] - 1):     # Normalization 꼭 해줄 필요가 없다.  정말 앞부분일 때에만 해준다
                    print ("\t\t\tMstep.py : condition = {}\tcheckall (strict) = {}\tcheckall (lenient) = {}\tnormalized ({})".format(condition, step.checkall_strict,  step.checkall_lenient, step.mixture.flatten () )) 
                else:
                    print ("\t\t\tMstep.py : condition = {}\tcheckall (strict) = {}\tcheckall (lenient) = {}\t({})".format(condition, step.checkall_strict,  step.checkall_lenient, step.mixture.flatten () )) 
                
        
            step.membership_p_normalize = np.zeros((NUM_MUTATION, step.membership_p.shape[1]), dtype="float64")
            for k in range(NUM_MUTATION):
                if k in step.fp_member_index:
                    step.membership_p_normalize[k] = np.zeros(step.membership_p_normalize.shape[1], dtype="float64")  # Set  1 (FP_index) 0 0 0 0 0    
                    step.membership_p_normalize[k][step.fp_index] = 1
                else:
                    step.membership_p_normalize[k] = np.round(np.power(10, step.membership_p[k]) / np.power(10, step.membership_p[k]).sum(axis=0, keepdims=1), 2)     # 분율을 따진다
                    if step.fp_index != -1: 
                        step.membership_p_normalize[k][step.fp_index] = 0

        sum_mixture = np.sum (step.mixture, axis = 1)

        if  kwargs["VISUALIZATION"] == True:
            if (kwargs["NUM_BLOCK"] == 1):
                visualizationeachstep.drawfigure_1d_soft(step, np_vaf, kwargs["CLEMENT_DIR"] + "/trial/clone" + str(kwargs["NUM_CLONE_ITER"]) + "." + str(kwargs["TRIAL"]) + "-" + str(kwargs["STEP_TOTAL"]) + "(soft)." + kwargs["IMAGE_FORMAT"], **kwargs)
            if (kwargs["NUM_BLOCK"] == 2):
                visualizationeachstep.drawfigure_2d_soft(step, np_vaf, kwargs["CLEMENT_DIR"] + "/trial/clone" + str(kwargs["NUM_CLONE_ITER"]) + "." + str(kwargs["TRIAL"]) + "-" + str(kwargs["STEP_TOTAL"]) + "(soft)." + kwargs["IMAGE_FORMAT"], **kwargs)
            if (kwargs["NUM_BLOCK"] == 3):
                #visualizationeachstep.drawfigure_3d(step, np_vaf, kwargs["CLEMENT_DIR"] + "/trial/clone" + str(kwargs["NUM_CLONE_ITER"]) + "." + str(kwargs["TRIAL"]) + "-" + str(kwargs["STEP_TOTAL"]) + "(soft)." + kwargs["IMAGE_FORMAT"], sum_mixture, **kwargs)
                visualizationeachstep.drawfigure_3d_SVD(step, np_vaf, kwargs["CLEMENT_DIR"] + "/trial/clone" + str(kwargs["NUM_CLONE_ITER"]) + "." + str(kwargs["TRIAL"]) + "-" + str(kwargs["STEP_TOTAL"]) + "(soft)." + kwargs["IMAGE_FORMAT"], sum_mixture, **kwargs)

    #############################################################################

    return step
