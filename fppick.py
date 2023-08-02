import numpy as np
import EMsoft

def main (mixture, membership, membership_p, df, **kwargs):
    VERBOSE = kwargs["VERBOSE"]
    if (VERBOSE in ["True", "T"]) | (int(str(VERBOSE)) >= 1):
        print ("조정 전\n{}\n{}\n\n".format(np.unique(membership , return_counts = True ), mixture))

    # i차원에서 몇 개가 non-zero로 나와야 false positive로 해석할 수 있을지는 불분명하다
    NUM_MUTATION = len(membership)
    NUM_CLONE =  mixture.shape[1]
    NUM_BLOCK = kwargs["NUM_BLOCK"]

    #1. clone을 space, axis로 구분하기  +  mixture, membership_p를 재배치하기 (뽑아오기)
    space_index, axis_index =[], []
    for j in range(NUM_CLONE):
        zero_axis_count = np.count_nonzero(mixture[:,j] < 0.01)
        if zero_axis_count >= int(kwargs["NUM_BLOCK"])/2:       # axis축 colne만 뽑음
            axis_index.append(j)
        else:
            space_index.append(j)
    mixture_space = mixture[:, space_index]
    membership_p_space = membership_p[:, space_index]
    membership_new = np.zeros(NUM_MUTATION, dtype = "int")

    #2. axis clone을 OUTLIER로 밀어버리기
    OUTLIER_NO = NUM_CLONE - 1
    for k in range (NUM_MUTATION):
        if membership[k] in axis_index:
            membership_new[k] = OUTLIER_NO
        else:
            membership_new[k] = membership[k]

    #3. membership을 재배치하기

    newname_dict = {space_index[i]:i for i in range(0,len(space_index))}       # [0,1,3,4] -> [0,1,2,3]
    for k in range(NUM_MUTATION):
        if membership_new[k] in space_index:        # axis index는 이미 바꿔줬으니 공란이다
            membership_new[k] = newname_dict[ membership_new[k] ]      # 번호를 앞에서부터 순서대로 바꿔주기

    #4. mixture를 재계산하기
    OUTLIER_NO = mixture_space.shape[1] - 1
    if (VERBOSE in ["True", "T"]) | (int(str(VERBOSE)) >= 0):
        print ("axis index : {}\nNEW_OUTILER_NO : {}, NUM_CLONE : {}".format(axis_index, OUTLIER_NO, OUTLIER_NO + 1))

    mixture_new = EMsoft.Mstep(membership_new, membership_p_space, df, mixture_space.shape[1], NUM_BLOCK, NUM_MUTATION,  "Last")           # 새 mixture구하기 (Outlier 제외)

    for i in range (kwargs["NUM_BLOCK"]):
        sum_depth, sum_alt = 0 , 0
        for k in range(NUM_MUTATION):
            if membership_new[k] == OUTLIER_NO:    # outlier 찾기
                sum_depth = sum_depth + df[k][i]["depth"]
                sum_alt = sum_alt + df[k][i]["alt"]

        if sum_depth == 0:      
            mixture_new[i][-1] = 0
        else:                                # outlier clone만을 생각한 이상적인 분율을 일단 assign
            mixture_new[i][-1] = round((sum_alt * 2) / sum_depth, 2)

    NUM_CLONE = mixture_new.shape[1]

    #5. membership_p_normalize를 재계산하기
    membership_p_normalize_space = np.zeros((NUM_MUTATION, membership_p_space.shape[1]), dtype="float64")
    for k in range(NUM_MUTATION):
        membership_p_normalize_space[k] = np.round(np.power(10, membership_p_space[k])/np.power(10, membership_p_space[k]).sum(axis=0, keepdims=1), 2)   # 로그를 취해으나 다시 지수를 취해준다

    if (VERBOSE in ["True", "T"]) | (int(str(VERBOSE)) >= 1):
        print ("\n조정 후\n{}\n{}".format(np.unique(membership_new , return_counts = True ), mixture_new ))

    return mixture_new, membership_new, membership_p_space, membership_p_normalize_space, axis_index, NUM_CLONE