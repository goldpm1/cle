# Interpolation 후 아래에 있는 것들은 
def outlier (df, np_vaf, membership, clone_no):
    df_new = []
    df_new_index= []
    for i in [i for i in range(len(membership)) if membership[i] == clone_no] :       # membership 번호가 맨 끝벙니면 outlier로 처리
        df_new.append(df[i])
        df_new_index.append(i)
    return df_new, df_new_index, np_vaf[df_new_index]


def interpolatedoutlier (df, np_vaf, membership, clone_no):
    df_new = []
    df_new_index= []
    for i in [i for i in range(len(membership)) if membership[i] == clone_no] :       # membership 번호가 맨 끝벙니면 outlier로 처리
        df_new.append(df[i])
        df_new_index.append(i)
    return df_new, df_new_index, np_vaf[df_new_index]