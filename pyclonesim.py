import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import palettable
import scoring
import visualizationpair

#def visualization(np_vaf, OUTPUT_FILENAME, b):



def main (INPUT_NPVAF, INPUT_FILENAME, OUTPUT_FILENAME,  mixture_answer, membership_answer, samplename_dict_input, samplename_dict_input_rev):   
    df = pd.read_csv (INPUT_FILENAME, sep = "\t")


    # pyclone의 clustering 정보만 가져온다
    b1 = df[df["sample_id"] == "block0"].iloc[:, 0:3]
    b2 = df[df["sample_id"] == "block1"].iloc[:, 0:3]
    b = pd.merge (b1, b2, left_on = "mutation_id", right_on = "mutation_id")

    # 각 mutation_id의 vaf 정보를 가져오기 위함
    np_vaf = pd.read_csv(INPUT_NPVAF, sep = "\t")
    np_vaf.rename(columns = {"Unnamed: 0":"mutation_id"}, inplace = True) 
    np_vaf.head()

    b = pd.merge (np_vaf, b,  left_on = "mutation_id", right_on = "mutation_id")

    # mixture_pyclone, membership_pyclone 정해주기
    membership_pyclone =  list(b["cluster_id_x"])
    mixture_pyclone = np.zeros ((mixture_answer.shape[0], len(set(membership_pyclone))), dtype = "float")
    for j in set(membership_pyclone):
        x_mean = np.mean (b[b["cluster_id_x"] == j]["block0"])
        y_mean = np.mean (b[b["cluster_id_x"] == j]["block1"])
        mixture_pyclone[:, j ] = [x_mean * 2, y_mean * 2]
        


    # #visualization(np_vaf, OUTPUT_FILENAME, b)
    # vivid_10 = palettable.cartocolors.qualitative.Vivid_10.mpl_colors
    # blue_seq = palettable.scientific.sequential.Devon_18.mpl_colors
    # colorlist = {key:value for key, value in enumerate(vivid_10)}
    

    # plt.figure (figsize = (7,7))
    # plt.suptitle ("PyClone : (n = {0})".format(np_vaf.shape[0]), fontsize = 30)
    # plt.axis ([0,  np.max(np.array(np_vaf.iloc[:,1:-1])) * 1.1,  0,  np.max(np.array(np_vaf.iloc[:,1:-1])) * 1.1])
    # plt.xlabel ("VAF of the sample 1")
    # plt.ylabel ("VAF of the sample 2")

    # sns.set_style("darkgrid", {"grid.color": ".6", "grid.linestyle": ":"})
    # sns.scatterplot(data = b, x = 'block0', y = 'block1', hue = "cluster_id_x", palette = colorlist, s = 50)
    # plt.savefig(OUTPUT_FILENAME)


    # 채점하기
    import scoring
    membership_pyclone =  list(b["cluster_id_x"])
    # maxmaxmax_membership =  list(b["cluster_id_x"])
    # membership_answer_old = list(b["membership_answer"])
    # membership_answer_new = scoring.Sorting(membership_answer_old )
    # membership_answer_max, score, sample_dict_rev = scoring.Scoring \
    #     (membership_answer_old, membership_answer_new, maxmaxmax_membership)

    score_df, score = \
        scoring.mixturebased(mixture_answer, mixture_pyclone, membership_answer, membership_pyclone, samplename_dict_input, samplename_dict_input_rev, "No")

    mixture_pyclone = np.round(mixture_pyclone, 2)
    return score_df, score, membership_pyclone, mixture_pyclone          # 100점 만점 score 반환




if __name__ == "__main__":
    import argparse
    from collections import Counter

    print ("main source가 수행됨")

    parser = argparse.ArgumentParser(description='Here is usage direction.')
    parser.add_argument('--INPUT_NPVAF', default="")
    parser.add_argument('--OUTPUT_FILENAME', default="")

    args = parser.parse_args()
    INPUT_NPVAF = args.INPUT_NPVAF
    INPUT_FILENAME = args.INPUT_FILENAME
    OUTPUT_FILENAME = args.OUTPUT_FILENAME

    main (INPUT_NPVAF, INPUT_FILENAME, OUTPUT_FILENAME)