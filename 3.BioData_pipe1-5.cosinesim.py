if __name__ == "__main__":
    import pandas as pd
    from sklearn.metrics.pairwise import cosine_similarity
    import numpy as np
    import argparse

    parser = argparse.ArgumentParser(description='The below is usage direction.')
    parser.add_argument("--SIGPROFILER_PATH", type = str, default = "")
    parser.add_argument("--OUTPUT_PATH", type = str, default = "")

    kwargs = {}
    args = parser.parse_args()

    kwargs["SIGPROFILER_PATH"] = args.SIGPROFILER_PATH
    kwargs["OUTPUT_PATH"] = args.OUTPUT_PATH

    COSMIC_PATH = "/home/goldpm1/resources/COSMIC_SBS/COSMIC_v3.3.1_SBS_GRCh38.txt"
    #SIGPROFILER_PATH = "/data/project/Alzheimer/CLEMENT/03.combinedoutput/3.BioData/Moore_1D/bronchus_epithelium/PD28690-H7/SigProfiler/output/DeNovo/SBS96/All_Solutions/SBS96_2_Signatures/Signatures/SBS96_S2_Signatures.txt"

    COSMIC_df = pd.read_csv (COSMIC_PATH, sep = "\t", index_col=0)
    SIGPROFILER_df = pd.read_csv ( kwargs["SIGPROFILER_PATH"], sep = "\t", index_col= 0)
    SIGPROFILER_df.index.name = 'Type'

    print ( np.all(COSMIC_df.index == SIGPROFILER_df.index))   # A[C>A]A 등으 ㅣ순서가 정확히 동일한지 확인


    # Compute cosine similarity
    cos_sim = cosine_similarity( COSMIC_df.values.transpose(), SIGPROFILER_df.values.transpose() )

    # 가장 비슷한 (= cosine similarity가 가장 높은 SBS) SBS 찾기

    output_file = open ( kwargs["OUTPUT_PATH"], "w")

    for signature_i in range ( cos_sim.shape[1] ):
        lst = sorted(  cos_sim[:, signature_i], reverse=True)
        lst_arg = []
        for i in lst:
            lst_arg.append ( list(cos_sim[:, signature_i]).index (i) )

        max = np.max ( cos_sim[:, signature_i] )
        max_arg = np.argmax ( cos_sim[:, signature_i] )
        #print ("S{}\tmax = {}\tmax_arg = {}".format(signature_i + 1, round (max, 2), COSMIC_df.columns[ max_arg ] )) 

        print ("{}\tlst[0:5] = {}\tlst_arg[0:5] = {}".format ( SIGPROFILER_df.columns [ signature_i ], np.round (lst[0:5], 2) ,  COSMIC_df.columns [ [i for i in  lst_arg [0:5]] ]  ) , file = output_file )

    output_file.close()