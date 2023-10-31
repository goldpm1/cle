def GenerateMatrix (Moore_VCF, Clone_no, **kwargs):
    import pandas as pd

    Moore_VCF = Moore_VCF.reset_index(drop = True)

    SigProfMatrix = pd.DataFrame (index = range ( len (Moore_VCF)),  columns = ["Project", "Sample", "ID", "Genome", "mut_type", "chrom",	"pos_start","pos_end",	"ref",	"alt"	,"Type"]) 
    SigProfMatrix["Project"] = "BioData"
    SigProfMatrix["Sample"] = str(Clone_no)
    SigProfMatrix["ID"] = Moore_VCF["pos"]
    SigProfMatrix["Genome"] = "GRCh37"
    SigProfMatrix["mut_type"] = "SNP"
    SigProfMatrix["chrom"] = Moore_VCF["chr"]
    SigProfMatrix["pos_start"] = Moore_VCF["start"]
    SigProfMatrix["pos_end"] = Moore_VCF["end"]
    SigProfMatrix["ref"] = Moore_VCF["ref"]
    SigProfMatrix["alt"] = Moore_VCF["alt"]
    SigProfMatrix["Type"] = "SOMATIC"

    SigProfMatrix = SigProfMatrix.drop_duplicates()

    SigProfMatrix.to_csv( kwargs["OUTPUT_PATH"], index=False, header = True,  sep="\t")
    print ("clone {} : {}".format (Clone_no, SigProfMatrix.shape) )



def extract (Moore_VCF, **kwargs):
    import numpy as np
    import pandas as pd

    pos_whole = []

    # Clone number 대로  mutation_id를 뽑아주기
    with open( kwargs ["DECISION_MEMBERSHIP_PATH"], "r") as input_file:
        membership = np.array( [int(line.strip()) for line in input_file] )

    npvaf = pd.read_csv ( kwargs ["NPVAF_PATH"], sep = "\t")

    for Clone_no in sorted( np.unique ( membership) ) :
        npvaf_Clone_no =  npvaf.iloc [np.where(membership == Clone_no)[0] ]
        pos_Clone_no =  list ( npvaf_Clone_no.iloc[:, 0] )
        pos_whole += pos_Clone_no


        print (pos_Clone_no)
        Moore_VCF_Clone_no = Moore_VCF [Moore_VCF['pos'].isin( pos_Clone_no ) == True]
        print (Moore_VCF_Clone_no)
        
        kwargs["OUTPUT_PATH"] = kwargs["OUTPUT_DIR"] + "/clone{}.txt".format(Clone_no)
        GenerateMatrix (Moore_VCF_Clone_no, Clone_no, **kwargs)

    # # whole
    # kwargs["OUTPUT_PATH"] = kwargs["OUTPUT_DIR"] + "/whole.txt"
    # GenerateMatrix ( Moore_VCF [Moore_VCF['pos'].isin( list ( npvaf.iloc[:, 0] ) ) == True] , "whole", **kwargs)
    # print ("whole : {}".format ( Moore_VCF [Moore_VCF['pos'].isin( list ( npvaf.iloc[:, 0] ) ) == True].shape) )

    # whole
    kwargs["OUTPUT_PATH"] = kwargs["OUTPUT_DIR"] + "/whole.txt"
    Moore_VCF_whole = Moore_VCF[Moore_VCF['pos'].isin( pos_whole ) == True]
    GenerateMatrix (Moore_VCF_whole, "whole", **kwargs)





if __name__ == "__main__":
    import pandas as pd
    import argparse

    parser = argparse.ArgumentParser(description='The below is usage direction.')
    parser.add_argument("--DECISION_MEMBERSHIP_PATH", type = str, default = "")
    parser.add_argument('--NPVAF_PATH', type = str, default = "")
    parser.add_argument('--DONOR', type = str, default = "")
    parser.add_argument('--TISSUE', type = str, default = "")
    parser.add_argument('--OUTPUT_DIR', type = str, default = "")

    kwargs = {}
    args = parser.parse_args()

    kwargs["DECISION_MEMBERSHIP_PATH"] = args.DECISION_MEMBERSHIP_PATH
    kwargs["NPVAF_PATH"] = args.NPVAF_PATH
    kwargs["DONOR"] = args.DONOR
    kwargs["TISSUE"] = args.TISSUE
    kwargs["OUTPUT_DIR"] = args.OUTPUT_DIR

    print ("DONOR = {}\tTISSUE = {}\tOUTPUT_DIR = {}".format (kwargs["DONOR"], kwargs["TISSUE"], kwargs["OUTPUT_DIR"]))

    Moore_VCF = pd.read_csv("/data/project/Alzheimer/CLEMENT/resource/paper/old/whole_info.txt", sep = "\t")

    #"/data/project/Alzheimer/CLEMENT/resource/paper/whole_info.txt"
    

    # column을 소문자로 바꾸기
    Moore_VCF.rename(columns=lambda x: x.lower(), inplace=True)

    # lexicographically 정렬하기
    custom_order = ["1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "X", "Y"]
    Moore_VCF['chr_categorical'] = pd.Categorical(Moore_VCF['chr'], categories=custom_order, ordered=True)
    Moore_VCF = Moore_VCF.sort_values(by= ['chr_categorical', 'start']).reset_index(drop=True).drop ("chr_categorical", axis = 1)

    # 앞에 "chr" 붙여주기
    Moore_VCF["pos"]  = "chr" + Moore_VCF["chr"] + "_" + Moore_VCF["start"].astype("str")


    extract (Moore_VCF, **kwargs)