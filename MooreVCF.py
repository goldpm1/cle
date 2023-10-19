import pandas as pd
import numpy as np
import glob, os, subprocess

kwargs = {}

kwargs ["INPUT_TSV"] = "/data/project/Alzheimer/CLEMENT/01.INPUT_TSV/2.CellData/CellData_1D/250x/M1-2_input.txt"
kwargs ["INPUT_MEMBERSHIP"] = "/data/project/Alzheimer/CLEMENT/03.combinedoutput/2.CellData/CellData_1D/n1000_250x/parent_0/fp_0.0/axis_-1/M1-2/0/0.input_membership_numerical.txt"
kwargs ["INPUT_MIXTURE"] = "/data/project/Alzheimer/CLEMENT/03.combinedoutput/2.CellData/CellData_1D/n1000_250x/parent_0/fp_0.0/axis_-1/M1-2/0/0.input_mixture.txt"

kwargs ["CLEMENT_MEMBERSHIP"] = "/data/project/Alzheimer/CLEMENT/03.combinedoutput/2.CellData/CellData_1D/n1000_250x/parent_0/fp_0.0/axis_-1/M1-2/0/result/CLEMENT_decision.membership.txt"
kwargs ["CLEMENT_MIXTURE"] = "/data/project/Alzheimer/CLEMENT/03.combinedoutput/2.CellData/CellData_1D/n1000_250x/parent_0/fp_0.0/axis_-1/M1-2/0/result/CLEMENT_decision.mixture.txt"

kwargs ["PYCLONEVI_MEMBERSHIP"] = "/data/project/Alzheimer/CLEMENT/03.combinedoutput/2.CellData/CellData_1D/n1000_250x/parent_0/fp_0.0/axis_-1/M1-2/0/result/pyclonevi.membership.txt"
kwargs ["PYCLONEVI_MIXTURE"] = "/data/project/Alzheimer/CLEMENT/03.combinedoutput/2.CellData/CellData_1D/n1000_250x/parent_0/fp_0.0/axis_-1/M1-2/0/result/pyclonevi.mixture.txt"

Moore_VCF = pd.read_csv("/data/project/Alzheimer/CLEMENT/resource/paper/whole_info.polyphen.sift.alpha.txt", sep = "\t")

# column을 소문자로 바꾸기
Moore_VCF.rename(columns=lambda x: x.lower(), inplace=True)

# lexicographically 정렬하기
custom_order = ["1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "X", "Y"]
Moore_VCF['chr_categorical'] = pd.Categorical(Moore_VCF['chr'], categories=custom_order, ordered=True)
Moore_VCF = Moore_VCF.sort_values(by= ['chr_categorical', 'start']).reset_index(drop=True).drop ("chr_categorical", axis = 1)

# 앞에 "chr" 붙여주기
Moore_VCF["pos"]  = "chr" + Moore_VCF["chr"] + "_" + Moore_VCF["start"].astype("str")
Moore_VCF["refalt"] =  Moore_VCF["ref"] + "/" + Moore_VCF["alt"]
Moore_VCF["gene_single"] = Moore_VCF['gene'].str.split(',')
print (Moore_VCF.shape)


from tqdm import tqdm

for k in tqdm ( range ( Moore_VCF.shape[0] ) ) :
    check = False
    for g in Moore_VCF.iloc[k]["gene_single"]:
        if ("OR" not in g)  & ("LOC" not in g) & ("NONE" not in g):
            check = True
            break
    if check == False:
        ii = [index for index, item in enumerate( Moore_VCF.iloc[k]["gene_single"] ) if 'LOC' in item]
        if len(ii) >= 1:
            g = Moore_VCF.iloc[k]["gene_single"][ii [0] ]
        else:
            g = "-"

    Moore_VCF.loc[k, "gene_single"] = g

    if k % 50000 == 0:
        print (k)
    #print (g , Moore_VCF.loc[k, "gene_single"])


Moore_VCF.to_csv ("/data/project/Alzheimer/CLEMENT/resource/paper/whole_info_singlegene.txt", sep = "\t", index = False)
Moore_VCF.loc[:, ["chr", "start", "end", "refalt", "depth_count", "alt_count", "vaf", "sample", "donorid", "tissue", "func", "gene_single"]].to_csv ("/data/project/Alzheimer/CLEMENT/resource/paper/whole_info.vepinput.txt", sep = "\t", index = False)