import os



# Remove SimData
for data in ["1.SimData/SimData_1D", "1.SimData/SimData_2D", "1.SimData/SimData_3D"]:
    for tool in  ["CLEMENT", "pyclone-vi", "quantumclone", "sciclone", "SIMPLE_KMEANS"]:
        os.system ("rm -rf " + "/data/project/Alzheimer/YSscript/cle/data/" + tool + "/" + data )
    os.system ("rm -rf /data/project/Alzheimer/CLEMENT/02.npvaf/" + data )
    os.system ("rm -rf /data/project/Alzheimer/CLEMENT/03.combinedoutput/" + data )
    os.system ("rm -rf /data/project/Alzheimer/YSscript/cle/log/" + data )
    print ("{} removed".format(data))


# Remove CellData
# for data in ["2.CellData/CellData_2D"]:
#     for tool in  ["CLEMENT", "pyclone-vi", "quantumclone", "sciclone", "SIMPLE_KMEANS"]:
#         os.system ("rm -rf " + "/data/project/Alzheimer/YSscript/cle/data/" + tool + "/" + data )
#     os.system ("rm -rf /data/project/Alzheimer/CLEMENT/02.npvaf/" + data )
#     os.system ("rm -rf /data/project/Alzheimer/CLEMENT/03.combinedoutput/" + data )
#     os.system ("rm -rf /data/project/Alzheimer/YSscript/cle/log/" + data )
#     print ("{} removed".format(data))


# Remove BioData
#for data in ["3.BioData/Brunner_2D", "3.BioData/Moore_1D", "3.BioData/Moore_1D_AG", "3.BioData/Moore_2D", "3.BioData/Moore_2D_AG"]:
# for data in ["3.BioData/Brunner_2D"]:
#     for tool in  ["CLEMENT", "pyclone-vi", "quantumclone", "sciclone", "SIMPLE_KMEANS"]:
#         os.system ("rm -rf " + "/data/project/Alzheimer/YSscript/cle/data/" + tool + "/" + data )
#     os.system ("rm -rf /data/project/Alzheimer/CLEMENT/02.npvaf/" + data )
#     os.system ("rm -rf /data/project/Alzheimer/CLEMENT/03.combinedoutput/" + data )
#     os.system ("rm -rf /data/project/Alzheimer/YSscript/cle/log/" + data )
#     print ("{} removed".format(data))

