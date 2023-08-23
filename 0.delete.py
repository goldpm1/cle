import os

# Remove CellData

for data in ["2.CellData/CellData_1D"]:
    for tool in  ["CLEMENT", "pyclone-vi", "quantumclone", "sciclone", "SIMPLE_KMEANS"]:
        os.system ("rm -rf " + "/data/project/Alzheimer/YSscript/cle/data/" + tool + "/" + data )
    os.system ("rm -rf /data/project/Alzheimer/CLEMENT/02.npvaf/" + data )
    os.system ("rm -rf /data/project/Alzheimer/CLEMENT/03.combinedoutput/" + data )
    os.system ("rm -rf /data/project/Alzheimer/YSscript/cle/log/" + data )
    print ("{} removed".format(data))


# Remove BioData
# for data in ["3.BioData/Moore_1D", "3.BioData/Moore_1D_AG"]:
#     for tool in  ["CLEMENT", "pyclone-vi", "quantumclone", "sciclone", "SIMPLE_KMEANS"]:
#         os.system ("rm -rf " + "/data/project/Alzheimer/YSscript/cle/data/" + tool + "/" + data )
#     os.system ("rm -rf /data/project/Alzheimer/CLEMENT/02.npvaf/" + data )
#     os.system ("rm -rf /data/project/Alzheimer/CLEMENT/03.combinedoutput/" + data )
#     os.system ("rm -rf /data/project/Alzheimer/YSscript/cle/log/" + data )
#     print ("{} removed".format(data))

