import os



# Remove SimData
for data in ["1.SimData/SimData_1D", "1.SimData/SimData_2D", "1.SimData/SimData_3D"]:
    for tool in  ["CLEMENT", "pyclone-vi", "quantumclone", "sciclone", "SIMPLE_KMEANS"]:
        os.system ("rm -rf " + "/data/project/Alzheimer/YSscript/cle/data/" + tool + "/" + data  )
    os.system ("rm -rf /data/project/Alzheimer/CLEMENT/01.INPUT_TSV/" + data  )
    os.system ("rm -rf /data/project/Alzheimer/CLEMENT/02.npvaf/" + data  )
    os.system ("rm -rf /data/project/Alzheimer/CLEMENT/03.combinedoutput/" + data )
    os.system ("rm -rf /data/project/Alzheimer/YSscript/cle/log/" + data  )
    print ("{} removed".format(data))


# Remove CellData
# for data in ["2.CellData/CellData_3D"]:
#     for tool in  ["CLEMENT", "pyclone-vi", "quantumclone", "sciclone", "SIMPLE_KMEANS"]:
#         os.system ("rm -rf " + "/data/project/Alzheimer/YSscript/cle/data/" + tool + "/" + data  + "/*")
#     os.system ("rm -rf /data/project/Alzheimer/CLEMENT/02.npvaf/" + data + "/*" )
#     os.system ("rm -rf /data/project/Alzheimer/CLEMENT/03.combinedoutput/" + data + "/*" )
#     os.system ("rm -rf /data/project/Alzheimer/YSscript/cle/log/" + data + "/*" )
#     print ("{} removed".format(data))

# for data in ["2.CellData/CellData_2D/n100_250x", "2.CellData/CellData_2D/n500_250x", "2.CellData/CellData_2D/n500_30x", "2.CellData/CellData_2D/n1000_30x", "2.CellData/CellData_2D/n1000_125x", "2.CellData/CellData_2D/n1000_250x"]:
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

