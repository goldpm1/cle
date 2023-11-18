class ResultClass:
    def __init__(self, toollist, **kwargs):
        self.score_record = np.zeros ( ( len (toollist), kwargs["BENCHMARK_END"]  + 1  ), dtype = "int" )
        self.Yindex_record = np.zeros ( (len (toollist)  , kwargs["BENCHMARK_END"]  + 1 ), dtype = "float" )
        self.ARI_record = np.zeros ( (len (toollist)  , kwargs["BENCHMARK_END"]  + 1 ), dtype = "float" )
        self.NUM_CLONE_record = np.zeros ( (len (toollist) , kwargs["BENCHMARK_END"]  + 1 ), dtype = "int" )
        self.runningtime_record = np.zeros ( (len (toollist) , kwargs["BENCHMARK_END"]  + 1 ), dtype = "float" )
        self.f1score_record = np.zeros ( (len (toollist) , kwargs["BENCHMARK_END"]  + 1 ), dtype = "float" )

    def acc (self, I, J, score, Yindex, ARI, NUM_CLONE_answer, runningtime, f1score):
        self.score_record[I , J] = score
        self.Yindex_record[I , J] = Yindex
        self.ARI_record[I , J] = ARI
        self.NUM_CLONE_record[I , J] = NUM_CLONE_answer
        self.runningtime_record[I , J] = runningtime
        self.f1score_record[I , J] = f1score

class FinalResult:
    def __init__(self, toollist, DIR_num ):
        self.score_record = np.zeros ( ( len (toollist), DIR_num  ), dtype = "int" )
        self.ARI_record = np.zeros ( (len (toollist)  , DIR_num ), dtype = "float" )
        self.RMSE_record = np.zeros ( (len (toollist)  , DIR_num ), dtype = "float" )

    def acc (self, I, J, median_score, median_ARI, RMSE):   # I : toollist J : sample
        self.score_record[I , J] = median_score
        self.ARI_record[I , J] = median_ARI
        self.RMSE_record[I , J] = RMSE

########################################################################################################################################################################


def drawfigure_MS (result, toollist, toollist_concise, ax, ax_row, ax_col, **kwargs):
    import palettable
    import matplotlib
    import numpy as np
    import pandas as pd
    import seaborn as sns
    from scipy.stats import ttest_ind
    import itertools

    tabl = palettable.tableau.Tableau_20.mpl_colors
    safe7 = palettable.cartocolors.qualitative.Vivid_10.mpl_colors
    Gr_10 = palettable.scientific.sequential.GrayC_20.mpl_colors

    colorlist = [i for i in tabl]
    colorlist = ["royalblue", "firebrick", "forestgreen", "darkorange", Gr_10[6], Gr_10[5], Gr_10[4]]
    sns.set_style("white")
    #sns.set_palette("tab10")
    sns.set_palette(sns.color_palette(colorlist))

    matplotlib.rcParams["font.family"] = 'arial'

    # Seaborn을 위해 df를 만들기
    df = pd.DataFrame (columns = ["tool", "score", "Y-index", "ARI", "NUM_CLONE_answer", "f1score", "runningtime"] )
    matrix = []
    for j, tool in enumerate( toollist ):
        for k in range (0, len(result.score_record[0]) ):
            matrix.append ( [tool, result.score_record[j][k], result.Yindex_record[j][k], result.ARI_record[j][k], result.NUM_CLONE_record[j][k], result.f1score_record[j][k], result.runningtime_record[j][k]] )

    df = (pd.DataFrame.from_records (matrix, columns = df.columns))
    print (df)


    sns.boxplot (data = df, x = "tool", y = "score",  ax = ax[ax_row][ax_col], linewidth = 1)
    
    ax[ax_row][ax_col].set_title ( kwargs["SAMPLENAME"], y = 0.98, fontweight = "semibold", fontsize = 15, loc = "center")

    ax[ax_row][ax_col].set_xlabel("")
    ax[ax_row][ax_col].xaxis.label.set_fontsize(11)
    ax[ax_row][ax_col].yaxis.label.set_fontsize(13)
    ax[ax_row][ax_col].tick_params(axis = 'x', rotation = 35)
    ax[ax_row][ax_col].set_xticks ( list ( range ( len(toollist)  ) ) )   # 숫자로 tick을 박아놓고
    ax[ax_row][ax_col].set_xticklabels( toollist_concise )  # 표시는 toollist_concise로

    xmin, xmax, ymin, ymax = ax[ax_row][ax_col].axis()
    ymax = (ymax - ymin) * 0.2 + ymax
    ax[ax_row][ax_col].axis ( [ xmin, xmax, ymin, ymax ]) 


########################################################################################################################################################################


def drawfigure_EC (result, toollist, toollist_concise, ax, ax_row, ax_col, NUM_CLONE_ans, **kwargs):
    import palettable, math
    import matplotlib
    import numpy as np
    import pandas as pd
    import seaborn as sns
    from scipy.stats import ttest_ind
    import itertools


    tabl = palettable.tableau.Tableau_20.mpl_colors
    safe7 = palettable.cartocolors.qualitative.Vivid_10.mpl_colors
    Gr_10 = palettable.scientific.sequential.GrayC_20.mpl_colors

    colorlist = [i for i in tabl]
    sns.set_style("white")
    #sns.set_palette("tab10")
    sns.set_palette(sns.color_palette(colorlist))

    matplotlib.rcParams["font.family"] = 'arial'


    # Seaborn을 위해 df를 만들기
    df = pd.DataFrame (columns = ["tool", "score", "Y-index", "ARI", "NUM_CLONE_answer", "f1score", "runningtime"] )
    matrix = []
    for j, tool in enumerate( toollist ):
        for k in range (0, len(result.score_record[0]) ):
            matrix.append ( [tool, result.score_record[j][k], result.Yindex_record[j][k], result.ARI_record[j][k], result.NUM_CLONE_record[j][k], result.f1score_record[j][k], result.runningtime_record[j][k]] )

    df = (pd.DataFrame.from_records (matrix, columns = df.columns))
    #print (df)


    #sns.swarmplot (data = df, x = "tool", y = "NUM_CLONE_answer", hue = "tool", size =5 ,  legend = False, ax = ax[ax_row][ax_col] )
    for j, tool in enumerate( toollist ):
        value_count = df [df ["tool"] == tool ] ["NUM_CLONE_answer"].value_counts()
        #print ("{} -> value_count : {}".format (tool, value_count))
        value_count_dict = {}      # {2:8, 3:69, 4:148, 5:112, 6:26, 7:17}
        for i in value_count.index:
            value_count_dict [i] = value_count.loc[i]
            ax[ax_row][ax_col].scatter ( x = j, y = i, s = value_count_dict[i] * 80, color = tabl[j])

    
    ax[ax_row][ax_col].set_title ( kwargs["SAMPLENAME"], y = 0.98, fontweight = "semibold", fontsize = 15, loc = "center")

    ax[ax_row][ax_col].set_xlabel("")
    ax[ax_row][ax_col].xaxis.label.set_fontsize(11)
    ax[ax_row][ax_col].yaxis.label.set_fontsize(13)
    ax[ax_row][ax_col].tick_params(axis = 'x', rotation = 35)
    ax[ax_row][ax_col].set_xticks ( list ( range ( len(toollist)  ) ) )   # 숫자로 tick을 박아놓고
    ax[ax_row][ax_col].set_xticklabels( toollist_concise )  # 표시는 toollist_concise로

    xmin, xmax, ymin, ymax = ax[ax_row][ax_col].axis()
    ymax = ymax + (ymax - ymin) * 0.2
    ymin =  ymin - (ymax - ymin) * 0.2 
    ax[ax_row][ax_col].axis ( [ xmin, xmax, ymin, ymax ]) 

    ax[ax_row][ax_col].set_yticks ( np.arange ( math.ceil(ymin), ymax + 0.01, 1 ) )   # 1간격으로 y tick을 설정
    
    

    # 사각형 그리기 
    rect = matplotlib.patches.Rectangle((-0.5, NUM_CLONE_ans - 0.5),                                # 사각형 꼭지점의 시작위치
                                                    ax[ax_row][ax_col].get_xlim()[1] - ax[ax_row][ax_col].get_xlim()[0] + 0.5, 1,        # x 길이, y 길이
                                                    linewidth=0.5, edgecolor='red', facecolor='black', alpha=0.3)
    ax[ax_row][ax_col].add_patch(rect)



def drawfigure_Final (final, toollist, toollist_concise, **kwargs):
    import palettable, math
    import matplotlib
    import numpy as np
    import pandas as pd
    import seaborn as sns

    output_final = open (kwargs["OUTPUT_FINAL_TABLE"], "w")
    print ("\n")
    for i, tool in enumerate( toollist ):
        print ("{}\t{}\t{}\t{}\t{}\t{}\t{}".format( tool, round ( np.mean (final.score_record [i , : ]) , 2) ,  round ( np.std (final.score_record [i , : ]) , 2), round ( np.mean (final.RMSE_record [i , : ]) , 2) ,  round ( np.std (final.RMSE_record [i , : ]) , 2) , round ( np.mean (final.ARI_record [i , : ]) , 2) ,  round ( np.std (final.ARI_record [i , : ]) , 2)  )  )
        print ("{}\t{}\t{}\t{}\t{}\t{}\t{}".format( tool, round ( np.mean (final.score_record [i , : ]) , 2) ,  round ( np.std (final.score_record [i , : ]) , 2), round ( np.mean (final.RMSE_record [i , : ]) , 2) ,  round ( np.std (final.RMSE_record [i , : ]) , 2), round ( np.mean (final.ARI_record [i , : ]) , 2) ,  round ( np.std (final.ARI_record [i , : ]) , 2)   ), file = output_final)
    print ("\n")
    output_final.close()



    ###################Visualization ###################
    tabl = palettable.tableau.Tableau_20.mpl_colors
    safe7 = palettable.cartocolors.qualitative.Vivid_10.mpl_colors
    Gr_10 = palettable.scientific.sequential.GrayC_20.mpl_colors

    colorlist = [i for i in tabl]
    sns.set_style("white")
    #sns.set_palette("tab10")
    sns.set_palette(sns.color_palette(colorlist))

    matplotlib.rcParams["font.family"] = 'arial'

    figFinal, axFinal = matplotlib.pyplot.subplots(nrows = 1,  ncols = 3, figsize = (16, 6 ))
    figFinal.subplots_adjust (wspace = 0.2, hspace = 0.5, bottom = 0.15, top = 0.8, left = 0.05, right = 0.95)
    figFinal.suptitle( "{}".format(kwargs["CONDITIONNAME"]), fontsize = 30, fontweight = "bold")

    flierprops = {"marker": "x"}
    boxprops = { "edgecolor": "black",  "linewidth": 2, "alpha": 0.5  }
    medianprops = {"color": "blue"}
    meanprops = {"marker":"s","markerfacecolor":"white", "markeredgecolor":"black", "markersize" : 10}

    # axFinal[0] : Membership score
    final_long = pd.melt( pd.DataFrame( final.score_record.T, columns = toollist, index = list ( range (tt) ) ) , var_name='toollist', value_name='value')
    sns.boxplot (data=final_long, x = "toollist", y = "value", hue = "toollist", flierprops = flierprops,   boxprops = boxprops,  medianprops = medianprops,  showmeans=True, meanprops = meanprops,  dodge = False, ax = axFinal[0])
    sns.stripplot(data=final_long, x = "toollist", y = "value", color = "black", jitter=True,  ax = axFinal[0] )
    axFinal[0].set_ylabel( "Membership score" )

    # axFinal[1] : ARI
    final_long = pd.melt( pd.DataFrame( final.ARI_record.T, columns = toollist, index = list ( range (tt) ) ) , var_name='toollist', value_name='value')
    sns.boxplot(data=final_long, x = "toollist", y = "value", hue = "toollist", flierprops = flierprops,   boxprops = boxprops,  medianprops = medianprops,  showmeans=True, meanprops = meanprops, dodge = False, ax = axFinal[1])
    sns.stripplot(data=final_long, x = "toollist", y = "value", color = "black", jitter=True,  ax = axFinal[1] )
    axFinal[1].set_ylabel( "ARI" )

    # axFinal[2] : Membership score
    final_long = pd.melt( pd.DataFrame( final.RMSE_record.T, columns = toollist, index = list ( range (tt) ) ) , var_name='toollist', value_name='value')
    sns.violinplot(data=final_long, x = "toollist", y = "value", hue = "toollist", flierprops = flierprops,   boxprops = boxprops,  medianprops=medianprops,  showmedians=True,  showmeans=True, meanprops = meanprops,  dodge = False,  ax = axFinal[2])
    #sns.stripplot(data=final_long, x = "toollist", y = "value", color = "black", jitter=True,  ax = axFinal[2] )
    axFinal[2].set_ylabel( "RMSE" )


    for k in range (3):
        axFinal[k].set_xlabel("")
        axFinal[k].yaxis.label.set_fontsize(13)
        axFinal[k].tick_params(axis = 'x', rotation = 35)
        axFinal[k].set_xticks ( list ( range ( len(toollist)  ) ) )   # 숫자로 tick을 박아놓고
        axFinal[k].set_xticklabels( toollist_concise )  # 표시는 toollist_concise로
        axFinal[k].set_title (chr (k+65) + "." , fontsize = 20, loc = "left")    # A. B. C.
        axFinal[k].get_legend().remove()

    figFinal.savefig(kwargs["OUTPUT_FINAL_JPG"])




########################################################################################################################################################################

def Mingle_intra(SAMPLENAME, NUM_BLOCK):
    global MRS_set
    #MRS_set = [['M1-2', 'M1-5', 'M1-6', 'M1-7', 'M1-8'] ]
    MRS_set = [['M1-2', 'M1-4', 'M1-6', 'M1-8'],  ['M2-2', 'M2-4', 'M2-6', 'M2-8'] ]


    if NUM_BLOCK == 1:
        for i in range(len(MRS_set)):
            if (SAMPLENAME.split("_")[0] in MRS_set[i]): 
                return True
        # if 'M3' in SAMPLENAME.split("_")[0]:   # M1, M2를 다 돌려보자
        #     return False
        # return True
    elif NUM_BLOCK == 2:
        for i in range(len(MRS_set)):
            if (SAMPLENAME.split("_")[0] in MRS_set[i]) & (SAMPLENAME.split("_")[1] in MRS_set[i]):
                return True
    elif NUM_BLOCK == 3:
        for i in range(len(MRS_set)):
            if (SAMPLENAME.split("_")[0] in MRS_set[i]) & (SAMPLENAME.split("_")[1] in MRS_set[i]) & (SAMPLENAME.split("_")[2] in MRS_set[i]):
                return True

    return False

def Mingle_inter (SAMPLENAME, NUM_BLOCK):
    global MRS_set
    #MRS_set = [['M1-1', 'M1-4', 'M1-5', 'M1-6', 'M1-8'],  ['M2-2', 'M2-4', 'M2-6', 'M2-10', 'M2-12'], ['M3-1', 'M3-2', 'M3-4', 'M3-5', 'M3-11']]
    MRS_set = [['M1-2', 'M1-4', 'M1-6', 'M1-8'],  ['M2-2', 'M2-4', 'M2-6', 'M2-8'] ]

    if NUM_BLOCK == 1:
        for i in range(len(MRS_set)):
            if (SAMPLENAME.split("_")[0] in MRS_set[i]): 
                return True
        return False
    elif NUM_BLOCK == 2:
        if (SAMPLENAME.split("_")[0] in MRS_set[0]) &  (SAMPLENAME.split("_")[1] in MRS_set[1]):
            return True
        else:
            return False
    elif NUM_BLOCK == 3:
        if (SAMPLENAME.split("_")[0] in MRS_set[0]) &  (SAMPLENAME.split("_")[1] in MRS_set[0]) &  (SAMPLENAME.split("_")[2] in MRS_set[1]):
            return True
        elif (SAMPLENAME.split("_")[0] in MRS_set[0]) &  (SAMPLENAME.split("_")[1] in MRS_set[1]) &  (SAMPLENAME.split("_")[2] in MRS_set[1]):
            return True
        else:
            return False

########################################################################################################################################################################



import numpy as np 
import pandas as pd
import os, glob, math
import matplotlib.pyplot
import seaborn as sns
import palettable
from sklearn.metrics import mean_squared_error
import argparse 

parser = argparse.ArgumentParser(description='The below is usage direction.')
parser.add_argument('--INPUT_DIR', type = str, default = "")
parser.add_argument('--CONDITIONNAME', type = str, default = "")
parser.add_argument('--BENCHMARK_START', type = int, default = 0)
parser.add_argument('--BENCHMARK_END', type = int, default = 10)
parser.add_argument('--OUTPUT_MS_JPG', type = str, default = "")
parser.add_argument('--OUTPUT_EC_JPG', type = str, default = "")
parser.add_argument('--OUTPUT_FINAL_JPG', type = str, default = "")
parser.add_argument('--OUTPUT_FINAL_TABLE', type = str, default = "")

kwargs = {}
args = parser.parse_args()

kwargs["INPUT_DIR"] = args.INPUT_DIR
kwargs["CONDITIONNAME"] = args.CONDITIONNAME
kwargs["BENCHMARK_START"] = args.BENCHMARK_START
kwargs["BENCHMARK_END"] = args.BENCHMARK_END
kwargs["OUTPUT_MS_JPG"] = args.OUTPUT_MS_JPG
kwargs["OUTPUT_EC_JPG"] = args.OUTPUT_EC_JPG
kwargs["OUTPUT_FINAL_JPG"] = args.OUTPUT_FINAL_JPG
kwargs["OUTPUT_FINAL_TABLE"] = args.OUTPUT_FINAL_TABLE


from natsort import natsorted
INPUT_DIR_LIST = natsorted (  glob.glob(kwargs["INPUT_DIR"] + "/*") ) 


toollist = ["CLEMENT_decision", "pyclonevi", "sciclone", "quantumclone", "simpleK_elbow", "simpleK_silhouette", "simpleK_gap"]
toollist_concise = ["CLEMENT", "pyclonevi", "sciclone", "qc", "simpleK_elb", "simpleK_sil", "simpleK_gap*"]


INPUT_DIR_LIST_SELECTED = []
for k in INPUT_DIR_LIST:
    if os.path.isdir (k) == True:
        SAMPLENAME = k.split("/")[-1]
        NUM_BLOCK = k.split("/")[-1].count("_") + 1
        if (Mingle_inter(SAMPLENAME, NUM_BLOCK) == True) | (Mingle_intra(SAMPLENAME, NUM_BLOCK) == True):
        #if (Mingle_inter(SAMPLENAME, NUM_BLOCK) == True) :
            INPUT_DIR_LIST_SELECTED.append (k)

INPUT_DIR_LIST = natsorted (  INPUT_DIR_LIST_SELECTED ) 
        


####################################################################################################

num_directories = len (INPUT_DIR_LIST)

if num_directories >= 49:
    NUM_ROW = 7
elif num_directories >= 42:
    NUM_ROW = 6
elif num_directories >= 30:
    NUM_ROW = 5
elif num_directories > 15:
    NUM_ROW = 4
elif num_directories > 9:
    NUM_ROW = 3
else:
    NUM_ROW = 2
print ("INPUT_DIR = {}\tn = {}\n".format( kwargs["INPUT_DIR"], num_directories ))

NUM_COL = math.ceil ( num_directories / NUM_ROW)
figMS, axMS = matplotlib.pyplot.subplots(nrows = NUM_ROW,  ncols = NUM_COL, figsize = (18,NUM_ROW * 3.5))
figMS.subplots_adjust (wspace = 0.4, hspace = 0.5, bottom = 0.1, top = 0.9, left = 0.05, right = 0.95)
figMS.suptitle( "{}".format(kwargs["CONDITIONNAME"]), fontsize = 30, fontweight = "bold")

figEC, axEC = matplotlib.pyplot.subplots(nrows = NUM_ROW,  ncols = NUM_COL, figsize = (18,NUM_ROW * 3.5))
figEC.subplots_adjust (wspace = 0.4, hspace = 0.5, bottom = 0.1, top = 0.9, left = 0.05, right = 0.95)
figEC.suptitle( "{}".format(kwargs["CONDITIONNAME"]), fontsize = 30, fontweight = "bold")

####################################################################################################


tt = 0
final = FinalResult ( toollist, len (INPUT_DIR_LIST) )


for DIR_index, DIR in enumerate( INPUT_DIR_LIST) :
    if os.path.isfile(DIR) == True:
        continue

    if "clone" in DIR.split("/")[-1]:   # simData
        NUM_CLONE_ans = int ( DIR.split("/")[-1].split("_")[-1] )

    else:
        if "_" not in DIR.split("/")[-1]:    # 1D
            NUM_CLONE_ans = 3 if  "M1" in DIR.split("/")[-1] else 4
        elif DIR.split("/")[-1].count("_") == 1:  # 2D
            if DIR.split("/")[-1].count("M1") == 2:   # "M1-4_M1-6"
                NUM_CLONE_ans = 3 
            elif (DIR.split("/")[-1].count("M2") == 2) | (DIR.split("/")[-1].count("M3") == 2) :    # "M2-4_M2-8"
                NUM_CLONE_ans = 4
            else:   # "M1-4_M2-6" 
                NUM_CLONE_ans = 5
        elif DIR.split("/")[-1].count("_") == 2:  #3D
            if DIR.split("/")[-1].count("M1") == 3:   # "M1-4_M1-6_M1-8"
                NUM_CLONE_ans = 3
            elif (DIR.split("/")[-1].count("M2") == 3) | (DIR.split("/")[-1].count("M3") == 3) :    # "M2-4_M2-8"
                NUM_CLONE_ans = 4 
            elif (DIR.split("/")[-1].count("M1") == 1) & (DIR.split("/")[-1].count("M2") == 1) & (DIR.split("/")[-1].count("M3") == 1) :    # "M1-1_M2-1_M3_1"
                NUM_CLONE_ans = 6
            else:
                NUM_CLONE_ans = 5

        NUM_CLONE_ans += int (DIR.split("/")[-4].split("_")[-1])
    
    print ("SAMPLENAME = {}\tNUM_CLONE_ans = {}".format (  DIR.split("/")[-1], NUM_CLONE_ans) ) 

    result = ResultClass(toollist, **kwargs)

    for j in range (kwargs["BENCHMARK_START"], kwargs["BENCHMARK_END"] + 1):
        break_check = False
        for i, tool in enumerate( toollist ):
            score, Yindex, ARI, NUM_CLONE_answer, runningtime, f1score = 0, 0, 0, 0, 0, 0

            if os.path.exists( DIR + "/" + str(j) + "/result/" + tool + ".results.txt" ) == True:
                #print ( DIR + "/" + str(j) + "/result/" + tool + ".results.txt" )
                inputdf = pd.read_csv ( DIR + "/" + str(j) + "/result/" + tool + ".results.txt", sep = "\t", header = None)

                for k in range (inputdf.shape[0]):
                    if inputdf.iloc[k][0] == "score":
                        score = int ( inputdf.iloc[k][1].split("/")[0])
                    if inputdf.iloc[k][0] == "Y-index":
                        Yindex = float ( inputdf.iloc[k][1] )
                    if inputdf.iloc[k][0] == "ARI":
                        ARI = float ( inputdf.iloc[k][1] )
                    if inputdf.iloc[k][0] == "NUM_CLONE":
                        NUM_CLONE_answer = int ( inputdf.iloc[k][1] )
                    if inputdf.iloc[k][0] == "runningtime":
                        runningtime = int ( inputdf.iloc[k][1] )
                    if inputdf.iloc[k][0] == "F1":
                        if inputdf.iloc[k][1] != "None":
                            f1score = float ( inputdf.iloc[k][1] )
                        else:
                            f1score = 0

            result.acc (i, j, score, Yindex, ARI, NUM_CLONE_answer, runningtime, f1score)
    
    
    for i, tool in enumerate( toollist ):
        y_pred =  ( result.NUM_CLONE_record  [i , : ]  )
        y_ans = np.array ( [NUM_CLONE_ans] * ( kwargs["BENCHMARK_END"] - kwargs ["BENCHMARK_START"] + 1 ) )
        RMSE = mean_squared_error(y_pred, y_ans)**0.5
        final.acc (i, tt, np.median (result.score_record [i, : ] ) , np.median (result.ARI_record [i, : ]) , RMSE )

    kwargs["SAMPLENAME"] = DIR.split("/")[-1]
    ax_col = int (tt % NUM_COL)
    ax_row  = int (tt / NUM_COL)

    #print ( "tt= {}\tDIR_index = {}\tax_row = {}\tax_col = {}".format (tt, DIR_index,ax_row, ax_col))
    tt += 1
    
    drawfigure_MS (result, toollist, toollist_concise, axMS, ax_row, ax_col, **kwargs)
    drawfigure_EC (result, toollist, toollist_concise, axEC, ax_row, ax_col, NUM_CLONE_ans, **kwargs)
    
figMS.savefig(kwargs["OUTPUT_MS_JPG"])
figEC.savefig(kwargs["OUTPUT_EC_JPG"])



final.score_record = final.score_record [: , 0:tt]
final.ARI_record = final.ARI_record [: , 0:tt]
final.RMSE_record = final.RMSE_record [: , 0:tt]

drawfigure_Final (final, toollist, toollist_concise, **kwargs)

#python3 /data/project/Alzheimer/YSscript/cle/2.CellData_pipe3_benchmark.py --INPUT_DIR /data/project/Alzheimer/CLEMENT/03.combinedoutput/2.CellData/CellData_1D/n500_125x/parent_0/fp_0.0/axis_-1 --SAMPLENAME M3-9 --CONDITIONNAME n500_125x/parent_0/fp_0.0/axis_-1 --BENCHMARK_START 0 --BENCHMARK_END 3 --OUTPUT_MS_JPG /data/project/Alzheimer/CLEMENT/03.combinedoutput/2.CellData/CellData_1D/n500_125x/parent_0/fp_0.0/axis_-1/BM_MS.jpg --OUTPUT_EC_JPG /data/project/Alzheimer/CLEMENT/03.combinedoutput/2.CellData/CellData_1D/n500_125x/parent_0/fp_0.0/axis_-1/BM_EC.jpg --OUTPUT_FINAL_JPG /data/project/Alzheimer/CLEMENT/03.combinedoutput/2.CellData/CellData_1D/n500_125x/parent_0/fp_0.0/axis_-1/BM_FINAL.jpg --OUTPUT_FINAL_TABLE /data/project/Alzheimer/CLEMENT/03.combinedoutput/2.CellData/CellData_1D/n500_125x/parent_0/fp_0.0/axis_-1/BM_FINAL.tsv



#######################################################################################################################################################################################################################################

