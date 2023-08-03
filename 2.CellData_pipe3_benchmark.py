def drawfigure (result, toollist, toollist_concise, ax, ax_row, ax_col, **kwargs):
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
    sns.set_style("white")
    #sns.set_palette("tab10")
    sns.set_palette(sns.color_palette(colorlist))

    # font_dir = "/home/goldpm1/miniconda3/envs/cnvpytor/lib/python3.7/site-packages/matplotlib/mpl-data/fonts/ttf/"
    # font_dirs = matplotlib.font_manager.findSystemFonts(fontpaths=font_dir, fontext='ttf')
    # for font in font_dirs:
    #     matplotlib.font_manager.fontManager.addfont(font)
    matplotlib.rcParams["font.family"] = 'arial'


    # Seaborn을 위해 df를 만들기
    df = pd.DataFrame (columns = ["tool", "score", "Y-index", "ARI", "NUM_CLONE_answer", "f1score", "runningtime"] )
    matrix = []
    for j, tool in enumerate( toollist ):
        for k in range (1, len(result.score_record[0]) ):
            matrix.append ( [tool, result.score_record[j][k], result.Yindex_record[j][k], result.ARI_record[j][k], result.NUM_CLONE_record[j][k], result.f1score_record[j][k], result.runningtime_record[j][k]] )

    df = (pd.DataFrame.from_records (matrix, columns = df.columns))

    #print (df)


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







import numpy as np 
import pandas as pd
import os, glob
import matplotlib.pyplot

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



########################################################################################

import argparse 

parser = argparse.ArgumentParser(description='The below is usage direction.')
parser.add_argument('--INPUT_DIR', type = str, default = "")
parser.add_argument('--SAMPLENAME', type = str, default = "")
parser.add_argument('--CONDITIONNAME', type = str, default = "")
parser.add_argument('--BENCHMARK_START', type = int, default = 0)
parser.add_argument('--BENCHMARK_END', type = int, default = 10)
parser.add_argument('--OUTPUT_JPG', type = str, default = "")

kwargs = {}
args = parser.parse_args()

kwargs["INPUT_DIR"] = args.INPUT_DIR
kwargs["SAMPLENAME"] = args.SAMPLENAME
kwargs["CONDITIONNAME"] = args.CONDITIONNAME
kwargs["BENCHMARK_START"] = args.BENCHMARK_START
kwargs["BENCHMARK_END"] = args.BENCHMARK_END
kwargs["OUTPUT_JPG"] = args.OUTPUT_JPG



INPUT_DIR_LIST = sorted (  glob.glob(kwargs["INPUT_DIR"] + "/*") ) 

toollist = ["CLEMENT_decision", "pyclonevi", "sciclone", "quantumclone", "simpleK_elbow", "simpleK_silhouette", "simpleK_gap"]
toollist_concise = ["CLEMENT", "pyclonevi", "sciclone", "qc", "simpleK_elb", "simpleK_sil", "simpleK_gap*"]


####################################################################################################
if len (INPUT_DIR_LIST) > 15:
    NUM_ROW = 4
elif len (INPUT_DIR_LIST) > 9:
    NUM_ROW = 3

fig, ax = matplotlib.pyplot.subplots(nrows = NUM_ROW,  ncols = int (len(INPUT_DIR_LIST) / NUM_ROW) + 1, figsize = (18,NUM_ROW * 3.5))
fig.subplots_adjust (wspace = 0.4, hspace = 0.5, bottom = 0.1, top = 0.9, left = 0.05, right = 0.95)
fig.suptitle( "{}".format(kwargs["CONDITIONNAME"]), fontsize = 30, fontweight = "bold")

####################################################################################################



for DIR_index, DIR in enumerate( INPUT_DIR_LIST) :
    if os.path.isfile(DIR) == True:
        continue
    
    print (DIR)

    result = ResultClass(toollist, **kwargs)

    for j in range (kwargs["BENCHMARK_START"], kwargs["BENCHMARK_END"] + 1):
        break_check = False
        for i, tool in enumerate( toollist ):
            score, Yindex, ARI, NUM_CLONE_answer, runningtime, f1score = 0, 0, 0, 0, 0, 0

            if os.path.exists( DIR + "/" + str(j) + "/result/" + tool + ".results.txt" ) == True:
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


    ax_row  = int (DIR_index % NUM_ROW)
    ax_col = int (DIR_index / NUM_ROW)
    
    drawfigure (result, toollist, ax, ax_row, ax_col, **kwargs)
    
fig.savefig(kwargs["OUTPUT_JPG"])

