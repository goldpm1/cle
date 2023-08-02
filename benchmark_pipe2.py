import pandas as pd
import numpy as np
import argparse
import os

class ResultClass:
    def __init__(self, toollist, **kwargs):
        self.score_record = np.zeros ( ( len (toollist), kwargs["BENCHMARK_NO"] + 1), dtype = "int" )
        self.Yindex_record = np.zeros ( (len (toollist)  , kwargs["BENCHMARK_NO"] + 1), dtype = "float" )
        self.NUM_CLONE_record = np.zeros ( (len (toollist) , kwargs["BENCHMARK_NO"] + 1), dtype = "int" )
        self.runningtime_record = np.zeros ( (len (toollist) , kwargs["BENCHMARK_NO"] + 1), dtype = "float" )
        self.f1score_record = np.zeros ( (len (toollist) , kwargs["BENCHMARK_NO"] + 1), dtype = "float" )

    def acc (self, I, J, score, Yindex, NUM_CLONE_answer, runningtime, f1score):
        self.score_record[J , I] = score
        self.Yindex_record[J , I] = Yindex
        self.NUM_CLONE_record[J , I] = NUM_CLONE_answer
        self.runningtime_record[J , I] = runningtime
        self.f1score_record[J , I] = f1score





def drawfigure (result, toollist, **kwargs):
    import palettable
    import matplotlib
    import numpy as np
    import seaborn as sns
    import copy

    tabl = palettable.tableau.Tableau_20.mpl_colors
    safe7 = palettable.cartocolors.qualitative.Vivid_10.mpl_colors
    Gr_10 = palettable.scientific.sequential.GrayC_20.mpl_colors

    colorlist = [i for i in safe7]
    sns.set_style("whitegrid")
    #sns.set_palette("tab10")
    sns.set_palette(sns.color_palette(colorlist))

    font_dir = "/home/goldpm1/miniconda3/envs/cnvpytor/lib/python3.7/site-packages/matplotlib/mpl-data/fonts/ttf/"
    font_dirs = matplotlib.font_manager.findSystemFonts(fontpaths=font_dir, fontext='ttf')
    for font in font_dirs:
        matplotlib.font_manager.fontManager.addfont(font)
    matplotlib.rcParams["font.family"] = 'arial'


    # Seaborn을 위해 df를 만들기
    df = pd.DataFrame (columns = ["tool", "score", "Y-index", "NUM_CLONE_answer", "f1score", "runningtime"] )
    matrix = []
    for j, tool in enumerate( toollist ):
        for k in range (1, len(result.score_record[0])):
            matrix.append ( [tool, result.score_record[j][k], result.Yindex_record[j][k], result.NUM_CLONE_record[j][k], result.f1score_record[j][k], result.runningtime_record[j][k]] )

    df = (pd.DataFrame.from_records (matrix, columns = df.columns))


    fig, ax = matplotlib.pyplot.subplots(1, 5, figsize = (18,6))
    fig.subplots_adjust (wspace = 0.5, hspace = 1, bottom = 0.15, top = 0.85)
    fig.text (x = 0.5, y = 0, ha = "center", s = "NUM_PARENT = {}, FP_RATIO = {}, AXIS_RATIO ={}".format(kwargs["NUM_PARENT"], kwargs["FP_RATIO"], kwargs["AXIS_RATIO"]), backgroundcolor = Gr_10[2], fontsize = 12)


    sns.boxplot (data = df, x = "tool", y = "score",  ax = ax[0], linewidth = 1)
    sns.swarmplot (data = df, x = "tool", y = "score", color = ".25", ax = ax[0])

    sns.boxplot (data = df, x = "tool", y = "Y-index",  ax = ax[1], linewidth = 1)
    sns.stripplot (data = df, x = "tool", y = "Y-index", color = ".15", ax = ax[1])

    sns.swarmplot (data = df, x = "tool", y = "NUM_CLONE_answer", ax = ax[2])
    ax[2].set_yticks ( range (np.min (df["NUM_CLONE_answer"]) , np.max (df["NUM_CLONE_answer"]) + 1  ))

    sns.violinplot (data = df, x = "tool", y = "f1score", ax = ax[3], linewidth = 1, edgecolor = "black")
    ax[3].set_ylabel ("f1 score for FP")

    sns.barplot (data = df, x = "tool", y = "runningtime", ax = ax[4], linewidth = 1, edgecolor = "black")
    ax[4].set_ylabel ("runningtime (s)")
    ax[4].set_yticks ( range ( int(ax[4].get_ylim() [0]) , int ( ax[4].get_ylim() [1]  ) + 1, 60  ))


    for i, ax_individual in enumerate(ax):
        ax_individual.set_xlabel("")
        ax_individual.tick_params(axis = 'x', rotation = 45)
        ax_individual.xaxis.label.set_fontsize(13)
        ax_individual.yaxis.label.set_fontsize(13)
        ax_individual.set_title (chr (i+65) + "." , fontsize = 15, loc = "left")
    matplotlib.pyplot.suptitle (kwargs["SAMPLENAME"], fontsize = 30, y = 1, fontweight = "semibold")


    fig.savefig (kwargs["OUTPUT_RESULT"])






if __name__ == "__main__":    
    import argparse 

    parser = argparse.ArgumentParser(description='The below is usage direction.')
    parser.add_argument('--INPUT_DIR', type = str)
    parser.add_argument('--SAMPLENAME', type = str)
    parser.add_argument('--BENCHMARK_NO', type = int)
    parser.add_argument('--RANDOM_PICK', type = int)
    parser.add_argument('--NUM_PARENT', type = int)
    parser.add_argument('--FP_RATIO', type = float)
    parser.add_argument('--AXIS_RATIO', type = float)
    parser.add_argument('--OUTPUT_RESULT', type = str)


    kwargs = {}
    args = parser.parse_args()

    kwargs["INPUT_DIR"] = args.INPUT_DIR
    kwargs["SAMPLENAME"] = args.SAMPLENAME
    kwargs["BENCHMARK_NO"] = int(args.BENCHMARK_NO)
    kwargs["RANDOM_PICK"] = int(args.RANDOM_PICK)
    kwargs["NUM_PARENT"] = int(args.NUM_PARENT)
    kwargs["FP_RATIO"] = float(args.FP_RATIO)
    kwargs["AXIS_RATIO"] = float(args.AXIS_RATIO)
    kwargs["OUTPUT_RESULT"] = args.OUTPUT_RESULT

    # os.system ("rm -rf " + str (kwargs["OUTPUT_RESULT"].split("/")[:-1]))
    # os.system ("mkdir -p " +str ( kwargs["OUTPUT_RESULT"].split("/")[:-1]) )
    

    #toollist = ["TGIL_hard_1st", "TGIL_hard_2nd" , "TGIL_soft", "pyclonevi", "sciclone"]
    toollist = ["TGIL_hard_1st", "pyclonevi", "sciclone"]

    result = ResultClass(toollist, **kwargs)

    for i in range (1, kwargs["BENCHMARK_NO"] + 1):
        for j, tool in enumerate( toollist ):
            inputdf = pd.read_csv (kwargs["INPUT_DIR"] + "/" + kwargs["SAMPLENAME"] + "/" + str(kwargs["RANDOM_PICK"]) + "_" + str(kwargs["NUM_PARENT"]) + "_" + str(kwargs["FP_RATIO"]) + "_" + str(kwargs["AXIS_RATIO"]) + "/" + str(i) + "/" + tool + ".results.txt", sep = "\t", header = None)

            for k in range (inputdf.shape[0]):
                if inputdf.iloc[k][0] == "score":
                    score = int ( inputdf.iloc[k][1].split("/")[0])
                if inputdf.iloc[k][0] == "Y-index":
                    Yindex = float ( inputdf.iloc[k][1] )
                if inputdf.iloc[k][0] == "NUM_CLONE":
                    NUM_CLONE_answer = int ( inputdf.iloc[k][1] )
                if inputdf.iloc[k][0] == "runningtime":
                    runningtime = int ( inputdf.iloc[k][1] )
                if inputdf.iloc[k][0] == "F1":
                    if inputdf.iloc[k][1] != "None":
                        f1score = float ( inputdf.iloc[k][1] )
                    else:
                        f1score = 0

            result.acc (i, j, score, Yindex, NUM_CLONE_answer, runningtime, f1score)


    drawfigure (result, toollist, **kwargs)
        
