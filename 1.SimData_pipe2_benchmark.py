import pandas as pd
import numpy as np
import argparse
import os

class ResultClass:
    def __init__(self, toollist, **kwargs):
        self.score_record = np.zeros ( ( len (toollist), kwargs["BENCHMARK_END"] - kwargs["BENCHMARK_START"] + 1  ), dtype = "int" )
        self.Yindex_record = np.zeros ( (len (toollist)  , kwargs["BENCHMARK_END"] - kwargs["BENCHMARK_START"] + 1 ), dtype = "float" )
        self.ARI_record = np.zeros ( (len (toollist)  , kwargs["BENCHMARK_END"] - kwargs["BENCHMARK_START"] + 1 ), dtype = "float" )
        self.NUM_CLONE_record = np.zeros ( (len (toollist) , kwargs["BENCHMARK_END"] - kwargs["BENCHMARK_START"] + 1 ), dtype = "int" )
        self.runningtime_record = np.zeros ( (len (toollist) , kwargs["BENCHMARK_END"] - kwargs["BENCHMARK_START"] + 1 ), dtype = "float" )
        self.f1score_record = np.zeros ( (len (toollist) , kwargs["BENCHMARK_END"] - kwargs["BENCHMARK_START"] + 1 ), dtype = "float" )

    def acc (self, I, J, score, Yindex, ARI, NUM_CLONE_answer, runningtime, f1score):
        self.score_record[J , I] = score
        self.Yindex_record[J , I] = Yindex
        self.ARI_record[J , I] = ARI
        self.NUM_CLONE_record[J , I] = NUM_CLONE_answer
        self.runningtime_record[J , I] = runningtime
        self.f1score_record[J , I] = f1score





def drawfigure (result, toollist, i_real, **kwargs):
    import palettable
    import matplotlib
    import numpy as np
    import seaborn as sns
    from scipy.stats import ttest_ind
    import itertools

    tabl = palettable.tableau.Tableau_20.mpl_colors
    safe7 = palettable.cartocolors.qualitative.Vivid_10.mpl_colors
    Gr_10 = palettable.scientific.sequential.GrayC_20.mpl_colors

    colorlist = [i for i in safe7]
    sns.set_style("white")
    #sns.set_palette("tab10")
    sns.set_palette(sns.color_palette(colorlist))

    font_dir = "/home/goldpm1/miniconda3/envs/cnvpytor/lib/python3.7/site-packages/matplotlib/mpl-data/fonts/ttf/"
    font_dirs = matplotlib.font_manager.findSystemFonts(fontpaths=font_dir, fontext='ttf')
    for font in font_dirs:
        matplotlib.font_manager.fontManager.addfont(font)
    matplotlib.rcParams["font.family"] = 'arial'


    # Seaborn을 위해 df를 만들기
    df = pd.DataFrame (columns = ["tool", "score", "Y-index", "ARI", "NUM_CLONE_answer", "f1score", "runningtime"] )
    matrix = []
    for j, tool in enumerate( toollist ):
        for k in range (1, i_real ):
            matrix.append ( [tool, result.score_record[j][k], result.Yindex_record[j][k], result.ARI_record[j][k], result.NUM_CLONE_record[j][k], result.f1score_record[j][k], result.runningtime_record[j][k]] )

    df = (pd.DataFrame.from_records (matrix, columns = df.columns))


    # t test 검정하기
    with open (kwargs["OUTPUT_TTEST"] , "w", encoding = "utf8" ) as output_file_ttest:
        combi = list(itertools.combinations(range( len(toollist)), 2))
        for a, b in combi:
            print ("{} - {}".format( toollist[a], toollist[b] ), file = output_file_ttest )
            for col in df.columns[1:]:       # score ~ runningtime 까지 서로 비교
                aa =  list (df [ df ["tool"] == toollist[a] ][col])
                bb = list (df [ df ["tool"] == toollist[b] ][col])
            
                statistics, pvalue = ttest_ind(aa, bb)
                print ("\t{} ->  pvalue = {}".format( col,  round(pvalue,5) ) , file = output_file_ttest)

    NUM_CLONE_answer = int ( kwargs["COMBINED_OUTPUT_DIR"].split("_")[-1] )
    print ("NUM_CLONE_answer : {}".format(NUM_CLONE_answer))


    fig, ax = matplotlib.pyplot.subplots(1, 3, figsize = (14,6))
    fig.subplots_adjust (wspace = 0.5, hspace = 1, bottom = 0.15, top = 0.85)
    #fig.text (x = 0.5, y = 0.02, ha = "center", s = "NUM_PARENT = {}, FP_RATIO = {}, AXIS_RATIO ={}".format(kwargs["NUM_PARENT"], kwargs["FP_RATIO"], kwargs["AXIS_RATIO"]), backgroundcolor = Gr_10[2], fontsize = 12)


    sns.boxplot (data = df, x = "tool", y = "score",  ax = ax[0], linewidth = 1)
    sns.swarmplot (data = df, x = "tool", y = "score", color = ".25", ax = ax[0])

    sns.boxplot (data = df, x = "tool", y = "ARI",  ax = ax[1], linewidth = 1)
    sns.stripplot (data = df, x = "tool", y = "ARI", color = ".15", ax = ax[1])

    # sns.swarmplot (data = df, x = "tool", y = "NUM_CLONE_answer", ax = ax[2])
    # ax[2].set_yticks ( range (np.min (df["NUM_CLONE_answer"]) , np.max (df["NUM_CLONE_answer"]) + 1  ))

    for j, tool in enumerate( toollist ):
        value_count = df [df ["tool"] == tool ] ["NUM_CLONE_answer"].value_counts()

        #print ("{} -> value_count : {}".format (tool, value_count))

        value_count_dict = {}      # {2:8, 3:69, 4:148, 5:112, 6:26, 7:17}
        for i in value_count.index:
            value_count_dict [i] = value_count.loc[i]
            ax[2].scatter ( x = NUM_CLONE_answer + (j * 0.1), y = i, s = value_count_dict[i] * 3, c = tabl[j])

    for i, ax_individual in enumerate(ax):
        ax_individual.set_xlabel("")
        ax_individual.tick_params(axis = 'x', rotation = 35)
        ax_individual.xaxis.label.set_fontsize(13)
        ax_individual.yaxis.label.set_fontsize(13)
        #ylim_tuple = ax_individual.get_ylim
        xmin, xmax, ymin, ymax = ax_individual.axis()
        ymax = (ymax - ymin) * 0.2 + ymax
        ax_individual.axis ( [ xmin, xmax, ymin, ymax ]) 

        #ax_individual.set_title (chr (i+65) + "." , fontsize = 15, loc = "left")
    matplotlib.pyplot.suptitle (kwargs["SAMPLENAME"] , fontsize = 30, y = 0.98, fontweight = "semibold")

    fig.savefig (kwargs["OUTPUT_JPG"])



if __name__ == "__main__":    
    import argparse 

    parser = argparse.ArgumentParser(description='The below is usage direction.')
    parser.add_argument('--COMBINED_OUTPUT_DIR', type = str)
    parser.add_argument('--SAMPLENAME', type = str)
    parser.add_argument('--BENCHMARK_START', type = int)
    parser.add_argument('--BENCHMARK_END', type = int)
    parser.add_argument('--OUTPUT_TTEST', type = str)
    parser.add_argument('--OUTPUT_JPG', type = str)

    kwargs = {}
    args = parser.parse_args()


    kwargs["COMBINED_OUTPUT_DIR"] = args.COMBINED_OUTPUT_DIR
    kwargs["SAMPLENAME"] = args.SAMPLENAME
    kwargs["BENCHMARK_START"] = int(args.BENCHMARK_START)
    kwargs["BENCHMARK_END"] = int(args.BENCHMARK_END)
    kwargs["OUTPUT_TTEST"] = args.OUTPUT_TTEST
    kwargs["OUTPUT_JPG"] = args.OUTPUT_JPG


    toollist = ["CLEMENT_hard_1st", "pyclonevi", "sciclone", "quantumclone"]

    result = ResultClass(toollist, **kwargs)

    i_real = 0
    for i in range (kwargs["BENCHMARK_START"], kwargs["BENCHMARK_END"] + 1):
        for j, tool in enumerate( toollist ):
            try:
                inputdf = pd.read_csv ( kwargs["COMBINED_OUTPUT_DIR"] + "/" + str(i) + "/" + tool + ".results.txt", sep = "\t", header = None)

                Yindex, runningtime, f1score = 0, 0, 0
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

                result.acc (i_real, j, score, Yindex, ARI, NUM_CLONE_answer, runningtime, f1score)

            except:
                i_real = i_real - 1          # 뭔가 오류가 나서 result가 안 나오는 것도 있다
                continue
        i_real = i_real + 1

    drawfigure (result, toollist, i_real, **kwargs)
        


    # bash 1.SimData_pipe2_benchmark.sh --COMBINED_OUTPUT_DIR /data/project/Alzheimer/YSscript/EM_MRS/data/combinedoutput/simulation_2D/clone_4 --SAMPLENAME simulation_2D_clone_4 --BENCHMARK_NO 100 --OUTPUT_TTEST /data/project/Alzheimer/YSscript/EM_MRS/data/combinedoutput/simulation_2D/clone_4/ttest.txt --OUTPUT_JPG /data/project/Alzheimer/YSscript/EM_MRS/data/combinedoutput/simulation_2D/clone_4/benchmark.pdf
    # bash 1.SimData_pipe2_benchmark.sh --COMBINED_OUTPUT_DIR /data/project/Alzheimer/YSscript/EM_MRS/data/combinedoutput/simulation_2D/clone_5 --SAMPLENAME simulation_2D_clone_5 --BENCHMARK_NO 100 --OUTPUT_TTEST /data/project/Alzheimer/YSscript/EM_MRS/data/combinedoutput/simulation_2D/clone_5/ttest.txt --OUTPUT_JPG /data/project/Alzheimer/YSscript/EM_MRS/data/combinedoutput/simulation_2D/clone_5/benchmark.pdf