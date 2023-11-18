import pandas as pd
import numpy as np
import argparse
import os

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





def drawfigure (result, toollist, toollist_concise, **kwargs):
    import palettable
    import matplotlib
    import numpy as np
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

    # font_dir = "/home/goldpm1/miniconda3/envs/cnvpytor/lib/python3.7/site-packages/matplotlib/mpl-data/fonts/ttf/"
    # font_dirs = matplotlib.font_manager.findSystemFonts(fontpaths=font_dir, fontext='ttf')
    # for font in font_dirs:
    #     matplotlib.font_manager.fontManager.addfont(font)
    matplotlib.rcParams["font.family"] = 'arial'


    output_result = open (kwargs["OUTPUT_RESULT_TABLE"], "w")
    print ("\n")
    for i, tool in enumerate( toollist ):
        print ("{}\t{}\t{}\t{}\t{}\t{}\t{}".format( tool, round ( np.mean (result.score_record [i , : ]) , 2) ,  round ( np.std (result.score_record [i , : ]) , 2), round ( np.mean (result.NUM_CLONE_record [i , : ]) , 2) ,  round ( np.std (result.NUM_CLONE_record [i , : ]) , 2) , round ( np.mean (result.ARI_record [i , : ]) , 2) ,  round ( np.std (result.ARI_record [i , : ]) , 2)  )  )
        print ("{}\t{}\t{}\t{}\t{}\t{}\t{}".format( tool, round ( np.mean (result.score_record [i , : ]) , 2) ,  round ( np.std (result.score_record [i , : ]) , 2), round ( np.mean (result.NUM_CLONE_record [i , : ]) , 2) ,  round ( np.std (result.NUM_CLONE_record [i , : ]) , 2), round ( np.mean (result.ARI_record [i , : ]) , 2) ,  round ( np.std (result.ARI_record [i , : ]) , 2)   ), file = output_result)
    print ("\n")

    # Seaborn을 위해 df를 만들기
    df = pd.DataFrame (columns = ["tool", "score",  "ARI", "NUM_CLONE_answer", "runningtime"] )
    matrix = []
    for i, tool in enumerate( toollist ):
        for j in range ( kwargs["BENCHMARK_START"]  , kwargs["BENCHMARK_END"]  + 1 ):
            matrix.append ( [tool, result.score_record[i][j], result.ARI_record[i][j], result.NUM_CLONE_record[i][j], result.runningtime_record[i][j]] )

    df = (pd.DataFrame.from_records (matrix, columns = df.columns))

    print (df)

    # # t test 검정하기
    # with open (kwargs["OUTPUT_TTEST"] , "w", encoding = "utf8" ) as output_file_ttest:
    #     combi = list(itertools.combinations(range( len(toollist)), 2))
    #     for a, b in combi:
    #         print ("{} - {}".format( toollist[a], toollist[b] ), file = output_file_ttest )
    #         for col in df.columns[1:]:       # score ~ runningtime 까지 서로 비교
    #             aa =  list (df [ df ["tool"] == toollist[a] ][col])
    #             bb = list (df [ df ["tool"] == toollist[b] ][col])
            
    #             statistics, pvalue = ttest_ind(aa, bb)
    #             print ("\t{} ->  pvalue = {}".format( col,  round(pvalue,5) ) , file = output_file_ttest)
    
    


    fig, ax = matplotlib.pyplot.subplots(1, 3, figsize = (14,6))
    fig.subplots_adjust (wspace = 0.2, hspace = 0.8, bottom = 0.15, top = 0.85, left = 0.05, right = 0.95)
    matplotlib.pyplot.suptitle ( kwargs["CONDITIONNAME"] + "  " + kwargs["SAMPLENAME"] , fontsize = 30, y = 0.98, fontweight = "semibold")
    #matplotlib.pyplot.title ( kwargs["CONDITIONNAME"] , fontsize = 16, y = 0.88, fontweight = "semibold")
    

    #1.
    sns.boxplot (data = df, x = "tool", y = "score",  ax = ax[0],  linewidth = 1)
    sns.swarmplot (data = df, x = "tool", y = "score",  color = ".25", ax = ax[0 ])

    #2.
    sns.boxplot (data = df, x = "tool", y = "ARI",   ax = ax[1], linewidth = 1)
    sns.stripplot (data = df, x = "tool", y = "ARI",  color = ".15", ax = ax[1])

    #3.
    #sns.swarmplot (data = df, x = "tool", y = "NUM_CLONE_answer", hue = "tool", ax = ax[2])
    for j, tool in enumerate( toollist ):
        value_count = df [df ["tool"] == tool ] ["NUM_CLONE_answer"].value_counts()
        #print ("{} -> value_count : {}".format (tool, value_count))
        value_count_dict = {}      # {2:8, 3:69, 4:148, 5:112, 6:26, 7:17}
        for i in value_count.index:
            value_count_dict [i] = value_count.loc[i]
            ax[2].scatter ( x = j, y = i, s = value_count_dict[i] * 80, color = colorlist[j])

    # 사각형 그리기
    if "clone" in kwargs["SAMPLENAME"]:   # simData
        NUM_CLONE_ans = int ( kwargs["SAMPLENAME"].split("_")[-1] )

    else:
        if kwargs["SAMPLENAME"].count("_") == 0:  # 1D
            NUM_CLONE_ans = 3 if  "M1" in kwargs["SAMPLENAME"] else 4
        elif kwargs["SAMPLENAME"].count("_") == 1:  # 2D
            if kwargs["SAMPLENAME"].count("M1") == 2:
                NUM_CLONE_ans = 3
            if kwargs["SAMPLENAME"].count("M2") == 2:
                NUM_CLONE_ans = 4
            if (kwargs["SAMPLENAME"].count("M1") == 1) & (kwargs["SAMPLENAME"].count("M2") == 1):
                NUM_CLONE_ans = 5
        elif kwargs["SAMPLENAME"].count("_") == 2:  #3D
            if kwargs["SAMPLENAME"].count("M1") == 3:   # "M1-4_M1-6_M1-8"
                NUM_CLONE_ans = 3
            elif (kwargs["SAMPLENAME"].count("M2") == 3) | (kwargs["SAMPLENAME"].count("M3") == 3) :    # "M2-2_M2-4_M2-8"
                NUM_CLONE_ans = 4 
            elif (kwargs["SAMPLENAME"].count("M1") == 1) & (kwargs["SAMPLENAME"].count("M2") == 1) & (kwargs["SAMPLENAME"].count("M3") == 1) :    # "M1-1_M2-1_M3_1"
                NUM_CLONE_ans = 6
            else:
                NUM_CLONE_ans = 5
        NUM_CLONE_ans += int (kwargs["CONDITIONNAME"].split("/")[-3].split("_")[-1])      # parent 개수를 더해주기

    rect = matplotlib.patches.Rectangle((-0.5, NUM_CLONE_ans - 0.5),                                # 사각형 꼭지점의 시작위치
                                                    ax[2].get_xlim()[1] - ax[2].get_xlim()[0] + 0.5, 1,        # x 길이, y 길이
                                                    linewidth=0.5, edgecolor='red', facecolor='black', alpha=0.3)
    ax[2].add_patch(rect)


    #4.
    #sns.violinplot (data = df, x = "tool", y = "runningtime", hue = "tool", ax = ax[3], edgecolor = "black")


    for i, ax_individual in enumerate(ax):
        ax_individual.set_xlabel("")
        ax_individual.tick_params(axis = 'x', rotation = 35)
        ax_individual.xaxis.label.set_fontsize(13)
        ax_individual.yaxis.label.set_fontsize(13)
        #ylim_tuple = ax_individual.get_ylim
        xmin, xmax, ymin, ymax = ax_individual.axis()
        ymax = (ymax - ymin) * 0.2 + ymax
        ax_individual.axis ( [ xmin, xmax, ymin, ymax ]) 
        ax_individual.set_title (chr (i+65) + "." , fontsize = 15, loc = "left", fontweight = "semibold")    # A. B. C.
        ax_individual.set_xticks ( list ( range ( len(toollist)  ) ) )   # 숫자로 tick을 박아놓고
        ax_individual.set_xticklabels( toollist_concise )  # 표시는 toollist_concise로
        ax_individual.legend().set_visible(False)   # Legend는 안 보이게 한다


    ax[2].set_yticks ( list ( range ( 1, int(ax[2].get_xlim()[1]) + 1) ) )  # 정수만 표시해준다
    # ax[3].set_ylabel ("runningtime (s)")
    # ax[3].set_yticks ( range ( 0, int ( ax[3].get_ylim() [1]  ) + 1, 60  ))

    

    fig.savefig (kwargs["OUTPUT_JPG"])






if __name__ == "__main__":    
    import argparse, subprocess, os

    parser = argparse.ArgumentParser(description='The below is usage direction.')
    parser.add_argument('--INPUT_DIR', type = str)
    parser.add_argument('--LOG_DIR', type = str)
    parser.add_argument('--SAMPLENAME', type = str)
    parser.add_argument('--CONDITIONNAME', type = str)
    parser.add_argument('--BENCHMARK_START', type = int)
    parser.add_argument('--BENCHMARK_END', type = int)
    parser.add_argument('--OUTPUT_JPG', type = str)
    parser.add_argument('--OUTPUT_TTEST', type = str)
    parser.add_argument('--OUTPUT_RESULT_TABLE', type = str, default = "")

    kwargs = {}
    args = parser.parse_args()

    kwargs["INPUT_DIR"] = args.INPUT_DIR
    kwargs["LOG_DIR"] = args.LOG_DIR
    kwargs["SAMPLENAME"] = args.SAMPLENAME
    kwargs["CONDITIONNAME"] = args.CONDITIONNAME
    kwargs["BENCHMARK_START"] = int(args.BENCHMARK_START)
    kwargs["BENCHMARK_END"] = int(args.BENCHMARK_END)
    kwargs["OUTPUT_TTEST"] = args.OUTPUT_TTEST
    kwargs["OUTPUT_JPG"] = args.OUTPUT_JPG
    kwargs["OUTPUT_RESULT_TABLE"] = args.OUTPUT_RESULT_TABLE


    # if float(args.FP_RATIO) == 0:
    #     kwargs["FP_RATIO"] = int (0)

    print ("LOG_DIR = {}".format ( kwargs["LOG_DIR"] ))

    toollist = ["CLEMENT_decision", "pyclonevi", "sciclone", "quantumclone", "simpleK_elbow", "simpleK_silhouette", "simpleK_gap"]
    toollist_concise = ["CLEMENT", "pyclonevi", "sciclone", "qc", "simpleK_elb", "simpleK_sil", "simpleK_gap*"]

    result = ResultClass(toollist, **kwargs)

    for j in range (kwargs["BENCHMARK_START"], kwargs["BENCHMARK_END"] + 1):
        for i, tool in enumerate( toollist ):            
            score, Yindex, ARI, NUM_CLONE_answer, runningtime, f1score = 0, 0, 0, 0, 0, 0

            #print (  kwargs["INPUT_DIR"] + "/" + str(j) + "/result/" + tool + ".results.txt"  )

            if os.path.exists( kwargs["INPUT_DIR"] + "/" + str(j) + "/result/" + tool + ".results.txt" ) == True:
                inputdf = pd.read_csv ( kwargs["INPUT_DIR"] + "/" + str(j) + "/result/" + tool + ".results.txt", sep = "\t", header = None)
                
                subprocess.run(["mkdir", "-p",  kwargs["LOG_DIR"] + "/" + str(j) , kwargs["INPUT_DIR"] + "/" + str(j) + "/log"], shell=False)
                subprocess.run(["cp", "-rf",  kwargs["LOG_DIR"] + "/" + str(j) , kwargs["INPUT_DIR"] + "/" + str(j) + "/log"], shell=False)
                

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

    drawfigure (result, toollist, toollist_concise, **kwargs)
        
