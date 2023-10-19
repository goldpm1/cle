import palettable
import matplotlib
import extract
import numpy as np
import seaborn as sns
import pandas as pd
from scipy.stats import kde
from sklearn.decomposition import TruncatedSVD, PCA


def drawfigure_1d_hard(bunch, np_vaf, output_filename, sum_mixture_beforenormalization, **kwargs):
    
    if (bunch.makeone_index == []) | (bunch.makeone_index == None) | (bunch.checkall_strict == False):  
        matplotlib.pyplot.style.use("Solarize_Light2")
    else:
        matplotlib.pyplot.style.use("seaborn-white")


    tabl = palettable.tableau.Tableau_20.mpl_colors
    Gr_10 = palettable.scientific.sequential.GrayC_20.mpl_colors
    colorlist = [i for i in tabl]
    if (bunch.includefp == True) & (bunch.fp_index != -1):
        colorlist[bunch.fp_index] = Gr_10[8]


    fig, ax = matplotlib.pyplot.subplots(ncols=1, nrows = 1, figsize=(6, 6))
    #ax.set_axis( [0,  np.max(np_vaf[:, :]) * 2.1, 0,  np.max(np_vaf[:, :]) * 2.1])

    matplotlib.rcParams["font.family"] = kwargs["FONT_FAMILY"]

    max_y = 0

    x = np.linspace(0, 2, 200)

    for k in sorted(list(set(bunch.membership))):     
        np_vaf_new_index, np_vaf_new = extract.npvaf(np_vaf, bunch.membership, k)           # extract membership1  == clone_num

        try:   
            kde_np_vaf_new = kde.gaussian_kde(np_vaf_new[:, 0] * 2)
            weight = len(np_vaf_new) / kwargs["NUM_MUTATION"]  
            y = kde_np_vaf_new(x) * weight

            if max_y < np.max(y):
                max_y = np.max(y)

            if k in bunch.makeone_index:
                ax.plot(x, y, color=colorlist[k], linewidth=5, label="clone {}".format(k))
            else:
                ax.plot(x, y, color=colorlist[k], linewidth=2, label="clone {}".format(k), linestyle="-.")
            ax.text(bunch.mixture[0, k], np.max(y) * 1.1, "{} (n={})".format(bunch.mixture[0,  k],  np.bincount(bunch.membership)[k]), ha="center", fontdict = {"fontsize": 16, "fontweight" : "bold"})

        except:
            continue

    

    if (kwargs["STEP"] <= (kwargs["COMPULSORY_NORMALIZATION"] - 1)) :
        import matplotlib.patches as patches
        fig_width, fig_height = fig.get_size_inches()
        rect = patches.Rectangle((0, 0), fig_width * fig.dpi, fig_height * fig.dpi, linewidth=8, edgecolor="#C23373", facecolor='none')
        fig.patches.append(rect)
        # Set the layout to tight to ensure the rectangle is not cut off
        matplotlib.pyplot.tight_layout()
        matplotlib.pyplot.suptitle("clone{}.{}-{} (postnormalization) : {}".format(kwargs["NUM_CLONE_ITER"], kwargs["TRIAL"], kwargs["STEP"], round(bunch.likelihood)), fontsize=20)
    else:
        matplotlib.pyplot.suptitle("clone{}.{}-{} (no normalization) : {}".format(kwargs["NUM_CLONE_ITER"], kwargs["TRIAL"], kwargs["STEP"], round(bunch.likelihood)), fontsize=20)
    matplotlib.pyplot.title("Sum_child = {}\tprenorm = {}".format(list(np.round(np.sum(bunch.mixture[:, bunch.makeone_index], axis=1), 2)), sum_mixture_beforenormalization), fontsize=12)

    ax.set_xlim( 0, np.max(np_vaf[:, :]) * 2.1 )
    ax.set_ylim( 0, max_y * 1.3 )

    ax.legend()
    ax.set_xlabel("VAF x 2 of the Sample 1", fontdict={"fontsize": 14})
    ax.set_ylabel("Density", fontdict={"fontsize": 14})

    if output_filename != "NotSave":
        matplotlib.pyplot.savefig(output_filename)


################################################################################################################################################################################################


def drawfigure_1d_soft(bunch, np_vaf, output_filename, **kwargs):
    membership = bunch.membership
    mixture = bunch.mixture
    membership_p_normalize = bunch.membership_p_normalize

    if (bunch.makeone_index == []) | (bunch.makeone_index == None) | (bunch.checkall_strict == False): 
        matplotlib.pyplot.style.use("Solarize_Light2")
    else:
        matplotlib.pyplot.style.use("seaborn-white")

    tabl = palettable.tableau.Tableau_20.mpl_colors
    Gr_10 = palettable.scientific.sequential.GrayC_20.mpl_colors
    colorlist = [i for i in tabl]
    if (bunch.includefp == True) & (bunch.fp_index != -1):
        colorlist[bunch.fp_index] = Gr_10[8]

    # font_dir = "/home/goldpm1/miniconda3/envs/cnvpytor/lib/python3.7/site-packages/matplotlib/mpl-data/fonts/ttf/"
    # font_dirs = matplotlib.font_manager.findSystemFonts( fontpaths=font_dir, fontext='ttf')
    # for font in font_dirs:
    #     matplotlib.font_manager.fontManager.addfont(font)
    matplotlib.rcParams["font.family"] = kwargs["FONT_FAMILY"]

    matplotlib.pyplot.figure(figsize=(6, 6))
    matplotlib.pyplot.xlim(0, 1)
    matplotlib.pyplot.xlabel("VAF x 2 (= Mixture)")
    matplotlib.pyplot.ylabel("count (weighted)")

    if (kwargs["STEP"] <= (kwargs["COMPULSORY_NORMALIZATION"] - 1)) :
        matplotlib.pyplot.suptitle("clone{}.{}-{} (postnormalization) : {}".format(kwargs["NUM_CLONE_ITER"], kwargs["TRIAL"], kwargs["STEP_TOTAL"], round(bunch.likelihood)), fontsize=20)
    else:
        matplotlib.pyplot.suptitle("clone{}.{}-{} (no normalization) : {}".format(kwargs["NUM_CLONE_ITER"], kwargs["TRIAL"], kwargs["STEP_TOTAL"], round(bunch.likelihood)), fontsize=20)
    matplotlib.pyplot.title("sum_child = {}".format(list(np.round(np.sum(bunch.mixture[:, bunch.makeone_index], axis=1), 2)), bunch.makeone_index ), fontsize=12)

    for k in sorted(list(set(bunch.membership))):
        if k in bunch.makeone_index:
            sns.distplot(pd.DataFrame(np_vaf[:, 0] * 2, columns=["vaf"])["vaf"], hist_kws={"weights": membership_p_normalize[:, k], "linewidth": 1.4, "edgecolor": "black"}, kde_kws={
                         "linewidth": 5, "color": "gray"}, color=colorlist[k], kde=False, bins=50, label="cluster {}  (mixture = {})".format(k,  str(round((mixture[0][k]), 2))))
        else:
            sns.distplot(pd.DataFrame(np_vaf[:, 0] * 2, columns=["vaf"])["vaf"], hist_kws={"weights": membership_p_normalize[:, k], "rwidth": 0.8}, color=colorlist[k],
                         kde=False, bins=50, label="cluster {}  (mixture = {})".format(k,  str(round((mixture[0][k]), 2))))
    matplotlib.pyplot. legend()

    matplotlib.pyplot.xlabel("VAF * 2 (= Mixture)")
    matplotlib.pyplot.ylabel("count (weighted)")

    if output_filename != "NotSave":
        matplotlib.pyplot.savefig(output_filename)



################################################################################################################################################################################################



def drawfigure_2d(bunch, np_vaf, output_filename, sum_mixture_beforenormalization, **kwargs):
    samplename_dict = {k: k for k in range(0, bunch.mixture.shape[1])}

    tabl = palettable.tableau.Tableau_20.mpl_colors
    Gr_10 = palettable.scientific.sequential.GrayC_20.mpl_colors
    colorlist = [i for i in tabl]

    matplotlib.rcParams["font.family"] = kwargs["FONT_FAMILY"]

    fig, ax = matplotlib.pyplot.subplots(ncols=1, nrows = 1, figsize=(6, 6))
    #ax.set_axis( [0,  np.max(np_vaf[:, :]) * 2.1, 0,  np.max(np_vaf[:, :]) * 2.1])
    ax.set_xlim(0, np.max(np_vaf[:, :]) * 2.1)
    ax.set_ylim(0, np.max(np_vaf[:, :]) * 2.1)

    #matplotlib.pyplot.figure(figsize=(6, 6))
    #matplotlib.pyplot.axis( [0,  np.max(np_vaf[:, :]) * 2.1, 0,  np.max(np_vaf[:, :]) * 2.1])
    # matplotlib.pyplot.xlabel("VAF x 2 of the Sample 1", fontdict={"fontsize": 14})
    # matplotlib.pyplot.ylabel("VAF x 2 of the Sample 2", fontdict={"fontsize": 14})

    if (kwargs["STEP"] <= (kwargs["COMPULSORY_NORMALIZATION"] - 1)) :
        matplotlib.pyplot.suptitle("clone{}.{}-{} (postnormalization) : {}".format(kwargs["NUM_CLONE_ITER"], kwargs["TRIAL"], kwargs["STEP"], round(bunch.likelihood)), fontsize=20)
    else:
        matplotlib.pyplot.suptitle("clone{}.{}-{} (no normalization) : {}".format(kwargs["NUM_CLONE_ITER"], kwargs["TRIAL"], kwargs["STEP"], round(bunch.likelihood)), fontsize=20)
    matplotlib.pyplot.title("sum_child = {}\tprenorm = {}".format( list (np.round( np.sum (bunch.mixture[:, bunch.makeone_index], axis = 1), 3) ), sum_mixture_beforenormalization), fontsize = 12)


    if (bunch.makeone_index == []) | (bunch.makeone_index == None) | (bunch.checkall_strict == False):  
        matplotlib.pyplot.style.use("Solarize_Light2")
    # elif (bunch.checkall_strict == False):
    #     matplotlib.pyplot.style.use("bmh")
    else:
        matplotlib.pyplot.style.use("seaborn-white")

    if (kwargs["STEP"] <= (kwargs["COMPULSORY_NORMALIZATION"] - 1)) :
        import matplotlib.patches as patches
        fig_width, fig_height = fig.get_size_inches()
        rect = patches.Rectangle((0, 0), fig_width * fig.dpi, fig_height * fig.dpi, linewidth=8, edgecolor="#C23373", facecolor='none')
        fig.patches.append(rect)
        # Set the layout to tight to ensure the rectangle is not cut off
        matplotlib.pyplot.tight_layout()

    
    #matplotlib.pyplot.scatter(np_vaf[:, 0] * 2, np_vaf[:, 1] * 2, color=[colorlist[samplename_dict[k]] for k in bunch.membership])
    for k in range (len ( bunch.membership)):
        if (bunch.fp_index == bunch.membership[k]):
            ax.scatter(np_vaf[k, 0] * 2, np_vaf[k, 1] * 2, color = Gr_10[10]  )
        else:
            ax.scatter(np_vaf[k, 0] * 2, np_vaf[k, 1] * 2, color = colorlist[ samplename_dict [ bunch.membership[k] ]]  )


    for sample_index, sample in enumerate(samplename_dict):
        if (sample_index >= bunch.mixture.shape[1]):
            continue

        # Drawe the centroid based on mixture
        x_mean = round(bunch.mixture[0][sample_index], 2)
        y_mean = round(bunch.mixture[1][sample_index], 2)
        if x_mean < 0.02:
            ax.text(x_mean, y_mean, "{0}".format( [y_mean]), verticalalignment='top', horizontalalignment='right', fontdict = {"fontsize": 14, "fontweight" : "bold"} )    
        elif y_mean < 0.02:
            ax.text(x_mean, -0.02, "{0}".format( [x_mean]), verticalalignment='top', horizontalalignment='center',  fontdict = {"fontsize": 14, "fontweight" : "bold"} )    
        else:
            ax.text(x_mean, y_mean, "{0}".format( [x_mean, y_mean]), verticalalignment='top', horizontalalignment='center', fontdict = {"fontsize": 16, "fontweight" : "bold"} )


        if (bunch.fp_index == sample_index):
            ax.scatter(x_mean, y_mean, marker='s', color= Gr_10[10], edgecolor='black', s=200, label="cluster" + str(sample_index) + " : " + str(list(bunch.membership).count(sample_index)))            
        elif (bunch.makeone_index != []) & (bunch.makeone_index != None):
            if sample_index in bunch.makeone_index:
                ax.scatter(x_mean, y_mean, marker='*', color=colorlist[sample_index], edgecolor='black', s=400, label="cluster" + str(sample_index) + " : " + str(list(bunch.membership).count(sample_index)))
            else:
                ax.scatter(x_mean, y_mean, marker='s', color=colorlist[sample_index], edgecolor='black', s=200, label="cluster" + str(sample_index) + " : " + str(list(bunch.membership).count(sample_index)))
        else:
            ax.scatter(x_mean, y_mean, marker='+', color=colorlist[sample_index], edgecolor='black', s=200, label="cluster" + str(sample_index) + " : " + str(list(bunch.membership).count(sample_index)))

        ax.legend()

    if output_filename != "NotSave":
        matplotlib.pyplot.savefig(output_filename)



################################################################################################################################################################################################

def drawfigure_2d_soft(bunch, np_vaf, output_filename, **kwargs):
    if kwargs["OPTION"] in ["Hard", "hard"]:
        matplotlib.pyplot.style.use("seaborn-white")
    elif kwargs["OPTION"] in ["Soft", "soft"]:
        matplotlib.pyplot.style.use("seaborn-white") 
    elif kwargs["OPTION"] in ["fp", "fp"]:
        matplotlib.pyplot.style.use("Solarize_Light2")

    NUM_MUTATION = len(bunch.membership)
    NUM_CLONE = kwargs["NUM_CLONE"]
    NUM_BLOCK = kwargs["NUM_BLOCK"]

    samplename_dict = {k: k for k in range(0, bunch.mixture.shape[1])}

    tabl = palettable.tableau.Tableau_20.mpl_colors
    Gr_10 = palettable.scientific.sequential.GrayC_20.mpl_colors
    colorlist = [i for i in tabl]

    # font_dir = "/home/goldpm1/miniconda3/envs/cnvpytor/lib/python3.7/site-packages/matplotlib/mpl-data/fonts/ttf/"
    # font_dirs = matplotlib.font_manager.findSystemFonts(fontpaths=font_dir, fontext='ttf')
    # for font in font_dirs:
    #     matplotlib.font_manager.fontManager.addfont(font)
    matplotlib.rcParams["font.family"] = kwargs["FONT_FAMILY"]

    if bunch.includefp == True:
        fp_color_num = samplename_dict[bunch.fp_index]
        colorlist[fp_color_num] = Gr_10[10]

    matplotlib.pyplot.figure(figsize=(6, 6))

    dimensionreduction = ""
    if bunch.mixture.shape[0] > 2:  # More than 3 samples
        dimensionreduction = "SVD"

    if dimensionreduction == "SVD":
        print("SVD → 2D")
        tsvd = TruncatedSVD(n_components=2)
        tsvd.fit(np.concatenate([np_vaf, bunch.mixture]))
        np_vaf = tsvd.transform(np.concatenate(  [np_vaf, bunch.mixture]))[:-NUM_BLOCK]
        mixture = tsvd.transform(np.concatenate( [np_vaf, bunch.mixture]))[-NUM_BLOCK:]
        matplotlib.pyplot.axis([np.min(np_vaf[:, 0]) * 2.1,  np.max(np_vaf[:, 0]) * 2.1,  np.min(np_vaf[:, 1]) * 2.1,  np.max(np_vaf[:, 1]) * 2.1])
        matplotlib.pyplot.xlabel("SVD1")
        matplotlib.pyplot.ylabel("SVD2")
    elif dimensionreduction == "PCA":
        print("PCA → 2D")
        pca = PCA(n_components=2)
        pca.fit(np.concatenate([np_vaf, bunch.mixture]))
        np_vaf = pca.transform(np.concatenate( [np_vaf, bunch.mixture]))[:-NUM_BLOCK]
        mixture = pca.transform(np.concatenate(  [np_vaf, bunch.mixture]))[-NUM_BLOCK:]
        matplotlib.pyplot.axis([np.min(np_vaf[:, 0]) * 2.1,  np.max(np_vaf[:, 0])  * 2.1,  np.min(np_vaf[:, 1]) * 2.1,  np.max(np_vaf[:, 1]) * 2.1])
        matplotlib.pyplot.xlabel("PC1")
        matplotlib.pyplot.ylabel("PC2")
    else:
        matplotlib.pyplot.axis( [0,  np.max(np_vaf[:, :]) * 2.1,  0,  np.max(np_vaf[:, :]) * 2.1] )
        matplotlib.pyplot.xlabel("Feature 1 : VAF x 2 of the Sample 1")
        matplotlib.pyplot.ylabel("Feature 2 : VAF x 2 of the Sample 2")


    matplotlib.pyplot.suptitle("clone{}.{}-{} (postnormalization) : {}".format(kwargs["NUM_CLONE_ITER"], kwargs["TRIAL"], kwargs["STEP_TOTAL"], round(bunch.likelihood)), fontsize=20)
    matplotlib.pyplot.title("sum_child = {}".format( list ( np.round( np.sum (bunch.mixture[:, bunch.makeone_index], axis = 1), 2) ) ), fontsize = 12)

    if bunch.includefp == "Yes":
        fp_color_num = samplename_dict[bunch.fp_index]
        colorlist[fp_color_num] = Gr_10[10]

    
    
    for j in range(0, NUM_CLONE):
        if j == bunch.fp_index:
            for k in bunch.fp_member_index:
                matplotlib.pyplot.scatter(np_vaf[k, 0] * 2, np_vaf[k, 1] * 2, alpha=1, color=Gr_10[10])
        else:
            for k in range(NUM_MUTATION):
                if bunch.membership_p_normalize[k, j] > 0.8:
                    matplotlib.pyplot.scatter( np_vaf[k, 0] * 2, np_vaf[k, 1] * 2, alpha=1, color=[colorlist[samplename_dict[j]]])
                elif bunch.membership_p_normalize[k, j] > 0.1:
                    matplotlib.pyplot.scatter( np_vaf[k, 0] * 2, np_vaf[k, 1] * 2, alpha=bunch.membership_p_normalize[k, j], color=[colorlist[samplename_dict[j]]])
                else: 
                    continue

    for sample_index, sample_key in enumerate(samplename_dict):
        if (sample_index >= bunch.mixture.shape[1]):
            continue

        x_mean = round(bunch.mixture[0][sample_index], 2)
        y_mean = round(bunch.mixture[1][sample_index], 2)

        if x_mean < 0.02:
            matplotlib.pyplot.text(x_mean, y_mean, "{0}".format( [y_mean]), verticalalignment='top', horizontalalignment='right', fontdict = {"fontsize": 14, "fontweight" : "bold"} )    
        elif y_mean < 0.02:
            matplotlib.pyplot.text(x_mean, -0.02, "{0}".format( [x_mean]), verticalalignment='top', horizontalalignment='center',  fontdict = {"fontsize": 14, "fontweight" : "bold"} )    
        else:
            matplotlib.pyplot.text(x_mean, y_mean, "{0}".format( [x_mean, y_mean]), verticalalignment='top', horizontalalignment='center', fontdict = {"fontsize": 16, "fontweight" : "bold"} )

        if (bunch.fp_index == sample_index):
            matplotlib.pyplot.scatter(x_mean, y_mean, marker='s', color= Gr_10[10], edgecolor='black', s=200, label="cluster" + str(sample_index) + " : " + str(list(bunch.membership).count(sample_index)))            
        elif (bunch.makeone_index != []) & (bunch.makeone_index != None):
            if sample_index in bunch.makeone_index:
                matplotlib.pyplot.scatter(x_mean, y_mean, marker='*', color=colorlist[sample_index], edgecolor='black', s=200, label="cluster" + str(sample_index) + " : " + str(list(bunch.membership).count(sample_index)))
            else:
                matplotlib.pyplot.scatter(x_mean, y_mean, marker='s', color=colorlist[sample_index], edgecolor='black', s=100, label="cluster" + str(sample_index) + " : " + str(list(bunch.membership).count(sample_index)))
        else:
            matplotlib.pyplot.scatter(x_mean, y_mean, marker='+', color=colorlist[sample_index], edgecolor='black', s=200, label="cluster" + str(sample_index) + " : " + str(list(bunch.membership).count(sample_index)))

        matplotlib.pyplot.legend()

    if output_filename != "NotSave":
        matplotlib.pyplot.savefig(output_filename)


################################################################################################################################################################################################

def drawfigure_3d(bunch, np_vaf, output_filename, sum_mixture_beforenormalization, **kwargs):
    from mpl_toolkits.mplot3d import Axes3D
    samplename_dict = {k: k for k in range(0, bunch.mixture.shape[1])}

    tabl = palettable.tableau.Tableau_20.mpl_colors
    Gr_10 = palettable.scientific.sequential.GrayC_20.mpl_colors
    colorlist = [i for i in tabl]

    matplotlib.rcParams["font.family"] = kwargs["FONT_FAMILY"]

    #fig, ax = matplotlib.pyplot.subplots(ncols=1, nrows = 1, figsize=(6, 6))
    fig = matplotlib.pyplot.figure( figsize = (6, 6) )
    ax = fig.add_subplot(111, projection='3d')
    #ax.set_axis( [0,  np.max(np_vaf[:, :]) * 2.1, 0,  np.max(np_vaf[:, :]) * 2.1])
    ax.set_xlim(0, np.max(np_vaf[:, :]) * 2.1)
    ax.set_ylim(0, np.max(np_vaf[:, :]) * 2.1)
    ax.set_zlim(0, np.max(np_vaf[:, :]) * 2.1)


    if (kwargs["STEP"] <= (kwargs["COMPULSORY_NORMALIZATION"] - 1)) :
        matplotlib.pyplot.suptitle("clone{}.{}-{} (postnormalization) : {}".format(kwargs["NUM_CLONE_ITER"], kwargs["TRIAL"], kwargs["STEP_TOTAL"], round(bunch.likelihood)), fontsize=20)
    else:
        matplotlib.pyplot.suptitle("clone{}.{}-{} (no normalization) : {}".format(kwargs["NUM_CLONE_ITER"], kwargs["TRIAL"], kwargs["STEP_TOTAL"], round(bunch.likelihood)), fontsize=20)
    matplotlib.pyplot.title("sum_child = {}\tprenorm = {}".format( list (np.round( np.sum (bunch.mixture[:, bunch.makeone_index], axis = 1), 3) ), sum_mixture_beforenormalization), fontsize = 12)


    if (bunch.makeone_index == []) | (bunch.makeone_index == None) | (bunch.checkall_strict == False):  
        matplotlib.pyplot.style.use("Solarize_Light2")
    # elif (bunch.checkall_strict == False):
    #     matplotlib.pyplot.style.use("bmh")
    else:
        matplotlib.pyplot.style.use("seaborn-white")

    if (kwargs["STEP"] <= (kwargs["COMPULSORY_NORMALIZATION"] - 1)) :
        import matplotlib.patches as patches
        fig_width, fig_height = fig.get_size_inches()
        rect = patches.Rectangle((0, 0), fig_width * fig.dpi, fig_height * fig.dpi, linewidth=8, edgecolor="#C23373", facecolor='none')
        fig.patches.append(rect)
        # Set the layout to tight to ensure the rectangle is not cut off
        matplotlib.pyplot.tight_layout()

    
    #matplotlib.pyplot.scatter(np_vaf[:, 0] * 2, np_vaf[:, 1] * 2, color=[colorlist[samplename_dict[k]] for k in bunch.membership])
    for k in range (len ( bunch.membership)):
        if (bunch.fp_index == bunch.membership[k]):
            ax.scatter(np_vaf[k, 0] * 2, np_vaf[k, 1] * 2, np_vaf[k, 2] * 2, color = Gr_10[10]  )
        else:
            ax.scatter(np_vaf[k, 0] * 2, np_vaf[k, 1] * 2, np_vaf[k, 2] * 2, color = colorlist[ samplename_dict [ bunch.membership[k] ]]  )


    for sample_index, sample in enumerate(samplename_dict):
        if (sample_index >= bunch.mixture.shape[1]):
            continue

        # Drawe the centroid based on mixture
        x_mean = round(bunch.mixture[0][sample_index], 2)
        y_mean = round(bunch.mixture[1][sample_index], 2)
        z_mean = round(bunch.mixture[2][sample_index], 2)

        ax.text(x_mean, y_mean, z_mean, "{0}".format( [x_mean, y_mean, z_mean]), verticalalignment='top', horizontalalignment='center', fontdict = {"fontsize": 16, "fontweight" : "bold"} )

        if (bunch.fp_index == sample_index):
            ax.scatter(x_mean, y_mean, z_mean, marker='s', color= Gr_10[10], edgecolor='black', s=200, label="cluster" + str(sample_index) + " : " + str(list(bunch.membership).count(sample_index)))            
        elif (bunch.makeone_index != []) & (bunch.makeone_index != None):
            if sample_index in bunch.makeone_index:
                ax.scatter(x_mean, y_mean, z_mean, marker='*', color=colorlist[sample_index], edgecolor='black', s=400, label="cluster" + str(sample_index) + " : " + str(list(bunch.membership).count(sample_index)))
            else:
                ax.scatter(x_mean, y_mean, z_mean, marker='s', color=colorlist[sample_index], edgecolor='black', s=200, label="cluster" + str(sample_index) + " : " + str(list(bunch.membership).count(sample_index)))
        else:
            ax.scatter(x_mean, y_mean, z_mean, marker='+', color=colorlist[sample_index], edgecolor='black', s=200, label="cluster" + str(sample_index) + " : " + str(list(bunch.membership).count(sample_index)))

        ax.legend()

    if output_filename != "NotSave":
        matplotlib.pyplot.savefig(output_filename)





def drawfigure_3d_SVD(bunch, np_vaf, output_filename, sum_mixture_beforenormalization, **kwargs):
    from mpl_toolkits.mplot3d import Axes3D
    samplename_dict = {k: k for k in range(0, bunch.mixture.shape[1])}

    tabl = palettable.tableau.Tableau_20.mpl_colors
    Gr_10 = palettable.scientific.sequential.GrayC_20.mpl_colors
    colorlist = [i for i in tabl]

    matplotlib.rcParams["font.family"] = kwargs["FONT_FAMILY"]

    fig, ax = matplotlib.pyplot.subplots(ncols=1, nrows = 1, figsize=(6, 6))
    #ax.set_axis( [0,  np.max(np_vaf[:, :]) * 2.1, 0,  np.max(np_vaf[:, :]) * 2.1])
    ax.set_xlim(0, np.max(np_vaf[:, :]) * 2.1)
    ax.set_ylim(0, np.max(np_vaf[:, :]) * 2.1)

    tsvd = TruncatedSVD(n_components=2)
    tsvd.fit( np_vaf )

    np_vaf, mixture = tsvd.transform( np.concatenate(  [np_vaf, bunch.mixture.T]) )[:-kwargs["NUM_CLONE"]], tsvd.transform(np.concatenate( [np_vaf, bunch.mixture.T]) )[-kwargs["NUM_CLONE"]:].T
    #np_vaf = tsvd.transform( np.concatenate(  [np_vaf, bunch.mixture.T]) )
    matplotlib.pyplot.axis([np.min(np_vaf[:, 0]) * 2.1,  np.max(np_vaf[:, 0]) * 2.1,  np.min(np_vaf[:, 1]) * 2.1,  np.max(np_vaf[:, 1]) * 2.1])



    if (kwargs["STEP"] <= (kwargs["COMPULSORY_NORMALIZATION"] - 1)) :
        matplotlib.pyplot.suptitle("clone{}.{}-{} (postnormalization) : {}".format(kwargs["NUM_CLONE_ITER"], kwargs["TRIAL"], kwargs["STEP_TOTAL"], round(bunch.likelihood)), fontsize=20)
    else:
        matplotlib.pyplot.suptitle("clone{}.{}-{} (no normalization) : {}".format(kwargs["NUM_CLONE_ITER"], kwargs["TRIAL"], kwargs["STEP_TOTAL"], round(bunch.likelihood)), fontsize=20)
    matplotlib.pyplot.title("sum_child = {}\tprenorm = {}".format( list (np.round( np.sum (bunch.mixture[:, bunch.makeone_index], axis = 1), 3) ), sum_mixture_beforenormalization), fontsize = 12)


    if (bunch.makeone_index == []) | (bunch.makeone_index == None) | (bunch.checkall_strict == False):  
        matplotlib.pyplot.style.use("Solarize_Light2")
    # elif (bunch.checkall_strict == False):
    #     matplotlib.pyplot.style.use("bmh")
    else:
        matplotlib.pyplot.style.use("seaborn-white")

    if (kwargs["STEP"] <= (kwargs["COMPULSORY_NORMALIZATION"] - 1)) :
        import matplotlib.patches as patches
        fig_width, fig_height = fig.get_size_inches()
        rect = patches.Rectangle((0, 0), fig_width * fig.dpi, fig_height * fig.dpi, linewidth=8, edgecolor="#C23373", facecolor='none')
        fig.patches.append(rect)
        # Set the layout to tight to ensure the rectangle is not cut off
        matplotlib.pyplot.tight_layout()

    
    #matplotlib.pyplot.scatter(np_vaf[:, 0] * 2, np_vaf[:, 1] * 2, color=[colorlist[samplename_dict[k]] for k in bunch.membership])
    for k in range (len ( bunch.membership)):
        if (bunch.fp_index == bunch.membership[k]):
            ax.scatter(np_vaf[k, 0] * 2, np_vaf[k, 1] * 2, color = Gr_10[10]  )
        else:
            ax.scatter(np_vaf[k, 0] * 2, np_vaf[k, 1] * 2, color = colorlist[ samplename_dict [ bunch.membership[k] ]]  )


    for sample_index, sample in enumerate(samplename_dict):
        if (sample_index >= mixture.shape[1]):
            continue

        # Drawe the centroid based on mixture
        x_mean = round(mixture[0][sample_index], 2)
        y_mean = round(mixture[1][sample_index], 2)

        ax.text(x_mean, y_mean, "{0}".format( [ round(bunch.mixture[0][sample_index],2), round(bunch.mixture[1][sample_index],2),  round(bunch.mixture[2][sample_index],2)  ]), verticalalignment='top', horizontalalignment='center', fontdict = {"fontsize": 16, "fontweight" : "bold"} )

        if (bunch.fp_index == sample_index):
            ax.scatter(x_mean, y_mean, marker='s', color= Gr_10[10], edgecolor='black', s=200, label="cluster" + str(sample_index) + " : " + str(list(bunch.membership).count(sample_index)))            
        elif (bunch.makeone_index != []) & (bunch.makeone_index != None):
            if sample_index in bunch.makeone_index:
                ax.scatter(x_mean, y_mean, marker='*', color=colorlist[sample_index], edgecolor='black', s=400, label="cluster" + str(sample_index) + " : " + str(list(bunch.membership).count(sample_index)))
            else:
                ax.scatter(x_mean, y_mean, marker='s', color=colorlist[sample_index], edgecolor='black', s=200, label="cluster" + str(sample_index) + " : " + str(list(bunch.membership).count(sample_index)))
        else:
            ax.scatter(x_mean, y_mean, marker='+', color=colorlist[sample_index], edgecolor='black', s=200, label="cluster" + str(sample_index) + " : " + str(list(bunch.membership).count(sample_index)))

        ax.legend()

    if output_filename != "NotSave":
        matplotlib.pyplot.savefig(output_filename)
