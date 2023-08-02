class simpleKmeansObject:
    def __init__(self, **kwargs):
        self.likelihood_record = [[]] * (kwargs["NUM_CLONE_TRIAL_END"] + 1)
        self.mixture_record = [[]] * (kwargs["NUM_CLONE_TRIAL_END"] + 1)
        self.membership_record = [[]] * (kwargs["NUM_CLONE_TRIAL_END"] + 1)
        self.elbow_K = kwargs["NUM_CLONE_TRIAL_END"]
        self.silhouette_K = kwargs["NUM_CLONE_TRIAL_END"]
        self.gap_K = kwargs["NUM_CLONE_TRIAL_END"]

    def find_max_likelihood (self, start, end):
        import numpy as np
        i = np.argmax(self.likelihood_record [ start : end + 1]) + start
        return i




def decision_elbow ( cluster, **kwargs):
    import numpy as np
    from kneed import KneeLocator
    import matplotlib.pyplot as plt

    if kwargs["NUM_CLONE_TRIAL_START"] == kwargs["NUM_CLONE_TRIAL_END"]:       # 1개밖에 없으면  kneedle이 오류난다
        threshold_x = kwargs["NUM_CLONE_TRIAL_END"]
    else:
        x = np.arange(start = kwargs["NUM_CLONE_TRIAL_START"], stop = kwargs["NUM_CLONE_TRIAL_END"] + 1)
        y = np.array( cluster.likelihood_record  [ kwargs["NUM_CLONE_TRIAL_START"] :  kwargs["NUM_CLONE_TRIAL_END"] + 1 ]  )

        for S in list(np.arange(0.1, 3.0, 0.5)):
            kneedle = KneeLocator(x, y, S=S, curve="convex", direction="decreasing")
            if (kneedle.knee != None) & (kwargs["VERBOSE"] >= 0):
                kneedle.plot_knee()
                plt.savefig ( kwargs["SIMPLE_KMEANS_DIR"] + "/elbow/kneedle." + kwargs["IMAGE_FORMAT"] )
                break
        if kneedle.knee != None:             # 가장 정상적인 경우
            threshold_x = round(kneedle.knee)
        else:                                             # elbow 법에서 제대로 못 찾았을 경우  단순히 가장 높은 값으로 한다
            threshold_x = cluster.find_max_likelihood( kwargs["NUM_CLONE_TRIAL_START"], kwargs["NUM_CLONE_TRIAL_END"] )
            print ("최적의 clone 값을 elbow법으로 찾지 못함 : {}".format(threshold_x))

    return threshold_x   # maxmaxmax_NUM_CLONE


def decision_silhouette (cluster, np_vaf, **kwargs):
    import numpy as np
    from sklearn.metrics import silhouette_samples, silhouette_score

    Silhouette_list = np.zeros (kwargs["NUM_CLONE_TRIAL_END"] + 1, dtype ="float")

    for k in range (kwargs["NUM_CLONE_TRIAL_START"], kwargs["NUM_CLONE_TRIAL_END"] + 1):
        silhouette_score_alldata = silhouette_samples(np_vaf , cluster.membership_record [k] )   # (500,1 )  500
        Silhouette_list [k] = np.mean (silhouette_score_alldata)

    arg_list = [ list(Silhouette_list).index(i) for i in sorted(Silhouette_list, reverse=True)][:2]
    #print ("\tSilhouette_list : {}\narg_list : {}".format  (Silhouette_list, arg_list) )

    return arg_list[0]



def decision_gapstatistics (cluster, np_vaf, **kwargs):
    import pandas as pd
    import numpy as np
    import math, scipy, miscellaneous
    from sklearn.cluster import KMeans

    Input_B = 20
    Gap_list, Std_list, S_list = np.zeros (kwargs["NUM_CLONE_TRIAL_END"] + 1, dtype ="float"), np.zeros (kwargs["NUM_CLONE_TRIAL_END"] + 1, dtype ="float"), np.zeros (kwargs["NUM_CLONE_TRIAL_END"] + 1, dtype ="float")
    Gap_list = np.array ([ float(12345) for i in Gap_list])    
    maxmaxmax_NUM_CLONE  = 0

    
    for NUM_CLONE in range (kwargs["NUM_CLONE_TRIAL_START"], kwargs["NUM_CLONE_TRIAL_END"] + 1):
        if maxmaxmax_NUM_CLONE == 0:
            maxmaxmax_NUM_CLONE = NUM_CLONE
            
        membership = cluster.membership_record [NUM_CLONE]
        mixture = cluster.mixture_record [NUM_CLONE]

        #1. Intra cluster variation (Wk)
        Wk = 0
        for k in range(kwargs["NUM_MUTATION"]):
            j = membership [k]    
            Wk = Wk + math.pow (  scipy.spatial.distance.euclidean(np_vaf[k] * 2, mixture[:, j] ),  2)   # Sum of square 
        Wk = round(math.log10(Wk), 3)
        # if  (kwargs["VERBOSE"] >= 1):
        #     print ("\tMy Clustering\tWk  : {}" .format(Wk))


        #2. Random generation & ICC (Wkb)
        Wkb_list = []
        for b in range (Input_B):
            reference_multiply = miscellaneous.multiply_npvaf ( kwargs["NUM_MUTATION"] , kwargs["NUM_BLOCK"] , np_vaf, 
                                                                                    sorted (set( range(0, kwargs["NUM_MUTATION"]) ))  , b ) 
            reference_np = reference_multiply

            
            kmeans = KMeans(n_clusters=NUM_CLONE, init = cluster.mixture_record [NUM_CLONE].transpose() , max_iter = 10, random_state = 0)  
            kmeans.fit(reference_np)  # nparray
            Wkb_list.append ( round (math.log10(kmeans.inertia_), 3) )
            miscellaneous.drawfigure (reference_np, kmeans.labels_,  kwargs["SIMPLE_KMEANS_DIR"] + "/gap/Kmeans.clone"  +   str (NUM_CLONE) + "." + str(b) + "." + kwargs["IMAGE_FORMAT"], **kwargs)

        Gap_list [NUM_CLONE] = round ( np.mean(Wkb_list) - Wk, 3)
        Std_list [NUM_CLONE] = round ( np.std (Wkb_list), 3)
        S_list [NUM_CLONE] = round (Std_list[NUM_CLONE] * math.sqrt(1 + Input_B) , 3 )

    #     if (kwargs["VERBOSE"] in ["True", "T"]) | (kwargs["VERBOSE"] >= 1):
    #         print ("\tRandom noise (B = {}) : \tmean Wkb = {}\tsdk = {}\tsk (sdk * sqrt ({})) = {}\n\tGap (Wkb - Wk) = {}\n\tPosterior = {}".format (Input_B, round( np.mean(Wkb_list), 3) , Std_list[NUM_CLONE], Input_B + 1, S_list[NUM_CLONE]  ,  Gap_list[NUM_CLONE], round (cluster.likelihood_record[NUM_CLONE]) ))


    # if (kwargs["VERBOSE"] in ["True", "T"]) | (kwargs["VERBOSE"] >= 1):
    #     print ("Gap list : {}\nS list : {}\n".format(Gap_list [kwargs["NUM_CLONE_TRIAL_START"] : kwargs["NUM_CLONE_TRIAL_END"] + 1], S_list [kwargs["NUM_CLONE_TRIAL_START"] : kwargs["NUM_CLONE_TRIAL_END"] + 1] ))


    # Max Gap Number
    Gap_list_index = []
    maxmaxmax_NUM_CLONE  = np.argmax ( Gap_list [kwargs["NUM_CLONE_TRIAL_START"] : kwargs["NUM_CLONE_TRIAL_END"] + 1]  ) + kwargs["NUM_CLONE_TRIAL_START"]
    
    Gap_list_df = pd.DataFrame ( Gap_list  ).sort_values ( by = 0, ascending = False)
    
    Gap_list_df = Gap_list_df [ Gap_list_df[0] != float(12345)]
    Gap_list = Gap_list_df[0]   # because of 12345

    for i in range( len(Gap_list_df) ):
        #print ("Gap statistics method (max Gap): {}th optimal NUM_CLONE = {}".format(i  + 1, Gap_list_df.index[i]  ))
        Gap_list_index.append ( Gap_list_df.index[i]  )

    return Gap_list_index [0]


            




def clustering (np_vaf, **kwargs):
    from sklearn.cluster import KMeans
    import numpy as np

    simpleK =  simpleKmeansObject ( **kwargs)

    for k in range (kwargs ["NUM_CLONE_TRIAL_START"], kwargs["NUM_CLONE_TRIAL_END"] + 1):
        kmeans = KMeans(n_clusters = k, init='k-means++', max_iter=100,random_state=0)  # model generation
        kmeans.fit (np_vaf)  
        simpleK.membership_record [k] = kmeans.labels_ 
        simpleK.likelihood_record [k] = kmeans.inertia_ 
        simpleK.mixture_record [k] = kmeans.cluster_centers_.T * 2

    simpleK.elbow_K = decision_elbow (simpleK, **kwargs)
    simpleK.silhouette_K = decision_silhouette (simpleK, np_vaf, **kwargs)
    simpleK.gap_K  =  decision_gapstatistics (simpleK, np_vaf, **kwargs)

    return (kwargs, simpleK)


def visualization (simpleK, np_vaf, **kwargs):
    import visualizationsingle
    import numpy as np

    if kwargs ["NUM_BLOCK"] == 1:
        visualizationsingle.drawfigure_1d ( membership = simpleK.membership_record [simpleK.elbow_K],
                                                                output_suptitle = "simpleKmeans_Elbow",
                                                                output_filename = kwargs["SIMPLE_KMEANS_DIR"] + "/elbow/simpleKmeans_elbow." + kwargs["IMAGE_FORMAT"],
                                                                np_vaf = np_vaf,
                                                                samplename_dict = {k:"clone {}".format(k) for k in range(0, np.max( simpleK.membership_record [simpleK.elbow_K] ) + 1)},
                                                                includefp = False,
                                                                fp_index = -1,
                                                                makeone_index = [],
                                                                **kwargs)
    elif kwargs ["NUM_BLOCK"] == 2:
        visualizationsingle.drawfigure_2d ( membership = simpleK.membership_record [simpleK.elbow_K],
                                                                output_suptitle = "simpleKmeans_Elbow",
                                                                output_filename = kwargs["SIMPLE_KMEANS_DIR"] + "/elbow/simpleKmeans_elbow." + kwargs["IMAGE_FORMAT"],
                                                                np_vaf = np_vaf,
                                                                samplename_dict = {k:"clone {}".format(k) for k in range(0, np.max( simpleK.membership_record [simpleK.elbow_K] ) + 1)},
                                                                includefp = False,
                                                                fp_index = -1,
                                                                dimensionreduction = "None"
                                                                **kwargs)
    elif kwargs ["NUM_BLOCK"] >= 3:
        visualizationsingle.drawfigure_2d ( membership = simpleK.membership_record [simpleK.elbow_K],
                                                                output_suptitle = "simpleKmeans_Elbow",
                                                                output_filename = kwargs["SIMPLE_KMEANS_DIR"] + "/elbow/simpleKmeans_elbow." + kwargs["IMAGE_FORMAT"],
                                                                np_vaf = np_vaf,
                                                                samplename_dict = {k:"clone {}".format(k) for k in range(0, np.max( simpleK.membership_record [simpleK.elbow_K] ) + 1)},
                                                                includefp = False,
                                                                fp_index = -1,
                                                                dimensionreduction = "SVD"
                                                                **kwargs)



def scoring (membership_answer, membership_answer_numerical, simpleK, **kwargs):
    import scoring

    simpleK.elbow_K_score, sample_dict_PtoA, sample_dict_AtoP  = scoring.Scoring ( membership_answer, membership_answer_numerical,
                                                                                                                                            simpleK.membership_record [simpleK.elbow_K], -1 , [] ) # fp를 designate 하지 못하니까 무조건 fp_index는 -1, parent_index는 []
    simpleK.silhouette_K_score, sample_dict_PtoA, sample_dict_AtoP  = scoring.Scoring ( membership_answer, membership_answer_numerical,
                                                                                                                                            simpleK.membership_record [simpleK.silhouette_K], -1 , [] ) # fp를 designate 하지 못하니까 무조건 fp_index는 -1, parent_index는 []
    simpleK.gap_K_score, sample_dict_PtoA, sample_dict_AtoP  = scoring.Scoring ( membership_answer, membership_answer_numerical,
                                                                                                                                            simpleK.membership_record [simpleK.gap_K], -1 , [] ) # fp를 designate 하지 못하니까 무조건 fp_index는 -1, parent_index는 []

    return simpleK