def np_vaf_extract(df):
    import numpy as np

    NUM_BLOCK = len(df[0])
    NUM_MUTATION = len(df)

    np_vaf = np.zeros((NUM_MUTATION, NUM_BLOCK), dtype = 'float')
    for row in range (NUM_MUTATION):
        for col in range (NUM_BLOCK):
            if df[row][col]["depth"] == 0:
                np_vaf[row][col] = 0
            else:    
                np_vaf[row][col] = round (df[row][col]["alt"] / df[row][col]["depth"] , 3)

    return np_vaf

def iszerocolumn (step, **kwargs):
    import numpy as np 

    for j in range ( kwargs ["NUM_CLONE_ITER"] ) :
        if  np.array_equal (step.mixture[: , j],  np.zeros ( kwargs["NUM_BLOCK"]) )  == True :
            #print ("\t\tEmpty clone found.\t{}".format (  step.mixture  ))
            return True
    return False

##################################################################################################################################################################################################################################################

def checkall (step, condition, **kwargs):   # condition = "strict" or "lenient"
    import numpy as np

    if type (step) == type ([]):     # from isparent.py
        sum_mixture = step

    else:  # from Mstep.py
        sum_mixture = np.zeros ( kwargs["NUM_BLOCK"], dtype = "float")
        for i in range (kwargs["NUM_BLOCK"]):
            for j in range (kwargs ["NUM_CLONE_ITER"]):
                if j in step.makeone_index:
                    sum_mixture[i] += step.mixture[i][j]


    if condition == "lenient":       # 이때는 좀 널널하게 잡는다. 그래도 안되면 중간에 stop하게 시킴
        if kwargs["MAKEONE_STRICT"] == 1:  # SimData, CellData
            makeone_standard = np.array ( [ [0.84, 1.15], [0.82, 1.3], [0.82, 1.3]  ],dtype = float)   # 1st: 1D,  2nd : 2D, 3D
        elif kwargs["MAKEONE_STRICT"] == 2: #  SimData, CellData Low depth
            makeone_standard = np.array ( [ [0.84, 1.3], [0.82, 1.3], [0.82, 1.3] ],dtype = float)
        elif kwargs["MAKEONE_STRICT"] == 3: #  BioData
            if kwargs["NUM_CLONE_ITER"] == 1:
                makeone_standard = np.array ( [ [0.80, 1.27], [0.77, 1.3], [0.77, 1.3] ],dtype = float)      # monoclonal (liver, clone) 이면  좀더 널널하게 잡음
            elif kwargs["NUM_CLONE_ITER"]  > 1:
                makeone_standard = np.array ( [ [0.92, 1.27], [0.77, 1.3], [0.77, 1.3] ],dtype = float)      # clone number가 좀더 많으면 빡빡하게 잡음
    elif condition == "strict":   # 이때는 조금 더 빡빡하게 잡음
        if kwargs["MAKEONE_STRICT"] == 1:    # SimData, CellData
            makeone_standard = np.array ( [ [0.95, 1.05], [0.95, 1.06], [0.95, 1.06] ],dtype = float)  # 1st: 1D,  2nd : 2D, 3rd : 3D
        elif kwargs["MAKEONE_STRICT"] == 2:   # SimData, CellData Low depth
            makeone_standard = np.array ( [ [0.95, 1.20], [0.94, 1.3], [0.93, 1.3] ],dtype = float)  # 1st: 1D,  2nd : 2D, 3D
        elif kwargs["MAKEONE_STRICT"] == 3:   # BioData
            makeone_standard = np.array ( [ [0.85, 1.18], [0.85, 1.3], [0.85, 1.3] ],dtype = float)  # 1st: 1D,  2nd : 2D, 3D


    if type (step) != type ([]):  # from Mstep.py
        if len(step.makeone_index) == 1:  # If monoclonal, set extremely lenient condition (due to the homologous variant contam).
            makeone_standard = np.array ( [ [0.7, 1.3], [0.7, 1.3] ],dtype = float)  # 1st: 1D,  2nd : 2D, 3D


    if kwargs["NUM_BLOCK"] == 1:      # 1D
        if (sum_mixture[0] < makeone_standard[0][0]) | (sum_mixture[0] > makeone_standard[0][1]):
            return False, sum_mixture
        else:
            return True, sum_mixture
    elif kwargs["NUM_BLOCK"] == 2:      # 2D
        for i in range( kwargs["NUM_BLOCK"] ):  
            if (sum_mixture[i] < makeone_standard[1][0]) | (sum_mixture[i] > makeone_standard[1][1]): 
                return False, sum_mixture
        return True, sum_mixture
    elif kwargs["NUM_BLOCK"] == 3:      # 3D
        for i in range( kwargs["NUM_BLOCK"] ):  
            if (sum_mixture[i] < makeone_standard[2][0]) | (sum_mixture[i] > makeone_standard[2][1]): 
                return False, sum_mixture
        return True, sum_mixture
    



##################################################################################################################################################################################################################################################

def DeleteCentroid (membership_kmeans, mixture_kmeans, **kwargs):
    import numpy as np

    ################Delete n (membership) less than MIN_CLUSTER_SIZE################
    mask = []
    t =   np.unique (membership_kmeans, return_counts = True ) 
    print ( "np.unique (return_counts = True)  :  {}".format( t ) )
    for j in range ( len( t[0] ) ):
        if t[1][j] >= kwargs["MIN_CLUSTER_SIZE"]:        # over MIN_CLUSTER_SIZE
            mask.append (j)
    mixture_kmeans = mixture_kmeans [:, mask]
    ###########################################################################


    ######################### Delete centroid which value is over 1 ####################
    mask = []
    for j in range ( mixture_kmeans.shape[1] ):
        if ( np.any ( mixture_kmeans [ : , j] > 1)  == False):   # under 1
            mask.append (j)
    mixture_kmeans =  mixture_kmeans [ :, mask ]
    ###########################################################################

    return mixture_kmeans



def initial_kmeans (input_containpos, df, np_vaf, np_BQ, OUTPUT_FILENAME, **kwargs):
    import numpy as np
    import matplotlib.pyplot as plt
    from sklearn.cluster import KMeans
    import palettable, itertools, random, re, math, copy, subprocess
    import seaborn as sns
    from scipy.stats import kde
    import clementrec

    NUM_BLOCK = np_vaf.shape[1]
    NUM_MUTATION = np_vaf.shape[0]


    
    ##############################  Add axis centroid ##############################    

    if kwargs["NUM_BLOCK"] == 1:
        kmeans = KMeans(n_clusters=kwargs["KMEANS_CLUSTERNO"], init = 'k-means++', max_iter=100,random_state=0)  # model generation
        kmeans.fit(np_vaf)  
        membership_kmeans = kmeans.labels_         
        
        if kwargs["VERBOSE"] >= 1:
            print ( np.unique (membership_kmeans, return_counts=True)[1]  )
            print ( np.where ( np.unique (membership_kmeans, return_counts=True)[1] != 0 )[0] )

        #mixture_kmeans = np.zeros (( NUM_BLOCK, kwargs["KMEANS_CLUSTERNO"]), dtype = "float")
        mixture_kmeans = np.zeros ( ( NUM_BLOCK, len ( np.where ( np.unique (membership_kmeans, return_counts=True)[1] != 0 )[0] )  ) , dtype = "float")

        for j_index, j in enumerate ( np.where ( np.unique (membership_kmeans, return_counts=True)[1] != 0 )[0] ):
            for i in range (NUM_BLOCK):
                mixture_kmeans[i][j_index] = round(np.mean(np_vaf[membership_kmeans == (j)][:, i] * 2), 3)



    if kwargs["NUM_BLOCK"] == 2:
        combi =  itertools.chain(*map(lambda x: itertools.combinations( list (range (kwargs["NUM_BLOCK"])), x), range(1, len( list (range (  kwargs["NUM_BLOCK"])) ) )))

        # 빈 mixture 만들어주기
        mixture_kmeans = np.zeros ( (kwargs["NUM_BLOCK"], 0) , dtype = "float") 

        # 축 혹은 평면에서 CLEMENT 돌리기 (clementrec)
        for subset in combi:    # 0인 평면 (3D) 혹은 축 (3D, 2D)
            subdim = list (subset)       # [0,]  [1,]  [2,] [ 0, 1] [ 0, 2] [1, 2]
            nonzero_dim = sorted(list ( set (range (kwargs["NUM_BLOCK"])) - set(subdim) ))
            
            index_interest = set ( range (len(np_vaf)))
            for zero_dim in subdim:   # 그 평면 혹은 축에 있는 index만 고름
                index_interest = set(index_interest) & set ( np.where ( np_vaf[: , zero_dim] == 0 )[0] )
                index_interest = sorted ( list (index_interest))

            if len(index_interest) == 0:  # 축 상에 없으면 넘어가도 됨
                print ( "\n➜➜➜\nsubdim = {} ( {} = 0 ) \tnonzero_dim = {}\tPASS".format ( subdim,  [chr (k+120)   for k in subdim], nonzero_dim ))
                continue
            

            
            index_interest_nonzero =  np.where ( np.all ( np_vaf[: , nonzero_dim] != 0, axis = 1 )  ) [0]    # 그 축이 아닌 다른 축에서 0이라면 뺴자
            #index_interest_nonzero = np.arange ( kwargs["NUM_MUTATION"] )   # 그냥 싹 다 사영 내리는 것

            # 최대 몇 개 뽑을지 결정하기
            tt = []
            for i in nonzero_dim : 
                t = round ( np.mean (  np_vaf[ index_interest_nonzero , i]) , 2)
                if t <= 0.13:
                    tt.append (4)
                else:
                    tt.append (  math.ceil ( 0.5 / t  )  )   # 그 축에서 몇 개까지의 clone이 가능할 것인가 (올림하자)
            max_t = np.max(tt)
            select_t = np.min ( [max_t, math.ceil ( len(index_interest) / kwargs["MIN_CLUSTER_SIZE"]) ] )   # 한 clustert당 이만큼은 들어가야 하니..
            if (select_t == 0) | ( len(index_interest) < kwargs["MIN_CLUSTER_SIZE"] ):     # 그 축에 cluster가 있는게 말이 안 될때에는 그냥 넘어가자
                print ( "\n\n\n\n➜➜➜\nsubdim = {} ( {} = 0 ) \tnonzero_dim = {}\tnum_mutation on the plane = {}\tnum_mut not exclusive on the plane = {}\tmean VAF*2 = {} -> select_t = {}\tso PASS".format (  subdim, [chr (k+120)   for k in subdim], nonzero_dim, len ( index_interest), len(index_interest_nonzero), t * 2, select_t ))
                continue
            else:
                print ( "\n\n\n\n➜➜➜\nsubdim = {} ( {} = 0 ) \tnonzero_dim = {}\tnum_mutation on the plane = {}\tnum_mut not exclusive on the plane = {}\tmean VAF*2 = {} -> select_t = {}".format (  subdim, [chr (k+120)   for k in subdim], nonzero_dim, len ( index_interest), len(index_interest_nonzero), t * 2, select_t ))


            # zero plane만 따로 뽑아서 CLEMENT를 recursive하게 돌려주기
            kwargs_transfer = copy.deepcopy ( kwargs ) 
            kwargs_transfer["NUM_CLONE_TRIAL_START"] = np.min ([select_t, 1])          
            kwargs_transfer["NUM_CLONE_TRIAL_END"] = np.min ([select_t, 4])
            kwargs_transfer["FP_RATIO"] = 0
            kwargs_transfer["NUM_PARENT"] = 0
            kwargs_transfer["NUM_BLOCK"] = len(nonzero_dim)
            kwargs_transfer["TRIAL_NO"] = 3
            kwargs_transfer["VERBOSE"] = 0
            input_containpos = input_containpos.reset_index(drop = True) 
            

            # 정사영 내린 것만 보기 위해 dimension selection
            NUM_CLONE_recursive, mixture_recursive = clementrec.recursive  ( input_containpos = input_containpos , 
                                                    df = [ [row[col] for col in nonzero_dim] for row in df  ],
                                                    np_vaf = np_vaf[:,  nonzero_dim ],
                                                    np_BQ = np_BQ[:, nonzero_dim ],
                                                    subdim  = subdim, nonzero_dim = nonzero_dim, 
                                                    index_interest = index_interest, index_interest_nonzero = index_interest_nonzero,
                                                    kwargs = kwargs_transfer)
    

            # 원래 차원으로 회복시키고 mixture를 합쳐주기
            if np.all(mixture_recursive == 0) != True:     # 전체가 (0,0)이면 무시하자 (FP는 어차피 나중에 한번에 붙여줌)
                mixture_subdim = np.zeros ( (kwargs["NUM_BLOCK"], NUM_CLONE_recursive) , dtype = "float") 
                j = 0
                for i in range ( kwargs["NUM_BLOCK"]) :
                    if i in nonzero_dim:
                        mixture_subdim[ i, : ] = mixture_recursive [ j , :]
                        j = j + 1
                
                #print ("\tmixture_subdim = {}".format ( mixture_subdim ))
                mixture_kmeans = np.hstack(( mixture_kmeans, mixture_subdim )) 
                print ("\tmixture_kmeans = {}".format (  ", ".join(str(np.round(row, 2)) for row in mixture_kmeans  )  ))

            # 흔적 지워주기
            subprocess.run (["rm -rf " + kwargs["CLEMENT_DIR"] + "/trial/*"  ], shell = True)
            subprocess.run (["rm -rf " + kwargs["CLEMENT_DIR"] + "/candidate/*"  ], shell = True)



        # 공간 + 평면에서 K means 돌리기
        if kwargs["KMEANS_CLUSTERNO"] - mixture_kmeans.shape[1] <= 1:
            kwargs["KMEANS_CLUSTERNO"] += 1      # 공간상에 1개밖에 안남았다면, +1 시켜주자

        print ("\n\n\n➜➜➜➜➜➜➜➜ SPACE k means")    

        # 0 없는 평면에서 K means 돌리기
        non_zero_rows = np.all( np_vaf != 0, axis = 1 )  # 아예 space
        non_zero_np_vaf = np_vaf [non_zero_rows, : ]

        if (kwargs["KMEANS_CLUSTERNO"] - mixture_kmeans.shape[1] + 1) * kwargs ["MIN_CLUSTER_SIZE"] >= non_zero_np_vaf.shape[0]:     # 남은  cluster 개수가 과하다고 생각될 때
            kwargs["KMEANS_CLUSTERNO"] = mixture_kmeans.shape[1] + math.floor (non_zero_np_vaf.shape[0] / kwargs ["MIN_CLUSTER_SIZE"]  ) - 1
        print ( "전체 KMEANS_CLUSTERNO = {}\t현재 차지하고 있는 centroid = {}\tMIN_CLUSTER_SIZE = {}\tspace 상의 variant 개수 = {}".format ( kwargs["KMEANS_CLUSTERNO"], mixture_kmeans.shape[1] , kwargs ["MIN_CLUSTER_SIZE"], non_zero_np_vaf.shape[0]   )  )

        kmeans = KMeans(n_clusters = np.max ( [ kwargs["KMEANS_CLUSTERNO"] - mixture_kmeans.shape[1], 1 ] )  , init='k-means++', max_iter=100, random_state=0)  # model generation
        kmeans.fit ( non_zero_np_vaf )  
        membership_kmeans = kmeans.labels_     

        non_zeroplane_mixture = kmeans.cluster_centers_.T * 2
        print ("non_zeroplane_mixture (before removal) : \n{}".format ( np.round(non_zeroplane_mixture, 2) ))
        non_zeroplane_mixture = DeleteCentroid (membership_kmeans, non_zeroplane_mixture, **kwargs)     # space에서 기준 미달들만 제거
        mixture_kmeans = np.hstack(( mixture_kmeans, non_zeroplane_mixture ))    

        print ("non_zeroplane_mixture (after removal): \n{}".format ( np.round(non_zeroplane_mixture, 2) ))
        
        #mixture_kmeans = np.array ( [[0.18, 0,   0.86, 0.66, 0,   0.19, 0.2,  0.96], [0,   0.09, 0.12, 0.13, 0.2,  0.53, 0.72, 0.81] ] )
        #mixture_kmeans = np.array ( [ [0.06, 0.,   0.,   0.18, 0.,   0.74, 0.9 ], [0.,   0.08, 0.16, 0.18, 0.22, 0.53, 0.84] ]) 


        # 둘을 합쳐주기
        kwargs["KMEANS_CLUSTERNO"] = mixture_kmeans.shape[1]
        print ("\n\n\n\n\n")



    elif kwargs["NUM_BLOCK"] == 3:  # 3D 일때    (231014 : 불편하니까 일단 축만 뽑자)
        combi =  itertools.chain(*map(lambda x: itertools.combinations( list (range (3)), x), [2]  ))

        # 빈 mixture 만들어주기
        mixture_kmeans = np.zeros ( (kwargs["NUM_BLOCK"], 0) , dtype = "float") 

        for subset in combi:    # 0인 축 (3D)
            subdim = list (subset)       # [0, 1], [0, 2], [1, 2]
            nonzero_dim = sorted(list ( set (range (kwargs["NUM_BLOCK"])) - set(subdim) ))   # [2], [1], [0]

            index_interest = set ( range (len(np_vaf)))
            for zero_dim in subdim:   # 그 평면 혹은 축에 있는 index만 고름
                index_interest = set(index_interest) & set ( np.where ( np_vaf[: , zero_dim] == 0 )[0] )
                index_interest = sorted ( list (index_interest))

            if len(index_interest) == 0:  # 축 상에 하나도 없으면 넘어가도 됨
                print ( "\n➜➜➜\nsubdim = {} ( {} = 0 ) \tnonzero_dim = {}\tPASS".format ( subdim,  [chr (k+120)   for k in subdim], nonzero_dim ))
                continue
            
            index_interest_nonzero =  np.where ( np.all ( np_vaf[: , nonzero_dim] != 0, axis = 1 )  ) [0]    # 그 축이 아닌 다른 축에서 0이라면 뺴자
            #index_interest_nonzero = np.arange ( kwargs["NUM_MUTATION"] )   # 그냥 싹 다 사영 내리는 것

            # 최대 몇 개 뽑을지 결정하기
            tt = []
            for i in nonzero_dim : 
                t = round ( np.mean (  np_vaf[ index_interest_nonzero , i]) , 2)
                if t <= 0.13:
                    tt.append (4)
                else:
                    tt.append (  math.ceil ( 0.5 / t  )  )   # 그 축에서 몇 개까지의 clone이 가능할 것인가 (올림하자)
            max_t = np.max(tt)
            select_t = np.min ( [max_t, math.ceil ( len(index_interest) / kwargs["MIN_CLUSTER_SIZE"]) ] )   # 한 clustert당 이만큼은 들어가야 하니..
            if (select_t == 0) | ( len(index_interest) < kwargs["MIN_CLUSTER_SIZE"]):     # 그 축에 cluster가 있는게 말이 안 될때에는 그냥 넘어가자
                print ( "\n\n\n\n➜➜➜\nsubdim = {} ( {} = 0 ) \tnonzero_dim = {}\tnum_mutation on the plane = {}\tnum_mut not exclusive on the plane = {}\tmean VAF*2 = {} -> select_t = {}\tso PASS".format (  subdim, [chr (k+120)   for k in subdim], nonzero_dim, len ( index_interest), len(index_interest_nonzero), t * 2, select_t ))
                continue
            else:
                print ( "\n\n\n\n➜➜➜\nsubdim = {} ( {} = 0 ) \tnonzero_dim = {}\tnum_mutation on the plane = {}\tnum_mut not exclusive on the plane = {}\tmean VAF*2 = {} -> select_t = {}".format (  subdim, [chr (k+120)   for k in subdim], nonzero_dim, len ( index_interest), len(index_interest_nonzero), t * 2, select_t ))


            # zero plane만 따로 뽑아서 CLEMENT를 recursive하게 돌려주기
            kwargs_transfer = copy.deepcopy ( kwargs ) 
            kwargs_transfer["NUM_CLONE_TRIAL_START"] = np.min ([select_t, 1])          
            kwargs_transfer["NUM_CLONE_TRIAL_END"] = 4
            kwargs_transfer["FP_RATIO"] = 0
            kwargs_transfer["NUM_PARENT"] = 0
            kwargs_transfer["NUM_BLOCK"] = len(nonzero_dim)
            kwargs_transfer["TRIAL_NO"] = 3
            kwargs_transfer["VERBOSE"] = 0
            input_containpos = input_containpos.reset_index(drop = True) 
            

            # 정사영 내린 것만 보기 위해 dimension selection
            NUM_CLONE_recursive, mixture_recursive = clementrec.recursive  ( input_containpos = input_containpos , 
                                                    df = [ [row[col] for col in nonzero_dim] for row in df  ],
                                                    np_vaf = np_vaf[:,  nonzero_dim ],
                                                    np_BQ = np_BQ[:, nonzero_dim ],
                                                    subdim  = subdim, nonzero_dim = nonzero_dim, 
                                                    index_interest = index_interest, index_interest_nonzero = index_interest_nonzero,
                                                    kwargs = kwargs_transfer)
    

            # 원래 차원으로 회복시키고 mixture를 합쳐주기
            if np.all(mixture_recursive == 0) != True:     # 전체가 (0,0)이면 무시하자
                mixture_subdim = np.zeros ( (kwargs["NUM_BLOCK"], NUM_CLONE_recursive) , dtype = "float") 
                j = 0
                for i in range ( kwargs["NUM_BLOCK"]) :
                    if i in nonzero_dim:
                        mixture_subdim[ i, : ] = mixture_recursive [ j , :]
                        j = j + 1
                
                #print ("\tmixture_subdim = {}".format ( mixture_subdim ))
                mixture_kmeans = np.hstack(( mixture_kmeans, mixture_subdim )) 
                print ("\t(final) mixture_kmeans = {}".format (  ", ".join(str(np.round(row, 2)) for row in mixture_kmeans  )  ))


            # 흔적 지워주기
            subprocess.run (["rm -rf " + kwargs["CLEMENT_DIR"] + "/trial/*"  ], shell = True)
            subprocess.run (["rm -rf " + kwargs["CLEMENT_DIR"] + "/candidate/*"  ], shell = True)
            

        # 공간 + 평면에서 K means 돌리기
        if kwargs["KMEANS_CLUSTERNO"] - mixture_kmeans.shape[1] <= 1:
            kwargs["KMEANS_CLUSTERNO"] += 1      # 공간상에 1개밖에 안남았다면, +1 시켜주자

        print ("\n\n\n➜➜➜➜➜➜➜➜ SPACE k means")    

        kmeans = KMeans(n_clusters = np.max ( [ kwargs["KMEANS_CLUSTERNO"] - mixture_kmeans.shape[1], 1 ] )  , init='k-means++', max_iter=100, random_state=0)  # model generation
        #non_zero_rows = np.all( np_vaf != 0, axis = 1 )  # 아예 space
        non_zero_rows = np.sum(np_vaf != 0, axis=1) >= 2  # space + 평면
        non_zero_np_vaf = np_vaf [non_zero_rows, : ]
        print ("len(non_zero_rows) = {}".format ( non_zero_rows.sum() ) )
        kmeans.fit ( non_zero_np_vaf )  
        membership_kmeans = kmeans.labels_     

        non_zeroplane_mixture = kmeans.cluster_centers_.T * 2
        print ("non_zeroplane_mixture (before removal) : \n{}".format ( ", ".join(str(np.round(row, 2)) for row in non_zeroplane_mixture  ) ) )

        non_zeroplane_mixture = DeleteCentroid (membership_kmeans, non_zeroplane_mixture, **kwargs)     # space에서 기준 미달들만 제거
        mixture_kmeans = np.hstack(( mixture_kmeans, non_zeroplane_mixture ))    
        print ("non_zeroplane_mixture (after removal) : \n{}".format ( ", ".join(str(np.round(row, 2)) for row in non_zeroplane_mixture  ) ) )
        
        #mixture_kmeans = np.array ( [[0.18, 0,   0.86, 0.66, 0,   0.19, 0.2,  0.96], [0,   0.09, 0.12, 0.13, 0.2,  0.53, 0.72, 0.81] ] )
        #mixture_kmeans = np.array ( [ [0.06, 0.,   0.,   0.18, 0.,   0.74, 0.9 ], [0.,   0.08, 0.16, 0.18, 0.22, 0.53, 0.84] ]) 


        # 둘을 합쳐주기
        kwargs["KMEANS_CLUSTERNO"] = mixture_kmeans.shape[1]
        print ("\n# of KMEANS_CLUSTERNO = {}\n\n\n\n\n".format (kwargs["KMEANS_CLUSTERNO"]))




    ############################################################################
    # Use lexsort to sort the array by the third row, then the second row, and finally the first row
    sorted_indices = np.lexsort(tuple(mixture_kmeans))
    # Use the sorted indices to rearrange the original array
    mixture_kmeans = mixture_kmeans[:, sorted_indices]
    ############################################################################


    if kwargs["VERBOSE"] >= 1:
        print ("\n► All starting points (by K means)" )
        print("\n".join(str(np.round(row, 2)) for row in mixture_kmeans ))





    tabl = palettable.tableau.Tableau_20.mpl_colors
    Gr_10 = palettable.scientific.sequential.GrayC_20.mpl_colors
    colorlist = [i for i in tabl]
    plt.figure(figsize=(6, 6))
            
    if NUM_BLOCK == 1:
        sns.kdeplot ( np_vaf [:, 0] * 2, shade = True)
        kde_np_vaf = kde.gaussian_kde(np_vaf[:, 0] * 2)
        plt.xlabel("VAF x 2 of the Sample 1",  fontdict={"fontsize": 14})
        plt.ylabel("Density", fontdict={"fontsize": 14})
        plt.xlim(0, 1)
        for j in range( mixture_kmeans.shape[1] ):
            plt.axvline ( x = mixture_kmeans[0][j], color = colorlist [j  % 20] ) 
            plt.text  ( mixture_kmeans[0][j], kde_np_vaf(mixture_kmeans[0][j]) * 1.08, "{}".format( np.round ( mixture_kmeans[0][j] , 2)), verticalalignment='top', ha = "center", fontdict = {"fontsize": 16, "fontweight" : "bold"}  )
        y_min, y_max = plt.ylim()
        plt.suptitle ("INITIAL_KMEANS", fontsize =  26, fontweight='semibold' )
        plt.text ( 0.5, y_max - 0.05 * (y_max - y_min), "KMEANS_CLUSTERNO = {}".format ( mixture_kmeans.shape[1] ), verticalalignment='top', ha = "center", fontdict = {"fontsize": 16, "fontweight" : "bold"}  )
        plt.text ( 0.5, y_max - 0.1 * (y_max - y_min), "AXIS_RATIO = {}".format (round (kwargs["AXIS_RATIO"], 2) ), verticalalignment='top', ha = "center", fontdict = {"fontsize": 16, "fontweight" : "bold"}  )
        plt.savefig ( OUTPUT_FILENAME )
    elif NUM_BLOCK == 2:
        plt.xlabel("VAF x 2 of the Sample 1",  fontdict={"fontsize": 14})
        plt.ylabel("VAF x 2 of the Sample 2", fontdict={"fontsize": 14})
        plt.xlim(-0.02, 1);   plt.ylim(-0.02, 1)
        for k in range ( np_vaf.shape[0] ):
            plt.scatter ( x = np_vaf [k, 0] * 2, y = np_vaf [k, 1] * 2, s = 30, color = "#EAC696" , alpha = 0.8)
        plt.scatter ( x = 0, y = 0, s = 30, color = Gr_10[10] , alpha = 0.8)        
        for j in range(kwargs["KMEANS_CLUSTERNO"]):     # Centroid 찍기
            plt.scatter ( mixture_kmeans[0][j], mixture_kmeans[1][j], marker = '*', color = colorlist[j % 20], edgecolor = "black", s = 500, label = "clone {} : {}".format (j, mixture_kmeans[:, j] ) )
            plt.text  ( mixture_kmeans[0][j], mixture_kmeans[1][j] -0.04, "[{},{}]".format( np.round (mixture_kmeans[0][j] , 2), np.round (mixture_kmeans[1][j], 2)), verticalalignment='top', ha = "center", fontdict = {"fontsize": 16, "fontweight" : "bold"}  )
        y_min, y_max = plt.ylim()
        plt.suptitle ("INITIAL_KMEANS", fontsize =  26, fontweight='semibold' )
        plt.text ( 0.5, y_max - 0.05 * (y_max - y_min), "KMEANS_CLUSTERNO = {}".format (kwargs["KMEANS_CLUSTERNO"]), verticalalignment='top', ha = "center", fontdict = {"fontsize": 16, "fontweight" : "bold"}  )
        plt.text ( 0.5, y_max - 0.1 * (y_max - y_min), "AXIS_RATIO = {}".format ( round(kwargs["AXIS_RATIO"], 2) ), verticalalignment='top', ha = "center", fontdict = {"fontsize": 16, "fontweight" : "bold"}  )
        plt.savefig ( OUTPUT_FILENAME )
    elif NUM_BLOCK == 3:
        # from mpl_toolkits.mplot3d import Axes3D 
        # fig = plt.figure( figsize=(6, 6) )
        # ax = fig.add_subplot(111, projection='3d')
        # ax.set_xlabel("VAF x 2 of the Sample 1",  fontdict={"fontsize": 14})
        # ax.set_ylabel("VAF x 2 of the Sample 2", fontdict={"fontsize": 14})
        # ax.set_zlabel("VAF x 2 of the Sample 3", fontdict={"fontsize": 14})
        # ax.set_xlim(-0.02, 1); ax.set_ylim(-0.02, 1); ax.set_zlim(-0.02, 1)
        # for k in range ( np_vaf.shape[0] ):
        #     ax.scatter ( np_vaf [k, 0] * 2, np_vaf [k, 1] * 2, np_vaf [k, 2] * 2, s = 30, color = "#EAC696" , alpha = 0.8)
        # ax.scatter ( 0, 0, 0, s = 30, color = Gr_10[10] , alpha = 0.8)        
        # for j in range(kwargs["KMEANS_CLUSTERNO"]):     # Centroid 찍기
        #     ax.scatter ( mixture_kmeans[0][j], mixture_kmeans[1][j], mixture_kmeans[2][j], marker = '*', color = colorlist[j % 20], edgecolor = "black", s = 500, label = "clone " + str(j))  
        #     ax.text  ( mixture_kmeans[0][j], mixture_kmeans[1][j] , mixture_kmeans[2][j] , "[{},{},{}]".format( np.round (mixture_kmeans[0][j] , 2), np.round (mixture_kmeans[1][j], 2), np.round (mixture_kmeans[2][j], 2) ), verticalalignment='top', ha = "center", fontdict = {"fontsize": 16, "fontweight" : "bold"}  )
        #plt.savefig ( OUTPUT_FILENAME )

        from sklearn.decomposition import TruncatedSVD, PCA
        fig, ax = plt.subplots(ncols=1, nrows = 1, figsize=(6, 6))
        plt.suptitle ("INITIAL_KMEANS", fontsize =  26, fontweight='semibold' )

        tsvd = TruncatedSVD(n_components=2)
        tsvd.fit( np_vaf )

        np_vaf_transform, mixture_transform = tsvd.transform( np.concatenate(  [np_vaf, mixture_kmeans.T]) )[:-mixture_kmeans.shape[1]], tsvd.transform(np.concatenate( [np_vaf, mixture_kmeans.T]) )[-mixture_kmeans.shape[1]:].T

        plt.axis([np.min(np_vaf_transform[:, 0]) * 2.1,  np.max(np_vaf_transform[:, 0]) * 2.1,  np.min(np_vaf_transform[:, 1]) * 2.1,  np.max(np_vaf_transform[:, 1]) * 2.1])
        for k in range ( np_vaf_transform.shape[0] ):
            ax.scatter ( np_vaf_transform [k, 0] * 2, np_vaf_transform [k, 1] * 2, s = 30, color = "#EAC696" , alpha = 0.6)
        ax.scatter ( 0, 0, s = 30, color = Gr_10[10] , alpha = 0.8)        
        for j in range ( mixture_transform.shape[1] ):    # centroid찍 기
            ax.scatter( mixture_transform[0][j], mixture_transform[1][j], marker='*', color=colorlist[j], edgecolor='black', s=400, label="cluster" + str(j)  )
            ax.text( mixture_transform[0][j], mixture_transform[1][j], "[{}, {}, {}]".format ( round(mixture_kmeans[0][j], 2) , round (mixture_kmeans[1][j], 2) , round(mixture_kmeans[2][j] , 2) ), verticalalignment='top', horizontalalignment='center', fontdict = {"fontsize": 16, "fontweight" : "bold"}   )
        plt.savefig ( OUTPUT_FILENAME )


    return mixture_kmeans, kwargs




def set_initial_parameter(np_vaf, mixture_kmeans, OUTPUT_FILENAME, step, trial, **kwargs):  
    import random
    import numpy as np
    import matplotlib.pyplot as plt
    import seaborn as sns
    from scipy.stats import kde
    
    random.seed ( kwargs["TRIAL"] )
    initial_mixture = np.zeros((kwargs["NUM_BLOCK"], kwargs["NUM_CLONE_ITER"]), dtype="float")
    fail_num = 0

    while True:
        if (kwargs["MAXIMUM_NUM_PARENT"] == 0) & (kwargs["TRIAL"] % 2 == 1):  # Half of the cases (1, 3, 5, ...) :  pick kwargs["NUM_CLONE"] - 1
            if kwargs["NUM_CLONE_ITER"]  >= 2:
                initial_randomsample = sorted ( random.sample ( range( mixture_kmeans.shape[1] ), kwargs["NUM_CLONE_ITER"] - 1 ) )
                initial_mixture = mixture_kmeans [:, initial_randomsample ]
                initial_mixture = np.concatenate((initial_mixture, np.array (1 - np.sum(initial_mixture, axis=1)).reshape(-1, 1) ), axis=1) 
            elif kwargs["NUM_CLONE_ITER"]  == 1:
                initial_randomsample = random.sample ( range( mixture_kmeans.shape[1] ), 1 )
                initial_mixture = mixture_kmeans [:, initial_randomsample ]
        else:  # ,  The other half of the cases : pick kwargs["NUM_CLONE_ITER"] 
            initial_randomsample = sorted ( random.sample ( range( mixture_kmeans.shape[1] ), kwargs["NUM_CLONE_ITER"]  ) )
            initial_mixture = mixture_kmeans [:, initial_randomsample ]

        
        if initial_randomsample not in trial.initial_randomsample_record:   # 그동안 중복이 없어야 탈출가능
            trial.initial_randomsample = trial.initial_randomsample_record  [kwargs["TRIAL"]] = initial_randomsample
            break
        else:
            fail_num += 1
            if fail_num < 15:
                if kwargs ["VERBOSE"] >= 1:
                    fail_num = fail_num  # 그냥 할거 없으니까
                    #print ("{} : 겹침 -> 다시 돌리자".format ( initial_randomsample, np.where (trial.initial_randomsample_record == initial_randomsample) ))
            else:
                break
    
    if fail_num == 15:
        if kwargs ["VERBOSE"] >= 1:
            print ("\t\tCant' make new combination based on this Initial K-means")
            print (trial.initial_randomsample_record )
        step.initial_sampling = False
        return step, kwargs

    


    if np.any(initial_mixture < 0) == True:   # If at least one negative values are, turn it into 0
        if kwargs["NUM_BLOCK"] == 1:
            initial_mixture[:, -1][initial_mixture[:, -1] < 0] = 0.01        # 0? 0.01
        else:
            initial_mixture[:, -1][initial_mixture[:, -1] < 0] = 0        # 0? 0.01
    
    # Use lexsort to sort the array by the third row, then the second row, and finally the first row
    sorted_indices = np.lexsort(tuple(initial_mixture))
    # Use the sorted indices to rearrange the original array
    initial_mixture = initial_mixture[:, sorted_indices]

    #################  fp_index & tn_index  ###############################################
    step.mixture  = np.hstack(( initial_mixture,  np.zeros( (kwargs["NUM_BLOCK"], 1),  dtype = "float")  ))    
    step.fp_index = kwargs["NUM_CLONE_ITER"]
    kwargs["NUM_CLONE"] = kwargs["NUM_CLONE_ITER"] + 1

    #print ("\t Initial mixture : {}".format (step.mixture) )

    for j in range ( initial_mixture.shape[1] ):
        if ( np.all( np.round(initial_mixture[:, j], 2)  != 0) == False ) :  # 하나라도 0이 있다면  (231011_0.003 이런거는 0으로 그냥 보자)
            step.tn_index.append (j)
    
    #print (step.tn_index)
    #########################################################################################

    if kwargs["VERBOSE"] >= 1:
        print ("\t\tinitial_parameter : {}".format( ", ".join(str(np.round(row, 2)) for row in initial_mixture )  ) )
        

    import palettable
    tabl = palettable.tableau.Tableau_20.mpl_colors
    Gr_10 = palettable.scientific.sequential.GrayC_20.mpl_colors
    colorlist = [i for i in tabl]
    
               
    if kwargs["NUM_BLOCK"] == 1:
        plt.figure(figsize=(6, 6))
        plt.suptitle ("Step 0. Initial Selection", fontsize = 20)
        sns.kdeplot ( np_vaf [:, 0] * 2, shade = True)
        kde_np_vaf = kde.gaussian_kde(np_vaf[:, 0] * 2)
        #plt.xlabel("VAF x 2 of the Sample 1",  fontdict={"fontsize": 14})
        plt.ylabel("Density", fontdict={"fontsize": 14})
        plt.xlim(0, 1)
        x_mean_list = []
        for j in range( step.mixture.shape[1] ):
            plt.axvline ( x = step.mixture[0][j], color = colorlist [j  % 20] ) 
            plt.text  ( step.mixture[0][j], kde_np_vaf(step.mixture[0][j]) * 1.08, "{}".format( np.round (step.mixture[0][j] , 2)), verticalalignment='top', ha = "center", fontdict = {"fontsize": 16, "fontweight" : "bold"}  )
            x_mean_list.append ( step.mixture[0][j] )
        plt.title ( "sum = [{}]".format( round( np.sum ( np.array(x_mean_list) ) , 2)   ) ,  ha = 'center', va = 'top', fontdict = {"fontsize" : 12} )
    elif kwargs["NUM_BLOCK"] == 2:
        plt.figure(figsize=(6, 6))
        plt.suptitle ("Step 0. Initial Selection", fontsize = 20)
        plt.scatter ( x = np_vaf [:, 0] * 2, y = np_vaf [:, 1] * 2, s = 30, color =  "#EAC696", alpha = 0.8)
        plt.xlim(-0.02, 1)
        plt.ylim(-0.02, 1)
        x_mean_list, y_mean_list = [], []    
        for j in range( step.mixture.shape[1] ):
            if j == step.mixture.shape[1] - 1:  # FP
                plt.scatter (step.mixture[0][j], step.mixture[1][j], marker = '*', color = Gr_10[10], edgecolor = "black", s = 500, label = "clone {} ({},{})".format (str(j), step.mixture[0][j], step.mixture[1][j] )  ) 
            else:
                plt.scatter (step.mixture[0][j], step.mixture[1][j], marker = '*', color = colorlist[j % 20], edgecolor = "black", s = 500, label = "clone {} ({},{})".format (str(j), step.mixture[0][j], step.mixture[1][j] ) )                 
            #print ( "({},{})".format ( step.mixture[0][j], step.mixture[1][j] ))
            plt.text  ( step.mixture[0][j], step.mixture[1][j] -0.04, "[{},{}]".format( np.round (step.mixture[0][j] , 2), np.round (step.mixture[1][j]  , 2)), verticalalignment = 'top', ha = "center", fontdict = {"fontsize": 16, "fontweight" : "bold"}  )
            x_mean_list.append ( step.mixture[0][j] )
            y_mean_list.append ( step.mixture[1][j] )
        plt.title ( "sum = [{},{}]".format( round( np.sum ( np.array(x_mean_list) ) , 2) , round( np.sum( np.array(y_mean_list) ), 2)  ) ,  ha = 'center', va = 'top', fontdict = {"fontsize" : 12} )
        plt.legend()
    elif kwargs["NUM_BLOCK"] == 3:
        # from mpl_toolkits.mplot3d import Axes3D
        # fig = plt.figure( figsize = (6, 6) )
        # plt.suptitle ("Step 0. Initial Selection", fontsize = 20)
        # ax = fig.add_subplot(111, projection='3d')
        # ax.scatter ( np_vaf [:, 0] * 2, np_vaf [:, 1] * 2, np_vaf [:, 2] * 2 , s = 30, color =  "#EAC696", alpha = 0.8)
        # ax.set_xlim(-0.02, 1);   ax.set_ylim(-0.02, 1);   ax.set_zlim(-0.02, 1)
        # x_mean_list, y_mean_list, z_mean_list = [], [], []
        # for j in range( step.mixture.shape[1] ):
        #     if j == step.mixture.shape[1] - 1:  # FP
        #         ax.scatter (step.mixture[0][j], step.mixture[1][j], step.mixture[2][j], marker = '*', color = Gr_10[10], edgecolor = "black", s = 500, label = "clone {} ({},{},{})".format (str(j), round (step.mixture[0][j], 2) , round( step.mixture[1][j], 2) , round(step.mixture[2][j], 2)  )) 
        #     else:
        #         ax.scatter (step.mixture[0][j], step.mixture[1][j], step.mixture[2][j], marker = '*', color = colorlist[j % 20], edgecolor = "black", s = 500, label = "clone {} ({},{},{})".format (str(j), round(step.mixture[0][j], 2), round(step.mixture[1][j], 2), round (step.mixture[2][j], 2)  )) 
                
        #     #print ( "({},{})".format ( step.mixture[0][j], step.mixture[1][j] ))
        #     ax.text  ( step.mixture[0][j], step.mixture[1][j], step.mixture[2][j], "[{},{},{}]".format( np.round (step.mixture[0][j] , 2), np.round (step.mixture[1][j]  , 2), np.round (step.mixture[2][j]  , 2)), verticalalignment='top', ha = "center", fontdict = {"fontsize": 16, "fontweight" : "bold"}  )
        #     x_mean_list.append ( step.mixture[0][j] )
        #     y_mean_list.append ( step.mixture[1][j] )
        #     z_mean_list.append ( step.mixture[2][j] )
        # plt.title ( "sum = [{},{},{}]".format( round( np.sum ( np.array(x_mean_list) ) , 2) , round( np.sum( np.array(y_mean_list) ), 2), round( np.sum( np.array(z_mean_list) ), 2)  ) ,  ha = 'center', va = 'top', fontdict = {"fontsize" : 12} )
        # plt.legend()

        from sklearn.decomposition import TruncatedSVD, PCA
        fig, ax = plt.subplots(ncols=1, nrows = 1, figsize=(6, 6))
        plt.suptitle ("Step 0. Initial Selection", fontsize = 20)
        tsvd = TruncatedSVD(n_components=2)
        tsvd.fit( np_vaf )

        np_vaf_transform, mixture_transform = tsvd.transform( np.concatenate(  [np_vaf, step.mixture.T]) )[:-step.mixture.shape[1]], tsvd.transform(np.concatenate( [np_vaf, step.mixture.T]) )[-step.mixture.shape[1]:].T
        plt.axis([np.min(np_vaf_transform[:, 0]) * 2.1,  np.max(np_vaf_transform[:, 0]) * 2.1,  np.min(np_vaf_transform[:, 1]) * 2.1,  np.max(np_vaf_transform[:, 1]) * 2.1])
        for k in range ( np_vaf_transform.shape[0] ):
            ax.scatter ( np_vaf_transform [k, 0] * 2, np_vaf_transform [k, 1] * 2, s = 30, color = "#EAC696" , alpha = 0.6)
        
        for j in range ( mixture_transform.shape[1]  - 1 ):    # centroid찍기
            ax.scatter( mixture_transform[0][j], mixture_transform[1][j], marker='*', color=colorlist[j], edgecolor='black', s=400, label="cluster" + str(j)  )
            ax.text( mixture_transform[0][j], mixture_transform[1][j], "[{}, {}, {}]".format ( round(step.mixture[0][j] ,2) , round(step.mixture[1][j], 2), round (step.mixture[2][j], 2) ), verticalalignment='top', horizontalalignment='center', fontdict = {"fontsize": 16, "fontweight" : "bold"}   )
        ax.scatter( 0, 0, marker='s', color = "#413F42", edgecolor='black', s=300, label="FP"  )
        
        plt.legend()

        
    plt.savefig ( OUTPUT_FILENAME )


    return step, kwargs

##################################################################################################################################################################################################################################################


def normalize ( cluster, i ):
    import copy
    import numpy as np

    NUM_BLOCK = cluster.mixture_record[i].shape[0]
    NUM_CLONE = cluster.mixture_record[i].shape[1]

    makeone_index = copy.deepcopy ( cluster.makeone_index_record[i] )
    cluster_normalized_mixture = copy.deepcopy ( cluster.mixture_record[i]  )

    for i in range(NUM_BLOCK):     # Normalization 
        sum = 0
        for j in range(NUM_CLONE):
            if j in makeone_index :   
                sum = sum + cluster_normalized_mixture[i][j]
        cluster_normalized_mixture[i] = np.round( cluster_normalized_mixture[i] / sum, 2) if sum != 0 else 0   # If sum = 0, let mixture = 0
    
    return cluster_normalized_mixture



def movedcolumn ( cluster_hard, cluster_soft, i ):
    import scipy.spatial
    import math, copy
    import numpy as np

    #print ( cluster_hard.mixture_record[i] )

    if cluster_hard.mixture_record[i].shape[1] != cluster_soft.mixture_record[i].shape[1]:
        return []
    
    NUM_CLONE = cluster_hard.mixture_record[i].shape[1]
    NUM_BLOCK = cluster_hard.mixture_record[i].shape[0]

    #   Normalize를 하는게 아무래도 비교가 공평하다
    cluster_hard_normalized_mixture = normalize ( cluster_hard, i   )
    cluster_soft_normalized_mixture = normalize ( cluster_soft, i   )

    #print ("cluster_hard.mixture = {}\ncluster_soft.mixture = {}".format ( cluster_hard.mixture_record[i], cluster_soft.mixture_record[i]  ))
    #print ("cluster_hard_normalized_mixture = {}\ncluster_soft_normalized_mixture = {}".format ( cluster_hard_normalized_mixture, cluster_soft_normalized_mixture  ))

    col_list = []
    for j in range ( NUM_CLONE ):
        # 차원을 고려해줌
        NUM_BLOCK_EFFECTIVE = np.count_nonzero( cluster_hard_normalized_mixture[ : , j]  != 0)
        threshold_dist =  (math.sqrt ( NUM_BLOCK_EFFECTIVE ) * 0.05)
      
        if float(scipy.spatial.distance.euclidean( cluster_hard_normalized_mixture [ : , j ] , cluster_soft_normalized_mixture [ : , j ]  )) > threshold_dist:    # 이것보다 더 많이 움직이면 움직였다고 봄
            if (cluster_hard.fp_index_record [i] == j) :       
                continue
            col_list.append(j)
            print ( "j = {}\tNUM_BLOCK_EFFECTIVE = {}\tthreshold_dist = {}\tdistance = {}".format (j, NUM_BLOCK_EFFECTIVE, round(threshold_dist,3), float(scipy.spatial.distance.euclidean( cluster_hard_normalized_mixture [ : , j ] , cluster_soft_normalized_mixture [ : , j ]  ))))

    return col_list, cluster_hard_normalized_mixture, cluster_soft_normalized_mixture, threshold_dist


def std_movedcolumn ( mixture_matrix , moved_col_list ):
    import numpy as np

    tt =  mixture_matrix [ : , moved_col_list ] 

    #print ( "dimension 고려 없는 std : {}".format (np.std (tt) ))

    weight = np.zeros ( tt.shape[0] )
    std_by_block = np.zeros ( tt.shape[0] )

    for i in range ( tt.shape[0 ]):
        weight[i] = np.count_nonzero  ( tt [i, :] )
        std_by_block[i] = np.std ( tt [i, : ] )

    #print ( "std_by_block = {}\nweight = {}\n".format (std_by_block, weight))   
    #print ( "dimension 고려한 std : {}".format ( round ( np.average (std_by_block, weights = weight), 2)  ))

    return ( round ( np.average (std_by_block, weights = weight), 2)  )


##################################################################################################################################################################################################################################################

def GoStop(step, **kwargs):
    import numpy as np
    previous_No = 5

    if kwargs["STEP"] > previous_No:
        for p in range( kwargs["STEP_TOTAL"] - previous_No,  kwargs["STEP_TOTAL"] ):       # Past 5 records
            if np.sum(np.equal(step.membership, step.membership_record[p])) >= int(kwargs["NUM_MUTATION"] * 0.995):
                if kwargs["VERBOSE"] >= 1:
                    print ( "\t\t\t▶ Stop because the membership is nearly the same with #{0} th".format(p) )
                return "Stop"

            if (np.round(step.mixture, 2) == np.round( step.mixture_record[p] , 2)).all() == True:
                if kwargs["VERBOSE"] >= 1:
                    print ("\t\t\t▶ Stop because the mixture is nearly the same with #{0}".format(p))
                return "Stop"

        # Stopping condition
        i =  step.find_max_likelihood_step (0, kwargs["STEP"]) 
        if (step.likelihood) < 0.99 * step.likelihood_record [ i ]:
            if step.likelihood < -9999990:
                if kwargs["VERBOSE"] >= 1:
                    print ("\t\t\t▶ Stop due to unavailable to make 1")    
                return "Stop"
            else:
                if kwargs["VERBOSE"] >= 1:
                    print ("\t\t\t▶ Stop due to failing to increase the likelihood")
                return "Stop"
    if step.likelihood < -9999990:
        if kwargs["VERBOSE"] >= 1:
            print ("\t\t\t▶ Stop due to unavailable to make 1")    
        return "Stop"

    return "Go"

##################################################################################################################################################################################################################################################


def VAFdensitogram (np_vaf, output_suptitle, output_filename,  **kwargs):
    import matplotlib
    import seaborn as sns
    from scipy.stats import kde
    import numpy as np
    import palettable

    tabl = palettable.tableau.Tableau_20.mpl_colors
    colorlist = [i for i in tabl]

    fig, ax = matplotlib.pyplot.subplots (figsize = (6, 6))
    matplotlib.pyplot.rcParams["font.family"] = 'arial'
    matplotlib.pyplot.suptitle ("{}".format(output_suptitle), fontsize = 20)
    ax.set_xlabel("VAF", fontdict = {"fontsize" : 14})
    ax.set_ylabel("Density", fontdict = {"fontsize" : 14})
    max_y = 0

    for i in range(kwargs["NUM_BLOCK"]):
        x = np.linspace(0,1,101)
        kde_np_vaf_new = kde.gaussian_kde(np_vaf[:,i])
        y = kde_np_vaf_new(x)
        if max_y < np.max(y):
            max_y = np.max(y)
        ax.plot(x, y, color = colorlist[i], label = "sample {0}".format(i))
        
        x_median = np.median (np_vaf[:, i])
        y_median = y [ round (x_median * 100) ]
        if x_median > 0.4:
            ax.text( x_median, max_y * 1.15, "{0}\n→ Monoclonal".format(x_median), verticalalignment='top', ha = "center", fontsize = 15, color = "g")
        elif x_median > 0.25:
            ax.text( x_median, max_y * 1.15, "{0}\n→ Biclonal".format(x_median), verticalalignment='top', ha = "center", fontsize = 15, color = "g")
        else:
            ax.text( x_median, max_y * 1.15, "{}\n→ Polyclonal".format(x_median), verticalalignment='top', ha = "center", fontsize = 15, color = "g")
        ax.vlines( x = x_median, ymin = 0, ymax = y_median , color="g", linestyle = "--", label = "median VAF")

    ax.axis ([0,  1.01,  0,  max_y * 1.3])
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["bottom"].set_linewidth (3)
    ax.legend()



    if output_filename != "NotSave":
        matplotlib.pyplot.savefig(output_filename)
    
    return x_median



def drawfigure (np_vaf, membership, mixture, output_filename, **kwargs):
    import palettable
    import matplotlib
    import numpy as np
    import seaborn as sns

    matplotlib.rcParams["font.family"] =  kwargs["FONT_FAMILY"]

    tabl = palettable.tableau.Tableau_20.mpl_colors
    Gr_10 = palettable.scientific.sequential.GrayC_20.mpl_colors
    colorlist = [i for i in tabl]

    if np_vaf.shape[1] == 3:
        from sklearn.decomposition import TruncatedSVD, PCA
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(ncols=1, nrows = 1, figsize=(6, 6))
        mixture = mixture * 2
        
        tsvd = TruncatedSVD(n_components=2)
        tsvd.fit( np_vaf )

        np_vaf_transform, mixture_transform = tsvd.transform( np.concatenate(  [np_vaf, mixture.T]) )[:-mixture.shape[1]], tsvd.transform(np.concatenate( [np_vaf, mixture.T]) )[-mixture.shape[1]:].T

        plt.axis([np.min(np_vaf_transform[:, 0]) * 1.1,  np.max(np_vaf_transform[:, 0]) * 1.1,  np.min(np_vaf_transform[:, 1]) * 1.1,  np.max(np_vaf_transform[:, 1]) * 1.1])
        for k in range ( np_vaf_transform.shape[0] ):
            ax.scatter ( np_vaf_transform [k, 0] , np_vaf_transform [k, 1] , s = 30, color = "#EAC696" , alpha = 0.6)
        ax.scatter ( 0, 0, s = 30, color = Gr_10[10] , alpha = 0.8)        
        for j in range ( mixture_transform.shape[1] ):    # centroid찍 기
            ax.scatter( mixture_transform[0][j], mixture_transform[1][j], marker='*', color=colorlist[j], edgecolor='black', s=400, label="cluster" + str(j)  )
            ax.text( mixture_transform[0][j], mixture_transform[1][j], "[{}, {}, {}]".format ( round(mixture[0][j], 2) , round (mixture[1][j], 2) , round(mixture[2][j] , 2) ), verticalalignment='top', horizontalalignment='center', fontdict = {"fontsize": 16, "fontweight" : "bold"}   )
        plt.savefig ( output_filename )


        from mpl_toolkits.mplot3d import Axes3D
        fig = plt.figure( figsize = (6, 6) )
        ax = fig.add_subplot(111, projection='3d')
        ax.scatter ( np_vaf [:, 0] , np_vaf [:, 1] , np_vaf [:, 2]  , s = 30, color =  "#EAC696", alpha = 0.8)
        x_mean_list, y_mean_list, z_mean_list = [], [], []
        for j in range( mixture.shape[1] ):
            if j == mixture.shape[1] - 1:  # FP
                ax.scatter (mixture[0][j], mixture[1][j], mixture[2][j], marker = '*', color = Gr_10[10], edgecolor = "black", s = 500, label = "clone {} ({},{},{})".format (str(j), round (mixture[0][j], 2) , round( mixture[1][j], 2) , round(mixture[2][j], 2)  )) 
            else:
                ax.scatter (mixture[0][j], mixture[1][j], mixture[2][j], marker = '*', color = colorlist[j % 20], edgecolor = "black", s = 500, label = "clone {} ({},{},{})".format (str(j), round(mixture[0][j], 2), round(mixture[1][j], 2), round (mixture[2][j], 2)  )) 
                
            #print ( "({},{})".format ( mixture[0][j], mixture[1][j] ))
            ax.text  ( mixture[0][j], mixture[1][j], mixture[2][j], "[{},{},{}]".format( np.round (mixture[0][j] , 2), np.round (mixture[1][j]  , 2), np.round (mixture[2][j]  , 2)), verticalalignment='top', ha = "center", fontdict = {"fontsize": 16, "fontweight" : "bold"}  )
            x_mean_list.append ( mixture[0][j] )
            y_mean_list.append ( mixture[1][j] )
            z_mean_list.append ( mixture[2][j] )
        plt.title ( "sum = [{},{},{}]".format( round( np.sum ( np.array(x_mean_list) ) , 2) , round( np.sum( np.array(y_mean_list) ), 2), round( np.sum( np.array(z_mean_list) ), 2)  ) ,  ha = 'center', va = 'top', fontdict = {"fontsize" : 12} )
        plt.savefig ( output_filename.replace("jpg", "3D.jpg") )
        plt.legend()

    elif np_vaf.shape[1] == 2:
        matplotlib.pyplot.figure(figsize=(6, 6))
        #matplotlib.pyplot.axis([0,  1,  0,  1])
        for k in set (membership):
            matplotlib.pyplot.scatter (x = np_vaf [np.where (membership == k)[0] ,0], y = np_vaf [np.where (membership == k)[0], 1] , color = colorlist[k], 
                                                    label = "cluster{} : {}".format(k, np.unique (membership, return_counts=True)[1][k] )  )
        for j in range ( mixture.shape[1] ) :  # Centroid
            matplotlib.pyplot.scatter (x = mixture [0 , j] * 2, y = mixture [1, j] * 2 , color = colorlist [j] , edgecolor='black', marker = '*', s = 400  )
            matplotlib.pyplot.text( mixture [0, j] * 2 , mixture [1, j] * 2, "[{}, {}]".format( round(mixture [0, j] * 2, 2), round( mixture [1, j] * 2, 2) ), verticalalignment='top', horizontalalignment='right', fontdict = {"fontsize": 14, "fontweight" : "bold"} )    
        matplotlib.pyplot.legend()
        matplotlib.pyplot.savefig(output_filename)
    elif np_vaf.shape[1] == 1:
        matplotlib.pyplot.figure(figsize=(6, 2))
        #matplotlib.pyplot.xlim(0,  1)
        matplotlib.pyplot.yticks ( [] )
        for k in set (membership):
            matplotlib.pyplot.scatter (x = np_vaf [ np.where (membership == k) , 0] , y = [0.5] * len ( np.where ( membership == k ) [0] ),  color = colorlist[k],
                                                    label = "cluster{} : {}".format(k, np.unique (membership, return_counts=True)[1][k] )  )
        #matplotlib.pyplot.scatter (x = np_vaf [: ,0] , y = [0.5] * np_vaf.shape[0],  color=[colorlist[k] for k in membership]  )
        matplotlib.pyplot.savefig(output_filename)



# def f_distribution (NUM_MUTATION_local, NUM_BLOCK_local):
#     import numpy as np
#     import random

#     a = np.random.f (3, 10, (NUM_MUTATION_local * 5, NUM_BLOCK_local))

#     f_max = np.max (a)
#     a_list = []
#     for i in range ( len(a) ):
#         if np.all ( a [ i ,  : ] < ( f_max / 4 ) ) == True:
#             a_list.append(i) 

#     a = a [np.array (random.sample (a_list, NUM_MUTATION_local))]
#     f_max = np.max(a)
#     a = a / f_max

#     return a


def TN_prior_cal(x):
    from scipy.special import expit
    return (1 - expit( 4.2*x - 2.5)) * 2

def multiply_npvaf ( NUM_MUTATION, NUM_BLOCK, np_vaf, nonfp_member_index , b): 
    import numpy as np
    import random

    random.seed (b)

    rrr = np.array ( [80, 100, 100, 120] ) 

    model_np_vaf = np.zeros ( (NUM_MUTATION * len(rrr) , NUM_BLOCK),  dtype = "float")

    for k in range (NUM_MUTATION):
        for i in range (NUM_BLOCK):
            #rrr_weight = rrr * TN_prior_cal ( np_vaf[ nonfp_member_index[k] ][i] * 2 )   # 각 값에 따라 움직이는 범위에 weight를 준다. 낮은 VAF일수록 좀 더 요동치게 해서 background clustering을 방해하려는 목적
            rrr_weight = rrr
            for r_index in range ( len (rrr) ):
                while True:
                    r =  random.randint( int(rrr_weight[r_index] * 0.9), int (rrr_weight[r_index] * 1.1) ) / 100       # 0.72 ~ 0.88, 0.91 ~ 0.110, 0.91 ~ 0.110, 1.08 ~ 1.32 중 하나를 뽑음
                    model_np_vaf [NUM_MUTATION * r_index +  k][i] = np_vaf[ nonfp_member_index[k] ][i] * r    # E step을 다시 한 번 더 돌거면 곱하기 2해서 1 넘어가면 안된다
                    if np_vaf[ nonfp_member_index[k] ][i] * r * 2 <= 3:       
                        break
            # if np_vaf[k][0] < 0.04:
            #     print (model_np_vaf [NUM_MUTATION * 1 +  k], model_np_vaf [NUM_MUTATION * 2 +  k])

    model_sel_np_vaf = model_np_vaf [ random.sample( range(len(model_np_vaf)), NUM_MUTATION ) , : ]         # 쭉 늘어세운 후 NUM_MUTATION 개만큼만 랜덤으로 뽑는다

    return model_sel_np_vaf


def calc_likelihood( np_vaf, mixture, membership, bb, **kwargs ):
    import math, random
    import numpy as np
    import scipy

    NUM_MUTATION, NUM_BLOCK, NUM_CLONE = len(membership), mixture.shape[0], mixture.shape[1]
    #print ("NUM_MUTATION = {}\tmixture = {}".format (NUM_MUTATION, mixture.shape))

    prob = np.zeros(  NUM_MUTATION, dtype="float64")    

    #debug_k = np.where ( ((np_vaf[:, 0] > 0.39) & (np_vaf[:, 0] < 0.41) & (np_vaf [:, 1] > 0.39 ) & (np_vaf [:, 1] < 0.41 )) )[0]
    #print (debug_k)
    
    for k in range (NUM_MUTATION):
        j = membership[k]

        for i in range( NUM_BLOCK ):
            alt_obs, depth_obs = int ( np_vaf [k][i] * 100), 100
            alt_calc, depth_calc = int ( mixture[i][j] * 100), 100
            a, b = alt_calc, depth_obs - alt_calc            

            if alt_obs >= depth_obs:
                continue

            if alt_obs == 0:    # TN or FN?   (합치면 1이 되어야 한다))
                if mixture [i][j] == 0 :  # TN
                    p = math.log10 ( kwargs ["TN_CONFIDENTIALITY"] )

                else: # FN
                    p1 = math.log10 ( 1 - kwargs ["TN_CONFIDENTIALITY"] )
                    p1 = 0

                    try:
                        p2 = math.log10 ( scipy.stats.betabinom.pmf(alt_obs, depth_obs, a+1, b+1) )   # 분자
                    except:
                        p2 = float("-inf")
                
                    try:
                        p = p1 + p2
                    except:
                        p = p1 - 400


            else:   # TP or FP?
                SEQ_ERROR = 0.03

                if mixture [i][j] == 0: # FP
                    try:
                        p = math.log10(scipy.stats.binom.pmf(n = depth_obs, p = SEQ_ERROR, k = alt_obs))
                    except:
                        p = -400
                else:  # TP
                    try:
                        p1 = math.log10(1 - scipy.stats.binom.pmf(n = depth_obs, p = SEQ_ERROR, k = alt_obs))   # Not FP
                    except:
                        p1 = -400
                    p1 = 0
                    try:
                        p2 = math.log10 ( scipy.stats.betabinom.pmf(alt_obs, depth_obs, a+1, b+1) )
                    except:
                        print (alt_obs, depth_obs, a+1, b+1)
                        p2 = float ("-inf")

                    try:
                        p = p1 + p2
                    except:
                        p = p1 - 400

            # if (bb == 0) & (k  in debug_k):
            #     np.set_printoptions(suppress=True) 
            #     print ("k = {}, i = {}, p1 = {}, p2 = {}, p = {}\t\talt_obs = {},depth_obs = {}\talt_calc = {},depth_calc = {}".format (k, i, p1 , round(p2, 3) , round( p, 3), alt_obs, depth_obs, alt_calc, depth_calc  ))
        prob[k] += p



    # Calculate sum for each group
    if kwargs["VERBOSE"] >= 1:
        t = np.unique(membership, return_counts=True)
        if bb == 0:
            for j in range( len(t[0]) ):
                group_sum = np.sum(prob [membership == t[0][j] ])
                print("\t\tSum of cluster{} : {} (n = {})".format( t[0][j], round( group_sum, 2), t[1][j]  ) )

    if kwargs["VERBOSE"] >= 1:
        print ( "\t\tprob = {}\tprob(sum) = {}".format( prob [0:2], np.sum(prob) ))
    return np.sum (prob)
    




def decision_gapstatistics (cluster, np_vaf, **kwargs):
    import pandas as pd
    import numpy as np
    import math, scipy, subprocess
    from sklearn.cluster import KMeans

    NUM_MUTATION, NUM_BLOCK = kwargs["NUM_MUTATION"], kwargs["NUM_BLOCK"]

    Input_B = 5
    Gap_list, Std_list, S_list = np.zeros (kwargs["NUM_CLONE_TRIAL_END"] + 1, dtype ="float"), np.zeros (kwargs["NUM_CLONE_TRIAL_END"] + 1, dtype ="float"), np.zeros (kwargs["NUM_CLONE_TRIAL_END"] + 1, dtype ="float")
    Gap_list = np.array ([ float(12345) for i in Gap_list])    
    maxmaxmax_NUM_CLONE  = 0

    
    for NUM_CLONE_ITER in range (kwargs["NUM_CLONE_TRIAL_START"], kwargs["NUM_CLONE_TRIAL_END"] + 1):
        NUM_CLONE = NUM_CLONE_ITER if cluster.includefp_record [NUM_CLONE_ITER] == False else NUM_CLONE_ITER + 1

        if  (kwargs["VERBOSE"] >= 1):
            print ("NUM_CLONE_ITER == {}".format(NUM_CLONE_ITER))
            print ("\tincludefp = {} (fp_index = {})\t\t-> {}".format( cluster.includefp_record [NUM_CLONE_ITER], cluster.fp_index_record [NUM_CLONE_ITER],  NUM_MUTATION - len(cluster.fp_member_index_record [NUM_CLONE_ITER])   ))

        if (cluster.stepindex_record[NUM_CLONE_ITER] == 0) | (cluster.likelihood_record[NUM_CLONE_ITER] < -9999990): 
            Gap_list [NUM_CLONE_ITER] = float("-inf")

        else:
            if maxmaxmax_NUM_CLONE == 0:
                maxmaxmax_NUM_CLONE = NUM_CLONE_ITER
                
            membership = cluster.membership_record[NUM_CLONE_ITER]
            mixture = cluster.mixture_record[NUM_CLONE_ITER]

            #1. Intra cluster variation (Wk)
            Wk = 0
            for k in range(NUM_MUTATION):
                j = membership [k]
                if (cluster.includefp_record [NUM_CLONE_ITER] == True)  &  (k in cluster.fp_member_index_record [NUM_CLONE_ITER] ):   #  Exclude FP from the calculation
                    continue
                Wk = Wk + math.pow (  scipy.spatial.distance.euclidean(np_vaf[k] * 2, mixture[:, j]),  2)   # Sum of square 
            Wk = round(math.log10(Wk), 3)

            Wk = round (cluster.likelihood_record[NUM_CLONE_ITER])

            Wk = calc_likelihood ( np_vaf, mixture / 2, membership, 0, **kwargs)         # 그냥 여기서 다시 E step calculation을 돌자

            if  (kwargs["VERBOSE"] >= 1):
                print ("\tMy Clustering\tWk  : {}" .format( round(Wk, 1)))


            #2. Random generation & ICC (Wkb)
            Wkb_list = []
            for b in range (Input_B):    # Background Kmeans
                # if kwargs ["NUM_BLOCK"] >= 2:
                #     reference_multiply = multiply_npvaf ( NUM_MUTATION - len(cluster.fp_member_index_record [NUM_CLONE_ITER]) , NUM_BLOCK , np_vaf, 
                #                                                                         sorted (set( range(0, NUM_MUTATION) ) - set(cluster.fp_member_index_record [NUM_CLONE_ITER])) , b )
                # elif kwargs ["NUM_BLOCK"] == 1:
                #     reference_multiply = multiply_npvaf ( NUM_MUTATION , NUM_BLOCK , np_vaf, 
                #                                                                         sorted (set( range(0, NUM_MUTATION) ))  , b ) 

                reference_multiply = multiply_npvaf ( NUM_MUTATION , NUM_BLOCK , np_vaf, sorted (set( range(0, NUM_MUTATION) ))  , b ) 
                reference_np = reference_multiply

                kmeans = KMeans(n_clusters=NUM_CLONE_ITER, 
                                init  = 'k-means++',      # 알아서 초기 centroid 정하는 방법
                                #init = np.delete(cluster.mixture_record [NUM_CLONE_ITER].transpose(), cluster.fp_index_record [NUM_CLONE_ITER], axis = 0),   # 내가 원하는 initial centroid 정하는 방법
                                max_iter = 10, random_state = 0)  
                kmeans.fit( reference_np )  # nparray
                #Wkb_list.append ( round (math.log10(kmeans.inertia_), 3) )        # inertia는 제대로 된 TP, TN, FP, FN 의 확률을 구현할 수 없다
                Wkb_list.append ( calc_likelihood ( reference_np, kmeans.cluster_centers_.T , kmeans.labels_ , b, **kwargs) )        # inertia는 제대로 된 TP, TN, FP, FN 의 확률을 구현할 수 없다

                drawfigure (reference_np * 2, kmeans.labels_,  kmeans.cluster_centers_.T , kwargs["CLEMENT_DIR"] + "/Kmeans/Kmeans.clone"  +   str (NUM_CLONE_ITER) + "." + str(b) + ".jpg", **kwargs)
                subprocess.run (["cp " + kwargs["CLEMENT_DIR"] + "/Kmeans/Kmeans.clone"  +   str (NUM_CLONE_ITER) + "." + str(b) + ".jpg " +  kwargs["COMBINED_OUTPUT_DIR"] + "/Kmeans/Kmeans.clone"  +   str (NUM_CLONE_ITER) + "." + str(b) + ".jpg"], shell = True)

            #Gap_list [NUM_CLONE_ITER] = round ( np.mean(Wkb_list) - Wk, 3)       # 원래대로라면 np.mean(Wkb_list) - Wk 이지만, log10 씌웠다면 Wk - np.mean(Wkb_list)로 하는 것이 맞다
            Gap_list [NUM_CLONE_ITER] = round ( Wk - np.median(Wkb_list) , 3)       # 원래대로라면 np.mean(Wkb_list) - Wk 이지만, log10 씌웠다면 Wk - np.mean(Wkb_list)로 하는 것이 맞다
            Std_list [NUM_CLONE_ITER] = round ( np.std (Wkb_list), 3)
            S_list [NUM_CLONE_ITER] = round (Std_list[NUM_CLONE_ITER] * math.sqrt(1 + Input_B) , 3 )
            

            if (kwargs["VERBOSE"] >= 1):
                #print ("Gap_list : {}\nStd_list : {}\nS_list : {}".format ( Gap_list, Std_list, S_list) )
                #print ("Wkb_list : {}".format ( Wkb_list ) )
                print ("\tRandom noise (B = {}) : \tmedian Wkb = {}\tsdk = {}\tsk (sdk * sqrt ({})) = {}\n\tGap (Wk - Wkb) = {}\n\tPosterior = {}".format (Input_B, round( np.median(Wkb_list), 3) , Std_list[NUM_CLONE_ITER], Input_B + 1, S_list[NUM_CLONE_ITER]  ,  Gap_list[NUM_CLONE_ITER], round (cluster.likelihood_record[NUM_CLONE_ITER]) ))

    if (kwargs["VERBOSE"] >= 1):
        print ("Gap list : {}\nS list : {}\n".format( np.round ( Gap_list [kwargs["NUM_CLONE_TRIAL_START"] : kwargs["NUM_CLONE_TRIAL_END"] + 1]), S_list [kwargs["NUM_CLONE_TRIAL_START"] : kwargs["NUM_CLONE_TRIAL_END"] + 1] ))


    # Max Gap Number
    Gap_list_index = []
    maxmaxmax_NUM_CLONE  = np.argmax ( Gap_list [kwargs["NUM_CLONE_TRIAL_START"] : kwargs["NUM_CLONE_TRIAL_END"] + 1]  ) + kwargs["NUM_CLONE_TRIAL_START"]
    
    Gap_list_df = pd.DataFrame ( Gap_list  ).sort_values ( by = 0, ascending = False)
    
    Gap_list_df = Gap_list_df [ Gap_list_df[0] != float(12345)]
    Gap_list = Gap_list_df[0]   # because of 12345
    
    for i in range( len(Gap_list_df) ):
        if cluster.checkall_strict_record  [ Gap_list_df.index[i] ] == True:
            Gap_list_index.append ( Gap_list_df.index[i]  )
            if kwargs["VERBOSE"] >= 1:
                print ("Gap statistics method (max Gap): {}th ( checkall_strict == True ) optimal NUM_CLONE = {}".format( len(Gap_list_index), Gap_list_df.index[i]  ))
    for i in range( len(Gap_list_df) ):
        if cluster.checkall_strict_record  [ Gap_list_df.index[i] ] == False:
            Gap_list_index.append ( Gap_list_df.index[i]  )
            if kwargs["VERBOSE"] >= 1:
                print ("Gap statistics method (max Gap): {}th ( checkall_strict == False ) optimal NUM_CLONE = {}".format( len(Gap_list_index), Gap_list_df.index[i]  ))


    return Gap_list_index    

# cluster.makeone_compulsory_record

    # Original Gap Number 
    #  
    #     for NUM_CLONE in range (kwargs["NUM_CLONE_TRIAL_START"], kwargs["NUM_CLONE_TRIAL_END"]):
    #         if Gap_list [NUM_CLONE] != float("-inf"):
    #             if Gap_list[NUM_CLONE] >= Gap_list[NUM_CLONE + 1] - S_list [NUM_CLONE +1]:
    #                 maxmaxmax_NUM_CLONE  = NUM_CLONE
    #                 if (kwargs["VERBOSE"] in ["True", "T"]) | (kwargs["VERBOSE"] >= 1):
    #                     print ("Original gap statistics method : optimal NUM_CLONE = {}".format(NUM_CLONE))
    #                 return [maxmaxmax_NUM_CLONE]

    #     print ("Original Gap statistics method : optimal NUM_CLONE = {}".format( maxmaxmax_NUM_CLONE  ))
    #     return [ maxmaxmax_NUM_CLONE  ]     



def decision_silhouette (cluster, np_vaf, **kwargs):
    import pandas as pd
    import numpy as np
    import math, scipy
    from sklearn.metrics import silhouette_samples, silhouette_score
    
    NUM_MUTATION, NUM_BLOCK = kwargs["NUM_MUTATION"], kwargs["NUM_BLOCK"]
    Silhouette_list = np.zeros (kwargs["NUM_CLONE_TRIAL_END"] + 1, dtype ="float")
    
    print ("\n\n------------------- Silhouette Method -------------------\n")
    
    for NUM_CLONE in range (kwargs["NUM_CLONE_TRIAL_START"], kwargs["NUM_CLONE_TRIAL_END"] + 1):
        if (kwargs["VERBOSE"] in ["True", "T"]) | (kwargs["VERBOSE"] >= 0):
            print ("NUM_CLONE == {}".format(NUM_CLONE))
            print ("\tincludefp = {} (fp_index = {})\t\t-> {}".format( cluster.includefp_record [NUM_CLONE], cluster.fp_index_record [NUM_CLONE],  NUM_MUTATION - len(cluster.fp_member_index_record [NUM_CLONE])   ))

        if (cluster.stepindex_record[NUM_CLONE] == 0) | (cluster.likelihood_record[NUM_CLONE] < -9999990):  # 그 어떤 답도 찾을 수 없어서 stepindex = 0인 경우
            Silhouette_list [NUM_CLONE] = 0
            
        else:
            membership = cluster.membership_record[NUM_CLONE]
            membership_fpexclude = np.array ([membership[i]   for i in set (range (0, len(membership)) ) - set( cluster.fp_member_index_record [NUM_CLONE] ) ])
            np_vaf_fpexclude = np.array ([np_vaf[i]   for i in set (range (0, len(np_vaf)) ) - set( cluster.fp_member_index_record [NUM_CLONE] ) ])
                
            silhouette_score_alldata = silhouette_samples(np_vaf_fpexclude , membership_fpexclude)   # (500,1 )  500
            Silhouette_list [NUM_CLONE] = np.mean (silhouette_score_alldata)


    arg_list = [ list(Silhouette_list).index(i) for i in sorted(Silhouette_list, reverse=True)][:2]
    print ("Silhouette_list : {}\narg_list : {}".format  (Silhouette_list, arg_list) )
    
    return arg_list





def XieBeni_calcualtor (cluster, np_vaf, **kwargs):
    import math, scipy
    from sklearn.cluster import KMeans

    if kwargs["NUM_CLONE"] == 1:
        return 1111111111111


    min = float("inf")
    for j1 in range (kwargs["NUM_CLONE"]): 
        for j2 in range (j1 + 1, kwargs["NUM_CLONE"]) :
            if (j1 == cluster.fp_index_record [ kwargs["NUM_CLONE"] ]) | (j2 == cluster.fp_index_record [ kwargs["NUM_CLONE"] ]) :
                continue

            if min > math.pow (scipy.spatial.distance.euclidean(cluster.mixture_record[kwargs["NUM_CLONE"]][:, j1], cluster.mixture_record[kwargs["NUM_CLONE"]][:, j2]), 2):
                min = math.pow (scipy.spatial.distance.euclidean(cluster.mixture_record[kwargs["NUM_CLONE"]][:, j1], cluster.mixture_record[kwargs["NUM_CLONE"]][:, j2]), 2)

    v = 0
    for j in range (kwargs["NUM_CLONE"])  :  
        if (j == cluster.fp_index_record [ kwargs["NUM_CLONE"] ]) :
            continue

        for k in range (kwargs["NUM_MUTATION"]):
            if k in cluster.fp_member_index_record[kwargs["NUM_CLONE"]]:
                continue
            v1 = math.pow( cluster.membership_p_record[kwargs["NUM_CLONE"]][k, j], 2)
            v2 = math.pow ( (scipy.spatial.distance.euclidean(cluster.mixture_record[kwargs["NUM_CLONE"]][:, j], np_vaf[k] * 2)) , 2)
            v = v + (v1 * v2)

    if min == 0:  
        return 1111111111111 - kwargs["NUM_CLONE"]
    else:
        return (v / (min * ( kwargs["NUM_MUTATION"] - len(cluster.fp_member_index_record[kwargs["NUM_CLONE"]]) )))



def decision_XieBeni(cluster, np_vaf, **kwargs):
    import pandas as pd
    import numpy as np
    import math, scipy
    from sklearn.cluster import KMeans


    NUM_MUTATION, NUM_BLOCK = kwargs["NUM_MUTATION"], kwargs["NUM_BLOCK"]
 
    XieBeni_list = np.zeros (kwargs["NUM_CLONE_TRIAL_END"] + 1, dtype ="float")

    for NUM_CLONE in range (kwargs["NUM_CLONE_TRIAL_START"], kwargs["NUM_CLONE_TRIAL_END"] + 1):
        if (kwargs["VERBOSE"] >= 1):
            print ("NUM_CLONE == {}".format(NUM_CLONE))
            print ("\tincludefp = {} (fp_index = {})\t\t-> {}".format( cluster.includefp_record [NUM_CLONE], cluster.fp_index_record [NUM_CLONE],  NUM_MUTATION - len(cluster.fp_member_index_record [NUM_CLONE])   ))

        if (cluster.stepindex_record[NUM_CLONE] == 0) | (cluster.likelihood_record[NUM_CLONE] < -9999990):  
            XieBeni_list [NUM_CLONE] = float("inf")
        else:
            kwargs ["NUM_CLONE"] = NUM_CLONE
            XieBeni_list[NUM_CLONE] = XieBeni_calcualtor (cluster, np_vaf, **kwargs)

        if (kwargs["VERBOSE"] >= 1):
            print ("\tXie-Beni index : {}".format (XieBeni_list[NUM_CLONE]))


    # Max XieBeni Number
    XieBeni_list_df = pd.DataFrame ( XieBeni_list ).sort_values ( by = 0, ascending = True)
    XieBeni_list_df = XieBeni_list_df [ XieBeni_list_df[0] != 0]
    XieBeni_list_index = []

    for i in range( len(XieBeni_list_df) ):
        if kwargs["VERBOSE"] >= 1:
            print ("XieBeni method (min): {}th optimal NUM_CLONE = {}".format(i  + 1, XieBeni_list_df.index[i]  ))
        XieBeni_list_index.append ( XieBeni_list_df.index[i]  )

    return XieBeni_list_index 



def decision_max(cluster, np_vaf, **kwargs):
    import pandas as pd
    import numpy as np
    import math, scipy
    from sklearn.cluster import KMeans


    NUM_MUTATION, NUM_BLOCK = kwargs["NUM_MUTATION"], kwargs["NUM_BLOCK"]

    likelihood_list = np.array ( cluster.likelihood_record )
    num_child_list = np.zeros ( kwargs["NUM_CLONE_TRIAL_END"] + 1, dtype = "int" )
    for j in range (0, kwargs["NUM_CLONE_TRIAL_END"] + 1):
        num_child_list [j] = len (cluster.makeone_index_record [j])

    candidate_df = pd.DataFrame ( np.concatenate ( ([likelihood_list] ,[num_child_list]) , axis = 0) ).transpose().sort_values ( by = 0, ascending = False)
    candidate_df.columns = ["likelihood", "NUM_CHILD"]
    candidate_df = candidate_df.astype ({"NUM_CHILD" : "int"})
    candidate_df ["NUM_CLONE"] = candidate_df.index
    candidate_df = candidate_df.reset_index(drop = True)


    candidate_df_after = pd.DataFrame (columns = candidate_df.columns)
    check = []
    check_index = []

    for k in range ( candidate_df.shape [0] ) :
        if candidate_df.loc[k, "NUM_CHILD"] not in set(check):
            if candidate_df.loc[k, "likelihood"] == float("-inf"):
                candidate_df_after = pd.concat( [ candidate_df_after, candidate_df.loc[ list ( set(list(range( candidate_df.shape[0] ))) - set(check_index) ) ] ])    
                candidate_df_after = candidate_df_after.reset_index(drop = True)
                break

            temp_df = candidate_df [ candidate_df["NUM_CHILD"] ==candidate_df.loc[k, "NUM_CHILD"] ].sort_values ( by = "NUM_CLONE", ascending = True) 
            check_index.append(temp_df.index [0])
            temp_df = temp_df.reset_index(drop = True)
            
            candidate_df_after = pd.concat( [ candidate_df_after, temp_df.loc[ [0] ] ])
            check.append  (candidate_df.loc[k, "NUM_CHILD"])
        
    candidate_df_after = candidate_df_after [ candidate_df_after["NUM_CLONE"] != 0 ]   
    candidate_df_after = candidate_df_after.astype ({"NUM_CHILD" : "int",  "NUM_CLONE" : "int"})
    candidate_df_after = candidate_df_after.reset_index(drop = True)


    for i in range( len(candidate_df_after) ):
        if (kwargs["VERBOSE"] >= 1):
            print ("likelihood = {} : {}th optimal NUM_CLONE = {}\tNUM_CHILD = {}".format( round(candidate_df_after["likelihood"].iloc[i], 2) , i  + 1, candidate_df_after["NUM_CLONE"].iloc [i],  candidate_df_after["NUM_CHILD"].iloc [i]  ))

    return list (candidate_df_after ["NUM_CLONE"] )