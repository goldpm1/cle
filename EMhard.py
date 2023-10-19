import os
import numpy as np
import  Estep, Mstep, Bunch, miscellaneous
import warnings
warnings.simplefilter (action = 'ignore', category = FutureWarning)
warnings.filterwarnings("ignore")


def whether_trial_acc (max_step, step_index, step, trial, **kwargs):
    if step.likelihood_record [max_step] > trial.likelihood_record [ kwargs["TRIAL"]] : 
        print ("\t\t\t✓ max_step : #{}th step\t\tcheckall_strict = {}\t\tstep.likelihood_record [max_step] = {}".format( max_step , step.checkall_strict_record[max_step], np.round (step.likelihood_record [max_step] , 2) ))    
        trial.acc ( step.mixture_record [max_step],  step.membership_record [max_step], step.likelihood_record [max_step], step.membership_p_record [max_step], step.membership_p_normalize_record [max_step], 
                        step.makeone_index_record[max_step], step.tn_index_record[max_step],   step.fp_index_record[max_step],  step_index, step.fp_member_index_record[max_step], step.includefp_record[max_step], step.fp_involuntary_record[max_step], step.checkall_strict_record[max_step], step.checkall_lenient_record[max_step], max_step, kwargs["TRIAL"] )
    
    return step, trial


def main (input_containpos, df, np_vaf, np_BQ, mixture_kmeans, **kwargs):
    NUM_BLOCK, kwargs["NUM_BLOCK"]= len(df[0]), len(df[0])
    NUM_MUTATION =  kwargs["RANDOM_PICK"]

    cluster = Bunch.Bunch2(**kwargs)

    for NUM_CLONE in range(kwargs["NUM_CLONE_TRIAL_START"], kwargs["NUM_CLONE_TRIAL_END"] + 1):
        kwargs["NUM_CLONE"] = kwargs["NUM_CLONE_ITER"] =  NUM_CLONE
        if kwargs["VERBOSE"] >= 1:
            print("\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n\nNUM_CLONE = {0}".format(NUM_CLONE))
        trial = Bunch.Bunch1(NUM_MUTATION , NUM_BLOCK, NUM_CLONE, kwargs["TRIAL_NO"])
        

        if kwargs["KMEANS_CLUSTERNO"] < kwargs["NUM_CLONE"]:  
            continue

        else:  # Most of the cases
            trial_index, failure_num = 0, 0
            while trial_index < kwargs["TRIAL_NO"]:
                kwargs["TRIAL"] = trial_index

                # if NUM_CLONE ==3 :
                #     if kwargs["TRIAL"] != 10:
                #         trial_index = trial_index + 1
                #         continue

                # if NUM_CLONE == 4:
                #     trial_index = trial_index + 1
                #     continue

                # if NUM_CLONE == 5:
                #     if kwargs["TRIAL"] != 3:
                #         trial_index = trial_index + 1
                #         continue


                step = Bunch.Bunch1(NUM_MUTATION , NUM_BLOCK, NUM_CLONE, kwargs["STEP_NO"])

                step, kwargs = miscellaneous.set_initial_parameter(np_vaf, mixture_kmeans, 
                                                                kwargs["CLEMENT_DIR"] + "/trial/clone" + str(NUM_CLONE) + "." + str(kwargs["TRIAL"]) + "-0.initial_kmeans(hard)." + kwargs["IMAGE_FORMAT"] ,
                                                                step, trial, **kwargs)

                if kwargs["VERBOSE"] >= 1:
                    print("\tTrial #{}\t{}".format(trial_index, trial.initial_randomsample_record  [kwargs["TRIAL"]]  ))

                if step.initial_sampling == False:
                    break


                for step_index in range(0, kwargs["STEP_NO"]):
                    kwargs["STEP"], kwargs["STEP_TOTAL"] = step_index, step_index
                    kwargs["OPTION"] = "hard"
                    if (step_index == (kwargs["STEP_NO"] - 1)):  # 맨 뒤까지 오면 종료
                        trial_index += 1
                        continue
                    if kwargs["VERBOSE"] >= 1:
                        print ("\t\tStep #{}".format(step_index))


                    # if ( kwargs ["NUM_CLONE_ITER"] == 3 )  & (kwargs["TRIAL"] in [0] ) :
                    #     kwargs["DEBUG"]  = True
                    # else:
                    #     kwargs["DEBUG"] = False

                    step = Estep.main(input_containpos, df, np_vaf, np_BQ, step, **kwargs)  
                    
                    
                    ################################ Early terminating condition  (NUM_PARENT, MIN_CLUSTER_SIZE) ################################################
                    if (  kwargs["NUM_CLONE"] -  len (step.makeone_index) - 1  > kwargs["MAXIMUM_NUM_PARENT"] ):     #  1st early terminating condition
                        failure_num = failure_num + 1
                        if kwargs["VERBOSE"] >= 1:
                            print ("\t\t\t♣ STOP:  {}th step,  because in E step →  NUM_CHILD = {}\tNUM_PARENT = {}".format( step_index, len (step.makeone_index), kwargs["NUM_CLONE"] -  len (step.makeone_index) - 1  ))    
                        if kwargs["STEP"] >= 1:
                            print ("kwargs[STEP] = {}\t{}".format(kwargs["STEP"], step.likelihood) )
                            max_step =  step.find_max_likelihood_step(0, kwargs["STEP"] - 1)   # 이상, 이하
                            step, trial = whether_trial_acc (max_step, step_index, step, trial, **kwargs) 
                        break
                
                    #print ( np.unique(step.membership, return_counts=True), np.min( np.unique(step.membership, return_counts=True)[1][np.arange(len(set(step.membership))) != step.fp_index] ), len (step.membership) ) 
                    if  ( np.min( np.unique(step.membership, return_counts=True)[1][np.arange(len(set(step.membership))) != step.fp_index] ) < kwargs["MIN_CLUSTER_SIZE"]  )  :          #  2nd early terminating condition  (except for FP index (the last index))
                        if step.less_than_min_cluster_size == True:     # If it has previous history 
                            failure_num = failure_num + 1
                            extincted_clone_index = np.argmin( np.unique(step.membership, return_counts=True)[1][np.arange(len(set(step.membership))) != step.fp_index] )
                            extincted_clone_count = np.min( np.unique(step.membership, return_counts=True)[1][np.arange(len(set(step.membership))) != step.fp_index] )

                            if kwargs["VERBOSE"] >= 1:
                                print ("\t\t\t♣ STOP: {}th step, because in E step →  The number of variants in clone {}  is {}개 ( < {}). ({})".format(step_index, extincted_clone_index, extincted_clone_count,  kwargs["MIN_CLUSTER_SIZE"], np.unique(step.membership, return_counts=True)[1] ))
                            max_step, trial.checkall_strict_record[ kwargs["TRIAL"] ]  =  step.find_max_likelihood_step(0, kwargs["STEP"] - 2)            # 2번 용서해줬으니 그 전까지 봐야 함
                            step, trial = whether_trial_acc (max_step, step_index, step, trial, **kwargs)
                            break

                        else:    # Just give a pardon in the first
                            step.less_than_min_cluster_size = True
                            

                    if np.all(step.mixture[:, range(step.mixture.shape[1] - 1)] == 0, axis=0).any():      # centroid가 (0) or (0, 0)  or (0, 0, 0) 이 나오는 경우
                        failure_num = failure_num + 1
                        zero_column_index = np.where( np.all( step.mixture[:, range (step.mixture.shape[1] - 1) ] == 0, axis = 0) )[0]
                        if kwargs["VERBOSE"] >= 1:
                            print ("\t\t\t♣ STOP: {}th step, because before M step →  clone {}  is no other than zero point ".format(step_index, zero_column_index  ))
                        max_step, max_step_bool =  step.find_max_likelihood_step(0, kwargs["STEP"] - 1)          
                        step, trial = whether_trial_acc (max_step, step_index, step, trial, **kwargs)
                        break
                    ###########################################################################################
                    
                    
                    step = Mstep.main(input_containpos, df, np_vaf, np_BQ, step, "Hard", **kwargs)   # M step  (Draw figure + Select makeone )

                    if kwargs["STEP"] >= kwargs["COMPULSORY_NORMALIZATION"]:
                        if miscellaneous.checkall (step, "strict", **kwargs)[0] == False:
                            if kwargs["VERBOSE"] >= 2:
                                print ("\t\t➨ STEP {} : sum of mixture is unqualified to checkall_strict\t{}".format( kwargs["STEP"], step.mixture))
                            step.likelihood = -9999999

                    if step.likelihood > -9999990: #  Most of the cases, just accumulate the results
                        step.acc(step.mixture, step.membership, step.likelihood, step.membership_p, step.membership_p_normalize, step.makeone_index, step.tn_index, step.fp_index, step_index, step.fp_member_index, step.includefp, step.fp_involuntary, step.checkall_strict, step.checkall_lenient, kwargs["STEP"], kwargs["STEP"]) 
                        # if step.fp_involuntary == True:
                        #     if kwargs["VERBOSE"] >= 1:
                        #         print ("\t\t\t▶ makeone_index : {}\tparent_index : {}\tfp_index : {}\tfp_involuntary : {}\ttn_index : {}".format(step.makeone_index , sorted( list ( set( list (range(0, kwargs["NUM_CLONE"] )) ) - set( step.makeone_index ) - set ( [step.fp_index] ) )) ,  step.fp_index, step.fp_involuntary, step.tn_index ) )
                        # else:  # Most of the cases
                        if kwargs["VERBOSE"] >= 1:
                            print ("\t\t\t\t▶ makeone_index : {}\tparent_index : {}\tfp_index : {}".format( step.makeone_index , sorted ( list ( set( list (range(0, kwargs["NUM_CLONE"] )) ) - set( step.makeone_index ) - set ( [step.fp_index] ) )),  step.fp_index, step.checkall_lenient ) )

                    

                    if miscellaneous.GoStop(step, **kwargs) == "Stop":
                        max_step, trial.checkall_strict_record[ kwargs["TRIAL"] ] =  step.find_max_likelihood_step (0, kwargs["STEP"])     # Including 0th step,  Whether this trial's best step is prenormalized or not
                        if kwargs["VERBOSE"] >= 1:
                            print ("\t\t✓ max_step : #{}th step\t\tcheckall_strict = {}\t\tstep.likelihood_record [max_step] = {}".format( max_step , trial.checkall_strict_record[ kwargs["TRIAL"] ]  , round (step.likelihood_record [max_step] , 2 )  ))
                        
                        trial.acc ( step.mixture_record [max_step],  step.membership_record [max_step], step.likelihood_record [max_step], step.membership_p_record [max_step], step.membership_p_normalize_record [max_step], step.makeone_index_record[max_step], step.tn_index_record[max_step],  step.fp_index_record[max_step],  step_index + 1, step.fp_member_index_record[max_step], step.includefp_record[max_step], step.fp_involuntary_record[max_step], step.checkall_strict_record[max_step], step.checkall_lenient_record[max_step], max_step, kwargs["TRIAL"] )
                        trial_index = trial_index + 1
                        failure_num = 0
                        break

                if failure_num >= 1: 
                    if kwargs["VERBOSE"] >= 3:
                        print ("\t\t\t\tfailure_num = 1  → Give up and pass to the next trial")
                    
                    failure_num = 0
                    trial_index = trial_index + 1
                

            i, cluster.checkall_strict_record [ kwargs["NUM_CLONE_ITER"] ] =  trial.find_max_likelihood_trial ( 0, kwargs["TRIAL_NO"])             # Best trial을 찾되, 그것이 prenormalization = True인지 False인지 저장한다
        
            if kwargs["VERBOSE"] >= 1:
                print ("\n\n\tIn NUM_CLONE = {}, we chose {}th trial, {}th step\n\t(trial.likelihood_record : {}, checkall_strict : {})\n\tFP_index : {}\tlen(fp_member_index) : {}".format(kwargs["NUM_CLONE_ITER"], i, trial.max_step_index_record[i], np.round ( trial.likelihood_record ), cluster.checkall_strict_record[ kwargs["NUM_CLONE_ITER"] ] ,trial.fp_index_record[i],  len (trial.fp_member_index_record[i] ) ) )

            if trial.max_step_index_record [i]  != -1:   # If available in this trial
                os.system ("cp " + kwargs["CLEMENT_DIR"] + "/trial/clone" + str(kwargs["NUM_CLONE_ITER"]) + "." + str( i ) + "-"  + str(  trial.max_step_index_record [i]  ) + "\(hard\)." + kwargs["IMAGE_FORMAT"] + " " + 
                    kwargs["CLEMENT_DIR"] + "/candidate/clone" + str(kwargs["NUM_CLONE_ITER"]) + ".\(hard\)." + kwargs["IMAGE_FORMAT"]  ) 
            
            cluster.acc ( trial.mixture_record [i], trial.membership_record [i], trial.likelihood_record [i], trial.membership_p_record [i], trial.membership_p_normalize_record [i], trial.stepindex_record [i], i, trial.max_step_index_record [i], trial.makeone_index_record[i], trial.tn_index_record[i],  trial.fp_index_record[i], trial.includefp_record[i], trial.fp_involuntary_record[i], trial.fp_member_index_record [i], **kwargs )  

    return cluster

    #print ("cluster_hard.makeone_index_record : {}".format(cluster_hard.makeone_index_record))
