if __name__ == "__main__":
    from SigProfilerAssignment import Analyzer as Analyze
    import argparse

    parser = argparse.ArgumentParser(description='The below is usage direction.')
    parser.add_argument("--OUTPUT_SBS96", type = str, default = "")
    parser.add_argument("--ASSIGNMENT_DIR", type = str, default = "")

    kwargs = {}
    args = parser.parse_args()

    kwargs["OUTPUT_SBS96"] = args.OUTPUT_SBS96
    kwargs["ASSIGNMENT_DIR"] = args.ASSIGNMENT_DIR

    #help ( Analyze )
    Analyze.cosmic_fit ( kwargs["OUTPUT_SBS96"], 
                                    kwargs["ASSIGNMENT_DIR"], 
                                    input_type="matrix", 
                                    cosmic_version=3.3, 
                                    exome=False,
                                    genome_build="GRCh37", 
                                    signature_database=None,
                                    exclude_signature_subgroups = [ 'Colibactin_signatures', 'Unknown_signatures', 'APOBEC_signatures', 'HR_deficiency_signatures' , 'MMR_deficiency_signatures',   'BER_deficiency_signatures'],   #  'remove' 가 붙은건 cnvpytor, 없는건 master
                                    export_probabilities=False,
                                    export_probabilities_per_mutation=False, 
                                    make_plots=True,
                                    verbose=False)
    

    # exclude_signature_subgroups = ['MMR_deficiency_signatures',
    #                            'POL_deficiency_signatures',
    #                            'HR_deficiency_signatures' ,
    #                            'BER_deficiency_signatures',
    #                            'Chemotherapy_signatures',
    #                            'Immunosuppressants_signatures'
    #                            'Treatment_signatures'
    #                            'APOBEC_signatures',
    #                            'Tobacco_signatures',
    #                            'UV_signatures',
    #                            'AA_signatures',
    #                            'Colibactin_signatures',
    #                            'Artifact_signatures',
    #                            'Lymphoid_signatures']
    

#/opt/Yonsei/python/3.8.1/lib/python3.8/site-packages/SigProfilerAssignment/decomposition.py
#/home/goldpm1/miniconda3/envs/cnvpytor/lib/python3.7/site-packages/SigProfilerAssignment/decomposition.py

#'Unknown_signatures' :{'SBS':['8','12','16', '17a','17b','19','23','28','33','34','37','39','40','41', '89','91','93','94'], 'DBS':[], 'ID':[]} 