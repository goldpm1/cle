if __name__ == "__main__":
    from SigProfilerExtractor import sigpro as sig
    import argparse

    parser = argparse.ArgumentParser(description='The below is usage direction.')
    parser.add_argument("--OUTPUT_SBS96", type = str, default = "")
    parser.add_argument("--EXTRACTOR_DIR", type = str, default = "")

    kwargs = {}
    args = parser.parse_args()


    kwargs["OUTPUT_SBS96"] = args.OUTPUT_SBS96
    kwargs["EXTRACTOR_DIR"] = args.EXTRACTOR_DIR


    sig.sigProfilerExtractor( "matrix", 
                                     kwargs["EXTRACTOR_DIR"], 
                                    kwargs["OUTPUT_SBS96"],
                                     reference_genome="GRCh37", opportunity_genome = "GRCh37", context_type = "96", exome = False, 
                                    minimum_signatures=2, maximum_signatures=5, 
                                    nmf_replicates=30, resample = True, batch_size=1, 
                                    cpu= 20, gpu=False, 
                                    nmf_init="random", precision= "single", matrix_normalization= "gmm", seeds= "random", 
                                    min_nmf_iterations= 10000, max_nmf_iterations=1000000, nmf_test_conv= 10000, nmf_tolerance= 1e-15, get_all_signature_matrices = True)