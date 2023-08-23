if __name__ == "__main__":
    from SigProfilerMatrixGenerator.scripts import SigProfilerMatrixGeneratorFunc as matGen
    import argparse

    parser = argparse.ArgumentParser(description='The below is usage direction.')
    parser.add_argument("--OUTPUT_DIR", type = str, default = "")

    kwargs = {}
    args = parser.parse_args()

    kwargs["OUTPUT_DIR"] = args.OUTPUT_DIR

    matrices = matGen.SigProfilerMatrixGeneratorFunc ( "BioData",  # project
                                                                                        "GRCh37",   # reference_genome
                                                                                        kwargs["OUTPUT_DIR"],  # path_to_input_files 
                                                                                        exome=False, 
                                                                                        bed_file=None, 
                                                                                        chrom_based=False, 
                                                                                        plot=True, 
                                                                                        tsb_stat=False, 
                                                                                        seqInfo=True)