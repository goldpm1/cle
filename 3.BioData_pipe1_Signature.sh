#!/bin/bash
#$ -cwd
#$ -S /bin/bash


if ! options=$(getopt -o h --long SCRIPT_DIR:,DECISION_MEMBERSHIP_PATH:,NPVAF_PATH:,DONOR:,TISSUE:,OUTPUT_DIR:, -- "$@")
then
    echo "ERROR: invalid options"
    exit 1
fi

eval set -- $options

while true; do
    case "$1" in
        -h|--help)
            echo "Usage"
        shift ;;
        --SCRIPT_DIR)
            SCRIPT_DIR=$2
        shift 2 ;;
        --DECISION_MEMBERSHIP_PATH)
            DECISION_MEMBERSHIP_PATH=$2
        shift 2 ;;
        --NPVAF_PATH)
            NPVAF_PATH=$2
        shift 2 ;;
        --DONOR)
            DONOR=$2
        shift 2 ;;
        --TISSUE)
            TISSUE=$2
        shift 2 ;;
        --OUTPUT_DIR)
            OUTPUT_DIR=$2
        shift 2 ;;
        --)
            shift
            break
    esac
done

rm -rf ${OUTPUT_DIR}
mkdir -p ${OUTPUT_DIR}

echo -e "------------------------------------------------------------------------------------------------ #1. MatrixFormation ------------------------------------------------------------------------------------------------"
echo -e python3 ${SCRIPT_DIR}"/3.BioData_pipe1-1.MatrixFormation.py" \
 --DECISION_MEMBERSHIP_PATH ${DECISION_MEMBERSHIP_PATH} \
 --NPVAF_PATH ${NPVAF_PATH} \
 --DONOR ${DONOR} \
 --TISSUE ${TISSUE} \
 --OUTPUT_DIR ${OUTPUT_DIR}

python3 ${SCRIPT_DIR}"/3.BioData_pipe1-1.MatrixFormation.py" \
 --DECISION_MEMBERSHIP_PATH ${DECISION_MEMBERSHIP_PATH} \
 --NPVAF_PATH ${NPVAF_PATH} \
 --DONOR ${DONOR} \
 --TISSUE ${TISSUE} \
 --OUTPUT_DIR ${OUTPUT_DIR}
sleep 5s



# SigProfiler : MatrixGenerator
echo -e "-\n----------------------------------------------------------------------------------------------- #2. SigProfiler : MatrixGenerator ------------------------------------------------------------------------------------------------"
echo -e  python3 ${SCRIPT_DIR}"/3.BioData_pipe1-2.MatrixGenerator.py" \
 --OUTPUT_DIR ${OUTPUT_DIR}


python3 ${SCRIPT_DIR}"/3.BioData_pipe1-2.MatrixGenerator.py" \
 --OUTPUT_DIR ${OUTPUT_DIR}



# SigProfiler : Assignment
echo -e "\n------------------------------------------------------------------------------------------------ #3. SigProfiler : Assignment-----------------------------------------------------------------------------------------------"
echo -e  python3 ${SCRIPT_DIR}"/3.BioData_pipe1-3.Assignment.py" \
 --OUTPUT_SBS96 ${OUTPUT_DIR}"/output/SBS/BioData.SBS96.all" \
 --ASSIGNMENT_DIR ${OUTPUT_DIR}"/output/Assignment" 


python3 ${SCRIPT_DIR}"/3.BioData_pipe1-3.Assignment.py" \
 --OUTPUT_SBS96 ${OUTPUT_DIR}"/output/SBS/BioData.SBS96.all" \
 --ASSIGNMENT_DIR ${OUTPUT_DIR}"/output/Assignment" 

cp ${OUTPUT_DIR}"/output/Assignment/Assignment_Solution/Activities/Assignment_Solution_Activity_Plots.pdf" ${OUTPUT_DIR%/*}

#/data/project/Alzheimer/CLEMENT/03.combinedoutput/3.BioData/Moore_1D/adrenal_gland_zona_glomerulosa/PD28690-L2/SigProfiler/output/Assignment/Assignment_Solution/Activities/Assignment_Solution_Activity_Plots.pdf


# SigProfiler : Extractor
echo -e "\n------------------------------------------------------------------------------------------------ #4. SigProfiler : Extractor-----------------------------------------------------------------------------------------------"
echo -e  python3 ${SCRIPT_DIR}"/3.BioData_pipe1-4.Extractor.py" \
 --OUTPUT_SBS96 ${OUTPUT_DIR}"/output/SBS/BioData.SBS96.all" \
 --EXTRACTOR_DIR ${OUTPUT_DIR}"/output/DeNovo" 

python3 ${SCRIPT_DIR}"/3.BioData_pipe1-4.Extractor.py" \
 --OUTPUT_SBS96 ${OUTPUT_DIR}"/output/SBS/BioData.SBS96.all" \
 --EXTRACTOR_DIR ${OUTPUT_DIR}"/output/DeNovo" 


# Cosine similiarity
echo -e "\n------------------------------------------------------------------------------------------------ #5. Cosine Similiarity -----------------------------------------------------------------------------------------------"
echo -e python3 ${SCRIPT_DIR}"/3.BioData_pipe1-5.cosinesim.py" \
 --SIGPROFILER_PATH ${OUTPUT_DIR}"/output//DeNovo/SBS96/All_Solutions/SBS96_3_Signatures/Signatures/SBS96_S3_Signatures.txt" \
 --OUTPUT_PATH ${OUTPUT_DIR}"/output/DeNovo/Cosine_sim.txt" 

python3 ${SCRIPT_DIR}"/3.BioData_pipe1-5.cosinesim.py" \
 --SIGPROFILER_PATH ${OUTPUT_DIR}"/output//DeNovo/SBS96/All_Solutions/SBS96_3_Signatures/Signatures/SBS96_S3_Signatures.txt" \
 --OUTPUT_PATH ${OUTPUT_DIR}"/output/DeNovo/Cosine_sim.txt" 
