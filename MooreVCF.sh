#!/bin/bash
#$ -cwd
#$ -S /bin/bash

#python3 MooreVCF.py

REF_hg19="/home/goldpm1/reference/hg19/hg19.fa"

/opt/Yonsei/ensembl-vep/101.0/vep \
-v -assembly "GRCh37" --everything --terms "SO" --fork "24" \
-i "/data/project/Alzheimer/CLEMENT/resource/paper/whole_info.vepinput.noheader.txt" \
-o "/data/project/Alzheimer/CLEMENT/resource/paper/whole_info.vepoutput.txt" \
--force --no_stats \
--fasta ${REF_hg19} \
--cache_version 101 \
--dir_cache /data/public/VEP/101 \
--offline \
--cache
#--dir_plugins /data/public/VEP/101/Plugins \
#--plugin dbNSFP,/home/goldpm1/tools/ANNOVAR/humandb/dbnsfp4.2a/dbNSFP4.0a.gz,ALL \
#--plugin SpliceAI,snv=/data/project/DC_WGS/Resource/spliceAI/spliceai_scores.raw.snv.hg38.vcf.gz,indel=/data/project/DC_WGS/Resource/spliceAI/spliceai_scores.raw.indel.hg38.vcf.gz \
#--vcf \