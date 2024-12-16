#!/bin/bash
#SBATCH --chdir=./
#SBATCH --job-name=index
#SBATCH --partition nocona
#SBATCH --nodes=1 --ntasks=20
#SBATCH --time=48:00:00
#SBATCH --mem-per-cpu=4G

cd 

cd references

source activate bcftools


# index references


# Swamp Sparrow
# GCA_028018845.1_bMelGeo1.pri_genomic.fna
bwa-mem2 index GCA_028018845.1_bMelGeo1.pri_genomic.fna

samtools faidx GCA_028018845.1_bMelGeo1.pri_genomic.fna

java -jar picard.jar CreateSequenceDictionary \
R=/home/jmanthey/references/GCA_028018845.1_bMelGeo1.pri_genomic.fna \
O=/home/jmanthey/references/GCA_028018845.1_bMelGeo1.pri_genomic.dict


# Myrtle Warbler
# GCA_001746935.2_mywa_2.1_genomic.fna
bwa-mem2 index GCA_001746935.2_mywa_2.1_genomic.fna

samtools faidx GCA_001746935.2_mywa_2.1_genomic.fna

java -jar picard.jar CreateSequenceDictionary \
R=/home/jmanthey/references/GCA_001746935.2_mywa_2.1_genomic.fna \
O=/home/jmanthey/references/GCA_001746935.2_mywa_2.1_genomic.dict


# Galapagos flycatcher
# GCA_038380795.1_CAS-CDF_Pyrnan_1.0_genomic.fna
bwa-mem2 index GCA_038380795.1_CAS-CDF_Pyrnan_1.0_genomic.fna

samtools faidx GCA_038380795.1_CAS-CDF_Pyrnan_1.0_genomic.fna

java -jar picard.jar CreateSequenceDictionary \
R=/home/jmanthey/references/GCA_038380795.1_CAS-CDF_Pyrnan_1.0_genomic.fna \
O=/home/jmanthey/references/GCA_038380795.1_CAS-CDF_Pyrnan_1.0_genomic.dict


# Swainson's Thrush
# GCF_009819885.2_bCatUst1.pri.v2_genomic.fna
bwa-mem2 index GCF_009819885.2_bCatUst1.pri.v2_genomic.fna

samtools faidx GCF_009819885.2_bCatUst1.pri.v2_genomic.fna

java -jar picard.jar CreateSequenceDictionary \
R=/home/jmanthey/references/GCF_009819885.2_bCatUst1.pri.v2_genomic.fna \
O=/home/jmanthey/references/GCF_009819885.2_bCatUst1.pri.v2_genomic.dict









