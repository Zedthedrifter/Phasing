#!/bin/bash

function set_up_env {

ENV=shapeit5-env

CONDA_SUBDIR=osx-64 conda create -y -n $ENV python=3.9
#use specific env
source $HOME/miniforge3/bin/activate $ENV
echo $CONDA_DEFAULT_ENV #double check the env is correct
#install
conda config --env --set subdir osx-64
#conda install -c bioconda eagle2 -y #for human, not a good option, uninstalled
conda install -c bioconda shapeit5 -y -n $ENV

#MANUAL
}

#convert plink to vcf, per chr
function plink_to_vcf {

INDIR=$1
OUTDIR=$2
PREFIX1=$3
PREFIX2=$4

source $HOME/miniforge3/bin/activate eagle_intel
echo "USING ENV" $CONDA_DEFAULT_ENV

#for i in 1; do
for i in $(seq 1 3); do
    echo "extracting chr $i vcf.bgz file"
    plink2 \
        --pfile $INDIR/$PREFIX1 \
        --chr $i \
        --export vcf bgz id-paste=iid \
        --out $OUTDIR/${PREFIX2}.chr${i}.tmp
    
    # strip INFO right after export (focal is already GT-only in FORMAT)
    #echo "bcftools annotating chr $i vcf.bgz file: add AN AC field"
    #AN AC ARE REQUIRED FOR SHAPEIT5|UNLESS IT'S ALREADY THERE
    #bcftools +fill-tags $OUTDIR/${PREFIX2}.chr${i}.tmp.vcf.gz -Oz -o $OUTDIR/${PREFIX2}.chr${i}.vcf.bgz -- -t AN,AC
    # Reindex the renamed file
    mv $OUTDIR/${PREFIX2}.chr${i}.tmp.vcf.gz $OUTDIR/${PREFIX2}.chr${i}.vcf.bgz
    echo "indexing chr $i vcf.bgz file" #SHAPEIT ALSO REQUIRES INDEX
    bcftools index -f "$OUTDIR/${PREFIX2}.chr${i}.vcf.bgz"
    #RID OF THE FILE WITH INFO
    rm -f $OUTDIR/${PREFIX2}.chr${i}.tmp.vcf.gz $OUTDIR/${PREFIX2}.chr${i}.tmp.log
done
}

#PHASE COMMON ALLELES
function phase_common {

VCF=$1
CHR=$2
RANGE=$3
OUTDIR=$4

#CHOOSE ENV
source $HOME/miniforge3/bin/activate shapeit5-env
echo $CONDA_DEFAULT_ENV #double check the env is correct

# require indexed vcf
SHAPEIT5_phase_common \
    --input $VCF \
    --region ${CHR}:${RANGE} \
    --map ${RIS_reference_panel}/ReferencePanel_4.8M_SNP_chr${CHR}updated.map  \
    --output ${OUTDIR}/${CHR}_${RANGE}.bcf \
    --thread 8

}

function ligate {

INDIR=$1
CHR=$2
OUTDIR=$3

ls -1v ${INDIR}/${CHR}_*.bcf > $OUTDIR/${CHR}_bcf_file_list.txt.tmp


SHAPEIT5_ligate \
    --input $OUTDIR/${CHR}_bcf_file_list.txt.tmp \
    --output ${OUTDIR}/chr_${CHR}_scaffold.bcf \
    --thread 4 \
    --index

}

#---------------------execute------------------------

function main {  

WORKDIR=/Users/zchen13/SCRATCH/Phasing
SAMPLE=/Users/zchen13/data/oryza_19K_RGP
COMMON=/Users/zchen13/SCRATCH/COMMON

#
RESULT1=$WORKDIR/01_SHAPEIT
#RESULT2=$WORKDIR/02_est-sfs
#RESULT3=$WORKDIR/03_ancify #AA INFERENCE USING ANCIFY


#CONSTANTS
VCF_ALL=/Users/zchen13/SCRATCH/ancestral_allele/01_est-sfs/genome1.focal.chr1-12.vcf.gz
RIS_reference_panel=/Users/zchen13/miniforge3/manual/impute2/RIS_reference_panel

#set_up_env

function setup_dir {

mkdir -p $COMMON
mkdir -p $WORKDIR
mkdir -p $RESULT1
#mkdir -p $RESULT2
#mkdir -p $RESULT3
#mkdir -p $RESULT5
#mkdir -p $RESULT5

}
setup_dir

#CONVERT PLINK TO BZG PER CHR: 1-3
#plink_to_vcf $SAMPLE $COMMON genome1.biallelic.base.genomewide.SNPs.withID.PlinkFormat genome1.focal

#PHASE-COMMON--A SMALL REGION
phase_common $COMMON/genome1.focal.chr2.vcf.bgz 1 1-500000 $RESULT1 #CHR 1



}


main