#!/bin/bash

function set_up_intel_env {

ENV=eagle_intel

CONDA_SUBDIR=osx-64 conda create -y -n $ENV python=3.9
#use specific env
source $HOME/miniforge3/bin/activate $ENV
echo $CONDA_DEFAULT_ENV #double check the env is correct
#install
conda config --env --set subdir osx-64
#conda install -c bioconda eagle2 -y #for human, not a good option, uninstalled
conda install -c bioconda plink -y -n $ENV
conda install -c bioconda plink2 -y -n $ENV
conda install -c bioconda impute2 -y -n $ENV

#MANUAL
}

function convert_plink {

INDIR=$1
OUTDIR=$2
PREFIX1=$3
PREFIX2=$4

source $HOME/miniforge3/bin/activate eagle_intel
echo "USING ENV" $CONDA_DEFAULT_ENV

for chr in {1..12}; do
    plink2 --pfile $INDIR/$PREFIX1 \
           --chr Chr${chr} \
           --make-pgen \
           --out $OUTDIR/${PREFIX2}_Chr${chr}
    
    plink2 --pfile $OUTDIR/${PREFIX2}_Chr${chr} \
        --max-alleles 2 \
       --make-bed \
       --out $OUTDIR/${PREFIX2}_Chr${chr}
#convert to gen format for impute2
    plink --bfile $OUTDIR/${PREFIX2}_Chr${chr} \
        --recode oxford \
        --out $OUTDIR/${PREFIX2}_Chr${chr}
done

}

#impute2 from Purdue University and the International Rice Research Institute

#run in chunks of 5 Mb with 250 kb buffer on each side.
function run_impute2 {

INDIR=$1
OUTDIR=$2
PREFIX=$3
REFDIR=$4
CHUNK_SIZE=1000000  # 1 Mb chunks + 250 kb buffer (default from IMPUTE2 manual)
BUFFER=250000 #optional, e.g., 100kb

source $HOME/miniforge3/bin/activate eagle_intel
echo "USING ENV" $CONDA_DEFAULT_ENV

for chr in {1..2}; do
    start=$(cut -f 3 -d ' ' $INDIR/${PREFIX}_Chr${chr}.gen | head -1)
    end=$(cut -f 3 -d ' ' $INDIR/${PREFIX}_Chr${chr}.gen| tail -1)
    
    # Add a small buffer (optional, e.g., 100kb) to ensure full coverage
    echo $start $end
    start_buf=$((start - $BUFFER))
    end_buf=$((end + $BUFFER))
    #
    echo "Chromosome $chr range: $start to $end"
    
    # Calculate chunk boundaries starting from the chromosome start
    chunk_start=$((start - (start % CHUNK_SIZE)))
    # If chunk_start < start, adjust to ensure we cover from the beginning
    if [ $chunk_start -lt $start ]; then
        chunk_start=$start
    fi
    
    # Process each 5 Mb chunk
    chunk_count=0
    while [ $chunk_start -le $end ]; do
        chunk_count=$((chunk_count + 1))
        chunk_end=$((chunk_start + CHUNK_SIZE - 1))
        
        # Ensure we don't exceed chromosome end
        if [ $chunk_end -gt $end ]; then
            chunk_end=$end
        fi
        
        echo "  Chunk $chunk_count: $chunk_start - $chunk_end (int: $int_start - $int_end)"
        impute2 -m $REFDIR/ReferencePanel_4.8M_SNP_chr${chr}updated.map \
                -h $REFDIR/GORP_phased_chr${chr}.phased.hap \
                -l $REFDIR/GORP_phased_chr${chr}.phased.newID.legend \
                -g $INDIR/${PREFIX}_Chr${chr}.gen \
                -int $chunk_start $chunk_end \
                -Ne 10000 \
                -k_hap 100 \
                --threads 4 \
                -o $OUTDIR/${PREFIX}_Chr${chr}_chunk${chunk_count}    
    done
done

}
#---------------------execute------------------------

function main {

DATA=/Users/zchen13/data/oryza_3K
WORKDIR=/Users/zchen13/SCRATCH/oryza_3K_RGP
RESULT1=$WORKDIR/01_plink
RESULT2=$WORKDIR/02_impute2
IMPUTE2_REF=/Users/zchen13/miniforge3/manual/impute2/RIS_reference_panel
GENOME=IRGSP-1_0

#set_up_intel_env
function setup_dir {

mkdir -p $WORKDIR
mkdir -p $RESULT1
mkdir -p $RESULT2


}
setup_dir

#CONVERT PLINK FORMAT
#convert_plink $DATA $RESULT1 genome1.genomewide.SNPs.withID.PlinkFormat $GENOME
#RUN IMPUTE2
run_impute2 $RESULT1 $RESULT2 $GENOME $IMPUTE2_REF
}

main
