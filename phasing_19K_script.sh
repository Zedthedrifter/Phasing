#!/bin/bash

function set_up_intel_env {

CONDA_SUBDIR=osx-64 conda create -y -n eagle_intel python=3.9
#use specific env
source $HOME/miniforge3/bin/activate eagle_intel
echo $CONDA_DEFAULT_ENV #double check the env is correct
#install
conda config --env --set subdir osx-64
#conda install -c bioconda eagle2 -y #for human, not a good option, uninstalled
conda install -c bioconda plink -y
conda install -c bioconda plink2 -y
conda install -c bioconda impute2 -y

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
       --make-bed \
       --out $OUTDIR/${PREFIX2}_Chr${chr}
#convert to gen format for impute2
    plink --bfile $OUTDIR/${PREFIX2}_Chr${chr} \
        --recode oxford \
        --out $OUTDIR/${PREFIX2}_Chr${chr}
done

}

#eagle is for human. this is useless for rice. 
function run_eagle {

INDIR=$1
PREFIX=$2

source $HOME/miniforge3/bin/activate eagle_intel
echo "USING ENV" $CONDA_DEFAULT_ENV

genetic_map=/Users/zchen13/miniforge3/manual/Eagle_v2.4.1/tables/genetic_map_hg38_withX.txt.gz
eagle --bfile=$INDIR/${PREFIX} \
      --geneticMapFile=$genetic_map \
      --outPrefix=genome1_phased \
      --numThreads=4
}

#impute2 from Purdue University and the International Rice Research Institute

#run in chunks of 5 Mb with 250 kb buffer on each side.
function run_impute2 {

INDIR=$1
OUTDIR=$2
PREFIX=$3
REFDIR=$4
CHUNK_SIZE=10000  # 5 Mb chunks + 250 kb buffer (default from IMPUTE2 manual)

source $HOME/miniforge3/bin/activate eagle_intel
echo "USING ENV" $CONDA_DEFAULT_ENV

for chr in {1..12}; do
    start=$(cut -f 3 -d ' ' $INDIR/${PREFIX}_Chr${chr}.gen | head -1)
    end=$(cut -f 3 -d ' ' $INDIR/${PREFIX}_Chr${chr}.gen| tail -1)
    
    # Add a small buffer (optional, e.g., 100kb) to ensure full coverage
    echo $start $end
    start_buf=$((start - 100))
    end_buf=$((end + 100))
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
                -o $OUTDIR/${PREFIX}_Chr${chr}_chunk${chunk_count}    
    done
done

}
#---------------------execute------------------------

function main {

DATA=/Users/zchen13/data/oryza_19K_RGP
WORKDIR=/Users/zchen13/SCRATCH/oryza_19K_RGP
RESULT1=$WORKDIR/01_plink
RESULT2=$WORKDIR/02_impute2
IMPUTE2_REF=/Users/zchen13/miniforge3/manual/impute2/RIS_reference_panel
GENOME=IRGSP-1_0

#set_up_intel_env
function setup_dir {

mkdir -p $WORKDIR
mkdir -p $RESULT1
mkdir -p $RESULT2
#mkdir -p $RESULT3
#mkdir -p $RESULT4
#mkdir -p $RESULT5

}
setup_dir

#STEP 1: CONVERT PLINK FORMAT
#convert_plink $DATA $genome1.biallelic.base.genomewide.SNPs.withID.PlinkFormat $GENOME
#STEP 2: RUN IMPUTE2
run_impute2 $RESULT1 $RESULT2 $GENOME $IMPUTE2_REF
}

main
