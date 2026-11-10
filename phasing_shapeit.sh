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

function phase_common_chunks {
    VCF=$1
    CHR=$2
    START=$3
    END=$4
    STEP=$5
    OVERLAP=$6
    OUTDIR=$7

    source $HOME/miniforge3/bin/activate shapeit5-env
    echo "USING ENV: $CONDA_DEFAULT_ENV"

    mkdir -p "$OUTDIR"

    for ((chunk_start = START; chunk_start <= END; chunk_start += STEP )); do
        chunk_end=$((chunk_start + STEP + OVERLAP -1))
        if ((chunk_end > END)); then
            chunk_end=$END
        fi
        RANGE="${chunk_start}-${chunk_end}"
        echo "Phasing ${CHR}:${RANGE}"
        SHAPEIT5_phase_common \
            --input "$VCF" \
            --region "${CHR}:${RANGE}" \
            --map "${RIS_reference_panel}/ReferencePanel_4.8M_SNP_chr${CHR}updated.map" \
            --output "${OUTDIR}/${CHR}_${RANGE}.bcf" \
            --thread 8 \
            --filter-maf 0.001 \
            --hmm-ne 7500 
            #it's better to split indica and japonica here
    done
}

#BRINGING THE CHUNKS TOGETHER
function ligate {

INDIR=$1
CHR=$2
OUTDIR=$3

source $HOME/miniforge3/bin/activate shapeit5-env
echo "USING ENV: $CONDA_DEFAULT_ENV"
echo "ligating chunks of chr $CHR"

#the coordinates/files must be in ascending order
ls -1vd ${INDIR}/${CHR}_*-*.bcf | awk -F'/' '{      
    n = split($NF, a, "_");
    print a[2] "\t" $0
}' | sort -k1,1n | cut -f2-  > $OUTDIR/${CHR}_bcf_file_list.txt.tmp

cat $OUTDIR/${CHR}_bcf_file_list.txt.tmp
SHAPEIT5_ligate \
    --input $OUTDIR/${CHR}_bcf_file_list.txt.tmp \
    --output ${OUTDIR}/chr_${CHR}_scaffold.bcf \
    --thread 4 

#bcftools index -f "${OUTDIR}/chr_${CHR}_scaffold.bcf"
#rm -f $OUTDIR/${CHR}_bcf_file_list.txt.tmp
}

function phase_rare {

    VCF=$1
    CHR=$2
    START=$3
    END=$4
    STEP=$5
    OVERLAP=$6
    SCAFFOLD=$7
    OUTDIR=$8

    source $HOME/miniforge3/bin/activate shapeit5-env
    echo "USING ENV: $CONDA_DEFAULT_ENV"
    mkdir -p "$OUTDIR"

    # Verify the scaffold file is indexed
    if [ ! -f "${SCAFFOLD}.csi" ]; then
        echo "Indexing scaffold file..."
        bcftools index -f "$SCAFFOLD"
    fi

    core_start=$START
    core_end=$((core_start + STEP - 1))
        if ((core_end > END)); then core_end=$END; fi

        chunk_end=$((core_end + OVERLAP))
        chunk_start=$((core_start - OVERLAP))
        if ((chunk_end > END)); then
            chunk_end=$END
        fi
        if ((chunk_start < START)); then
            chunk_start=$START
        fi

        # Define ranges and output name using the core (non-overlapping) coordinates
        SCAFFOLD_RANGE="${CHR}:${chunk_start}-${chunk_end}"
        CORE_RANGE="${CHR}:${core_start}-${core_end}"
        OUT_NAME="${OUTDIR}/${CHR}_${core_start}-${core_end}.rare.bcf"

        echo "Phasing rare variants for ${CORE_RANGE} with scaffold range : ${SCAFFOLD_RANGE}"
        SHAPEIT5_phase_rare \
            --input "$VCF" \
            --input-region "$CORE_RANGE" \
            --scaffold "$SCAFFOLD" \
            --scaffold-region "$SCAFFOLD_RANGE" \
            --map "${RIS_reference_panel}/ReferencePanel_4.8M_SNP_chr${CHR}updated.map" \
            --output "$OUT_NAME" \
            --thread 8
 

}

function phase_rare_chunks {
    VCF=$1
    CHR=$2
    START=$3
    END=$4
    STEP=$5
    OVERLAP=$6
    SCAFFOLD=$7
    OUTDIR=$8

    source $HOME/miniforge3/bin/activate shapeit5-env
    echo "USING ENV: $CONDA_DEFAULT_ENV"
    mkdir -p "$OUTDIR"

    # Verify the scaffold file is indexed
    if [ ! -f "${SCAFFOLD}.csi" ]; then
        echo "Indexing scaffold file..."
        bcftools index -f "$SCAFFOLD"
    fi

    for ((core_start = START; core_start <= END; core_start += STEP)); do
        core_end=$((core_start + STEP - 1))
        if ((core_end > END)); then core_end=$END; fi

        chunk_end=$((core_end + OVERLAP))
        chunk_start=$((core_start - OVERLAP))
        if ((chunk_end > END)); then
            chunk_end=$END
        fi
        if ((chunk_start < START)); then
            chunk_start=$START
        fi

        # Define ranges and output name using the core (non-overlapping) coordinates
        SCAFFOLD_RANGE="${CHR}:${chunk_start}-${chunk_end}"
        CORE_RANGE="${CHR}:${core_start}-${core_end}"
        OUT_NAME="${OUTDIR}/${CHR}_${core_start}-${core_end}.rare.bcf"

        echo "Phasing rare variants for ${CORE_RANGE} with scaffold range : ${SCAFFOLD_RANGE}"
        SHAPEIT5_phase_rare \
            --input "$VCF" \
            --input-region "$CORE_RANGE" \
            --scaffold "$SCAFFOLD" \
            --scaffold-region "$SCAFFOLD_RANGE" \
            --map "${RIS_reference_panel}/ReferencePanel_4.8M_SNP_chr${CHR}updated.map" \
            --output "$OUT_NAME" \
            --thread 8
    done
}

#CONCATENATE CHUNKS OF A CHROMOSOME
function concatenate {

INDIR=$1
CHR=$2
OUTDIR=$3

source $HOME/miniforge3/bin/activate est-sfs
echo "USING ENV: $CONDA_DEFAULT_ENV"
echo "CONCATENATE ALL chunks of chr $CHR"

#the coordinates/files must be in ascending order
ls -1vd ${INDIR}/${CHR}_*-*.bcf | awk -F'/' '{      
    n = split($NF, a, "_");
    print a[2] "\t" $0
}' | sort -k1,1n | cut -f2- > $OUTDIR/${CHR}_bcf_file_list.txt.tmp

cat $OUTDIR/${CHR}_bcf_file_list.txt.tmp

bcftools concat \
    --naive \
    -f $OUTDIR/${CHR}_bcf_file_list.txt.tmp \
    -o ${OUTDIR}/chr_${CHR}_scaffold.bcf \
    --threads 8

bcftools index \
    -f ${OUTDIR}/chr_${CHR}_scaffold.bcf
}

#CLEAN ANCESTRAL ALLELE INFO

function clean_AA {

source $HOME/miniforge3/bin/activate est-sfs
echo "USING ENV" $CONDA_DEFAULT_ENV

AA=$1
NAME=$2
OUTDIR=$3

echo "Preparing AA TSV (converting NA to .)..."
    
    # 1. Convert 'NA' to '.' using awk
    # Assuming your TSV has no header and NA is in the 3rd column
    head $AA
    awk -F'\t' 'BEGIN{OFS="\t"} {if ($3 == "NA" || $3 == "") $3="."; print}' \
        "$AA" > "$OUTDIR/${NAME}_cleaned.tsv"
    
    echo "Sorting and compressing cleaned TSV..."
    
    # 2. Sort, bgzip, and index the TSV for bcftools
    # Note: bcftools requires the annotation file to be sorted
    sort -k1,1 -k2,2n "$OUTDIR/${NAME}_cleaned.tsv" | bgzip > "$OUTDIR/${NAME}_cleaned.tsv.gz"
    tabix -s1 -b2 -e2 "$OUTDIR/${NAME}_cleaned.tsv.gz"
}
#ADD ANCESTRAL ALLELE INFO
function add_AA {

source $HOME/miniforge3/bin/activate est-sfs
echo "USING ENV" $CONDA_DEFAULT_ENV

AA=$1
VCF=$2
NAME=$3
OUTDIR=$4
    
    echo "Annotating VCF..."
    echo '##INFO=<ID=AA,Number=1,Type=String,Description="Ancestral Allele">' > $OUTDIR/aa_header.txt
    # 3. Annotate the VCF
    # -c CHROM,POS,AA: matches on the first two columns, adds the AA tag
    bcftools annotate \
        -a "$AA" \
        -h "$OUTDIR/aa_header.txt" \
        -c CHROM,POS,AA \
        -Oz -o "$OUTDIR/${NAME}_wAA.vcf.gz" \
        "$VCF"
    
    echo "Indexing final VCF..."
    bcftools index -f "$OUTDIR/${NAME}_wAA.vcf.gz"
    
    echo "Output: $OUTDIR/${NAME}_wAA.vcf.gz"
}

#CONVERT TO VCZ
function convert_vcz {

INFILE=$1
OUTDIR=$2
NAME=$3

source $HOME/miniforge3/bin/activate ts
echo "USING ENV" $CONDA_DEFAULT_ENV

echo "convert to vcz: ${NAME}"

vcf2zarr convert $INFILE ${OUTDIR}/${NAME}.vcz

}

#---------------------execute------------------------

function main {  

WORKDIR=/Users/zchen13/SCRATCH/Phasing
SAMPLE=/Users/zchen13/data/oryza_19K_RGP
COMMON=/Users/zchen13/SCRATCH/COMMON
AA=/Users/zchen13/SCRATCH/ancestral_allele/02_est-sfs/AA_Obrth-mrdonls-95-85.tsv

#
RESULT1=$WORKDIR/01_S5-COMMON
RESULT2=$WORKDIR/02_S5-LIGATE
RESULT3=$WORKDIR/03_S5-RARE #AA INFERENCE USING ANCIFY
RESULT4=$WORKDIR/04_S5-FINAL
RESULT5=$WORKDIR/05_ADD_AA
RESULT6=$WORKDIR/06_VCZ

#CONSTANTS
RIS_reference_panel=/Users/zchen13/miniforge3/manual/impute2/RIS_reference_panel

#set_up_env

function setup_dir {

mkdir -p $COMMON
mkdir -p $WORKDIR
mkdir -p $RESULT1
mkdir -p $RESULT2
mkdir -p $RESULT3
mkdir -p $RESULT4
mkdir -p $RESULT5
mkdir -p $RESULT6

}
setup_dir

#CONVERT PLINK TO BZG PER CHR: 1-3
#plink_to_vcf $SAMPLE $COMMON genome1.biallelic.base.genomewide.SNPs.withID.PlinkFormat genome1.focal

#PHASE-COMMON--A SMALL REGION
#phase_common $COMMON/genome1.focal.chr1.vcf.bgz 1 1-500000 $RESULT1 #CHR 1
#phase_common_chunks $COMMON/genome1.focal.chr1.vcf.bgz 1 1 200000 50000 20000 $RESULT1
#phase_common_chunks $COMMON/genome1.focal.chr1.vcf.bgz 1 1 43270923 500000 20000 $RESULT1 #total length of chr1: 43270923

#LIGATE CHUNKS
#ligate $RESULT1 1 $RESULT2

#PHASE RARE ALLELES
#phase_rare_chunks $COMMON/genome1.focal.chr1.vcf.bgz 1 1000000 43270923 50000 20000 $RESULT2/chr_1_scaffold.bcf $RESULT3 #use smaller steps! otherwise it will hit error and abort

#PHASE WINDOWS THAT FAILED IN BULK (expand, shrink, or adjust window size)
#phase_rare $COMMON/genome1.focal.chr1.vcf.bgz 1 10500000 10650000 150000 20000 $RESULT2/chr_1_scaffold.bcf $RESULT3 #use smaller steps! otherwise it will hit error and abort
#phase_rare_chunks $COMMON/genome1.focal.chr1.vcf.bgz 1 37200000 37250000 10000 5000 $RESULT2/chr_1_scaffold.bcf $RESULT3 #use smaller steps! otherwise it will hit error and abort
#phase_rare $COMMON/genome1.focal.chr1.vcf.bgz 1 43200000 43270923 70923 20000 $RESULT2/chr_1_scaffold.bcf $RESULT3 #use smaller steps! otherwise it will hit error and abort

#REMOVE OVERLAPPED
#cd $RESULT3
#rm -f 1_43200000-43249999.rare.bcf* 1_10500000-10549999.rare.bcf* 1_10600000-10649999.rare.bcf*
#cd -

#CONCATENATE
#concatenate $RESULT3 1 $RESULT4

#NOW 01 02 03 CAN ALL BE CLEARED
#BETTER TO MAKE A FOLDER FOR EACH CHR DURING THE PROCESS

#ADD ANCESTRAL ALLELES
#clean_AA $AA AA_Obrth-mrdonls-95-85 $RESULT5 #only need to run once for all alleles
#A SMALL CHUNK FOR EXPERIMENTING
add_AA $RESULT5/AA_Obrth-mrdonls-95-85_cleaned.tsv.gz $RESULT3/1_500000-549999.rare.bcf  1_500000-549999 $RESULT5
#add_AA $RESULT5/AA_Obrth-mrdonls-95-85_cleaned.tsv.gz $RESULT4/chr_1_scaffold.bcf  chr_1_scaffold $RESULT5

#CONVERT TO VCZ
convert_vcz $RESULT5/1_500000-549999_wAA.vcf.gz $RESULT6 1_500000-549999
#convert_vcz $RESULT5/chr_1_scaffold_wAA.vcf.gz $RESULT6 chr_1_scaffold
}


main