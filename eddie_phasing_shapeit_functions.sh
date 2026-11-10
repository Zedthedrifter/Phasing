#!/bin/bash

# Initialise the environment modules
. /etc/profile.d/modules.sh
module load miniforge
module load roslin/python/3.8.10 
#module load igmm/apps/bcftools/1.20 

#convert plink to vcf, per chr
function plink_to_vcf {

CHR=$1
INDIR=$2
OUTDIR=$3
PREFIX1=$4
PREFIX2=$5

i=$CHR

conda activate shapeit5-env
echo "USING ENV" $CONDA_DEFAULT_ENV

    echo "extracting chr $i vcf.bgz file"
    plink2 \
        --pfile $INDIR/$PREFIX1 \
        --chr $i \
        --export vcf bgz id-paste=iid \
        --out $OUTDIR/${PREFIX2}.chr${i}
    #the output is already vcf.bgz format
    echo "indexing chr $i vcf.bgz file" #SHAPEIT ALSO REQUIRES INDEX
    bcftools index -f "$OUTDIR/${PREFIX2}.chr${i}.vcf.bgz"
    #RID OF THE FILE WITH INFO
    rm -f $OUTDIR/${PREFIX2}.chr${i}.tmp.log
}

#PHASE COMMON 
function phase_common {

    chunk_start=$1
    VCF=$2
    CHR=$3
    END=$4
    STEP=$5
    OVERLAP=$6
    OUTDIR=$7

#CHOOSE ENV
conda activate shapeit5-env
echo "USING ENV" $CONDA_DEFAULT_ENV3

RIS_reference_panel=/exports/eddie/scratch/zchen13/DATA/RIS_reference_panel

#MAKE CHUNK
chunk_end=$((chunk_start + STEP + OVERLAP -1))
        if ((chunk_end > END)); then
            chunk_end=$END
        fi
        RANGE="${chunk_start}-${chunk_end}"
        echo "Phasing ${CHR}:${RANGE}"

# require indexed vcf
SHAPEIT5_phase_common \
    --input $VCF \
    --region ${CHR}:${RANGE} \
    --map ${RIS_reference_panel}/ReferencePanel_4.8M_SNP_chr${CHR}updated.map  \
    --output ${OUTDIR}/${CHR}_${RANGE}.bcf \
    --thread 1 \
    --filter-maf 0.001 \
    --hmm-ne 7500 

}

#BRINGING THE CHUNKS TOGETHER
function ligate {

CHR=$1
INDIR=$2
OUTDIR=$3

#CHOOSE ENV
conda activate shapeit5-env
echo "USING ENV" $CONDA_DEFAULT_ENV3

echo "ligating chunks of chr $CHR"

#the coordinates/files must be in ascending order
ls -1vd ${INDIR}/CHR_${CHR}/${CHR}_*-*.bcf | awk -F'/' '{      
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

    core_start=$1
    VCF=$2
    SCAFFOLD=$3
    CHR=$4
    END=$5
    STEP=$6
    OVERLAP=$7
    OUTDIR=$8

#CHOOSE ENV
conda activate shapeit5-env
echo "USING ENV" $CONDA_DEFAULT_ENV3

RIS_reference_panel=/exports/eddie/scratch/zchen13/DATA/RIS_reference_panel

# Verify the scaffold file is indexed
    if [ ! -f "${SCAFFOLD}.csi" ]; then
        echo "Indexing scaffold file..."
        bcftools index -f "$SCAFFOLD"
    fi

    #MAKE CHUNK
    core_end=$((core_start + STEP - 1))
    if ((core_end > END)); then core_end=$END; fi

    chunk_start=$((core_start - OVERLAP))
    chunk_end=$((core_end + OVERLAP))
    
    if ((chunk_start < START)); then
        chunk_start=$START
    fi
    if ((chunk_end > END)); then
        chunk_end=$END
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

#CONCATENATE CHUNKS OF A CHROMOSOME
function concatenate {

CHR=$1
INDIR=$2
OUTDIR=$3

#CHOOSE ENV
conda activate shapeit5-env
echo "USING ENV" $CONDA_DEFAULT_ENV3
echo "CONCATENATE ALL chunks of chr $CHR"

#the coordinates/files must be in ascending order
ls -1vd ${INDIR}/CHR_${CHR}/${CHR}_*-*.bcf | awk -F'/' '{      
    n = split($NF, a, "_");
    print a[2] "\t" $0
}' | sort -k1,1n | cut -f2- > $OUTDIR/${CHR}_bcf_file_list.txt.tmp

cat $OUTDIR/${CHR}_bcf_file_list.txt.tmp

bcftools concat \
    --naive \
    -f $OUTDIR/${CHR}_bcf_file_list.txt.tmp \
    -o ${OUTDIR}/chr_${CHR}.conca.bcf \
    --threads 8
echo "index the final bcf"
bcftools index \
    -f ${OUTDIR}/chr_${CHR}.conca.bcf
}

function add_AA {

conda activate shapeit5-env
echo "USING ENV" $CONDA_DEFAULT_ENV

CHR=$1
AA=$2
SURFIX=$3
INDIR=$4
OUTDIR=$5
    
    #secho '##INFO=<ID=AA,Number=1,Type=String,Description="Ancestral Allele">' > $OUTDIR/aa_header.txt
    # 3. Annotate the VCF
    # -c CHROM,POS,AA: matches on the first two columns, adds the AA tag
    VCF=$INDIR/chr_${CHR}.${SURFIX}
    echo "Annotating $VCF..."
    bcftools annotate \
        -a "$AA" \
        -h "$OUTDIR/aa_header.txt" \
        -c CHROM,POS,AA \
        -Oz -o "$OUTDIR/chr_${CHR}.${SURFIX}" \
        "$VCF"
    
    echo "Indexing final VCF..."
    bcftools index -f "$OUTDIR/chr_${CHR}.${SURFIX}"
    
    echo "Output: $OUTDIR/chr_${CHR}.${SURFIX}"
}

#CONVERT TO VCZ
function convert_vcz {

CHR=$1
SURFIX=$2
INDIR=$3
OUTDIR=$4


conda activate ts-env
echo "USING ENV" $CONDA_DEFAULT_ENV

INFILE=$INDIR/chr_${CHR}.${SURFIX}
echo "convert to vcz: ${INFILE}"

vcf2zarr convert $INFILE ${OUTDIR}/chr_${CHR}.vcz

}

function python_tsinfer {

CHR=$1
INDIR=$2
OUTDIR=$3

conda activate ts-env
echo "USING ENV" $CONDA_DEFAULT_ENV

echo "INFERRING ARG FOR: chr_${CHR}"
#using python3 of the env called
python3 $HOME/Phasing/eddie_ts_tree.py run_tsinfer --chr 1 \
                                                   --vcz $INDIR/chr_${CHR}.vcz \
                                                   --outfile $OUTDIR/chr_${CHR}.tree #\
#                                                   --start 100 \
#                                                   --end 10100   
}
#================================================================================
#USER INPUTS
MY_FUNCTION=$1 #YOU don't need this line. just a reminder of what the first variable means
echo "Processing GE-ID $SGE_TASK_ID"
# Check command-line argument and call the function
case "$1" in
    plink_to_vcf)
        plink_to_vcf $SGE_TASK_ID $2 $3 $4 $5 # 
        ;;
    phase_common)
        phase_common $SGE_TASK_ID $2 $3 $4 $5 $6 $7 # 
        ;;
    ligate)
        ligate $SGE_TASK_ID $2 $3 # 
        ;;
    phase_rare)
        phase_rare $SGE_TASK_ID $2 $3 $4 $5 $6 $7 $8 # 
        ;;
    concatenate)
        concatenate $SGE_TASK_ID $2 $3 # 
        ;;
    add_AA)
        add_AA $SGE_TASK_ID $2 $3 $4 $5 # 
        ;;
    convert_vcz)
        convert_vcz $SGE_TASK_ID $2 $3 $4 # 
        ;;
    python_tsinfer)
        python_tsinfer $SGE_TASK_ID $2 $3 # 
        ;;
    *)
        echo "Unknown function: $1: Usage: $0 {plink_to_vcf} [args]"
        exit 1
        ;;
esac