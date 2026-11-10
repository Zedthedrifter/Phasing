#!/bin/sh
# Grid Engine options (lines prefixed with #$)



WORKDIR=/exports/eddie/scratch/zchen13/Phasing
SAMPLE=/exports/eddie/scratch/zchen13/DATA/oryza_19K_RGP/


#
RESULT0=$WORKDIR/00_PLINK_CONVERT
RESULT1=$WORKDIR/01_S5-COMMON
RESULT2=$WORKDIR/02_S5-LIGATE
RESULT3=$WORKDIR/03_S5-RARE #AA INFERENCE USING ANCIFY
RESULT4=$WORKDIR/04_S5-FINAL
RESULT5=$WORKDIR/05_ADD_AA
RESULT6=$WORKDIR/06_VCZ
RESULT7=$WORKDIR/07_TS

#CONSTANTS
RIS_reference_panel=/exports/eddie/scratch/zchen13/DATA/RIS_reference_panel
AA=/exports/eddie/scratch/zchen13/Phasing/05_ADD_AA/AA_Obrth-mrdonls-95-85_cleaned.tsv.gz #copied from local

function setup_dir {

mkdir -p $WORKDIR
mkdir -p $RESULT0
mkdir -p $RESULT1
mkdir -p $RESULT2
mkdir -p $RESULT3
mkdir -p $RESULT4
mkdir -p $RESULT5
mkdir -p $RESULT6
mkdir -p $RESULT7
mkdir -p $WORKDIR/run_common
mkdir -p $WORKDIR/run_rare #if you are running a lot of jobs, put the runfiles in scratch so they do not fill up home

}
setup_dir

# Initialise the environment modules
#. /etc/profile.d/modules.sh
#module load miniforge

#CONVERT PLINK TO VCZ
function convert_plink {
    # Grid Engine uses qsub. Resource requests are passed as -l (lowercase L) options.
    # Memory (h_vmem) is per-core, so for 3G total with 8 cores, we request 3G per core (or adjust as needed).
    # -pe sharedmem 8 requests the 8 CPUs.
    JOB=$(qsub -t 2-12 \
               -pe sharedmem 1 \
               -l h_vmem=4G \
               -N "convert_plink_vcf" \
               -M zchen13@ed.ac.uk     \
               -m e \
               $HOME/Phasing/eddie_phasing_shapeit_functions.sh plink_to_vcf $SAMPLE \
                                                                             $RESULT0 \
                                                                             genome1.biallelic.base.genomewide.SNPs.withID.PlinkFormat genome1.focal
                                                                             )

    # qsub returns the job ID directly (often with a cluster suffix like 12345.1)
    echo "$JOB"
}
#PHASE COMMON 
function phase_common {

    CHR=$1
    START=$2
    END=$3
    STEP=$4
    OVERLAP=$5
    MEM=$6
    JOB=$(qsub -t ${START}-$END:$STEP \
               -pe sharedmem 1 \
               -l h_vmem=${MEM} \
               -N "com_${CHR}" \
               $HOME/Phasing/eddie_phasing_shapeit_functions.sh phase_common $RESULT0/genome1.focal.chr${CHR}.vcf.bgz \
                                                                             $CHR \
                                                                             $END \
                                                                             $STEP \
                                                                             $OVERLAP \
                                                                             $RESULT1/CHR_${CHR}
    )                                          
   echo "$JOB"
}

#PHASE COMMON ALL
function phase_common_iter {

module load igmm/apps/bcftools/1.20 
    STEP=500000
    OVERLAP=20000
    MEM=16G #16 IS THE MINI SIZE OF A CPU.IT WILL NOT GIVE YOU 'HALP A CPU' ANYWAY
    for CHR in {4..12}; do

        END=$(bcftools view $RESULT0/genome1.focal.chr${CHR}.vcf.bgz |head -50|grep '^##contig='|cut -d ',' -f 2 |cut -d '=' -f 2)
        mkdir $RESULT1/CHR_${CHR}
        echo "lengh of CHR${CHR} = ${END}"
        phase_common $CHR 1 $END $STEP $OVERLAP $MEM                          

    done

}

function phase_common_fix {

STEP=500000
OVERLAP=20000

#for f in phase_common.e*; do     if grep -qE "bad_alloc|Segmentation fault|Aborted" "$f" 2>/dev/null; then         echo "FAILED: $f";     fi; done

CHR=4
#phase_common $CHR 1000001 $((3500001+$STEP -1)) $STEP $OVERLAP 6G
#phase_common $CHR 8000001 $((8000001+$STEP -1)) $STEP $OVERLAP 6G 
CHR=6 
#phase_common $CHR 22500001 $((22500001+$STEP -1)) $STEP $OVERLAP 6G  
CHR=7
#phase_common $CHR 14000001 $((14000001+$STEP -1)) $STEP $OVERLAP 6G  
CHR=11
#phase_common $CHR 6500001 $((6500001+$STEP -1)) $STEP $OVERLAP 6G  
#phase_common $CHR 28000001 $((28000001+$STEP -1)) $STEP $OVERLAP 6G  
CHR=12
#phase_common $CHR 13000001 $((13000001+$STEP -1)) $STEP $OVERLAP 6G
#phase_common $CHR 16000001 $((16000001+$STEP -1)) $STEP $OVERLAP 6G 
#phase_common $CHR 16500001 $((16500001+$STEP -1)) $STEP $OVERLAP 6G 
#phase_common $CHR 20000001 $((20000001+$STEP -1)) $STEP $OVERLAP 6G 
phase_common $CHR 6000001 $((6000001+$STEP -1)) $STEP $OVERLAP 12G     
#phase_common $CHR 8000001 $((8500001+$STEP -1)) $STEP $OVERLAP 6G     
}
#LIGATE CHUNKS
function ligate_chunks {

JOB=$(qsub -t 11-12 \
               -pe sharedmem 1 \
               -l h_vmem=4G \
               -N "ligate" \
               -l h='!(node1e06|node1f07)'  \
               -o $HOME/Phasing \
               -e $HOME/Phasing \
               $HOME/Phasing/eddie_phasing_shapeit_functions.sh ligate $RESULT1 $RESULT2
    ) 
    echo "$JOB"
}

#PHASE RARE
function phase_rare {

    CHR=$1
    START=$2
    END=$3
    STEP=$4
    OVERLAP=$5
    MEM=$6

    echo "PHASE RARE CHR ${CHR}: ${START}-${END} WITH STEP=${STEP} OVERLAP=${OVERLAP}"
    JOB=$(qsub -t ${START}-$END:$STEP \
               -pe sharedmem 1 \
               -l h_vmem=${MEM} \
               -l h='!(node1e06|node1e07|node1e09)' \
               -N "rare_${CHR}" \
               -o $WORKDIR/run_rare \
               -e $WORKDIR/run_rare \
               $HOME/Phasing/eddie_phasing_shapeit_functions.sh phase_rare $RESULT0/genome1.focal.chr${CHR}.vcf.bgz \
                                                                           $RESULT2/chr_${CHR}_scaffold.bcf  \
                                                                           $CHR \
                                                                           $END \
                                                                           $STEP \
                                                                           $OVERLAP \
                                                                           $RESULT3/CHR_${CHR}
    )                                          
    echo "$JOB"
    #-l h='(node1d01|node1d02|node1d03|node1d10|node1e01|node1e02|node1e03|node1e04|node1e05|node1f01|node1f03|node1f04|node1f16)'  \
    #-l h='!(node1e06|node1e07|node1e09)' \
}

#PAHSE ALL RARE ALLELES IN BULK
function phase_rare_iter {

    module load igmm/apps/bcftools/1.20 
    STEP=50000
    OVERLAP=10000
    MEM=16G
#need to do 10 11 12
    for CHR in 10 11 12  ; do
        END=$(bcftools view $RESULT0/genome1.focal.chr${CHR}.vcf.bgz |head -50|grep '^##contig='|cut -d ',' -f 2 |cut -d '=' -f 2)
        mkdir $RESULT3/CHR_${CHR}
        phase_rare $CHR 1 $END $STEP $OVERLAP $MEM
    done   
}

#MANUALLY FIX A FEW REGIONS
function fix_chr1_rare {

#for f in phase_com.e*; do     if grep -qE "command not found| bad_alloc|Segmentation fault|Aborted" "$f" 2>/dev/null; then         echo "FAILED: $f";     fi; done
#for f in  /exports/eddie/scratch/zchen13/Phasing/run_rare/rare_8.e*; do if grep -qE "void genotype_set" "$f" 2>/dev/null; then         echo "$f";     fi; done
#chr 9 passed without error
    phase_rare 4 200000 300000 5000 2000
}

function concatenate {

JOB=$(qsub -t 5 \
               -pe sharedmem 1 \
               -l h_vmem=32G \
               -N "concatenate" \
               -l h='!(node1e06|node1f07|node1f09)'  \
               -o $HOME/Phasing \
               -e $HOME/Phasing \
               $HOME/Phasing/eddie_phasing_shapeit_functions.sh concatenate $RESULT3 $RESULT4
    ) 
    echo "$JOB"
}

function add_AA {

    echo '##INFO=<ID=AA,Number=1,Type=String,Description="Ancestral Allele">' > $RESULT5/aa_header.txt
    JOB=$(qsub -t 9 \
               -pe sharedmem 1 \
               -l h_vmem=16G \
               -N "ADD_AA" \
               -l h='!(node1e06|node1f07|node1f09)'  \
               -o $HOME/Phasing \
               -e $HOME/Phasing \
               $HOME/Phasing/eddie_phasing_shapeit_functions.sh add_AA $AA conca.bcf $RESULT4 $RESULT5
    ) 
    echo "$JOB"
}

#CONVERT TO VCZ FORMAT FOR tsinfer/tskit
function convert_vcz {

    echo '##INFO=<ID=AA,Number=1,Type=String,Description="Ancestral Allele">' > $RESULT5/aa_header.txt
    JOB=$(qsub -t 9 \
               -pe sharedmem 16 \
               -l h_vmem=20G \
               -N "VCZ_CONVERT" \
               -l h='!(node1e06|node1f07|node1f09)'  \
               -o $HOME/Phasing \
               -e $HOME/Phasing \
               $HOME/Phasing/eddie_phasing_shapeit_functions.sh convert_vcz conca.bcf $RESULT5 $RESULT6
    ) 
    echo "$JOB"
}

function tsinfer {

JOB=$(qsub -t 9 \
               -pe sharedmem 64 \
               -l h_vmem=40G,h_rss=40G \
               -N "ts_infer_WG" \
               -l h='!(node1e06|node1f07|node1f09)'  \
               -o $HOME/Phasing \
               -e $HOME/Phasing \
               $HOME/Phasing/eddie_phasing_shapeit_functions.sh python_tsinfer $RESULT6 $RESULT7
    ) 
    echo "$JOB"

}
#================================

function main {

#CONVERT-PLINK 
#convert_plink
#
#phase_common_iter
#phase_common_fix
#ligate_chunks
#phase_rare_iter
#concatenate
#add_AA
#convert_vcz
#
tsinfer

}

main
