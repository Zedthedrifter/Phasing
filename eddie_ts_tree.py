import msprime
import zarr
import tsinfer
import tskit 
from contextlib import redirect_stdout
import subprocess as sbp
import sys
import argparse
import csv
from Bio import SeqIO
from Bio.Seq import Seq
from Bio import Phylo
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np

import zarr
import numpy as np
import tsinfer

RateMap_dir='/exports/eddie/scratch/zchen13/DATA/RIS_reference_panel'

chr_length={1:43270923,
            2:35937250,
            3:36413819,
            4:35502694,
            5:29958434,
            6:31248787,
            7:29697621,
            8:28443022,
            9:23012720,
            10:23207287,
            11:29021106,
            12:27531856}

def run_tsinfer(vcz, outfile, chr, start=0,end=None):
    #MAKE RECOMBINATION MAP
    map=f"{RateMap_dir}/ReferencePanel_4.8M_SNP_chr{chr}updated.map"
    # 1. Read your map file (tab-separated)
    print(f"msprime read recombination map CHR {chr} for tsinfer")
    df = pd.read_csv(map, sep='\t')

    # 2. Extract positions and rates
    # RateMap requires position length = rate length + 1. Your file's positions are the start of each interval.
    # You need to append an end coordinate (e.g., chromosome length) as the final position.
    positions = df['position'].values
    positions = np.concatenate(([0], positions))
    rates = df['COMBINED_rate.cM.Mb.'].values
    print(chr, chr_length[chr],len(positions),len(rates))
    # Assuming chromosome 1 ends at the last position (or a known length)
    # Note: the rate unit is cM/Mb. msprime expects rate in "per base per generation" (Morgans/bp).
    # Conversion: 1 cM/Mb = 0.01 Morgan / 1,000,000 bp = 1e-8 per bp.
    rates_in_morgans = rates * 1e-8

    # 3. Construct RateMap
    # Ensure length(positions) == length(rates_in_morgans) + 1
    # If your positions do not include the chromosome start point 0, you may need to manually prepend 0.
    rate_map = msprime.RateMap(position=positions, rate=rates_in_morgans)
    #
    #INFER TREE
    vcf_zarr = zarr.open(vcz)
    aa_array = vcf_zarr["variant_AA"][:]
    valid_bases = {'A', 'C', 'G', 'T'}
    cleaned_aa = np.array([
        base if base in valid_bases else '.'
        for base in aa_array
    ])
    print("Cleaned AA set:", set(cleaned_aa))
    n_sites = vcf_zarr["variant_position"].shape[0]
    print(f"CHR {chr} contains {n_sites} variants")
    if end is not None:
        print(f"make vdata for variants {start}:{end}")
        site_mask = np.ones(n_sites, dtype=bool)
        site_mask[start:end] = False
        vdata = tsinfer.VariantData(
            vcz,
            cleaned_aa[start:end],
            site_mask=site_mask
        )
    else:
        vdata = tsinfer.VariantData(vcz, cleaned_aa)
    print(type(vdata))
    ts = tsinfer.infer(vdata, 
                       recombination_rate=rate_map,
                       num_threads=64,
                       path_compression=True,
                       post_process=True,
                       progress_monitor=True)
    #mmap_temp_dir="/exports/eddie/scratch/zchen13/tmp",
    ts.dump(outfile)
    print(f"Saved tree sequence to {outfile}")
    print(f"Tree sequence contains {ts.num_trees} trees, {ts.num_sites} sites, "
          f"{ts.num_samples} samples")
    return vdata, ts

#=============================
def main():
    count=0
    parser = argparse.ArgumentParser(description="Python script with multiple executable functions")
    subparsers = parser.add_subparsers(dest='command', help='Available commands')
        
    # Common arguments that will be reused
    common_args = argparse.ArgumentParser(add_help=False)
    common_args.add_argument('-v','--vcf', help='vcf file with GT and DP info',metavar='vcf')
    common_args.add_argument('-o','--outfile', help='full path for outputfile',metavar='outfile')
    common_args.add_argument('--start', type=int,default=0,help='use from n nt of the chromosome',metavar='start')
    common_args.add_argument('--end', type=int,default=None,help='use up to n nt of the chromosome',metavar='end')
    common_args.add_argument('--csv', help='csv file, high resolution region, species specific SNPs, etc',metavar='csv')
    common_args.add_argument('--sp', help='species to be simulated',metavar='sp')
    common_args.add_argument('--model', help='model to use for the simulation',metavar='MODEL')
    common_args.add_argument('--chr', type=int,help='chromosome',metavar='chr')
    common_args.add_argument('--vcz', help='vcz file (real or from simulation)',metavar='vcz')

    #functions
    #tsinfer
    run_tsinfer_parser=subparsers.add_parser('run_tsinfer', parents=[common_args],help='run_tsinfer', 
                                            usage = './ts.py run_tsinfer--vcz <vcz file> --outfile <svg output>')
    run_tsinfer_parser.set_defaults(func=run_tsinfer)
    #TSKIT

    #parse arguments
    args = parser.parse_args()

    if not args.command:
        parser.print_help()
        sys.exit(1)
        
    # Prepare common kwargs for the function call
    kwargs = {k: v for k, v in vars(args).items() 
                if k in ['vcf','csv', 'sp','MODEL','chr','outfile','vcz', 'start','end'] 
                and v is not None}
        
    # Call the appropriate function
    args.func(**kwargs)
    
if __name__ == '__main__': 
    main()

