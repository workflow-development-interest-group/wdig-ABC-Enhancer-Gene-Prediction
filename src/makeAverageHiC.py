from hic import *
import numpy as np
from functools import reduce
import argparse
import sys, os, os.path
from tools import write_params

def parseargs():
    class formatter(argparse.ArgumentDefaultsHelpFormatter, argparse.RawTextHelpFormatter):
        pass

    epilog = ("")
    parser = argparse.ArgumentParser(description='Make average HiC dataset',
                                     epilog=epilog,
                                     formatter_class=formatter)
    readable = argparse.FileType('r')
    parser.add_argument('--celltypes', required=True, help="Comma delimitted list of cell types")
    parser.add_argument('--basedir', required=True, help="Basedir")
    parser.add_argument('--outDir', required=True, help="Output directory")
    parser.add_argument('--resolution', default=5000, type=int, help="Resolution of hic dataset (in bp)")
    parser.add_argument('--ref_scale', default=-1.91, type=float, help="Resolution of hic dataset (in bp)")
    parser.add_argument('--ref_gamma', default=-.876, type=float, help="Resolution of hic dataset (in bp)")

    args = parser.parse_args()
    return(args)

def main():
    args = parseargs()
    os.makedirs(args.outDir, exist_ok=True)

    #Write params file
    write_params(args, os.path.join(args.outDir, "params.txt"))

    #Parse cell types
    cell_types = args.celltypes.split(",")

    chromosomes = ['chr' + str(x) for x in range(1,23)] + ['chrX'] 
    chromosomes = ['chr22']

    for chromosome in chromosomes:
        hic_list = [process_chr(cell_type, chromosome, args.basedir, args.resolution, args.ref_scale, args.ref_gamma) for cell_type in cell_types]

        #Make average
        #TO DO: Enforce minimum number of cell types required for averaging
        all_hic = reduce(lambda x, y: pd.merge(x, y, on = ['bin1', 'bin2'], how = 'outer'), hic_list)
        all_hic.fillna(value=0, inplace=True)
        cols_for_avg = list(filter(lambda x:'hic_kr_' in x, all_hic.columns))
        all_hic['avg_hic'] = all_hic[cols_for_avg].mean(axis=1)
        all_hic.drop(cols_for_avg, inplace=True, axis=1)

        all_hic['bin1'] = all_hic['bin1'] * args.resolution
        all_hic['bin2'] = all_hic['bin2'] * args.resolution

        os.makedirs(os.path.join(args.outDir, chromosome), exist_ok=True)
        all_hic.to_csv(os.path.join(args.outDir, chromosome, chromosome + ".KRobserved.gz"), sep="\t", header=False, index=False, compression="gzip")

def scale_hic_with_powerlaw(hic, resolution, scale_ref, gamma_ref, scale, gamma):

    #get_powerlaw_at_distance expects positive gamma
    gamma_ref = -1 * gamma_ref
    gamma = -1 * gamma

    dists = (hic['bin2'] - hic['bin1']) * resolution
    pl_ref = get_powerlaw_at_distance(dists, gamma_ref, scale_ref) 
    pl = get_powerlaw_at_distance(dists, gamma, scale)

    hic['hic_kr'] = hic['hic_kr'] * (pl_ref / pl)

    return hic

def process_chr(cell_type, chromosome, basedir, resolution, scale_ref, gamma_ref):

    # import pdb
    # pdb.set_trace()

    hic_file = os.path.join(basedir, cell_type, "5kb_resolution_intra", chromosome, chromosome + ".KRobserved")

    #Load gamma and scale
    pl_summary = pd.read_csv(os.path.join(basedir, cell_type, "5kb_resolution_intra/powerlaw/hic.powerlaw.txt"), sep="\t")

    #Read in and normalize to make DS
    hic = load_hic(hic_file, 
                    hic_type="juicebox", 
                    hic_resolution = resolution, 
                    tss_hic_contribution=np.NaN, 
                    window = np.Inf, 
                    min_window=0, 
                    gamma = -1*pl_summary['pl_gamma'].values[0],
                    apply_diagonal_bin_correction=False)

    #power law scale 
    hic = scale_hic_with_powerlaw(hic, resolution, scale_ref, gamma_ref, scale = pl_summary['pl_scale'].values[0], gamma = pl_summary['pl_gamma'].values[0])

    return hic


if __name__ == '__main__':
    main()



    #Loop over each chromosome
    #Reach HiC file from each cell type for this chromosome
    #KR Normalize
    #Powerlaw scale
    #Average together
    #Write average file
