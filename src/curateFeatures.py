import pandas as pd
import numpy as np
import argparse
import os
from peaks import *
import traceback
from itertools import chain

#TO DO
#1. Assumes genome build is hg19

pd.options.display.max_colwidth = 100000 #is this necessary... (really?)

def parseargs(required_args=True):
    class formatter(argparse.ArgumentDefaultsHelpFormatter, argparse.RawTextHelpFormatter):
        pass

    epilog = ("")
    parser = argparse.ArgumentParser(description='Make peaks file for a given cell type',
                                     epilog=epilog,
                                     formatter_class=formatter)
    readable = argparse.FileType('r')
    parser.add_argument('--cellType', required=required_args, help="Name of cell type")
    parser.add_argument('--params_file', required=required_args, help="Parameters file")
    parser.add_argument('--genome', required=required_args, help="File listing genome annotations for each species/build")
    parser.add_argument('--outDir', required=required_args)
    
    parser.add_argument('--pval_cutoff', default=.1, help="pvalue cutoff for MACS2")
    parser.add_argument('--nStrongestPeaks', default=175000, help="Number of peaks to use for defining candidate regions")
    parser.add_argument('--peakExtendFromSummit', default=250, help="Number of base pairs to extend each preak from its summit")

    parser.add_argument('--regions_whitelist', default="", help="Bed file of regions to forcibly include in candidate enhancers. Overrides regions_blacklist")
    parser.add_argument('--regions_blacklist', default="", help="Bed file of regions to forcibly exclude from candidate enhancers")
    
    args = parser.parse_args()
    return(args)

def processCellType(cellType, args):
	os.makedirs(os.path.join(args.outDir), exist_ok=True)
	write_params(args, os.path.join(args.outDir, "params.txt"))

	params = pd.read_csv(args.params_file, sep="\t")
	assert(cellType in params["cell_type"].tolist()), "Celltype not found in parameters file"
	params = params.loc[params["cell_type"] == cellType, :]

	genome_params = pd.read_csv(args.genome, sep="\t").set_index("name").T.to_dict()
	genome = genome_params[params['genome'].values[0]]

	#Note: Will run both DHS and ATAC files if they both exist
	d1 = {filename : "DHS" for filename in params['feature_DHS'].to_string(index=False).split(",") if filename != 'NaN'}
	d2 = {filename : "ATAC" for filename in params['feature_ATAC'].to_string(index=False).split(",") if filename != 'NaN'}
	d3 = {filename : "H3K27ac" for filename in params['feature_H3K27ac'].to_string(index=False).split(",") if filename != 'NaN'}
	file_dict = dict(chain.from_iterable(d.items() for d in (d1, d2, d3)))

	qc_list = []
	for this_file, feature_type in file_dict.items():
		qc_stats = pd.DataFrame({'cellType' : [cellType]})
		qc_stats['feature'] = feature_type
		qc_stats['default_accessibility_feature'] = params["default_accessibility_feature"].values[0]
		qc_stats['file'] = this_file

		#MACS
		if feature_type != "H3K27ac":
			macs_format = get_macs_format(feature_type) #assumes ATAC is paired end and DNase is not

			peaks_file_prefix = os.path.basename(this_file.replace(".tagAlign.gz", ".macs2").replace(".bam", ".macs2"))
			run_macs2_callpeak(this_file, args.outDir, peaks_file_prefix, macs_format, args.pval_cutoff, genome['sizes'], force=False)	
			qc_stats['peak_file'] = os.path.join(args.outDir, peaks_file_prefix + '_peaks.narrowPeak')
			qc_stats['num_peaks'] = compute_macs_stats(os.path.join(args.outDir, peaks_file_prefix + '_peaks.narrowPeak'))

			#Make candidate regions
			make_candidate_regions_from_summits(macs_peaks = qc_stats['peak_file'].values[0], 
												accessibility_file = this_file, 
												genome_sizes = genome['sizes'], 
												regions_whitelist = args.regions_whitelist,
												regions_blacklist = args.regions_blacklist,
												n_enhancers = args.nStrongestPeaks, 
												peak_extend = args.peakExtendFromSummit, 
												outdir = args.outDir)
			qc_stats['candidate_region_file'] = os.path.join(args.outDir, os.path.basename(qc_stats['peak_file'].values[0]) + ".candidateRegions.bed")

		#Count Reads 
		try:
			qc_stats['read_count'] = count_total(this_file)
		except Exception as e:
			print(e)
			qc_stats['read_count'] = np.nan

		qc_list.append(qc_stats)

	all_qc_stats = pd.concat(qc_list)
	all_qc_stats.to_csv(os.path.join(args.outDir, "feature.stats.txt"), sep="\t", index=False)
	print(all_qc_stats)
	
def write_params(args, file):
    with open(file, 'w') as outfile:
        for arg in vars(args):
            outfile.write("--" + arg + " " + str(getattr(args, arg)) + " ")

def main(args):
    processCellType(args.cellType, args)

if __name__ == '__main__':
    args = parseargs()
    main(args)




