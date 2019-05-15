import argparse
import numpy as np
import pandas as pd

def parseargs(required_args=True):
    class formatter(argparse.ArgumentDefaultsHelpFormatter, argparse.RawTextHelpFormatter):
        pass

    epilog = ("")
    parser = argparse.ArgumentParser(description='Run neighborhood for a given cell type',
                                     epilog=epilog,
                                     formatter_class=formatter)

    parser.add_argument('--enhancers', required=required_args, help="EnhancerList file")
    parser.add_argument('--outfile', default="", help="columns for which to compute qnorm reference")
    parser.add_argument('--cols', default="DHS.RPM,H3K27ac.RPM,ATAC.RPM", help="columns for which to compute qnorm reference")

    args = parser.parse_args()
    return args


def makeQnorm(args):
	enhancers = pd.read_csv(args.enhancers, sep="\t")
	ref = pd.DataFrame({'quantile' : np.concatenate([np.linspace(0,.99,100), np.linspace(.991,.999,9)])})
	ref['rank'] = (ref['quantile'] * enhancers.shape[0]).round(decimals=0)

	cols = set(set(vars(args)['cols'].split(",")) & set(enhancers.columns))
	for col in cols:
		ref[col] = np.percentile(enhancers[col], ref['quantile']*100)

	ref.to_csv(args.outfile, sep="\t", index=False, float_format="%.3f")

def main(args):
    makeQnorm(args)

if __name__ == '__main__':
    args = parseargs()
    main(args)