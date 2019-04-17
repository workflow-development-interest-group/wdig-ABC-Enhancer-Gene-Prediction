import argparse
import progressbar as pb
from predictor import Predictor
from tools import *
import pandas as pd
import numpy as np
import sys, traceback, os, os.path

def get_model_argument_parser():
    class formatter(argparse.ArgumentDefaultsHelpFormatter, argparse.RawTextHelpFormatter):
        pass

    epilog = ("model is:\n"
              "\tactivity = (H3K27ac * DHS ^ Ratio) ^ (1 / (1 + Ratio))\n"
              "\tscore = (activity ^ A.beta) * (proximity ^ P.beta)\n"
              "\trelative effect = score / (sum of scores in neighborhood)")


    parser = argparse.ArgumentParser(description='Predict enhancer relative effects.',
                                     epilog=epilog,
                                     formatter_class=formatter)
    readable = argparse.FileType('r')

    #Basic parameters
    parser.add_argument('--cellType', required=False, help="Name of cell type")
    parser.add_argument('--nbhd_directory', help="Directory with neighborhoods files")
    parser.add_argument('--outdir', required=True, help="output directory")
    parser.add_argument('--params_file', help="Parameters file")
    parser.add_argument('--genes', type=readable, required=False, help="Table of genes for which predictions should be made. Overrides GeneList.txt in neighborhoods directory")
    parser.add_argument('--HiC_directory_listing', default="", help="File mapping cell names to hic directories")
    parser.add_argument('--window', type=int, default=5000000, help="Make predictions for all enhancers within this distance of the gene's TSS")
    parser.add_argument('--threshold', type=float, required=True, default=None, help="Threshold on ABC Score to call a predicted positive")

    #qnorm
    parser.add_argument('--qnorm', default='', help="json file containing to quantile normalize epigenetica data to")

    #hic
    parser.add_argument('--HiCdir', help="Directory with hic bedgraphs")
    parser.add_argument('--hic_cap', type=float, default=100, help="HiC cap (in normalized units 0-100)")
    parser.add_argument('--tss_hic_contribution', type=float, default=100, help="Weighting of diagonal bin of hic matrix as a percentage of its neighbors")
    parser.add_argument('--hic_pseudocount_distance', type=int, default=1e6, help="A pseudocount is added equal to the powerlaw fit at this distance")

    #Power law
    parser.add_argument('--scale_hic_using_powerlaw', action="store_true", help="Scale Hi-C values using powerlaw relationship")
    parser.add_argument('--hic_gamma', type=float, default=1, help="powerlaw exponent of hic_cell_type. Must be positive")
    parser.add_argument('--hic_gamma_reference', type=float, default=1, help="powerlaw exponent to scale to. Must be positive")

    #Genes to run through model
    parser.add_argument('--run_all_genes', action='store_true', help="Do not check for gene expression, make predictions for all genes")
    parser.add_argument('--expression_cutoff', type=float, default=1, help="Make predictions for genes with expression higher than this value")
    parser.add_argument('--promoter_activity_quantile_cutoff', type=float, default=.4, help="Quantile cutoff on promoter activity")

    #Output formatting
    parser.add_argument('--skip_gene_files', action="store_true", help="Do not make individual gene files")
    parser.add_argument('--skinny_gene_files', action="store_true", help="Use subset of columns for genes files")
    parser.add_argument('--make_all_putative', action="store_true", help="Make big file with concatenation of all genes file")

    #Other
    parser.add_argument('--tss_slop', type=int, default=500, help="Distance from tss to search for self-promoters")
    parser.add_argument('--include_chrY', '-y', action='store_true', help="Include Y chromosome")

    return parser


def get_predict_argument_parser():
    parser = get_model_argument_parser()
    return parser


def parse_cell_type_args(args, cellType):
    file_params = pd.read_table(args.params_file, sep='\t')
    file_params = file_params.loc[file_params["cell_type"] == cellType, ]
    args.DHS_column = file_params['default_accessibility_feature'].values[0] + ".RPM"

    args.enhancers = os.path.join(args.nbhd_directory, "EnhancerList.txt")   
    if args.genes is None:
        args.genes = os.path.join(args.nbhd_directory, "GeneList.txt")
     
    return args

def main():
    parser = get_predict_argument_parser()
    args = parser.parse_args()
    args = parse_cell_type_args(args, args.cellType)

    if not os.path.exists(args.outdir):
        os.makedirs(args.outdir)

    preddir = os.path.join(args.outdir, "genes")
    if not os.path.exists(preddir):
        os.makedirs(preddir)

    write_prediction_params(args, os.path.join(args.outdir, "parameters.predict.txt"))
    
    print("reading genes")
    genes = read_genes(args.genes)

    print("reading enhancers")
    enhancers = read_enhancers(args.enhancers)

    print("building predictor")
    predictor = Predictor(enhancers, **vars(args))

    print("applying qnorm")
    predictor.add_normalized_data_to_enhancers(enhancers)

    chromosomes = predictor.chromosomes()
    print("data loaded for chromosomes: {}".format(" ".join(sorted(chromosomes))))

    #Initialize Prediction files
    pred_file = os.path.join(args.outdir, "EnhancerPredictions.txt")
    all_pred_file = os.path.join(args.outdir, "EnhancerPredictionsAllPutative.txt.gz")
    pred_bed = os.path.join(args.outdir, "EnhancerPredictions.bed")

    args.score_column = "ABC.Score"

    all_positive_list = []
    all_putative_list = []
    gene_stats = []
    failed_genes = []
    pbar = pb.ProgressBar(max_value=len(genes), redirect_stdout=True)
    for idx, gene in pbar(genes.iterrows()):
        if gene.chr == 'chrY' and not args.include_chrY:
            continue
        if gene.chr not in chromosomes:
            print("\nNo data for {}".format(gene.chr))
            continue
        print("\nPredicting {} with {} {} TSS".format(gene["name"], gene["chr"], gene["tss"]))

        try:
            nearby_enhancers = enhancers.within_range(gene.chr, gene.tss - args.window, gene.tss + args.window)
            predictor.predict_from_normalized_to_enhancers(nearby_enhancers, gene, args.window, tss_slop=args.tss_slop)
            
            col_names=['chr','start','end','TargetGene','TargetGeneTSS','class','Score.Fraction','Score','distance','hic.distance','hic.distance.adj','estimatedCP','estimatedCP.adj','normalized_dhs','normalized_h3k27ac','TargetGeneExpression','TargetGeneTSSActivityQuantile']
            if not args.skip_gene_files:
                if not args.skinny_gene_files:
                    write_scores(preddir, gene, nearby_enhancers)
                else:
                    write_scores(preddir, gene, nearby_enhancers[col_names])

            if args.make_all_putative:
                all_putative_list.append(nearby_enhancers[col_names])
            
            gene_is_expressed_proxy = check_gene_for_runnability(gene, args.expression_cutoff, args.promoter_activity_quantile_cutoff)
            if args.run_all_genes or gene_is_expressed_proxy:
                positives = nearby_enhancers.ix[nearby_enhancers[args.score_column] >= args.threshold,:]
                print("{} enhancers predicted for {}".format(positives.shape[0], gene["name"]))
                all_positive_list.append(positives)

            #Add gene to gene summary file
            if nearby_enhancers.shape[0] > 0:
                stats = predictor.get_gene_prediction_stats(args, nearby_enhancers)
                stats['prediction_file'] = get_score_filename(gene)
                stats['gene_is_expressed_proxy'] = gene_is_expressed_proxy
                gene_stats.append(stats)

        except:
            failed_genes.append(gene["chr"] + "\t" + gene["name"])
            print("Failed on " + gene["name"] + " ... skipping. Traceback:")
            traceback.print_exc(file=sys.stdout)
            continue

    if args.score_column is not None:
        all_positive = pd.concat(all_positive_list)
        all_positive.to_csv(pred_file, sep="\t", index=False, header=True, float_format="%.4f")
        write_connections_bedpe_format(all_positive.loc[all_positive["class"] != "promoter"], outfile=os.path.join(args.outdir, "Predictions_nopromoters.bedpe"), score_column=args.score_column)

    gene_stats = pd.concat(gene_stats, axis=1).T
    gene_stats.to_csv(os.path.join(args.outdir, "GenePredictionStats.txt"), sep="\t", index=False)

    with open(os.path.join(args.outdir, "FailedGenes.txt"), 'w') as failed_file:
        for gene in failed_genes:
            failed_file.write(gene + "\n")

    if args.make_all_putative:
        all_putative = pd.concat(all_putative_list)
        all_putative.to_csv(all_pred_file, sep="\t", index=False, header=True, compression="gzip", float_format="%.4f", na_rep="NaN", chunksize=100000)

def write_prediction_params(args, file):
    with open(file, 'w') as outfile:
        for arg in vars(args):
            outfile.write("--" + arg + " " + str(getattr(args, arg)) + " ")


if __name__ == '__main__':
    main()
