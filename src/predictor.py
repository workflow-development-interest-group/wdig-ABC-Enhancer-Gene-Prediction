import numpy as np
import pandas as pd
from tools import get_gene_name
import sys
import time
import pyranges as pr
from hic import *

def get_hic_file(chromosome, hic_dir):
    hic_file = os.path.join(hic_dir, chromosome, chromosome + ".KRobserved")

    assert(len(hic_file) == 1)

    return(hic_file)

def make_predictions(chromosome, enhancers, genes, hic_file, args):
    pred = make_pred_table(chromosome, enhancers, genes, args)
    pred = add_hic_to_enh_gene_table(enhancers, genes, pred, hic_file, chromosome, args)
    pred = annotate_predictions(pred)

    pred = compute_score(pred, [pred['activity_base'], pred['hic_kr_pl_scaled_adj']], "ABC")
    #pred = compute_score(pred, [pred['activity_base'], pred['estimatedCP.adj']], "powerlaw")

    return(pred)

def make_pred_table(chromosome, enh, genes, args):
    print('Making putative predictions table...')
    t = time.time()
 

    enh['enh_midpoint'] = (enh['start'] + enh['end'])/2
    enh['enh_idx'] = enh.index
    genes['gene_idx'] = genes.index
    enh_pr = df_to_pyranges(enh)
    genes_pr = df_to_pyranges(genes, start_col = 'TargetGeneTSS', end_col = 'TargetGeneTSS', start_slop=args.window, end_slop = args.window)

    pred = enh_pr.join(genes_pr).df.drop(['Start_b','End_b','chr_b','Chromosome','Start','End'], axis = 1)
    #pred['enh_midpoint'] = (pred['start'] + pred['end'])/2
    pred['distance'] = abs(pred['enh_midpoint'] - pred['TargetGeneTSS'])
    pred = pred.loc[pred['distance'] < args.window,:] #for backwards compatability

    # else:
    #     enh['temp_merge_key'] = 0
    #     genes['temp_merge_key'] = 0

    #     #Make cartesian product and then subset to EG pairs within window. 
    #     #TO DO: Replace with pyranges equivalent of bedtools intersect or GRanges overlaps 
    #     pred = pd.merge(enh, genes, on = 'temp_merge_key')

    #     pred['enh_midpoint'] = (pred['start'] + pred['end'])/2
    #     pred['distance'] = abs(pred['enh_midpoint'] - pred['TargetGeneTSS'])
    #     pred = pred.loc[pred['distance'] < args.window,:]

    #     print('Done. There are {} putative enhancers for chromosome {}'.format(pred.shape[0], chromosome))
    #     print('Elapsed time: {}'.format(time.time() - t))

    return pred
    
def df_to_pyranges(df, start_col='start', end_col='end', chr_col='chr', start_slop=0, end_slop=0):
    df['Chromosome'] = df[chr_col]
    df['Start'] = df[start_col] - start_slop
    df['End'] = df[end_col] + end_slop

    return(pr.PyRanges(df))

def add_hic_to_enh_gene_table(enh, genes, pred, hic_file, chromosome, args):
    print('Begin HiC')
    HiC = load_hic(hic_file = hic_file, hic_type = args.hic_type, hic_resolution = args.hic_resolution, tss_hic_contribution = args.tss_hic_contribution, window = args.window, min_window = 0, gamma = args.hic_gamma)

    # import pdb
    # pdb.set_trace()

    #Add hic to pred table
    #At this point we have a table where each row is an enhancer/gene pair. 
    #We need to add the corresponding HiC matrix entry.
    #If the HiC is provided in juicebox format (ie constant resolution), then we can just merge using the indices
    #But more generally we do not want to assume constant resolution. In this case hic should be provided in bedpe format

    t = time.time()
    if args.hic_type == "bedpe":
        #Use pyranges to compute overlaps between enhancers/genes and hic bedpe table
        #Consider each range of the hic matrix separately - and merge each range into both enhancers and genes. 
        #Then remerge on hic index

        HiC['hic_idx'] = HiC.index
        hic1 = df_to_pyranges(HiC, start_col='x1', end_col='x2')
        hic2 = df_to_pyranges(HiC, start_col='y1', end_col='y2')

        #Overlap in one direction
        enh_hic1 = df_to_pyranges(enh, start_col = 'enh_midpoint', end_col = 'enh_midpoint', end_slop = 1).join(hic1).df
        genes_hic2 = df_to_pyranges(genes, start_col = 'TargetGeneTSS', end_col = 'TargetGeneTSS', end_slop = 1).join(hic2).df
        ovl12 = enh_hic1[['enh_idx','hic_idx','hic_kr']].merge(genes_hic2[['gene_idx', 'hic_idx']], on = 'hic_idx')

        #Overlap in the other direction
        enh_hic2 = df_to_pyranges(enh, start_col = 'enh_midpoint', end_col = 'enh_midpoint', end_slop = 1).join(hic2).df
        genes_hic1 = df_to_pyranges(genes, start_col = 'TargetGeneTSS', end_col = 'TargetGeneTSS', end_slop = 1).join(hic1).df
        ovl21 = enh_hic2[['enh_idx','hic_idx','hic_kr']].merge(genes_hic1[['gene_idx', 'hic_idx']], on = ['hic_idx'])

        #Concatenate both directions and merge into preditions
        ovl = pd.concat([ovl12, ovl21]).drop_duplicates() #.drop('hic_idx', axis = 1)
        pred = pred.merge(ovl, on = ['enh_idx', 'gene_idx'], how = 'left')
        pred.fillna(value={'hic_kr' : 0}, inplace=True)
    elif args.hic_type == "juicebox":
        #Merge directly using indices
        #Could also do this by indexing into the sparse matrix (instead of merge) but this seems to be slower
        #Index into sparse matrix
        #pred['hic_kr'] = [HiC[i,j] for (i,j) in pred[['enh_bin','tss_bin']].values.tolist()]
        pred['enh_bin'] = np.floor(pred['enh_midpoint'] / args.hic_resolution).astype(int)
        pred['tss_bin'] = np.floor(pred['TargetGeneTSS'] / args.hic_resolution).astype(int)
        pred['bin1'] = np.amin(pred[['enh_bin', 'tss_bin']], axis = 1)
        pred['bin2'] = np.amax(pred[['enh_bin', 'tss_bin']], axis = 1)
        pred = pred.merge(HiC, how = 'left', on = ['bin1','bin2'])
        pred.fillna(value={'hic_kr' : 0}, inplace=True)


    pred.drop(['x1','x2','y1','y2','bin1','bin2','enh_idx','gene_idx','hic_idx','enh_midpoint','tss_bin','enh_bin'], inplace=True, axis = 1, errors='ignore')
        
    print('HiC added to predictions table. Elapsed time: {}'.format(time.time() - t))

    # Add powerlaw scaling
    pred = scale_with_powerlaw(pred, args)

    #Add pseudocount
    pred = add_hic_pseudocount(pred, args)

    print("HiC Complete")
    #print('Elapsed time: {}'.format(time.time() - t))

    return(pred)

def scale_with_powerlaw(pred, args):

    #TO DO: is this np.exp right? Shouldn't it be dependant on distance?
    if not args.scale_hic_using_powerlaw:
        pred['hic_kr_pl_scaled'] = pred['hic_kr']
    else:
        powerlaw_estimate = get_powerlaw_at_distance(pred['distance'].values, args.hic_gamma)
        powerlaw_estimate_reference = get_powerlaw_at_distance(pred['distance'].values, args.hic_gamma_reference)
        pred['powerlaw_contact'] = powerlaw_estimate
        pred['powerlaw_contact_reference'] = powerlaw_estimate_reference
        pred['hic_kr_pl_scaled'] = pred['hic_kr'] * (powerlaw_estimate_reference / powerlaw_estimate)

    return(pred)

def add_hic_pseudocount(pred, args):

    #TO DO: Include Hi-C scale here - or deal with this constant issue. The pseudocount is too big
    powerlaw_fit = get_powerlaw_at_distance(pred['distance'].values, args.hic_gamma)
    powerlaw_fit_at_ref = get_powerlaw_at_distance(args.hic_pseudocount_distance, args.hic_gamma)
    
    pseudocount = np.amin(pd.DataFrame({'a' : powerlaw_fit, 'b' : powerlaw_fit_at_ref}), axis = 1)
    pred['hic_pseudocount'] = pseudocount
    pred['hic_kr_pl_scaled_adj'] = pred['hic_kr_pl_scaled'] + pseudocount

    return(pred)

def compute_score(enhancers, product_terms, prefix):

    scores = np.column_stack(product_terms).prod(axis = 1)

    enhancers[prefix + '.Score.Numerator'] = scores
    enhancers[prefix + '.Score'] = enhancers[prefix + '.Score.Numerator'] / enhancers.groupby('TargetGene')[prefix + '.Score.Numerator'].transform('sum')

    return(enhancers)

def annotate_predictions(pred):
    #Add is self promoter etc

    return(pred)

def make_gene_prediction_stats(pred, args):
    summ1 = pred.groupby(['chr','TargetGene','TargetGeneTSS']).agg({'TargetGeneIsExpressed' : lambda x: set(x).pop(), 'ABC.Score' : lambda x: all(np.isnan(x)) ,  'name' : 'count'})
    summ1.columns = ['geneIsExpressed', 'geneFailed','nEnhancersConsidered']

    summ2 = pred.loc[pred['class'] != 'promoter',:].groupby(['chr','TargetGene','TargetGeneTSS']).agg({args.score_column : lambda x: sum(x > args.threshold)})
    summ2.columns = ['nDistalEnhancersPredicted']
    summ1 = summ1.merge(summ2, left_index=True, right_index=True)

    summ1.to_csv(os.path.join(args.outdir, "GenePredictionStats.txt"), sep="\t", index=False)
