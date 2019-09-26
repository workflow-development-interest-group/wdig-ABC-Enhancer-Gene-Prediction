#from weakref import WeakValueDictionary

import numpy as np
import scipy.sparse as ssp
import pandas as pd
import time

def load_hic(hic_file, hic_type, hic_resolution, tss_hic_contribution, window, min_window):
    print("Loading HiC")

    if hic_type == 'juicebox':
        HiC_sparse_mat = hic_to_sparse(hic_file, hic_resolution)
        HiC = process_hic(HiC_sparse_mat, hic_resolution, tss_hic_contribution, window, min_window)
        #HiC = juicebox_to_bedpe(HiC, chromosome, args)
    elif hic_type == 'bedpe':
        HiC = pd.read_csv(hic_file, sep="\t")

    return(HiC)

# def juicebox_to_bedpe(hic, chromosome, resolution):
#     hic['chr'] = chromosome
#     hic['x1'] = hic['bin1'] * resolution
#     hic['x2'] = (hic['bin1'] + 1) * resolution
#     hic['y1'] = hic['bin2'] * resolution
#     hic['y2'] = (hic['bin2'] + 1) * resolution

#     return(hic)

def process_hic(hic_mat, resolution, tss_hic_contribution, window, min_window=0, hic_is_doubly_stochastic=False, apply_diagonal_bin_correction=True):
    #Make doubly stochastic.
    #Juicer produces a matrix with constant row/column sums. But sum is not 1 and is variable across chromosomes
    t = time.time()

    if not hic_is_doubly_stochastic:
        sums = hic_mat.sum(axis = 0)
        assert(np.max(sums)/np.min(sums[sums > 0]) < 1.001)
        mean_sum = np.mean(sums[sums > 0])
        print('HiC Matrix has row sums of {}, making doubly stochastic...'.format(mean_sum))
        hic_mat = hic_mat.multiply(1/mean_sum)

    #Slow version. Its a constant scalar so don't need to to the matrix multiplication
    # kr_vec = np.repeat(np.sqrt(mean_sum), sums.shape[1])
    # norm_mat = ssp.dia_matrix((1.0 / kr_vec, [0]), (sums.shape[1], sums.shape[1]))
    # hic_mat = norm_mat * hic_mat * norm_mat

    #Adjust diagonal of matrix based on neighboring bins
    # hic_mat = hic_mat.tolil(copy=False)
    # for ii in range(1, sums.shape[1] - 1):
    #     hic_mat[ii,ii] = max(hic_mat[ii,ii-1], hic_mat[ii,ii+1]) * args.tss_hic_contribution / 100

    #Adjust diagonal of matrix based on neighboring bins
    #First and last rows need to be treated differently
    if apply_diagonal_bin_correction:
        last_idx = hic_mat.shape[0] - 1
        nonzero_diag = hic_mat.nonzero()[0][hic_mat.nonzero()[0] == hic_mat.nonzero()[1]]
        nonzero_diag = list(set(nonzero_diag) - set(np.array([last_idx])) - set(np.array([0])))

        for ii in nonzero_diag:
            hic_mat[ii,ii] = max(hic_mat[ii,ii-1], hic_mat[ii,ii+1]) * tss_hic_contribution / 100

        if hic_mat[0,0] != 0:
            hic_mat[0, 0] = hic_mat[0,1] * tss_hic_contribution / 100

        if hic_mat[last_idx, last_idx] != 0:
            hic_mat[last_idx, last_idx] = hic_mat[last_idx, last_idx - 1] * tss_hic_contribution / 100

    #Remove lower triangle
    hic_mat = ssp.triu(hic_mat)

    #Turn into dataframe
    hic_mat = hic_mat.tocoo(copy=False)
    hic_df = pd.DataFrame({'bin1': hic_mat.row, 'bin2': hic_mat.col, 'hic_kr': hic_mat.data})

    #Prune to window
    hic_df = hic_df.loc[np.logical_and(abs(hic_df['bin1'] - hic_df['bin2']) <= window/resolution, abs(hic_df['bin1'] - hic_df['bin2']) >= min_window/resolution)]
    print("HiC has {} rows after windowing between {} and {}".format(hic_df.shape[0], min_window, window))

    print('process.hic: Elapsed time: {}'.format(time.time() - t))


    return(hic_df)

def hic_to_sparse(filename, resolution, hic_is_doubly_stochastic=False):
    t = time.time()
    HiC = pd.read_table(filename, names=["bin1", "bin2", "hic_kr"],
                        header=None, engine='c', memory_map=True)

    # verify our assumptions
    assert np.all(HiC.bin1 <= HiC.bin2)

    # find largest entry
    max_pos = max(HiC.bin1.max(), HiC.bin2.max())
    hic_size = max_pos // resolution + 1

    # TO DO: Interpolate NAN
    # drop NaNs from hic
    HiC = HiC.loc[~np.isnan(HiC['hic_kr']),:]
    print("HiC has {} rows after dropping NaNs".format(HiC.shape[0]))

    # convert to sparse matrix in CSR (compressed sparse row) format, chopping
    # down to HiC bin size.  note that conversion to scipy sparse matrices
    # accumulates repeated indices, so this will do the right thing.
    row = np.floor(HiC.bin1.values / resolution).astype(int)
    col = np.floor(HiC.bin2.values / resolution).astype(int)
    dat = HiC.hic_kr.values

    #JN: Need both triangles in order to compute row/column sums to make double stochastic.
    #If juicebox is upgraded to return DS matrices, then can remove one triangle
    #TO DO: Remove one triangle when juicebox is updated.
    # we want a symmetric matrix.  Easiest to do that during creation, but have to be careful of diagonal
    if not hic_is_doubly_stochastic:
        mask = (row != col)  # off-diagonal
        row2 = col[mask]  # note the row/col swap
        col2 = row[mask]
        dat2 = dat[mask]

        # concat and create
        row = np.hstack((row, row2))
        col = np.hstack((col, col2))
        dat = np.hstack((dat, dat2))

    print('hic.to.sparse: Elapsed time: {}'.format(time.time() - t))

    return ssp.csr_matrix((dat, (row, col)), (hic_size, hic_size))
