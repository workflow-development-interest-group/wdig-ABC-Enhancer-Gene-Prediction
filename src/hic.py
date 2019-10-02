import numpy as np
import scipy.sparse as ssp
import pandas as pd
import time

def load_hic(hic_file, hic_type, hic_resolution, tss_hic_contribution, window, min_window, gamma, interpolate_nan=True, apply_diagonal_bin_correction=True):
    print("Loading HiC")

    if hic_type == 'juicebox':
        HiC_sparse_mat = hic_to_sparse(hic_file, hic_resolution)
        HiC = process_hic(hic_mat = HiC_sparse_mat, 
                            resolution = hic_resolution, 
                            tss_hic_contribution = tss_hic_contribution, 
                            window = window, 
                            min_window = min_window, 
                            gamma = gamma,
                            interpolate_nan = interpolate_nan)
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

def process_hic(hic_mat, resolution, tss_hic_contribution, window, min_window=0, hic_is_doubly_stochastic=False, apply_diagonal_bin_correction=True, interpolate_nan=True, gamma=None):
    #Make doubly stochastic.
    #Juicer produces a matrix with constant row/column sums. But sum is not 1 and is variable across chromosomes
    t = time.time()

    if not hic_is_doubly_stochastic:
        sums = hic_mat.sum(axis = 0)
        assert(np.max(sums[sums > 0])/np.min(sums[sums > 0]) < 1.001)
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

    #Fill NaN
    #NaN in the KR normalized matrix are not zeros. They are entries where the KR algorithm did not converge
    #So need to fill these. Use powerlaw. 
    #Not ideal obviously but the scipy interpolation algos are either very slow or don't work since the nan structure implies that not all nans are interpolated
    if interpolate_nan:
        nan_loc = np.isnan(hic_df['hic_kr'])
        hic_df.loc[nan_loc,'hic_kr'] = get_powerlaw_at_distance(abs(hic_df.loc[nan_loc,'bin1'] - hic_df.loc[nan_loc,'bin2']), gamma)

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

    # drop NaNs from hic
    # If loading KR Norm, then don't want to remove NaN
    #HiC = HiC.loc[~np.isnan(HiC['hic_kr']),:]
    #print("HiC has {} rows after dropping NaNs".format(HiC.shape[0]))

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

def get_powerlaw_at_distance(distances, gamma, scale=None):
    assert(gamma > 0)
    log_dists = np.log(distances + 1)

    #Determine scale parameter
    #A powerlaw distribution has two parameters: the exponent and the minimum domain value 
    #In our case, the minimum domain value is always constant (equal to 1 HiC bin) so there should only be 1 parameter
    #The current fitting approach does a linear regression in log-log space which produces both a slope (gamma) and a intercept (scale)
    #Empirically there is a linear relationship between these parameters (which makes sense since we expect only a single parameter distribution)
    #It should be possible to analytically solve for scale using gamma. But this doesn't quite work since the hic data does not actually follow a power-law
    #So could pass in the scale parameter explicity here. Or just kludge it as I'm doing now
    #TO DO: Eventually the pseudocount should be replaced with a more appropriate smoothing procedure.

    #4.80 and 11.63 come from a linear regression of scale on gamma across 20 hic cell types at 5kb resolution. Do the params change across resolutions?
    if scale is None:
        scale = -4.80 + 11.63 * gamma

    powerlaw_contact = np.exp(scale + -1*gamma * log_dists)

    return(powerlaw_contact)
