#########################################################
# Libraries
#########################################################

import numpy as np
import pandas as pd

from datetime import datetime
from multiprocessing import Pool

#########################################################
# Functions
#########################################################

def gen_W(df):
    '''
    Objective: Generate the design matrix W for calculating the projection matrix P_w
                and the vector of the expected values of the data (eta)
    
    Parameters:
        df: dataframe of the observed values that have the assigned bins

    Returns: Design matrix W
    '''
    
    W = np.zeros((len(df), len(df.bin.unique())))
    i = 0
    for j in df.bin.unique():
        W[:, i] = (df.bin == j)
        i += 1
        
    return(W)

def pdm_aov_fun2(data, eta, p, n_iter, iter_cols, bin_width):

    df = data.copy()
    df['phase_scale'] = (df.obs_times%p)/p
    df['bin'] = np.floor(df.phase_scale/bin_width)
    
    N = len(df)
    r = len(df.groupby('bin').count())
    
    try:
    
        W = gen_W(df)
        p1 = (1/N)*np.outer(np.ones(N), np.ones(N))
        pw = np.dot(W, np.dot(np.linalg.pinv(np.dot(W.T, W)), W.T))
        phi = np.dot(eta, np.dot(pw - p1, eta))
        
        
        bin_counts = df.groupby('bin', as_index = False).obs_times.count()
        bin_counts = bin_counts.rename(columns = {'obs_times': 'bin_count'})
        
        all_group_means = df.groupby('bin', as_index = False).mean().copy()
        all_group_means = all_group_means.drop(columns=['obs_times', 'phase_scale']).copy()
        
        colnames = ['bmean' + str(i) for i in range(n_iter)]
        colnames.insert(0, 'bin')
        all_group_means.columns = colnames
        
        bins = bin_counts.merge(all_group_means).copy()
        
        temp = df.merge(all_group_means)
        
        all_means = df.mean(axis = 0)
        
        bmean_cols = ['bmean' + str(i) for i in range(n_iter)]
        
        s0 = [np.linalg.norm(df[x] - all_means[x])**2/(N - 1) for x in iter_cols]
        s2 = [np.linalg.norm(temp[x] - temp[y])**2/(N - r) for x,y in zip(iter_cols, bmean_cols)]
        s1 = [sum(bins.bin_count*(bins[x] - all_means[y])**2)/(r-1) for x,y in zip(bmean_cols, iter_cols)]
        
        pdm = [(N - r)/(N-1)*(x/y) for x,y in zip(s2, s0)] # scaled pdm
        aov = [x/y for x,y in zip(s1, s2)]
        
    except ZeroDivisionError:
        pdm = np.empty(n_iter)*np.nan
        aov = np.empty(n_iter)*np.nan
        phi = np.nan
    
    return(np.array(pdm), np.array(aov), phi)

def calc_sl(ordered_phase, ordered_obs):
    
    lagged_phase = np.roll(ordered_phase, 1)
    lagged_phase[0] = lagged_phase[0] - 1
    
    lagged_x = np.roll(ordered_obs, 1)

    sl = (np.sqrt((ordered_phase - lagged_phase)**2 + (ordered_obs - lagged_x)**2)).sum()
    
    return sl

#########################################################
# Calculate the AOV and PDM statistics from the data
#########################################################

def aov_pdm_wrapper(p):
    
    niter = 300
    bw = 0.1
    colnames = data.columns.tolist()
    colnames.remove('obs_times')
    
    #path = 'Documents/Astrostatistics/PD_review_paper/new_code/'
    path = ''
    
    flux_avg = float(np.load(path + 'inter/flux_avg.npy'))
    eta = flux_avg*np.ones(len(data))
    
    pdm, aov, phi = pdm_aov_fun2(data, eta, p, niter, colnames, bw)
    
    return(np.array(pdm), np.array(aov), phi)
    
    
def calc_sl_wrapper(p):
    
    df = data.copy()
    df['phase_scale'] = (df.obs_times%p)/p
    
    # calculate the string-length metric
    df_ordered = df.sort_values(by = ['phase_scale']).copy()
    sims_ordered = df_ordered.drop(['obs_times', 'phase_scale'], axis = 1).copy()
    
    sl_all_sims = np.array(sims_ordered.apply(lambda x: calc_sl(df_ordered.phase_scale, x), axis = 0))
    
    return(sl_all_sims)
    




if __name__ == '__main__':
    
    # Import libraries
    import argparse

    # Allow arguments from the command-line
    parser = argparse.ArgumentParser(description='Which dataset to use?')

    parser.add_argument('-dname', '--dname', type = str, help='Data name')
    args = parser.parse_args()
    
    data = pd.read_pickle('inter/' + args.dname + '.pkl')
    #data = pd.read_pickle('Documents/Astrostatistics/PD_review_paper/new_code/inter/null_even_missing.pkl')
    
    f_min = 0.005
    f_max = 15
    f_seq = np.arange(f_min, f_max + 10**(-5), 10**(-5))
    p_seq = 1/f_seq
    
    # append the true period for plotting purposes
    p_seq = np.append(p_seq, 150)
    
    # make it manageable for the purposes of just getting graphs and stuff together
    #p_seq = np.arange(20, 40 + 10**(-2), 10**(-2))
    
    
    periods = p_seq.tolist()
    
    start = datetime.now()
    
    with Pool() as pool:
        result_pdm_aov = pool.map(aov_pdm_wrapper, periods)
        
        
    # scale the columns for the string-length calculation
    sims = data.drop(['obs_times'], axis = 1).copy()
    sims_scaled = sims.apply(lambda x: (x - x.min())/(x.max() - x.min()))
    data1 = sims_scaled.copy()
    data1.insert(0, 'obs_times', data.obs_times.copy())
    data = data1.copy()
        
    with Pool() as pool:
        result_sl = pool.map(calc_sl_wrapper, periods)
        
        
    # compile the results
    pdm_all = np.zeros((len(p_seq), len(data.columns)))
    pdm_all[:,0] = p_seq
    
    aov_all = np.zeros((len(p_seq), len(data.columns)))
    aov_all[:,0] = p_seq
    
    phi_all = np.zeros(len(p_seq))
    
    
    for i in range(len(p_seq)):
        pdm_all[i,1:] = result_pdm_aov[i][0] 
        aov_all[i,1:] = result_pdm_aov[i][1]
        phi_all[i] = result_pdm_aov[i][2]

    
    # make results into dataframes
    pdm_all_df = pd.DataFrame(pdm_all)
    aov_all_df = pd.DataFrame(aov_all)
    phi_all_df = pd.DataFrame(phi_all)
    
    sl_df = pd.DataFrame(result_sl)
    sl_df.insert(0, 'p', p_seq)
        
    # save results
    pdm_all_df.to_pickle('output/sims_data/pdm_all' + args.dname + '.pkl')
    aov_all_df.to_pickle('output/sims_data/aov_all' + args.dname + '.pkl')
    phi_all_df.to_pickle('output/sims_data/ncp_all' + args.dname + '.pkl')
    sl_df.to_pickle('output/sims_data/sl_all' + args.dname + '.pkl')
    
    
    end = datetime.now()
    
    rt = end - start
    print(rt)
















































