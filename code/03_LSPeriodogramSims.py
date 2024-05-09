#######################################################
# Libraries
#######################################################

import numpy as np
import pandas as pd
from astropy.timeseries import LombScargle

#######################################################
# Functions
#######################################################

    

def calc_ls_power(data, iter_cols, frequency_grid, path, name):
    
    niter = len(iter_cols)


    power_df = pd.DataFrame(np.empty((len(f_seq), niter)))
    
    power_df_scaled = pd.DataFrame(np.empty((len(f_seq), niter)))
    
    
    for i in range(1,niter):
        
        var_data = data.iloc[:, i].var()
        
        one_sim_power = LombScargle(data.obs_times, data[iter_cols[i]],
                    fit_mean=False, center_data=True,
                    normalization='psd').power(frequency_grid)
        
        power_df.iloc[:, i-1] = one_sim_power
        
        power_df_scaled.iloc[:, i-1] = 2*one_sim_power/var_data
        
        
    # Add in a column of the trial periods (using periods for plotting purposes)
    power_df.insert(0, 'p', 1/f_seq)
    power_df_scaled.insert(0, 'p', 1/f_seq)
        
    
    # Save the dataframes
    power_df.to_pickle(path + 'ls_power_all_' + name + '.pkl')
    power_df_scaled.to_pickle(path + 'ls_power_scaled_all_' + name + '.pkl')


#########################################################
# Read in the data
#########################################################

#path = 'Documents/Astrostatistics/PD_review_paper/new_code/' 
path = ''

# Read in the simulated data from Gen_Data.py
null_even = pd.read_pickle(path + 'inter/null_even.pkl')
null_uneven = pd.read_pickle(path + 'inter/null_uneven.pkl')
null_hetero_trans = pd.read_pickle(path + 'inter/null_hetero_trans.pkl')
null_hetero_sim = pd.read_pickle(path + 'inter/null_hetero_sim.pkl')
h1 = pd.read_pickle(path + 'inter/h1.pkl')

#########################################################
# Calculate the periodogram for each dataset
#########################################################


f_min = 0.005
f_max = 15
f_seq = np.arange(f_min, f_max + 10**(-5), 10**(-5))
    # append the true period for plotting purposes
f_seq = np.append(f_seq, 1/150)


# Calculate the Lomb Scargle periodogram for each sim
iter_cols = null_even.columns.copy()
iter_cols = iter_cols.drop('obs_times')

calc_ls_power(null_even, iter_cols, f_seq, path + 'output/sims_data/', 'null_even')
calc_ls_power(null_uneven, iter_cols, f_seq, path + 'output/sims_data/', 'null_uneven')
calc_ls_power(null_hetero_trans, iter_cols, f_seq, path + 'output/sims_data/', 'null_hetero_trans')
calc_ls_power(null_hetero_sim, iter_cols, f_seq, path + 'output/sims_data/', 'null_hetero_sim')
calc_ls_power(h1, iter_cols, f_seq, path + 'output/sims_data/', 'h1')































