##################################################################################################################
# Generate data for PDM, AOV, String-Length, and Lomb-Scarge Periodogram Analysis
##################################################################################################################

# Null Hypothesis: p = \infty (constant function with normally distributed errors)

# Alternative Hypothesis: p \neq \infty (non-constant, periodic function with normally distributed errors)
    # options: 1) Well-detached eccentric system using two non-overlapping Gaussians
            #  2) Tight circular system with two overlapping Gaussians
            #  3) Tight circular system with an out-of-eclipse ellipsoidal variation
            #  4) Sinusoidal variation like RV measurements?
        # See figure 3 of Mowlavi 2022
        
###############################################################################
# Libraries and functions
###############################################################################

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os
from math import floor

plt.rcParams.update({'font.size': 15})

def gauss1(phi, d, mu, sigma):
    '''
    Objective: Generate the kernel of Gaussian density
    
    Parameters:
        phi: array of observed phase values
        d: ampitude of the Gaussian
        mu: mean of the Gaussian
        sigma: standard devation of the Gaussian

    Returns: Kernel of Gaussian density
    '''
    
    mag_gauss = d*np.exp(-((phi - mu)**2)/(2*sigma**2))
    
    return(mag_gauss)

def str_to_array(array_as_string):
    list_str = array_as_string.strip('][').split(',')
    list_float = [float(i) for i in list_str]
    array = np.asarray(list_float)
    return(array)

###############################################################################
# Paths
###############################################################################

#path = 'Documents/Astrostatistics/PD_review_paper/new_code/'
path = ''
raw_eb = path + 'raw/gdr3_vari_eclipsing_binary/'
raw_photo = path + 'raw/gdr3_epoch_photometry/'

inter_eb = path + 'inter/gdr3_vari_eclipsing_binary/'
inter_photo = path + 'inter/gdr3_epoch_photometry/'


###############################################################################
# Explore Mowlavi data to make sure you generate it correctly
###############################################################################

    
########################    
# Get the average number of observations for the data used in Mowlavi 2022
########################

# Get a list of the files produced by Mowlavi 2022
mowlavi_files_all = os.listdir(raw_eb)
    # remove the text file that was used to get the CSVs
mowlavi_files = [i for i in mowlavi_files_all if 'txt' not in i] 

# Get the source ID's that Mowlavi used
    # Get the number of sources
n_files = len(mowlavi_files)
n_sources = np.zeros(n_files)

for i in range(n_files):
    file_name = mowlavi_files[i]
    eb_gauss_df = pd.read_csv(raw_eb + file_name, header = 267)
    n_sources[i] = len(eb_gauss_df)
    
    # Save the file for later
    eb_gauss_df.to_pickle(inter_eb + file_name.replace('csv', 'pkl'))
    

    # Get the total number of sources
n_sources = [int(x) for x in n_sources]
np.save(inter_eb + 'n_sources', n_sources)
n_sources = np.load(inter_eb + 'n_sources.npy')
total_n_sources = np.sum(n_sources)


    # Get the data we need from Mowlavi
mowlavi_sol_ids = np.zeros(total_n_sources, dtype = float)
mowlavi_source_ids = np.zeros(total_n_sources, dtype = float)
mowlavi_b = np.zeros(total_n_sources, dtype = float)
mowlavi_g1_d = np.zeros(total_n_sources, dtype = float)
mowlavi_g1_s = np.zeros(total_n_sources, dtype = float)
mowlavi_g1_mu = np.zeros(total_n_sources, dtype = float)
mowlavi_g2_d = np.zeros(total_n_sources, dtype = float)
mowlavi_g2_s = np.zeros(total_n_sources, dtype = float)
mowlavi_g2_mu = np.zeros(total_n_sources, dtype = float)


j = 0
for i in range(n_files):
    file_name = mowlavi_files[i].replace('csv', 'pkl')
    eb_gauss_df = pd.read_pickle(inter_eb + file_name) 
    mowlavi_sol_ids[j:(j + n_sources[i])] = eb_gauss_df.solution_id.copy()
    mowlavi_source_ids[j:(j + n_sources[i])] = eb_gauss_df.source_id.copy()
    mowlavi_b[j:(j + n_sources[i])] = eb_gauss_df.geom_model_reference_level.copy()
    mowlavi_g1_d[j:(j + n_sources[i])] = eb_gauss_df.geom_model_gaussian1_depth.copy()
    mowlavi_g1_s[j:(j + n_sources[i])] = eb_gauss_df.geom_model_gaussian1_sigma.copy()
    mowlavi_g1_mu[j:(j + n_sources[i])] = eb_gauss_df.geom_model_gaussian1_phase.copy()
    mowlavi_g2_d[j:(j + n_sources[i])] = eb_gauss_df.geom_model_gaussian2_depth.copy()
    mowlavi_g2_s[j:(j + n_sources[i])] = eb_gauss_df.geom_model_gaussian2_sigma.copy()
    mowlavi_g2_mu[j:(j + n_sources[i])] = eb_gauss_df.geom_model_gaussian2_phase.copy()
    j = j + n_sources[i]
    
    
# Make pandas dataframes of Mowlavi data
mowlavi_df = pd.DataFrame({'solution_id': mowlavi_sol_ids,
                           'source_id': mowlavi_source_ids,
                           'b': mowlavi_b,
                           'd1': mowlavi_g1_d,
                           's1': mowlavi_g1_s,
                           'mu1': mowlavi_g1_mu,
                           'd2': mowlavi_g2_d,
                           's2': mowlavi_g2_s,
                           'mu2': mowlavi_g2_mu})

mowlavi_df.to_pickle(inter_eb + '/mowlavi.pkl')


# We're going to generate data under option 1 of the alternative hypothesis; what's the range of magnitudes that we'll encouter
eb_gauss_df = pd.read_pickle(inter_eb + '/VariEclipsingBinary_460242-460258.pkl')
i = 1
b = eb_gauss_df.geom_model_reference_level[i].copy()
d1 = eb_gauss_df.geom_model_gaussian1_depth[i].copy()
d2 = eb_gauss_df.geom_model_gaussian2_depth[i].copy()

# b - d1 = 18.3 and b - d2 = 19.2






# Get the source IDs, flux values, flux errors, and number of observations for Gaia objects
photo_files_all = os.listdir(raw_photo)
    # remove the text file that was used to get the CSVs
photo_files = [i for i in photo_files_all if 'txt' not in i]

    # Get the number of sources
n_files = len(photo_files)
n_sources = np.zeros(n_files)

for i in range(n_files):
    file_name = photo_files[i]
    photo_df = pd.read_csv(raw_photo + file_name, header = 365)
    n_sources[i] = len(photo_df)
    
    # Save the file for later
    photo_df.to_pickle(inter_photo + file_name.replace('csv', 'pkl'))
    
n_sources = [int(x) for x in n_sources]
np.save(inter_photo + 'n_sources', n_sources)
n_sources = np.load(inter_photo + 'n_sources.npy')
total_n_sources = np.sum(n_sources)

    # Get the source id's, magnitudes, fluxes, flux errors, and number of transits
        # for observations with observed magnitudes btw 18.3 and 19.5
photo_source_ids = np.zeros(total_n_sources, dtype = float)
photo_transit_times = [[]]*total_n_sources
photo_mag = [[]]*total_n_sources
photo_flux = [[]]*total_n_sources
photo_flux_err_ratio = [[]]*total_n_sources
n_transits = np.zeros(total_n_sources, dtype = float)

j = 0
for i in range(n_files):
    
    file_name = photo_files[i].replace('csv', 'pkl')
    photo_df = pd.read_pickle(inter_photo + file_name) 
       
    # Keep only fluxes and flux errors with observed magnitudes btw 18.3 and 19.5
    temp_gmag = photo_df.g_transit_mag.copy()
    
    temp_n = len(temp_gmag)
    
    for k in range(temp_n):
        # Get the indices of the magnitudes between 18.3 and 19.5
        idx = (18.3 < str_to_array(temp_gmag[k])) & (str_to_array(temp_gmag[k]) < 19.5)
        
        # Get the corresponding fluxes, flux errors, and G magnitude
        flux_keep = str_to_array(photo_df.g_transit_flux[k])[idx].copy()
        flux_err_keep = str_to_array(photo_df.g_transit_flux_error[k])[idx].copy()
        flux_error_ratio = flux_err_keep/flux_keep
        
        photo_mag[j+k] = str_to_array(photo_df.g_transit_mag[k])[idx].copy()
        photo_flux[j+k] = flux_keep.copy()
        photo_flux_err_ratio[j+k] = flux_error_ratio.copy()
        
    photo_source_ids[j:(j + temp_n)] = photo_df.source_id.copy()
    photo_transit_times[j:(j + temp_n)] = photo_df.g_transit_time.copy()
    n_transits[j:(j + temp_n)] = photo_df.n_transits.copy()
    j = j + temp_n

     
# Drop indices we didn't use
idx = [i for i,x in enumerate(photo_flux) if len(x) != 0]
photo_source_ids_keep = photo_source_ids[idx].copy()
photo_transit_times_keep = [x for x,y in zip(photo_transit_times, photo_flux) if len(y) != 0]
photo_mag_keep = [x for x in photo_flux if len(x) != 0]
photo_flux_keep = [x for x in photo_flux if len(x) != 0]
photo_flux_err_ratio_keep = [x for x in photo_flux_err_ratio if len(x) != 0]
n_transits_keep = n_transits[idx].copy()


# Make pandas dataframes of only the Gaia data that we need
gaia_df_need = pd.DataFrame({'source_id': photo_source_ids_keep,
                             'transit_times': photo_transit_times_keep,
                             'g_mag': photo_mag_keep,
                             'flux': photo_flux_keep,
                             'flux_err_ratio': photo_flux_err_ratio_keep,
                             'n_transits': n_transits_keep})
gaia_df_need.to_pickle(inter_photo + '/gaia_df_need.pkl')



# Merge the gaia and the mowlavi datasets
mg = pd.merge(gaia_df_need, mowlavi_df)
mg.to_pickle(path + 'inter/mg.pkl')


###############################################################################
# Set simulation parameters
###############################################################################


# Get the average number of transits in the combined Mowlavi and Gaia dataset
mg = pd.read_pickle(path + 'inter/mg.pkl')

mg_ntransits = floor(mg.n_transits.mean())
print(mg_ntransits)

# Set observation times using actual observation times from
    # an eclipsing binary with the same number of observations
    # as the mean number of observations of all eclipsing binaries

    # select a source
potential_sources = mg.loc[mg.n_transits == mg_ntransits, 'source_id'].copy()

sample_source = potential_sources.sample(1, random_state = 23897)
sample_source = float(sample_source.iloc[0])


obs_times = mg.loc[mg.source_id == sample_source, 'transit_times'].copy()
    # clean up the observation times to be an array
idx = obs_times.index[0]
obs_times_clean = str_to_array(obs_times[idx])

n_sims = 300


# Need to set a homoscedastic variance for the flux
        # 1) Calculate the average flux for observed magnitudes btw 18.3 and 19.5
        # 2) Calculate ratio of flux error to flux for observed magnitudes btw 18.3 and 19.5
        # 3) Calculate average error ratio
        # 4) Use (average error ratio)*(flux average) as the homoscedastic variance
flux_sum = 0
flux_err_ratio_sum = 0
flux_err_ratio_min = 1
flux_err_ratio_max = 0
flux_n = 0
for i in range(len(mg)):
    flux_sum += mg.flux[i].sum()
    flux_n += len(mg.flux[i])
    flux_err_ratio_sum += mg.flux_err_ratio[i].sum()
    min_err_ratio = mg.flux_err_ratio[i].min()
    max_err_ratio = mg.flux_err_ratio[i].max()
    
    flux_err_ratio_min = np.min([flux_err_ratio_min, min_err_ratio])
    flux_err_ratio_max = np.max([flux_err_ratio_max, max_err_ratio])
    

flux_avg = flux_sum/flux_n
flux_err_ratio_avg = flux_err_ratio_sum/flux_n
new_flux_error = flux_err_ratio_avg*flux_avg

np.save(path + 'inter/flux_avg', flux_avg)

np.random.seed(332489)



###############################################################################
# Generate data under the null hypothesis (no periodic component)
###############################################################################

####################################################
# Generate all of the obseved values that we'll use
####################################################
colnames = ['sim' + str(x) for x in range(n_sims)]
colnames.insert(0, 'obs_times')

null = pd.DataFrame(index = range(mg_ntransits), columns = colnames)

hetero = pd.DataFrame(index = range(mg_ntransits), columns = colnames)

for i in range(n_sims):
    
    new_flux = np.random.normal(loc = flux_avg, scale = new_flux_error, size = mg_ntransits)
    new_mag = 25.6874 + -2.5*np.log10(new_flux) # From Carrasco (2016) and Riello (2021)
    new_mag_sd = np.sqrt(((2.5*new_flux_error)/(new_flux*np.log(10)))**2 + 0.0028**2)
    null['sim' + str(i)] = new_mag/new_mag_sd
    hetero['sim' + str(i)] = new_mag


####################################################
# Just change the observation times
####################################################

# Evenly spaced time points
    # make a grid of evenly spaced points
step = (np.max(obs_times_clean) - np.min(obs_times_clean))/mg_ntransits
obs_times_even = np.arange(np.min(obs_times_clean), np.max(obs_times_clean), 
                           step)

null.obs_times = obs_times_even
null.to_pickle(path + 'inter/null_even.pkl')
null.to_csv(path + 'inter/null_even.csv', index = False)




# Evenly spaced time points; missing every 5th observation
step = (np.max(obs_times_clean) - np.min(obs_times_clean))/mg_ntransits
obs_times_extra = np.arange(np.min(obs_times_clean), np.max(obs_times_clean) + 11*step, 
                           step)
obs_times_even_drop = np.delete(obs_times_extra, np.arange(0, len(obs_times_extra), 5))

    # make sure we have the same number of observations
len(obs_times_even_drop)

    # Get a list of the dropped times
dropped_times = obs_times_extra[np.arange(0, len(obs_times_extra), 5)]
dropped_times_trunc = np.trunc(dropped_times)


null.obs_times = obs_times_even_drop
null.to_pickle(path + 'inter/null_even_missing.pkl')
null.to_csv(path + 'inter/null_even_missing.csv', index = False)




# Unevenly spaced time points
null.obs_times = obs_times_clean
null.to_pickle(path + 'inter/null_uneven.pkl')
null.to_csv(path + 'inter/null_uneven.csv', index = False)



# Heteroscedasticity from transformation
hetero.obs_times = obs_times_clean
hetero.to_pickle(path + 'inter/null_hetero_trans.pkl')
hetero.to_csv(path + 'inter/null_hetero_trans.csv', index = False)



# Unevenly spaced time points; missing obs
obs_times_drop1 = np.random.uniform(np.min(obs_times_extra), np.max(obs_times_extra), size = mg_ntransits + 10)
obs_times_drop2 = [x for x in obs_times_drop1 if np.trunc(x) not in dropped_times_trunc]
obs_times_drop = np.sort(np.random.choice(obs_times_drop2, mg_ntransits, replace = False))

null.obs_times = obs_times_drop
null.to_pickle(path + 'inter/null_uneven_missing.pkl')
null.to_csv(path + 'inter/null_uneven_missing.csv', index = False)


####################################################
# Generate data under the null hypothesis (no periodic component)
    # Unevenly spaced time points; heteroscedastic variance forced
####################################################

colnames = ['sim' + str(x) for x in range(n_sims)]
colnames.insert(0, 'obs_times')

null = pd.DataFrame(index = range(mg_ntransits), columns = colnames)
null.obs_times = obs_times_clean

hetero_error = np.arange(flux_err_ratio_min, flux_err_ratio_max + 0.1,
          step = (flux_err_ratio_max - flux_err_ratio_min)/(mg_ntransits-1))

for i in range(n_sims):
    
    new_flux = np.random.normal(loc = flux_avg, scale = hetero_error, size = mg_ntransits)
    new_mag = 25.6874 + -2.5*np.log10(new_flux) # From Carrasco (2016) and Riello (2021)
    null['sim' + str(i)] = new_mag

null.to_pickle(path + 'inter/null_hetero_sim.pkl')
null.to_csv(path + 'inter/null_hetero_sim.csv', index = False)


####################################################
# Generate data under option 1 of the alternative hypothesis
####################################################

# Arbitrarily set a period for the observation times
p = 150
phi = np.sort((obs_times_clean % p)/p)


# Extract data from one of the fitted models in Mowlavi 2022
eb_gauss_df = pd.read_pickle(inter_eb + '/VariEclipsingBinary_460242-460258.pkl')

# Find a fitted model without a cosine component and where the Gaussians do not overlap
    # (option 1 of the alternative hypothesis)
i = 1

b = eb_gauss_df.geom_model_reference_level[i].copy()
d1 = eb_gauss_df.geom_model_gaussian1_depth[i].copy()
sigma1 = eb_gauss_df.geom_model_gaussian1_sigma[i].copy()
mu1_orig = eb_gauss_df.geom_model_gaussian1_phase[i].copy()

mu1 = mu1_orig*p + p*11
sigma1 = sigma1*p

d2 = eb_gauss_df.geom_model_gaussian2_depth[i].copy()
sigma2 = eb_gauss_df.geom_model_gaussian2_sigma[i].copy()
mu2_orig = eb_gauss_df.geom_model_gaussian2_phase[i].copy()

mu2 = mu2_orig*p + p*11
sigma2 = sigma2*p

x = np.linspace(np.min(obs_times_clean), np.max(obs_times_clean), num = 1000)

mag_gauss = b + gauss1(x, d1, mu1, sigma1) + gauss1(x, d2, mu2, sigma2) +\
    gauss1(x, d1, mu1 + p, sigma1) + gauss1(x, d2, mu2 + p, sigma2) +\
        gauss1(x, d1, mu1 + 2*p, sigma1) + gauss1(x, d2, mu2 + 2*p, sigma2) +\
           gauss1(x, d1, mu1 + 3*p, sigma1) + gauss1(x, d2, mu2 + 3*p, sigma2)  +\
               gauss1(x, d1, mu1 + 4*p, sigma1) + gauss1(x, d2, mu2 + 4*p, sigma2) +\
                   gauss1(x, d1, mu1 + 5*p, sigma1) + gauss1(x, d2, mu2 + 5*p, sigma2) +\
                       gauss1(x, d1, mu1 + 6*p, sigma1) + gauss1(x, d2, mu2 + 6*p, sigma2) 
                       
flux_mean = 10**((mag_gauss - 25.6874)/(-2.5))
new_mag_sd = np.sqrt(((2.5*new_flux_error)/(flux_mean*np.log(10)))**2 + 0.0028**2)

plt.scatter(x, mag_gauss/new_mag_sd, c = 'black', s = 10)
plt.ylabel('Scaled G-Magnitude')
plt.xlabel('Barycentric Julian Date from 1/1/2010')
plt.savefig(path + 'output/sims_plots/h1_function.pdf', bbox_inches = 'tight',
            dpi = 600)
plt.close()



x = obs_times_clean.copy()

mag_gauss = b + gauss1(x, d1, mu1, sigma1) + gauss1(x, d2, mu2, sigma2) +\
    gauss1(x, d1, mu1 + p, sigma1) + gauss1(x, d2, mu2 + p, sigma2) +\
        gauss1(x, d1, mu1 + 2*p, sigma1) + gauss1(x, d2, mu2 + 2*p, sigma2) +\
           gauss1(x, d1, mu1 + 3*p, sigma1) + gauss1(x, d2, mu2 + 3*p, sigma2)  +\
               gauss1(x, d1, mu1 + 4*p, sigma1) + gauss1(x, d2, mu2 + 4*p, sigma2) +\
                   gauss1(x, d1, mu1 + 5*p, sigma1) + gauss1(x, d2, mu2 + 5*p, sigma2) +\
                       gauss1(x, d1, mu1 + 6*p, sigma1) + gauss1(x, d2, mu2 + 6*p, sigma2) 
        
plt.scatter(x, -mag_gauss)

np.save(path + 'inter/mag_gauss', mag_gauss)

# convert f_guass to flux value
flux_mean = 10**((mag_gauss - 25.6874)/(-2.5))



# Generate data under the alternative hypothesis for option 1;
    # no error associated with Gmag; simulate it
    # https://gea.esac.esa.int/archive/documentation/GDR3/Data_analysis/chap_cu4sso/sec_cu4sso_processingsteps/ssec_cu4sso_photometryprocessing.html

h1 = pd.DataFrame(index = range(mg_ntransits), columns = colnames)
h1.obs_times = obs_times_clean

for i in range(n_sims):
    
    # generate observed flux
    new_flux = np.random.normal(loc = flux_mean, scale = new_flux_error, size = mg_ntransits)
    
    # convert back to magnitude
    h1_mag = 25.6874 + -2.5*np.log10(new_flux)
    new_mag_sd = np.sqrt(((2.5*new_flux_error)/(new_flux*np.log(10)))**2 + 0.0028**2)
    
    h1['sim' + str(i)] = h1_mag/new_mag_sd


h1.to_pickle(path + 'inter/h1.pkl')
h1.to_csv(path + 'inter/h1.csv', index = False)

























