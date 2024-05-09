# -*- coding: utf-8 -*-
"""
Created on Sat Nov 13 14:02:18 2021

@author: ngierty

# Purpose: Validate the Lomb-Scargle Periodogram extremum bounds

# References: https://docs.astropy.org/en/stable/timeseries/lombscargle.html#id35
"""

#######################################################
# Libraries
#######################################################

import numpy as np
import random
from math import floor, sqrt, pi
import matplotlib.pyplot as plt
from astropy.timeseries import LombScargle

plt.rcParams.update({'font.size': 15})

#######################################################
# Parameters
# T: Amount of time spanned; 
# n: Number of observations
# p
# n0
#######################################################

np.random.seed(49853)

# Parameters
n = 100
T = 24*5
n0 = 5
p = 0.01
niter = 10000

np.random.RandomState(348)
random.seed(4567)

path = 'Documents/Astrostatistics/PD_review_paper/new_code/output/sims_plots/'

#######################################################
# Simulations
#######################################################


# Generate the random observation times
# Following VanderPlas 2018 to get a precise Nyquist frequency for the trial frequencies
ni = np.sort(random.sample(range(T), n))

ts = ni * p

f_max = 24*10

frequency_grid = np.linspace(1/(n0 * T), f_max, n0 * T)

# Save the Lomb Scargle Periodogram values for each frequency
direct_ls_power = np.zeros((len(frequency_grid), niter))


for i in range(niter):
    
    np.random.seed(i)

    # Generate the X values
    xs = np.random.normal(size = n)
    
    # Calculate the Lomb Scargle periodogram, using the default set of frequencies
    power = LombScargle(ts, xs, fit_mean=False, center_data=True,
                        normalization='psd').power(frequency_grid)
    
    # Save the first set of periodogram values
    direct_ls_power[:, i] = power


# Check that we didn't accidentally generate the exact same dataset for every dataset
plt.plot(frequency_grid, direct_ls_power[:, 0])
plt.plot(frequency_grid, direct_ls_power[:, 1])
plt.show()


# For the last simulated dataset, bootstrap niter times 
    # to get bootstrap estimates of the FAP
boots_direct_ls_power = np.zeros((len(frequency_grid), niter))

for i in range(niter):
    
    # bootstrap the sample
    np.random.seed(i)
    resample = np.random.randint(0, n, n)
    xs_boot = xs[resample].copy() 
    
    # calculate the Lomb-Scargle power of the bootstrapped sample
    power = LombScargle(ts, xs_boot, fit_mean=False, center_data=True,
                        normalization='psd').power(frequency_grid)
    
    boots_direct_ls_power[:, i] = power



# Calculate the maximum values for each of the periodograms and save
direct_ls_power_max = np.zeros(niter)
boots_ls_power_max = np.zeros(niter)

for i in range(niter):
    direct_ls_power_max[i] = direct_ls_power[:, i].max()
    boots_ls_power_max[i] = boots_direct_ls_power[:, i].max()
    
    

plt.hist(direct_ls_power_max, bins = floor(sqrt(niter)),
         density = True, color = 'lightgray')
plt.show()

plt.hist(boots_ls_power_max, bins = floor(sqrt(niter)),
         density = True, color = 'lightgray')
plt.show()

# For a grid of u values calculate the probability that the max of Lomb-Scargle
# Periodogram is bigger than u and calculate the Baluev (2008) bound
u_grid_n = 100

u_min = max(direct_ls_power_max.min(), boots_ls_power_max.min())
u_max = min(direct_ls_power_max.max(), boots_ls_power_max.max())

u_grid = np.linspace(u_min, u_max, u_grid_n)
direct = np.zeros(u_grid_n)
boots_direct = np.zeros(u_grid_n)


for i in range(u_grid_n):

    u = u_grid[i]
    direct[i] = sum(direct_ls_power_max >= u)/niter
    boots_direct[i] = sum(boots_ls_power_max >= u)/niter
    
    
    
var_t = ts.var()
Teff = sqrt(4*pi*var_t)
W = f_max*Teff
tau_davies = W*np.exp(-u_grid)*np.sqrt(u_grid)
davies = np.exp(-u_grid) + tau_davies
baluev = 1 - (1 - np.exp(-u_grid))*np.exp(-tau_davies)


plt.plot(u_grid, direct, c = 'black', label = 'Direct')
plt.plot(u_grid, boots_direct, c = 'red', linestyle = ':', label = 'Bootstrap')
#plt.plot(u_grid, davies, c = 'green', linestyle = '-.', label = 'Davies Bound')
plt.plot(u_grid, baluev, c = 'blue', linestyle = '--', label = 'Baluev Bound')
plt.xlabel('Maximum Value of Lomb-Scargle Periodogram')
plt.ylabel('False Alarm Probability')
plt.ylim(0,1)
plt.legend()
#plt.show()
plt.savefig(path + 'fap_bound.png', bbox_inches = 'tight', dpi = 600)
plt.close()


###############
# Change the units on the y-axis to match Baluev and Vanderplas
###############

log10_direct = np.log10(direct)
log10_boot_direct = np.log10(boots_direct)
log10_davies = np.log10(davies)
log10_baleuv = np.log10(baluev)

# Get reasonable minimum for plotting purposes
all_probs = np.concatenate((direct, boots_direct, davies, baluev),axis=0)
y_min = floor(np.log10(all_probs.min()))

# Set yticks and labels
yticks = np.arange(y_min, 1, 1.0)
ylabels = 10**yticks


plt.plot(u_grid, log10_direct, c = 'black', label = 'Direct')
plt.plot(u_grid, log10_boot_direct, c = 'red', linestyle = ':', label = 'Bootstrap')
#plt.plot(u_grid, log10_davies, c = 'green', linestyle = '-.', label = 'Davies Bound')
plt.plot(u_grid, log10_baleuv, c = 'blue', linestyle = '--', label = 'Baluev Bound')
plt.xlabel('Largest LS Power')
plt.ylabel('False Alarm Probability')
plt.yticks(ticks = yticks, labels = ylabels)
plt.ylim(y_min,0.5)
plt.legend()
#plt.show()
plt.savefig(path + 'fap_bound_log.pdf', bbox_inches = 'tight', dpi = 600)
plt.close()

















