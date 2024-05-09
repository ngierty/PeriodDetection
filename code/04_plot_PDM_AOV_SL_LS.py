###############################################################################
# Libraries
###############################################################################

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from math import floor, sqrt
import scipy.stats as ss

plt.rcParams['agg.path.chunksize'] = 100000

###############################################################################
# Functions
###############################################################################

# scatter plot of the observed data
def plot_obs_data(data, ylab, path, name):
    
    plt.scatter(data.obs_times, data.sim0, c = 'black', s = 10)
    plt.ylabel(ylab)
    plt.xlabel('Barycentric Julian Date from 1/1/2010')
    plt.savefig(path + 'obs_data_' + name + '.pdf', bbox_inches = 'tight',
                dpi = 600)
    plt.close()


# scatter plot of the observed data in-phase
def plot_phase(data, p, n_iter, iter_cols, bin_width, ax = None):

    df = data.copy()
    df['phase_scale'] = (df.obs_times%p)/p
    df['bin'] = np.floor(df.phase_scale/bin_width)
    
    
    bin_counts = df.groupby('bin', as_index = False).obs_times.count()
    bin_counts = bin_counts.rename(columns = {'obs_times': 'bin_count'})
    
    all_group_means = df.groupby('bin', as_index = False).mean().copy()
    all_group_means = all_group_means.drop(columns=['obs_times', 'phase_scale']).copy()
    
    colnames = ['bmean' + str(i) for i in range(n_iter)]
    colnames.insert(0, 'bin')
    all_group_means.columns = colnames
    
    bins = bin_counts.merge(all_group_means).copy()
    
    # Plot the observations of each bin and the group mean for one of the iterations
    bin_labs = df.bin.unique()
    bin_labs.sort()
    
    # Plot the phase-folded data
    for i in bin_labs:
        ax.scatter(df.phase_scale[df.bin == i], df[iter_cols[0]][df.bin == i],
                    c = 'black', alpha = (i+1)*0.07)
            
        ax.plot([bin_width*i, bin_width*(i+1)], [bins.bmean0[bins.bin == i], bins.bmean0[bins.bin == i]],
                 c = 'black', linewidth = 3)

    
    ax.set_xlabel('Phase-folded Time')
    
    return ax
    
    
# histogram of the AOV statistics
def plot_aov(aov, data, p, bin_width, max_aov, ax = None):
    
    df = data.copy()
    df['phase_scale'] = (df.obs_times%p)/p
    df['bin'] = np.floor(df.phase_scale/bin_width)
    
    N = len(df)
    r = len(df.groupby('bin').count())
    
    n_iter = len(aov)

    # Plot histograms of the AOV statistics
    ax.hist(aov, bins = floor(sqrt(n_iter)),
             density = True, color = 'lightgray')
    df1 = r - 1
    df2 = N - r
    x = np.linspace(0, max_aov, 1000)
    pdf = ss.ncf.pdf(x, df1, df2, 0)
    ax.plot(x, pdf, c = 'black', ls = '-')
    
    ax.set_xlabel('AOV Statistic')
    
    return ax


# histogram of the PDM statistics
def plot_pdm(pdm, data, p, bin_width, min_pdm, ax = None):
    
    df = data.copy()
    df['phase_scale'] = (df.obs_times%p)/p
    df['bin'] = np.floor(df.phase_scale/bin_width)
    
    N = len(df)
    r = len(df.groupby('bin').count())
    
    n_iter = len(pdm)
    
    # Plot histograms of the PDM statistics
    ax.hist(pdm, bins = floor(sqrt(n_iter)),
             density = True, color = 'lightgray')
    x = np.linspace(min_pdm, 1, 1000)
    const = (N - r)/(r - 1)
    pdf = const*(x**(-2))*ss.ncf.pdf(const*((1/x)-1), r - 1, N - r, 0)
    ax.plot(x, pdf, c = 'black', ls = '-')
    ax.set_xlabel('Scaled PDM Statistic')
    
    return ax



# histogram of the SL statistics
def plot_sl(sl, xlab, ax = None):
    
    n_iter = len(sl)
    
    ax.hist(sl, bins = floor(sqrt(n_iter)),
             density = True, color = 'lightgray')
    ax.set_xlabel('SL Statistic')

    return ax


# Histograms of the Lomb-Scargle Power statistics
def plot_lsp(lsp, max_lsp, ax = None):

    n_iter = len(lsp)    

    ax.hist(lsp, bins = floor(sqrt(n_iter)),
             density = True, color = 'lightgray')
    x = np.linspace(0, max_lsp, 1000)
    pdf = ss.chi2.pdf(x, 2)
    ax.plot(x, pdf, c = 'black', ls = '-')
    ax.set_xlabel('Scaled Lomb-Scargle Power')
    
    return ax
    
    
    


# final figure of the scatter plot of the phase-folded observations,
    # AOV, PDM, and SL statistics for a set of trial periods
def plot_stat_hists(scatter_data, aov_data, pdm_data, sl_data, lsp_data,
                    trial_periods, iter_cols, bin_width, ylab,
                    output_path, name):
    
    
    n_iter = scatter_data.shape[1] - 1


    fig = plt.figure()
    fig.set_figheight(30)
    fig.set_figwidth(22)
    widths = [4, 4, 4]
    heights = [1, 1, 1, 1, 1]
    gs = fig.add_gridspec(nrows = 5, ncols = 3, wspace=0.1, hspace = 0.3,
                          width_ratios=widths, height_ratios=heights)
    axs = gs.subplots(sharex='row', sharey='row')
    
    # Get just the rows we need
    row_idxs = [[]]*len(trial_periods)
    for i in range(len(trial_periods)):
        tp = trial_periods[i]
        row_idxs[i] = np.where(aov_data.p == tp)[0][0].copy()
    
    
    # Get the limits for the density graphs
    aov_max = np.max(aov_data.iloc[row_idxs, 1:])
    lsp_max = np.max(lsp_data.iloc[row_idxs, 1:])
    pdm_min = np.min(aov_data.iloc[row_idxs, 1:])
    
    
        # Plot the observed values in phase for several choices of trial period
    for i in range(len(trial_periods)):
        
        tp = trial_periods[i]
        
        # Plot the simulated data in phase
        plot_phase(scatter_data, tp, n_iter, iter_cols, bin_width, axs[0,i])
        
        
        # Plot a histogram of the AOV statistics
        idx = row_idxs[i]
        plot_aov(aov_data.iloc[idx, 1:], scatter_data, tp, bin_width, aov_max, axs[1,i])
        
        
        # Plot a histogram of the PDM statistics
        plot_pdm(pdm_data.iloc[idx, 1:], scatter_data, tp, bin_width, pdm_min, axs[2,i])
        
    
        # Plot a histogram of the SL statistics
        plot_sl(sl_data.iloc[idx, 1:], 'SL Statistic', axs[3,i])
        
        
        # Plot a histogram of the Lomb-Scargle Power statistics
        plot_lsp(lsp_data.iloc[idx, 1:], lsp_max, axs[4,i])
        
        
        # Add in the column titles
        axs[0,i].set_title('Trial Period ' + '{:.3f}'.format(tp))
        
   
    # Set the y labels
    axs[0, 0].set_ylabel(ylab)
    axs[1, 0].set_ylabel('Density')
    axs[2, 0].set_ylabel('Density')
    axs[3, 0].set_ylabel('Density')
    axs[4, 0].set_ylabel('Density')
    fig.align_labels()
    
    plt.savefig(output_path + name + '.pdf', bbox_inches = 'tight', dpi = 600)
    plt.close()



def plot_max_stat(data, nsim, ax = None):
    
    max_stat = data.iloc[:, 1:].max(axis = 0).copy()
    # Plot the histogram of the extreme value
    ax.hist(max_stat, bins = floor(sqrt(nsim)),
             density = True, color = 'lightgray')
    
    return ax


def plot_min_stat(data, nsim, ax = None):
    
    min_stat = data.iloc[:, 1:].min(axis = 0).copy()
    
    # Plot the histogram of the extreme value
    ax.hist(min_stat, bins = floor(sqrt(nsim)),
             density = True, color = 'lightgray')
    
    return ax



def plot_h0_extremevals(aov1, lsp1, pdm1, sl1, aov2, lsp2, pdm2, sl2, output_path, name):
    
    nsims = aov1.shape[1]
    

    fig = plt.figure()
    fig.set_figheight(15)
    fig.set_figwidth(40)
    widths = [4, 4, 4, 4]
    heights = [1, 1]
    gs = fig.add_gridspec(nrows = 2, ncols = 4, wspace=0.22, hspace = 0.1,
                          width_ratios=widths, height_ratios=heights,
                          figure = fig)
    
    
    ax_aov1 = fig.add_subplot(gs[0, 0])
    ax_aov2 = fig.add_subplot(gs[1, 0], sharex = ax_aov1, sharey = ax_aov1)
    
    ax_lsp1 = fig.add_subplot(gs[0, 1])
    ax_lsp2 = fig.add_subplot(gs[1, 1], sharex = ax_lsp1, sharey = ax_lsp1)
    
    ax_pdm1 = fig.add_subplot(gs[0, 2])
    ax_pdm2 = fig.add_subplot(gs[1, 2], sharex = ax_pdm1, sharey = ax_pdm1)
    
    ax_sl1 = fig.add_subplot(gs[0, 3])
    ax_sl2 = fig.add_subplot(gs[1, 3], sharex = ax_sl1, sharey = ax_sl1)
    
    
    
        # Set 1 Observations
    plot_max_stat(aov1, nsims, ax_aov1)
    plot_max_stat(lsp1, nsims, ax_lsp1)
    plot_min_stat(pdm1, nsims, ax_pdm1)
    plot_min_stat(sl1, nsims, ax_sl1)
    
    
        # Set 2 Observations
    plot_max_stat(aov2, nsims, ax_aov2)
    plot_max_stat(lsp2, nsims, ax_lsp2)
    plot_min_stat(pdm2, nsims, ax_pdm2)
    plot_min_stat(sl2, nsims, ax_sl2)
    
    
    ax_aov1.set_title('Maximum Observed AOV')
    ax_lsp1.set_title('Maximum Observed LS Power')
    ax_pdm1.set_title('Minimum Observed PDM')
    ax_sl1.set_title('Minimum Observed SL')
    
    ax_aov1.set_ylabel('Density')
    ax_aov2.set_ylabel('Density')
    
    ax_aov2.set_xlabel('AOV Statistic')
    ax_lsp2.set_xlabel('LS Power')
    ax_pdm2.set_xlabel('PDM Statistic')
    ax_sl2.set_xlabel('SL Statistic')
    
    fig.align_labels()
    
    plt.savefig(output_path + name + '.pdf', bbox_inches = 'tight', dpi = 600)
    plt.close()


def plot_h1_extremevals(aov, lsp, pdm, sl, output_path, name):
    
    nsims = aov.shape[1]
    

    fig = plt.figure()
    fig.set_figheight(8)
    fig.set_figwidth(40)
    widths = [4, 4, 4, 4]
    heights = [1]
    gs = fig.add_gridspec(nrows = 1, ncols = 4, wspace=0.22,
                          width_ratios=widths, height_ratios=heights,
                          figure = fig)
    
    
    ax_aov = fig.add_subplot(gs[0])
    ax_lsp = fig.add_subplot(gs[1])
    ax_pdm = fig.add_subplot(gs[2])
    ax_sl = fig.add_subplot(gs[3])
    
    
    plot_max_stat(aov, nsims, ax_aov)
    plot_max_stat(lsp, nsims, ax_lsp)
    plot_min_stat(pdm, nsims, ax_pdm)
    plot_min_stat(sl, nsims, ax_sl)
    
    
    ax_aov.set_title('Maximum Observed AOV')
    ax_lsp.set_title('Maximum Observed LS Power')
    ax_pdm.set_title('Minimum Observed PDM')
    ax_sl.set_title('Minimum Observed SL')
    
    ax_aov.set_ylabel('Density')
    
    ax_aov.set_xlabel('AOV Statistic')
    ax_lsp.set_xlabel('LS Power')
    ax_pdm.set_xlabel('PDM Statistic')
    ax_sl.set_xlabel('SL Statistic')
    
    fig.align_labels()
    
    plt.savefig(output_path + name + '.pdf', bbox_inches = 'tight', dpi = 600)
    plt.close()
    
    
def plot_periodogram(df, vline = False, p = None, ax = None):
    
    # Calculate the mean and standard deviation of the statistic
        # for a specific period
    stat_mean = df.iloc[:, 1:].mean(axis = 1)
    stat_sd = df.iloc[:, 1:].std(axis = 1)
    
    # Put into a data frame for plotting purposes
    periodogram_data = pd.DataFrame({'p': df.p, 'stat_mean': stat_mean, 'stat_sd': stat_sd})
    
    # Add in the monte carlo 1 sd error bounds
    periodogram_data['l'] = periodogram_data.stat_mean - periodogram_data.stat_sd
    periodogram_data['u'] = periodogram_data.stat_mean + periodogram_data.stat_sd

    # Sort for plotting purposes
    periodogram_data_sorted = periodogram_data.sort_values(by = 'p').copy()
    
    # Plot
    ax.plot(periodogram_data_sorted.p, periodogram_data_sorted.stat_mean, c = 'black')
    #ax.fill_between(periodogram_data_sorted.p,
    #                periodogram_data_sorted.l,
    #                periodogram_data_sorted.u, color = 'black', alpha=0.2)
    ax.plot(periodogram_data_sorted.p, periodogram_data_sorted.l, c = 'blue')
    ax.plot(periodogram_data_sorted.p, periodogram_data_sorted.u, c = 'blue')
    if vline:
        ax.axvline(p, c = 'red')
    
    return ax


def plot_h0_periodograms(aov1, lsp1, pdm1, sl1, aov2, lsp2, pdm2, sl2, output_path, name):

    fig = plt.figure()
    fig.set_figheight(12)
    fig.set_figwidth(35)
    widths = [4, 4, 4, 4]
    heights = [1, 1]
    gs = fig.add_gridspec(nrows = 2, ncols = 4, wspace=0.25, hspace = 0.1,
                          width_ratios=widths, height_ratios=heights)
    axs = gs.subplots(sharex='col', sharey='col')
    
        # Evenly Spaced Observations
    plot_periodogram(aov1, ax = axs[0,0])
    plot_periodogram(lsp1, ax = axs[0,1])
    plot_periodogram(pdm1, ax = axs[0,2])
    plot_periodogram(sl1, ax = axs[0,3])
    
    
        # Unevenly Spaced Observations
    plot_periodogram(aov2, ax = axs[1,0])
    plot_periodogram(lsp2, ax = axs[1,1])
    plot_periodogram(pdm2, ax = axs[1,2])
    plot_periodogram(sl2, ax = axs[1,3])
    
    
    axs[0, 0].set_title('AOV Periodogram')
    axs[0, 1].set_title('LS Power Periodogram')
    axs[0, 2].set_title('PDM Periodogram')
    axs[0, 3].set_title('SL Periodogram')
    
    axs[1, 0].set_xlabel('Trial P')
    axs[1, 1].set_xlabel('Trial P')
    axs[1, 2].set_xlabel('Trial P')
    axs[1, 3].set_xlabel('Trial P')
    fig.align_labels()
    
    
    plt.savefig(output_path + name + '.pdf', bbox_inches = 'tight', dpi = 600)
    plt.close()
    
    
def plot_h1_periodograms(aov, lsp, pdm, sl, output_path, name):

    fig = plt.figure()
    fig.set_figheight(8)
    fig.set_figwidth(35)
    widths = [4, 4, 4, 4]
    heights = [1]
    gs = fig.add_gridspec(nrows = 1, ncols = 4, wspace=0.25,
                          width_ratios=widths, height_ratios=heights)
    axs = gs.subplots(sharex='col', sharey='col')
    
    plot_periodogram(aov, True, 150, ax = axs[0])
    plot_periodogram(lsp, True, 150, ax = axs[1])
    plot_periodogram(pdm, True, 150, ax = axs[2])
    plot_periodogram(sl, True, 150, ax = axs[3])
    
    
    axs[0].set_title('AOV Periodogram')
    axs[1].set_title('LS Power Periodogram')
    axs[2].set_title('PDM Periodogram')
    axs[3].set_title('SL Periodogram')
    
    axs[0].set_xlabel('Trial P')
    axs[1].set_xlabel('Trial P')
    axs[2].set_xlabel('Trial P')
    axs[3].set_xlabel('Trial P')
    fig.align_labels()
    
    
    plt.savefig(output_path + name + '.pdf', bbox_inches = 'tight', dpi = 600)
    plt.close()


###############################################################################
# Read in the data
###############################################################################

#path = 'Documents/Astrostatistics/PD_review_paper/new_code/' 
path = ''

# Read in the simulated data from Gen_Data.py
null_even = pd.read_pickle(path + 'inter/null_even.pkl')
null_uneven = pd.read_pickle(path + 'inter/null_uneven.pkl')
null_hetero_trans = pd.read_pickle(path + 'inter/null_hetero_trans.pkl')
null_hetero_sim = pd.read_pickle(path + 'inter/null_hetero_sim.pkl')
h1 = pd.read_pickle(path + 'inter/h1.pkl')

# Read in the statistics
    # AOV
null_even_aov = pd.read_pickle(path + 'output/sims_data/aov_allnull_even.pkl')
null_uneven_aov = pd.read_pickle(path + 'output/sims_data/aov_allnull_uneven.pkl')
null_hetero_trans_aov = pd.read_pickle(path + 'output/sims_data/aov_allnull_hetero_trans.pkl')
null_hetero_sim_aov = pd.read_pickle(path + 'output/sims_data/aov_allnull_hetero_sim.pkl')
h1_aov = pd.read_pickle(path + 'output/sims_data/aov_allh1.pkl')

    # PDM
null_even_pdm = pd.read_pickle(path + 'output/sims_data/pdm_allnull_even.pkl')
null_uneven_pdm = pd.read_pickle(path + 'output/sims_data/pdm_allnull_uneven.pkl')
null_hetero_trans_pdm = pd.read_pickle(path + 'output/sims_data/pdm_allnull_hetero_trans.pkl')
null_hetero_sim_pdm = pd.read_pickle(path + 'output/sims_data/pdm_allnull_hetero_sim.pkl')
h1_pdm = pd.read_pickle(path + 'output/sims_data/pdm_allh1.pkl')


    # SL
null_even_sl = pd.read_pickle(path + 'output/sims_data/sl_allnull_even.pkl')
null_uneven_sl = pd.read_pickle(path + 'output/sims_data/sl_allnull_uneven.pkl')
null_hetero_trans_sl = pd.read_pickle(path + 'output/sims_data/sl_allnull_hetero_trans.pkl')
null_hetero_sim_sl = pd.read_pickle(path + 'output/sims_data/sl_allnull_hetero_sim.pkl')
h1_sl = pd.read_pickle(path + 'output/sims_data/sl_allh1.pkl')


    # Lomb-Scargle Power
null_even_ls_power = pd.read_pickle(path + 'output/sims_data/ls_power_scaled_all_null_even.pkl')
null_uneven_ls_power = pd.read_pickle(path + 'output/sims_data/ls_power_scaled_all_null_uneven.pkl')
null_hetero_trans_ls_power = pd.read_pickle(path + 'output/sims_data/ls_power_scaled_all_null_hetero_trans.pkl')
null_hetero_sim_ls_power = pd.read_pickle(path + 'output/sims_data/ls_power_scaled_all_null_hetero_sim.pkl')
h1_ls_power = pd.read_pickle(path + 'output/sims_data/ls_power_scaled_all_h1.pkl')



# Change the column names of the statistics for iterating over later
colnames = null_even.columns.tolist()
colnames.remove('obs_times')
colnames.insert(0, 'p')

null_even_aov.columns = colnames
null_uneven_aov.columns = colnames
null_hetero_trans_aov.columns = colnames
null_hetero_sim_aov.columns = colnames
h1_aov.columns = colnames

null_even_pdm.columns = colnames
null_uneven_pdm.columns = colnames
null_hetero_trans_pdm.columns = colnames
null_hetero_sim_pdm.columns = colnames
h1_pdm.columns = colnames


null_even_sl.columns = colnames
null_uneven_sl.columns = colnames
null_hetero_trans_sl.columns = colnames
null_hetero_sim_sl.columns = colnames
h1_sl.columns = colnames


null_even_ls_power.columns = colnames
null_uneven_ls_power.columns = colnames
null_hetero_trans_ls_power.columns = colnames
null_hetero_sim_ls_power.columns = colnames
h1_ls_power.columns = colnames

    # save the columns we'll iterate over
icols = colnames.copy()
icols.remove('p')


# Save the sequence of trial p-values
p_seq = null_even_aov.p.copy()



###############################################################################
# Plot the data
###############################################################################

plots_output = path + 'output/sims_plots/'


############################################
# Plot the scatter plots of the observed data
############################################

plt.rcParams.update({'font.size': 15})

# Plot one of the sims without phase-folding
plot_obs_data(null_even, 'Scaled G-Magnitude', plots_output, 'null_even')
plot_obs_data(null_uneven, 'Scaled G-Magnitude', plots_output, 'null_uneven')
plot_obs_data(null_hetero_trans, 'G-Magnitude', plots_output, 'null_hetero_trans')
plot_obs_data(null_hetero_sim, 'G-Magnitude', plots_output, 'null_hetero_sim')
plot_obs_data(h1, 'Scaled G-Magnitude', plots_output, 'h1')



############################################
# Plot the figures of the phase-folded data with the AOV, PDM, SL, and LS statistics
############################################

plt.rcParams.update({'font.size': 25})

# Make the set of trial periods
    # period closest to 75
i = np.where(np.abs(p_seq - 75) == np.min(np.abs(p_seq - 75)))[0][0].copy()
p_75 = p_seq[i]

    # period near 83
p_83 = p_seq[702]

trial_ps = [p_75, p_83, 150]

# Set the bin width
bw = 0.1


# H0, even, homoscedastic
plot_stat_hists(null_even, null_even_aov, null_even_pdm, null_even_sl, null_even_ls_power,
                trial_ps, icols, bw, 'Scaled G-Magnitude', 
                plots_output, 'null_even_stat_hists')


# H0, uneven, homoscedastic
plot_stat_hists(null_uneven, null_uneven_aov, null_uneven_pdm, null_uneven_sl, null_uneven_ls_power,
                trial_ps, icols, bw, 'Scaled G-Magnitude', 
                plots_output, 'null_uneven_stat_hists')


# H0, uneven, heteroscedastic noise from not scaling
plot_stat_hists(null_hetero_trans, null_hetero_trans_aov, null_hetero_trans_pdm,
                null_hetero_trans_sl, null_hetero_trans_ls_power,
                trial_ps, icols, bw, 'G-Magnitude', 
                plots_output, 'null_hetero_trans_stat_hists')


# H0, uneven, heteroscedastic increasing noise
plot_stat_hists(null_hetero_sim, null_hetero_sim_aov, null_hetero_sim_pdm,
                null_hetero_sim_sl, null_hetero_sim_ls_power,
                trial_ps, icols, bw, 'G-Magnitude', 
                plots_output, 'null_hetero_sim_stat_hists')


# H1
plot_stat_hists(h1, h1_aov, h1_pdm, h1_sl, h1_ls_power,
                trial_ps, icols, bw, 'Scaled G-Magnitude', 
                plots_output, 'h1_stat_hists')



############################################
# Plot the figures of the maximum AOV and Lomb-Scargle Power and the minimum PDM and SL statistics
############################################

plt.rcParams.update({'font.size': 30})


# Extreme values histograms for evenly- and unevenly-spaced observations, homoscedastic noise
plot_h0_extremevals(null_even_aov, null_even_ls_power, null_even_pdm, null_even_sl,
                    null_uneven_aov, null_uneven_ls_power, null_uneven_pdm, null_uneven_sl,
                    plots_output, 'h0_extreme_homoscedastic_hists')



# Extreme values histograms for unevenly-spaced observations, heteroscedastic noise
plot_h0_extremevals(null_hetero_trans_aov, null_hetero_trans_ls_power, null_hetero_trans_pdm, null_hetero_trans_sl,
                    null_hetero_sim_aov, null_hetero_sim_ls_power, null_hetero_sim_pdm, null_hetero_sim_sl,
                    plots_output, 'h0_extreme_heteroscedastic_hists')


# Extreme values histograms for observations under Ha
plot_h1_extremevals(h1_aov, h1_ls_power, h1_pdm, h1_sl, plots_output, 'h1_extreme_hists')




############################################
# Plot the periodograms of the AOV, PDM, and SL statistics and LS power
############################################

plt.rcParams.update({'font.size': 30})


# Periograms for evenly- and unevenly-spaced observations, homoscedastic noise
plot_h0_periodograms(null_even_aov, null_even_ls_power, null_even_pdm, null_even_sl,
                     null_uneven_aov, null_uneven_ls_power, null_uneven_pdm, null_uneven_sl,
                     plots_output, 'h0_homoscedastic_periodograms')


# Periograms for unevenly-spaced observations, heteroscedastic noise
plot_h0_periodograms(null_hetero_trans_aov, null_hetero_trans_ls_power, null_hetero_trans_pdm, null_hetero_trans_sl,
                     null_hetero_sim_aov, null_hetero_sim_ls_power, null_hetero_sim_pdm, null_hetero_sim_sl,
                     plots_output, 'h0_heteroscedastic_periodograms')



# Periodograms for observations under Ha
plot_h1_periodograms(h1_aov, h1_ls_power, h1_pdm, h1_sl, 
                     plots_output, 'h1_periodograms')





