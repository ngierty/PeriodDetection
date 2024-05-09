#####################################################
# Paths and Libraries
#####################################################

import matplotlib.pyplot as plt
import numpy as np
from math import pi, sin
import pandas as pd

#####################################################
### Phase-folding Example
#####################################################

SMALL_SIZE = 30
MEDIUM_SIZE = 30

plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=SMALL_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels

x = np.linspace(0, 6*pi, num = 100)
x = x[x < 6*pi] # removes the one awkard point at 6pi when plotting
sinx = [sin(t) for t in x]
y = np.random.normal(sinx, 0.2)
phasedf = pd.DataFrame({'x': x, 'y': y})
phasedf['color_label'] = np.floor(x/(2*pi))
groups = phasedf.groupby('color_label')

markers = ['o', '^', 's']


# Set up the figure
fig = plt.figure()
fig.set_figheight(8)
fig.set_figwidth(30)
widths = [4, 4, 4]
heights = [1]
gs = fig.add_gridspec(nrows = 1, ncols = 3, wspace=0.05,
                      width_ratios=widths, height_ratios=heights,
                      figure = fig)
axs = gs.subplots(sharex=None, sharey='row')



# Plot without phase folding
axs[0].margins(0.05) # Optional, just adds 5% padding to the autoscaling
k = 0
for name, group in groups:
    axs[0].plot(group.x, group.y, marker=markers[k], linestyle='', label=name)
    k += 1
axs[0].set_xlabel('Time (Days)')
axs[0].set_ylabel('Observed Value')


# Plot with phase folding at incorrect period
phasedf['phase'] = (x % (pi/3)) / (pi/3)
axs[1].margins(0.05) # Optional, just adds 5% padding to the autoscaling
k = 0
for name, group in groups:
    axs[1].plot(group.phase, group.y, marker=markers[k], linestyle='', label=name)
    k += 1
axs[1].set_xlabel('Phase-Folded Time')


# Plot with phase folding at correct period
phasedf['phase'] = (x % (2*pi)) / (2*pi)
axs[2].margins(0.05) # Optional, just adds 5% padding to the autoscaling
k = 0
for name, group in groups:
    axs[2].plot(group.phase, group.y, marker=markers[k], linestyle='', label=name)
    k += 1
axs[2].set_xlabel('Phase-Folded Time')



path = 'Documents/Astrostatistics/PD_review_paper/new_code/' 

plt.savefig(path + 'output/phase_folding.pdf', bbox_inches = 'tight', dpi = 600)
plt.close()











































