##############
## PLOTTING ##
##############

import matplotlib 
from matplotlib import pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import gaussian_kde


initialData = pd.read_csv('calcPSamide_oligomers.csv', sep=',')


pointDensity = np.vstack([initialData['calcNPR1'], initialData['calcNPR2']])
height = gaussian_kde(pointDensity)(pointDensity)

kwargs = dict(xlabel = 'NPR1')

plt.figure(figsize=(10,7), dpi= 300)
matplotlib.pyplot.ylim(0.5,1)
matplotlib.pyplot.xlim(0,1)
graph = plt.scatter(initialData['calcNPR1'], initialData['calcNPR2'], c=height, s=40, edgecolor='')

#plt.scatter(initialData['calcNPR1'], initialData['calcNPR2'])
#matplotlib.pyplot.ylim(0.5,1)
#matplotlib.pyplot.xlim(0,1)
