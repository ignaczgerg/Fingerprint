import numpy as np
import matplotlib
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

'''
initialData_1 = pd.read_csv('calcFDA.csv')
initialData_2 = pd.read_csv('calcOdyes.csv')
initialData_3 = pd.read_csv('calcPS_oligomers.csv')
initialData_4 = pd.read_csv('calcUNPD_1.csv')
plt.figure(figsize=(10,7), dpi= 300)
initialData_4['calcTPSA'].hist(bins = 150)
initialData_1['calcTPSA'].hist(bins = 150)
#initialData_1['calcTPSA'].hist(bins = 100)
#initialData_3['calcTPSA'].hist(bins = 13)
'''
'''
################################################################
### This is a normal plotting for a single frequency diagram ###
################################################################
initialData = pd.read_csv('calcUNPD_1.csv')
plt.rcParams.update({'figure.figsize':(10,7), 'figure.dpi':100})

molweight, bins = np.histogram(initialData['Molweight'], bins=150, range=[0, 1500])
plt.hist(molweight, bins=100)
plt.gca().set(title='Frequency Histogram', ylabel='Frequency');
'''



##################################################
### This is a Seaborn histogram of two dataset ###
##################################################

sns.set_style('white')
initialData_1 = pd.read_csv('calcPS_oligomers.csv')
#initialData_2 = pd.read_csv('calcOdyes.csv')
#initialData_3 = pd.read_csv('calcPS_oligomers.csv')
#initialData_4 = pd.read_csv('calcUNPD_1.csv')
#initialData_5 = pd.read_csv('calcPEG.csv')
#initialData_6 = pd.read_csv('calcPSamide_mod.csv')

kwargs = dict(hist_kws={'alpha':0.5}, kde_kws={'linewidth':2})


plt.figure(figsize=(10,7), dpi= 300)
sns.distplot(initialData_1['calcSpherocityIndex'], color="dodgerblue", label="Drug molecules", **kwargs, kde=True, hist=False)
#sns.distplot(initialData_2['calcSpherocityIndex'], color="orange", label="Organic dyes", **kwargs, bins=15, kde=True, hist=False)
#sns.distplot(initialData_3['Molweight'], color="pink", label="PS markers", **kwargs, bins=10, kde=False)
#sns.distplot(initialData_4['calcSpherocityIndex'], color="lightgreen", label="UNPD", **kwargs, kde=True, bins=200, hist=False)
#sns.distplot(initialData_5['Molweight'], color="yellow", label="PEG markers", **kwargs, bins=15, kde=False)
#sns.distplot(initialData_6['calcSpherocityIndex'], color="grey", label="PS mod", **kwargs, kde=True, hist=False, bins=120)

plt.legend()

