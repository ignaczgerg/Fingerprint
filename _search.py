import rdkit
from rdkit import Chem
from rdkit import DataStructs
import numpy as np
import pandas as pd
import matplotlib
import seaborn 

initialData = pd.read_csv('UNPD1to11.csv', sep=';')
df = pd.DataFrame(initialData, columns = ['Column1', 'Column30'])
df['Column30'].replace('.' or '', np.nan, inplace=True)
df.dropna(subset=['Column30'], inplace=True)

ms = []
fps = []
finger = []

for i in range(len(df)):
    ms.append(df.iloc[i,1])
    fps.append(Chem.MolFromSmiles(ms[i]))
    finger.append(Chem.RDKFingerprint(fps[i]))

