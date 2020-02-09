import rdkit
from rdkit import Chem
from rdkit.Chem import AllChem
import pandas as pd
import numpy as np


initialData = pd.read_csv('dimer.csv', sep=',')
df = pd.DataFrame(initialData, columns = ['PMA_ID', 'SMILES'])
df['SMILES'].replace('.' or '', np.nan, inplace=True)
df.dropna(subset=['SMILES'], inplace=True)

optimizedStructures  = []
embedStructures = []


def StructureOpt(mol):
    mol = Chem.MolFromSmiles(mol)
    mol = Chem.AddHs(mol)
    if rdkit.Chem.rdDistGeom.EmbedMolecule(mol, maxAttempts=10, randomSeed = 0, useRandomCoords=True) == -1:
        mol = np.nan
    else:
        AllChem.MMFFOptimizeMolecule(mol, maxIters = 300)
    return(mol)    


def Iteration(df):
    for i in range(len(df)):
        optimizedStructures.append(StructureOpt(df.iloc[i,1]))
    return(optimizedStructures)

Iteration(df)

w = Chem.SDWriter('dimer.sdf') 
m = optimizedStructures.count(np.nan)      
while m > 0:
    optimizedStructures.remove(np.nan)
    m = m-1
for n in range(len(df)): w.write(optimizedStructures[n])

