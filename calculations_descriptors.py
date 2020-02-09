import rdkit
from rdkit import Chem
from rdkit.Chem import Descriptors3D
from rdkit.Chem import Descriptors
from rdkit.Chem import Crippen
from rdkit.Chem import Lipinski
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem.rdmolfiles import SDMolSupplier
import pandas as pd
import numpy as np



## Molweight
molweight = []
## 3D Descriptors
calcAsphericity = []
calcEccentricity = []
calcInertialShapeFactor = []
calcRadiusOfGyration = []
calcSpherocityIndex = []
calcPMI1  = []
calcPMI2  = []
calcPMI3  = []
calcNPR1  = []
calcNPR2  = []
## Crippen module
calcCrippenLopP = []
calcCrippenMR = []
## Lipinski module
calcLipinskiFractionCSP3 = []
calcLipinskiHeavyAtomCount = []
calcLipinskiNHOHCount = []
calcLipinskiNOCount = []
calcLipinskiNumAliphaticCarbocycles = []
calcLipinskiNumAliphaticHeterocycles = []
calcLipinskiNumAliphaticRings = []
calcLipinskiNumAromaticCarbocycles = []
calcLipinskiNumAromaticHeterocycles = []
calcLipinskiNumAromaticRings = []
calcLipinskiNumHAcceptors = []
calcLipinskiNumHDonors = []
calcLipinskiNumHeteroatoms = []
calcLipinskiNumRotatableBonds = []
calcLipinskiNumSaturatedCarbocycles = []
calcLipinskiNumSaturatedHeterocycles = []
calcLipinskiNumSaturatedRings = []
calcLipinskiRingCount = []
## SMILES
backSmiles = []
## ID
IDNumber = []
## Free surface area
calcCalcSASA = [] # Compute the Solvent Accessible Surface Area 
calcTPSA = []

suppl = SDMolSupplier("dimer.sdf")
df = pd.DataFrame(columns = ['Name', 'SMILES'])


for mol in suppl:
    IDNumber.append(mol)
    molweight.append(rdMolDescriptors.CalcExactMolWt(mol))
    ##
    calcPMI1.append(rdkit.Chem.Descriptors3D.PMI1(mol))
    calcPMI2.append(rdkit.Chem.Descriptors3D.PMI2(mol))
    calcPMI3.append(rdkit.Chem.Descriptors3D.PMI3(mol))
    calcNPR1.append(rdkit.Chem.Descriptors3D.NPR1(mol))
    calcNPR2.append(rdkit.Chem.Descriptors3D.NPR2(mol))
    ##
    calcCrippenLopP.append(rdkit.Chem.Crippen.MolLogP(mol))
    calcCrippenMR.append(rdkit.Chem.Crippen.MolMR(mol))
    ##
    calcLipinskiFractionCSP3.append(rdkit.Chem.Lipinski.FractionCSP3(mol))
    calcLipinskiHeavyAtomCount.append(rdkit.Chem.Lipinski.HeavyAtomCount(mol))
    calcLipinskiNHOHCount.append(rdkit.Chem.Lipinski.NHOHCount(mol))
    calcLipinskiNOCount.append(rdkit.Chem.Lipinski.NOCount(mol))
    calcLipinskiNumAliphaticCarbocycles.append(rdkit.Chem.Lipinski.NumAliphaticCarbocycles(mol))
    calcLipinskiNumAliphaticHeterocycles.append(rdkit.Chem.Lipinski.NumAliphaticHeterocycles(mol))
    calcLipinskiNumAliphaticRings.append(rdkit.Chem.Lipinski.NumAliphaticRings(mol))
    calcLipinskiNumAromaticCarbocycles.append(rdkit.Chem.Lipinski.NumAromaticCarbocycles(mol))
    calcLipinskiNumAromaticHeterocycles.append(rdkit.Chem.Lipinski.NumAromaticHeterocycles(mol))
    calcLipinskiNumAromaticRings.append(rdkit.Chem.Lipinski.NumAromaticRings(mol))
    calcLipinskiNumHAcceptors.append(rdkit.Chem.Lipinski.NumHAcceptors(mol))
    calcLipinskiNumHDonors.append(rdkit.Chem.Lipinski.NumHDonors(mol))
    calcLipinskiNumHeteroatoms.append(rdkit.Chem.Lipinski.NumHeteroatoms(mol))
    calcLipinskiNumRotatableBonds.append(rdkit.Chem.Lipinski.NumRotatableBonds(mol))
    calcLipinskiNumSaturatedCarbocycles.append(rdkit.Chem.Lipinski.NumSaturatedCarbocycles(mol))
    calcLipinskiNumSaturatedHeterocycles.append(rdkit.Chem.Lipinski.NumSaturatedHeterocycles(mol))
    calcLipinskiNumSaturatedRings.append(rdkit.Chem.Lipinski.NumSaturatedRings(mol))
    calcLipinskiRingCount.append(rdkit.Chem.Lipinski.RingCount(mol))
    ##
    calcAsphericity.append(rdkit.Chem.Descriptors3D.Asphericity(mol))
    calcEccentricity.append(rdkit.Chem.Descriptors3D.Eccentricity(mol))
    calcInertialShapeFactor.append(rdkit.Chem.Descriptors3D.InertialShapeFactor(mol))
    calcRadiusOfGyration.append(rdkit.Chem.Descriptors3D.RadiusOfGyration(mol))
    calcSpherocityIndex.append(rdkit.Chem.Descriptors3D.SpherocityIndex(mol))
    #calcCalcSASA.append(rdkit.Chem.rdFreeSASA.CalcSASA(mol))
    calcTPSA.append(rdkit.Chem.Descriptors.TPSA(mol))
    backSmiles.append(Chem.MolToSmiles(mol))

df['Name'] = IDNumber
df['SMILES'] = backSmiles
df['Molweight'] = molweight
df['calcPM1'] = calcPMI1
df['calcPM2'] = calcPMI2
df['calcPM3'] = calcPMI3
df['calcNPR1'] = calcNPR1
df['calcNPR2'] = calcNPR2
df['calcCrippenLopP'] = calcCrippenLopP
df['calcCrippenMR'] = calcCrippenMR
df['calcLipinskiFractionCSP3'] = calcLipinskiFractionCSP3
df['calcLipinskiHeavyAtomCount'] = calcLipinskiHeavyAtomCount
df['calcLipinskiNHOHCount'] = calcLipinskiNHOHCount
df['calcLipinskiNOCount'] = calcLipinskiNOCount
df['calcLipinskiNumAliphaticCarbocycles'] = calcLipinskiNumAliphaticCarbocycles
df['calcLipinskiNumAliphaticHeterocycles'] = calcLipinskiNumAliphaticHeterocycles
df['calcLipinskiNumAliphaticRings'] = calcLipinskiNumAliphaticRings
df['calcLipinskiNumAromaticCarbocycles'] = calcLipinskiNumAromaticCarbocycles
df['calcLipinskiNumAromaticHeterocycles'] = calcLipinskiNumAromaticHeterocycles
df['calcLipinskiNumAromaticRings'] = calcLipinskiNumAromaticRings
df['calcLipinskiNumHAcceptors'] = calcLipinskiNumHAcceptors
df['calcLipinskiNumHDonors'] = calcLipinskiNumHDonors
df['calcLipinskiNumHeteroatoms'] = calcLipinskiNumHeteroatoms
df['calcLipinskiNumRotatableBonds'] = calcLipinskiNumRotatableBonds
df['calcLipinskiNumSaturatedCarbocycles'] = calcLipinskiNumSaturatedCarbocycles
df['calcLipinskiNumSaturatedHeterocycles'] = calcLipinskiNumSaturatedHeterocycles
df['calcLipinskiNumSaturatedRings'] = calcLipinskiNumSaturatedRings
df['calcLipinskiRingCount'] = calcLipinskiRingCount
df['calcAsphericity'] = calcAsphericity
df['calcEccentricity'] = calcEccentricity
df['calcInertialShapeFactor'] = calcInertialShapeFactor
df['calcRadiusOfGyration'] = calcRadiusOfGyration
df['calcSpherocityIndex'] = calcSpherocityIndex
#df['calcCalcSASA'] = calcCalcSASA
df['calcTPSA'] = calcTPSA



df.to_csv('calcPSamide_mod.csv')
