# Original source-code from CIME's repository:
# Author: Thomas Wolf
# https://github.com/jku-vds-lab/cime/blob/main/Examples/CimeExample.ipynb


#####################################
### import the used functions
#####################################
#for umap you need the umap-learn package
# import umap.umap_ as umap
from rdkit import Chem
import numpy as np
from rdkit.Chem import AllChem, DataStructs
import pandas as pd

import sys
sys.path.append('./')

#######################################
###define the used functions
#######################################

#define a function that trnasforms smiles RDKit Mol objects
def smiles_to_mols(smis):
    """Convert a list of smiles to a list of RDKit Mol objects"""
    if not isinstance(smis, list):
        raise TypeError("Please provide smiles as a list")
    mols = []
    successful_inds = []
    # print("smis", smis)
    for ind, smi in enumerate(smis):
        # RDKit will sometimes 'fail' but instead of throwing an
        # error it will (sometimes!) print a warning and return None.
        # print('Smiles:', ind, smi)
        m = Chem.MolFromSmiles(smi, sanitize=True)
        if m is None:
            print("Mol generation failed for", smi, "(index", ind, ")\n")
        else:
            mols.append(m)
            successful_inds.append(ind)
    return mols, successful_inds


# define a function to calculate count based morgan fingerprints
def calc_morgan_counts_for_mol(mol, radius=1, nBits=2048):
    try:
        fpCounts = AllChem.GetHashedMorganFingerprint(
            mol, radius, nBits=nBits, includeRedundantEnvironments=False)
        array = np.zeros((0,), dtype=np.int8)
        DataStructs.ConvertToNumpyArray(fpCounts, array)
    except:
        print("Failure in calc_morgan_counts_for_mol()")
        array = np.repeat(np.nan, nBits)

    return array


#get the substructure indices for a specfied radius around an atomID
def getSubstructIndex(mol,atomID,radius):
    if radius>0:
        env = Chem.FindAtomEnvironmentOfRadiusN(mol,radius,atomID)
        atomsToUse=[]
        for b in env:
            atomsToUse.append(mol.GetBondWithIdx(b).GetBeginAtomIdx())
            atomsToUse.append(mol.GetBondWithIdx(b).GetEndAtomIdx())
        atomsToUse = list(set(atomsToUse))
    else:
        atomsToUse = [atomID]
        env=None
    return(atomsToUse)

#get the substructure indices for each fingerprint bit
def calc_indices_for_fingerprints(smiles,radius = 1,nBits=2048):
   mol = Chem.MolFromSmiles(smiles,sanitize = True)
   bi = {}
   fp = AllChem.GetHashedMorganFingerprint(mol, radius = radius, nBits=nBits, includeRedundantEnvironments=False,bitInfo=bi)
   fp = fp.GetNonzeroElements()
   #weights = [0  for i in mol.GetAtoms()]
   #weights[slice(mol.GetSubstructMatch(Chem.MolFromSmiles('CNC(=O)C(C)(C)C',sanitize = True)))] += 1


   fp_indices = {}
   for fragment_id, count in fp.items():
      sub_indices = []
      for i in range(count):
         #print(i)
         root, radius = bi[fragment_id][i]
         subindex = getSubstructIndex(mol,root,radius)
         sub_indices = sub_indices + subindex

      sub_indices = list(set(sub_indices))

      fp_indices[fragment_id] = sub_indices
   fp_indices = {"fp_morgan_counts_" +  str(k): v for k, v in fp_indices.items() }

   return(fp_indices)



#we also need a function that calculates the attribution of each atom to the final 
#prediction, for direction you can choose between "up": all positive attributions
#"down": all negative attributions and recommended "both": summing up over all attributions
def atom_attributions(mol, indices, shap_vals, featureNames, direction="both"):
    """
    Sum up the attributions to the model predictions (Shap values) for each atom in mol.

    Parameters:
        mol: rdkit molecule object
        indices (dict): a dictionary with keys = feature names and values = list of atom positions in mol (output of calc_indices_for_fingerprints())
        shap_vals: a numpy array of dimensions (len(featureNames))
        featureNames (np.array): numpy array containing the feature names for shap_vals
        direction (str): whether to consider positive (up), negative (down) or both shap_vals contributions

    Returns:
        np.array: each element = 1 atom of mol and the values are the summed contributions from shap_vals
    """
    ## TO DO: add support for MACCS keys

    # get for each fingerprint bit the atom numbers mapping to that bit and initialise weights to 0
    curWeights = np.array([0 for i in range(mol.GetNumAtoms())], 'float64')

    # step through each index (dictionary)
    for j in indices:
        # retrieve index to subset into the featureNames and the corresponding shap values array
        curFPIndex = [i for i, s in enumerate(featureNames) if s == j]
        if len(curFPIndex) == 0:
            continue
        if len(curFPIndex) > 1:
            raise IndexError("Found several matches for feature " + j)
        # do something only if (i) the fingerprint bit was in the selected features
        # and (ii) only if the retrieved shap value is positive/negative (direction argument)
        if (direction == "up" and shap_vals[curFPIndex[0]] > 0) or (direction == "down" and shap_vals[curFPIndex[0]] < 0) or (direction == "both"):
            # use the atom numbers ('positions') to update the curWeights array, where
            # each atom of the current molecule is represented
            # one atom may be part of several fingerprints so the summing has to be done across
            # all fingerprints
            for pos in indices[j]:
                #value = (curWeights[pos] + shap_vals[curFPIndex[0]])
                value = (curWeights[pos] + (shap_vals[curFPIndex[0]]/len(indices[j])))
                curWeights[pos] = value

    return curWeights

