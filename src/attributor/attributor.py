# Source-code Author: Linlin Zhao, 2022
# if you use this code, please cite our articles:
# xBCF - ...
# XSMILES - ...
#

import os
import re
import numpy as np
import pandas as pd
from random import randint
# import torch
# import torch.nn as nn
# import torch.nn.functional as F
#from tinydb import TinyDB, Query
from rdkit import Chem
from rdkit.Chem.Draw import SimilarityMaps
#import logomaker
import sys
# sys.path.append('./cddd/')
# import matplotlib.pyplot as plt
from cddd.inference import InferenceModel
# torch.manual_seed(124)

# use_gpu = torch.cuda.is_available()
use_gpu = False
# os.environ['CUDA_VISIBLE_DEVICES'] = '0'

char_dict = {
    0: '</s>',
    1: '#',
    2: '%',
    3: ')',
    4: '(',
    5: '+',
    6: '-',
    7: '1',
    8: '0',
    9: '3',
    10: '2',
    11: '5',
    12: '4',
    13: '7',
    14: '6',
    15: '9',
    16: '8',
    17: ':',
    18: '=',
    19: '@',
    20: 'C',
    21: 'B',
    22: 'F',
    23: 'I',
    24: 'H',
    25: 'O',
    26: 'N',
    27: 'P',
    28: 'S',
    29: '[',
    30: ']',
    31: 'c',
    32: 'i',
    33: 'o',
    34: 'n',
    35: 'p',
    36: 's',
    37: 'Cl',
    38: 'Br',
    39: '<s>'
    }

#cddd512_model = "cddd_default_model/"
#db = TinyDB('tinydb.json')

class Attributor:
    """The generalized attributor for any downstream models based on CDDD
    """
    def __init__(self, qsar_mdl, cddd_model, cpu_threads=2):
        self.cddd_model = cddd_model
        self.qsar_mdl = qsar_mdl

    def smiles_attribution(self, smiles, ids=None, method="perturb", plot=True):
        """The main interface for providing atom-wise and char-wise attributions


        Parameters
        ----------
        smiles : str or list of str
            The input SMILES strings (single or list), if it is a single string, the results will be
            printed or plotted to the screen. If it is a list, the results will be return as a
            dict for being stored as json files, which can also be viewed by SMILESdrawer.
        ids : list
            A list of identifiers (e.g. the name of the compounds) when smiles is a list. If
            provided, it should have the same length as smiles
        method : str, optional
            either "perturb" or "mask", by default "perturb"
        plot : bool, optional
            If True, it plots the atom attributions.
 
        Returns
        -------
        pandas dataframe
            return the attribution scores for each character in the input SMILES
        """
 
        # define char and atom lists for getting attributions
 
 
        if isinstance(smiles, str):
            # REGEX_SML = r'Cl|Br|[#%\)\(\+\-1032547698:=@CBFIHONPS\[\]cionps]'
            # char_list = re.findall(REGEX_SML, smiles)
            # char_list = re.findall(REGEX_SML, smiles)
 
            self._smiles_parser(smiles)
            attributions = self._get_attributions(smiles, method=method)
            chars = list(char_dict.values())[1:-1]
            chars_df = pd.concat(
                [
                    pd.DataFrame(self.char_list, columns=['Char']),
                    pd.DataFrame(attributions["bcf_attributions"], columns=['BCF_score']),
                    pd.DataFrame(attributions["lod_attributions"], columns=['LOD_score']),
                    pd.DataFrame(np.array(attributions["bcf_attributions"]) - 0.6*np.array(attributions["lod_attributions"]), columns=['Diff'])
                ], axis=1)
 
            # atom_chars = ["Cl","Br","C","B","F","I","H","O","N","P","S","c","i","o","n","p","s"]
            # atom_char_list = [i for i in range(len(char_list)) if char_list[i] in atom_chars]
            # print(f"atom_char_list={atom_char_list}")
            # true_char_list = [self.char_list[i] for i in self.atom_char_list]
            atom_bcf_scores = np.array(attributions["bcf_attributions"])[self.atom_char_list]
            atom_lod_scores = np.array(attributions["lod_attributions"])[self.atom_char_list]
            atom_diff = atom_bcf_scores - 0.6*atom_lod_scores
 
            atom_df = pd.concat(
                [
                    pd.DataFrame(self.true_char_list, columns=['Atom']),
                    pd.DataFrame(atom_bcf_scores, columns=['bcf_scores']),
                    pd.DataFrame(atom_lod_scores, columns=['lod_scores']),
                    pd.DataFrame(atom_diff, columns=['bcf_diff_lod'])
                ], axis=1)
 
            special_atom_dict = self._get_special_dict(attributions)
            if special_atom_dict:
                special_atom_df = pd.DataFrame.from_dict(special_atom_dict)
                print(f"Special atoms dataframe: {special_atom_df}")
 
            if plot:
                print(f"atom_id: {atom_df}")
                print(f"char_id: {chars_df}")
 
                #atom_logo(smiles, attributions)
                mol = Chem.MolFromSmiles(smiles)
                print("BCF attributions: ")
                fig = SimilarityMaps.GetSimilarityMapFromWeights(mol, atom_bcf_scores)
                print("LOD attributions:")
                fig = SimilarityMaps.GetSimilarityMapFromWeights(mol, atom_lod_scores)
                print("Score difference (BCF - LOD):")
                fig = SimilarityMaps.GetSimilarityMapFromWeights(mol, atom_diff)
            return chars_df
 
        elif isinstance(smiles, list):
 
            sml_attr = []
            special_atom_dfs = []
            for i, sml in enumerate(smiles):
                print(f"Processing {sml}")
                self._smiles_parser(sml)
                attributions = self._get_attributions(sml, method=method)

                attr = attributions['attributions']

                attributes_dict =  {
                    "attributes":
                    {
                        "predicted_value": attributions['predicted_value'],
                        f"max_attr_{method}": np.max(np.abs(attr)).item()
                    }
                }

                special_atom_dict = self._get_special_dict(attributions)

                # IF ids are given then update the attribute dict
                if ids is not None:
                    assert len(ids) == len(smiles)
                    attributes_dict['attributes'].update({"Compound ID": ids[i]})
                    if special_atom_dict:
                        special_atom_dict.update({"Compound ID": ids[i]})
                        special_atom_df = pd.DataFrame.from_dict(special_atom_dict)
                        special_atom_dfs.append(special_atom_df)
                res_dict = {
                        "string": sml,
                        "sequence": attributions['chars'],
                        "methods": [
                            {"name": f"Attribution {method}", "scores": attr}
                        ],
                        "attributes": attributes_dict['attributes']
                        }
                sml_attr.append(res_dict)
                #db.insert({'res': res_dict, 'smiles': sml, 'count': i})
            if special_atom_dfs:
                merge_special_df = pd.concat(special_atom_dfs, axis=0)
                return sml_attr, merge_special_df
            else:
                return sml_attr, "No special atoms"

        else:
            raise ValueError("The input SMILES must be a string or a list of strings.")

    def _smiles_parser(self, sml):
        """Parse the input smiles to get list of chars or atom
        """
        REGEX_SML = r'Cl|Br|[#%\)\(\+\-1032547698:=@CBFIHONPS\[\]cionps]'
        self.atom_chars = ["Cl","Br","C","B","F","I","H","O","N","P","S","c","i","o","n","p","s"]
        self.char_list = re.findall(REGEX_SML, sml)
        self.atom_char_list = [i for i in range(len(self.char_list)) if self.char_list[i] in self.atom_chars]
        self.true_char_list = [self.char_list[i] for i in self.atom_char_list]

    def _get_special_dict(self, attributions):
        subset_key = ['special_char', 'attr_special_lod',
                      'attr_special_bcf', 'special_char_idx',
                      'special_atom_idx', 'SMILES']
        special_atom_dict = {k: attributions[k] for k in subset_key if k in attributions}
        return special_atom_dict

    def _get_attributions(self, sml, method='perturb'):
        """Get character attributions for the given SMILES
        Two ways of getting attributions:
            1. mask: mask the on-position value
                    --> need to test dummy atom
            2. perturb: flip the off-position values one by one and get the average
 
        The flipping idea is what is the expected prediction change if right atom
        is replaced with any other atom.
 
        1. Compute original BCF values
        2. for each atom:
            a. Generate fake smiles strings for every character in the vocab
            b. Get embeddings for those smiles
            c. Compute BCF values according to embeddings
            d. Average BCF values
            g. Compute the difference as the atom attribution
            h. Return atom attributions
 
        The masking method is to estimate the effect of absence of certain atoms
 
        1. Compute original logBCF
        2. for each atom:
            a. Pass its index to cddd model creator for masking the atom
            b. Get the CDDD embeddings
            c. Pass embeddings to pytorch model
 
 
        Paramters
        -----------
        method : str
            Specify which method to perturb SMILES
 
        Returns
        -------
        1-d array
            An array of scores for all characters in the input SMILES
        """
        emb = self.cddd_model.seq_to_emb(sml)
        orig_y = self.qsar_mdl.predict(emb)
 
        char_vocab = list(char_dict.values())[1:-1]
 
        chars = self.char_list
 
        attributions = {
            "attributions": [],
            "chars": chars,
            "predicted_value": orig_y.item(),
        }
 
        if method == "perturb":
            for i, sml_char in enumerate(chars):
                sml_copy = chars
                # mutated_smls = [f"{sml_copy[:i]}{w}{sml_copy[i+1:]}" for w in chars]
                mutated_smls = ["".join(sml_copy[:i] + [w] + sml_copy[i+1:]) for w in char_vocab]
                # mutated_smls.append(sml)
                mutated_emb = self.cddd_model.seq_to_emb(mutated_smls)
                # x = torch.tensor(mutated_emb).float()
                mutated_y = self.qsar_mdl.predict(mutated_emb)
                # print(y)
                # orig_bcf = y[0][-1]
                # orig_lod = y[1][-1]
                # exclude the value for original sml
                # attr_dict[i]
                attributions['attributions'].append((orig_y - np.mean(mutated_y)).item())
                # attributions['chars'].append(sml_char)
                # attributions['pred_logBCF'].append(orig_bcf.item())
                # attributions['pred_logD'].append(orig_lod.item())
 
 
            # print(f"The predicted value={orig_y.item()}")
 
            # print(f"The ground truth logBCF={self.bcf}")
        elif method == "mask":
            for i, sml_char in enumerate(chars):
                sml_copy = chars
                # print(sml_copy)
                sml_copy = "".join(sml_copy[:i] + ["A"] + sml_copy[i+1:])
                # print(sml_copy)
                masked_emb = self.cddd_model.seq_to_emb("".join(sml_copy))
                # masked_x = torch.tensor(masked_emb).float()
                # model.eval()
                masked_y = self.qsar_mdl.predict(masked_emb)
                # orig_y = self.qsar_mdl(torch.tensor(emb).float())
                attributions['attributions'].append((orig_y - masked_y).item())
        # Get the attributions values for the special chars for comparing with correction factors in EPWIN
        special_chars = ['Br', 'Cl', 'I', 'n']
        special_char_idx = [i for i in range(len(chars)) if chars[i] in special_chars]
        # update the attributions dict when there is special atoms
        if special_char_idx:
            special_char_list = [chars[i] for i in special_char_idx]
            special_atom_idx = [i for i in range(len(self.true_char_list)) if self.true_char_list[i] in special_chars]
            attr_special = np.array(attributions['attributions'])[special_char_idx]
            attributions.update({
                'special_char': special_char_list,
                'attr_special': attr_special,
                'special_char_idx': special_char_idx,
                'special_atom_idx': special_atom_idx,
                'SMILES': sml
            })
        return attributions

    def predict(self, sml, return_emb=False, ids=None):
        """
        Interface for predicting logBCF for the given SMILES

        Parameters
        ----------
        sml: str or list of str
            The SMILES of interest for prediction
        return_emb: bool
            If return the generated embedding for the given sml
        """
        # if isinstance(sml, str):
        # the code should work for both str and list
        new_emb = self.cddd_model.seq_to_emb(sml)
        # x = torch.tensor(new_emb).float()
        val = self.qsar_mdl.predict(new_emb)
        # bcf = val[0].detach().numpy().tolist()
        # lod = val[1].detach().numpy().tolist()

        if ids is not None:
            res = {
                    "Compound ID": ids,
                    'prediction': val
                    }
        else:
            res = {
                    'prediction': val
                    }
        if return_emb:
            return res, new_emb
        else:
            return res