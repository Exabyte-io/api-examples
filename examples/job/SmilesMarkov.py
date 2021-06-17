#!/usr/bin/env python
# coding: utf-8

# # SMILES Markov Chains

# In[]:


get_ipython().system('pip install optuna rdkit_pypi deepchem transformers matplotlib tqdm hmmlearn xgboost')


# In[]:


import os
from typing import List
import random
import urllib
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import tqdm
import rdkit.Chem.Descriptors
from hmmlearn.hmm import MultinomialHMM
import xgboost

import rdkit
rdkit.RDLogger.DisableLog('rdApp.*')
from rdkit.Chem.Draw import IPythonConsole
from rdkit.Chem.Draw.MolDrawing import MolDrawing
from rdkit.Chem import Lipinski, Descriptors
import functools


# # Read in our corpus of SMILES strings

# In[]:


dataset_path = '../assets/Molecules_Dataset_Collection/latest'
files = [i for i in os.listdir(dataset_path) if i.endswith('csv')]

all_smiles = set()

for file in files:
    df = pd.read_csv(os.path.join(dataset_path,file))
    all_smiles.update(df.smiles.dropna())
    
all_smiles = set(smile for smile in all_smiles if ':' not in smile)
print(f"total of {len(all_smiles)} smiles")
print("Sample smiles:")
for count,smile in enumerate(all_smiles):
    print(smile)
    if count > 10:
        break


# # Import the SMILES Tokenizer from DeepChem

# In[]:


# Need to patch deepchem first per https://github.com/deepchem/deepchem/issues/2519
import re, importlib
deepchem_root_path=importlib.util.find_spec('deepchem').origin
smiles_tokenizer_path = os.path.join(os.path.dirname(deepchem_root_path), "feat", "smiles_tokenizer.py")
with open(smiles_tokenizer_path, "r") as inp:
  file_contents = inp.read()
with open(smiles_tokenizer_path, "w") as outp:
  for line in file_contents.split("\n"):
    outp.write(re.sub("max_len(\s|$)", "model_max_length", line) + "\n")
    
from deepchem.feat.smiles_tokenizer import SmilesTokenizer


# # Tokenize the Strings

# In[]:


class SmilesEncoder():
    def __init__(self):
        # Download Deepchem's vocabulary file if we don't already have it
        vocab_url = "https://raw.githubusercontent.com/deepchem/deepchem/master/deepchem/feat/tests/data/vocab.txt"
        self.vocab_path = "../assets/deepchem_smiles_vocab.txt"
        if not os.path.exists(self.vocab_path):
            with open(vocab_path, "wb") as outp:
                download = urllib.request.urlopen(vocab_url)
                outp.write(download.read())
                
                # Also, make sure colons are in our vocabulary
                outp.write(b"\n:")
                
        self.tokenizer = SmilesTokenizer(self.vocab_path)
        self.inversion_dict = {value:key for key,value in self.tokenizer.vocab.items()}
        
        self.ignored_range = list(range(self.tokenizer.vocab['[PAD]'], self.tokenizer.vocab['[MASK]']))
        
    def transform(self, X:str) -> List[int]:
        # Don't care about the start/end
        tokenized = self.tokenizer.encode(X)
        return tokenized
    
    def inverse_transform(self, X:List[int]) -> str:
        sep = self.tokenizer.vocab['[SEP]']
        if sep in X:
            end_index = np.random.choice(np.where(X==sep)[0])
            if end_index <= 1:
                end_index = -1
        else:
            end_index = -1
        result = "".join(self.inversion_dict[i] for i in X[1:end_index] if i not in self.ignored_range)
        return result   
enc = SmilesEncoder()


# In[]:


encoded = list(map(enc.transform, all_smiles))
np.random.shuffle(encoded)


# # Tune the Hidden Markov Model

# In[]:


# Here, we're just going to use grid search, assesing the model accuracy via random subsampling from the dataset

def objective(n_states: int) -> int:
    """
    Optimize the number of hidden states. Train on a random subsample from the dataset, since this is rather expensive to run.
    """
    n_subsamples = 16
    subsample_size = 128
    
    total = 0
    for n in tqdm.trange(n_subsamples):
        model = MultinomialHMM(n_states, verbose=False, n_iter=64)
        x=np.concatenate([random.choice(encoded) for _ in range(subsample_size)]).reshape(-1,1)
        model.fit(x)
        total += model.score(x)
    mean = total / n_subsamples
    return mean

vals = []
n_max = 4
for n in range(1,n_max+1):
    print(f"{n} states", flush=True)
    vals.append(objective(n))


# In[]:


plt.plot(list(range(1,len(vals)+1)), vals)
plt.xlabel("Number of Hidden States")
plt.ylabel("Log Loss")
plt.show()
plt.close()


# # Train the Hidden Markov Model

# In[]:


best_n, best_loss = max(enumerate(vals, 1), key=lambda i: i[0])
print(f"Training with {best_n} hidden states")

model = MultinomialHMM(best_n, verbose=True, n_iter=1024)

# This is still a fairly expensive model to train, so we'll restrict ourselves to 10000 molecules in the dataset
model.fit(np.concatenate(encoded).reshape(-1,1))


# In[]:


# Now that we have the model, let's save it to disk (just in case)
import pickle
with open("../assets/smiles_markov_model.pkl", "wb") as outp:
    pickle.dump(model, outp)


# # Create a Molecule Pipeline

# In[]:


def generate_molecules(min_heavy_atoms:int, max_heavy_atoms:int, num_mols:int, verbose:bool=True) -> List[rdkit.Chem.Mol]:
    """
    A pipeline that generates our molecules, from our generator!
    """
    pred_mols = []
    while len(pred_mols) < num_mols:
        sample, _ = model.sample(max_heavy_atoms)
        sample = sample.T[0]
        pred = enc.inverse_transform(sample)

        try:
            mol = rdkit.Chem.MolFromSmiles(pred)
        except:
            mol = None
        if mol and mol.GetNumAtoms() >= min_heavy_atoms:
            pred_mols.append(mol)
            if verbose:
                print(f"Found {len(pred_mols)} molecules. New molecule has weight {np.round((rdkit.Chem.Descriptors.MolWt(mol)),3)} amu")
    return pred_mols


# In[]:


molecules = generate_molecules(min_heavy_atoms=4, max_heavy_atoms=64, num_mols=8, verbose=True)
img = rdkit.Chem.Draw.MolsToGridImage(molecules, subImgSize=[256,256], molsPerRow=4)
img


# # Screen for Bitterness

# In[]:


with open("../assets/bittersweet_xgboost.pkl", "rb") as inp:
    pipeline = pickle.load(inp)
    
# Create a function to do the featurization and prediction
def predict_bitter_from_smiles(mol, model=pipeline):
    # The XGBoost classifier we trained in the last webinar takes in these features in this exact order
    descriptors = (Lipinski.NumAromaticRings,
                   Lipinski.NumHAcceptors,
                   Lipinski.NumHDonors,
                   Descriptors.MolWt,
                   Descriptors.MinPartialCharge,
                   Descriptors.MaxPartialCharge,
                   Descriptors.FpDensityMorgan1,
                   Descriptors.FpDensityMorgan2,
                   Descriptors.FpDensityMorgan3,
                  )
    
    features = [descriptor(mol) for descriptor in descriptors]
    result = model.predict([features])
    return result


# In[]:


bitter_candidates = []
n_desired_candidates = 4
while len(bitter_candidates) < n_desired_candidates:
    bitter_candidates += list(filter(predict_bitter_from_smiles, generate_molecules(min_heavy_atoms=4,
                                                          max_heavy_atoms=30,
                                                          num_mols=16,
                                                          verbose=False)
          ))
img = rdkit.Chem.Draw.MolsToGridImage(bitter_candidates[:n_desired_candidates], subImgSize=[256,256], molsPerRow=4)
img


# Of interest is that the pipeline generated a molecule that does indeed have a bitter taste. The first molecule in the above image is 2-propanol, which has a slightly bitter taste (per the MSDS as https://www2.atmos.umd.edu/~russ/MSDS/isopropanol_fisher.html). 
