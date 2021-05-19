#!/usr/bin/env python
# coding: utf-8

# In[]:


import abc 
import numpy as np
import pandas as pd
import sklearn
import sklearn.pipeline
import sklearn.linear_model
import sklearn.decomposition
import sklearn.feature_selection
import inspect
import re
import joblib
import optuna

from rdkit.Chem import MolFromSmiles
import rdkit.Chem.Descriptors as Descriptors
import rdkit.Chem.Crippen as Crippen
import rdkit.Chem.Lipinski as Lipinski
import rdkit.Chem.MolSurf as MolSurf

import matplotlib.pyplot as plt

print(sklearn.__file__)


# In[]:


data = pd.read_csv("logP_dataset.csv",header=None,names=["SMILES", "LogP"])
data.head()


# In[]:


# RDKit supports an older model of LogP. It's computationally-cheap and will probably correlate well with the true LogP.
# Wildman, S. A. and Crippen, G. M. Prediction of Physicochemical Parameters by Atomic Contributions. J. Chem. Inf. Comput. Sci. 1999, 39 (5) 868-873.
# https://doi.org/10.1021/ci990307l
data["crippenLogP"] = data.SMILES.apply(lambda smiles: Crippen.MolLogP(MolFromSmiles(smiles)))
crippen_rmse = np.sqrt(sklearn.metrics.mean_squared_error(data.LogP, data.crippenLogP))

# Plot it
plt.scatter(data.crippenLogP, data.LogP)
plt.title(f"Crippen LogP Model, RMSE = {round(crippen_rmse, 2)}")
plt.xlabel("Crippen LogP")
plt.ylabel("True LogP")

# Make both axes equal
lims = (min(data.crippenLogP.min(), data.LogP.min()), max(data.crippenLogP.max(), data.LogP.min()))
plt.xlim(lims)
plt.ylim(lims)

# Show the plot
plt.show()
plt.close()


# In[]:


# Select a few intuitive descriptors for our dataset
descriptors = [(Descriptors.MolWt, "MolecularWeight"),
               (Descriptors.MinPartialCharge, "MinPartialCharge"),
               (Descriptors.MaxPartialCharge, "MaxPartialCharge"),
               (Descriptors.FpDensityMorgan1, "MorganFingerprint Density1"),
               (Descriptors.FpDensityMorgan2, "MorganFingerprint Density2"),
               (Descriptors.FpDensityMorgan3, "MorganFingerprint Density3"),
               
               (MolSurf.LabuteASA, "SurfaceArea"),
               (Lipinski.FractionCSP3, "FractionSP3Carbons"),
               (Lipinski.HeavyAtomCount, "HeavyAtomCount"),
               (Lipinski.NHOHCount, "NumberNHorOH"),
               (Lipinski.NOCount, "NumberNorO"),
               (Lipinski.NumAliphaticCarbocycles, "NumAliphaticCarbocycles"),
               (Lipinski.NumAliphaticHeterocycles, "NumAliphaticHeterocycles"),
               (Lipinski.NumAliphaticRings, "NumAliphaticRings"),
               (Lipinski.NumAromaticCarbocycles, "NumAromaticCarbocycles"),
               (Lipinski.NumAromaticRings, "NumAromaticRings"),
               (Lipinski.NumHAcceptors, "HBondAcceptors"),
               (Lipinski.NumHDonors, "HBondDonors"),
               (Lipinski.NumHeteroatoms, "NumHeteroatoms"),
               (Lipinski.NumRotatableBonds, "NumRotatableBonds"),
               (Lipinski.NumSaturatedCarbocycles, "NumSaturatedCarbocycles"),
               (Lipinski.NumSaturatedHeterocycles, "NumSaturatedHeterocycles"),
               (Lipinski.NumSaturatedRings, "NumSaturatedRings"),
               (Lipinski.RingCount, "RingCount")
               ]
for descriptor, name in descriptors:
    if name not in data.columns:
        data[name] = data.SMILES.apply(lambda smiles: descriptor(MolFromSmiles(smiles)))
    


# In[]:


data["SmilesSTringLen"] = data.SMILES.apply(lambda smiles: len(smiles))

# Get all the VDW Surface Area descriptors. Since solvation involves this type of long range interaciton, these will likely be important descriptors.
for name, descriptor in filter(lambda result: re.search("PEOE_VSA\d+", result[0]), inspect.getmembers(MolSurf, inspect.isfunction)):
    data[name] = data.SMILES.apply(lambda smiles: descriptor(MolFromSmiles(smiles)))
data.head()


# In[]:


data.describe()


# In[]:


data = data.dropna()
data.describe()


# In[]:


# Let's plot the correlations
plt.close()
cax = plt.matshow(data.corr(), cmap="plasma")
plt.gca().set_xticklabels(data.columns, fontsize=4)
plt.setp(plt.gca().get_xticklabels(), rotation=90)
plt.gca().set_yticklabels(data.columns, fontsize=4)
plt.xticks(range(len(data.columns)))
plt.yticks(range(len(data.columns)))
plt.colorbar(label="Correlation")
plt.rcParams['figure.dpi'] = 400


# In[]:



data_train, data_test = sklearn.model_selection.train_test_split(data, test_size=0.9)

x = data_train.drop(columns=["SMILES", "LogP"]).to_numpy()
n_features = len(data_train.columns) - 2
y = (data_train.LogP).to_numpy()


# In[]:


def create_pipeline(params):
    pipeline = sklearn.pipeline.Pipeline([
        ("Standard Scaler", sklearn.preprocessing.StandardScaler()),
        ("Feature Selection", sklearn.feature_selection.SelectKBest(score_func=sklearn.feature_selection.f_regression)),
        ("Kernelized PCA", sklearn.decomposition.KernelPCA()),
        ("Bayesian Ridge Regression", sklearn.linear_model.BayesianRidge())
    ])
    pipeline.set_params(**params)
    return pipeline
    
def objective(trial: optuna.Trial):
    with joblib.parallel_backend("loky"):
        params = {
            "Feature Selection__k": trial.suggest_int("Feature Selection__k", 3, n_features),
            "Kernelized PCA__n_components": trial.suggest_int("Kernelized PCA__n_components", 1, 3),
            "Kernelized PCA__kernel": trial.suggest_categorical("Kernelized PCA__kernel", ["linear", "poly", "rbf", "sigmoid", "cosine"]),
            "Bayesian Ridge Regression__alpha_1": trial.suggest_loguniform("Bayesian Ridge Regression__alpha_1", 1e-12, 1),
            "Bayesian Ridge Regression__alpha_2": trial.suggest_loguniform("Bayesian Ridge Regression__alpha_2", 1e-12, 1),
            "Bayesian Ridge Regression__lambda_1": trial.suggest_loguniform("Bayesian Ridge Regression__lambda_1", 1e-12, 1),
            "Bayesian Ridge Regression__lambda_2": trial.suggest_loguniform("Bayesian Ridge Regression__lambda_2", 1e-12, 1)
        }

        pipeline = create_pipeline(params)

        score = sklearn.model_selection.cross_val_score(pipeline, x, y, scoring="neg_root_mean_squared_error", n_jobs=1, cv=5).mean()
        return score

study = optuna.create_study(direction="maximize")


# In[]:


study.optimize(objective, n_trials=100, n_jobs=1)


# In[]:


study.best_params
pipeline = create_pipeline(study.best_params)
pipeline.fit(x,y)
error = study.best_value

# Predict Data
data["PredictedLogP"], data["PredictedLogPStd"] = pipeline.predict(data.drop(columns=["SMILES", "LogP"]).to_numpy(), return_std=True)


# In[]:


data = data.sort_values("PredictedLogP")
plt.scatter(data.PredictedLogP, data.LogP)
plt.gca().fill_between(data.PredictedLogP, data.PredictedLogP - data.PredictedLogPStd, data.PredictedLogP + data.PredictedLogPStd,
                        alpha=0.5, color="#FFAAAA")

plt.title(f"LogP, RMSE = {round(-study.best_value, 2)}")
plt.xlabel("Predicted LogP")
plt.ylabel("True LogP")

# Show the plot
plt.show()
plt.close()


# In[]:


data = data.sort_values("crippenLogP")
plt.scatter(data.crippenLogP, data.PredictedLogP)
plt.gca().fill_between(data.crippenLogP, data.PredictedLogP - data.PredictedLogPStd, data.PredictedLogP + data.PredictedLogPStd,
                       alpha=0.5, color="#FFAAAA")

plt.title(f"LogP Model, RMSE = {round(study.best_value, 2)}")
plt.xlabel("Crippen LogP")
plt.ylabel("True LogP")


# Show the plot
plt.show()
plt.close()

