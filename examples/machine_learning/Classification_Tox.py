#!/usr/bin/env python
# coding: utf-8

# In[]:


# From https://www.kaggle.com/katyaarnold/bittersweet


# In[]:


import pandas as pd
import sklearn.preprocessing
import sklearn.pipeline
import xgboost 
import optuna
import joblib
import numpy as np
import matplotlib
import matplotlib.pyplot as plt

import functools
import argparse
import warnings

from rdkit.Chem import MolFromSmiles, Lipinski
import rdkit.Chem.Descriptors as Descriptors

import inspect
import re


# In[]:


raw_data = pd.read_csv("../assets/bitter-sweet.csv")


# In[]:


raw_data.Taste.value_counts()


# In[]:


# Drop a few columns we don't want
data = raw_data[["Name", "SMILES"]].copy()

# Let's take a 1 v all approach to this, and just encode "Bitter" as 1 and everything else as 0.
classes = ("Sweet", "Bitter", "Tasteless", "Non-Bitter")
binarized = sklearn.preprocessing.label_binarize(raw_data.Taste.to_numpy(), classes=classes)
bitter = binarized[:,1]
data["Bitter"] = bitter
data = data.dropna()
data


# In[]:


# Let's get a bunch of descriptors into the dataset
descriptors = ((Descriptors.FpDensityMorgan1, "MorganFingerprint Density1"),
               (Descriptors.FpDensityMorgan2, "MorganFingerprint Density2"),
               (Descriptors.FpDensityMorgan3, "MorganFingerprint Density3"),
               (Lipinski.NumAromaticRings, "NumAromaticRings"),
               (Lipinski.NumHAcceptors, "HBondAcceptors"),
               (Lipinski.NumHDonors, "HBondDonors"),
               (Descriptors.MolWt, "MolecularWeight"),
               (Descriptors.MinPartialCharge, "MinPartialCharge"),
               (Descriptors.MaxPartialCharge, "MaxPartialCharge"),
              )

# The dataset is a little unclean, with some badly-formatted SMILEs strings. We'll skip errors there, and just return None instead.
def get_descriptor(smiles, fun):
    try:
        result = fun(MolFromSmiles(smiles))
    except:
        result = None
    return result

for descriptor, name in descriptors:
    data[name] = data.SMILES.apply(functools.partial(get_descriptor, fun=descriptor))

# And finally, we'll drop rows that didn't make it through the featurization
data = data.dropna()


# In[]:


data.head()


# In[]:


data_train, data_test = sklearn.model_selection.train_test_split(data, test_size=0.2)

x = data_train.drop(columns=["Name", "SMILES", "Bitter"]).to_numpy()
n_features = x.shape[1]
y = (data_train.Bitter).to_numpy()


# In[]:


def create_pipeline(params):
    pipeline = sklearn.pipeline.Pipeline([
        ("Standard Scaler", sklearn.preprocessing.StandardScaler()),
        ("XGBoost Classifier", xgboost.XGBClassifier(eval_metric="auc",
                                                    use_label_encoder=False))
    ])
    pipeline.set_params(**params)
    return pipeline
    
def objective(trial: optuna.Trial):
    with joblib.parallel_backend("loky"):
        params = {
            "XGBoost Classifier__eta": trial.suggest_loguniform("XGBoost Classifier__eta", 0.01, 0.4),
            "XGBoost Classifier__min_child_weight": trial.suggest_int("XGBoost Classifier__min_child_weight", 1,100),
            "XGBoost Classifier__max_depth": trial.suggest_int("XGBoost Classifier__max_depth", 3,50),
            "XGBoost Classifier__gamma": trial.suggest_uniform("XGBoost Classifier__gamma", 0, 1),
            "XGBoost Classifier__subsample": trial.suggest_uniform("XGBoost Classifier__subsample", 0, 1),
            "XGBoost Classifier__lambda": trial.suggest_loguniform("XGBoost Classifier__lambda", 0.1, 1000),
            "XGBoost Classifier__alpha": trial.suggest_loguniform("XGBoost Classifier__alpha", 0.1, 1000),
        }

        pipeline = create_pipeline(params)

        score = sklearn.model_selection.cross_val_score(pipeline, x, y, scoring="roc_auc", n_jobs=1, cv=5).mean()
        return score

study = optuna.create_study(direction="maximize")


# In[ ]:


study.optimize(objective, n_trials=100)


# In[]:


study.best_params
pipeline = create_pipeline(study.best_params)
pipeline.fit(x,y)
train_probs = pipeline.predict_proba(x)

error = study.best_value
print(f"5-CV ROC AUC: {error}")


# In[]:


test_preds = pipeline.predict(data_test.drop(columns=["Name", "SMILES", "Bitter"]).to_numpy())
test_probs = pipeline.predict_proba(data_test.drop(columns=["Name", "SMILES", "Bitter"]).to_numpy())
test_score = sklearn.metrics.roc_auc_score(y_true=data_test.Bitter, y_score=test_probs[:,1])
print(f"Test ROC AUC: {test_score}")


# In[]:


fpr, tpr, thresholds = sklearn.metrics.roc_curve(y_true=data_test.Bitter, y_score=test_probs[:,1])
threshes[0] -= 1


# In[]:


# Plot the curve
fig, ax = plt.subplots(dpi=300)
points = np.array([fpr, tpr]).T.reshape(-1, 1, 2)
segments = np.concatenate([points[:-1], points[1:]], axis=1)
norm = plt.Normalize(thresholds.min(), thresholds.max())
lc = matplotlib.collections.LineCollection(segments, cmap='jet', norm=norm, linewidths=2)
lc.set_array(thresholds)
line = ax.add_collection(lc)
fig.colorbar(line, ax=ax).set_label('Threshold')

# Padding to ensure we see the line
ax.margins(0.01)

# plt.plot(false_positive_rate, true_positive_rate, c=colors, label=f"ROC2 Cure, AUC={roc_auc}")
plt.title(f"ROC curve, AUC={np.round(error,3)}")
plt.xlabel("False Positive Rate")
plt.ylabel("True Positive Rate")


# In[ ]:





# In[ ]:




