#!/usr/bin/env python
# coding: utf-8

# # Background
# 
# This notebook walks the user though the design and tuning of a neural network to predict the adsorption energetics of $\cdot{CH_3}$, $\cdot{OH}$, and $CO$ (several key intermediates in a variety of catalytic processes including $CO_2$ hydrogenation) to nanoparticles made of Cu, Ag, and Au. 
# 
# ## Adsorption Energetics
# 
# Broadly, this notebook re-implements a recent adsorption model found by Dean et al and published in [Science Advances](https://advances.sciencemag.org/content/5/9/eaax5101.abstract). This work was motivated by the need to predict the adsorption energetics of molecular species to heterogeneous catalysts. Understanding this physical property is important, as catalytic activity is [often a function](https://en.wikipedia.org/wiki/Sabatier_principle) of the adsorption strength of key intermediates, and tuning a catalysts's adsorption strength can be a key part of optimizing its catalytic activity. Typically, one calculates adsorption energy by performing three separate Density-Functional Theory (DFT) calculations:
# 
# 1. The adsorbate (for example, CO)
# 2. The catalyst (for example, a nanoparticle)
# 3. The adsorbed state (for example, CO adsorbed to a particular adsorption site on a nanoparticle)
# 
# The adsorption energy is then calculated as follows, where $E_{ads}$ is the adsorption energy, $E_{adsorbate}$ is the energy of the adsorbate (calculation #1 above), $E_{catalyst}$ is the energy of the catalyst (calculation #2 above), and $E_{adsorbate-catalyst}$ is the energy of the complex between the two (calculation #3 above):
# 
# $E_{ads} = E_{adsorbate-catalyst} - E_{adsorbate} - E_{catalyst}$
# 
# In simpler terms, this equations calculates the energy required to begin at the adsorbed state, and separate the adsorbate and catalyst at an infinite distance from one-another.
# 
# This can become computationally expensive, as it requires three different DFT calculations for one material. Hence, this motivates the model created by [Dean et al](https://advances.sciencemag.org/content/5/9/eaax5101.abstract), as the creation of a computationally-inexpensive method of predicting adsorption energetics can increase the throughput of catalyst screening, thus accelerating catalyst discovery.
# 
# ## Descriptors Used
# 
# In the work, three desriptors were found: $CE_{local}$, $\mu_{adsorbate}$, and MADs. They were chosen for their rapid determination, either being tabulated or computationally-inexpensive to determine.
# 
# - $CE_{local}$ is a descriptor which represents the stability of the binding site. Specifically, if we consider a metal atom bound to a nanoparticle, $CE_{local}$ is equal to the sum of bond energies for every bond it forms with the cluster. For example, if atom A was bound to atoms B, C, and D in a nanoparticle, $CE_{local}$ would be equal to the sum of the A-B bond energy, the A-C bond energy, and the A-D bond energy. Implementation-wise, to avoid requiring DFT calculations, Dean et al utilized a variation of the [Bond-Centric Model](https://pubs.acs.org/doi/abs/10.1021/acs.nanolett.8b00670) (BCM) of nanoparticle stability, which enables the rapid prediction of bond energies in a metallic system.
# 
# - $\mu_{adsorbate}$ is the chemical potential of the adsorbate, which can be approximated under Hard/Soft Acid/Base (HSAB) theory via a first-order central differencing approximation using the Electron Affinity (EA) and Ionization Potential (IP) of the molecule. This results to the negative average of the IP and EA. For many relevant chemical species, the IP and EA area readily tabulated on resources such as the NIST WebBook or CRC Handbook of Chemistry.
# 
# - MADs functions as a descriptor of the intrinsic tendency of the <u>M</u>etal and <u>ADS</u>orbate to bind. It is the binding energy between the adsorbate and a single metal atom of the relevant metal type (for example, for CO adsorbing to Cu, it would be the binding strength of a single CO molecule adsorbing to a single atom of Cu). Although this is performed via DFT, there are few-enough atoms that this calculation is highly tractable. Moreover, because it is only performed once for a given adsorbate-metal pair, it does not present a bottleneck to a high-throughput screening approach.
# 
# ## Neural Networks
# 
# In the original paper, a linear regression is found using the above descriptors with a 10-fold cross validated Root Mean Squared Error (RMSE) of 0.179 eV:
# 
# $E_{ads} \approx a*CE_{local} + b*\mu_{adsorbate} + c*MADs + d$
# 
# In this notebook, we replicate the results of this linear regression using a neural network. We proceed through a -step process:
# 
# 1. We import the adsorption data along with the above descriptors descriptor (available in the supporting information of [the paper](https://advances.sciencemag.org/content/5/9/eaax5101.abstract))
# 2. We create a function to programatically create a neural network using TensorFlow
# 3. We optimze several neural network hyperparameters using Scikit-Learn's cross-validated grid search implementation
# 4. We train a final neural network using the tuned hyperparameters, and evaluate it by:
#     a. Learning curves for the train / test set
#     b. A parity plot of predicted versus ground-truth adsorption energy

# # Setting Up

# In[]:


import random
import logging
import os
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import sklearn.model_selection
import subprocess
import tensorflow as tf
import IPython
from tensorflow.keras.layers import Dense, BatchNormalization, InputLayer
from tensorflow.keras.wrappers.scikit_learn import KerasRegressor


# In[]:


# Check if GraphViz is installed; used to plot the NN architecture later on.
rc = subprocess.call(['dot', '-v'])
if rc == 0:
    dot_exists = True
else:
    dot_exists = False


# In[]:


# Random seeds for reproducibility
np.random.seed(1234)
random.seed(1234)


# Next, we'll get everything into a Pandas dataframe.

# In[]:


# Load data into a Pandas dataframe for convenience
df = pd.read_csv("Data.csv", sep=",")
y_cols = ["PBE_BE_eV"]
x_cols = ["CE_Local_eV", "ChemPot_eV", "MADS_eV"]

train_x = df[x_cols]
train_y = df[y_cols]

IPython.display.display(df)


# ## Hyperparameter Optimization
# 
# In this example, we use Scaled-Exponential Linear Unit activation functions, which [provably](https://arxiv.org/pdf/1706.02515.pdf) result in a self-normalizing effect as information flows through a neural network's layers, provided the network is sufficiently deep. One assumption made in the proof is the use of LeCunn's initialization strategy for weights and biases.

# In[]:


# Weight initialization strategy
initializer = tf.keras.initializers.lecun_normal() 


# Next, we create a function to automatically create a TF feedforward neural network with a given width an depth. We then use this function alongside TF's wrapper for the Sklearn regressor API to produce a build function for use in the grid search.

# In[]:


def create_model(hidden_layers_width: int, hidden_layers_depth: int) -> tf.keras.Sequential:
    model = tf.keras.Sequential()
    # Input layer of the model
    model.add(InputLayer(input_shape=len(x_cols),name="Input_Layer"))
    # Batch normalization. Not strictly required since SELU networks self-normalize with sufficient depth. But, we're looking at shallow networks here, so it might be a good idea.
    model.add(BatchNormalization(name="Batch_Normalization"))
    for layer in range(hidden_layers_depth):
        model.add(Dense(hidden_layers_width, activation="selu",
                        kernel_initializer=initializer,
                        name=f"Hidden_Layer_{layer+1}"))
    # Output layer
    model.add(Dense(1, "linear",
                    kernel_initializer=initializer,
                    name="Output_Layer"))
    
    # Compile the model and return it
    model.compile(optimizer="adam", loss="mse")
    return model

# Use TF's wrapper for compliance the Sklearn regressor API
nn_model = KerasRegressor(build_fn=create_model,
                          verbose=0,
                          epochs=256,
                          batch_size=16
                          )


# We now conduct a grid search over several possible architectures.

# In[]:


# Number of folds for cross-validation
folds = 5

# How many threads to use for the grid search
n_threads = 4

# Hidden layer widths to check
hidden_layer_widths = [3]

# Depths to check
hidden_layer_depths = [0, 1, 2, 4, 8]


# In[]:


print(f"Training on {len(train_x)} samples across {len(train_x.columns)} features")
print(f"{folds}-fold CV")

tuning_parameters = {"hidden_layers_width": hidden_layer_widths,
                     "hidden_layers_depth": hidden_layer_depths}
grid = sklearn.model_selection.GridSearchCV(estimator=nn_model,
                                            param_grid=tuning_parameters,
                                            scoring="neg_mean_squared_error",
                                            verbose=1,
                                            cv=folds,
                                            n_jobs=n_threads
                                            )
logger = tf.get_logger()
logger.setLevel(logging.ERROR)

grid.fit(X=train_x, y=train_y)


# Now, let's print out some information about the best architecture found:

# In[]:


# Print the RMSE
rmse = np.round(np.sqrt(abs(grid.best_score_)),3)
print(f"RMSE = {rmse}")

# Print the architecture
best_params = grid.best_params_
print("\nBest parameters are:")
for key, value in best_params.items():
    print(f"\t{key} : {value}")

# Compile the best model, and print information about it
best_model = create_model(hidden_layers_width=best_params["hidden_layers_width"],
                          hidden_layers_depth=best_params["hidden_layers_depth"])

print("\nDetailed Information:")
print(best_model.summary())

print("\nArchitecture:")
if dot_exists:
    arch_filename = "architecture.png"
    tf.keras.utils.plot_model(best_model, show_shapes=True, to_file=arch_filename)
    image = IPython.display.Image(filename=arch_filename)
    IPython.display.display(image)
else:
    print("dot is not found in the System path, but was needed to draw the architecture.")


# ## Final Training
# 
# Now that we have identified a good model architecture, we train it.

# In[]:


# Train the final network
history = best_model.fit(train_x, train_y, epochs=100000, batch_size=16,
                         verbose=0,
                         validation_split=0.1,
                         callbacks=[
                             tf.keras.callbacks.EarlyStopping(
                                 monitor="val_loss",
                                 min_delta=0,
                                 patience=50,
                                 mode="min",
                                 restore_best_weights=True)
                         ]
                         ).history


# ## Final Evaluation
# 
# Now that we have a trained model, we'll plot its learning curve and create a parity plot.

# In[]:


plt.rcParams['figure.figsize'] = 5,5
plt.rcParams['figure.dpi'] = 200


# In[]:


# Plot the training history
loss = np.sqrt(history["loss"])
val_loss = np.sqrt(history["val_loss"])
plt.plot(loss, label="Training Set Loss")
plt.plot(val_loss, label="Validation Set Loss")
plt.legend()
plt.ylabel("RMSE Loss (eV)")
plt.xlabel("Epoch Number")
fig = plt.gcf()
fig.savefig("learning_curve.png")
plt.show()
plt.close()


# In[]:


# Parity plot
df["predicted_adsorption_eV"] = best_model.predict(train_x)

plt.scatter(x=df["predicted_adsorption_eV"], y=train_y,
            label=f"Predictions (RMSE = {np.round(rmse, 2)} eV)",
            c="#1c2957")

lims_min = min(min(df["predicted_adsorption_eV"]),
               float(train_y.min())
               )
lims_max = max(max(df["predicted_adsorption_eV"]),
                float(train_y.max())
               )
lims = (lims_min, lims_max)

plt.xlim(lims)
plt.ylim(lims)
plt.plot([lims[0], lims[1]],
         [lims[0], lims[1]],
         label="Parity", c="black")
plt.legend()
plt.ylabel("DFT Adsorption Energy (eV)")
plt.xlabel("Predicted Adsorption Energy (eV)")

fig = plt.gcf()
fig.savefig("parity.png")
plt.show()
plt.close()


# ## Save Model
# 
# And last, we will save the final trained model.

# In[]:


# Save the final model
print(best_params)
best_model.save("best_model.mdl")

