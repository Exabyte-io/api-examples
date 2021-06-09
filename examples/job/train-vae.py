#!/usr/bin/env python
# coding: utf-8

# In[]:


get_ipython().system('pip install tensorflow rdkit_pypi pymonad tqdm deepchem')


# In[]:


import urllib.request
import itertools
import os

import tensorflow as tf
import pandas as pd
import numpy as np

from deepchem.feat.smiles_tokenizer import SmilesTokenizer

import rdkit.Chem
import tqdm

RANDOM_SEED = 42
tqdm.tqdm.pandas()


# # Read the Dataset

# In[]:


data = pd.read_csv("assets/curated-solubility-dataset.csv")

# Drop duplicates based on the InChI representation
nondupes = data.drop_duplicates(subset="InChI")

smiles = nondupes.SMILES
smiles


# # Data Augmentation

# In[ ]:


augmented_data = list(smiles.copy())
print(f"Starting with a dataset of length {len(augmented_data)}")
total_augmentations = 100
print(f"Augmenting with {total_augmentations} per string and dropping duplicates/malformed SMILES", flush=True)
def augment_smiles(smiles_string: str):
    try:
        mol = rdkit.Chem.MolFromSmiles(smiles_string)
    except (RuntimeError, TypeError):
        return None
    results = []
    try:
        results += [rdkit.Chem.MolToRandomSmilesVect(mol, total_augmentations, RANDOM_SEED)]
    except (RuntimeError, TypeError):
        pass
    return tuple(*results)

with tqdm.tqdm(total=total_augmentations*len(augmented_data)) as pbar:
    for smile in smiles:
        if augmentation:=augment_smiles(smile):
            augmented_data += augmentation
        pbar.update(total_augmentations)

augmented_data = pd.Series(list(set(augmented_data)))
print(f"Finishing with a dataset of length {len(augmented_data)}")


# # Tokenization with DeepChem

# In[ ]:


# DeepChem's SmileTokenizer uses the WordPiece transformer by HuggingFace (https://huggingface.co/transformers/tokenizer_summary.html), with the regular expression SMILES tokenization strategy developed by Schwaller, P. et al in https://doi.org/10.1039/c8sc02339e
from deepchem.feat.smiles_tokenizer import SmilesTokenizer

# Download Deepchem's vocabulary file if we don't already have it
vocab_url = "https://raw.githubusercontent.com/deepchem/deepchem/master/deepchem/feat/tests/data/vocab.txt"
vocab_path = "assets/deepchem_smiles_vocab.txt"
if not os.path.exists(vocab_path):
    with open(vocab_path, "wb") as outp:
        download = urllib.request.urlopen(vocab_url)
        outp.write(download.read())
        
tokenizer = SmilesTokenizer(vocab_path)
tokenized_smiles = augmented_data.progress_map(tokenizer.encode)
tokenized_smiles


# In[ ]:


# Next up, we'll remove SMILES with any unknown characters, since we don't want our generator putting those tokens in the output
tokenized_smiles = tokenized_smiles[tokenized_smiles.apply(lambda smile: tokenizer.vocab["[UNK]"] not in smile)]

# And we'll 0-pad the strings out to a constant length
maxlen = tokenized_smiles.apply(len).max()
maxlen += maxlen%2
tokenized_smiles = tokenized_smiles.apply(lambda smile: smile + [tokenizer.vocab["[PAD]"]]*(maxlen-len(smile)))
tokenized_smiles.apply(len).describe()


# # Create the VAE

# In[ ]:


from tensorflow.keras.layers import Input, Dense, Conv1D, Layer, Flatten, Reshape, Conv1DTranspose


# In[ ]:


class Sampling(Layer):
    """
    Taken from https://keras.io/examples/generative/vae/
    """
    def call(self, inputs):
        z_mean, z_log_var = inputs
        batch = tf.shape(z_mean)[0]
        dim = tf.shape(z_mean)[1]
        epsilon = tf.keras.backend.random_normal(shape=(batch, dim))
        return z_mean + tf.exp(0.5 * z_log_var) * epsilon


# In[ ]:


latent_dim = 1

encoder_inputs = Input(shape=(maxlen, 1))
h = Conv1D(32, 3, activation='relu', strides=2, padding='same')(encoder_inputs)
h = Conv1D(64, 3, activation='relu', strides=2, padding='same')(h)
h = Flatten()(h)
h = Dense(16, activation='relu')(h)
z_mean = Dense(latent_dim, name='z_mean')(h)
z_log_var = Dense(latent_dim, name='z_log_var')(h)
z = Sampling()([z_mean, z_log_var])

encoder = tf.keras.Model(encoder_inputs, [z_mean, z_log_var, z], name='encoder')
encoder.summary()


# In[ ]:


latent_inputs = Input(shape=(latent_dim,))
h = Dense(78*64)(latent_inputs)
h = Reshape((78,64))(h)
h = Conv1DTranspose(64, 2, activation="relu", strides=2, padding="same")(h)
h = Conv1DTranspose(32, 2, activation="relu", strides=2, padding="same")(h)
decoder_outputs = Conv1DTranspose(1,3,activation="sigmoid", padding="same")(h)

decoder = tf.keras.Model(latent_inputs, decoder_outputs, name="decoder")
decoder.summary()


# In[ ]:


class VAE(tf.keras.Model):
    def __init__(self, encoder, decoder, **kwargs):
        super(VAE, self).__init__(**kwargs)
        self.encoder = encoder
        self.decoder = decoder
        self.total_loss_tracker = tf.keras.metrics.Mean(name="total_loss")
        self.reconstruction_loss_tracker = tf.keras.metrics.Mean(
            name="reconstruction_loss"
        )
        self.kl_loss_tracker = tf.keras.metrics.Mean(name="kl_loss")

    @property
    def metrics(self):
        return [
            self.total_loss_tracker,
            self.reconstruction_loss_tracker,
            self.kl_loss_tracker,
        ]

    def train_step(self, data):
        with tf.GradientTape() as tape:
            z_mean, z_log_var, z = self.encoder(data)
            reconstruction = self.decoder(z)
            reconstruction_loss = tf.reduce_mean(
                tf.reduce_sum(
                    tf.keras.losses.binary_crossentropy(data, reconstruction), axis=[1]
                )
            )
            kl_loss = -0.5 * (1 + z_log_var - tf.square(z_mean) - tf.exp(z_log_var))
            kl_loss = tf.reduce_mean(tf.reduce_sum(kl_loss, axis=0))
            total_loss = reconstruction_loss + kl_loss
        grads = tape.gradient(total_loss, self.trainable_weights)
        self.optimizer.apply_gradients(zip(grads, self.trainable_weights))
        self.total_loss_tracker.update_state(total_loss)
        self.reconstruction_loss_tracker.update_state(reconstruction_loss)
        self.kl_loss_tracker.update_state(kl_loss)
        return {
            "loss": self.total_loss_tracker.result(),
            "reconstruction_loss": self.reconstruction_loss_tracker.result(),
            "kl_loss": self.kl_loss_tracker.result(),
        }


# # Train the VAE

# In[ ]:


# Reshape the data as needed, scale between 0 and 1
train_data = np.array([i for i in tokenized_smiles])
train_data = np.expand_dims(train_data, -1).astype("float32") / maxlen
train_data.reshape(1,312,-1).shape


# In[ ]:


vae = VAE(encoder, decoder)
vae.compile(optimizer=tf.keras.optimizers.Adam())
vae.fit(train_data, epochs=1000, batch_size=64)


# In[ ]:





# # Detokenization
# Todo: Move this after the VAE stuff

# In[ ]:


# Create a method to detokenize the encoded SMILES string
detokenize_dict = {value:key for key,value in tokenizer.vocab.items()}
def invert_tokenization(tokens):
    inverted = [detokenize_dict[token] for token in tokens]
    partial_smiles = "".join(inverted)
    
    # Remove things we can't inverse transform
    ignored = [detokenize_dict[i] for i in range(tokenizer.vocab['[PAD]'], tokenizer.vocab['[MASK]']+1)]
    for to_remove in ignored:
        partial_smiles = partial_smiles.replace(to_remove, "")
    return partial_smiles


# In[ ]:


predictions = (vae.decoder.predict([np.random.uniform(low=-10, high=10, size=100)])[:,:,0] * maxlen).astype(int)
list(map(lambda tokenized: invert_tokenization(tokenized), predictions))


# In[ ]:




