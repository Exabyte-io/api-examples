# ML Datasets

This folder contains datasets which are used in our [machine learning examples](../machine_learning).

- `bitter-sweet.csv` contains a slightly modified version of the dataset from
  https://www.kaggle.com/katyaarnold/bittersweet. We have dropped the `In Bitter Domain` and `Bitter` columns and
  retained the `Taste` column.
- `gamma_alumina_digne_et_al.poscar` contains a structure of gamma alumina published in https://doi.org/10.1016/j.jcat.2004.04.020

## Sci_Adv_Dean_Data.csv

This dataset comes from a recent publication by Dean et al in [Science Advances](https://advances.sciencemag.org/content/5/9/eaax5101.abstract). This work was motivated by the need to predict the adsorption energetics of molecular species to heterogeneous catalysts. Understanding this physical property is important, as catalytic activity is [often a function](https://en.wikipedia.org/wiki/Sabatier_principle) of the adsorption strength of key intermediates, and tuning a catalysts's adsorption strength can be a key part of optimizing its catalytic activity. Typically, one calculates adsorption energy by performing three separate Density-Functional Theory (DFT) calculations:

1. The adsorbate (for example, CO)
2. The catalyst (for example, a nanoparticle)
3. The adsorbed state (for example, CO adsorbed to a particular adsorption site on a nanoparticle)

The adsorption energy is then calculated as follows, where $E_{ads}$ is the adsorption energy, $E_{adsorbate}$ is the energy of the adsorbate (calculation #1 above), $E_{catalyst}$ is the energy of the catalyst (calculation #2 above), and $E_{adsorbate-catalyst}$ is the energy of the complex between the two (calculation #3 above):

$E_{ads} = E_{adsorbate-catalyst} - E_{adsorbate} - E_{catalyst}$

In simpler terms, this equations calculates the energy required to begin at the adsorbed state, and separate the adsorbate and catalyst at an infinite distance from one-another.

This can become computationally expensive, as it requires three different DFT calculations for one material. Hence, this motivates the model created by [Dean et al](https://advances.sciencemag.org/content/5/9/eaax5101.abstract), as the creation of a computationally-inexpensive method of predicting adsorption energetics can increase the throughput of catalyst screening, thus accelerating catalyst discovery.

## The Training Set

The training set consists of the adsorption energetics to sites on various nanoparticles made of copper, silver, and gold (metals which are commonly investigated for catalytic applications). Adsorbates investigated were $\cdot{CH_3}$, $\cdot{OH}$, and $CO$, which were chosen because of their relevance as intermediates in various catalytic cycles. Data was generated using [CP2K](https://www.cp2k.org/)'s implementation of DFT. A Double-Zeta Valence Polarized (DZVP) basis set was used in conjunction with the pseudopotentials of Goedecker, Teter, and Hutter (with a 500Ry cutoff) and the Perdew-Burke-Ernzerhof (PBE) functional. Overall, this results in a robust methodology which is generally good at capturing the adsorption interaction of catalytically-relevant species to transition metals.

## Descriptors Used

In the work, three desriptors were found: $CE_{local}$, $\mu_{adsorbate}$, and MADs. They were chosen for their rapid determination, either being tabulated or computationally-inexpensive to determine.

- $CE_{local}$ is a descriptor which represents the stability of the binding site. Specifically, if we consider a metal atom bound to a nanoparticle, $CE_{local}$ is equal to the sum of bond energies for every bond it forms with the cluster. For example, if atom A was bound to atoms B, C, and D in a nanoparticle, $CE_{local}$ would be equal to the sum of the A-B bond energy, the A-C bond energy, and the A-D bond energy. Implementation-wise, to avoid requiring DFT calculations, Dean et al utilized a variation of the [Bond-Centric Model](https://pubs.acs.org/doi/abs/10.1021/acs.nanolett.8b00670) (BCM) of nanoparticle stability, which enables the rapid prediction of bond energies in a metallic system.

- $\mu_{adsorbate}$ is the chemical potential of the adsorbate, which can be approximated under Hard/Soft Acid/Base (HSAB) theory via a first-order central differencing approximation using the Electron Affinity (EA) and Ionization Potential (IP) of the molecule. This results to the negative average of the IP and EA. For many relevant chemical species, the IP and EA area readily tabulated on resources such as the NIST WebBook or CRC Handbook of Chemistry.

- MADs functions as a descriptor of the intrinsic tendency of the <u>M</u>etal and <u>ADS</u>orbate to bind. It is the binding energy between the adsorbate and a single metal atom of the relevant metal type (for example, for CO adsorbing to Cu, it would be the binding strength of a single CO molecule adsorbing to a single atom of Cu). Although this is performed via DFT, there are few-enough atoms that this calculation is highly tractable. Moreover, because it is only performed once for a given adsorbate-metal pair, it does not present a bottleneck to a high-throughput screening approach.
