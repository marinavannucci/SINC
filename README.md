# Simultaneous Inference of Networks and their Covariates (SINC)
SINC Algorithm from the paper "Latent Network Estimation and Variable Selection for Compositional Data via Variational EM"

Python code and simulations used in the above mentioned paper. The following items can be found in this repository

## SINC_functions.py

This file contains the SINC function and all other dependent functions to run SINC. These are dependent on the python modules numpy, scipy, multiprocess, and multiprocessing.

## small_example

A folder containing a small simulation tutorial on how the function may be used. In this folder you will find small simulation.ipynb, a notebook showing the models use, with supplementing synthetic data saved as csv files.

## Simulations

The simulations reported in section 4 of the paper can be found in this folder, labeled as Band, Cluster, Hub, and Random, corresponding to the graph structure.