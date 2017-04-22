Project 1: Learning Galaxy Formation from Large-scale Structure and Weak Lensing
======================

Project Lead: Ying Zu
--------


# Program per session (2 hours each , total 8 hours):

## 1. Fitting Halo Mass to Weak Lensing Profiles.

* Introduce the weak lensing measurement of Mandelbaum+2016 on the average halo mass of galaxies at fixed
stellar mass split by colors, <M_h|M_*r>_red vs. <M_h|M_*>_blue. (~20 mins)

* Build a model for NFW profiles at fixed M_h and concentration. (~30 mins)

* Vary parameters in the NFW profile on a grid and overplot the results against the measurements
from i. (~20 mins)

* Fit individual weak lensing profiles using simple curve fitting. (~20 mins)

* Fit individual weak lensing profiles using MCMC sampling, and compare the results to iv. (~20 mins)

## 2. Predicting Halo Mass from Theory.

* Interpret the results from session 1 and describe the steps for predicting them from theory. (~20 mins)

* Measure halo mass function (HMF) and stellar-to-halo mass relation (SHMR) from simulated halos and mock
galaxies (provided by me) (~30 mins)

* Predict HMF and SHMR from theory (plenty of HMF code online!). (~30 mins)

* Convolve HMF and SMHR to compute <M_h|M_*>. (~30 mins)

## 3. Constructing Galaxy Quenching Models.

*  Introduce the halo quenching model of Zu & Mandelbaum (2016). (~20 mins)

* Implement a quenching model into the <M_h|M_*> code from session #2, and predict <M_h|M_*r>_red vs.
    <M_h|M_*>_blue. (~30 mins)

* Measure <M_h|M_*r>_red vs. <M_h|M_*>_blue from the mock catalog. (~20 mins)

* Vary parameters in the quenching model on a grid and overplot the results against the measurements
    from iii. (~30 mins)


## 4. Model inference.

* Write the likelihood functions and run MCMC inferences (~ 30 mins).

* Try different quenching models and repeat the analysis. (~ 30 mins)

* Comment on model selection, quenching of satellite population, extra constraints from correlation
    functions, etc.

# Requirements (software, bibliography):

## software:

* Python (including scipy and numpy)
* h5py
* your favorite HMF and cosmology code (type "halo mass function" into github searchbar).

## biblio:

* Mandelbaum et al. 2016 (http://adsabs.harvard.edu/abs/2016MNRAS.457.3200M)
* Zu and Mandelbaum 2015 (http://adsabs.harvard.edu/abs/2015MNRAS.454.1161Z)
* Zu and Mandelbaum 2016 (http://adsabs.harvard.edu/abs/2016MNRAS.457.4360Z)

## Are there any concepts (cosmology/statistics) that you would like the students review previous to your sessions?
Basic MCMC concepts, basic large scale structure

* None. Just ask questions during the sessions.


