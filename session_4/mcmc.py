# 05 May 2017 03:11:46

# import sys
# sys.path.insert(1, '../session_3/')
import emcee
import numpy as np
import matplotlib.pyplot as plt

"""Running MCMC Sampling"""

def read_zm16_data(use_red=True):
    """read mock data"""
    if use_red:
        fname = '../data/ZM16/halo_mass_red.dat'
    else:
        fname = '../data/ZM16/halo_mass_blue.dat'
    lgms, lgmh, errlgmh = np.genfromtxt(fname, unpack=True)
    return(lgms, lgmh, errlgmh)

def predict_halomass_bimodality(p):
    """p: input quenching parameters"""
    pass

def loglikelihood():
    pass

def logposteriror():
    pass

def run_mcmchammer():
    pass


if __name__ == "__main__":
    pass
