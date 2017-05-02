# 02 May 2017 16:39:54
import sys
sys.path.insert(1, '../session_2/')
try:
    from fred import RedFractionGenerator
    from zypy.zycosmo import get_halofuncs, get_cosmo, CosmoParams
    from hod.shmr import SHMR_Leauthaud, HSMR, HSMR_split
    from hod.predict import get_f_sigma_lnMs
    has_fred = True
except ImportError:
    has_fred = False
import numpy as np
import matplotlib.pyplot as plt
from read_mock import read_mock



"""Implement Halo Quenching."""


def is_red(lgms, gcolor):
    cut = 0.8*(lgms/10.5)**0.6
    isred = gcolor >= cut
    return(isred)

def test_mock_hsmr_split(mockfile):
    """Check the color-split version of halo to stellar mass relations in the mock."""
    galrec = read_mock(mockfile)
    iscen = galrec['lg_halo_mass'] > 1
    lgmh = galrec['lg_halo_mass'][iscen]
    lgms = galrec['lg_stellar_mass'][iscen]
    gcolor = galrec['g-r:'][iscen]
    lgms_bins = np.linspace(9.8, 12.0, 9)
    lgms_cens = (lgms_bins[1:] + lgms_bins[:-1]) / 2.0
    lgmh_cens_red = np.empty_like(lgms_cens)
    lgmh_err_red = np.empty_like(lgms_cens)
    lgmh_cens_blue = np.empty_like(lgms_cens)
    lgmh_err_blue = np.empty_like(lgms_cens)
    for i in xrange(lgms_cens.size):
        sel = (lgms >= lgms_bins[i]) & (lgms < lgms_bins[i+1])
        isred = is_red(lgms[sel], gcolor[sel])
        nred = float(np.sum(isred))
        isblue = ~isred
        nblue = float(np.sum(isblue))
        lgmh_cens_red[i] = np.mean(lgmh[sel][isred])
        lgmh_err_red[i] = np.std(lgmh[sel][isred])/np.sqrt(nred)
        lgmh_cens_blue[i] = np.mean(lgmh[sel][isblue])
        lgmh_err_blue[i] = np.std(lgmh[sel][isblue])/np.sqrt(nblue)
    if has_fred:
        pass
        #
    plt.errorbar(lgms_cens, lgmh_cens_red,  yerr=lgmh_err_red, marker="o", ms=5, color="r")
    plt.errorbar(lgms_cens, lgmh_cens_blue,  yerr=lgmh_err_blue, marker="s", ms=5, color="b")
    plt.xlabel(r"$M_*\;[M_\odot/h^2]$")
    plt.ylabel(r"$M_h\;[M_\odot/h]$")
    plt.show()

if __name__ == "__main__":
    mockfile = '/Users/ying/Dropbox/Public/iHODcatalog_bolshoi.h5'
    test_mock_hsmr_split(mockfile)

