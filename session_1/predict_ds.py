# 01 May 2017 18:12:10
try:
    from zypy.zycosmo import Density, NFW, get_cosmo
except ImportError:
    print("Cosmology/NFW code is required.")
from read_m16 import read_m16
import matplotlib.pyplot as plt
import numpy as np


def predict_ds(rp, mass=1e13, c=5, rho_mean=52700242579.0, DELTA_HALO=200, h=0.7):
    nfw = NFW(mass, c, rho_mean, DELTA_HALO)
    ds = nfw.DelSig(rp/h)
    # [Msun/Mpc^2] to [h Msun/pc^2]
    ds = ds / 1e12 / h
    return(ds)

def get_density(z=0.1):
    cosmo = get_cosmo('PLANCK2')
    den = Density(cosmo)
    # Density unit in solar masses per cubic (physical) Mpc.
    rho_mean = den.rho_mean_z(z)
    return(rho_mean)

def test_predict_ds(mass_bin='11.0_15.0'):
    for use_red in [True, False]:
        if use_red:
            label = "M16 Red"
            color = 'red'
            marker = "o"
        else:
            label = "M16 Blue"
            color = 'blue'
            marker = "s"
        rp, ds, ds_err = read_m16(use_red=use_red, mass_bin=mass_bin)
        plt.errorbar(rp, ds, yerr=ds_err, marker=marker, ms=5, ls="--", lw=0.5, color=color, label=label)
    _rp = np.logspace(-1.5, 0.8, 20)
    lgmass_list = [11.5, 12, 12.5, 13, 13.5, 14]
    nlgmass = len(lgmass_list) * 1.0
    for i, lgmass in enumerate(lgmass_list):
        _ds =  predict_ds(_rp, 10**lgmass)
        plt.plot(_rp, _ds, ls="-", label=format(lgmass, '3.1f'), color=format((nlgmass-float(i))/(nlgmass+1.0), '2.1f'))
    plt.legend(loc=1)
    plt.xlabel(r"$R\;[Mpc/h]$")
    plt.ylabel(r"$\Delta\Sigma\;[h M_\odot/pc^2]$")
    plt.xscale('log')
    plt.yscale('log')
    plt.show()

if __name__ == "__main__":
    # print get_density()
    test_predict_ds(mass_bin='11.0_15.0')
    test_predict_ds(mass_bin='10.4_10.7')
    test_predict_ds(mass_bin='10.7_11.0')
